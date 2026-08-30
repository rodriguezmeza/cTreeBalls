#!/usr/bin/env python3
"""Compare the cTreeBalls radial and three-dimensional Ly-alpha estimators.

The default synthetic catalog is strictly collinear.  In that geometry the
3D 2PCF can be projected by summing its numerator and denominator over the
transverse bins, and the 3D 3PCF can be projected onto the two signed radial
lags.  Those projected histograms must agree with the native 1D estimators.

For a non-collinear catalog the estimators answer different questions: the 1D
methods deliberately ignore transverse separation, while the 3D methods do
not.  This example therefore rejects non-collinear input instead of presenting
that expected physical difference as a numerical error.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import asdict, dataclass
import json
import math
import os
from pathlib import Path
import subprocess
import time
from typing import Dict, Iterable, Optional

import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parents[1]


@dataclass
class ComparisonConfig:
    output_dir: Path
    cballs_executable: Path = PROJECT_ROOT / "cballs"
    catalog_path: Optional[Path] = None
    n_pixels: int = 160
    pixels_per_forest: int = 4
    radial_min: float = 50.0
    radial_max: float = 80.0
    threads: int = min(4, os.cpu_count() or 1)
    seed: int = 20260830
    rp_max: float = 12.0
    rp_bins: int = 12
    rt_max: float = 2.0
    rt_bins: int = 4
    r3_max: float = 8.0
    r3_bins: int = 8
    theta_bins: int = 4
    mu_bins: int = 4
    include_3pcf: bool = True
    include_combined: bool = True

    def validate(self) -> None:
        self.output_dir = Path(self.output_dir).expanduser().resolve()
        self.cballs_executable = Path(
            self.cballs_executable
        ).expanduser().resolve()
        if self.catalog_path is not None:
            self.catalog_path = Path(self.catalog_path).expanduser().resolve()
        if not self.cballs_executable.is_file():
            raise ValueError(f"cballs executable not found: {self.cballs_executable}")
        if self.catalog_path is not None and not self.catalog_path.is_file():
            raise ValueError(f"catalog not found: {self.catalog_path}")
        if self.n_pixels < 3 or self.pixels_per_forest < 1:
            raise ValueError("n_pixels must be at least 3 and pixels_per_forest positive")
        if self.threads < 1:
            raise ValueError("threads must be positive")
        if not math.isfinite(self.radial_min) or not math.isfinite(self.radial_max):
            raise ValueError("radial limits must be finite")
        if self.radial_min <= 0.0 or self.radial_max <= self.radial_min:
            raise ValueError("require 0 < radial_min < radial_max")
        maxima = (self.rp_max, self.rt_max, self.r3_max)
        if not all(math.isfinite(value) and value > 0.0 for value in maxima):
            raise ValueError("all separation maxima must be finite and positive")
        bins = (self.rp_bins, self.rt_bins, self.r3_bins, self.mu_bins)
        if any(value < 1 for value in bins):
            raise ValueError("all histogram bin counts must be positive")
        if self.include_3pcf and self.theta_bins < 2:
            raise ValueError("theta_bins must be at least 2 to distinguish signed lags")


@dataclass
class Histogram:
    numerator: np.ndarray
    denominator: np.ndarray

    @property
    def value(self) -> np.ndarray:
        result = np.zeros_like(self.numerator, dtype=np.float64)
        np.divide(
            self.numerator,
            self.denominator,
            out=result,
            where=self.denominator != 0.0,
        )
        return result


def generate_collinear_catalog(config: ComparisonConfig) -> np.ndarray:
    """Return ``x y z delta weight forest_id`` on one line of sight."""
    rng = np.random.default_rng(config.seed)
    radii = rng.uniform(config.radial_min, config.radial_max, config.n_pixels)
    order = rng.permutation(config.n_pixels)
    radii = radii[order]
    delta = (
        0.22 * np.sin(0.43 * radii)
        + 0.11 * np.cos(0.17 * radii)
        + rng.normal(0.0, 0.08, config.n_pixels)
    )
    weights = rng.uniform(0.65, 1.35, config.n_pixels)
    forest_ids = np.arange(config.n_pixels) // config.pixels_per_forest
    forest_ids = forest_ids[order]
    return np.column_stack(
        (radii, np.zeros(config.n_pixels), np.zeros(config.n_pixels),
         delta, weights, forest_ids)
    )


def load_catalog(path: Path) -> np.ndarray:
    catalog = np.loadtxt(path, dtype=np.float64, comments="#", ndmin=2)
    if catalog.shape[1] != 6:
        raise ValueError(
            f"{path} must contain x y z delta weight forest_id; "
            f"found {catalog.shape[1]} columns"
        )
    return catalog


def validate_collinear_catalog(catalog: np.ndarray, tolerance: float = 2e-12) -> None:
    if catalog.ndim != 2 or catalog.shape[1] != 6 or len(catalog) < 3:
        raise ValueError("catalog must have at least three rows and six columns")
    if not np.all(np.isfinite(catalog)):
        raise ValueError("catalog contains a non-finite value")
    positions = catalog[:, :3]
    radii = np.linalg.norm(positions, axis=1)
    if np.any(radii <= 0.0):
        raise ValueError("all catalog positions must be nonzero")
    directions = positions / radii[:, None]
    alignment = directions @ directions[0]
    if np.max(np.abs(alignment - 1.0)) > tolerance:
        raise ValueError(
            "exact 1D/3D projection requires every point on the same "
            "observer-centered ray"
        )
    if np.any(catalog[:, 4] <= 0.0):
        raise ValueError("weights must be positive")
    if np.max(np.abs(catalog[:, 5] - np.rint(catalog[:, 5]))) > tolerance:
        raise ValueError("forest_id values must be integers")


def validate_3pcf_projection_bins(
    catalog: np.ndarray, config: ComparisonConfig, tolerance: float = 2e-12
) -> None:
    """Reject the measure-zero cases not recoverable from a binned 3D file."""
    radii = np.linalg.norm(catalog[:, :3], axis=1)
    separations = np.abs(radii[:, None] - radii[None, :])
    forest_ids = np.rint(catalog[:, 5]).astype(np.int64)
    upper = np.triu_indices(len(radii), k=1)
    separations = separations[upper]
    distinct_forests = forest_ids[upper[0]] != forest_ids[upper[1]]
    active = distinct_forests & (separations < config.r3_max)
    scaled = separations[active] * config.r3_bins / config.r3_max
    if np.any(np.isclose(scaled, np.rint(scaled), rtol=0.0, atol=tolerance)):
        raise ValueError(
            "a 3PCF radial separation lies on a bin boundary; the signed 1D "
            "bin cannot be reconstructed uniquely from the binned 3D output"
        )


def write_catalog(path: Path, catalog: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(
        path,
        catalog,
        fmt=("%.17g",) * 5 + ("%.0f",),
        header="x y z delta weight forest_id",
    )


def thread_environment(threads: int) -> Dict[str, str]:
    env = os.environ.copy()
    env.update(
        {
            "OMP_NUM_THREADS": str(threads),
            "OMP_DYNAMIC": "FALSE",
            "OPENBLAS_NUM_THREADS": "1",
            "MKL_NUM_THREADS": "1",
            "VECLIB_MAXIMUM_THREADS": "1",
            "NUMEXPR_NUM_THREADS": "1",
        }
    )
    return env


def run_method(
    config: ComparisonConfig, catalog_path: Path, method: str
) -> tuple[Path, float]:
    root = config.output_dir / "runs" / method
    root.mkdir(parents=True, exist_ok=True)
    for stale in root.glob("hist*lya*.txt"):
        stale.unlink()
    command = [
        os.fspath(config.cballs_executable),
        f"search={method}",
        f"infile={catalog_path}",
        "infileformat=lya-ascii",
        "iCatalogs=1",
        f"rootDir={root}",
        f"numberThreads={config.threads}",
        "usePeriodic=false",
        "useLogHist=false",
        f"rangeN={max(config.rp_max, config.r3_max) * 2.5:.17g}",
        "rminHist=0.1",
        f"sizeHistN={config.rp_bins}",
        f"lya2RpMax={config.rp_max:.17g}",
        f"lya2RtMax={config.rt_max:.17g}",
        f"lya2RpBins={config.rp_bins}",
        f"lya2RtBins={config.rt_bins}",
        f"lya3RMax={config.r3_max:.17g}",
        f"lya3RBins={config.r3_bins}",
        f"lya3ThetaBins={config.theta_bins}",
        f"lya3MuBins={config.mu_bins}",
        "verbose=0",
        "verbose_log=0",
        "options=",
    ]
    start = time.perf_counter()
    completed = subprocess.run(
        command,
        cwd=PROJECT_ROOT,
        env=thread_environment(config.threads),
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    elapsed = time.perf_counter() - start
    (root / "run.log").write_text(completed.stdout, encoding="utf-8")
    if completed.returncode:
        raise RuntimeError(
            f"{method} failed with status {completed.returncode}; "
            f"see {root / 'run.log'}\n{completed.stdout[-3000:]}"
        )
    return root, elapsed


def data_rows(path: Path) -> Iterable[list[str]]:
    if not path.is_file():
        raise FileNotFoundError(f"expected output was not produced: {path}")
    for line in path.read_text(encoding="utf-8").splitlines():
        if line and not line.startswith("#"):
            yield line.split()


def read_1d_2pcf(path: Path, bins: int) -> Histogram:
    numerator = np.zeros(bins)
    denominator = np.zeros(bins)
    for fields in data_rows(path):
        index = int(fields[0])
        numerator[index] = float(fields[3])
        denominator[index] = float(fields[4])
    return Histogram(numerator, denominator)


def read_projected_3d_2pcf(path: Path, rp_bins: int, rt_bins: int) -> Histogram:
    numerator = np.zeros((rp_bins, rt_bins))
    denominator = np.zeros((rp_bins, rt_bins))
    for fields in data_rows(path):
        bp, bt = map(int, fields[:2])
        numerator[bp, bt] = float(fields[5])
        denominator[bp, bt] = float(fields[6])
    return Histogram(numerator.sum(axis=1), denominator.sum(axis=1))


def read_1d_3pcf(path: Path, bins_per_side: int) -> Histogram:
    total_bins = 2 * bins_per_side
    numerator = np.zeros((total_bins, total_bins))
    denominator = np.zeros((total_bins, total_bins))
    for fields in data_rows(path):
        b1, b2 = map(int, fields[:2])
        numerator[b1, b2] = float(fields[5])
        denominator[b1, b2] = float(fields[6])
    return Histogram(numerator, denominator)


def signed_bin_from_3d(radial_bin: int, theta_bin: int,
                       bins_per_side: int, theta_bins: int) -> int:
    if theta_bin == 0:
        return bins_per_side + radial_bin
    if theta_bin == theta_bins - 1:
        return bins_per_side - 1 - radial_bin
    raise ValueError(
        "occupied intermediate theta bin: the input is not collinear at "
        "the requested numerical precision"
    )


def read_projected_3d_3pcf(
    path: Path, bins_per_side: int, theta_bins: int, mu_bins: int
) -> Histogram:
    total_bins = 2 * bins_per_side
    numerator = np.zeros((total_bins, total_bins))
    denominator = np.zeros((total_bins, total_bins))
    for fields in data_rows(path):
        b1, b2, t1, t2, bmu = map(int, fields[:5])
        row_numerator = float(fields[11])
        row_denominator = float(fields[12])
        if row_denominator == 0.0:
            continue
        s1 = signed_bin_from_3d(b1, t1, bins_per_side, theta_bins)
        s2 = signed_bin_from_3d(b2, t2, bins_per_side, theta_bins)
        expected_mu = mu_bins - 1 if (s1 < bins_per_side) == (
            s2 < bins_per_side
        ) else 0
        if bmu != expected_mu:
            raise ValueError(
                "occupied opening-angle bin is inconsistent with a collinear catalog"
            )
        numerator[s1, s2] += row_numerator
        denominator[s1, s2] += row_denominator
    return Histogram(numerator, denominator)


def relative_difference(reference: np.ndarray, candidate: np.ndarray) -> np.ndarray:
    scale = np.maximum.reduce(
        (np.abs(reference), np.abs(candidate), np.full(reference.shape, 1e-14))
    )
    return np.abs(candidate - reference) / scale


def comparison_metrics(reference: Histogram, candidate: Histogram) -> dict:
    occupied = (reference.denominator != 0.0) | (candidate.denominator != 0.0)
    if not np.any(occupied):
        return {
            "occupied_bins": 0,
            "max_absolute": 0.0,
            "max_relative": 0.0,
            "max_numerator_relative": 0.0,
            "max_denominator_relative": 0.0,
        }
    value_difference = np.abs(candidate.value - reference.value)
    return {
        "occupied_bins": int(np.count_nonzero(occupied)),
        "max_absolute": float(np.max(value_difference[occupied])),
        "max_relative": float(
            np.max(relative_difference(reference.value, candidate.value)[occupied])
        ),
        "max_numerator_relative": float(
            np.max(
                relative_difference(reference.numerator, candidate.numerator)[occupied]
            )
        ),
        "max_denominator_relative": float(
            np.max(
                relative_difference(reference.denominator, candidate.denominator)[occupied]
            )
        ),
    }


def write_2pcf_csv(
    path: Path, config: ComparisonConfig, projected: Histogram,
    scan: Histogram, tree: Histogram
) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            [
                "bin", "rp_center", "xi_3d_projected", "xi_1d_scan",
                "xi_1d_tree", "relative_scan_vs_3d", "relative_tree_vs_3d",
                "numerator_3d_projected", "denominator_3d_projected",
                "numerator_1d_scan", "denominator_1d_scan",
                "numerator_1d_tree", "denominator_1d_tree",
            ]
        )
        scan_relative = relative_difference(projected.value, scan.value)
        tree_relative = relative_difference(projected.value, tree.value)
        for index in range(config.rp_bins):
            writer.writerow(
                [
                    index, (index + 0.5) * config.rp_max / config.rp_bins,
                    projected.value[index], scan.value[index], tree.value[index],
                    scan_relative[index], tree_relative[index],
                    projected.numerator[index], projected.denominator[index],
                    scan.numerator[index], scan.denominator[index],
                    tree.numerator[index], tree.denominator[index],
                ]
            )


def write_3pcf_csv(
    path: Path, config: ComparisonConfig, projected: Histogram,
    radial: Histogram
) -> None:
    total_bins = 2 * config.r3_bins
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            [
                "b1", "b2", "lag1_center", "lag2_center",
                "zeta_3d_projected", "zeta_1d", "relative_1d_vs_3d",
                "numerator_3d_projected", "denominator_3d_projected",
                "numerator_1d", "denominator_1d",
            ]
        )
        relative = relative_difference(projected.value, radial.value)
        for b1 in range(total_bins):
            lag1 = -config.r3_max + (b1 + 0.5) * config.r3_max / config.r3_bins
            for b2 in range(total_bins):
                lag2 = (
                    -config.r3_max
                    + (b2 + 0.5) * config.r3_max / config.r3_bins
                )
                writer.writerow(
                    [
                        b1, b2, lag1, lag2, projected.value[b1, b2],
                        radial.value[b1, b2], relative[b1, b2],
                        projected.numerator[b1, b2], projected.denominator[b1, b2],
                        radial.numerator[b1, b2], radial.denominator[b1, b2],
                    ]
                )


def write_timings_csv(path: Path, timings: Dict[str, float]) -> None:
    references = {
        "lya-1d-2pcf-omp": "lya-2pcf-omp",
        "lya-1d-tree-2pcf-omp": "lya-2pcf-omp",
        "lya-1d-3pcf-omp": "lya-3pcf-omp",
        "lya-1d-2pcf-3pcf-omp": "lya-2pcf-3pcf-omp",
    }
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            ["method", "wall_seconds", "matching_3d_method", "speedup_vs_3d"]
        )
        for method, elapsed in timings.items():
            reference = references.get(method, method)
            writer.writerow([method, elapsed, reference, timings[reference] / elapsed])


def make_plot(
    path: Path, config: ComparisonConfig, projected2: Histogram,
    scan2: Histogram, tree2: Histogram, timings: Dict[str, float],
    projected3: Optional[Histogram] = None,
    radial3: Optional[Histogram] = None,
) -> None:
    mpl_config = config.output_dir / ".matplotlib"
    mpl_config.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", os.fspath(mpl_config))
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise RuntimeError("plotting requires matplotlib") from exc

    centers = (np.arange(config.rp_bins) + 0.5) * config.rp_max / config.rp_bins
    figure, axes = plt.subplots(2, 3, figsize=(15, 8), constrained_layout=True)
    axes[0, 0].plot(centers, projected2.value, "o-", label="3D projected")
    axes[0, 0].plot(centers, scan2.value, "s--", label="1D scan")
    axes[0, 0].plot(centers, tree2.value, "^:", label="1D tree")
    axes[0, 0].set(xlabel=r"$r_\parallel$", ylabel=r"$\xi$", title="2PCF")
    axes[0, 0].legend()
    active2 = projected2.denominator != 0.0
    scan_relative = relative_difference(projected2.value, scan2.value)
    tree_relative = relative_difference(projected2.value, tree2.value)
    axes[1, 0].semilogy(
        centers[active2], np.maximum(scan_relative[active2], 1e-17),
        "s--", label="1D scan",
    )
    axes[1, 0].semilogy(
        centers[active2], np.maximum(tree_relative[active2], 1e-17),
        "^:", label="1D tree",
    )
    axes[1, 0].set(
        xlabel=r"$r_\parallel$", ylabel="symmetric relative difference",
        title="2PCF numerical agreement",
    )
    axes[1, 0].legend()

    if projected3 is not None and radial3 is not None:
        extent = (-config.r3_max, config.r3_max, -config.r3_max, config.r3_max)
        active3 = (projected3.denominator != 0.0) | (radial3.denominator != 0.0)
        if np.any(active3):
            amplitude = max(
                float(np.max(np.abs(projected3.value[active3]))),
                float(np.max(np.abs(radial3.value[active3]))),
                1e-14,
            )
        else:
            amplitude = 1.0
        image1 = axes[0, 1].imshow(
            np.ma.masked_where(projected3.denominator == 0.0, projected3.value),
            origin="lower", extent=extent, cmap="coolwarm",
            vmin=-amplitude, vmax=amplitude, aspect="auto",
        )
        axes[0, 1].set(title="3D 3PCF projected", xlabel="lag 2", ylabel="lag 1")
        figure.colorbar(image1, ax=axes[0, 1], label=r"$\zeta$")
        image2 = axes[0, 2].imshow(
            np.ma.masked_where(radial3.denominator == 0.0, radial3.value),
            origin="lower", extent=extent, cmap="coolwarm",
            vmin=-amplitude, vmax=amplitude, aspect="auto",
        )
        axes[0, 2].set(title="Native 1D 3PCF", xlabel="lag 2", ylabel="lag 1")
        figure.colorbar(image2, ax=axes[0, 2], label=r"$\zeta$")
        relative3 = relative_difference(projected3.value, radial3.value)
        relative_image = np.ma.masked_where(
            ~active3, np.log10(np.maximum(relative3, 1e-17))
        )
        image3 = axes[1, 1].imshow(
            relative_image, origin="lower", extent=extent, cmap="viridis",
            vmin=-17.0, vmax=0.0, aspect="auto",
        )
        axes[1, 1].set(
            title="3PCF numerical agreement", xlabel="lag 2", ylabel="lag 1"
        )
        figure.colorbar(image3, ax=axes[1, 1], label="log10 relative difference")
    else:
        for axis in (axes[0, 1], axes[0, 2], axes[1, 1]):
            axis.set_axis_off()

    labels = list(timings)
    values = [timings[label] for label in labels]
    short_labels = [label.replace("lya-", "").replace("-omp", "") for label in labels]
    axes[1, 2].barh(short_labels, values, color="#4472a8")
    axes[1, 2].invert_yaxis()
    axes[1, 2].set(xlabel="wall time [s]", title=f"Runtime ({config.threads} threads)")
    figure.suptitle("cTreeBalls Ly-alpha: native 1D versus projected 3D")
    figure.savefig(path, dpi=170)
    plt.close(figure)


def run_comparison(config: ComparisonConfig) -> dict:
    config.validate()
    config.output_dir.mkdir(parents=True, exist_ok=True)
    if config.catalog_path is None:
        catalog = generate_collinear_catalog(config)
        catalog_path = config.output_dir / "collinear_catalog.txt"
        write_catalog(catalog_path, catalog)
    else:
        catalog = load_catalog(config.catalog_path)
        catalog_path = config.catalog_path
    validate_collinear_catalog(catalog)
    if config.include_3pcf:
        validate_3pcf_projection_bins(catalog, config)

    methods = ["lya-2pcf-omp", "lya-1d-2pcf-omp", "lya-1d-tree-2pcf-omp"]
    if config.include_3pcf:
        methods.extend(("lya-3pcf-omp", "lya-1d-3pcf-omp"))
    if config.include_combined and config.include_3pcf:
        methods.extend(("lya-1d-2pcf-3pcf-omp", "lya-2pcf-3pcf-omp"))

    roots: Dict[str, Path] = {}
    timings: Dict[str, float] = {}
    for method in methods:
        print(f"Running {method} ...", flush=True)
        roots[method], timings[method] = run_method(config, catalog_path, method)

    projected2 = read_projected_3d_2pcf(
        roots["lya-2pcf-omp"] / "histXi2pcf_lya.txt",
        config.rp_bins,
        config.rt_bins,
    )
    scan2 = read_1d_2pcf(
        roots["lya-1d-2pcf-omp"] / "histXi2pcf_lya1d.txt", config.rp_bins
    )
    tree2 = read_1d_2pcf(
        roots["lya-1d-tree-2pcf-omp"] / "histXi2pcf_lya1d.txt",
        config.rp_bins,
    )
    numerical = {
        "2pcf_scan_vs_projected_3d": comparison_metrics(projected2, scan2),
        "2pcf_tree_vs_projected_3d": comparison_metrics(projected2, tree2),
    }

    projected3 = None
    radial3 = None
    if config.include_3pcf:
        projected3 = read_projected_3d_3pcf(
            roots["lya-3pcf-omp"] / "histZetaM_lya5d.txt",
            config.r3_bins,
            config.theta_bins,
            config.mu_bins,
        )
        radial3 = read_1d_3pcf(
            roots["lya-1d-3pcf-omp"] / "histZetaM_lya1d.txt",
            config.r3_bins,
        )
        numerical["3pcf_1d_vs_projected_3d"] = comparison_metrics(
            projected3, radial3
        )

    if config.include_combined and config.include_3pcf:
        combined1d2 = read_1d_2pcf(
            roots["lya-1d-2pcf-3pcf-omp"] / "histXi2pcf_lya1d.txt",
            config.rp_bins,
        )
        numerical["2pcf_1d_combined_vs_dedicated"] = comparison_metrics(
            scan2, combined1d2
        )
        combined1d3 = read_1d_3pcf(
            roots["lya-1d-2pcf-3pcf-omp"] / "histZetaM_lya1d.txt",
            config.r3_bins,
        )
        combined3d2 = read_projected_3d_2pcf(
            roots["lya-2pcf-3pcf-omp"] / "histXi2pcf_lya.txt",
            config.rp_bins,
            config.rt_bins,
        )
        combined3d3 = read_projected_3d_3pcf(
            roots["lya-2pcf-3pcf-omp"] / "histZetaM_lya5d.txt",
            config.r3_bins,
            config.theta_bins,
            config.mu_bins,
        )
        numerical["3pcf_1d_combined_vs_dedicated"] = comparison_metrics(
            radial3, combined1d3
        )
        numerical["2pcf_3d_combined_vs_dedicated"] = comparison_metrics(
            projected2, combined3d2
        )
        numerical["3pcf_3d_combined_vs_dedicated"] = comparison_metrics(
            projected3, combined3d3
        )

    write_2pcf_csv(
        config.output_dir / "comparison_2pcf.csv",
        config,
        projected2,
        scan2,
        tree2,
    )
    if projected3 is not None and radial3 is not None:
        write_3pcf_csv(
            config.output_dir / "comparison_3pcf.csv",
            config,
            projected3,
            radial3,
        )
    write_timings_csv(config.output_dir / "timings.csv", timings)
    plot_path = config.output_dir / "lya_1d_vs_3d.png"
    make_plot(
        plot_path, config, projected2, scan2, tree2, timings,
        projected3, radial3,
    )

    summary = {
        "config": {
            **asdict(config),
            "output_dir": os.fspath(config.output_dir),
            "cballs_executable": os.fspath(config.cballs_executable),
            "catalog_path": (
                os.fspath(config.catalog_path) if config.catalog_path else None
            ),
        },
        "catalog": os.fspath(catalog_path),
        "timings_seconds": timings,
        "numerical_comparisons": numerical,
        "plot": os.fspath(plot_path),
    }
    (config.output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    return summary


def parse_arguments() -> ComparisonConfig:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir", type=Path,
        default=PROJECT_ROOT / "examples" / "output_lya_1d_vs_3d",
    )
    parser.add_argument("--cballs", type=Path, default=PROJECT_ROOT / "cballs")
    parser.add_argument(
        "--catalog", type=Path,
        help="optional collinear x y z delta weight forest_id catalog",
    )
    parser.add_argument("--n-pixels", type=int, default=160)
    parser.add_argument("--pixels-per-forest", type=int, default=4)
    parser.add_argument("--radial-min", type=float, default=50.0)
    parser.add_argument("--radial-max", type=float, default=80.0)
    parser.add_argument("--threads", type=int, default=min(4, os.cpu_count() or 1))
    parser.add_argument("--seed", type=int, default=20260830)
    parser.add_argument("--rp-max", type=float, default=12.0)
    parser.add_argument("--rp-bins", type=int, default=12)
    parser.add_argument("--rt-max", type=float, default=2.0)
    parser.add_argument("--rt-bins", type=int, default=4)
    parser.add_argument("--r3-max", type=float, default=8.0)
    parser.add_argument("--r3-bins", type=int, default=8)
    parser.add_argument("--theta-bins", type=int, default=4)
    parser.add_argument("--mu-bins", type=int, default=4)
    parser.add_argument("--skip-3pcf", action="store_true")
    parser.add_argument("--skip-combined", action="store_true")
    arguments = parser.parse_args()
    return ComparisonConfig(
        output_dir=arguments.output_dir,
        cballs_executable=arguments.cballs,
        catalog_path=arguments.catalog,
        n_pixels=arguments.n_pixels,
        pixels_per_forest=arguments.pixels_per_forest,
        radial_min=arguments.radial_min,
        radial_max=arguments.radial_max,
        threads=arguments.threads,
        seed=arguments.seed,
        rp_max=arguments.rp_max,
        rp_bins=arguments.rp_bins,
        rt_max=arguments.rt_max,
        rt_bins=arguments.rt_bins,
        r3_max=arguments.r3_max,
        r3_bins=arguments.r3_bins,
        theta_bins=arguments.theta_bins,
        mu_bins=arguments.mu_bins,
        include_3pcf=not arguments.skip_3pcf,
        include_combined=not arguments.skip_combined,
    )


def main() -> None:
    try:
        summary = run_comparison(parse_arguments())
    except (OSError, RuntimeError, ValueError) as exc:
        raise SystemExit(f"ERROR: {exc}") from exc
    print("\nNumerical comparisons:")
    for name, metrics in summary["numerical_comparisons"].items():
        print(
            f"  {name}: max relative={metrics['max_relative']:.3e}, "
            f"max absolute={metrics['max_absolute']:.3e}"
        )
    print("\nWall times:")
    for method, elapsed in summary["timings_seconds"].items():
        print(f"  {method}: {elapsed:.6f} s")
    print(f"\nResults: {summary['config']['output_dir']}")
    print(f"Plot: {summary['plot']}")


if __name__ == "__main__":
    main()
