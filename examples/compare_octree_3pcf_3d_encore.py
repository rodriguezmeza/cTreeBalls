#!/usr/bin/env python3
"""Compare the cTreeBalls 3D survey estimator with ENCORE.

The example generates a deterministic non-periodic survey, runs
``octree-3pcf-3d-omp`` on data and random catalogs, runs ENCORE on the
matching ``D-alpha*R`` and ``alpha*R`` catalogs, applies ENCORE's supplied
3PCF coupling matrix, and compares raw and edge-corrected results.

The bundled ENCORE executable may have been compiled for another platform.
Unless ``--encore-executable`` is supplied, this script compiles a small CPU
executable from the ENCORE source tree with matching NBIN and ORDER values.
The ENCORE source tree is never modified.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import asdict, dataclass
import json
import math
import os
from pathlib import Path
import re
import shutil
import shlex
import signal
import subprocess
import sys
import time
from typing import Dict, Iterable, Optional, Tuple

import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_ENCORE_SOURCE = Path(
    os.environ.get(
        "ENCORE_SOURCE",
        "/Users/mar/Documents/Codex/lyman_alpha/Eladio/encore_2026-08-27",
    )
)


@dataclass
class ComparisonConfig:
    output_dir: Path
    search_method: str = "octree-3pcf-3d-omp"
    mpi_ranks: int = 2
    mpiexec: str = "mpiexec"
    mpi_extra_arg: tuple[str, ...] = ()
    timeout: float = 3600.0
    encore_source: Path = DEFAULT_ENCORE_SOURCE
    encore_executable: Optional[Path] = None
    cballs_executable: Path = PROJECT_ROOT / "cballs"
    backend: str = "cython"
    cxx: str = "c++"
    n_data: int = 90
    n_random: int = 300
    nbins: int = 6
    measured_lmax: int = 3
    rmin: float = 0.05
    rmax: float = 0.80
    threads: int = 4
    encore_nside: int = 12
    seed: int = 8675309

    @property
    def reported_lmax(self) -> int:
        return self.measured_lmax - 1

    def validate(self) -> None:
        if self.search_method not in {"octree-3pcf-3d-omp", "octree-3pcf-3d-mpi"}:
            raise ValueError("select the OMP or MPI physical-3D engine")
        if self.search_method.endswith("-mpi") and self.backend != "cli":
            raise ValueError("MPI examples require backend='cli'")
        if self.mpi_ranks < 1 or not math.isfinite(self.timeout) or self.timeout <= 0:
            raise ValueError("mpi-ranks and timeout must be positive")
        self.output_dir = Path(self.output_dir).expanduser().resolve()
        self.encore_source = Path(self.encore_source).expanduser().resolve()
        self.cballs_executable = Path(self.cballs_executable).expanduser().resolve()
        if self.encore_executable is not None:
            self.encore_executable = Path(
                self.encore_executable
            ).expanduser().resolve()
        if self.backend not in {"cython", "cli"}:
            raise ValueError("backend must be 'cython' or 'cli'")
        if self.n_data < 3 or self.n_random < 3:
            raise ValueError("both catalogs need at least three points")
        if self.nbins < 2:
            raise ValueError("ENCORE 3PCF needs at least two radial bins")
        if not 1 <= self.measured_lmax <= 10:
            raise ValueError("measured_lmax must be in 1..10")
        if not math.isfinite(self.rmin) or not math.isfinite(self.rmax):
            raise ValueError("rmin and rmax must be finite")
        if self.rmin < 0.0 or self.rmax <= self.rmin:
            raise ValueError("require 0 <= rmin < rmax")
        if self.threads < 1 or self.encore_nside < 1:
            raise ValueError("threads and encore_nside must be positive")


def _sample_window(rng: np.random.Generator, count: int) -> np.ndarray:
    """Sample an angular/radial selection with a boundary and a small hole."""
    accepted = []
    while sum(chunk.shape[0] for chunk in accepted) < count:
        points = rng.uniform(-0.72, 0.72, size=(max(4 * count, 256), 3))
        radius = np.linalg.norm(points, axis=1)
        angular_cut = points[:, 2] > -0.48 + 0.22 * points[:, 0]
        radial_cut = (radius > 0.12) & (radius < 0.92)
        hole = (
            (points[:, 0] - 0.18) ** 2
            + (points[:, 1] + 0.16) ** 2
            < 0.055**2
        )
        selected = points[angular_cut & radial_cut & ~hole]
        if selected.size:
            accepted.append(selected)
    return np.concatenate(accepted, axis=0)[:count]


def _sample_clustered_data(
    rng: np.random.Generator, count: int
) -> np.ndarray:
    uniform_count = count // 2
    result = [_sample_window(rng, uniform_count)]
    centers = np.asarray(
        [[-0.22, 0.08, 0.10], [0.30, -0.12, 0.24]], dtype=np.float64
    )
    clustered = []
    needed = count - uniform_count
    while sum(chunk.shape[0] for chunk in clustered) < needed:
        center_ids = rng.integers(0, len(centers), size=max(4 * needed, 128))
        points = centers[center_ids] + rng.normal(
            scale=(0.11, 0.13, 0.10), size=(len(center_ids), 3)
        )
        radius = np.linalg.norm(points, axis=1)
        angular_cut = points[:, 2] > -0.48 + 0.22 * points[:, 0]
        radial_cut = (radius > 0.12) & (radius < 0.92)
        hole = (
            (points[:, 0] - 0.18) ** 2
            + (points[:, 1] + 0.16) ** 2
            < 0.055**2
        )
        selected = points[angular_cut & radial_cut & ~hole]
        if selected.size:
            clustered.append(selected)
    result.append(np.concatenate(clustered, axis=0)[:needed])
    data = np.concatenate(result, axis=0)
    rng.shuffle(data)
    return data


def generate_catalogs(
    config: ComparisonConfig,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return data positions/weights and random positions/weights."""
    rng = np.random.default_rng(config.seed)
    data_positions = _sample_clustered_data(rng, config.n_data)
    random_positions = _sample_window(rng, config.n_random)
    data_weights = 0.85 + 0.20 * (
        1.0 + data_positions[:, 2] + 0.25 * data_positions[:, 0]
    )
    random_weights = 0.80 + 0.18 * (
        1.0 + random_positions[:, 2] + 0.25 * random_positions[:, 0]
    )
    return data_positions, data_weights, random_positions, random_weights


def write_catalog(path: Path, positions: np.ndarray, weights: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    values = np.column_stack((positions, weights))
    np.savetxt(path, values, fmt="%.17e", header="x y z weight")


def write_ctree_catalog(
    path: Path, positions: np.ndarray, weights: np.ndarray
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    values = np.column_stack(
        (positions, np.ones(len(positions), dtype=np.float64), weights)
    )
    np.savetxt(path, values, fmt="%.17e", header="x y z value weight")


def _run_checked(
    command: Iterable[str], cwd: Path, log_path: Path, env=None, timeout=3600.0
) -> float:
    start = time.perf_counter()
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w") as stream:
        process = subprocess.Popen(
            [str(token) for token in command], cwd=cwd, env=env, stdout=stream,
            stderr=subprocess.STDOUT, start_new_session=True,
        )
        try:
            code = process.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            os.killpg(process.pid, signal.SIGKILL)
            process.wait()
            raise RuntimeError(f"command exceeded {timeout}s; see {log_path}") from None
    elapsed = time.perf_counter() - start
    if code != 0:
        raise RuntimeError(
            f"command failed with status {code}: "
            f"{' '.join(str(token) for token in command)}\n"
            f"See {log_path}\n{log_path.read_text()[-3000:]}"
        )
    return elapsed


def thread_environment(threads: int) -> Dict[str, str]:
    """Return a controlled OpenMP environment for one benchmark process."""
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


def _ctree_settings(config: ComparisonConfig, root: Path) -> Dict[str, object]:
    return {
        "searchMethod": "octree-3pcf-3d-omp",
        "rangeN": config.rmax,
        "rminHist": config.rmin,
        "sizeHistN": config.nbins,
        "sizeHistPhi": 8,
        "mChebyshev": config.measured_lmax,
        "lengthBox": 2.0,
        "useLogHist": False,
        "usePeriodic": False,
        "numberThreads": config.threads,
        "scanLevel": 2,
        "stepState": 1000000,
        "verbose": 0,
        "verbose_log": 0,
        "rootDir": str(root),
        "iCatalogs": "1,2",
        "options": (
            "survey-estimator-3d,compute-2pcf-3d,compute-3pcf-3d,"
            "survey-keep-top-multipole"
        ),
    }


def run_ctreeballs_cython(
    config: ComparisonConfig,
    data_positions: np.ndarray,
    data_weights: np.ndarray,
    random_positions: np.ndarray,
    random_weights: np.ndarray,
) -> float:
    try:
        from cyballs import cballs
    except ImportError as exc:
        raise RuntimeError(
            "cyballs is not importable; install it or use --backend cli"
        ) from exc

    root = config.output_dir / "ctreeballs"
    root.mkdir(parents=True, exist_ok=True)
    balls = cballs()
    balls.set(_ctree_settings(config, root))
    balls.set_catalog(data_positions, weights=data_weights, catalog=0)
    balls.set_catalog(random_positions, weights=random_weights, catalog=1)
    start = time.perf_counter()
    try:
        balls.Run(level=["MainLoop"])
    finally:
        balls.struct_cleanup()
    return time.perf_counter() - start


def run_ctreeballs_cli(
    config: ComparisonConfig,
    data_positions: np.ndarray,
    data_weights: np.ndarray,
    random_positions: np.ndarray,
    random_weights: np.ndarray,
) -> float:
    if not config.cballs_executable.is_file():
        raise FileNotFoundError(
            f"cballs executable not found: {config.cballs_executable}"
        )
    root = config.output_dir / "ctreeballs"
    inputs = config.output_dir / "catalogs"
    root.mkdir(parents=True, exist_ok=True)
    data_file = inputs / "data_ctreeballs.txt"
    random_file = inputs / "random_ctreeballs.txt"
    write_ctree_catalog(data_file, data_positions, data_weights)
    write_ctree_catalog(random_file, random_positions, random_weights)
    command = [
        str(config.cballs_executable),
        f"search={config.search_method}",
        f"in={data_file},{random_file}",
        "infmt=multi-columns-ascii,multi-columns-ascii",
        "columns=1,2,3,4,5",
        "iCatalogs=1,2",
        f"rootDir={root}",
        f"numberThreads={config.threads}",
        "usePeriodic=false",
        "useLogHist=false",
        f"rangeN={config.rmax:.17g}",
        f"rminHist={config.rmin:.17g}",
        f"sizeHistN={config.nbins}",
        f"mChebyshev={config.measured_lmax}",
        "stepState=1000000",
        "verbose=0",
        "verbose_log=0",
        (
            "options=pos-and-convergence-weight,survey-estimator-3d,"
            "compute-2pcf-3d,compute-3pcf-3d,"
            "survey-keep-top-multipole"
        ),
    ]
    if config.search_method.endswith("-mpi"):
        command = [*shlex.split(config.mpiexec), *config.mpi_extra_arg,
                   "-n", str(config.mpi_ranks), *command]
    return _run_checked(
        command,
        PROJECT_ROOT,
        root / "run.log",
        env=thread_environment(config.threads),
        timeout=config.timeout,
    )


def _replace_define(source: str, name: str, value: int) -> str:
    pattern = rf"(?m)^#define\s+{re.escape(name)}\s+\d+\s*$"
    replacement = f"#define {name} {value}"
    updated, count = re.subn(pattern, replacement, source)
    if count != 1:
        raise RuntimeError(
            f"expected one '#define {name}' in ENCORE source, found {count}"
        )
    return updated


def prepare_encore(config: ComparisonConfig) -> Tuple[Path, Path, float]:
    """Return executable, runtime directory, and native build time."""
    source = config.encore_source
    if not (source / "encore.cpp").is_file():
        raise FileNotFoundError(
            f"ENCORE source not found at {source}; use --encore-source"
        )
    coupling_source = source / "coupling_matrices"
    coupling_file = coupling_source / (
        f"edge_correction_matrix_3pcf_LMAX{config.measured_lmax}.npy"
    )
    if not coupling_file.is_file():
        raise FileNotFoundError(
            f"ENCORE coupling matrix not found: {coupling_file}"
        )

    runtime = config.output_dir / "encore"
    runtime.mkdir(parents=True, exist_ok=True)
    shutil.copytree(
        coupling_source, runtime / "coupling_matrices", dirs_exist_ok=True
    )
    (runtime / "output").mkdir(exist_ok=True)

    if config.encore_executable is not None:
        if not config.encore_executable.is_file():
            raise FileNotFoundError(
                f"ENCORE executable not found: {config.encore_executable}"
            )
        return config.encore_executable, runtime, 0.0

    configured_source = runtime / "encore_configured.cpp"
    source_text = (source / "encore.cpp").read_text(encoding="utf-8")
    source_text = _replace_define(source_text, "NBIN", config.nbins)
    source_text = _replace_define(
        source_text, "ORDER", config.measured_lmax
    )
    source_text = _replace_define(source_text, "MAXTHREAD", 1)
    configured_source.write_text(source_text, encoding="utf-8")
    executable = runtime / "encore_native"
    command = [
        config.cxx,
        "-std=c++11",
        "-O3",
        f"-I{source}",
        str(configured_source),
        "-o",
        str(executable),
    ]
    elapsed = _run_checked(command, runtime, runtime / "build.log")
    return executable, runtime, elapsed


def run_encore(
    config: ComparisonConfig,
    data_positions: np.ndarray,
    data_weights: np.ndarray,
    random_positions: np.ndarray,
    random_weights: np.ndarray,
) -> Dict[str, float]:
    executable, runtime, build_seconds = prepare_encore(config)
    alpha = float(np.sum(data_weights) / np.sum(random_weights))
    dmr_positions = np.concatenate((data_positions, random_positions), axis=0)
    dmr_weights = np.concatenate(
        (data_weights, -alpha * random_weights), axis=0
    )
    normalized_random_weights = alpha * random_weights
    input_dir = config.output_dir / "catalogs"
    dmr_file = input_dir / "data_minus_random_encore.txt"
    random_file = input_dir / "normalized_random_encore.txt"
    write_catalog(dmr_file, dmr_positions, dmr_weights)
    write_catalog(random_file, random_positions, normalized_random_weights)

    common = [
        str(executable),
        "-rmin",
        f"{config.rmin:.17g}",
        "-rmax",
        f"{config.rmax:.17g}",
        "-nside",
        str(config.encore_nside),
    ]
    env = thread_environment(config.threads)
    random_seconds = _run_checked(
        common + ["-in", str(random_file), "-outstr", "comparison.r"],
        runtime,
        runtime / "random.log",
        env=env,
    )
    numerator_seconds = _run_checked(
        common + ["-in", str(dmr_file), "-outstr", "comparison.n00"],
        runtime,
        runtime / "data_minus_random.log",
        env=env,
    )
    return {
        "alpha": alpha,
        "build_seconds": build_seconds,
        "random_seconds": random_seconds,
        "numerator_seconds": numerator_seconds,
    }


def _load_numeric_rows(path: Path) -> np.ndarray:
    rows = [
        [float(token) for token in line.split()]
        for line in path.read_text(encoding="ascii").splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]
    return np.asarray(rows, dtype=np.float64)


def read_ctreeballs(config: ComparisonConfig) -> Dict[str, np.ndarray]:
    root = config.output_dir / "ctreeballs"
    xi = np.atleast_2d(
        _load_numeric_rows(root / "histXi2pcf_3d_survey.txt")
    )
    zeta = np.atleast_2d(
        _load_numeric_rows(root / "histZetaM_3d_survey.txt")
    )
    if xi.shape[1] != 6 or zeta.shape[1] != 14:
        raise RuntimeError("unexpected cTreeBalls survey output format")
    return {"xi": xi, "zeta": zeta}


def read_encore_2pcf(path: Path, nbins: int) -> np.ndarray:
    rows = _load_numeric_rows(path)
    if rows.shape != (2, nbins):
        raise RuntimeError(f"unexpected ENCORE 2PCF format in {path}")
    expected_bins = np.arange(nbins, dtype=np.float64)
    if not np.array_equal(rows[0], expected_bins):
        raise RuntimeError("ENCORE 2PCF radial-bin labels do not match")
    return rows[1]


def read_encore_3pcf(
    path: Path, nbins: int, measured_lmax: int
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    lines = [
        line.split()
        for line in path.read_text(encoding="ascii").splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]
    pair_count = nbins * (nbins - 1) // 2
    if len(lines) != measured_lmax + 3:
        raise RuntimeError(f"unexpected ENCORE 3PCF row count in {path}")
    bin1 = np.asarray(lines[0], dtype=np.int64)
    bin2 = np.asarray(lines[1], dtype=np.int64)
    if len(bin1) != pair_count or len(bin2) != pair_count:
        raise RuntimeError("unexpected ENCORE radial-pair count")
    values = np.empty((measured_lmax + 1, pair_count), dtype=np.float64)
    for expected_ell, tokens in enumerate(lines[2:]):
        if len(tokens) != pair_count + 1 or int(tokens[0]) != expected_ell:
            raise RuntimeError("unexpected ENCORE multipole layout")
        values[expected_ell] = np.asarray(tokens[1:], dtype=np.float64)
    return bin1, bin2, values


def read_and_correct_encore(
    config: ComparisonConfig,
) -> Dict[str, np.ndarray]:
    output = config.output_dir / "encore" / "output"
    numerator_2 = read_encore_2pcf(
        output / "comparison.n00_2pcf.txt", config.nbins
    )
    random_2 = read_encore_2pcf(
        output / "comparison.r_2pcf.txt", config.nbins
    )
    bin1, bin2, numerator_3 = read_encore_3pcf(
        output / "comparison.n00_3pcf.txt",
        config.nbins,
        config.measured_lmax,
    )
    random_bin1, random_bin2, random_3 = read_encore_3pcf(
        output / "comparison.r_3pcf.txt",
        config.nbins,
        config.measured_lmax,
    )
    if not np.array_equal(bin1, random_bin1) or not np.array_equal(
        bin2, random_bin2
    ):
        raise RuntimeError("ENCORE numerator/random radial bins differ")

    xi = np.divide(
        numerator_2,
        random_2,
        out=np.zeros_like(numerator_2),
        where=random_2 != 0.0,
    )
    coupling_path = config.encore_source / "coupling_matrices" / (
        f"edge_correction_matrix_3pcf_LMAX{config.measured_lmax}.npy"
    )
    coupling = np.load(coupling_path)
    mode_count = config.measured_lmax + 1
    if coupling.shape != (mode_count, mode_count, mode_count):
        raise RuntimeError(
            f"unexpected ENCORE coupling shape {coupling.shape}"
        )
    corrected = np.zeros_like(numerator_3)
    valid = np.zeros(numerator_3.shape[1], dtype=bool)
    condition = np.full(numerator_3.shape[1], np.inf)
    for pair in range(numerator_3.shape[1]):
        r0 = random_3[0, pair]
        if not np.isfinite(r0) or r0 == 0.0:
            continue
        window = random_3[:, pair] / r0
        matrix = np.einsum("ijk,k->ij", coupling, window)
        try:
            corrected[:, pair] = np.linalg.solve(
                matrix, numerator_3[:, pair] / r0
            )
        except np.linalg.LinAlgError:
            continue
        condition[pair] = np.linalg.cond(matrix)
        valid[pair] = np.isfinite(corrected[:, pair]).all()
        if not valid[pair]:
            corrected[:, pair] = 0.0
    return {
        "xi": xi,
        "numerator_2": numerator_2,
        "random_2": random_2,
        "bin1": bin1,
        "bin2": bin2,
        "numerator_3": numerator_3,
        "random_3": random_3,
        "corrected_3": corrected,
        "valid_3": valid,
        "condition_3": condition,
    }


def _relative_difference(
    left: np.ndarray, right: np.ndarray
) -> Tuple[np.ndarray, float]:
    reference_scale = float(
        max(np.max(np.abs(left)), np.max(np.abs(right)), 1.0)
    )
    floor = 1.0e-12 * reference_scale
    denominator = np.maximum(np.abs(right), floor)
    return (left - right) / denominator, floor


def align_results(
    config: ComparisonConfig,
    ctree: Dict[str, np.ndarray],
    encore: Dict[str, np.ndarray],
) -> Dict[str, np.ndarray]:
    centers = config.rmin + (
        np.arange(config.nbins, dtype=np.float64) + 0.5
    ) * (config.rmax - config.rmin) / config.nbins
    ctree_xi = ctree["xi"][:, 2]
    if len(ctree_xi) != config.nbins:
        raise RuntimeError("cTreeBalls 2PCF bin count differs from ENCORE")
    xi_relative, xi_floor = _relative_difference(ctree_xi, encore["xi"])
    xi_numerator_relative, _ = _relative_difference(
        ctree["xi"][:, 3], encore["numerator_2"]
    )
    xi_random_relative, _ = _relative_difference(
        ctree["xi"][:, 4], encore["random_2"]
    )

    rows_by_key = {
        (int(row[0]), int(row[1]) - 1, int(row[2]) - 1): row
        for row in ctree["zeta"]
    }
    pair_count = len(encore["bin1"])
    ctree_corrected = np.empty(
        (config.measured_lmax + 1, pair_count), dtype=np.float64
    )
    ctree_numerator = np.empty_like(ctree_corrected)
    ctree_random = np.empty_like(ctree_corrected)
    ctree_valid = np.empty_like(ctree_corrected, dtype=bool)
    for ell in range(config.measured_lmax + 1):
        for pair, (first, second) in enumerate(
            zip(encore["bin1"], encore["bin2"])
        ):
            key = (ell, int(first), int(second))
            if key not in rows_by_key:
                raise RuntimeError(f"missing cTreeBalls row {key}")
            row = rows_by_key[key]
            ctree_corrected[ell, pair] = row[6]
            ctree_numerator[ell, pair] = row[9]
            ctree_random[ell, pair] = row[10]
            ctree_valid[ell, pair] = bool(int(row[13]))

    zeta_relative, zeta_floor = _relative_difference(
        ctree_corrected[: config.reported_lmax + 1],
        encore["corrected_3"][: config.reported_lmax + 1],
    )
    numerator_relative, numerator_floor = _relative_difference(
        ctree_numerator, encore["numerator_3"]
    )
    random_relative, random_floor = _relative_difference(
        ctree_random, encore["random_3"]
    )
    return {
        "centers": centers,
        "ctree_xi": ctree_xi,
        "encore_xi": encore["xi"],
        "xi_relative": xi_relative,
        "xi_numerator_relative": xi_numerator_relative,
        "xi_random_relative": xi_random_relative,
        "ctree_xi_numerator": ctree["xi"][:, 3],
        "ctree_xi_random": ctree["xi"][:, 4],
        "encore_xi_numerator": encore["numerator_2"],
        "encore_xi_random": encore["random_2"],
        "bin1": encore["bin1"],
        "bin2": encore["bin2"],
        "ctree_zeta": ctree_corrected,
        "encore_zeta": encore["corrected_3"],
        "zeta_relative": zeta_relative,
        "ctree_numerator_3": ctree_numerator,
        "encore_numerator_3": encore["numerator_3"],
        "numerator_relative": numerator_relative,
        "ctree_random_3": ctree_random,
        "encore_random_3": encore["random_3"],
        "random_relative": random_relative,
        "ctree_valid_3": ctree_valid,
        "encore_valid_3": encore["valid_3"],
        "encore_condition_3": encore["condition_3"],
        "xi_relative_floor": np.asarray(xi_floor),
        "zeta_relative_floor": np.asarray(zeta_floor),
        "numerator_relative_floor": np.asarray(numerator_floor),
        "random_relative_floor": np.asarray(random_floor),
    }


def write_tables(config: ComparisonConfig, aligned: Dict[str, np.ndarray]) -> None:
    with (config.output_dir / "comparison_2pcf.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.writer(stream)
        writer.writerow(
            [
                "bin",
                "r_center",
                "ctreeballs_xi",
                "encore_xi",
                "absolute_difference",
                "relative_difference",
                "ctreeballs_numerator",
                "encore_numerator",
                "ctreeballs_random",
                "encore_random",
            ]
        )
        for radial_bin in range(config.nbins):
            left = aligned["ctree_xi"][radial_bin]
            right = aligned["encore_xi"][radial_bin]
            writer.writerow(
                [
                    radial_bin,
                    aligned["centers"][radial_bin],
                    left,
                    right,
                    left - right,
                    aligned["xi_relative"][radial_bin],
                    aligned["ctree_xi_numerator"][radial_bin],
                    aligned["encore_xi_numerator"][radial_bin],
                    aligned["ctree_xi_random"][radial_bin],
                    aligned["encore_xi_random"][radial_bin],
                ]
            )

    with (config.output_dir / "comparison_3pcf.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.writer(stream)
        writer.writerow(
            [
                "ell",
                "bin1",
                "bin2",
                "r1_center",
                "r2_center",
                "ctreeballs_zeta_encore_basis",
                "encore_zeta",
                "absolute_difference",
                "relative_difference",
                "ctreeballs_numerator",
                "encore_numerator",
                "ctreeballs_random",
                "encore_random",
                "ctreeballs_valid",
                "encore_valid",
                "encore_condition",
            ]
        )
        for ell in range(config.reported_lmax + 1):
            for pair, (first, second) in enumerate(
                zip(aligned["bin1"], aligned["bin2"])
            ):
                left = aligned["ctree_zeta"][ell, pair]
                right = aligned["encore_zeta"][ell, pair]
                writer.writerow(
                    [
                        ell,
                        int(first),
                        int(second),
                        aligned["centers"][int(first)],
                        aligned["centers"][int(second)],
                        left,
                        right,
                        left - right,
                        aligned["zeta_relative"][ell, pair],
                        aligned["ctree_numerator_3"][ell, pair],
                        aligned["encore_numerator_3"][ell, pair],
                        aligned["ctree_random_3"][ell, pair],
                        aligned["encore_random_3"][ell, pair],
                        int(aligned["ctree_valid_3"][ell, pair]),
                        int(aligned["encore_valid_3"][pair]),
                        aligned["encore_condition_3"][pair],
                    ]
                )


def make_plots(config: ComparisonConfig, aligned: Dict[str, np.ndarray]) -> None:
    matplotlib_cache = config.output_dir / ".matplotlib"
    matplotlib_cache.mkdir(exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(matplotlib_cache))
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise RuntimeError("matplotlib is required to create the plots") from exc

    figure, axes = plt.subplots(
        2, 1, figsize=(8.2, 6.4), sharex=True,
        gridspec_kw={"height_ratios": (2.2, 1.0)},
    )
    axes[0].plot(
        aligned["centers"], aligned["ctree_xi"], "o-", label="cTreeBalls"
    )
    axes[0].plot(
        aligned["centers"], aligned["encore_xi"], "s--", label="ENCORE"
    )
    axes[0].set_ylabel(r"$\xi(r)$")
    axes[0].grid(alpha=0.25)
    axes[0].legend()
    axes[1].axhline(0.0, color="black", linewidth=0.8)
    axes[1].plot(aligned["centers"], aligned["xi_relative"], "o-")
    axes[1].set_xlabel("Radial separation")
    axes[1].set_ylabel("Relative diff.")
    axes[1].grid(alpha=0.25)
    figure.suptitle("Survey 2PCF: cTreeBalls versus ENCORE")
    figure.tight_layout()
    figure.savefig(config.output_dir / "comparison_2pcf.png", dpi=170)
    plt.close(figure)

    pair_index = np.arange(len(aligned["bin1"]))
    figure, axes = plt.subplots(
        2, 1, figsize=(10.0, 7.0), sharex=True,
        gridspec_kw={"height_ratios": (2.4, 1.1)},
    )
    for ell in range(config.reported_lmax + 1):
        axes[0].plot(
            pair_index,
            aligned["ctree_zeta"][ell],
            marker="o",
            linewidth=1.2,
            label=fr"cTreeBalls $\ell={ell}$",
        )
        axes[0].plot(
            pair_index,
            aligned["encore_zeta"][ell],
            linestyle="none",
            marker="x",
            markersize=6,
            label=fr"ENCORE $\ell={ell}$",
        )
        axes[1].plot(
            pair_index,
            aligned["zeta_relative"][ell],
            marker="o",
            linewidth=1.1,
            label=fr"$\ell={ell}$",
        )
    axes[0].set_ylabel(r"Edge-corrected $\zeta_\ell$")
    axes[0].grid(alpha=0.25)
    axes[0].legend(ncol=2, fontsize=8)
    axes[1].axhline(0.0, color="black", linewidth=0.8)
    axes[1].set_xlabel("Radial-bin pair index")
    axes[1].set_ylabel("Relative diff.")
    axes[1].grid(alpha=0.25)
    figure.suptitle("Survey 3PCF multipoles: cTreeBalls versus ENCORE")
    figure.tight_layout()
    figure.savefig(config.output_dir / "comparison_3pcf.png", dpi=170)
    plt.close(figure)


def summarize(
    config: ComparisonConfig,
    aligned: Dict[str, np.ndarray],
    ctree_seconds: float,
    encore_timings: Dict[str, float],
) -> Dict[str, object]:
    science_slice = slice(0, config.reported_lmax + 1)
    summary = {
        "configuration": {
            **asdict(config),
            "output_dir": str(config.output_dir),
            "encore_source": str(config.encore_source),
            "encore_executable": (
                str(config.encore_executable)
                if config.encore_executable is not None
                else None
            ),
            "cballs_executable": str(config.cballs_executable),
        },
        "normalization": {
            "alpha": encore_timings["alpha"],
            "reported_lmax": config.reported_lmax,
            "measured_lmax": config.measured_lmax,
            "relative_difference_definition": (
                "(cTreeBalls-ENCORE)/max(abs(ENCORE), floor)"
            ),
            "xi_floor": float(aligned["xi_relative_floor"]),
            "zeta_floor": float(aligned["zeta_relative_floor"]),
            "encore_text_output_precision": (
                "ENCORE writes raw counts with %le (six digits after the "
                "decimal), limiting file-based comparisons to about 1e-6"
            ),
        },
        "differences": {
            "xi_max_absolute": float(
                np.max(np.abs(aligned["ctree_xi"] - aligned["encore_xi"]))
            ),
            "xi_max_relative": float(np.max(np.abs(aligned["xi_relative"]))),
            "raw_xi_numerator_max_relative": float(
                np.max(np.abs(aligned["xi_numerator_relative"]))
            ),
            "raw_xi_random_max_relative": float(
                np.max(np.abs(aligned["xi_random_relative"]))
            ),
            "zeta_max_absolute": float(
                np.max(
                    np.abs(
                        aligned["ctree_zeta"][science_slice]
                        - aligned["encore_zeta"][science_slice]
                    )
                )
            ),
            "zeta_max_relative": float(
                np.max(np.abs(aligned["zeta_relative"]))
            ),
            "raw_numerator_max_relative": float(
                np.max(np.abs(aligned["numerator_relative"]))
            ),
            "raw_random_max_relative": float(
                np.max(np.abs(aligned["random_relative"]))
            ),
        },
        "timings_seconds": {
            "ctreeballs": ctree_seconds,
            **{
                key: value
                for key, value in encore_timings.items()
                if key.endswith("seconds")
            },
        },
    }
    summary["agreement_within_encore_text_precision"] = bool(
        summary["differences"]["xi_max_relative"] < 2.0e-5
        and summary["differences"]["zeta_max_relative"] < 2.0e-5
        and summary["differences"]["raw_numerator_max_relative"] < 2.0e-5
        and summary["differences"]["raw_random_max_relative"] < 2.0e-5
    )
    (config.output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    return summary


def run_comparison(config: ComparisonConfig) -> Dict[str, object]:
    config.validate()
    config.output_dir.mkdir(parents=True, exist_ok=True)
    data_positions, data_weights, random_positions, random_weights = (
        generate_catalogs(config)
    )
    catalogs = config.output_dir / "catalogs"
    write_catalog(catalogs / "data.txt", data_positions, data_weights)
    write_catalog(catalogs / "random.txt", random_positions, random_weights)

    if config.backend == "cython":
        ctree_seconds = run_ctreeballs_cython(
            config,
            data_positions,
            data_weights,
            random_positions,
            random_weights,
        )
    else:
        ctree_seconds = run_ctreeballs_cli(
            config,
            data_positions,
            data_weights,
            random_positions,
            random_weights,
        )
    encore_timings = run_encore(
        config,
        data_positions,
        data_weights,
        random_positions,
        random_weights,
    )
    ctree = read_ctreeballs(config)
    correction_start = time.perf_counter()
    encore = read_and_correct_encore(config)
    encore_timings["correction_seconds"] = (
        time.perf_counter() - correction_start
    )
    aligned = align_results(config, ctree, encore)
    write_tables(config, aligned)
    make_plots(config, aligned)
    summary = summarize(config, aligned, ctree_seconds, encore_timings)
    return {
        "config": config,
        "catalogs": {
            "data_positions": data_positions,
            "data_weights": data_weights,
            "random_positions": random_positions,
            "random_weights": random_weights,
        },
        "ctreeballs": ctree,
        "encore": encore,
        "aligned": aligned,
        "summary": summary,
    }


def print_summary(summary: Dict[str, object]) -> None:
    differences = summary["differences"]
    timings = summary["timings_seconds"]
    print("cTreeBalls versus ENCORE survey comparison")
    print(f"  xi max absolute difference: {differences['xi_max_absolute']:.6e}")
    print(f"  xi max relative difference: {differences['xi_max_relative']:.6e}")
    print(
        "  zeta max absolute difference: "
        f"{differences['zeta_max_absolute']:.6e}"
    )
    print(
        "  zeta max relative difference: "
        f"{differences['zeta_max_relative']:.6e}"
    )
    print(
        "  raw numerator max relative difference: "
        f"{differences['raw_numerator_max_relative']:.6e}"
    )
    print(
        "  raw random max relative difference: "
        f"{differences['raw_random_max_relative']:.6e}"
    )
    print(
        "  agreement within ENCORE text precision: "
        f"{summary['agreement_within_encore_text_precision']}"
    )
    print(f"  cTreeBalls run time: {timings['ctreeballs']:.3f} s")
    encore_count = timings["random_seconds"] + timings["numerator_seconds"]
    encore_total = encore_count + timings.get("correction_seconds", 0.0)
    print(f"  ENCORE count time: {encore_count:.3f} s")
    print(f"  ENCORE edge-correction time: {timings.get('correction_seconds', 0.0):.3f} s")
    print(f"  ENCORE complete estimator time: {encore_total:.3f} s")


def parse_args(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("octree_3pcf_3d_encore_output"),
    )
    parser.add_argument(
        "--encore-source", type=Path, default=DEFAULT_ENCORE_SOURCE
    )
    parser.add_argument("--encore-executable", type=Path)
    parser.add_argument(
        "--cballs", type=Path, default=PROJECT_ROOT / "cballs"
    )
    parser.add_argument("--backend", choices=("cython", "cli"), default="cython")
    parser.add_argument("--search-method", choices=("octree-3pcf-3d-omp", "octree-3pcf-3d-mpi"),
                        default="octree-3pcf-3d-omp")
    parser.add_argument("--mpi-ranks", type=int, default=2)
    parser.add_argument("--mpiexec", default="mpiexec")
    parser.add_argument("--mpi-extra-arg", action="append", default=[])
    parser.add_argument("--timeout", type=float, default=3600.0)
    parser.add_argument("--cxx", default="c++")
    parser.add_argument("--n-data", type=int, default=90)
    parser.add_argument("--n-random", type=int, default=300)
    parser.add_argument("--nbins", type=int, default=6)
    parser.add_argument("--lmax", type=int, default=3)
    parser.add_argument("--rmin", type=float, default=0.05)
    parser.add_argument("--rmax", type=float, default=0.80)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--encore-nside", type=int, default=12)
    parser.add_argument("--seed", type=int, default=8675309)
    return parser.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    config = ComparisonConfig(
        output_dir=args.output,
        encore_source=args.encore_source,
        encore_executable=args.encore_executable,
        cballs_executable=args.cballs,
        backend=args.backend,
        search_method=args.search_method, mpi_ranks=args.mpi_ranks,
        mpiexec=args.mpiexec, mpi_extra_arg=args.mpi_extra_arg, timeout=args.timeout,
        cxx=args.cxx,
        n_data=args.n_data,
        n_random=args.n_random,
        nbins=args.nbins,
        measured_lmax=args.lmax,
        rmin=args.rmin,
        rmax=args.rmax,
        threads=args.threads,
        encore_nside=args.encore_nside,
        seed=args.seed,
    )
    result = run_comparison(config)
    print_summary(result["summary"])
    print(f"Outputs: {config.output_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
