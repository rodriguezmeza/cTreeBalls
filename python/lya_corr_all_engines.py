#!/usr/bin/env python3
"""Run the native cTreeBalls Ly-alpha forest engines on one retained catalog.

DESI/PICCA FITS, six-column lya-ascii, NPZ and synthetic input are supported.
MPI rank zero reads once, broadcasts arrays once, and all ranks register a
forest catalog with cyballs. Radial and anisotropic 3D estimators are distinct.
See python/README_lya_corr_all_engines.md for commands and scientific limitations.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass, field
import glob
import json
import math
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import time

import numpy as np

from kappa_corr_all_engines import (
    broadcast_array, discover_cython_methods, get_mpi_comm,
    mpi_environment_size,
)

PROJECT_ROOT = Path(__file__).resolve().parents[1]
MPI_CHILD = "CTREEBALLS_LYA_MPI_CHILD"
DESI_BASE = ("https://data.desi.lbl.gov/public/dr1/vac/dr1/lya-deltas/"
             "v1.0/delta-lya-0-0/Delta/")
DESI_EXAMPLE = DESI_BASE + "delta-1019.fits.gz"
LYA_WAVELENGTH = 1215.67


@dataclass(frozen=True)
class EngineSpec:
    radial: bool
    orders: tuple[int, ...]
    mpi: bool


LYA_ENGINES = {
    name + "-" + parallel: EngineSpec(radial, orders, parallel == "mpi")
    for parallel in ("omp", "mpi")
    for name, radial, orders in (
        ("lya-2pcf", False, (2,)),
        ("lya-3pcf", False, (3,)),
        ("lya-2pcf-3pcf", False, (2, 3)),
        ("lya-1d-2pcf", True, (2,)),
        ("lya-1d-3pcf", True, (3,)),
        ("lya-1d-2pcf-3pcf", True, (2, 3)),
        ("lya-1d-tree-2pcf", True, (2,)),
    )
}
# These scalar kernels do not enforce three distinct forests. They must not
# silently substitute for the forest estimators, even with exclude-same-los.
EXCLUDED = (
    "Angular convergence/shear, point-count and octree-3pcf-3d engines are not "
    "forest estimators: their geometry or same-quasar triplet exclusions differ."
)
PRODUCTS = {
    "3d_2pcf": ("histXi2pcf_lya.txt", 2, 7),
    "3d_3pcf": ("histZetaM_lya5d.txt", 5, 13),
    "1d_2pcf": ("histXi2pcf_lya1d.txt", 1, 5),
    "1d_3pcf": ("histZetaM_lya1d.txt", 2, 7),
}


@dataclass
class ForestCatalog:
    positions: np.ndarray
    delta: np.ndarray
    weights: np.ndarray
    forest_ids: np.ndarray
    metadata: dict = field(default_factory=dict)

    def normalized(self):
        pos = np.ascontiguousarray(self.positions, dtype=np.float64)
        delta = np.ascontiguousarray(self.delta, dtype=np.float64)
        weights = np.ascontiguousarray(self.weights, dtype=np.float64)
        ids = np.asarray(self.forest_ids)
        if pos.ndim != 2 or pos.shape[1] != 3 or len(pos) < 3:
            raise ValueError("positions must have shape (N, 3), N >= 3")
        if any(v.ndim != 1 or len(v) != len(pos) for v in (delta, weights, ids)):
            raise ValueError("delta, weights and forest_ids must have shape (N,)")
        if ids.dtype.kind not in "iu":
            raise ValueError("forest_ids must be integers, not rounded floating-point IDs")
        if ids.dtype.kind == "u" and np.any(ids > np.iinfo(np.int64).max):
            raise ValueError("forest_ids exceed signed int64")
        if not all(np.all(np.isfinite(v)) for v in (pos, delta, weights)):
            raise ValueError("catalog contains non-finite values")
        if np.any(weights < 0) or not np.any(weights > 0):
            raise ValueError("weights must be non-negative with at least one positive value")
        distances = np.hypot(np.hypot(pos[:, 0], pos[:, 1]), pos[:, 2])
        if not np.all(np.isfinite(distances) & (distances > 0)):
            raise ValueError("positions need positive finite observer distances")
        return ForestCatalog(pos, delta, weights, np.ascontiguousarray(ids, dtype=np.int64),
                             dict(self.metadata))

    @property
    def nbody(self):
        return len(self.positions)


@dataclass
class RunConfig:
    engines: tuple[str, ...]
    output_dir: Path
    threads: int = 1
    rp_max: float = 200.0
    rt_max: float = 200.0
    rp_bins: int = 50
    rt_bins: int = 50
    r3_max: float = 20.0
    r3_bins: int = 4
    theta_bins: int = 4
    mu_bins: int = 4
    max_hist_mib: float = 1024.0
    plots: bool = True

    def validate(self):
        if not self.engines or any(e not in LYA_ENGINES for e in self.engines):
            raise ValueError("select at least one registered Ly-alpha engine")
        for name in ("threads", "rp_bins", "rt_bins", "r3_bins", "theta_bins", "mu_bins"):
            value = getattr(self, name)
            if not isinstance(value, int) or not 1 <= value <= 100000:
                raise ValueError(f"{name} must be an integer in [1, 100000]")
        for name in ("rp_max", "rt_max", "r3_max", "max_hist_mib"):
            if not math.isfinite(getattr(self, name)) or getattr(self, name) <= 0:
                raise ValueError(f"{name} must be positive and finite")
        if not math.isfinite(math.hypot(self.rp_max, self.rt_max)):
            raise ValueError("separation limits overflow")
        for engine in self.engines:
            spec = LYA_ENGINES[engine]
            bins2 = self.rp_bins * (1 if spec.radial else self.rt_bins)
            bins3 = ((2 * self.r3_bins)**2 if spec.radial else
                     self.r3_bins**2 * self.theta_bins**2 * self.mu_bins)
            bins = (bins2 if 2 in spec.orders else 0) + (bins3 if 3 in spec.orders else 0)
            # Raw sums/counts plus reduction scratch; trees/catalogs are extra.
            mib = bins * 64 * (self.threads + 2) / 2**20
            if mib > self.max_hist_mib:
                raise ValueError(f"{engine}: histogram estimate {mib:.1f} MiB/rank "
                                 f"exceeds --max-hist-mib={self.max_hist_mib}")
        return self


def resolve_engines(tokens, available, statistics="2pcf"):
    requested = [s.strip() for token in tokens for s in token.split(",") if s.strip()]
    selected = []
    for token in requested:
        if token in ("all", "all-omp", "all-mpi", "all-1d", "all-3d"):
            candidates = [name for name, spec in LYA_ENGINES.items()
                          if name in available
                          and (token != "all-omp" or not spec.mpi)
                          and (token != "all-mpi" or spec.mpi)
                          and (token != "all-1d" or spec.radial)
                          and (token != "all-3d" or not spec.radial)]
            if statistics != "both":
                candidates = [n for n in candidates
                              if LYA_ENGINES[n].orders == (int(statistics[0]),)]
        else:
            if token not in LYA_ENGINES:
                raise ValueError(f"{token}: not a native forest engine. {EXCLUDED}")
            if token not in available:
                raise ValueError(f"{token} is not compiled in the imported cyballs")
            if statistics != "both" and LYA_ENGINES[token].orders != (int(statistics[0]),):
                raise ValueError(f"{token} needs --statistics=both or its matching order")
            candidates = [token]
        selected.extend(n for n in candidates if n not in selected)
    if not selected:
        raise ValueError("no requested engines are compiled; enable LYAFORESTOMPON/"
                         "LYAFORESTMPION and rebuild cyballs")
    return tuple(selected)


def expand_inputs(patterns):
    paths = []
    for pattern in patterns:
        matches = sorted(glob.glob(os.path.expanduser(str(pattern))))
        if not matches:
            raise FileNotFoundError(f"no files match {pattern}")
        for item in matches:
            path = Path(item).resolve()
            if path in paths:
                raise ValueError(f"duplicate input file: {path}")
            paths.append(path)
    return paths


def read_desi(paths, *, omega_m=0.315, h=0.674, z_min=0.0, z_max=10.0,
              max_forests=None, pixel_stride=1, delta_field="auto"):
    """Read the DESI DR1 image layout; preserve delta and supplied weights.

    Flat LambdaCDM with Tcmb0=0 converts absorption redshift to Mpc/h. This is
    an explicit configurable fiducial model, not a claim of DESI pipeline parity.
    Limits/stride make demonstration runs small; stride is not pixel rebinning.
    """
    from astropy.cosmology import FlatLambdaCDM
    from astropy.io import fits
    from astropy import units as u

    if not (0 < omega_m < 1 and math.isfinite(h) and h > 0):
        raise ValueError("require 0 < omega_m < 1 and finite h > 0")
    if not (math.isfinite(z_min) and math.isfinite(z_max) and 0 <= z_min < z_max):
        raise ValueError("require finite 0 <= z_min < z_max")
    if pixel_stride < 1 or (max_forests is not None and max_forests < 1):
        raise ValueError("pixel_stride and max_forests must be positive")
    cosmology = FlatLambdaCDM(H0=100*h, Om0=omega_m, Tcmb0=0)
    pieces, seen, provenance = [], set(), []
    for path in expand_inputs(paths):
        if max_forests is not None and len(pieces) >= max_forests:
            break
        with fits.open(path, memmap=False) as hdus:
            names = {hdu.name for hdu in hdus}
            if not {"LAMBDA", "METADATA", "WEIGHT"} <= names:
                raise ValueError(f"{path}: expected DESI/PICCA image-layout FITS")
            field_name = delta_field
            if field_name == "auto":
                field_name = "DELTA_BLIND" if "DELTA_BLIND" in names else "DELTA"
            if field_name not in names:
                raise ValueError(f"{path}: missing {field_name}")
            wave = np.asarray(hdus["LAMBDA"].data, dtype=float)
            meta = hdus["METADATA"]
            required = {"LOS_ID", "RA", "DEC"}
            if not required <= set(meta.columns.names):
                raise ValueError(f"{path}: missing metadata columns {required}")
            if meta.data["LOS_ID"].dtype.kind not in "iu":
                raise ValueError(f"{path}: LOS_ID must be stored as integers")
            waveunit = hdus["LAMBDA"].header.get("BUNIT",
                       hdus["LAMBDA"].header.get("BUNITS", "Angstrom"))
            wave = (wave * u.Unit(waveunit)).to_value(u.AA)
            if wave.ndim != 1 or not np.all(np.isfinite(wave) & (wave > 0)):
                raise ValueError(f"{path}: invalid linear LAMBDA grid")
            values, weights = hdus[field_name].data, hdus["WEIGHT"].data
            expected = (len(meta.data), len(wave))
            if values.shape != expected or weights.shape != expected:
                raise ValueError(f"{path}: expected delta/weight shape {expected}")
            z_abs = wave / LYA_WAVELENGTH - 1
            good_wave = (z_abs > 0) & (z_abs >= z_min) & (z_abs < z_max)
            chi = np.zeros(len(wave))
            chi[good_wave] = cosmology.comoving_distance(z_abs[good_wave]).value * h
            ra_unit = u.Unit(meta.columns["RA"].unit or "rad")
            dec_unit = u.Unit(meta.columns["DEC"].unit or "rad")
            loaded = 0
            for row, data, weight in zip(meta.data, values, weights):
                if max_forests is not None and len(pieces) >= max_forests:
                    break
                identifier = int(row["LOS_ID"])
                if identifier in seen:
                    raise ValueError(f"{path}: duplicate LOS_ID {identifier}; do not "
                                     "combine overlapping releases/delta regions")
                seen.add(identifier)
                ra = (float(row["RA"]) * ra_unit).to_value(u.rad)
                dec = (float(row["DEC"]) * dec_unit).to_value(u.rad)
                if not (math.isfinite(ra) and math.isfinite(dec) and abs(dec) <= np.pi/2):
                    raise ValueError(f"{path}: invalid sky coordinates for {identifier}")
                valid = good_wave & np.isfinite(data) & np.isfinite(weight) & (weight > 0)
                indices = np.flatnonzero(valid)[::pixel_stride]
                if not len(indices):
                    continue
                los = np.array([np.cos(dec)*np.cos(ra), np.cos(dec)*np.sin(ra), np.sin(dec)])
                pieces.append((chi[indices, None]*los, np.array(data[indices]),
                               np.array(weight[indices]),
                               np.full(len(indices), identifier, dtype=np.int64)))
                loaded += 1
            provenance.append(dict(path=str(path), delta_field=field_name,
                                   blinding=str(meta.header.get("BLINDING", "unknown")),
                                   accepted_forests=loaded))
    if not pieces:
        raise ValueError("no valid forest pixels survived the selection")
    arrays = [np.concatenate([piece[i] for piece in pieces]) for i in range(4)]
    return ForestCatalog(*arrays, metadata=dict(
        input_format="DESI/PICCA images", files=provenance, distance_unit="Mpc/h",
        cosmology=dict(model="FlatLambdaCDM", omega_m=omega_m, h=h, Tcmb0=0),
        wavelength_lya_angstrom=LYA_WAVELENGTH, z_min=z_min, z_max=z_max,
        max_forests=max_forests, pixel_stride=pixel_stride,
        weights="input WEIGHT unchanged; zero/negative/nonfinite pixels excluded",
    )).normalized()


def read_ascii(path):
    dtype = [(n, "f8") for n in ("x", "y", "z", "delta", "weights")] + [("ids", "i8")]
    rows = np.loadtxt(path, dtype=dtype, ndmin=1)
    return ForestCatalog(np.column_stack([rows[n] for n in ("x", "y", "z")]),
                         rows["delta"], rows["weights"], rows["ids"],
                         dict(input_format="lya-ascii", path=str(path))).normalized()


def read_npz(path):
    with np.load(path, allow_pickle=False) as data:
        metadata = json.loads(str(data["metadata"])) if "metadata" in data else {}
        return ForestCatalog(*(data[n] for n in ("positions", "delta", "weights", "forest_ids")),
                             metadata=metadata).normalized()


def save_catalog(path, catalog):
    np.savez_compressed(path, positions=catalog.positions, delta=catalog.delta,
                        weights=catalog.weights, forest_ids=catalog.forest_ids,
                        metadata=json.dumps(catalog.metadata))


def synthetic_catalog(forests=8, pixels=12, seed=1234):
    if forests < 3 or pixels < 1:
        raise ValueError("synthetic input needs >= 3 forests and >= 1 pixel/forest")
    rng = np.random.default_rng(seed)
    ra, dec = rng.uniform(-.008, .008, (2, forests))
    los = np.column_stack((np.cos(dec)*np.cos(ra), np.cos(dec)*np.sin(ra), np.sin(dec)))
    chi = 4000 + np.arange(pixels)*2.1 + rng.uniform(-.4, .4, (forests, 1))
    pos = (chi[..., None]*los[:, None, :]).reshape(-1, 3)
    return ForestCatalog(pos, rng.normal(0, .2, len(pos)), rng.uniform(.5, 2, len(pos)),
                         np.repeat(np.arange(forests, dtype=np.int64), pixels),
                         dict(input_format="synthetic", seed=seed,
                              distance_unit="Mpc/h")).normalized()


def collective(comm, function, *, root_only=False):
    """Propagate Python-stage errors before entering the next C collective."""
    result, error = None, None
    if not root_only or comm.rank == 0:
        try:
            result = function()
        except Exception as exc:
            error = f"{type(exc).__name__}: {exc}"
    errors = [e for e in comm.allgather(error) if e]
    if errors:
        raise RuntimeError("; ".join(dict.fromkeys(errors)))
    return result


def broadcast_catalog(comm, catalog):
    meta = comm.bcast(catalog.metadata if comm.rank == 0 else None, root=0)
    arrays = [broadcast_array(comm, getattr(catalog, name) if comm.rank == 0 else None)
              for name in ("positions", "delta", "weights", "forest_ids")]
    return collective(comm, lambda: ForestCatalog(*arrays, metadata=meta).normalized())


def engine_parameters(config, engine, root):
    return dict(searchMethod=engine, rootDir=str(root), iCatalogs="1",
                numberThreads=config.threads, usePeriodic=False, useLogHist=False,
                rangeN=max(config.rp_max, config.rt_max, config.r3_max),
                rminHist=0.0, sizeHistN=4, verbose=0, verbose_log=0, options="",
                lya2RpMax=config.rp_max, lya2RtMax=config.rt_max,
                lya2RpBins=config.rp_bins, lya2RtBins=config.rt_bins,
                lya3RMax=config.r3_max, lya3RBins=config.r3_bins,
                lya3ThetaBins=config.theta_bins, lya3MuBins=config.mu_bins)


def read_products(root, engine):
    spec, result = LYA_ENGINES[engine], {}
    for order in spec.orders:
        key = f"{'1d' if spec.radial else '3d'}_{order}pcf"
        filename, index_count, columns = PRODUCTS[key]
        path = root / filename
        with path.open() as stream:
            has_rows = any(line.strip() and not line.startswith("#") for line in stream)
        table = np.loadtxt(path, ndmin=2) if has_rows else np.empty((0, columns))
        if table.shape[1] != columns or not np.all(np.isfinite(table)):
            raise ValueError(f"{path}: invalid result columns or non-finite results")
        if len(table):
            expected = np.divide(table[:, -2], table[:, -1],
                                 out=np.zeros(len(table)), where=table[:, -1] != 0)
            np.testing.assert_allclose(table[:, -3], expected, rtol=5e-12, atol=1e-14)
            if np.any(table[:, -1] < 0):
                raise ValueError(f"{path}: negative denominator")
        result[key] = table
    return result


def compare_products(reference, candidate, key):
    ndim = PRODUCTS[key][1]
    def keyed(table):
        return {tuple(row[:ndim].astype(int)): row[-3:] for row in table}
    left, right = keyed(reference), keyed(candidate)
    keys = sorted(set(left) | set(right))
    a = np.array([left.get(k, np.zeros(3)) for k in keys]).reshape(-1, 3)
    b = np.array([right.get(k, np.zeros(3)) for k in keys]).reshape(-1, 3)
    difference = b - a
    relative = np.divide(difference[:, 0], a[:, 0], out=np.full(len(keys), np.nan),
                         where=a[:, 0] != 0)
    finite = np.isfinite(relative)
    metrics = dict(max_abs_correlation=float(np.max(abs(difference[:, 0]), initial=0)),
                   max_abs_numerator=float(np.max(abs(difference[:, 1]), initial=0)),
                   max_abs_denominator=float(np.max(abs(difference[:, 2]), initial=0)),
                   max_abs_relative=float(np.max(abs(relative[finite]))) if finite.any() else None,
                   undefined_relative_bins=int((~finite).sum()),
                   occupied_bins=int(np.count_nonzero((a[:, 2] > 0) | (b[:, 2] > 0))))
    rows = np.column_stack((np.asarray(keys).reshape(-1, ndim), a[:, 0], b[:, 0],
                            difference[:, 0], relative))
    return metrics, rows


def projected_plot_data(table, key, config):
    if key == "1d_2pcf":
        shape = (config.rp_bins,)
    elif key == "3d_2pcf":
        shape = (config.rp_bins, config.rt_bins)
    elif key == "1d_3pcf":
        shape = (2*config.r3_bins, 2*config.r3_bins)
    else:
        shape = (config.r3_bins, config.r3_bins)
    num, den = np.zeros(shape), np.zeros(shape)
    if len(table):
        index = tuple(table[:, i].astype(int) for i in range(len(shape)))
        np.add.at(num, index, table[:, -2])
        np.add.at(den, index, table[:, -1])
    return np.divide(num, den, out=np.full(shape, np.nan), where=den > 0)


def make_plots(results, config):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    paths = []
    for key in PRODUCTS:
        members = [(name, r["products"][key]) for name, r in results.items()
                   if key in r["products"]]
        if not members:
            continue
        fig, axes = plt.subplots(1, len(members), figsize=(5*len(members), 4), squeeze=False)
        fields = [projected_plot_data(table, key, config) for _, table in members]
        finite = np.concatenate([v[np.isfinite(v)] for v in fields])
        limit = float(np.max(abs(finite), initial=1e-15))
        for ax, (name, _), values in zip(axes[0], members, fields):
            if values.ndim == 1:
                ax.plot((np.arange(config.rp_bins)+.5)*config.rp_max/config.rp_bins, values)
                ax.set(xlabel="|chi_j - chi_i| [Mpc/h]", ylabel="xi")
            else:
                if key == "3d_2pcf":
                    extent, labels = (0, config.rt_max, 0, config.rp_max), ("r transverse", "r parallel")
                elif key == "1d_3pcf":
                    extent = (-config.r3_max, config.r3_max)*2
                    labels = ("lag 2", "lag 1")
                else:
                    extent, labels = (0, config.r3_max)*2, ("r2", "r1")
                plot = ax.imshow(values, origin="lower", extent=extent, aspect="auto",
                                 cmap="RdBu_r", vmin=-limit, vmax=limit)
                ax.set(xlabel=labels[0]+" [Mpc/h]", ylabel=labels[1]+" [Mpc/h]")
                fig.colorbar(plot, ax=ax, label="xi" if "2pcf" in key else "zeta")
            ax.set_title(name, fontsize=10)
        title = key + (" (angularly summed numerator / denominator)" if key == "3d_3pcf" else "")
        fig.suptitle(title)
        fig.tight_layout()
        path = config.output_dir / f"{key}.png"
        fig.savefig(path, dpi=150)
        plt.close(fig)
        paths.append(str(path))
    return paths


def run_engine_suite(catalog, config, comm=None):
    """Retain one NumPy catalog per process; return root-only result tables."""
    from cyballs import cballs
    if comm is None:
        comm = get_mpi_comm(
            mpi_environment_size() > 1 or any(
                LYA_ENGINES[e].mpi for e in config.engines if e in LYA_ENGINES))
    config = collective(comm, config.validate)
    catalog = collective(comm, catalog.normalized)
    available = collective(comm, lambda: discover_cython_methods(config.engines))
    if any(names != list(config.engines) for names in comm.allgather(available)):
        raise ValueError("selected engines are not compiled identically on every rank")
    forest_count = len(np.unique(catalog.forest_ids))
    required = max(max(LYA_ENGINES[e].orders) for e in config.engines)
    if forest_count < required:
        raise ValueError(f"selected statistics need at least {required} distinct forests")
    config.output_dir = Path(config.output_dir).expanduser().resolve()
    collective(comm, lambda: config.output_dir.mkdir(parents=True, exist_ok=True), root_only=True)
    balls = collective(comm, cballs)
    results = {}
    try:
        if not hasattr(balls, "set_forest_catalog"):
            raise RuntimeError("rebuild cyballs: set_forest_catalog() is required")
        collective(comm, lambda: balls.set_forest_catalog(
            catalog.positions, catalog.delta, catalog.weights, catalog.forest_ids))
        for engine in config.engines:
            root = config.output_dir / engine
            # Never reuse stale histograms from an earlier failed run.
            def prepare():
                root.mkdir(parents=True, exist_ok=True)
                if any((root / spec[0]).exists() for spec in PRODUCTS.values()):
                    raise FileExistsError(f"{root}: results already exist; choose a new output directory")
            collective(comm, prepare, root_only=True)
            spec = LYA_ENGINES[engine]
            participates = spec.mpi or comm.rank == 0
            collective(comm, lambda: balls.set(engine_parameters(config, engine, root))
                       if participates else None)
            comm.barrier()
            if comm.rank == 0:
                print(f"Running {engine}: {catalog.nbody} pixels, {forest_count} forests, "
                      f"{comm.size if spec.mpi else 1} rank(s) x {config.threads} threads",
                      flush=True)
            started = time.perf_counter()
            try:
                collective(comm, lambda: balls.Run(level=["MainLoop"]) if participates else None)
            finally:
                collective(comm, lambda: balls.struct_cleanup() if participates else None)
            elapsed = max(comm.allgather(time.perf_counter() - started))
            product = collective(comm, lambda: read_products(root, engine), root_only=True)
            if comm.rank == 0:
                results[engine] = dict(wall_seconds=elapsed, products=product,
                                       ranks=comm.size if spec.mpi else 1)
        def finish():
            summary = dict(catalog={**catalog.metadata, "pixels": catalog.nbody,
                                    "forests": forest_count},
                           config={**config.__dict__, "output_dir": str(config.output_dir)},
                           catalog_registrations_per_rank=1,
                           timing="wall time: C startup, tree, search, output, cleanup; excludes FITS conversion",
                           engines={e: {k: v for k, v in r.items() if k != "products"}
                                    for e, r in results.items()}, comparisons={})
            for key in PRODUCTS:
                members = [(e, r["products"][key]) for e, r in results.items()
                           if key in r["products"]]
                if not members:
                    continue
                base_name, baseline = members[0]
                for name, table in members[1:]:
                    metrics, rows = compare_products(baseline, table, key)
                    tag = f"{key}__{name}__vs__{base_name}"
                    summary["comparisons"][tag] = metrics
                    np.savetxt(config.output_dir / (tag+".csv"), rows, delimiter=",",
                               header=",".join([f"bin{i}" for i in range(PRODUCTS[key][1])]
                                               + ["reference", "candidate", "difference", "relative"]),
                               comments="")
            summary["plots"] = make_plots(results, config) if config.plots else []
            with (config.output_dir / "summary.json").open("w") as stream:
                json.dump(summary, stream, indent=2, allow_nan=False)
                stream.write("\n")
            return summary
        collective(comm, finish, root_only=True)
    finally:
        collective(comm, balls.clear_catalogs)
    return results


def parse_arguments(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group()
    source.add_argument("--fits", nargs="+", help="DESI DR1 image-layout delta FITS files/globs")
    source.add_argument("--catalog", type=Path, help="NPZ with positions, delta, weights, forest_ids")
    source.add_argument("--ascii", type=Path, help="six columns: x y z delta weight forest_id")
    source.add_argument("--synthetic", action="store_true", help="small generated catalog (default)")
    parser.add_argument("--engine", "--engines", nargs="+", default=["lya-2pcf-omp"])
    parser.add_argument("--statistics", choices=["2pcf", "3pcf", "both"], default="2pcf")
    parser.add_argument("--list-engines", action="store_true")
    parser.add_argument("--output", type=Path, default=Path("Output_lya_all_engines"))
    parser.add_argument("--save-catalog", type=Path, help="cache the selected input as NPZ")
    parser.add_argument("--threads", type=int, default=1)
    for name, default in (("rp-max", 200.), ("rt-max", 200.), ("r3-max", 20.)):
        parser.add_argument("--"+name, type=float, default=default, help="Mpc/h")
    for name, default in (("rp-bins", 50), ("rt-bins", 50), ("r3-bins", 4),
                          ("theta-bins", 4), ("mu-bins", 4)):
        parser.add_argument("--"+name, type=int, default=default)
    parser.add_argument("--omega-m", type=float, default=.315)
    parser.add_argument("--h", type=float, default=.674)
    parser.add_argument("--z-min", type=float, default=0.)
    parser.add_argument("--z-max", type=float, default=10.)
    parser.add_argument("--max-forests", type=int, help="FITS demonstration subset")
    parser.add_argument("--pixel-stride", type=int, default=1, help="FITS subsampling, not rebinning")
    parser.add_argument("--delta-field", choices=["auto", "DELTA", "DELTA_BLIND"], default="auto")
    parser.add_argument("--synthetic-forests", type=int, default=8)
    parser.add_argument("--synthetic-pixels", type=int, default=12)
    parser.add_argument("--seed", type=int, default=1234)
    parser.add_argument("--max-hist-mib", type=float, default=1024.)
    parser.add_argument("--no-plots", action="store_true")
    parser.add_argument("--mpi-ranks", type=int, default=1)
    parser.add_argument("--mpiexec", help="launcher built against the same MPI as cyballs/mpi4py")
    parser.add_argument("--mpi-extra-arg", action="append", default=[])
    return parser.parse_args(argv)


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    args = parse_arguments(argv)
    if args.mpi_ranks < 1:
        raise ValueError("--mpi-ranks must be positive")
    if args.mpi_ranks > 1 and mpi_environment_size() == 1 and not os.environ.get(MPI_CHILD):
        launcher = args.mpiexec or shutil.which("mpiexec") or shutil.which("mpirun")
        if not launcher:
            raise RuntimeError("mpiexec not found; specify --mpiexec")
        command = [launcher]
        for extra in args.mpi_extra_arg:
            command.extend(shlex.split(extra))
        command += ["-n", str(args.mpi_ranks), sys.executable, str(Path(__file__).resolve()), *argv]
        return subprocess.run(command, env={**os.environ, MPI_CHILD: "1"}, check=False).returncode
    # Importing cyballs does not initialize MPI; multi-rank execution does.
    comm = get_mpi_comm(mpi_environment_size() > 1)
    available = collective(comm, lambda: discover_cython_methods(list(LYA_ENGINES)))
    if args.list_engines:
        if comm.rank == 0:
            for name, spec in LYA_ENGINES.items():
                print(f"{name:28s} {'available' if name in available else 'not built':10s} "
                      f"{'radial' if spec.radial else '3D':6s} {spec.orders}")
            print(EXCLUDED)
        return 0
    engines = collective(comm, lambda: resolve_engines(args.engine, available, args.statistics))
    if comm.size == 1 and any(LYA_ENGINES[e].mpi for e in engines):
        comm = get_mpi_comm(True)
    config = RunConfig(engines, args.output, threads=args.threads, rp_max=args.rp_max,
                       rt_max=args.rt_max, rp_bins=args.rp_bins, rt_bins=args.rt_bins,
                       r3_max=args.r3_max, r3_bins=args.r3_bins, theta_bins=args.theta_bins,
                       mu_bins=args.mu_bins, max_hist_mib=args.max_hist_mib,
                       plots=not args.no_plots)
    collective(comm, config.validate)
    def load():
        if not args.fits and (args.max_forests is not None or args.pixel_stride != 1
                              or args.z_min != 0 or args.z_max != 10
                              or args.omega_m != .315 or args.h != .674
                              or args.delta_field != "auto"):
            raise ValueError("FITS selection/cosmology options require --fits; "
                             "ASCII and NPZ coordinates are already comoving")
        if args.fits:
            catalog = read_desi(args.fits, omega_m=args.omega_m, h=args.h,
                                z_min=args.z_min, z_max=args.z_max, max_forests=args.max_forests,
                                pixel_stride=args.pixel_stride, delta_field=args.delta_field)
        elif args.catalog:
            catalog = read_npz(args.catalog)
        elif args.ascii:
            catalog = read_ascii(args.ascii)
        else:
            catalog = synthetic_catalog(args.synthetic_forests, args.synthetic_pixels, args.seed)
        if args.save_catalog:
            args.save_catalog.parent.mkdir(parents=True, exist_ok=True)
            save_catalog(args.save_catalog, catalog)
        return catalog
    catalog = collective(comm, load, root_only=True)
    catalog = broadcast_catalog(comm, catalog)
    if comm.rank == 0:
        print("Catalog loaded once; retained NumPy input will be reused by every engine.")
        for source in catalog.metadata.get("files", []):
            print(f"DESI field={source['delta_field']}, BLINDING={source['blinding']}")
        if any(3 in LYA_ENGINES[e].orders for e in engines):
            print("3PCF can be expensive on dense forests; start with --max-forests and --pixel-stride.")
    run_engine_suite(catalog, config, comm)
    if comm.rank == 0:
        print(f"Results: {config.output_dir / 'summary.json'}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
