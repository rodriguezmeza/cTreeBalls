#!/usr/bin/env python3
"""Run active cTreeBalls Ly-alpha forest engines on one retained catalog.

DESI and eBOSS/PICCA FITS, six-column lya-ascii, NPZ and synthetic input are supported.
MPI rank zero reads once, broadcasts arrays once, and all ranks register a
forest catalog with cyballs. Radial, anisotropic 3D, and multipole estimators
are distinct.
See tests/python/README_lya_corr_all_engines.md for commands and scientific limitations.
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
    flatten_radial_matrix, mpi_environment_size, aggregate_rank_timings,
    _timing_metadata, timing_summary, write_timing_report,
)

PROJECT_ROOT = Path(__file__).resolve().parents[2]
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
    tree: bool = False
    family: str = "lya"
    supports_smooth_pivot: bool = False


LYA_ENGINES = {
    name + "-" + parallel: EngineSpec(radial, orders, parallel == "mpi", tree)
    for parallel in ("omp", "mpi")
    for name, radial, orders, tree in (
        ("lya-2pcf", False, (2,), False),
        ("lya-3pcf", False, (3,), False),
        ("lya-2pcf-3pcf", False, (2, 3), False),
        ("lya-1d-2pcf", True, (2,), False),
        ("lya-1d-3pcf", True, (3,), False),
        ("lya-1d-2pcf-3pcf", True, (2, 3), False),
        ("lya-1d-tree-2pcf", True, (2,), True),
        ("lya-1d-tree-3pcf", True, (3,), True),
    )
}
LYA_ENGINES.update({
    "lya-los-tree-2pcf-omp": EngineSpec(False, (2,), False, tree=True),
    "lya-los-tree-3pcf-omp": EngineSpec(False, (3,), False, tree=True),
    "lya-los-tree-2pcf-3pcf-omp": EngineSpec(False, (2, 3), False, tree=True),
    "octree-3pcf-3d-omp": EngineSpec(False, (3,), False, family="multipole"),
    "octree-3pcf-3d-mpi": EngineSpec(False, (3,), True, family="multipole"),
    "lya-1d-tree-same-los-2pcf-omp": EngineSpec(
        True, (2,), False, tree=True, family="same_los"
    ),
})
EXCLUDED = (
    "Angular convergence/shear and point-count engines are not forest estimators: "
    "their geometry or field contracts differ."
)
FOREST_SELECTION_NOTE = (
    "Angular read-mask/edge-corrections and scalar dual-node frontier controls "
    "are not forest options; select forests/pixels before set_forest_catalog()."
)
INCOMPATIBLE_ENGINE_REASONS = {
    name: (
        "requires gamma1/gamma2 spin-2 data and spherical or flat-sky transport; "
        "use tests/python/shear_corr_all_engines.py"
    )
    for name in (
        "octree-shear-sphere-2balls-omp",
        "kdtree-shear-sphere-2balls-omp",
        "balltree-shear-sphere-2balls-omp",
    )
}
PRODUCTS = {
    "3d_2pcf": ("histXi2pcf_lya.txt", 2, 7),
    "3d_3pcf": ("histZetaM_lya5d.txt", 5, 13),
    "1d_2pcf": ("histXi2pcf_lya1d.txt", 1, 5),
    "1d_3pcf": ("histZetaM_lya1d.txt", 2, 7),
    "1d_same_los_2pcf": ("histXi2pcf_lya1d_same_los.txt", 1, 5),
    "multipole_3pcf": ("histZetaM_3d.txt", 3, 8),
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
    multipole_lmax: int = 3
    multipole_rmin: float = 0.0
    max_hist_mib: float = 1024.0
    plots: bool = True
    flatten_plots: bool = True
    lya2pcf_source: Path | None = None
    reference_covariance: bool = False
    reference_nside: int = 32
    reference_timeout: float = 0.0
    rtol: float = 1e-8
    atol: float = 1e-12
    relative_floor: float = 1e-12
    fail_on_mismatch: bool = False
    analysis: bool = True
    covariance: Path | None = None
    distortion_matrix: Path | None = None
    model_correlation: Path | None = None
    wedge_bins: int = 50
    wedge_mu_edges: tuple[float, ...] = (0., .5, .8, .95, 1.)
    wedge_subsamples: int = 10
    wedge_r_max: float | None = None

    def validate(self):
        if not self.engines or any(e not in LYA_ENGINES for e in self.engines):
            raise ValueError("select at least one registered Ly-alpha engine")
        for name in ("threads", "rp_bins", "rt_bins", "r3_bins", "theta_bins", "mu_bins"):
            value = getattr(self, name)
            if not isinstance(value, int) or not 1 <= value <= 100000:
                raise ValueError(f"{name} must be an integer in [1, 100000]")
        if not isinstance(self.multipole_lmax, int) or not 0 <= self.multipole_lmax <= 10:
            raise ValueError("multipole_lmax must be an integer in [0, 10]")
        for name in ("rp_max", "rt_max", "r3_max", "max_hist_mib"):
            if not math.isfinite(getattr(self, name)) or getattr(self, name) <= 0:
                raise ValueError(f"{name} must be positive and finite")
        if (not math.isfinite(self.multipole_rmin) or self.multipole_rmin < 0
                or self.multipole_rmin >= self.r3_max):
            raise ValueError("require 0 <= multipole_rmin < r3_max")
        if not math.isfinite(math.hypot(self.rp_max, self.rt_max)):
            raise ValueError("separation limits overflow")
        for name in ("rtol", "atol", "relative_floor", "reference_timeout"):
            if not math.isfinite(getattr(self, name)) or getattr(self, name) < 0:
                raise ValueError(f"{name} must be finite and non-negative")
        if self.reference_nside < 1 or self.reference_nside > 8192 or self.reference_nside & (self.reference_nside-1):
            raise ValueError("reference_nside must be a power of two in [1, 8192]")
        if not 1 <= self.wedge_bins <= 10000 or not 1 <= self.wedge_subsamples <= 100:
            raise ValueError("require wedge_bins in [1, 10000] and wedge_subsamples in [1, 100]")
        edges = np.asarray(self.wedge_mu_edges)
        if (edges.ndim != 1 or not 2 <= len(edges) <= 17 or not np.isfinite(edges).all()
                or edges[0] != 0 or edges[-1] != 1 or np.any(np.diff(edges) <= 0)):
            raise ValueError("wedge_mu_edges must increase strictly from 0 to 1 (at most 16 wedges)")
        if self.wedge_r_max is not None and (not math.isfinite(self.wedge_r_max) or self.wedge_r_max <= 0):
            raise ValueError("wedge_r_max must be positive and finite")
        has_2pcf = any(2 in LYA_ENGINES[e].orders for e in self.engines)
        has_3d_2pcf = any(2 in LYA_ENGINES[e].orders and not LYA_ENGINES[e].radial for e in self.engines)
        if self.lya2pcf_source:
            if not has_2pcf:
                raise ValueError("lya2pcf reference requires a selected 2PCF engine; no upstream 3PCF reference is available")
            self.lya2pcf_source = Path(self.lya2pcf_source).expanduser().resolve()
            for name in ("parameters.py", "correlation_procedures_cpu.py", "post_processing.py"):
                if not (self.lya2pcf_source/name).is_file():
                    raise ValueError(f"--lya2pcf-source is missing {name}")
        if self.reference_covariance and (not self.lya2pcf_source or not has_3d_2pcf):
            raise ValueError("reference covariance requires --lya2pcf-source and a 3D 2PCF engine")
        if bool(self.distortion_matrix) != bool(self.model_correlation):
            raise ValueError("use --distortion-matrix together with --model-correlation (forward modelling)")
        if (self.covariance or self.distortion_matrix or self.reference_covariance) and (not self.analysis or not has_3d_2pcf):
            raise ValueError("covariance/distortion analysis requires a 3D 2PCF engine and enabled analysis")
        if self.analysis and has_3d_2pcf:
            from lya_analysis import analysis_workspace_bytes
            if analysis_workspace_bytes(self) > self.max_hist_mib*2**20:
                raise ValueError("2PCF analysis workspace exceeds --max-hist-mib; reduce bins or increase budget")
        for name in ("covariance", "distortion_matrix", "model_correlation"):
            if getattr(self, name):
                path = Path(getattr(self, name)).expanduser().resolve()
                if not path.is_file():
                    raise ValueError(f"{name}: input file does not exist: {path}")
                setattr(self, name, path)
        for engine in self.engines:
            spec = LYA_ENGINES[engine]
            bins2 = self.rp_bins * (1 if spec.radial else self.rt_bins)
            if spec.family == "multipole":
                bins3 = (self.multipole_lmax + 1)*self.r3_bins**2
            else:
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
        if token in ("all", "all-omp", "all-mpi", "all-1d", "all-3d", "all-tree",
                     "all-multipole"):
            candidates = [name for name, spec in LYA_ENGINES.items()
                          if name in available
                          and (token != "all-omp" or not spec.mpi)
                          and (token != "all-mpi" or spec.mpi)
                          and (token != "all-1d" or spec.radial)
                          and (token != "all-3d" or not spec.radial)
                          and (token != "all-tree" or spec.tree)
                          and (token != "all-multipole" or spec.family == "multipole")]
            if statistics != "both":
                candidates = [n for n in candidates
                              if LYA_ENGINES[n].orders == (int(statistics[0]),)]
        else:
            if token not in LYA_ENGINES:
                if token in INCOMPATIBLE_ENGINE_REASONS:
                    raise ValueError(f"{token}: {INCOMPATIBLE_ENGINE_REASONS[token]}")
                raise ValueError(f"{token}: not a registered forest engine. {EXCLUDED}")
            if token not in available:
                raise ValueError(
                    f"{token} is unavailable: compile its cTreeBalls addon"
                )
            if statistics != "both" and LYA_ENGINES[token].orders != (int(statistics[0]),):
                raise ValueError(f"{token} needs --statistics=both or its matching order")
            candidates = [token]
        selected.extend(n for n in candidates if n not in selected)
    if not selected:
        raise ValueError("no requested engines are available; enable the corresponding "
                         "LYAFOREST/OCTREE3PCF3D addon")
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
              max_forests=None, pixel_stride=1, delta_field="auto", project_delta=False,
              redshift_weight_exponent=0., weight_z_ref=2.25):
    """Read the DESI DR1 image layout; preserve delta and supplied weights.

    Flat LambdaCDM with Tcmb0=0 converts absorption redshift to Mpc/h. This is
    an explicit configurable fiducial model, not a claim of DESI pipeline parity.
    Limits/stride make demonstration runs small; stride is not pixel rebinning.
    """
    from astropy.cosmology import FlatLambdaCDM
    from astropy.io import fits
    from astropy import units as u
    from lya_fits import preprocess_forest

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
                retained_delta, retained_weight = preprocess_forest(
                    data[indices], weight[indices], np.log10(wave[indices]),
                    project_delta=project_delta, redshift_weight_exponent=redshift_weight_exponent,
                    weight_z_ref=weight_z_ref)
                pieces.append((chi[indices, None]*los, retained_delta, retained_weight,
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
        project_delta=project_delta, projection_domain="retained pixels after cuts/stride",
        redshift_weight_exponent=redshift_weight_exponent, weight_z_ref=weight_z_ref,
        weights=("input WEIGHT unchanged; zero/negative/nonfinite pixels excluded" if redshift_weight_exponent == 0
                 else "input WEIGHT times ((1+z)/(1+z_ref))**exponent"),
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
    spec = LYA_ENGINES[engine]
    if spec.family == "multipole":
        return dict(searchMethod=engine, rootDir=str(root), iCatalogs="1",
                    numberThreads=config.threads, usePeriodic=False,
                    useLogHist=False, rangeN=config.r3_max,
                    rminHist=config.multipole_rmin,
                    sizeHistN=config.r3_bins, mChebyshev=config.multipole_lmax,
                    verbose=0, verbose_log=0,
                    options="only-3pcf-3d,exclude-all-same-los,no-smooth-pivot")
    return dict(searchMethod=engine, rootDir=str(root), iCatalogs="1",
                numberThreads=config.threads, usePeriodic=False, useLogHist=False,
                rangeN=max(config.rp_max, config.rt_max, config.r3_max),
                rminHist=0.0, sizeHistN=4, verbose=0, verbose_log=0,
                options="no-smooth-pivot",
                lya2RpMax=config.rp_max, lya2RtMax=config.rt_max,
                lya2RpBins=config.rp_bins, lya2RtBins=config.rt_bins,
                lya3RMax=config.r3_max, lya3RBins=config.r3_bins,
                lya3ThetaBins=config.theta_bins, lya3MuBins=config.mu_bins)


def read_products(root, engine):
    spec, result = LYA_ENGINES[engine], {}
    for order in spec.orders:
        key = (f"multipole_{order}pcf" if spec.family == "multipole"
               else (f"1d_same_los_{order}pcf" if spec.family == "same_los"
                     else f"{'1d' if spec.radial else '3d'}_{order}pcf"))
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
        if key == "multipole_3pcf":
            table = table[table[:, 1] < table[:, 2]]
        result[key] = table
    return result


def compare_products(reference, candidate, key, *, relative_floor=0., rtol=1e-8, atol=1e-12):
    ndim = PRODUCTS[key][1]
    def keyed(table):
        values = {tuple(row[:ndim].astype(int)): row[-3:] for row in table}
        if len(values) != len(table):
            raise ValueError(f"{key}: duplicate histogram bin indices")
        return values
    left, right = keyed(reference), keyed(candidate)
    keys = sorted(set(left) | set(right))
    a = np.array([left.get(k, np.zeros(3)) for k in keys]).reshape(-1, 3)
    b = np.array([right.get(k, np.zeros(3)) for k in keys]).reshape(-1, 3)
    difference = b - a
    relative = np.divide(difference[:, 0], a[:, 0], out=np.full(len(keys), np.nan),
                         where=np.abs(a[:, 0]) > relative_floor)
    finite = np.isfinite(relative)
    metrics = dict(max_abs_correlation=float(np.max(abs(difference[:, 0]), initial=0)),
                   max_abs_numerator=float(np.max(abs(difference[:, 1]), initial=0)),
                   max_abs_denominator=float(np.max(abs(difference[:, 2]), initial=0)),
                   max_abs_relative=float(np.max(abs(relative[finite]))) if finite.any() else None,
                   undefined_relative_bins=int((~finite).sum()),
                   occupied_bins=int(np.count_nonzero((a[:, 2] > 0) | (b[:, 2] > 0))),
                   occupancy_mismatches=int(np.count_nonzero((a[:, 2] > 0) != (b[:, 2] > 0))),
                   passed_correlation=bool(np.allclose(b[:, 0], a[:, 0], rtol=rtol, atol=atol)),
                   passed_raw_sums=bool(np.allclose(b[:, 1:], a[:, 1:], rtol=rtol, atol=atol)),
                   rtol=rtol, atol=atol, relative_floor=relative_floor)
    metrics["passed"] = (metrics["passed_correlation"] and metrics["passed_raw_sums"]
                         and metrics["occupancy_mismatches"] == 0)
    rows = np.column_stack((np.asarray(keys).reshape(-1, ndim), a[:, 0], b[:, 0],
                            difference[:, 0], relative))
    return metrics, rows


def projected_plot_data(table, key, config):
    if key == "multipole_3pcf":
        table = table[table[:, 0] == 0]
    if key in ("1d_2pcf", "1d_same_los_2pcf"):
        shape = (config.rp_bins,)
    elif key == "3d_2pcf":
        shape = (config.rp_bins, config.rt_bins)
    elif key == "1d_3pcf":
        shape = (2*config.r3_bins, 2*config.r3_bins)
    elif key == "multipole_3pcf":
        shape = (config.r3_bins, config.r3_bins)
    else:
        shape = (config.r3_bins, config.r3_bins)
    num, den = np.zeros(shape), np.zeros(shape)
    if len(table):
        if key == "multipole_3pcf":
            index = (table[:, 1].astype(int)-1, table[:, 2].astype(int)-1)
        else:
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
        columns = min(3, len(members))
        rows = math.ceil(len(members)/columns)
        fig, axes = plt.subplots(rows, columns, figsize=(5*columns, 4*rows), squeeze=False)
        fields = [projected_plot_data(table, key, config) for _, table in members]
        finite_parts = [v[np.isfinite(v)] for v in fields if np.any(np.isfinite(v))]
        finite = np.concatenate(finite_parts) if finite_parts else np.array([1e-15])
        limit = float(np.max(abs(finite), initial=1e-15))
        for ax, (name, _), values in zip(axes.flat, members, fields):
            if values.ndim == 1:
                ax.plot((np.arange(config.rp_bins)+.5)*config.rp_max/config.rp_bins, values)
                ax.set(xlabel="|chi_j - chi_i| [Mpc/h]", ylabel="xi")
            else:
                if key == "3d_2pcf":
                    extent, labels = (0, config.rt_max, 0, config.rp_max), ("r transverse", "r parallel")
                elif key == "1d_3pcf":
                    extent = (-config.r3_max, config.r3_max)*2
                    labels = ("lag 2", "lag 1")
                elif key == "multipole_3pcf":
                    extent = (config.multipole_rmin, config.r3_max)*2
                    labels = ("r2", "r1")
                else:
                    extent, labels = (0, config.r3_max)*2, ("r2", "r1")
                plot = ax.imshow(values, origin="lower", extent=extent, aspect="auto",
                                 cmap="RdBu_r", vmin=-limit, vmax=limit)
                ax.set(xlabel=labels[0]+" [Mpc/h]", ylabel=labels[1]+" [Mpc/h]")
                fig.colorbar(plot, ax=ax, label="xi" if "2pcf" in key else "zeta")
            ax.set_title(name, fontsize=10)
        for ax in list(axes.flat)[len(members):]:
            ax.set_visible(False)
        title = key + (" (angular-summed numerator / denominator)" if key == "3d_3pcf"
                       else " (ell=0 monopole)" if key == "multipole_3pcf" else "")
        fig.suptitle(title)
        fig.tight_layout()
        path = config.output_dir / f"{key}.png"
        fig.savefig(path, dpi=150)
        plt.close(fig)
        paths.append(str(path))
        if config.flatten_plots and "3pcf" in key and all(v.ndim == 2 for v in fields):
            fig, ax = plt.subplots(figsize=(9, 4.5))
            for (name, _), values in zip(members, fields):
                flattened = flatten_radial_matrix(values)
                ax.plot(np.arange(flattened.size), flattened, label=name)
            width = fields[0].shape[1]
            for boundary in range(width, fields[0].size, width):
                ax.axvline(boundary - 0.5, color="0.82", linewidth=0.55)
            ax.set(
                xlabel="flattened radial-bin index (bin 1 major, bin 2 minor)",
                ylabel="zeta", title=title+": flattened radial-bin matrix",
            )
            ax.grid(True, axis="y", linestyle=":", alpha=0.45)
            ax.legend(fontsize=8)
            fig.tight_layout()
            flat_path = config.output_dir / f"{key}_flattened.png"
            fig.savefig(flat_path, dpi=150)
            plt.close(fig)
            paths.append(str(flat_path))
    return paths


def run_engine_suite(catalog, config, comm=None):
    """Retain one NumPy catalog per process; return root-only result tables."""
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
    required = max(1 if LYA_ENGINES[e].family == "same_los"
                   else max(LYA_ENGINES[e].orders) for e in config.engines)
    if forest_count < required:
        raise ValueError(f"selected statistics need at least {required} distinct forests")
    config.output_dir = Path(config.output_dir).expanduser().resolve()
    collective(comm, lambda: config.output_dir.mkdir(parents=True, exist_ok=True), root_only=True)
    balls = None
    from cyballs import cballs
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
                if any(root.iterdir()):
                    raise FileExistsError(f"{root}: results already exist; choose a new output directory")
            collective(comm, prepare, root_only=True)
            spec = LYA_ENGINES[engine]
            participates = spec.mpi or comm.rank == 0
            setup_wall = setup_cpu = compute_wall = compute_cpu = native_cpu = 0.0
            native_timings = None
            def setup():
                nonlocal setup_wall, setup_cpu
                if participates:
                    wall, cpu = time.perf_counter(), time.process_time()
                    balls.set(engine_parameters(config, engine, root))
                    balls.Run(level=["SetNumberThreads"])
                    setup_wall, setup_cpu = time.perf_counter()-wall, time.process_time()-cpu
            collective(comm, setup)
            comm.barrier()
            if comm.rank == 0:
                print(f"Running {engine}: {catalog.nbody} pixels, {forest_count} forests, "
                      f"{comm.size if spec.mpi else 1} rank(s) x {config.threads} threads",
                      flush=True)
            started = time.perf_counter()
            def compute():
                nonlocal compute_wall, compute_cpu, native_cpu, native_timings
                if participates:
                    wall, cpu = time.perf_counter(), time.process_time()
                    balls.Run(level=["MainLoop"])
                    compute_wall, compute_cpu = time.perf_counter()-wall, time.process_time()-cpu
                    native_cpu = float(balls.getCPUTime())
                    native_timings = balls.getTimings()
            try:
                collective(comm, compute)
            finally:
                collective(comm, lambda: balls.struct_cleanup() if participates else None)
            elapsed = max(comm.allgather(time.perf_counter() - started))
            rank_timing = _timing_metadata(
                setup_wall, setup_cpu, compute_wall, compute_cpu,
                "cTreeBalls parameter/thread setup plus MainLoop including native output; "
                "catalog registration, Python result loading and cleanup excluded")
            rank_timing.update(rank=comm.rank, native_reported_cpu_time=native_cpu)
            if participates:
                rank_timing.update(
                    native_mainloop_wall_time=float(native_timings["wall_seconds"]),
                    native_mainloop_cpu_time=float(native_timings["process_cpu_seconds"]))
            rows = comm.allgather(rank_timing if participates else None)
            timings = aggregate_rank_timings([row for row in rows if row is not None])
            product = collective(comm, lambda: read_products(root, engine), root_only=True)
            if comm.rank == 0:
                results[engine] = dict(wall_seconds=elapsed, products=product,
                                       ranks=comm.size if spec.mpi else 1,
                                       smooth_pivot="unsupported",
                                       provenance=balls.getRunMetadata(),
                                       threads_per_rank=config.threads,
                                       parameters=engine_parameters(config, engine, root))
                results[engine].update(timings)
                print(f"  completed: compute wall {timings['compute_wall_time']:.6f}s, "
                      f"summed process CPU {timings['compute_cpu_time']:.6f}s", flush=True)
        def finish():
            from lya_analysis import write_comparisons, analyse_2pcf
            native_timings = timing_summary(results)
            write_timing_report(config.output_dir / "timing_report.txt", native_timings)
            if config.lya2pcf_source:
                from lya_reference import run_references
                families = {key for result in results.values() for key in result["products"]}
                results.update(run_references(catalog, config, families))
            config_json = {
                key: str(value) if isinstance(value, Path) else value
                for key, value in config.__dict__.items()
            }
            summary = dict(catalog={**catalog.metadata, "pixels": catalog.nbody,
                                    "forests": forest_count},
                           config=config_json,
                           catalog_registrations_per_rank=1,
                           timings=native_timings,
                           timing=("native wall_seconds includes MainLoop, output, synchronization and "
                                   "cleanup; FITS conversion is excluded. Reference wall_seconds is "
                                   "warmed CPU pair traversal/reduction; worker_wall_seconds includes "
                                   "imports, JIT, analysis and I/O. These are different scopes, not a speedup ratio."),
                           selection_contract=FOREST_SELECTION_NOTE,
                           engines={e: {k: v for k, v in r.items() if k != "products"}
                                    for e, r in results.items()}, comparisons={})
            summary["comparisons"], difference_paths = write_comparisons(results, config)
            summary["comparison_contract"] = (
                "All pairs within each estimator family only; 3D, radial-only, same-LOS and "
                "multipole products are not interchangeable. No upstream 3PCF reference is run.")
            summary["plots"] = (make_plots(results, config) if config.plots else []) + difference_paths
            if config.analysis:
                summary["analysis_2pcf"], analysis_paths = analyse_2pcf(results, config)
                summary["plots"].extend(analysis_paths)
            failed = [tag for tag, result in summary["comparisons"].items() if not result["passed"]]
            summary["validation"] = dict(passed=not failed, failed_comparisons=failed,
                                         comparisons=len(summary["comparisons"]))
            with (config.output_dir / "summary.json").open("w") as stream:
                json.dump(summary, stream, indent=2, allow_nan=False)
                stream.write("\n")
            print(f"Comparisons: {len(summary['comparisons'])-len(failed)}/{len(summary['comparisons'])} passed "
                  f"(rtol={config.rtol:g}, atol={config.atol:g})", flush=True)
            if failed and config.fail_on_mismatch:
                raise ValueError(f"{len(failed)} comparison(s) failed; inspect {config.output_dir/'summary.json'}")
            return summary
        collective(comm, finish, root_only=True)
    finally:
        collective(comm, balls.clear_catalogs)
    return results


def parse_arguments(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group()
    source.add_argument("--fits", nargs="+", help="DESI image-layout or eBOSS/PICCA forest-HDU delta FITS files/globs")
    source.add_argument("--catalog", type=Path, help="NPZ with positions, delta, weights, forest_ids")
    source.add_argument("--ascii", type=Path, help="six columns: x y z delta weight forest_id")
    source.add_argument("--synthetic", action="store_true", help="small generated catalog (default)")
    parser.add_argument(
        "--engine", "--engines", nargs="+", default=["lya-2pcf-omp"],
        help=("engine names or all/all-omp/all-mpi/all-1d/all-3d/all-tree/"
              "all-multipole"),
    )
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
    parser.add_argument("--multipole-lmax", type=int, default=3)
    parser.add_argument("--multipole-rmin", type=float, default=0.0, help="Mpc/h")
    parser.add_argument("--omega-m", type=float, default=.315)
    parser.add_argument("--h", type=float, default=.674)
    parser.add_argument("--z-min", type=float, default=0.)
    parser.add_argument("--z-max", type=float, default=10.)
    parser.add_argument("--max-forests", type=int, help="FITS demonstration subset")
    parser.add_argument("--pixel-stride", type=int, default=1, help="FITS subsampling, not rebinning")
    parser.add_argument("--delta-field", choices=["auto", "DELTA", "DELTA_BLIND"], default="auto")
    parser.add_argument("--fits-layout", choices=["auto", "desi", "eboss"], default="auto")
    parser.add_argument("--eboss-angle-unit", choices=["rad", "deg"], default="rad")
    parser.add_argument("--project-delta", action="store_true", help="remove weighted mean and log-wavelength slope on retained FITS pixels")
    parser.add_argument("--redshift-weight-exponent", type=float, default=0., help="optional WEIGHT multiplier ((1+z)/(1+z_ref))**exponent")
    parser.add_argument("--weight-z-ref", type=float, default=2.25)
    parser.add_argument("--lya2pcf-source", type=Path, help="enable external CPU pair-kernel references for every selected 2PCF family")
    parser.add_argument("--reference-covariance", action="store_true", help="use upstream weighted sky-subsampling covariance for 3D 2PCF")
    parser.add_argument("--reference-nside", type=int, default=32)
    parser.add_argument("--reference-timeout", type=float, default=0., help="seconds per reference worker; 0 means no timeout")
    parser.add_argument("--rtol", type=float, default=1e-8)
    parser.add_argument("--atol", type=float, default=1e-12)
    parser.add_argument("--relative-floor", type=float, default=1e-12, help="omit relative errors for smaller absolute reference correlations")
    parser.add_argument("--fail-on-mismatch", action="store_true")
    parser.add_argument("--no-analysis", action="store_true", help="omit 3D 2PCF archives, wedges and covariance/model analysis")
    parser.add_argument("--covariance", type=Path, help="shared 3D 2PCF covariance: NPY/NPZ or FITS CO")
    parser.add_argument("--distortion-matrix", type=Path, help="forward-model matrix: NPY/NPZ or FITS DM (not an inverse correction)")
    parser.add_argument("--model-correlation", type=Path, help="unprojected 3D model on matching bins: NPY/NPZ or FITS DA")
    parser.add_argument("--wedge-bins", type=int, default=50)
    parser.add_argument("--wedge-mu-edges", default="0,0.5,0.8,0.95,1")
    parser.add_argument("--wedge-subsamples", type=int, default=10, help="sub-bin samples per dimension; 100 reproduces notebook resolution")
    parser.add_argument("--wedge-r-max", type=float, help="Mpc/h; default min(rp_max, rt_max)")
    parser.add_argument("--synthetic-forests", type=int, default=8)
    parser.add_argument("--synthetic-pixels", type=int, default=12)
    parser.add_argument("--seed", type=int, default=1234)
    parser.add_argument("--max-hist-mib", type=float, default=1024.)
    parser.add_argument("--no-plots", action="store_true")
    parser.add_argument(
        "--no-flatten-plots", action="store_true",
        help="omit flattened radial-bin views for 3PCF products",
    )
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
    native_candidates = list(LYA_ENGINES)
    probe_candidates = [*native_candidates, *INCOMPATIBLE_ENGINE_REASONS]
    discovered = collective(comm, lambda: discover_cython_methods(probe_candidates))
    available = [name for name in native_candidates if name in discovered]
    if args.list_engines:
        if comm.rank == 0:
            for name, spec in LYA_ENGINES.items():
                geometry = ("1D-tree" if spec.tree else
                            ("radial" if spec.radial else "3D"))
                smooth = "unsupported"
                print(f"{name:28s} {'available' if name in available else 'not built':10s} "
                      f"{geometry:7s} {spec.orders} smooth={smooth} "
                      "mask=preselect edge=not-applicable")
            rejected = [name for name in INCOMPATIBLE_ENGINE_REASONS if name in discovered]
            if rejected:
                print("\nBuilt spin-2 engines intentionally rejected by the forest driver:")
                for name in rejected:
                    print(f"  {name:42s} {INCOMPATIBLE_ENGINE_REASONS[name]}")
            print(EXCLUDED)
            print(FOREST_SELECTION_NOTE)
        return 0
    engines = collective(comm, lambda: resolve_engines(args.engine, available, args.statistics))
    if comm.size == 1 and any(LYA_ENGINES[e].mpi for e in engines):
        comm = get_mpi_comm(True)
    config = RunConfig(engines, args.output, threads=args.threads, rp_max=args.rp_max,
                       rt_max=args.rt_max, rp_bins=args.rp_bins, rt_bins=args.rt_bins,
                       r3_max=args.r3_max, r3_bins=args.r3_bins, theta_bins=args.theta_bins,
                       mu_bins=args.mu_bins, multipole_lmax=args.multipole_lmax,
                       multipole_rmin=args.multipole_rmin,
                       max_hist_mib=args.max_hist_mib,
                       plots=not args.no_plots,
                       flatten_plots=not args.no_flatten_plots,
                       lya2pcf_source=args.lya2pcf_source, reference_covariance=args.reference_covariance,
                       reference_nside=args.reference_nside, reference_timeout=args.reference_timeout,
                       rtol=args.rtol, atol=args.atol, relative_floor=args.relative_floor,
                       fail_on_mismatch=args.fail_on_mismatch, analysis=not args.no_analysis,
                       covariance=args.covariance, distortion_matrix=args.distortion_matrix,
                       model_correlation=args.model_correlation, wedge_bins=args.wedge_bins,
                       wedge_mu_edges=tuple(float(x) for x in args.wedge_mu_edges.split(",")),
                       wedge_subsamples=args.wedge_subsamples, wedge_r_max=args.wedge_r_max)
    collective(comm, config.validate)
    def load():
        if not args.fits and (args.max_forests is not None or args.pixel_stride != 1
                              or args.z_min != 0 or args.z_max != 10
                              or args.omega_m != .315 or args.h != .674
                              or args.delta_field != "auto" or args.fits_layout != "auto"
                              or args.eboss_angle_unit != "rad" or args.project_delta
                              or args.redshift_weight_exponent != 0 or args.weight_z_ref != 2.25):
            raise ValueError("FITS selection/cosmology options require --fits; "
                             "ASCII and NPZ coordinates are already comoving")
        if args.fits:
            from lya_fits import read_fits
            catalog = read_fits(args.fits, omega_m=args.omega_m, h=args.h,
                                z_min=args.z_min, z_max=args.z_max, max_forests=args.max_forests,
                                pixel_stride=args.pixel_stride, delta_field=args.delta_field,
                                fits_layout=args.fits_layout, eboss_angle_unit=args.eboss_angle_unit,
                                project_delta=args.project_delta,
                                redshift_weight_exponent=args.redshift_weight_exponent, weight_z_ref=args.weight_z_ref)
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
        print("Smooth-pivot: unsupported by forest/multipole engines; native runs use no-smooth-pivot.")
        print(FOREST_SELECTION_NOTE)
        for source in catalog.metadata.get("files", []):
            print(f"FITS field={source['delta_field']}, BLINDING={source['blinding']}")
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
