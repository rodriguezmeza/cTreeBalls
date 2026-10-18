#!/usr/bin/env python3
"""Run active cTreeBalls Ly-alpha forest engines on one retained catalog.

DESI, eBOSS/PICCA, Cartesian and RA-DEC-R FITS, lya-ascii, NPZ and synthetic input are supported.
MPI rank zero reads once, broadcasts arrays once, and all ranks register a
forest catalog with cyballs. Radial, anisotropic 3D, and multipole estimators
are distinct.
See tests/python/README_lya_corr_all_engines.md for commands and scientific limitations.
"""
from __future__ import annotations

import sys
from pathlib import Path

# Locate sibling analysis and benchmark helpers for direct script execution.
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python"))


import argparse
from run_guidance import add_run_guidance
from dataclasses import dataclass, field
import glob
from importlib.machinery import EXTENSION_SUFFIXES
import json
import math
import os
import shlex
import shutil
import subprocess
import time

import numpy as np
from catalog_memory import (add_memory_arguments, register_forest_catalog, preflight_transfer, memory_metadata)

from kappa_corr_all_engines import (
    broadcast_array, discover_cython_methods, get_mpi_comm,
    flatten_radial_matrix, mpi_environment_size, aggregate_rank_timings,
    _timing_metadata, timing_summary, write_timing_report,
)

PROJECT_ROOT = Path(__file__).resolve().parents[2]
# Prefer a matching local rebuild, including when this driver is copied into
# run2/python. An explicit source override remains authoritative.
_source_override = os.environ.get("CTREEBALLS_SOURCE_ROOT")
if _source_override:
    sys.path.insert(0, str(Path(_source_override).expanduser().resolve()))
else:
    for _candidate in (PROJECT_ROOT, *(p / "cTreeBalls" for p in Path(__file__).resolve().parents)):
        if ((_candidate / "Makefile").is_file()
                and any((_candidate / ("cyballs" + suffix)).is_file()
                        for suffix in EXTENSION_SUFFIXES)):
            sys.path.insert(0, str(_candidate))
            break
MPI_CHILD = "CTREEBALLS_LYA_MPI_CHILD"
DESI_BASE = ("https://data.desi.lbl.gov/public/dr1/vac/dr1/lya-deltas/"
             "v1.0/delta-lya-0-0/Delta/")
DESI_EXAMPLE = DESI_BASE + "delta-1019.fits.gz"
LYA_WAVELENGTH = 1215.67
ANISOTROPIC_MULTIPOLE = "lya-anisotropic-multipole-3pcf-omp"


@dataclass(frozen=True)
class EngineSpec:
    radial: bool
    orders: tuple[int, ...]
    mpi: bool
    tree: bool = False
    family: str = "lya"
    supports_smooth_pivot: bool = False
    explicit_only: bool = False


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
    ANISOTROPIC_MULTIPOLE: EngineSpec(False, (3,), False, tree=True,
                                    family="anisotropic_multipole", explicit_only=True),
    "lya-rprt-quadtree-2pcf-omp": EngineSpec(False, (2,), False, tree=True),
    "lya-los-tree-2pcf-omp": EngineSpec(False, (2,), False, tree=True),
    "lya-los-tree-3pcf-omp": EngineSpec(False, (3,), False, tree=True),
    "lya-los-tree-2pcf-3pcf-omp": EngineSpec(False, (2, 3), False, tree=True),
    "lya-los-tree-2pcf-mpi": EngineSpec(False, (2,), True, tree=True),
    "lya-los-tree-3pcf-mpi": EngineSpec(False, (3,), True, tree=True),
    "lya-los-tree-2pcf-3pcf-mpi": EngineSpec(False, (2, 3), True, tree=True),

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
    log_lambda: np.ndarray | None = None
    distortion_domain: dict | None = None

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
        ids = np.ascontiguousarray(ids, dtype=np.int64)
        lam = None if self.log_lambda is None else np.ascontiguousarray(self.log_lambda,dtype=np.float64)
        if lam is not None and (lam.shape != (len(pos),) or not np.isfinite(lam).all()):
            raise ValueError("log_lambda must be a finite (N,) array")
        from lya_distortion import normalize_domain
        domain = normalize_domain(self.distortion_domain,pos,ids,lam)
        return ForestCatalog(pos, delta, weights, ids, dict(self.metadata),lam,domain)

    @property
    def nbody(self):
        return len(self.positions)


@dataclass
class RunConfig:
    engines: tuple[str, ...]
    output_dir: Path
    qualification_reference: Path | None = None
    qualification_rtol: float = .05
    qualification_atol: float = 1e-10
    accuracy_report_only: bool = False
    threads: int = 1
    verbose: int = 1
    verbose_log: int = 0
    rp_max: float = 200.0
    rt_max: float = 200.0
    rp_bins: int = 50
    rt_bins: int = 50
    r3_max: float = 20.0
    r3_bins: int = 4
    theta_bins: int = 4
    mu_bins: int = 4
    lya2_kernel: int = 0
    lya2_rp_slop: float = 0.0
    lya2_rt_slop: float = 0.0
    lya3_kernel: int = 0
    lya3_pivot_cell_max: int = 8
    lya3_lmax: int = 8
    lya3_mu_mode: str = "exact"
    lya3_multipole_kernel: int = 0
    multipole_lmax: int = 3
    multipole_rmin: float = 0.0
    max_hist_mib: float = 1024.0
    plots: bool = True
    flatten_plots: bool = True
    plot_multipole_orders: tuple[int, ...] | None = None
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
    compute_distortion: bool = False
    distortion_projection_order: int = 1
    distortion_max_mib: float = 512.
    wedge_bins: int = 50
    wedge_mu_edges: tuple[float, ...] = (0., .5, .8, .95, 1.)
    wedge_subsamples: int = 10
    wedge_r_max: float | None = None
    transfer_catalog_ownership: bool = False

    def validate(self):
        if not self.engines or any(e not in LYA_ENGINES for e in self.engines):
            raise ValueError("select at least one registered Ly-alpha engine")
        for name in ("threads", "rp_bins", "rt_bins", "r3_bins", "theta_bins", "mu_bins"):
            value = getattr(self, name)
            if not isinstance(value, int) or not 1 <= value <= 100000:
                raise ValueError(f"{name} must be an integer in [1, 100000]")
        if self.verbose not in range(4) or self.verbose_log not in range(4):
            raise ValueError("verbose and verbose_log must be in [0,3]")
        if self.lya2_kernel not in (0,1) or self.lya3_kernel not in range(6):
            raise ValueError("require lya2_kernel=0..1 and lya3_kernel=0..5")
        for name in ("lya2_rp_slop", "lya2_rt_slop"):
            if not math.isfinite(getattr(self, name)) or not 0 <= getattr(self, name) <= 1:
                raise ValueError(f"{name} must be finite in [0,1]")
        if self.lya2_rp_slop or self.lya2_rt_slop:
            pair_engines = [e for e in self.engines if LYA_ENGINES[e].family == "lya"
                            and not LYA_ENGINES[e].radial and 2 in LYA_ENGINES[e].orders]
            if not pair_engines:
                raise ValueError("pair bin slop requires a distinct-forest 3D 2PCF engine")
            if self.lya2_kernel != 1 and any(e != "lya-rprt-quadtree-2pcf-omp" for e in pair_engines):
                raise ValueError("pair bin slop requires --lya2-kernel 1 (automatic for the rp-rt quadtree)")
            if self.compute_distortion or self.distortion_matrix:
                raise ValueError("nonzero pair bin slop cannot use an exact hard-bin distortion matrix; use zero pair slops for distortion modelling")
        if not isinstance(self.lya3_pivot_cell_max,int) or not 1<=self.lya3_pivot_cell_max<=64:
            raise ValueError("lya3_pivot_cell_max must be in [1,64]")
        if self.plot_multipole_orders is not None:
            from lya_multipole_plots import selected_orders
            selected_orders(self, range(33))
            if self.plots:
                for engine in self.engines:
                    spec = LYA_ENGINES[engine]
                    limit = self.lya3_lmax if engine == ANISOTROPIC_MULTIPOLE else self.multipole_lmax
                    if spec.family in ('anisotropic_multipole', 'multipole') and not selected_orders(self, range(limit+1)):
                        raise ValueError(f'{engine}: none of --plot-multipole-orders are available')
        if not isinstance(self.lya3_lmax, int) or not 0 <= self.lya3_lmax <= 32:
            raise ValueError("lya3_lmax must be an integer in [0,32]")
        if self.lya3_mu_mode not in ("exact", "finite-multipole"):
            raise ValueError("lya3_mu_mode must be exact or finite-multipole")
        if self.lya3_multipole_kernel not in (0, 1, 2):
            raise ValueError("lya3_multipole_kernel must be 0, 1 or 2")
        if ANISOTROPIC_MULTIPOLE in self.engines:
            hard_3pcf = any(LYA_ENGINES[e].family == "lya" and not LYA_ENGINES[e].radial
                            and 3 in LYA_ENGINES[e].orders for e in self.engines)
            if self.lya3_kernel and not hard_3pcf:
                raise ValueError("--lya3-kernel applies to hard-bin forest engines; for anisotropic "
                                 "multipoles use --lya3-multipole-kernel 0, 1 or 2")
            if self.qualification_reference:
                raise ValueError("anisotropic moment packets cannot use --qualification-reference; "
                                 "use the calibration script for raw moments, or select an exact "
                                 "hard-bin engine alongside this engine to compare its exact bins")
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
        if self.compute_distortion and self.distortion_matrix:
            raise ValueError("choose --compute-distortion or --distortion-matrix")
        if self.distortion_projection_order not in (0,1):
            raise ValueError("distortion_projection_order must be 0 or 1")
        if not math.isfinite(self.distortion_max_mib) or self.distortion_max_mib <= 0:
            raise ValueError("distortion_max_mib must be positive and finite")
        if self.compute_distortion and not has_3d_2pcf:
            raise ValueError("--compute-distortion requires a distinct-forest 3D 2PCF engine")
        if not self.compute_distortion and bool(self.distortion_matrix) != bool(self.model_correlation):
            raise ValueError("use --distortion-matrix together with --model-correlation (forward modelling)")
        if self.model_correlation and not self.analysis:
            raise ValueError("model comparison requires enabled analysis")
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
            elif spec.family == "anisotropic_multipole":
                bins3 = self.r3_bins**2 * self.theta_bins**2 * (
                    self.lya3_lmax + 1 + (self.mu_bins if self.lya3_mu_mode == "exact" else 0))
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
                          and not spec.explicit_only
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
    """Resolve files/globs and immediate FITS children of directories deterministically.

    Cached NPZ files and unrelated sidecars in a delta directory are ignored.
    Overlapping inputs are rejected to prevent silently counting a forest twice.
    """
    paths = []
    for pattern in patterns:
        matches = sorted(glob.glob(os.path.expanduser(str(pattern))))
        if not matches:
            raise FileNotFoundError(f"no files match {pattern}")
        for item in matches:
            path = Path(item).resolve()
            candidates = ([p for p in sorted(path.iterdir()) if p.is_file()
                           and p.name.lower().endswith((".fits", ".fits.gz", ".fit", ".fit.gz", ".fits.fz"))]
                          if path.is_dir() else [path])
            if not candidates:
                raise FileNotFoundError(f"no FITS files in {path}")
            for candidate in candidates:
                candidate = candidate.resolve()
                if not candidate.is_file():
                    raise FileNotFoundError(f"not an input file: {candidate}")
                if candidate in paths:
                    raise ValueError(f"duplicate input file: {candidate}")
                paths.append(candidate)
    return paths


def read_desi(paths, *, omega_m=0.315, h=0.674, z_min=0.0, z_max=10.0,
              max_forests=None, pixel_stride=1, delta_field="auto", project_delta=False,
              redshift_weight_exponent=0., weight_z_ref=2.25, retain_projection_domain=False):
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
    pieces, seen, provenance, domains, wavelengths = [], set(), [], [], []
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
            fit_wave = z_abs > 0
            distance_wave = fit_wave if retain_projection_domain and not project_delta else good_wave
            chi[distance_wave] = cosmology.comoving_distance(z_abs[distance_wave]).value * h
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
                wavelengths.append(np.log10(wave[indices]))
                if retain_projection_domain and not project_delta:
                    full = np.flatnonzero(fit_wave & np.isfinite(data) & np.isfinite(weight) & (weight>0))
                    domains.append(dict(positions=chi[full,None]*los,
                        forest_ids=np.full(len(full),identifier,dtype=np.int64),
                        log_lambda=np.log10(wave[full]),projection_weights=weight[full],
                        data_indices=np.searchsorted(full,indices)))
                loaded += 1
            provenance.append(dict(path=str(path), delta_field=field_name,
                                   blinding=str(meta.header.get("BLINDING", "unknown")),
                                   accepted_forests=loaded))
    if not pieces:
        raise ValueError("no valid forest pixels survived the selection")
    pieces = [list(piece) for piece in pieces]
    arrays = []
    for i in range(4):
        arrays.append(np.concatenate([piece[i] for piece in pieces]))
        for piece in pieces:
            piece[i] = None
    pieces.clear()
    from lya_distortion import concatenate_domains
    return ForestCatalog(*arrays, metadata=dict(
        input_format="DESI/PICCA images", files=provenance, distance_unit="Mpc/h",
        cosmology=dict(model="FlatLambdaCDM", omega_m=omega_m, h=h, Tcmb0=0),
        wavelength_lya_angstrom=LYA_WAVELENGTH, z_min=z_min, z_max=z_max,
        max_forests=max_forests, pixel_stride=pixel_stride,
        project_delta=project_delta, projection_domain=("retained pixels after cuts/stride" if project_delta else "input projection unchanged; original valid forest assumed for response"),
        redshift_weight_exponent=redshift_weight_exponent, weight_z_ref=weight_z_ref,
        weights=("input WEIGHT unchanged; zero/negative/nonfinite pixels excluded" if redshift_weight_exponent == 0
                 else "input WEIGHT times ((1+z)/(1+z_ref))**exponent"),
    ), log_lambda=np.concatenate(wavelengths),
        distortion_domain=concatenate_domains(domains)).normalized()


def read_ascii(path):
    dtype = [(n, "f8") for n in ("x", "y", "z", "delta", "weights")] + [("ids", "i8")]
    rows = np.loadtxt(path, dtype=dtype, ndmin=1)
    return ForestCatalog(np.column_stack([rows[n] for n in ("x", "y", "z")]),
                         rows["delta"], rows["weights"], rows["ids"],
                         dict(input_format="lya-ascii", path=str(path))).normalized()


def read_npz(path):
    with np.load(path, allow_pickle=False) as data:
        metadata = json.loads(str(data["metadata"])) if "metadata" in data else {}
        from lya_distortion import DOMAIN_FIELDS
        keys = ['distortion_'+k for k in DOMAIN_FIELDS]
        if any(k in data for k in keys+['distortion_log_lambda']) and not all(k in data for k in keys):
            raise ValueError("incomplete continuum-fitting domain in NPZ cache")
        domain = {k: data['distortion_'+k] for k in DOMAIN_FIELDS} if all(k in data for k in keys) else None
        if domain is not None and 'distortion_log_lambda' in data:
            domain['log_lambda'] = data['distortion_log_lambda']
        return ForestCatalog(*(data[n] for n in ("positions", "delta", "weights", "forest_ids")),
                             metadata=metadata, log_lambda=data['log_lambda'] if 'log_lambda' in data else None,
                             distortion_domain=domain).normalized()


def save_catalog(path, catalog):
    extra = {} if catalog.log_lambda is None else dict(log_lambda=catalog.log_lambda)
    if catalog.distortion_domain is not None:
        extra.update({'distortion_'+k: v for k,v in catalog.distortion_domain.items()})
    np.savez_compressed(path, positions=catalog.positions, delta=catalog.delta,
                        weights=catalog.weights, forest_ids=catalog.forest_ids,
                        metadata=json.dumps(catalog.metadata),**extra)


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
    # The distortion calculation runs once on rank zero; its fitting metadata
    # need not be duplicated on the other MPI ranks.
    lam = catalog.log_lambda if comm.rank == 0 else None
    domain = catalog.distortion_domain if comm.rank == 0 else None
    return collective(comm, lambda: ForestCatalog(*arrays, metadata=meta,
                       log_lambda=lam,distortion_domain=domain).normalized())


def engine_parameters(config, engine, root):
    spec = LYA_ENGINES[engine]
    if spec.family == "multipole":
        return dict(searchMethod=engine, rootDir=str(root), iCatalogs="1",
                    numberThreads=config.threads, usePeriodic=False,
                    useLogHist=False, rangeN=config.r3_max,
                    rminHist=config.multipole_rmin,
                    sizeHistN=config.r3_bins, mChebyshev=config.multipole_lmax,
                    verbose=config.verbose, verbose_log=config.verbose_log,
                    lya3MuMode=0, lya2Kernel=0, lya2RpSlop=0.0, lya2RtSlop=0.0,
                    options="only-3pcf-3d,exclude-all-same-los,no-smooth-pivot")
    # set() merges overrides: reset the order after a rank-zero-only
    # multipole run so later MPI startup controls agree on every rank.
    return dict(searchMethod=engine, rootDir=str(root), iCatalogs="1",
                numberThreads=config.threads, usePeriodic=False, useLogHist=False,
                rangeN=(config.r3_max if engine == ANISOTROPIC_MULTIPOLE
                        else max(config.rp_max, config.rt_max, config.r3_max)),
                rminHist=0.0, sizeHistN=4, mChebyshev=config.multipole_lmax,
                verbose=config.verbose, verbose_log=config.verbose_log,
                options="no-smooth-pivot",
                lya2RpMax=config.rp_max, lya2RtMax=config.rt_max,
                lya2RpBins=config.rp_bins, lya2RtBins=config.rt_bins,
                lya3RMax=config.r3_max, lya3RBins=config.r3_bins,
                lya3ThetaBins=config.theta_bins, lya3MuBins=config.mu_bins,
                lya2Kernel=1 if engine == "lya-rprt-quadtree-2pcf-omp" else (config.lya2_kernel if not spec.radial and 2 in spec.orders else 0),
                lya2RpSlop=config.lya2_rp_slop if spec.family == "lya" and not spec.radial and 2 in spec.orders else 0.0,
                lya2RtSlop=config.lya2_rt_slop if spec.family == "lya" and not spec.radial and 2 in spec.orders else 0.0,
                lya3Kernel=(config.lya3_multipole_kernel if engine == ANISOTROPIC_MULTIPOLE
                            else config.lya3_kernel if not spec.radial and 3 in spec.orders else 0),
                lya3LMax=config.lya3_lmax,
                lya3MuMode=int(engine == ANISOTROPIC_MULTIPOLE and config.lya3_mu_mode == "exact"),
                lya3PivotCellMax=config.lya3_pivot_cell_max)


def read_products(root, engine, config=None):
    spec, result = LYA_ENGINES[engine], {}
    if engine == ANISOTROPIC_MULTIPOLE:
        if config is None:
            raise ValueError("anisotropic products require the requested RunConfig, including mu mode")
        if config.lya3_mu_mode == "finite-multipole":
            # Signed raw moments and reconstructed bins have distinct contracts;
            # neither is inserted into the exact-ratio comparison pipeline.
            return result
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
    from lya_multipole_plots import write_scalar_plots
    for name, result in results.items():
        if "multipole_3pcf" in result["products"]:
            paths.extend(write_scalar_plots(config.output_dir/name,
                                            result["products"]["multipole_3pcf"], config))
    return paths


def replot_multipoles(source, output=None, orders=None):
    """Read saved histograms without loading a catalog or running an engine."""
    source = Path(source).expanduser().resolve()
    summary_path = source/'summary.json'
    saved = json.loads(summary_path.read_text())
    settings = {key: value for key, value in saved['config'].items()
                if key in RunConfig.__dataclass_fields__}
    settings.update(output_dir=Path(output).expanduser().resolve() if output else source/'multipole_plots',
                    plots=True, plot_multipole_orders=tuple(orders) if orders is not None else None)
    config = RunConfig(**settings)
    if (config.r3_bins < 1 or config.theta_bins < 1 or config.mu_bins < 1
            or not math.isfinite(config.r3_max) or config.r3_max <= 0
            or not 0 <= config.lya3_lmax <= 32 or not 0 <= config.multipole_lmax <= 10):
        raise ValueError('saved summary has an invalid multipole grid')
    from lya_multipole_plots import selected_orders, write_scalar_plots
    selected_orders(config, range(33))
    report = dict(source_summary=str(summary_path), native_search_run=False, plots=[], engines={})
    for engine in saved['engines']:
        if engine not in LYA_ENGINES:
            continue
        target = config.output_dir/engine
        if engine == ANISOTROPIC_MULTIPOLE:
            from lya_anisotropic_multipoles import inspect_outputs
            exact = read_products(source/engine, engine, config).get('3d_3pcf')
            item = inspect_outputs(source/engine, config, exact, plot_root=target)
            report['engines'][engine] = item
            report['plots'].extend(item['plots'])
        elif LYA_ENGINES[engine].family == 'multipole':
            table = read_products(source/engine, engine, config)['multipole_3pcf']
            paths = write_scalar_plots(target, table, config)
            report['engines'][engine] = dict(plots=paths, normalization='native scalar numerator / denominator')
            report['plots'].extend(paths)
    if not report['engines']:
        raise ValueError('the saved summary contains no supported 3PCF multipole engine')
    config.output_dir.mkdir(parents=True, exist_ok=True)
    path = config.output_dir/'plot_summary.json'
    path.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    return path


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
    distortion_report = None
    if config.compute_distortion:
        def compute_response():
            from lya_distortion import compute_catalog_distortion, save_distortion
            print("Computing continuum distortion matrix (rank 0, OpenMP)",flush=True)
            return save_distortion(config.output_dir/'distortion',compute_catalog_distortion(catalog,config))
        distortion_report = collective(comm,compute_response,root_only=True)
    balls = None
    from cyballs import cballs
    balls = collective(comm, cballs)
    results = {}
    try:
        if not hasattr(balls, "set_forest_catalog"):
            raise RuntimeError("rebuild cyballs: set_forest_catalog() is required")
        collective(comm, lambda: register_forest_catalog(balls, config.transfer_catalog_ownership,
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
                if (config.lya2_rp_slop or config.lya2_rt_slop) and spec.family == "lya" and not spec.radial and 2 in spec.orders:
                    print(f"  Approximate 2PCF bin placement: rp slop={config.lya2_rp_slop:g}, "
                          f"rt slop={config.lya2_rt_slop:g} bin widths; outer cuts and forest exclusions remain exact.", flush=True)
                if engine == ANISOTROPIC_MULTIPOLE:
                    print(f"  Lmax={config.lya3_lmax}, moment kernel={config.lya3_multipole_kernel}, "
                          f"mu output={config.lya3_mu_mode}", flush=True)
            started = time.perf_counter()
            def compute():
                nonlocal compute_wall, compute_cpu, native_cpu, native_timings
                if participates:
                    wall, cpu = time.perf_counter(), time.process_time()
                    balls.Run(level=["MainLoop"])
                    compute_wall, compute_cpu = time.perf_counter()-wall, time.process_time()-cpu
                    native_cpu = float(balls.getCPUTime())
                    native_timings = balls.getTimings()
                    if comm.rank == 0:
                        from cyballs import publish_result_packet
                        publish_result_packet(balls,root/'qualification',
                            Path(config.qualification_reference)/engine/'qualification' if config.qualification_reference else None,
                            rtol=config.qualification_rtol,atol=config.qualification_atol,
                    fail_on_rejection=not config.accuracy_report_only)
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
            product = collective(comm, lambda: read_products(root, engine, config), root_only=True)
            def anisotropic_outputs():
                if engine != ANISOTROPIC_MULTIPOLE:
                    return None
                from lya_anisotropic_multipoles import inspect_outputs
                return inspect_outputs(root, config, product.get("3d_3pcf"))
            extra = collective(comm, anisotropic_outputs, root_only=True)
            if comm.rank == 0:
                results[engine] = dict(wall_seconds=elapsed, products=product,
                                       ranks=comm.size if spec.mpi else 1,
                                       smooth_pivot="unsupported",
                                       provenance=balls.getRunMetadata(),
                                       threads_per_rank=config.threads,
                                       parameters=engine_parameters(config, engine, root))
                results[engine].update(timings)
                if extra is not None:
                    results[engine]["anisotropic_multipole"] = extra
                    print(f"  anisotropic output: {config.lya3_mu_mode}; raw moments through L={config.lya3_lmax}; "
                          "finite-L reconstruction remains approximate", flush=True)
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
                           memory_policy=memory_metadata(config.transfer_catalog_ownership),
                           timings=native_timings,
                           timing=("native wall_seconds includes MainLoop, output, synchronization and "
                                   "cleanup; FITS conversion is excluded. Reference wall_seconds is "
                                   "warmed CPU pair traversal/reduction; worker_wall_seconds includes "
                                   "imports, JIT, analysis and I/O. These are different scopes, not a speedup ratio."),
                           selection_contract=FOREST_SELECTION_NOTE,
                           engines={e: {k: v for k, v in r.items() if k != "products"}
                                    for e, r in results.items()}, comparisons={})
            if distortion_report is not None:
                summary['computed_distortion'] = distortion_report
            summary["comparisons"], difference_paths = write_comparisons(results, config)
            summary["comparison_contract"] = (
                "All pairs within each estimator family only; 3D, radial-only, same-LOS and "
                "multipole products are not interchangeable. No upstream 3PCF reference is run.")
            summary["plots"] = (make_plots(results, config) if config.plots else []) + difference_paths
            for row in results.values():
                if "anisotropic_multipole" in row:
                    summary["plots"].extend(row["anisotropic_multipole"]["plots"])
            if config.analysis:
                from dataclasses import replace
                analysis_config = replace(config,distortion_matrix=config.output_dir/'distortion'/'distortion.npz') if config.compute_distortion and config.model_correlation else config
                summary["analysis_2pcf"], analysis_paths = analyse_2pcf(results, analysis_config)
                summary["plots"].extend(analysis_paths)
            failed = [tag for tag, result in summary["comparisons"].items() if not result["passed"]]
            summary["validation"] = dict(passed=not failed, failed_comparisons=failed,
                                         comparisons=len(summary["comparisons"]))
            if ANISOTROPIC_MULTIPOLE in results:
                summary["validation"]["scope"] = (
                    "Exact-ratio product comparisons only; raw anisotropic moments and approximate "
                    "finite-L reconstruction are excluded. Zero comparisons is not an accuracy certificate.")
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
    add_memory_arguments(parser, healpix=False)
    parser.add_argument('--qualification-reference','--accuracy-reference',type=Path,help='Exact driver output root containing ENGINE/qualification packets for this same catalog')
    parser.add_argument('--qualification-rtol',type=float,default=.05, help='relative L2 target (default 0.05); does not bound every coefficient')
    parser.add_argument('--qualification-atol',type=float,default=1e-10)
    parser.add_argument('--accuracy-report-only',action='store_true',help='retain results and report penalties even above the relative L2 target')
    source = parser.add_mutually_exclusive_group()
    source.add_argument("--fits", nargs="+", help="DESI, eBOSS/PICCA, Cartesian, or RA-DEC-R forest FITS files/globs/directories")
    source.add_argument("--catalog", type=Path, help="NPZ with positions, delta, weights, forest_ids")
    source.add_argument("--ascii", type=Path, help="six columns: x y z delta weight forest_id")
    source.add_argument("--plot-existing", type=Path, help="plot saved multipole histograms from an output root; no catalog loading or search")
    source.add_argument("--synthetic", action="store_true", help="small generated catalog (default)")
    parser.add_argument(
        "--engine", "--engines", nargs="+", default=["lya-2pcf-omp"],
        help=("engine names or all/all-omp/all-mpi/all-1d/all-3d/all-tree/"
              "all-multipole"),
    )
    parser.add_argument("--statistics", choices=["2pcf", "3pcf", "both"], default="2pcf")
    parser.add_argument("--list-engines", action="store_true")
    parser.add_argument("--output", "--outdir", type=Path, help="output root (default Output_lya_all_engines; plot-existing: SOURCE/multipole_plots)")
    parser.add_argument("--save-catalog", type=Path, help="cache the selected input as NPZ")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--verbose", type=int, choices=range(4), default=1,
                        help="native verbosity; 1 enables live counting, 0 silences it")
    parser.add_argument("--verbose-log", type=int, choices=range(4), default=0,
                        help="native log verbosity; 1 records live counting in cballs.log")
    for name, default in (("rp-max", 200.), ("rt-max", 200.), ("r3-max", 20.)):
        parser.add_argument("--"+name, type=float, default=default, help="Mpc/h")
    for name, default in (("rp-bins", 50), ("rt-bins", 50), ("r3-bins", 4),
                          ("theta-bins", 4), ("mu-bins", 4)):
        parser.add_argument("--"+name, type=int, default=default)
    parser.add_argument("--lya2-kernel",type=int,choices=(0,1),default=0,help="3D pair kernel: 0 pixel, 1 certified forest cells/range moments")
    parser.add_argument("--lya2-bin-slop", "--bin-slop", type=float, default=None,
                        help="set both 3D 2PCF bin slops in bin-width units [0,1]; default 0 (exact); per-axis options override")
    parser.add_argument("--lya2-rp-slop", type=float, default=None,
                        help="parallel-bin leakage [0,1]; overrides --bin-slop for rp")
    parser.add_argument("--lya2-rt-slop", type=float, default=None,
                        help="transverse-bin leakage [0,1]; overrides --bin-slop for rt")
    parser.add_argument("--lya3-kernel",type=int,choices=range(6),default=0,help="3D triple kernel: 5 adaptive hierarchical radial moments; 0 retains segments")
    parser.add_argument("--lya3-pivot-cell-max",type=int,default=8,help="maximum pivot-cell population for kernels 4/5 (1..64)")
    parser.add_argument("--lya3-lmax", type=int, default=8, help="anisotropic Legendre moment order, 0..32")
    parser.add_argument("--lya3-mu-mode", choices=("exact", "finite-multipole"), default="exact",
                        help="anisotropic engine: exact adds exact angular bins; finite-multipole writes moments/reconstruction only")
    parser.add_argument("--lya3-multipole-kernel", type=int, choices=(0, 1, 2), default=0,
                        help="anisotropic moments: 0 automatic hierarchy, 1 prefix, 2 forced forest groups")
    parser.add_argument("--multipole-lmax", type=int, default=3)
    parser.add_argument("--multipole-rmin", type=float, default=0.0, help="Mpc/h")
    parser.add_argument("--omega-m", type=float, default=.315)
    parser.add_argument("--h", type=float, default=.674)
    parser.add_argument("--z-min", type=float, default=0.)
    parser.add_argument("--z-max", type=float, default=10.)
    parser.add_argument("--max-forests", type=int, help="FITS demonstration subset")
    parser.add_argument("--pixel-stride", type=int, default=1, help="FITS subsampling, not rebinning")
    parser.add_argument("--delta-field", choices=["auto", "DELTA", "DELTA_BLIND"], default="auto")
    parser.add_argument("--fits-layout", choices=["auto", "desi", "eboss", "cartesian", "radec-r"], default="auto")
    parser.add_argument("--eboss-angle-unit", choices=["rad", "deg"], default="rad")
    parser.add_argument("--radec-r-angle-unit", choices=["deg", "rad"], default="deg",
                        help="RA/DEC unit fallback for unlabelled RA-DEC-R columns (default deg)")
    parser.add_argument("--radec-r-distance-unit", choices=["Mpc/h", "Mpc"], default="Mpc/h",
                        help="rcdis unit fallback for unlabelled RA-DEC-R columns; Mpc is multiplied by h")
    parser.add_argument("--radec-r-distance-source", choices=["input", "redshift"], default="input",
                        help="input: rcdis must be an absolute observer distance; redshift: "
                             "compute chi(zobs) using --omega-m/--h in Mpc/h")
    parser.add_argument("--project-delta", action="store_true", help="remove each forest weighted mean and log-wavelength slope on retained FITS pixels or an NPZ with log_lambda")
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
    parser.add_argument("--compute-distortion",action="store_true",help="Compute native continuum 2PCF response once and save NPZ/FITS; optionally apply to --model-correlation")
    parser.add_argument("--distortion-projection-order",type=int,choices=[0,1],default=1,help="0=weighted mean; 1=weighted mean+slope in true log wavelength (default)")
    parser.add_argument("--distortion-max-mib",type=float,default=512.,help="Dense native matrix workspace budget; fewer bins/threads reduce memory")
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
    parser.add_argument("--plot-multipole-orders", type=int, nargs="+", choices=range(33),
                        help="orders for multipole maps and radial slices (default: all available); CSVs retain all orders")
    parser.add_argument("--no-plots", action="store_true")
    parser.add_argument(
        "--no-flatten-plots", action="store_true",
        help="omit flattened radial-bin views for 3PCF products",
    )
    parser.add_argument("--mpi-ranks", type=int, default=1)
    parser.add_argument("--mpiexec", help="launcher built against the same MPI as cyballs/mpi4py")
    parser.add_argument("--mpi-extra-arg", action="append", default=[])
    parser.epilog = 'Native timings separate Python Run from C MainLoop: MPI wall=max participating ranks, CPU=sum participating ranks. CPU includes worker threads and may exceed wall time. --verbose 1 enables live completion counters; --verbose-log 1 records them. LOS-tree methods are included in all-tree; estimator families are compared separately. Angular shear/scalar pivot-reuse and count-window options do not apply to forest estimators.'
    add_run_guidance(parser, field='lya', benchmark=False, private=False)
    args = parser.parse_args(argv)
    common = args.lya2_bin_slop
    if common is not None and (not math.isfinite(common) or not 0 <= common <= 1):
        parser.error("--bin-slop/--lya2-bin-slop must be finite in [0,1]")
    if args.lya2_rp_slop is None: args.lya2_rp_slop = 0.0 if common is None else common
    if args.lya2_rt_slop is None: args.lya2_rt_slop = 0.0 if common is None else common
    return args


def catalog_field_diagnostics(catalog):
    """Summarize the registered field with bounded temporary memory."""
    scale = float(np.max(catalog.weights))
    total_weight, weighted_delta, total_delta = 0., 0., 0.
    radius_min, radius_max = math.inf, 0.
    for start in range(0, catalog.nbody, 262144):
        s = slice(start, start+262144)
        weight = catalog.weights[s]/scale
        total_weight += float(np.sum(weight))
        weighted_delta += float(np.dot(weight, catalog.delta[s]))
        total_delta += float(np.sum(catalog.delta[s]))
        pos = catalog.positions[s]
        radius = np.hypot(np.hypot(pos[:, 0], pos[:, 1]), pos[:, 2])
        radius_min = min(radius_min, float(radius.min()))
        radius_max = max(radius_max, float(radius.max()))
    mean = weighted_delta/total_weight
    return dict(scope="retained pixels as registered with the native engine",
                weighted_mean_delta=mean, unweighted_mean_delta=total_delta/catalog.nbody,
                squared_weighted_mean_delta=mean**2,
                min_observer_distance_mpc_h=radius_min, max_observer_distance_mpc_h=radius_max)


def load_catalog(args):
    """Read and validate one catalog; shared by the forest benchmark frontend."""
    if not args.fits and (args.max_forests is not None or args.pixel_stride != 1
                          or args.z_min != 0 or args.z_max != 10
                          or args.omega_m != .315 or args.h != .674
                          or args.delta_field != "auto" or args.fits_layout != "auto"
                          or args.eboss_angle_unit != "rad"
                          or getattr(args, "radec_r_angle_unit", "deg") != "deg"
                          or getattr(args, "radec_r_distance_unit", "Mpc/h") != "Mpc/h"
                          or getattr(args, "radec_r_distance_source", "input") != "input"
                          or args.redshift_weight_exponent != 0 or args.weight_z_ref != 2.25):
        raise ValueError("FITS selection/cosmology options require --fits; "
                         "ASCII and NPZ coordinates are already comoving")
    if args.fits:
        from lya_fits import read_fits
        catalog = read_fits(args.fits, omega_m=args.omega_m, h=args.h,
                            z_min=args.z_min, z_max=args.z_max, max_forests=args.max_forests,
                            pixel_stride=args.pixel_stride, delta_field=args.delta_field,
                            fits_layout=args.fits_layout, eboss_angle_unit=args.eboss_angle_unit,
                            radec_r_angle_unit=getattr(args, "radec_r_angle_unit", "deg"),
                            radec_r_distance_unit=getattr(args, "radec_r_distance_unit", "Mpc/h"),
                            radec_r_distance_source=getattr(args, "radec_r_distance_source", "input"),
                            project_delta=args.project_delta,
                            redshift_weight_exponent=args.redshift_weight_exponent, weight_z_ref=args.weight_z_ref,
                            retain_projection_domain=getattr(args,'compute_distortion',False))
    elif args.catalog:
        catalog = read_npz(args.catalog)
    elif args.ascii:
        catalog = read_ascii(args.ascii)
    else:
        catalog = synthetic_catalog(args.synthetic_forests, args.synthetic_pixels, args.seed)
    if args.project_delta and not args.fits:
        if not args.catalog or catalog.log_lambda is None:
            raise ValueError('--project-delta requires FITS wavelengths or an NPZ with true log_lambda; '
                             'it cannot be inferred from Cartesian distance')
        from lya_projection import project_lya_delta
        catalog.delta, report = project_lya_delta(catalog.delta, catalog.weights,
            catalog.forest_ids, catalog.log_lambda, return_diagnostics=True)
        catalog.metadata.update(project_delta=True,
            projection_domain='retained cached pixels with correlation weights',
            continuum_projection=report)
        # Reprojection fits the retained input; a cached earlier full-forest
        # response domain no longer describes this operator.
        catalog.distortion_domain = None
    if args.project_delta:
        from lya_projection import forest_projection_diagnostics
        verified = forest_projection_diagnostics(catalog.delta, catalog.weights,
                                                 catalog.forest_ids, catalog.log_lambda)
        catalog.metadata['projection_constraints'] = verified
    # Record the actual retained field, including any selection, stride or
    # projection. An unweighted zero mean does not imply a weighted zero mean.
    catalog.metadata["field_diagnostics"] = catalog_field_diagnostics(catalog)
    if args.save_catalog:
        args.save_catalog.parent.mkdir(parents=True, exist_ok=True)
        save_catalog(args.save_catalog, catalog)
    return catalog


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    args = parse_arguments(argv)
    if args.plot_existing is not None:
        if args.no_plots or args.mpi_ranks != 1 or mpi_environment_size() > 1:
            raise ValueError('--plot-existing requires plotting enabled and a single Python process')
        path = replot_multipoles(args.plot_existing, args.output, args.plot_multipole_orders)
        print(f'Multipole plots regenerated without a native search: {path}')
        return 0
    if args.output is None:
        args.output = Path('Output_lya_all_engines')
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
                geometry = (("1D-tree" if spec.radial else "3D-tree") if spec.tree else
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
    config = RunConfig(engines, args.output, threads=args.threads,
            transfer_catalog_ownership=args.transfer_catalog_ownership,
                       verbose=args.verbose, verbose_log=args.verbose_log,
            qualification_reference=args.qualification_reference,
            accuracy_report_only=args.accuracy_report_only,
            qualification_rtol=args.qualification_rtol,
            qualification_atol=args.qualification_atol, rp_max=args.rp_max,
                       rt_max=args.rt_max, rp_bins=args.rp_bins, rt_bins=args.rt_bins,
                       r3_max=args.r3_max, r3_bins=args.r3_bins, theta_bins=args.theta_bins,
                       mu_bins=args.mu_bins, multipole_lmax=args.multipole_lmax,
                       lya2_kernel=args.lya2_kernel,lya3_kernel=args.lya3_kernel,
                       lya2_rp_slop=args.lya2_rp_slop, lya2_rt_slop=args.lya2_rt_slop,
                       lya3_pivot_cell_max=args.lya3_pivot_cell_max,
                       lya3_lmax=args.lya3_lmax, lya3_mu_mode=args.lya3_mu_mode,
                       lya3_multipole_kernel=args.lya3_multipole_kernel,
                       multipole_rmin=args.multipole_rmin,
                       max_hist_mib=args.max_hist_mib,
                       plot_multipole_orders=tuple(args.plot_multipole_orders) if args.plot_multipole_orders is not None else None,
                       plots=not args.no_plots,
                       flatten_plots=not args.no_flatten_plots,
                       lya2pcf_source=args.lya2pcf_source, reference_covariance=args.reference_covariance,
                       reference_nside=args.reference_nside, reference_timeout=args.reference_timeout,
                       rtol=args.rtol, atol=args.atol, relative_floor=args.relative_floor,
                       fail_on_mismatch=args.fail_on_mismatch, analysis=not args.no_analysis,
                       covariance=args.covariance, distortion_matrix=args.distortion_matrix,
                       model_correlation=args.model_correlation, wedge_bins=args.wedge_bins,
                       compute_distortion=args.compute_distortion,
                       distortion_projection_order=args.distortion_projection_order,
                       distortion_max_mib=args.distortion_max_mib,
                       wedge_mu_edges=tuple(float(x) for x in args.wedge_mu_edges.split(",")),
                       wedge_subsamples=args.wedge_subsamples, wedge_r_max=args.wedge_r_max)
    collective(comm, config.validate)
    collective(comm, lambda: preflight_transfer(
        config.transfer_catalog_ownership and any(e in LYA_ENGINES for e in config.engines), forest=True))
    catalog = collective(comm, lambda: load_catalog(args), root_only=True)
    catalog = broadcast_catalog(comm, catalog)
    if comm.rank == 0:
        print("Catalog loaded once; retained NumPy input will be reused by every engine.")
        print("Smooth-pivot: unsupported by forest/multipole engines; native runs use no-smooth-pivot.")
        print(FOREST_SELECTION_NOTE)
        if catalog.metadata.get('project_delta'):
            check = catalog.metadata.get('projection_constraints')
            print('Continuum projection: per-forest weighted mean and log-wavelength slope removed before all engines.')
            if check:
                print(f"  maximum weighted mean residual={check['max_abs_weighted_mean_after']:.3g}; "
                      f"normalized first-moment residual={check['max_abs_normalized_first_moment_after']:.3g}")
        for source in catalog.metadata.get("files", []):
            print(f"FITS field={source['delta_field']}, BLINDING={source['blinding']}")
        if "distance_source" in catalog.metadata:
            print(f"Observer distances: {catalog.metadata['distance_source']}")
        if "field_diagnostics" in catalog.metadata:
            diagnostic = catalog.metadata["field_diagnostics"]
            print(f"Retained delta: weighted mean={diagnostic['weighted_mean_delta']:.8g}; "
                  f"unweighted mean={diagnostic['unweighted_mean_delta']:.8g}")
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
