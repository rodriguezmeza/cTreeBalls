#!/usr/bin/env python3
"""Run one in-memory convergence catalog through cTreeBalls search engines.

The command-line interface reads a HEALPix FITS map, NPZ catalog, or generated
test map once.  It converts the selected pixels to Cartesian unit vectors and
registers those NumPy arrays once with ``cyballs.set_catalog()``.  Subsequent
engines rebuild only their C-owned trees and histograms; they do not read the
catalog from disk again.

MPI engines use the same contract.  Rank 0 reads the catalog once and mpi4py
broadcasts the arrays before every rank calls ``set_catalog()``.  Use
``--mpi-ranks`` to let this script relaunch itself with mpiexec, or invoke the
script under mpiexec directly.  mpi4py must use the same MPI implementation as
the one used to compile cTreeBalls.

The module is also importable.  A Python script or notebook can construct a
``KappaCatalog`` and call ``run_engine_suite`` without any FITS input.
``--edge-corrections`` selects corrected complex 3PCF multipoles. Engine
groups then include only methods supporting that product and the input mask.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass, field
import json
import math
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import time
from typing import Any, Dict, Iterable, Optional, Sequence

import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CBALLS = PROJECT_ROOT / "cballs"
MPI_CHILD_ENV = "CTREEBALLS_KAPPA_MPI_CHILD"


@dataclass(frozen=True)
class EngineSpec:
    mpi: bool = False
    has_2pcf: bool = True
    has_3pcf: bool = True
    force_log_bins: bool = False
    supports_mask: bool = False
    supports_edge: bool = False
    note: str = ""


# Scalar angular/convergence methods accepted by this driver.  Availability is
# still determined from the active executable, so optional profiles do not
# appear in ``all`` unless they were compiled.
KAPPA_ENGINES: Dict[str, EngineSpec] = {
    "octree-sincos-omp": EngineSpec(),
    "balls-omp": EngineSpec(),
    "kdtree-omp": EngineSpec(),
    "balltree-omp": EngineSpec(),
    "balltree-2balls-omp": EngineSpec(supports_edge=True),
    "balltree-2balls-mpi": EngineSpec(mpi=True, supports_edge=True),
    "octree-2balls-omp": EngineSpec(supports_mask=True, supports_edge=True),
    "octree-2balls-mpi": EngineSpec(mpi=True, supports_mask=True, supports_edge=True),
    "balltree-2balls-omp_3pcf": EngineSpec(has_2pcf=False, supports_edge=True),
    "balltree-2balls-mpi_3pcf": EngineSpec(mpi=True, has_2pcf=False, supports_edge=True),
    "balltree-mpi": EngineSpec(mpi=True),
    "octree-kkk-omp": EngineSpec(),
    "octree-ggg-omp": EngineSpec(supports_mask=True, supports_edge=True),
    "octree-ggg-mpi": EngineSpec(mpi=True, supports_mask=True, supports_edge=True),
    "octree-sincos-omp-addons": EngineSpec(),
    "tree-omp-sincos": EngineSpec(),
    "direct-sincos": EngineSpec(),
    "direct-simple-sincos": EngineSpec(),
    "octree-ggg-omp-triangles": EngineSpec(has_2pcf=False),
    "octree-balls4-omp": EngineSpec(force_log_bins=True, supports_mask=True, supports_edge=True),
    "octree-balls4-mpi": EngineSpec(mpi=True, force_log_bins=True, supports_mask=True, supports_edge=True),
    "octree-kkk-balls4-omp-triangles": EngineSpec(
        has_2pcf=False, force_log_bins=True
    ),
    "octree-ggg": EngineSpec(),
    "direct-simple-sincos-loopId": EngineSpec(),
    "balls-omp-0357": EngineSpec(),
}


INCOMPATIBLE_ENGINE_REASONS = {
    "octree-shear-omp": "requires gamma1/gamma2 spin-2 data, not a kappa map",
    "octree-ggg-cross-omp": "requires multiple distinct field catalogs",
    "kdtree-cute-box": "is a periodic Cartesian-box estimator",
    "kdtree-box-omp": "is a periodic Cartesian-box estimator",
    "octree-box-omp": "is a periodic Cartesian-box estimator",
    "neighbor-boxes-omp": "is a periodic Cartesian-box estimator",
    "octree-3pcf-3d-omp": "computes a physical 3D statistic, not angular kappa",
    "octree-ggg-3d-omp": "alias of the physical 3D scalar estimator",
    "octree-3pcf-3d-mpi": "computes a physical 3D statistic, not angular kappa",
    "octree-ggg-3d-mpi": "alias of the physical 3D scalar MPI estimator",
    "lya-2pcf-omp": "requires Ly-alpha forest IDs and radial-pixel semantics",
    "lya-3pcf-omp": "requires Ly-alpha forest IDs and radial-pixel semantics",
    "lya-2pcf-3pcf-omp": "requires Ly-alpha forest IDs and radial-pixel semantics",
    "lya-1d-2pcf-omp": "requires Ly-alpha forest IDs and radial-pixel semantics",
    "lya-1d-3pcf-omp": "requires Ly-alpha forest IDs and radial-pixel semantics",
    "lya-1d-2pcf-3pcf-omp": "requires Ly-alpha forest IDs and radial-pixel semantics",
    "lya-1d-tree-2pcf-omp": "requires Ly-alpha forest IDs and radial-pixel semantics",
    "lya-2pcf-mpi": "requires Ly-alpha forest IDs and radial-pixel semantics",
    "lya-3pcf-mpi": "requires Ly-alpha forest IDs and radial-pixel semantics",
    "lya-2pcf-3pcf-mpi": "requires Ly-alpha forest IDs and radial-pixel semantics",
    "lya-1d-2pcf-mpi": "requires Ly-alpha forest IDs and radial-pixel semantics",
    "lya-1d-3pcf-mpi": "requires Ly-alpha forest IDs and radial-pixel semantics",
    "lya-1d-2pcf-3pcf-mpi": "requires Ly-alpha forest IDs and radial-pixel semantics",
    "lya-1d-tree-2pcf-mpi": "requires Ly-alpha forest IDs and radial-pixel semantics",
}

for _forest_method in tuple(INCOMPATIBLE_ENGINE_REASONS):
    if _forest_method.startswith("lya-"):
        INCOMPATIBLE_ENGINE_REASONS[_forest_method] = (
            "requires forest IDs and observer distances; use "
            "python/lya_corr_all_engines.py and set_forest_catalog()"
        )

METHOD_ALIASES = {
    "octree-ggg-3d-omp": "octree-3pcf-3d-omp",
    "octree-ggg-3d-mpi": "octree-3pcf-3d-mpi",
}


@dataclass
class KappaCatalog:
    positions: np.ndarray
    kappa: np.ndarray
    weights: Optional[np.ndarray] = None
    mask: Optional[np.ndarray] = None
    metadata: Dict[str, Any] = field(default_factory=dict)

    def normalized(self) -> "KappaCatalog":
        positions = np.ascontiguousarray(self.positions, dtype=np.float64)
        kappa = np.ascontiguousarray(self.kappa, dtype=np.float64)
        if positions.ndim != 2 or positions.shape[1] != 3:
            raise ValueError(
                f"positions must have shape (N, 3), got {positions.shape}"
            )
        count = positions.shape[0]
        if count < 3:
            raise ValueError("a convergence catalog needs at least three points")
        if kappa.ndim != 1 or kappa.shape[0] != count:
            raise ValueError(f"kappa must have shape ({count},), got {kappa.shape}")
        if not np.all(np.isfinite(positions)) or not np.all(np.isfinite(kappa)):
            raise ValueError("positions and kappa must contain only finite values")

        weights = _optional_vector(self.weights, count, "weights", np.float64)
        if weights is not None:
            if not np.all(np.isfinite(weights)) or np.any(weights < 0.0):
                raise ValueError("weights must be finite and non-negative")
        mask = _optional_vector(self.mask, count, "mask", np.uint8)
        if mask is not None and not np.all((mask == 0) | (mask == 1)):
            raise ValueError("mask values must be boolean or 0/1")
        return KappaCatalog(
            positions=positions,
            kappa=kappa,
            weights=weights,
            mask=mask,
            metadata=dict(self.metadata),
        )

    @property
    def nbody(self) -> int:
        return int(self.positions.shape[0])


@dataclass
class RunConfig:
    engines: Sequence[str]
    output_dir: Path
    theta_min: float = 0.12250467913471644
    theta_max: float = 3.6279974066581295
    theta_scale: str = "degree"
    bins: int = 20
    multipoles: int = 7
    threads: int = max(1, (os.cpu_count() or 2) - 1)
    use_log_bins: bool = True
    tree_theta: float = 1.0
    nsmooth: int = 16
    phi_left: float = 0.0
    phi_right: float = 90.0
    theta_left: float = 0.0
    theta_right: float = 90.0
    options: Sequence[str] = field(default_factory=tuple)
    result_type: str = "sincos"
    verbose: int = 1
    verbose_log: int = 1
    continue_on_error: bool = False
    edge_corrections: bool = False

    @property
    def wants_edge_corrections(self) -> bool:
        return (
            self.edge_corrections
            or self.result_type == "edge_effects"
            or "edge-corrections" in _split_options(self.options)
        )

    def normalized(self) -> "RunConfig":
        output_dir = Path(self.output_dir).expanduser().resolve()
        engines = tuple(self.engines)
        if not engines:
            raise ValueError("at least one engine is required")
        if self.bins < 4:
            raise ValueError("out-m-HistZeta requires bins >= 4")
        if self.threads < 1:
            raise ValueError("threads must be positive")
        if self.multipoles < 2:
            raise ValueError("cTreeBalls requires multipoles >= 2")
        if self.nsmooth < 1:
            raise ValueError("nsmooth must be positive")
        if self.theta_scale not in {"degree", "radian", "au"}:
            raise ValueError("theta_scale must be degree, radian, or au")
        if not math.isfinite(self.theta_min) or not math.isfinite(self.theta_max):
            raise ValueError("theta limits must be finite")
        if self.theta_min <= 0.0 or self.theta_max <= self.theta_min:
            raise ValueError("require 0 < theta_min < theta_max")
        if self.result_type not in {"sincos", "edge_effects"}:
            raise ValueError("result_type must be sincos or edge_effects")
        options = tuple(_split_options(self.options))
        edge = self.wants_edge_corrections
        if edge and "only-2pcf" in options:
            raise ValueError("edge corrections require 3PCF; remove only-2pcf")
        return RunConfig(
            **{
                **self.__dict__,
                "engines": engines,
                "output_dir": output_dir,
                "options": options,
                "edge_corrections": edge,
                "result_type": "edge_effects" if edge else self.result_type,
            }
        )


class SerialComm:
    rank = 0
    size = 1

    def bcast(self, value: Any, root: int = 0) -> Any:
        return value

    def Bcast(self, value: np.ndarray, root: int = 0) -> None:
        return None

    def allgather(self, value: Any) -> list[Any]:
        return [value]

    def barrier(self) -> None:
        return None


def _optional_vector(
    values: Optional[np.ndarray], count: int, name: str, dtype: np.dtype
) -> Optional[np.ndarray]:
    if values is None:
        return None
    result = np.ascontiguousarray(values, dtype=dtype)
    if result.ndim != 1 or result.shape[0] != count:
        raise ValueError(f"{name} must have shape ({count},), got {result.shape}")
    return result


def _split_options(values: Iterable[str] | str | None) -> list[str]:
    if values is None:
        return []
    if isinstance(values, str):
        values = [values]
    result: list[str] = []
    for value in values:
        for item in str(value).split(","):
            item = item.strip()
            if item and item.lower() != "none" and item not in result:
                result.append(item)
    return result


def discover_search_methods(executable: Path = DEFAULT_CBALLS) -> list[str]:
    executable = Path(executable).expanduser().resolve()
    if not executable.is_file():
        raise FileNotFoundError(f"cballs executable not found: {executable}")
    completed = subprocess.run(
        [os.fspath(executable), "options=print-search-methods"],
        cwd=executable.parent,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    methods = re.findall(r"^- ([^ ]+) \(id=-?\d+\)$", completed.stdout, re.MULTILINE)
    if not methods:
        raise RuntimeError(
            f"could not discover search methods from {executable}:\n"
            f"{completed.stdout[-2000:]}"
        )
    return methods


def discover_cython_methods(candidates: Sequence[str]) -> list[str]:
    try:
        from cyballs import search_method_id
    except ImportError as exc:
        raise RuntimeError(
            "the installed cyballs extension lacks search_method_id(); "
            "rebuild it with `make cyballs` from this source tree"
        ) from exc
    return [name for name in candidates if search_method_id(name) >= 0]


def print_engine_table(available: Sequence[str]) -> None:
    available_set = set(available)
    print("Kappa-compatible engines:")
    for name, spec in KAPPA_ENGINES.items():
        state = "available" if name in available_set else "not built"
        parallel = "MPI+OpenMP" if spec.mpi else "OpenMP/serial"
        products = "3PCF" if not spec.has_2pcf else (
            "2PCF" if not spec.has_3pcf else "2PCF+3PCF"
        )
        mask = "mask=yes" if spec.supports_mask else "mask=no"
        edge = "edge=yes" if spec.supports_edge else "edge=no"
        print(f"  {name:38s} {state:11s} {parallel:12s} {products:9s} {mask} {edge}")
    incompatible = [name for name in available if name in INCOMPATIBLE_ENGINE_REASONS]
    if incompatible:
        print("\nBuilt engines intentionally rejected by a kappa-map driver:")
        for name in incompatible:
            print(f"  {name:38s} {INCOMPATIBLE_ENGINE_REASONS[name]}")


def resolve_engines(
    tokens: Sequence[str], available: Sequence[str], *,
    edge_corrections: bool = False, masked: bool = False,
) -> list[str]:
    requested: list[str] = []
    for token in tokens or ("octree-ggg-omp",):
        requested.extend(item.strip() for item in token.split(",") if item.strip())
    available_set = set(available)
    selected: list[str] = []
    for name in requested:
        if name in {"all", "all-kappa"}:
            additions = [item for item in KAPPA_ENGINES if item in available_set]
        elif name == "all-omp":
            additions = [
                item for item, spec in KAPPA_ENGINES.items()
                if not spec.mpi and item in available_set
            ]
        elif name == "all-mpi":
            additions = [
                item for item, spec in KAPPA_ENGINES.items()
                if spec.mpi and item in available_set
            ]
        else:
            if name in INCOMPATIBLE_ENGINE_REASONS:
                raise ValueError(f"{name}: {INCOMPATIBLE_ENGINE_REASONS[name]}")
            if name not in KAPPA_ENGINES:
                raise ValueError(f"{name} is not registered as a kappa engine")
            if name not in available_set:
                raise ValueError(f"{name} is not present in the active cballs build")
            if edge_corrections and not KAPPA_ENGINES[name].supports_edge:
                raise ValueError(f"{name} does not support angular edge corrections")
            if masked and not KAPPA_ENGINES[name].supports_mask:
                raise ValueError(f"{name} does not support read-mask in this driver")
            additions = [name]
        if edge_corrections:
            additions = [item for item in additions if KAPPA_ENGINES[item].supports_edge]
        if masked:
            additions = [item for item in additions if KAPPA_ENGINES[item].supports_mask]
        for item in additions:
            if item not in selected:
                selected.append(item)
    if not selected:
        requirements = []
        if edge_corrections:
            requirements.append("edge corrections")
        if masked:
            requirements.append("read-mask")
        suffix = " supporting " + " and ".join(requirements) if requirements else ""
        raise ValueError("the selected build contains no requested kappa engines" + suffix)
    return selected


def angular_limits(config: RunConfig) -> tuple[float, float]:
    if config.theta_scale == "degree":
        return math.radians(config.theta_min), math.radians(config.theta_max)
    return config.theta_min, config.theta_max


def _healpix_valid(values: np.ndarray, hp: Any) -> np.ndarray:
    return np.isfinite(values) & ~hp.mask_bad(values)


def catalog_from_healpix(
    fits_path: Path,
    field_index: int = 0,
    nside_down: int = 0,
    mask_path: Optional[Path] = None,
    mask_threshold: float = 0.0,
    center_field: bool = True,
) -> KappaCatalog:
    try:
        import healpy as hp
    except ImportError as exc:
        raise RuntimeError("HEALPix FITS input requires healpy") from exc

    fits_path = Path(fits_path).expanduser().resolve()
    values = np.asarray(
        hp.read_map(fits_path, field=field_index, dtype=np.float64),
        dtype=np.float64,
    )
    nside_in = int(hp.get_nside(values))
    if nside_down:
        if not hp.isnsideok(nside_down):
            raise ValueError("nside_down must be a valid HEALPix NSIDE")
        if nside_down > nside_in:
            raise ValueError("nside_down cannot exceed the input NSIDE")
        if nside_down != nside_in:
            values = hp.ud_grade(
                values, nside_out=nside_down,
                order_in="RING", order_out="RING", power=0.0,
            )
    nside = int(hp.get_nside(values))
    valid = _healpix_valid(values, hp)

    mask_values = None
    if mask_path is not None:
        mask_path = Path(mask_path).expanduser().resolve()
        mask_values = np.asarray(
            hp.read_map(mask_path, field=0, dtype=np.float64), dtype=np.float64
        )
        mask_nside = int(hp.get_nside(mask_values))
        if mask_nside != nside:
            mask_values = hp.ud_grade(
                mask_values, nside_out=nside,
                order_in="RING", order_out="RING", power=0.0,
            )
        mask_valid = _healpix_valid(mask_values, hp)
        mask_full = mask_valid & (mask_values > mask_threshold)
    else:
        mask_full = np.ones(values.shape, dtype=bool)

    pixels = np.flatnonzero(valid)
    if pixels.size < 3:
        raise ValueError("the HEALPix map contains fewer than three valid pixels")
    x, y, z = hp.pix2vec(nside, pixels, nest=False)
    positions = np.ascontiguousarray(np.column_stack((x, y, z)), dtype=np.float64)
    kappa = np.ascontiguousarray(values[pixels], dtype=np.float64)
    selected_mask = np.ascontiguousarray(mask_full[pixels], dtype=np.uint8)
    if center_field:
        mean_selection = selected_mask.astype(bool)
        if not np.any(mean_selection):
            raise ValueError("the mask excludes every valid convergence pixel")
        kappa = kappa.copy()
        kappa -= np.mean(kappa[mean_selection], dtype=np.float64)
    return KappaCatalog(
        positions=positions,
        kappa=kappa,
        mask=selected_mask if mask_values is not None else None,
        metadata={
            "source": os.fspath(fits_path),
            "mask_source": os.fspath(mask_path) if mask_path is not None else None,
            "nside_input": nside_in,
            "nside": nside,
            "field": field_index,
            "centered": center_field,
        },
    ).normalized()


def synthetic_healpix_catalog(nside: int = 4, center_field: bool = True) -> KappaCatalog:
    try:
        import healpy as hp
    except ImportError as exc:
        raise RuntimeError("synthetic HEALPix input requires healpy") from exc
    if not hp.isnsideok(nside):
        raise ValueError("synthetic_nside must be a valid HEALPix NSIDE")
    pixels = np.arange(hp.nside2npix(nside))
    x, y, z = hp.pix2vec(nside, pixels, nest=False)
    positions = np.column_stack((x, y, z))
    kappa = 0.35 * x - 0.21 * y + 0.08 * z * z
    if center_field:
        kappa -= np.mean(kappa)
    weights = 0.9 + 0.1 * (1.0 + z)
    return KappaCatalog(
        positions=positions,
        kappa=kappa,
        weights=weights,
        metadata={"source": "synthetic", "nside": nside, "centered": center_field},
    ).normalized()


def catalog_from_npz(path: Path, center_field: bool = False) -> KappaCatalog:
    path = Path(path).expanduser().resolve()
    with np.load(path, allow_pickle=False) as archive:
        catalog = KappaCatalog(
            positions=archive["positions"],
            kappa=archive["kappa"],
            weights=archive["weights"] if "weights" in archive else None,
            mask=archive["mask"] if "mask" in archive else None,
            metadata={"source": os.fspath(path)},
        ).normalized()
    if center_field:
        selection = (
            catalog.mask.astype(bool)
            if catalog.mask is not None
            else np.ones(catalog.nbody, dtype=bool)
        )
        catalog.kappa = catalog.kappa.copy()
        catalog.kappa -= np.mean(catalog.kappa[selection], dtype=np.float64)
        catalog.metadata["centered"] = True
    return catalog


def save_catalog_npz(path: Path, catalog: KappaCatalog) -> None:
    payload: Dict[str, np.ndarray] = {
        "positions": catalog.positions,
        "kappa": catalog.kappa,
    }
    if catalog.weights is not None:
        payload["weights"] = catalog.weights
    if catalog.mask is not None:
        payload["mask"] = catalog.mask
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, **payload)


def get_mpi_comm(required: bool) -> Any:
    if not required:
        return SerialComm()
    try:
        from mpi4py import MPI
    except ImportError as exc:
        raise RuntimeError(
            "MPI engines (including one-rank suites) require mpi4py; install mpi4py "
            "against the same MPI used to build cTreeBalls"
        ) from exc
    return MPI.COMM_WORLD


def mpi_environment_size() -> int:
    """Best-effort detection before mpi4py initializes the Python MPI layer."""
    for name in (
        "OMPI_COMM_WORLD_SIZE", "PMI_SIZE", "PMIX_SIZE",
        "MV2_COMM_WORLD_SIZE", "MPI_LOCALNRANKS",
    ):
        value = os.environ.get(name)
        if value:
            try:
                return int(value)
            except ValueError:
                pass
    return 1


def broadcast_array(comm: Any, array: Optional[np.ndarray]) -> Optional[np.ndarray]:
    present = comm.bcast(array is not None if comm.rank == 0 else None, root=0)
    if not present:
        return None
    descriptor = comm.bcast(
        (array.shape, array.dtype.str) if comm.rank == 0 else None, root=0
    )
    if comm.rank != 0:
        array = np.empty(descriptor[0], dtype=np.dtype(descriptor[1]))
    comm.Bcast(array, root=0)
    return array


def broadcast_catalog(comm: Any, catalog: Optional[KappaCatalog]) -> KappaCatalog:
    if comm.size == 1:
        if catalog is None:
            raise ValueError("rank 0 did not provide a catalog")
        return catalog.normalized()
    metadata = comm.bcast(catalog.metadata if comm.rank == 0 else None, root=0)
    result = KappaCatalog(
        positions=broadcast_array(comm, catalog.positions if comm.rank == 0 else None),
        kappa=broadcast_array(comm, catalog.kappa if comm.rank == 0 else None),
        weights=broadcast_array(
            comm, catalog.weights if comm.rank == 0 and catalog is not None else None
        ),
        mask=broadcast_array(
            comm, catalog.mask if comm.rank == 0 and catalog is not None else None
        ),
        metadata=metadata,
    )
    return result.normalized()


def engine_parameters(config: RunConfig, engine: str, masked: bool) -> dict:
    spec = KAPPA_ENGINES[engine]
    rmin, rmax = angular_limits(config)
    angular_bins_needed = max(4, 2 * config.multipoles + 1)
    angular_bins = 1 << (angular_bins_needed - 1).bit_length()
    options = ["compute-HistN", "and-CF", "out-m-HistZeta", "KKKCorrelation"]
    options.extend(config.options)
    if masked and "read-mask" not in options:
        options.append("read-mask")
    options = _split_options(options)
    if engine in {"octree-balls4-omp", "octree-balls4-mpi"} and "smooth-pivot" in options:
        raise ValueError(f"{engine} does not support options=smooth-pivot")
    if config.wants_edge_corrections:
        if not spec.supports_edge:
            raise ValueError(f"{engine} does not support angular edge corrections")
        if "only-2pcf" in options:
            raise ValueError("edge corrections require 3PCF; remove only-2pcf")
        options = _split_options(
            [*options, "edge-corrections", "no-normalize-HistZeta"]
        )
    return {
        "searchMethod": engine,
        "rangeN": rmax,
        "rminHist": rmin,
        "sizeHistN": config.bins,
        "mChebyshev": config.multipoles,
        "sizeHistPhi": angular_bins,
        "numberThreads": config.threads,
        "useLogHist": bool(config.use_log_bins or spec.force_log_bins),
        "usePeriodic": False,
        "lengthBox": 2.0,
        "theta": config.tree_theta,
        "nsmooth": config.nsmooth,
        "phiL": math.radians(config.phi_left),
        "phiR": math.radians(config.phi_right),
        "thetaL": math.radians(config.theta_left),
        "thetaR": math.radians(config.theta_right),
        "iCatalogs": "1",
        "rootDir": os.fspath(config.output_dir / engine),
        "options": ",".join(options),
        "verbose": config.verbose,
        "verbose_log": config.verbose_log,
    }


def copy_engine_results(balls: Any, engine: str, config: RunConfig) -> dict:
    spec = KAPPA_ENGINES[engine]
    options = set(_split_options(config.options))
    edge = config.wants_edge_corrections
    want_2pcf = spec.has_2pcf and "only-3pcf" not in options
    want_3pcf = spec.has_3pcf and "only-2pcf" not in options
    result: Dict[str, Any] = {
        "engine": engine,
        "cpu_time": float(balls.getCPUTime()),
        "nbody": int(balls.getNBody()),
        "result_type": "edge_effects" if edge else config.result_type,
        "edge_corrections": edge,
        "warnings": [],
    }
    try:
        result["r"] = np.asarray(balls.getrBins()).copy()
    except Exception as exc:
        result["warnings"].append(f"r bins unavailable: {exc}")
    if want_2pcf:
        for key, getter in (
            ("xi", balls.getHistXi2pcf),
            ("nn", balls.getHistNN),
        ):
            try:
                result[key] = np.asarray(getter()).copy()
            except Exception as exc:
                result["warnings"].append(f"{key} unavailable: {exc}")
    if want_3pcf:
        try:
            multipoles = int(balls.getnMultipoles())
            result["multipoles"] = multipoles
            if not edge:
                for order in range(1, multipoles + 2):
                    components = []
                    for component, name in enumerate(
                        ("cos", "sin", "sincos", "cossin"), start=1
                    ):
                        value = np.asarray(
                            balls.getHistZetaMsincos(order, component)
                        ).copy()
                        result[f"zeta_{name}_{order}"] = value
                        components.append(value)
                    result[f"zeta_cos_plus_sin_{order}"] = (
                        components[0] + components[1]
                    )
            else:
                for order in range(1, multipoles + 2):
                    result[f"zeta_edge_{order}"] = np.asarray(
                        balls.getHistZetaM_EE(order)
                    ).copy()
                    result[f"zeta_edge_im_{order}"] = np.asarray(
                        balls.getHistZetaM_EE_Im(order)
                    ).copy()
                    for key in (f"zeta_edge_{order}", f"zeta_edge_im_{order}"):
                        value = result[key]
                        if value.shape != (config.bins, config.bins):
                            raise ValueError(f"{key} has unexpected shape {value.shape}")
                        if not np.all(np.isfinite(value)):
                            raise ValueError(f"{key} contains nonfinite values")
                    result[f"zeta_edge_complex_{order}"] = (
                        result[f"zeta_edge_{order}"]
                        + 1j * result[f"zeta_edge_im_{order}"]
                    )
        except Exception as exc:
            if edge:
                raise RuntimeError(
                    f"{engine}: requested edge-corrected 3PCF is unavailable: {exc}"
                ) from exc
            result["warnings"].append(f"3PCF multipoles unavailable: {exc}")
    return result


def save_engine_results(root: Path, result: dict) -> None:
    python_root = root / "python"
    python_root.mkdir(parents=True, exist_ok=True)
    arrays = {key: value for key, value in result.items() if isinstance(value, np.ndarray)}
    np.savez_compressed(python_root / "histograms.npz", **arrays)
    if "r" in result and "xi" in result:
        np.savetxt(
            python_root / "histXi2pcf.txt",
            np.column_stack((result["r"], result["xi"])),
            header="r xi",
        )
    for key, value in arrays.items():
        match = re.fullmatch(r"zeta_(?:cos_plus_sin|edge)_(\d+)", key)
        if match:
            np.savetxt(
                python_root / f"histZetaM_EE_{match.group(1)}.txt", value
            )
        match_im = re.fullmatch(r"zeta_edge_im_(\d+)", key)
        if match_im:
            np.savetxt(
                python_root / f"histZetaM_EE_Im_{match_im.group(1)}.txt", value
            )
    metadata = {
        key: value
        for key, value in result.items()
        if not isinstance(value, np.ndarray)
    }
    metadata["array_shapes"] = {key: list(value.shape) for key, value in arrays.items()}
    (python_root / "result.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )


def run_engine_suite(
    catalog: KappaCatalog,
    config: RunConfig,
    comm: Any = None,
) -> dict[str, dict]:
    """Run all selected engines while retaining one registered NumPy catalog."""
    from cyballs import cballs

    config = config.normalized()
    if comm is None:
        comm = get_mpi_comm(
            mpi_environment_size() > 1
            or any(KAPPA_ENGINES[e].mpi for e in config.engines if e in KAPPA_ENGINES)
        )
    catalog = catalog.normalized()
    masked = catalog.mask is not None
    for engine in config.engines:
        if engine not in KAPPA_ENGINES:
            raise ValueError(f"unknown kappa engine: {engine}")
        # Validate requested products before constructing an object on any rank.
        engine_parameters(config, engine, masked)
        if masked and not KAPPA_ENGINES[engine].supports_mask:
            supported = ", ".join(
                name for name, spec in KAPPA_ENGINES.items() if spec.supports_mask
            )
            raise ValueError(
                f"{engine} does not document read-mask semantics; use "
                f"one of {supported} for a masked kappa map"
            )

    config.output_dir.mkdir(parents=True, exist_ok=True)
    balls = cballs()
    balls.set_catalog(
        catalog.positions,
        kappa=catalog.kappa,
        weights=catalog.weights,
        mask=catalog.mask,
    )
    if balls.catalog_count != 1:
        raise RuntimeError("cyballs did not retain exactly one in-memory catalog")

    results: dict[str, dict] = {}
    failures: dict[str, str] = {}
    try:
        for engine in config.engines:
            spec = KAPPA_ENGINES[engine]
            comm.barrier()
            participates = spec.mpi or comm.rank == 0
            local_error = None
            result = None
            if participates:
                engine_root = config.output_dir / engine
                engine_root.mkdir(parents=True, exist_ok=True)
                balls.set(engine_parameters(config, engine, masked))
                if comm.rank == 0:
                    print(
                        f"Running {engine} with {config.threads} OpenMP thread(s)"
                        + (f" on {comm.size} MPI rank(s)" if spec.mpi else ""),
                        flush=True,
                    )
                started = time.perf_counter()
                try:
                    balls.Run(level=["MainLoop"])
                    if comm.rank == 0:
                        result = copy_engine_results(balls, engine, config)
                        result["wall_time"] = time.perf_counter() - started
                except Exception as exc:
                    local_error = f"{type(exc).__name__}: {exc}"
                finally:
                    try:
                        balls.struct_cleanup()
                    except Exception as exc:
                        cleanup_error = f"cleanup {type(exc).__name__}: {exc}"
                        local_error = (
                            f"{local_error}; {cleanup_error}"
                            if local_error else cleanup_error
                        )

            if spec.mpi:
                errors = comm.allgather(local_error)
            else:
                root_error = comm.bcast(local_error if comm.rank == 0 else None, root=0)
                errors = [root_error]
            errors = [error for error in errors if error]
            if errors:
                message = "; ".join(dict.fromkeys(errors))
                failures[engine] = message
                if comm.rank == 0:
                    print(f"FAILED {engine}: {message}", file=sys.stderr)
                if not config.continue_on_error:
                    raise RuntimeError(f"{engine} failed: {message}")
            elif comm.rank == 0 and result is not None:
                results[engine] = result
                save_engine_results(config.output_dir / engine, result)
            comm.barrier()
    finally:
        balls.clear_catalogs()

    if comm.rank == 0:
        summary = {
            "catalog": {**catalog.metadata, "nbody": catalog.nbody},
            "engines": {
                engine: {
                    "cpu_time": value.get("cpu_time"),
                    "wall_time": value.get("wall_time"),
                    "edge_corrections": value["edge_corrections"],
                    "result_type": value["result_type"],
                    "warnings": value.get("warnings", []),
                }
                for engine, value in results.items()
            },
            "failures": failures,
            "catalog_reads": 1,
            "set_catalog_calls_per_process": 1,
            "mpi_ranks": comm.size,
            "threads_per_rank": config.threads,
        }
        (config.output_dir / "summary.json").write_text(
            json.dumps(summary, indent=2) + "\n", encoding="utf-8"
        )
    return results


def spawn_mpi(args: argparse.Namespace) -> int:
    executable = args.mpiexec or shutil.which("mpiexec") or shutil.which("mpirun")
    if not executable:
        raise RuntimeError("mpiexec was not found; pass --mpiexec /path/to/mpiexec")
    command = [os.fspath(executable)]
    for value in args.mpi_extra_arg:
        command.extend(shlex.split(value))
    command.extend(
        ["-n", str(args.mpi_ranks), sys.executable, os.fspath(Path(__file__).resolve())]
    )
    command.extend(sys.argv[1:])
    environment = os.environ.copy()
    environment[MPI_CHILD_ENV] = "1"
    completed = subprocess.run(command, env=environment, check=False)
    return completed.returncode


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    sources = parser.add_mutually_exclusive_group()
    sources.add_argument("--fits", type=Path, help="HEALPix convergence FITS map")
    sources.add_argument("--catalog-npz", type=Path, help="positions/kappa NPZ catalog")
    sources.add_argument("--synthetic-nside", type=int, help="generated smoke-test map")
    parser.add_argument("--field", type=int, default=0)
    parser.add_argument("--mask", type=Path, help="optional HEALPix 0/1 mask")
    parser.add_argument("--mask-threshold", type=float, default=0.0)
    parser.add_argument("--nside-down", type=int, default=0)
    parser.add_argument("--no-center-field", action="store_true")
    parser.add_argument("--save-catalog-npz", type=Path)
    parser.add_argument(
        "--engine", "--engines", dest="engines", action="append", default=[],
        help="repeat or use comma lists; all, all-omp, and all-mpi select available "
             "methods compatible with the requested edge correction and mask",
    )
    parser.add_argument("--list-engines", action="store_true")
    parser.add_argument("--cballs", type=Path, default=DEFAULT_CBALLS)
    parser.add_argument("--outdir", type=Path, default=Path("Output_all_engines"))
    parser.add_argument("--threads", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    parser.add_argument("--mpi-ranks", type=int, default=1)
    parser.add_argument("--mpiexec", type=Path)
    parser.add_argument("--mpi-extra-arg", action="append", default=[])
    parser.add_argument("--theta-min", type=float, default=0.12250467913471644)
    parser.add_argument("--theta-max", type=float, default=3.6279974066581295)
    parser.add_argument(
        "--theta-scale", choices=("degree", "radian", "au"), default="degree"
    )
    parser.add_argument("--nbins", type=int, default=20)
    parser.add_argument("--multipoles", type=int, default=7)
    parser.add_argument("--tree-theta", type=float, default=1.0)
    parser.add_argument("--nsmooth", type=int, default=16)
    parser.add_argument("--linear-bins", action="store_true")
    parser.add_argument("--thetaL", type=float, default=0.0)
    parser.add_argument("--thetaR", type=float, default=90.0)
    parser.add_argument("--phiL", type=float, default=0.0)
    parser.add_argument("--phiR", type=float, default=90.0)
    parser.add_argument("--more-options", action="append", default=[])
    parser.add_argument(
        "--edge-corrections", action="store_true",
        help="compute and save complex edge-corrected 3PCF; equivalent to --type edge_effects",
    )
    parser.add_argument(
        "--type", dest="result_type", choices=("sincos", "edge_effects"),
        default="sincos",
        help="3PCF output: sincos components or edge-corrected complex multipoles",
    )
    parser.add_argument("--verbose", type=int, default=1)
    parser.add_argument("--verbose-log", type=int, default=1)
    parser.add_argument("--continue-on-error", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_arguments()
    if args.mpi_ranks < 1:
        raise SystemExit("ERROR: --mpi-ranks must be positive")

    try:
        executable_methods = discover_search_methods(args.cballs)
        cython_candidates = list(
            dict.fromkeys(
                executable_methods
                + list(KAPPA_ENGINES)
                + list(INCOMPATIBLE_ENGINE_REASONS)
            )
        )
        cython_methods = discover_cython_methods(cython_candidates)
        available = [name for name in executable_methods if name in cython_methods]
        executable_only = [
            name for name in executable_methods if name not in cython_methods
        ]
        cython_only = [
            name for name in cython_methods
            if name not in executable_methods and name not in METHOD_ALIASES
        ]
        if executable_only or cython_only:
            print(
                "WARNING: cballs and cyballs were built with different search "
                "profiles; only their intersection will be used.",
                file=sys.stderr,
            )
            if executable_only:
                print(
                    "  executable only: " + ", ".join(executable_only),
                    file=sys.stderr,
                )
            if cython_only:
                print(
                    "  cyballs only: " + ", ".join(cython_only),
                    file=sys.stderr,
                )
        if args.list_engines:
            print_engine_table(available)
            return 0
        config = RunConfig(
            engines=args.engines or ("octree-ggg-omp",),
            output_dir=args.outdir,
            theta_min=args.theta_min,
            theta_max=args.theta_max,
            theta_scale=args.theta_scale,
            bins=args.nbins,
            multipoles=args.multipoles,
            threads=args.threads,
            use_log_bins=not args.linear_bins,
            tree_theta=args.tree_theta,
            nsmooth=args.nsmooth,
            phi_left=args.phiL,
            phi_right=args.phiR,
            theta_left=args.thetaL,
            theta_right=args.thetaR,
            options=args.more_options,
            result_type=args.result_type,
            verbose=args.verbose,
            verbose_log=args.verbose_log,
            continue_on_error=args.continue_on_error,
            edge_corrections=args.edge_corrections,
        ).normalized()
        engines = resolve_engines(
            args.engines, available, edge_corrections=config.wants_edge_corrections,
            masked=args.mask is not None or "read-mask" in config.options,
        )
        needs_mpi = any(KAPPA_ENGINES[engine].mpi for engine in engines)
        if (needs_mpi and args.mpi_ranks > 1 and mpi_environment_size() == 1
                and os.environ.get(MPI_CHILD_ENV) != "1"):
            return spawn_mpi(args)

        # Python owns MPI for the entire suite, including one-rank MPI runs.
        # An engine cleanup must not finalize MPI before the next engine runs.
        python_mpi_required = needs_mpi or mpi_environment_size() > 1
        comm = get_mpi_comm(python_mpi_required)
        if os.environ.get(MPI_CHILD_ENV) == "1" and comm.size != args.mpi_ranks:
            raise RuntimeError(
                f"mpiexec started {comm.size} ranks, expected {args.mpi_ranks}"
            )

        def load_catalog():
            if args.fits is not None:
                catalog = catalog_from_healpix(
                    args.fits,
                    field_index=args.field,
                    nside_down=args.nside_down,
                    mask_path=args.mask,
                    mask_threshold=args.mask_threshold,
                    center_field=not args.no_center_field,
                )
            elif args.catalog_npz is not None:
                if args.mask is not None:
                    raise ValueError("--mask is only valid with --fits")
                catalog = catalog_from_npz(
                    args.catalog_npz, center_field=not args.no_center_field
                )
            elif args.synthetic_nside is not None:
                if args.mask is not None:
                    raise ValueError("--mask is only valid with --fits")
                catalog = synthetic_healpix_catalog(
                    args.synthetic_nside, center_field=not args.no_center_field
                )
            else:
                raise ValueError(
                    "select --fits, --catalog-npz, or --synthetic-nside"
                )
            if args.save_catalog_npz is not None:
                save_catalog_npz(args.save_catalog_npz, catalog)
            print(
                f"Catalog loaded once on rank 0: {catalog.nbody} bodies",
                flush=True,
            )
            return catalog
        catalog = None
        load_error = None
        if comm.rank == 0:
            try:
                catalog = load_catalog()
            except Exception as exc:
                load_error = f"{type(exc).__name__}: {exc}"
        load_error = comm.bcast(load_error, root=0)
        if load_error:
            raise RuntimeError(load_error)
        catalog = broadcast_catalog(comm, catalog)

        # NPZ masks are only known after the one catalog read and broadcast.
        config.engines = tuple(resolve_engines(
            args.engines, available, edge_corrections=config.wants_edge_corrections,
            masked=catalog.mask is not None or "read-mask" in config.options,
        ))
        if comm.rank == 0:
            print(
                "Selected " + ("edge-corrected " if config.wants_edge_corrections else "")
                + "engines: " + ", ".join(config.engines), flush=True,
            )
        run_engine_suite(catalog, config, comm=comm)
        if comm.rank == 0:
            print(f"Results written to {Path(args.outdir).expanduser().resolve()}")
        return 0
    except (OSError, RuntimeError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
