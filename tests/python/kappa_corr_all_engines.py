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
Native-resolution implicit HEALPix files are scanned in bounded memory;
``--max-points`` provides reproducible bottom-hash thinning for comparison
and scaling runs without constructing a full-sky Cartesian catalog.
"""

from __future__ import annotations

import argparse
from contextlib import ExitStack
from dataclasses import dataclass, field
import json
import math
import os
from collections.abc import Mapping
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
import time
from typing import Any, Dict, Iterable, Optional, Sequence

import numpy as np


# healpy imports matplotlib. Keep a persistent writable cache when the home
# cache is locked, rather than rebuilding matplotlib's font cache on every run.
if "MPLCONFIGDIR" not in os.environ:
    _matplotlib_home = Path.home() / ".matplotlib"
    if not _matplotlib_home.is_dir() or not os.access(_matplotlib_home, os.W_OK):
        _matplotlib_cache = (
            Path(tempfile.gettempdir()) / f"ctreeballs-matplotlib-{os.getuid()}"
        )
        _matplotlib_cache.mkdir(parents=True, exist_ok=True)
        os.environ["MPLCONFIGDIR"] = os.fspath(_matplotlib_cache)


PROJECT_ROOT = Path(__file__).resolve().parents[2]


def _default_cballs_executable() -> Path:
    configured = os.environ.get("CTREEBALLS_CBALLS") or os.environ.get("CBALLS")
    if configured:
        return Path(configured).expanduser()
    source_executable = PROJECT_ROOT / "cballs"
    if source_executable.is_file():
        return source_executable
    path_executable = shutil.which("cballs")
    return Path(path_executable) if path_executable else source_executable


DEFAULT_CBALLS = _default_cballs_executable()
MPI_CHILD_ENV = "CTREEBALLS_KAPPA_MPI_CHILD"

# A directly executed source-tree script otherwise searches ``python/`` before
# the repository root and may load a stale site-installed extension.
if any(PROJECT_ROOT.glob("cyballs*.so")) or any(PROJECT_ROOT.glob("cyballs*.pyd")):
    _project_root = os.fspath(PROJECT_ROOT)
    sys.path[:] = [entry for entry in sys.path if entry != _project_root]
    sys.path.insert(0, _project_root)


@dataclass(frozen=True)
class EngineSpec:
    mpi: bool = False
    has_2pcf: bool = True
    has_3pcf: bool = True
    force_log_bins: bool = False
    supports_mask: bool = False
    supports_edge: bool = False
    supports_dual_node_bin_slop: bool = False
    supports_smooth_pivot: bool = False
    note: str = ""


# Scalar angular/convergence methods accepted by this driver.  Availability is
# still determined from the active executable, so optional profiles do not
# appear in ``all`` unless they were compiled.
KAPPA_ENGINES: Dict[str, EngineSpec] = {
    "kdtree-2balls-omp": EngineSpec(
        supports_mask=True, supports_edge=True,
        supports_dual_node_bin_slop=True, supports_smooth_pivot=True,
    ),
    "kdtree-2balls-mpi": EngineSpec(
        mpi=True, supports_mask=True, supports_edge=True,
        supports_dual_node_bin_slop=True, supports_smooth_pivot=True,
    ),
    "balltree-2balls-omp": EngineSpec(
        supports_mask=True, supports_edge=True,
        supports_dual_node_bin_slop=True, supports_smooth_pivot=True,
    ),
    "balltree-2balls-mpi": EngineSpec(
        mpi=True, supports_mask=True, supports_edge=True,
        supports_dual_node_bin_slop=True, supports_smooth_pivot=True,
    ),
    "octree-2balls-omp": EngineSpec(
        supports_mask=True, supports_edge=True, supports_dual_node_bin_slop=True
    ),
    "octree-2balls-mpi": EngineSpec(
        mpi=True, supports_mask=True, supports_edge=True,
        supports_dual_node_bin_slop=True,
    ),
}


LYA_FOREST_ENGINES = tuple(
    f"{name}-{parallel}"
    for parallel in ("omp", "mpi")
    for name in (
        "lya-2pcf", "lya-3pcf", "lya-2pcf-3pcf",
        "lya-1d-2pcf", "lya-1d-3pcf", "lya-1d-2pcf-3pcf",
        "lya-1d-tree-2pcf", "lya-1d-tree-3pcf",
    )
)

SHEAR_FIELD_ENGINES = (
    "octree-shear-sphere-2balls-omp",
    "kdtree-shear-sphere-2balls-omp",
    "balltree-shear-sphere-2balls-omp",
)


INCOMPATIBLE_ENGINE_REASONS = {
    "kdtree-box-omp": "is a periodic Cartesian-box estimator",
    "neighbor-boxes-omp": "is a periodic Cartesian-box estimator",
    "octree-3pcf-3d-omp": "computes a physical 3D statistic, not angular kappa",
    "octree-3pcf-3d-mpi": "computes a physical 3D statistic, not angular kappa",
}

INCOMPATIBLE_ENGINE_REASONS.update({
    method: (
        "requires gamma1/gamma2 spin-2 data and its documented projection/transport; "
        "use tests/python/shear_corr_all_engines.py"
    )
    for method in SHEAR_FIELD_ENGINES
})

INCOMPATIBLE_ENGINE_REASONS.update({
    method: ("requires forest IDs and observer distances; use "
             "tests/python/lya_corr_all_engines.py and set_forest_catalog()")
    for method in LYA_FOREST_ENGINES
})

METHOD_ALIASES = {}


@dataclass(frozen=True)
class AngularPatch:
    """Open longitude/colatitude rectangle, with bounds in degrees."""

    phi_left: float = 0.0
    phi_right: float = 90.0
    theta_left: float = 0.0
    theta_right: float = 90.0

    def __post_init__(self) -> None:
        if not all(math.isfinite(value) for value in (
            self.phi_left, self.phi_right, self.theta_left, self.theta_right,
        )):
            raise ValueError("patch bounds must be finite degrees")
        if not 0.0 <= self.phi_left < self.phi_right <= 360.0:
            raise ValueError("patch requires 0 <= phiL < phiR <= 360 degrees")
        if not 0.0 <= self.theta_left < self.theta_right <= 180.0:
            raise ValueError("patch requires 0 <= thetaL < thetaR <= 180 degrees")

    def select(self, theta: np.ndarray, phi: np.ndarray) -> np.ndarray:
        # Match the native FITS patch predicate, including open boundaries.
        return (
            (theta > math.radians(self.theta_left))
            & (theta < math.radians(self.theta_right))
            & (phi > math.radians(self.phi_left))
            & (phi < math.radians(self.phi_right))
        )

    def metadata(self) -> dict:
        return {
            "phiL": self.phi_left, "phiR": self.phi_right,
            "thetaL": self.theta_left, "thetaR": self.theta_right,
            "units": "degree", "theta_convention": "colatitude",
            "boundaries": "exclusive",
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
    dual_node_bin_slop: bool = False
    smooth_pivot_compiled: Optional[bool] = None
    plots: bool = True
    flatten_plots: bool = True
    patch: bool = False

    @property
    def angular_patch(self) -> Optional[AngularPatch]:
        options = _split_options(self.options)
        if not self.patch and "patch" not in options:
            return None
        if "patch-with-all" in options:
            raise ValueError("patch-with-all is not supported with the shared patch filter")
        return AngularPatch(
            self.phi_left, self.phi_right, self.theta_left, self.theta_right,
        )

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
        if self.theta_scale == "degree" and self.theta_max > 180.0:
            raise ValueError("angular sphere separations cannot exceed 180 degrees")
        if self.theta_scale == "radian" and self.theta_max > math.pi:
            raise ValueError("angular sphere separations cannot exceed pi radians")
        if self.result_type not in {"sincos", "edge_effects"}:
            raise ValueError("result_type must be sincos or edge_effects")
        options = tuple(_split_options(self.options))
        patch = self.angular_patch
        edge = self.wants_edge_corrections
        if edge and "only-2pcf" in options:
            raise ValueError("edge corrections require 3PCF; remove only-2pcf")
        return RunConfig(
            **{
                **self.__dict__,
                "engines": engines,
                "output_dir": output_dir,
                "options": options,
                "patch": patch is not None,
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


def merge_statistics_options(
    statistics: str, values: Iterable[str] | str | None,
) -> tuple[str, ...]:
    """Translate the public statistics selector to native runtime options."""
    if statistics not in {"2pcf", "3pcf", "both"}:
        raise ValueError("statistics must be 2pcf, 3pcf, or both")
    options = _split_options(values)
    has_only_2pcf = "only-2pcf" in options
    has_only_3pcf = "only-3pcf" in options
    if has_only_2pcf and has_only_3pcf:
        raise ValueError("only-2pcf and only-3pcf cannot be combined")
    if statistics == "2pcf":
        if has_only_3pcf:
            raise ValueError("--statistics 2pcf conflicts with options=only-3pcf")
        if not has_only_2pcf:
            options.append("only-2pcf")
    elif statistics == "3pcf":
        if has_only_2pcf:
            raise ValueError("--statistics 3pcf conflicts with options=only-2pcf")
        if not has_only_3pcf:
            options.append("only-3pcf")
    return tuple(options)


def statistics_from_options(values: Iterable[str] | str | None) -> str:
    options = set(_split_options(values))
    if "only-2pcf" in options:
        return "2pcf"
    if "only-3pcf" in options:
        return "3pcf"
    return "both"


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


def discover_make_settings(executable: Path = DEFAULT_CBALLS) -> dict[str, str]:
    executable = Path(executable).expanduser().resolve()
    completed = subprocess.run(
        [os.fspath(executable), "options=make-info"],
        cwd=executable.parent,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    if completed.returncode not in (0, 1):
        raise RuntimeError(
            f"could not inspect the cballs build profile:\n{completed.stdout[-2000:]}"
        )
    return dict(re.findall(
        r"^\s*([A-Za-z_][A-Za-z0-9_]*)\s*=\s*(.*?)\s*$",
        completed.stdout,
        re.MULTILINE,
    ))


def discover_cython_methods(candidates: Sequence[str]) -> list[str]:
    try:
        from cyballs import search_method_id
    except (ImportError, OSError) as exc:
        raise RuntimeError(
            "the installed cyballs extension lacks search_method_id(); "
            "rebuild it with `make cyballs` from this source tree"
        ) from exc
    return [name for name in candidates if search_method_id(name) >= 0]


def print_engine_table(
    available: Sequence[str], smooth_pivot_compiled: Optional[bool] = None
) -> None:
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
        if not spec.supports_smooth_pivot:
            smooth = "smooth=unsupported"
        elif smooth_pivot_compiled is True:
            smooth = "smooth=default-on"
        elif smooth_pivot_compiled is False:
            smooth = "smooth=not-compiled"
        else:
            smooth = "smooth=build-dependent"
        print(
            f"  {name:38s} {state:11s} {parallel:12s} {products:9s} "
            f"{mask} {edge} {smooth} frontier=dual-node"
        )
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
    for token in tokens or ("octree-2balls-omp",):
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
                raise ValueError(
                    f"{name} is not available in the active cballs/cyballs build"
                )
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
        minimum = math.radians(config.theta_min)
        maximum = math.radians(config.theta_max)
    elif config.theta_scale == "radian":
        minimum, maximum = config.theta_min, config.theta_max
    else:
        return config.theta_min, config.theta_max
    return 2.0 * math.sin(0.5 * minimum), 2.0 * math.sin(0.5 * maximum)


def _healpix_valid(values: np.ndarray, hp: Any) -> np.ndarray:
    return np.isfinite(values) & ~hp.mask_bad(values)


def _implicit_healpix_column(
    stack: ExitStack, path: Path, field_index: int,
) -> tuple[np.ndarray, int, bool, list[str]]:
    """Return a memory-mapped implicit HEALPix column and its layout."""
    try:
        from astropy.io import fits
    except ImportError as exc:
        raise RuntimeError("sparse HEALPix FITS input requires astropy") from exc

    path = Path(path).expanduser().resolve()
    hdus = stack.enter_context(fits.open(path, memmap=True, lazy_load_hdus=True))
    tables = [
        hdu for hdu in hdus
        if hasattr(hdu, "columns") and hdu.data is not None
    ]
    if not tables:
        raise ValueError(f"no binary table was found in {path}")
    table = tables[0]
    scheme = str(table.header.get("INDXSCHM", "IMPLICIT")).strip().upper()
    if scheme != "IMPLICIT":
        raise ValueError(
            f"sparse loading requires an IMPLICIT HEALPix table, got {scheme!r}"
        )
    nside = int(table.header.get("NSIDE", 0))
    if nside <= 0:
        raise ValueError(f"missing or invalid NSIDE in {path}")
    ordering = str(table.header.get("ORDERING", "RING")).strip().upper()
    if ordering not in {"RING", "NESTED", "NEST"}:
        raise ValueError(f"unsupported HEALPix ORDERING={ordering!r} in {path}")
    try:
        column = np.asarray(table.data.field(field_index)).reshape(-1)
    except (IndexError, TypeError, ValueError) as exc:
        raise ValueError(
            f"field {field_index} is unavailable in {path}; "
            f"fields are {table.columns.names}"
        ) from exc
    expected = 12*nside*nside
    if column.size != expected:
        raise ValueError(
            f"HEALPix field has {column.size} pixels, expected {expected} "
            f"for NSIDE={nside}"
        )
    return column, nside, ordering != "RING", list(table.columns.names)


def _splitmix64(indices: np.ndarray, seed: int) -> np.ndarray:
    """Return reproducible, platform-independent priorities for pixel IDs."""
    with np.errstate(over="ignore"):
        value = np.asarray(indices, dtype=np.uint64) + np.uint64(seed)
        value += np.uint64(0x9E3779B97F4A7C15)
        value = (value ^ (value >> np.uint64(30))) * np.uint64(0xBF58476D1CE4E5B9)
        value = (value ^ (value >> np.uint64(27))) * np.uint64(0x94D049BB133111EB)
        return value ^ (value >> np.uint64(31))


def _retain_lowest_priorities(
    retained_pixels: np.ndarray,
    retained_priorities: np.ndarray,
    pixels: np.ndarray,
    priorities: np.ndarray,
    maximum: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Merge one chunk into a deterministic bottom-hash sample."""
    if pixels.size > maximum:
        keep = np.argpartition(priorities, maximum - 1)[:maximum]
        pixels, priorities = pixels[keep], priorities[keep]
    if retained_pixels.size:
        pixels = np.concatenate((retained_pixels, pixels))
        priorities = np.concatenate((retained_priorities, priorities))
    if pixels.size > maximum:
        keep = np.argpartition(priorities, maximum - 1)[:maximum]
        pixels, priorities = pixels[keep], priorities[keep]
    return pixels, priorities


def _native_resolution_healpix_catalog(
    fits_path: Path,
    field_index: int,
    mask_path: Optional[Path],
    mask_threshold: float,
    center_field: bool,
    max_points: int,
    sampling_seed: int,
    chunk_pixels: int,
    patch: Optional[AngularPatch] = None,
) -> KappaCatalog:
    """Build a catalog by scanning memory-mapped map columns in chunks."""
    try:
        import healpy as hp
    except ImportError as exc:
        raise RuntimeError("HEALPix FITS input requires healpy") from exc

    if max_points < 0 or 0 < max_points < 3:
        raise ValueError("max_points must be zero or at least three")
    if not 0 <= sampling_seed < 2**64:
        raise ValueError("sampling_seed must lie in [0, 2**64)")
    if chunk_pixels < 1024:
        raise ValueError("fits_chunk_pixels must be at least 1024")

    with ExitStack() as stack:
        values, nside, nested, fields = _implicit_healpix_column(
            stack, fits_path, field_index
        )
        mask_values = None
        resolved_mask = None
        if mask_path is not None:
            resolved_mask = Path(mask_path).expanduser().resolve()
            mask_values, mask_nside, mask_nested, _ = _implicit_healpix_column(
                stack, resolved_mask, 0
            )
            if mask_nside != nside or mask_nested != nested:
                raise ValueError(
                    "native-resolution sparse loading requires map and mask "
                    "to have the same NSIDE and ORDERING"
                )

        retained_pixels = np.empty(0, dtype=np.int64)
        retained_priorities = np.empty(0, dtype=np.uint64)
        pixel_chunks: list[np.ndarray] = []
        eligible_count = 0
        before_patch = 0
        for start in range(0, values.size, chunk_pixels):
            stop = min(start + chunk_pixels, values.size)
            selected = _healpix_valid(values[start:stop], hp)
            if mask_values is not None:
                mask_chunk = mask_values[start:stop]
                selected &= _healpix_valid(mask_chunk, hp)
                selected &= mask_chunk > mask_threshold
            local = np.flatnonzero(selected).astype(np.int64, copy=False)
            if not local.size:
                continue
            local += start
            before_patch += int(local.size)
            if patch is not None:
                theta, phi = hp.pix2ang(nside, local, nest=nested)
                local = local[patch.select(theta, phi)]
                if not local.size:
                    continue
            eligible_count += int(local.size)
            if max_points:
                priorities = _splitmix64(local, sampling_seed)
                retained_pixels, retained_priorities = _retain_lowest_priorities(
                    retained_pixels, retained_priorities,
                    local, priorities, max_points,
                )
            else:
                pixel_chunks.append(local)

        if max_points:
            pixels = retained_pixels
        elif pixel_chunks:
            pixels = np.concatenate(pixel_chunks)
        else:
            pixels = np.empty(0, dtype=np.int64)
        pixels.sort()
        if pixels.size < 3:
            raise ValueError("the map, mask, and patch select fewer than three valid pixels")
        kappa = np.ascontiguousarray(values[pixels], dtype=np.float64)

    x, y, z = hp.pix2vec(nside, pixels, nest=nested)
    positions = np.ascontiguousarray(np.column_stack((x, y, z)), dtype=np.float64)
    selected_mask = (
        np.ones(pixels.size, dtype=np.uint8) if mask_path is not None else None
    )
    if center_field:
        kappa = kappa.copy()
        kappa -= np.mean(kappa, dtype=np.float64)
    return KappaCatalog(
        positions=positions,
        kappa=kappa,
        mask=selected_mask,
        metadata={
            "source": os.fspath(Path(fits_path).expanduser().resolve()),
            "mask_source": os.fspath(resolved_mask) if resolved_mask else None,
            "nside_input": nside,
            "nside": nside,
            "ordering": "NESTED" if nested else "RING",
            "field": field_index,
            "source_fields": fields,
            "centered": center_field,
            "mask_preselected": mask_path is not None,
            "patch": patch.metadata() if patch is not None else None,
            "eligible_pixels_before_patch": before_patch,
            "eligible_pixels_before_thinning": eligible_count,
            "max_points": max_points,
            "sampling_seed": sampling_seed if max_points else None,
            "sampling_method": "splitmix64-bottom-k" if max_points else "none",
            "fits_chunk_pixels": chunk_pixels,
            "loader": "memory-mapped-native-resolution",
        },
    ).normalized()


def catalog_from_healpix(
    fits_path: Path,
    field_index: int = 0,
    nside_down: int = 0,
    mask_path: Optional[Path] = None,
    mask_threshold: float = 0.0,
    center_field: bool = True,
    max_points: int = 0,
    sampling_seed: int = 8675309,
    chunk_pixels: int = 1 << 20,
    patch: Optional[AngularPatch] = None,
) -> KappaCatalog:
    try:
        import healpy as hp
    except ImportError as exc:
        raise RuntimeError("HEALPix FITS input requires healpy") from exc
    if max_points < 0 or 0 < max_points < 3:
        raise ValueError("max_points must be zero or at least three")
    if not 0 <= sampling_seed < 2**64:
        raise ValueError("sampling_seed must lie in [0, 2**64)")
    if chunk_pixels < 1024:
        raise ValueError("fits_chunk_pixels must be at least 1024")

    fits_path = Path(fits_path).expanduser().resolve()
    values = None
    sparse_available = True
    try:
        with ExitStack() as stack:
            _, nside_in, _, _ = _implicit_healpix_column(
                stack, fits_path, field_index
            )
    except ValueError as exc:
        if "requires an IMPLICIT HEALPix table" not in str(exc):
            raise
        sparse_available = False
        values = np.asarray(
            hp.read_map(fits_path, field=field_index, dtype=np.float64),
            dtype=np.float64,
        )
        nside_in = int(hp.get_nside(values))

    target_nside = nside_down or nside_in
    if not hp.isnsideok(target_nside):
        raise ValueError("nside_down must be a valid HEALPix NSIDE")
    if target_nside > nside_in:
        raise ValueError("nside_down cannot exceed the input NSIDE")
    if target_nside == nside_in and sparse_available:
        try:
            return _native_resolution_healpix_catalog(
                fits_path, field_index, mask_path, mask_threshold, center_field,
                max_points, sampling_seed, chunk_pixels,
                patch,
            )
        except ValueError as exc:
            if "map and mask to have the same NSIDE and ORDERING" not in str(exc):
                raise

    if values is None:
        values = np.asarray(
            hp.read_map(fits_path, field=field_index, dtype=np.float64),
            dtype=np.float64,
        )
    if target_nside != nside_in:
        values = hp.ud_grade(
            values, nside_out=target_nside,
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

    eligible = valid & mask_full
    pixels = np.flatnonzero(eligible)
    before_patch = int(pixels.size)
    if patch is not None:
        keep_patch = np.empty(pixels.size, dtype=bool)
        for start in range(0, pixels.size, chunk_pixels):
            stop = min(start + chunk_pixels, pixels.size)
            theta, phi = hp.pix2ang(nside, pixels[start:stop], nest=False)
            keep_patch[start:stop] = patch.select(theta, phi)
        pixels = pixels[keep_patch]
    eligible_count = int(pixels.size)
    if max_points:
        if max_points < 3:
            raise ValueError("max_points must be zero or at least three")
        priorities = _splitmix64(pixels, sampling_seed)
        if pixels.size > max_points:
            keep = np.argpartition(priorities, max_points - 1)[:max_points]
            pixels = np.sort(pixels[keep])
    if pixels.size < 3:
        raise ValueError("the map, mask, and patch select fewer than three valid pixels")
    x, y, z = hp.pix2vec(nside, pixels, nest=False)
    positions = np.ascontiguousarray(np.column_stack((x, y, z)), dtype=np.float64)
    kappa = np.ascontiguousarray(values[pixels], dtype=np.float64)
    selected_mask = (
        np.ones(pixels.size, dtype=np.uint8) if mask_values is not None else None
    )
    if center_field:
        kappa = kappa.copy()
        kappa -= np.mean(kappa, dtype=np.float64)
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
            "mask_preselected": mask_values is not None,
            "patch": patch.metadata() if patch is not None else None,
            "eligible_pixels_before_patch": before_patch,
            "eligible_pixels_before_thinning": eligible_count,
            "max_points": max_points,
            "sampling_seed": sampling_seed if max_points else None,
            "sampling_method": "splitmix64-bottom-k" if max_points else "none",
            "loader": (
                "healpy-downgrade" if target_nside != nside_in
                else "healpy-full-map"
            ),
        },
    ).normalized()


def synthetic_healpix_catalog(
    nside: int = 4, center_field: bool = True,
    patch: Optional[AngularPatch] = None,
) -> KappaCatalog:
    try:
        import healpy as hp
    except ImportError as exc:
        raise RuntimeError("synthetic HEALPix input requires healpy") from exc
    if not hp.isnsideok(nside):
        raise ValueError("synthetic_nside must be a valid HEALPix NSIDE")
    pixels = np.arange(hp.nside2npix(nside))
    before_patch = int(pixels.size)
    if patch is not None:
        theta, phi = hp.pix2ang(nside, pixels, nest=False)
        pixels = pixels[patch.select(theta, phi)]
        if pixels.size < 3:
            raise ValueError("the patch selects fewer than three valid pixels")
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
        metadata={
            "source": "synthetic", "nside": nside, "centered": center_field,
            "patch": patch.metadata() if patch is not None else None,
            "eligible_pixels_before_patch": before_patch,
            "eligible_pixels_before_thinning": int(pixels.size),
        },
    ).normalized()


def filter_catalog_patch(
    catalog: KappaCatalog, patch: Optional[AngularPatch],
) -> KappaCatalog:
    """Filter already-loaded arrays without moving the observer or input rows."""
    if patch is None or catalog.metadata.get("patch") == patch.metadata():
        return catalog
    catalog = catalog.normalized()
    radius = np.linalg.norm(catalog.positions, axis=1)
    if np.any(radius == 0.0) or not np.all(np.isfinite(radius)):
        raise ValueError("angular patch selection requires finite nonzero position radii")
    theta = np.arccos(np.clip(catalog.positions[:, 2] / radius, -1.0, 1.0))
    phi = np.mod(np.arctan2(catalog.positions[:, 1], catalog.positions[:, 0]), 2*np.pi)
    selected = patch.select(theta, phi)
    active = (catalog.mask.astype(bool) if catalog.mask is not None
              else np.ones(catalog.nbody, dtype=bool))
    if np.count_nonzero(selected & active) < 3:
        raise ValueError("the patch and mask select fewer than three valid points")
    kappa = catalog.kappa[selected].copy()
    if catalog.metadata.get("centered", False):
        kappa -= np.mean(kappa[active[selected]], dtype=np.float64)
    return KappaCatalog(
        positions=catalog.positions[selected], kappa=kappa,
        weights=catalog.weights[selected] if catalog.weights is not None else None,
        mask=catalog.mask[selected] if catalog.mask is not None else None,
        metadata={
            **catalog.metadata, "patch": patch.metadata(),
            "eligible_pixels_before_patch": int(np.count_nonzero(active)),
            "eligible_pixels_before_thinning": int(np.count_nonzero(selected & active)),
        },
    ).normalized()


def catalog_from_npz(
    path: Path, center_field: bool = False,
    patch: Optional[AngularPatch] = None,
) -> KappaCatalog:
    path = Path(path).expanduser().resolve()
    with np.load(path, allow_pickle=False) as archive:
        catalog = KappaCatalog(
            positions=archive["positions"],
            kappa=archive["kappa"],
            weights=archive["weights"] if "weights" in archive else None,
            mask=archive["mask"] if "mask" in archive else None,
            metadata={"source": os.fspath(path)},
        ).normalized()
    catalog = filter_catalog_patch(catalog, patch)
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
    if engine not in KAPPA_ENGINES:
        raise ValueError(f"{engine} does not support this active scalar driver")
    spec = KAPPA_ENGINES[engine]
    rmin, rmax = angular_limits(config)
    angular_bins_needed = max(4, 2 * config.multipoles + 1)
    angular_bins = 1 << (angular_bins_needed - 1).bit_length()
    options = ["compute-HistN", "and-CF", "out-m-HistZeta", "KKKCorrelation"]
    # The shared Python catalog has already been cut; do not filter it again in C.
    options.extend(option for option in _split_options(config.options) if option != "patch")
    if masked and "read-mask" not in options:
        options.append("read-mask")
    options = _split_options(options)
    if "smooth-pivot" in options and not spec.supports_smooth_pivot:
        raise ValueError(f"{engine} does not support options=smooth-pivot")
    if "smooth-pivot" in options and config.smooth_pivot_compiled is False:
        raise ValueError("options=smooth-pivot requires SMOOTHPIVOTON=1")
    if "smooth-pivot" in options and "no-smooth-pivot" in options:
        raise ValueError("options=smooth-pivot and no-smooth-pivot conflict")
    if config.wants_edge_corrections:
        if not spec.supports_edge:
            raise ValueError(f"{engine} does not support angular edge corrections")
        if "only-2pcf" in options:
            raise ValueError("edge corrections require 3PCF; remove only-2pcf")
        options = _split_options(
            [*options, "edge-corrections", "no-normalize-HistZeta"]
        )
    if config.dual_node_bin_slop and spec.supports_dual_node_bin_slop:
        options = _split_options([*options, "dual-node-bin-slop"])
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


def smooth_pivot_mode(config: RunConfig, engine: str) -> str:
    spec = KAPPA_ENGINES[engine]
    options = set(_split_options(config.options))
    if not spec.supports_smooth_pivot:
        return "unsupported"
    if "no-smooth-pivot" in options:
        return "disabled-by-option"
    if config.smooth_pivot_compiled is False:
        return "not-compiled"
    if "smooth-pivot" in options:
        return "enabled-explicitly"
    return "enabled-by-build-default" if config.smooth_pivot_compiled else "build-default"


def settings_to_json(value: Any) -> Any:
    """Detach immutable run settings into JSON-compatible containers."""
    if isinstance(value, Mapping):
        return {key: settings_to_json(item) for key, item in value.items()}
    if isinstance(value, tuple):
        return [settings_to_json(item) for item in value]
    return value


def copy_engine_results(balls: Any, engine: str, config: RunConfig) -> dict:
    spec = KAPPA_ENGINES[engine]
    options = set(_split_options(config.options))
    edge = config.wants_edge_corrections
    want_2pcf = spec.has_2pcf and "only-3pcf" not in options
    want_3pcf = spec.has_3pcf and "only-2pcf" not in options
    result: Dict[str, Any] = {
        "engine": engine,
        "smooth_pivot": smooth_pivot_mode(config, engine),
        "cpu_time": float(balls.getCPUTime()),
        "nbody": int(balls.getNBody()),
        "result_type": "edge_effects" if edge else config.result_type,
        "edge_corrections": edge,
        "run_settings": settings_to_json(balls.run_settings),
        "warnings": [],
    }
    if want_3pcf:
        result["three_pcf_normalization"] = (
            "edge-corrected-mode-coupling" if edge else
            "raw-distinct-triplet-sum"
            if "no-normalize-HistZeta" in options else
            "per-bin-distinct-triplet-weight"
        )
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
                    result[f"zeta_m_{order - 1}"] = (
                        components[0] + components[1]
                        + 1j * (components[2] - components[3])
                    )
            else:
                diagnostics = balls.getScalarWindowDiagnostics()
                for name in ("status", "valid", "window_monopole", "pivot_ratio"):
                    result[f"scalar_window_{name}"] = diagnostics[name].copy()
                valid = diagnostics["valid"]
                if valid.shape != (config.bins, config.bins):
                    raise ValueError("scalar window validity has unexpected shape")
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
                        if not np.all(np.isfinite(value[valid])):
                            raise ValueError(f"{key} contains nonfinite valid estimates")
                        if not np.all(np.isnan(value[~valid])):
                            raise ValueError(f"{key} must mark unsupported bins as NaN")
                    result[f"zeta_edge_complex_{order}"] = (
                        result[f"zeta_edge_{order}"]
                        + 1j * result[f"zeta_edge_im_{order}"]
                    )
                    result[f"zeta_m_{order - 1}"] = result[
                        f"zeta_edge_complex_{order}"
                    ]
        except Exception as exc:
            if edge:
                raise RuntimeError(
                    f"{engine}: requested edge-corrected 3PCF is unavailable: {exc}"
                ) from exc
            result["warnings"].append(f"3PCF multipoles unavailable: {exc}")
    return result


def solve_scalar_mode_coupling(
    signal: np.ndarray, window: np.ndarray, max_n: int,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Solve C[ell,n] zeta[n] = signal[ell]/window[0] per radial bin."""
    signal = np.asarray(signal, dtype=np.complex128)
    window = np.asarray(window, dtype=np.complex128)
    if max_n < 0:
        raise ValueError("max_n must be non-negative")
    if signal.ndim != 3 or signal.shape[0] != signal.shape[1]:
        raise ValueError(f"signal must have shape (B, B, 2M+1), got {signal.shape}")
    bins = signal.shape[0]
    multipoles = 2 * max_n + 1
    if signal.shape != (bins, bins, multipoles):
        raise ValueError("inconsistent scalar signal multipole dimensions")
    if window.shape != (bins, bins, 4 * max_n + 1):
        raise ValueError("inconsistent scalar window multipole dimensions")
    if not np.all(np.isfinite(signal)) or not np.all(np.isfinite(window)):
        raise ValueError("scalar signal/window multipoles contain nonfinite values")

    orders = np.arange(-max_n, max_n + 1)
    corrected = np.full_like(signal, complex(np.nan, np.nan))
    status = np.full((bins, bins), 2, dtype=np.uint8)
    ratio = np.full((bins, bins), np.nan)
    diagnostics = {"empty_bins": 0, "singular_bins": 0, "nonfinite_bins": 0,
                   "status": status, "window_monopole": window[:, :, 2*max_n].real.copy(),
                   "pivot_ratio": ratio}
    for bin_1 in range(bins):
        for bin_2 in range(bins):
            n_zero = window[bin_1, bin_2, 2 * max_n]
            if n_zero.real <= 0 or abs(n_zero) <= np.finfo(float).tiny:
                diagnostics["empty_bins"] += 1
                continue
            matrix = np.empty((multipoles, multipoles), dtype=np.complex128)
            matrix_scale = 0.0
            for row, ell in enumerate(orders):
                for column, order in enumerate(orders):
                    value = window[bin_1, bin_2, ell - order + 2 * max_n]
                    matrix[row, column] = value / n_zero
                    matrix_scale = max(matrix_scale, abs(value) ** 2)
            rhs = signal[bin_1, bin_2, :] / n_zero
            tolerance = 128.0 * np.finfo(float).eps * (
                1.0 + math.sqrt(matrix_scale / abs(n_zero) ** 2)
            )
            singular = False
            nonfinite = not (np.all(np.isfinite(matrix)) and np.all(np.isfinite(rhs)))
            smallest_pivot, largest_pivot = math.inf, 0.0
            for column in range(multipoles):
                if nonfinite:
                    break
                pivot = column + int(np.argmax(np.abs(matrix[column:, column]) ** 2))
                pivot_size = abs(matrix[pivot, column])
                if not np.isfinite(pivot_size):
                    nonfinite = True
                    break
                if pivot_size ** 2 <= tolerance ** 2:
                    singular = True
                    break
                smallest_pivot = min(smallest_pivot, pivot_size)
                largest_pivot = max(largest_pivot, pivot_size)
                if pivot != column:
                    matrix[[column, pivot], :] = matrix[[pivot, column], :]
                    rhs[[column, pivot]] = rhs[[pivot, column]]
                divisor = matrix[column, column]
                matrix[column, :] /= divisor
                rhs[column] /= divisor
                for row in range(multipoles):
                    if row == column:
                        continue
                    factor = matrix[row, column]
                    if factor != 0.0:
                        matrix[row, :] -= factor * matrix[column, :]
                        rhs[row] -= factor * rhs[column]
            if nonfinite or not np.all(np.isfinite(rhs)):
                diagnostics["nonfinite_bins"] += 1
                status[bin_1, bin_2] = 4
            elif singular:
                diagnostics["singular_bins"] += 1
                status[bin_1, bin_2] = 3
                ratio[bin_1, bin_2] = 0.0
            else:
                status[bin_1, bin_2] = 1
                ratio[bin_1, bin_2] = smallest_pivot/largest_pivot
                corrected[bin_1, bin_2, :] = rhs
    diagnostics["valid"] = status == 1
    return corrected, diagnostics


def _timing_metadata(
    setup_wall: float, setup_cpu: float, compute_wall: float,
    compute_cpu: float, scope: str,
) -> dict[str, Any]:
    return {
        "setup_wall_time": float(setup_wall),
        "setup_cpu_time": float(setup_cpu),
        "compute_wall_time": float(compute_wall),
        "compute_cpu_time": float(compute_cpu),
        "total_wall_time": float(setup_wall + compute_wall),
        "total_cpu_time": float(setup_cpu + compute_cpu),
        "timing_scope": scope,
    }


def aggregate_rank_timings(rows: Sequence[dict[str, Any]]) -> dict[str, Any]:
    """Use critical-path wall time and consumed CPU time, never rank-zero alone."""
    if not rows:
        raise ValueError("at least one participating rank is required")
    result = {
        name: (max if "wall" in name else sum)(float(row[name]) for row in rows)
        for name in ("setup_wall_time", "setup_cpu_time", "compute_wall_time",
                     "compute_cpu_time", "total_wall_time", "total_cpu_time",
                     "native_reported_cpu_time")
    }
    result.update(ranks=len(rows), rank_timings=list(rows),
                  timing_scope=rows[0]["timing_scope"] +
                  "; wall=max(participating ranks), CPU=sum(participating ranks)")
    return result



def flatten_radial_matrix(values: np.ndarray) -> np.ndarray:
    """Flatten (radial-bin 1, radial-bin 2) in stable row-major order."""
    matrix = np.asarray(values)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError(f"3PCF radial matrix must be square, got {matrix.shape}")
    return matrix.reshape(-1, order="C")


def result_multipoles(result: dict) -> dict[int, np.ndarray]:
    modes: dict[int, np.ndarray] = {}
    for key, value in result.items():
        match = re.fullmatch(r"zeta_m_(\d+)", key)
        if match and isinstance(value, np.ndarray):
            modes[int(match.group(1))] = np.asarray(value)
    return modes


def _finite_limit(arrays: Sequence[np.ndarray]) -> float:
    finite = [
        np.abs(np.asarray(value)[np.isfinite(value)])
        for value in arrays if np.any(np.isfinite(value))
    ]
    return max((float(np.max(value)) for value in finite if value.size), default=1.0e-15)


def compare_result_arrays(results: dict[str, dict]) -> dict[str, dict[str, Any]]:
    """Compare common 2PCF/3PCF products against the first engine."""
    if not results:
        return {}
    reference_name = next(iter(results))
    reference = results[reference_name]
    comparisons: dict[str, dict[str, Any]] = {}
    for name, result in list(results.items())[1:]:
        keys = sorted(
            key for key in reference.keys() & result.keys()
            if key in {"nn", "xi"} or re.fullmatch(r"zeta_m_\d+", key)
        )
        for key in keys:
            left = np.asarray(reference[key])
            right = np.asarray(result[key])
            if left.shape != right.shape:
                continue
            absolute = np.abs(right - left)
            relative = 2.0 * absolute / np.maximum(
                np.abs(right) + np.abs(left), 1.0e-15
            )
            finite = np.isfinite(absolute) & np.isfinite(relative)
            tag = f"{name}__vs__{reference_name}__{key}"
            comparisons[tag] = {
                "max_absolute": float(np.max(absolute[finite])) if np.any(finite) else None,
                "rms_absolute": float(np.sqrt(np.mean(absolute[finite] ** 2)))
                if np.any(finite) else None,
                "max_symmetric_relative": float(
                    np.max(relative[finite])
                ) if np.any(finite) else None,
                "unsupported_reference_bins": int(np.count_nonzero(~np.isfinite(left))),
                "unsupported_candidate_bins": int(np.count_nonzero(~np.isfinite(right))),
                "compared_bins": int(np.count_nonzero(finite)),
            }
    return comparisons


def make_plots(results: dict[str, dict], config: RunConfig) -> list[str]:
    """Write 2PCF curves, 3PCF heat maps, and Figure-7-style flattened modes."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plot_root = config.output_dir / "plots"
    plot_root.mkdir(parents=True, exist_ok=True)
    paths: list[str] = []
    pair_members = [
        (name, value) for name, value in results.items()
        if "r" in value and "xi" in value
    ]
    if pair_members:
        rows = 2 if len(pair_members) > 1 else 1
        fig, axes = plt.subplots(rows, 1, figsize=(8.0, 3.5 * rows), squeeze=False)
        for name, value in pair_members:
            axes[0, 0].plot(value["r"], value["xi"], marker=".", label=name)
        axes[0, 0].set(ylabel="xi", title="Convergence 2PCF")
        axes[0, 0].legend(fontsize=8)
        axes[0, 0].grid(True, linestyle=":", alpha=0.5)
        if rows == 2:
            reference_name, reference = pair_members[0]
            for name, value in pair_members[1:]:
                if np.shape(value["xi"]) != np.shape(reference["xi"]):
                    continue
                denominator = np.maximum(
                    np.maximum(np.abs(value["xi"]), np.abs(reference["xi"])),
                    1.0e-15,
                )
                difference = (value["xi"] - reference["xi"]) / denominator
                axes[1, 0].plot(value["r"], difference, marker=".", label=name)
            axes[1, 0].set(ylabel="relative difference", title=f"relative to {reference_name}")
            axes[1, 0].grid(True, linestyle=":", alpha=0.5)
            axes[1, 0].legend(fontsize=8)
        axes[-1, 0].set_xlabel("Euclidean chord separation")
        fig.tight_layout()
        path = plot_root / "two_pcf.png"
        fig.savefig(path, dpi=160)
        plt.close(fig)
        paths.append(os.fspath(path))

    multipoles = {name: result_multipoles(value) for name, value in results.items()}
    modes = sorted({mode for values in multipoles.values() for mode in values})
    for mode in modes:
        members = [(name, values[mode]) for name, values in multipoles.items() if mode in values]
        for component, selector in (("real", np.real), ("imag", np.imag)):
            arrays = [selector(value) for _, value in members]
            if component == "imag" and not any(np.any(value != 0.0) for value in arrays):
                continue
            limit = _finite_limit(arrays)
            fig, axes = plt.subplots(
                1, len(members), figsize=(4.8 * len(members), 4.1), squeeze=False
            )
            for axis, (name, _), values in zip(axes[0], members, arrays):
                image = axis.imshow(
                    values, origin="lower", aspect="equal", cmap="RdBu_r",
                    vmin=-limit, vmax=limit,
                )
                axis.set(
                    title=name, xlabel="radial bin 2", ylabel="radial bin 1"
                )
                fig.colorbar(image, ax=axis, shrink=0.82)
            fig.suptitle(f"3PCF mode m={mode}, {component} component")
            fig.tight_layout()
            path = plot_root / f"three_pcf_m{mode}_{component}.png"
            fig.savefig(path, dpi=160)
            plt.close(fig)
            paths.append(os.fspath(path))
        if len(members) > 1:
            reference_name, reference = members[0]
            candidates = members[1:]
            fig, axes = plt.subplots(
                1, len(candidates), figsize=(4.8 * len(candidates), 4.1), squeeze=False
            )
            for axis, (name, values) in zip(axes[0], candidates):
                relative = 2.0 * np.abs(values - reference) / np.maximum(
                    np.abs(values) + np.abs(reference), 1.0e-15
                )
                image = axis.imshow(relative, origin="lower", aspect="equal", cmap="magma")
                axis.set(
                    title=f"{name} vs {reference_name}",
                    xlabel="radial bin 2", ylabel="radial bin 1",
                )
                fig.colorbar(image, ax=axis, shrink=0.82)
            fig.suptitle(f"3PCF mode m={mode}, symmetric relative difference")
            fig.tight_layout()
            path = plot_root / f"three_pcf_m{mode}_relative.png"
            fig.savefig(path, dpi=160)
            plt.close(fig)
            paths.append(os.fspath(path))

    if config.flatten_plots and modes:
        for component, selector in (("real", np.real), ("imag", np.imag)):
            has_component = component == "real" or any(
                np.any(selector(values[mode]) != 0.0)
                for values in multipoles.values() for mode in values
            )
            if not has_component:
                continue
            fig, axes = plt.subplots(
                len(modes), 1, figsize=(10.0, max(3.2, 2.6 * len(modes))),
                squeeze=False, sharex=True,
            )
            for axis, mode in zip(axes[:, 0], modes):
                for name, values in multipoles.items():
                    if mode not in values:
                        continue
                    flattened = flatten_radial_matrix(selector(values[mode]))
                    axis.plot(np.arange(flattened.size), flattened, label=name)
                for boundary in range(config.bins, config.bins ** 2, config.bins):
                    axis.axvline(boundary - 0.5, color="0.82", linewidth=0.55)
                axis.set_ylabel(f"m={mode}")
                axis.grid(True, axis="y", linestyle=":", alpha=0.45)
                axis.legend(fontsize=7, ncol=2)
            axes[-1, 0].set_xlabel(
                "flattened radial-bin index (bin 1 major, bin 2 minor)"
            )
            fig.suptitle(f"3PCF flattened radial-bin matrices, {component} component")
            fig.tight_layout()
            path = plot_root / f"three_pcf_flattened_{component}.png"
            fig.savefig(path, dpi=160)
            plt.close(fig)
            paths.append(os.fspath(path))
    return paths


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


def timing_summary(results: dict[str, dict]) -> dict[str, dict[str, Any]]:
    """Return comparable setup/compute wall and process-CPU measurements."""
    summary: dict[str, dict[str, Any]] = {}
    for engine, result in results.items():
        compute_wall = float(result.get("compute_wall_time", 0.0))
        compute_cpu = float(result.get("compute_cpu_time", 0.0))
        summary[engine] = {
            "backend": "ctreeballs",
            "setup_wall_time": float(result.get("setup_wall_time", 0.0)),
            "setup_cpu_time": float(result.get("setup_cpu_time", 0.0)),
            "compute_wall_time": compute_wall,
            "compute_cpu_time": compute_cpu,
            "total_wall_time": float(result.get("total_wall_time", compute_wall)),
            "total_cpu_time": float(result.get("total_cpu_time", compute_cpu)),
            "compute_cpu_wall_ratio": (
                compute_cpu / compute_wall if compute_wall > 0.0 else None
            ),
            "native_reported_cpu_time": result.get("native_reported_cpu_time"),
            "timing_scope": result.get("timing_scope", "unspecified"),
        }
    return summary


def write_timing_report(path: Path, timings: dict[str, dict[str, Any]]) -> None:
    """Write a human-readable timing table alongside the JSON summary."""
    columns = (
        ("engine", 35), ("backend", 10), ("setup_wall_s", 14),
        ("compute_wall_s", 16), ("total_wall_s", 14),
        ("setup_cpu_s", 13), ("compute_cpu_s", 15), ("total_cpu_s", 13),
        ("native_cpu_s", 14), ("cpu/wall", 10),
    )
    lines = [" ".join(name.ljust(width) for name, width in columns)]
    lines.append(" ".join("-" * width for _, width in columns))
    for engine, values in timings.items():
        ratio = values["compute_cpu_wall_ratio"]
        native_cpu = values["native_reported_cpu_time"]
        fields = (
            engine,
            values["backend"],
            f'{values["setup_wall_time"]:.6f}',
            f'{values["compute_wall_time"]:.6f}',
            f'{values["total_wall_time"]:.6f}',
            f'{values["setup_cpu_time"]:.6f}',
            f'{values["compute_cpu_time"]:.6f}',
            f'{values["total_cpu_time"]:.6f}',
            "n/a" if native_cpu is None else f"{native_cpu:.6f}",
            "n/a" if ratio is None else f"{ratio:.3f}",
        )
        lines.append(" ".join(str(value).ljust(width) for value, (_, width) in zip(fields, columns)))
    lines.extend((
        "",
        "CPU values are process CPU seconds; they may exceed wall time for threaded work.",
        "MPI CPU is summed across participating ranks; wall time is their maximum.",
        "Read timing_scope in summary.json: setup, output and cleanup scopes differ by driver.",
    ))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_engine_suite(
    catalog: KappaCatalog,
    config: RunConfig,
    comm: Any = None,
) -> dict[str, dict]:
    """Run all selected engines while retaining one registered NumPy catalog."""
    config = config.normalized()
    if comm is None:
        comm = get_mpi_comm(
            mpi_environment_size() > 1
            or any(KAPPA_ENGINES[e].mpi for e in config.engines if e in KAPPA_ENGINES)
        )
    catalog = catalog.normalized()
    catalog = filter_catalog_patch(catalog, config.angular_patch)
    masked = catalog.mask is not None
    native_engines = list(config.engines)
    for engine in config.engines:
        if engine not in KAPPA_ENGINES:
            raise ValueError(f"unknown kappa engine: {engine}")
        spec = KAPPA_ENGINES[engine]
        # Validate requested products before constructing an object on any rank.
        engine_parameters(config, engine, masked)
        if masked and not spec.supports_mask:
            supported = ", ".join(
                name for name, spec in KAPPA_ENGINES.items() if spec.supports_mask
            )
            raise ValueError(
                f"{engine} does not document read-mask semantics; use "
                f"one of {supported} for a masked kappa map"
            )

    config.output_dir.mkdir(parents=True, exist_ok=True)
    balls = None
    if native_engines:
        from cyballs import cballs
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
            rank_timing = None
            if participates:
                engine_root = config.output_dir / engine
                engine_root.mkdir(parents=True, exist_ok=True)
                if comm.rank == 0:
                    print(
                        f"Running {engine} with {config.threads} thread(s)"
                        + (f" on {comm.size} MPI rank(s)" if spec.mpi else ""),
                        flush=True,
                    )
                try:
                    setup_started = time.perf_counter()
                    setup_cpu_started = time.process_time()
                    balls.set(engine_parameters(config, engine, masked))
                    # Complete parameter parsing, catalog publication, and
                    # thread setup before timing tree construction + search.
                    balls.Run(level=["SetNumberThreads"])
                    setup_seconds = time.perf_counter() - setup_started
                    setup_cpu_seconds = time.process_time() - setup_cpu_started
                    started = time.perf_counter()
                    started_cpu = time.process_time()
                    balls.Run(level=["MainLoop"])
                    compute_seconds = time.perf_counter() - started
                    compute_cpu_seconds = time.process_time() - started_cpu
                    rank_timing = _timing_metadata(
                        setup_seconds, setup_cpu_seconds,
                        compute_seconds, compute_cpu_seconds,
                        "cTreeBalls parameter/thread setup plus MainLoop; "
                        "in-memory catalog registration and cleanup excluded",
                    )
                    rank_timing.update(rank=comm.rank,
                                       native_reported_cpu_time=float(balls.getCPUTime()))
                    if comm.rank == 0:
                        result = copy_engine_results(balls, engine, config)
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
            if not errors:
                rank_timings = comm.allgather(rank_timing) if spec.mpi else [rank_timing]
                if comm.rank == 0 and result is not None:
                    result.update(aggregate_rank_timings(rank_timings))
                    result["wall_time"] = result["compute_wall_time"]
                    result["cpu_time"] = result["native_reported_cpu_time"]
                    result["threads_per_rank"] = config.threads
                    result["parameters"] = engine_parameters(config, engine, masked)
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
                print(
                    f"Finished {engine}: compute wall "
                    f"{result['compute_wall_time']:.6g} s, summed process CPU "
                    f"{result['compute_cpu_time']:.6g} s",
                    flush=True,
                )
            comm.barrier()
    finally:
        if balls is not None:
            balls.clear_catalogs()

    if comm.rank == 0:
        timings = timing_summary(results)
        summary = {
            "catalog": {**catalog.metadata, "nbody": catalog.nbody},
            "requested_statistics": statistics_from_options(config.options),
            "engines": {
                engine: {
                    "cpu_time": value.get("cpu_time"),
                    "wall_time": value.get("wall_time"),
                    "setup_wall_time": value.get("setup_wall_time"),
                    "setup_cpu_time": value.get("setup_cpu_time"),
                    "compute_wall_time": value.get("compute_wall_time"),
                    "compute_cpu_time": value.get("compute_cpu_time"),
                    "total_wall_time": value.get("total_wall_time"),
                    "total_cpu_time": value.get("total_cpu_time"),
                    "native_reported_cpu_time": value.get("native_reported_cpu_time"),
                    "timing_scope": value.get("timing_scope"),
                    "rank_timings": value["rank_timings"],
                    "ranks": value["ranks"],
                    "parameters": value["parameters"],
                    "edge_corrections": value["edge_corrections"],
                    "result_type": value["result_type"],
                    "smooth_pivot": value["smooth_pivot"],
                    "warnings": value.get("warnings", []),
                }
                for engine, value in results.items()
            },
            "failures": failures,
            "catalog_reads": 1,
            "set_catalog_calls_per_process": 1 if native_engines else 0,
            "mpi_ranks": comm.size,
            "threads_per_rank": config.threads,
            "comparisons": compare_result_arrays(results),
            "timings": timings,
        }
        summary["plots"] = make_plots(results, config) if config.plots else []
        (config.output_dir / "summary.json").write_text(
            json.dumps(summary, indent=2) + "\n", encoding="utf-8"
        )
        write_timing_report(config.output_dir / "timing_report.txt", timings)
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


def parse_arguments(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    sources = parser.add_mutually_exclusive_group()
    sources.add_argument("--fits", type=Path, help="HEALPix convergence FITS map")
    sources.add_argument("--catalog-npz", type=Path, help="positions/kappa NPZ catalog")
    sources.add_argument("--synthetic-nside", type=int, help="generated smoke-test map")
    parser.add_argument("--field", type=int, default=0)
    parser.add_argument("--mask", type=Path, help="optional HEALPix 0/1 mask")
    parser.add_argument("--mask-threshold", type=float, default=0.0)
    parser.add_argument("--nside-down", type=int, default=0)
    parser.add_argument(
        "--max-points", type=int, default=0,
        help="deterministically retain at most this many valid pixels; 0 keeps all",
    )
    parser.add_argument(
        "--sampling-seed", type=int, default=8675309,
        help="seed for reproducible HEALPix bottom-hash thinning",
    )
    parser.add_argument(
        "--fits-chunk-pixels", type=int, default=1 << 20,
        help="pixels per memory-mapped FITS scan chunk",
    )
    parser.add_argument("--no-center-field", action="store_true")
    parser.add_argument("--save-catalog-npz", type=Path)
    parser.add_argument(
        "--engine", "--engines", dest="engines", action="append", default=[],
        help="repeat or use comma lists; all, all-omp, and all-mpi select available "
             "methods compatible with the requested edge correction and mask",
    )
    parser.add_argument("--list-engines", action="store_true")
    parser.add_argument(
        "--cballs", type=Path, default=DEFAULT_CBALLS,
        help=(
            "cballs executable; default lookup uses CTREEBALLS_CBALLS/CBALLS, "
            "the source tree, then PATH"
        ),
    )
    parser.add_argument("--outdir", type=Path, default=Path("Output_all_engines"))
    parser.add_argument(
        "--statistics", choices=("2pcf", "3pcf", "both"), default="both",
        help="requested work; maps to only-2pcf/only-3pcf for native engines",
    )
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
    parser.add_argument(
        "--dual-node-bin-slop", action="store_true",
        help="enable theta-sized approximate dual-node acceptance",
    )
    parser.add_argument("--nsmooth", type=int, default=16)
    parser.add_argument("--linear-bins", action="store_true")
    parser.add_argument(
        "--patch", action="store_true",
        help="select the same angular patch for every engine before thinning and centering; "
             "also enabled by --more-options patch",
    )
    parser.add_argument("--thetaL", type=float, default=0.0,
                        help="patch lower colatitude in degrees (0=north pole)")
    parser.add_argument("--thetaR", type=float, default=90.0,
                        help="patch upper colatitude in degrees (180=south pole)")
    parser.add_argument("--phiL", type=float, default=0.0,
                        help="patch lower longitude in degrees [0, 360]")
    parser.add_argument("--phiR", type=float, default=90.0,
                        help="patch upper longitude in degrees [0, 360]; no wraparound")
    parser.add_argument("--more-options", action="append", default=[])
    parser.add_argument(
        "--no-smooth-pivot", action="store_true",
        help="disable the SMOOTHPIVOTON default on engines that support it",
    )
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
    parser.add_argument("--no-plots", action="store_true")
    parser.add_argument(
        "--no-flatten-plots", action="store_true",
        help="omit Figure-7-style flattened 3PCF radial-bin plots",
    )
    return parser.parse_args(argv)


def main() -> int:
    args = parse_arguments()
    if args.mpi_ranks < 1:
        raise SystemExit("ERROR: --mpi-ranks must be positive")

    try:
        executable_methods = discover_search_methods(args.cballs)
        make_settings = discover_make_settings(args.cballs)
        smooth_setting = make_settings.get("SMOOTHPIVOTON")
        smooth_pivot_compiled = (
            smooth_setting == "1" if smooth_setting in {"0", "1"} else None
        )
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
            print_engine_table(available, smooth_pivot_compiled)
            return 0
        config = RunConfig(
            engines=args.engines or ("octree-2balls-omp",),
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
            patch=args.patch,
            options=merge_statistics_options(
                args.statistics,
                tuple(args.more_options) + (
                    ("no-smooth-pivot",) if args.no_smooth_pivot else ()
                ),
            ),
            result_type=args.result_type,
            verbose=args.verbose,
            verbose_log=args.verbose_log,
            continue_on_error=args.continue_on_error,
            edge_corrections=args.edge_corrections,
            dual_node_bin_slop=args.dual_node_bin_slop,
            smooth_pivot_compiled=smooth_pivot_compiled,
            plots=not args.no_plots,
            flatten_plots=not args.no_flatten_plots,
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
                    max_points=args.max_points,
                    sampling_seed=args.sampling_seed,
                    chunk_pixels=args.fits_chunk_pixels,
                    patch=config.angular_patch,
                )
            elif args.catalog_npz is not None:
                if args.mask is not None:
                    raise ValueError("--mask is only valid with --fits")
                catalog = catalog_from_npz(
                    args.catalog_npz, center_field=not args.no_center_field,
                    patch=config.angular_patch,
                )
            elif args.synthetic_nside is not None:
                if args.mask is not None:
                    raise ValueError("--mask is only valid with --fits")
                catalog = synthetic_healpix_catalog(
                    args.synthetic_nside, center_field=not args.no_center_field,
                    patch=config.angular_patch,
                )
            else:
                raise ValueError(
                    "select --fits, --catalog-npz, or --synthetic-nside"
                )
            if args.save_catalog_npz is not None:
                save_catalog_npz(args.save_catalog_npz, catalog)
            if config.angular_patch is not None:
                patch = config.angular_patch
                print(
                    f"Shared patch (degrees, open bounds): "
                    f"{patch.phi_left:g} < phi < {patch.phi_right:g}, "
                    f"{patch.theta_left:g} < colatitude < {patch.theta_right:g}; "
                    f"{catalog.metadata['eligible_pixels_before_patch']} eligible -> "
                    f"{catalog.metadata['eligible_pixels_before_thinning']} in patch -> "
                    f"{catalog.nbody} retained",
                    flush=True,
                )
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
            print(
                "Smooth-pivot modes: " + ", ".join(
                    f"{engine}={smooth_pivot_mode(config, engine)}"
                    for engine in config.engines
                ),
                flush=True,
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
