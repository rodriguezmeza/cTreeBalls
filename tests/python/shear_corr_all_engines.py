#!/usr/bin/env python3
"""Compare enabled full-sky shear dual-node engines on one retained catalog.

The input catalog is loaded once and passed directly to ``cyballs``. The
octree, median-KD-tree, and PCA-ball-tree engines compute transported spin-2
2PCF and accepted-node LogMultipole 3PCF products on spherical catalogs.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass, field
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time
from typing import Any, Iterable, Optional, Sequence

import numpy as np


# healpy imports matplotlib.  Keep its cache stable when a locked home
# directory would otherwise trigger a fresh font-cache build on every run.
if "MPLCONFIGDIR" not in os.environ:
    _matplotlib_home = Path.home() / ".matplotlib"
    if not _matplotlib_home.is_dir() or not os.access(_matplotlib_home, os.W_OK):
        _matplotlib_cache = Path(tempfile.gettempdir()) / f"ctreeballs-matplotlib-{os.getuid()}"
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
SHEAR_SPHERE_TWO_BALLS_ENGINE = "octree-shear-sphere-2balls-omp"
SHEAR_SPHERE_KDTREE_TWO_BALLS_ENGINE = "kdtree-shear-sphere-2balls-omp"
SHEAR_SPHERE_BALLTREE_TWO_BALLS_ENGINE = "balltree-shear-sphere-2balls-omp"
SHEAR_SPHERE_ENGINES = (
    SHEAR_SPHERE_TWO_BALLS_ENGINE,
    SHEAR_SPHERE_KDTREE_TWO_BALLS_ENGINE,
    SHEAR_SPHERE_BALLTREE_TWO_BALLS_ENGINE,
)
NATIVE_SHEAR_ENGINES = SHEAR_SPHERE_ENGINES
ENGINE_ORDER = NATIVE_SHEAR_ENGINES
SHEAR_ENGINE_SETTINGS = {
    SHEAR_SPHERE_TWO_BALLS_ENGINE: "OCTREESHEARSPHERE2BALLSOMPON",
    SHEAR_SPHERE_KDTREE_TWO_BALLS_ENGINE: "KDTREESHEARSPHERE2BALLSOMPON",
    SHEAR_SPHERE_BALLTREE_TWO_BALLS_ENGINE: "BALLTREESHEARSPHERE2BALLSOMPON",
}
COMPONENT_LABELS = ("Gamma0", "Gamma1", "Gamma2", "Gamma3")

# Prefer the extension built beside this source tree over a site installation.
if any(PROJECT_ROOT.glob("cyballs*.so")) or any(PROJECT_ROOT.glob("cyballs*.pyd")):
    _root = os.fspath(PROJECT_ROOT)
    sys.path[:] = [entry for entry in sys.path if entry != _root]
    sys.path.insert(0, _root)


@dataclass
class ShearCatalog:
    positions: np.ndarray
    gamma1: np.ndarray
    gamma2: np.ndarray
    weights: Optional[np.ndarray] = None
    geometry: str = "flat"
    metadata: dict[str, Any] = field(default_factory=dict)

    def normalized(self) -> "ShearCatalog":
        positions = np.ascontiguousarray(self.positions, dtype=np.float64)
        gamma1 = np.ascontiguousarray(self.gamma1, dtype=np.float64)
        gamma2 = np.ascontiguousarray(self.gamma2, dtype=np.float64)
        geometry = str(self.geometry).strip().lower()
        if geometry not in {"flat", "sphere"}:
            raise ValueError("geometry must be flat or sphere")
        if positions.ndim != 2 or positions.shape[1] not in (2, 3):
            raise ValueError(f"positions must have shape (N, 2) or (N, 3), got {positions.shape}")
        if positions.shape[1] == 2 and geometry == "sphere":
            raise ValueError("spherical shear positions require three coordinates")
        if positions.shape[1] == 2:
            positions = np.column_stack((positions, np.zeros(positions.shape[0])))
        count = positions.shape[0]
        if count < 3:
            raise ValueError("a shear catalog needs at least three points")
        for name, values in (("gamma1", gamma1), ("gamma2", gamma2)):
            if values.ndim != 1 or values.shape[0] != count:
                raise ValueError(f"{name} must have shape ({count},), got {values.shape}")
        if self.weights is None:
            weights = np.ones(count, dtype=np.float64)
        else:
            weights = np.ascontiguousarray(self.weights, dtype=np.float64)
            if weights.ndim != 1 or weights.shape[0] != count:
                raise ValueError(f"weights must have shape ({count},), got {weights.shape}")
        if not all(np.all(np.isfinite(values)) for values in
                   (positions, gamma1, gamma2, weights)):
            raise ValueError("positions, shear, and weights must be finite")
        if np.any(weights < 0.0):
            raise ValueError("weights must be non-negative")
        if geometry == "flat":
            if not np.allclose(positions[:, 2], positions[0, 2],
                               rtol=0.0, atol=1.0e-12):
                raise ValueError("flat-sky shear requires one tangent plane")
        else:
            norms = np.linalg.norm(positions, axis=1)
            if np.any(norms <= np.finfo(float).tiny):
                raise ValueError("spherical positions must be nonzero vectors")
            positions = np.ascontiguousarray(positions/norms[:, None])
        return ShearCatalog(
            positions=np.ascontiguousarray(positions),
            gamma1=gamma1,
            gamma2=gamma2,
            weights=weights,
            geometry=geometry,
            metadata=dict(self.metadata),
        )

    @property
    def nbody(self) -> int:
        return int(self.positions.shape[0])

    @property
    def gamma(self) -> np.ndarray:
        return self.gamma1 + 1j*self.gamma2


@dataclass
class RunConfig:
    statistics: str = "both"
    min_sep: float = 0.05
    max_sep: float = 2.0
    sep_units: str = "degree"
    bins: int = 6
    multipoles: int = 3
    phi_bins: int = 32
    threads: int = max(1, (os.cpu_count() or 2) - 1)
    tree_theta: float = 1.0
    smooth_radius: Optional[float] = None
    use_log_bins: bool = True
    options: Sequence[str] = field(default_factory=tuple)
    output_dir: Path = Path("Output_shear_all_engines")
    plots: bool = True
    plot_order: int = 0
    verbose: int = 1
    verbose_log: int = 0
    continue_on_error: bool = False

    def normalized(self) -> "RunConfig":
        if self.statistics not in {"2pcf", "3pcf", "both"}:
            raise ValueError("statistics must be 2pcf, 3pcf, or both")
        if self.sep_units not in {"degree", "arcmin", "radian", "projected"}:
            raise ValueError(
                "sep_units must be degree, arcmin, radian, or projected"
            )
        if not math.isfinite(self.min_sep) or not math.isfinite(self.max_sep):
            raise ValueError("separation limits must be finite")
        if self.min_sep <= 0.0 or self.max_sep <= self.min_sep:
            raise ValueError("require 0 < min_sep < max_sep")
        if self.sep_units == "degree" and self.max_sep >= 180.0:
            raise ValueError("degree separations must be below 180")
        if self.sep_units == "arcmin" and self.max_sep >= 180.0*60.0:
            raise ValueError("arcmin separations must be below 10800")
        if self.sep_units == "radian" and self.max_sep >= math.pi:
            raise ValueError("radian separations must be below pi")
        if self.bins < 4:
            raise ValueError("the active cTreeBalls parameter contract requires bins >= 4")
        if self.multipoles < 2:
            raise ValueError("cTreeBalls requires multipoles >= 2")
        if self.phi_bins < 4:
            raise ValueError("phi_bins must be at least 4")
        if self.threads < 1:
            raise ValueError("threads must be positive")
        if (self.smooth_radius is not None
                and (not math.isfinite(self.smooth_radius)
                     or self.smooth_radius < 0.0)):
            raise ValueError("smooth_radius must be finite and non-negative")
        if not -self.multipoles <= self.plot_order <= self.multipoles:
            raise ValueError("plot_order must lie on the requested multipole axis")
        return RunConfig(
            **{
                **self.__dict__,
                "options": tuple(split_options(self.options)),
                "output_dir": Path(self.output_dir).expanduser().resolve(),
            }
        )

    def native_limits(self, geometry: str = "flat") -> tuple[float, float]:
        if self.sep_units == "projected":
            if geometry == "sphere":
                raise ValueError("sep_units=projected is unavailable on the sphere")
            return self.min_sep, self.max_sep
        angular_scale = math.pi/180.0
        if self.sep_units == "arcmin":
            angular_scale /= 60.0
        minimum = self.min_sep*angular_scale \
            if self.sep_units in {"degree", "arcmin"} else self.min_sep
        maximum = self.max_sep*angular_scale \
            if self.sep_units in {"degree", "arcmin"} else self.max_sep
        transform = math.sin if geometry == "sphere" else math.tan
        return 2.0*transform(0.5*minimum), 2.0*transform(0.5*maximum)

    def projected_limits(self) -> tuple[float, float]:
        return self.native_limits("flat")


def split_options(values: Iterable[str] | str | None) -> list[str]:
    if values is None:
        return []
    if isinstance(values, str):
        values = (values,)
    result: list[str] = []
    for value in values:
        for item in str(value).split(","):
            item = item.strip()
            if item and item.lower() != "none" and item not in result:
                result.append(item)
    return result


def validate_spherical_smooth_radius(config: RunConfig,
                                     options: Sequence[str]) -> None:
    """Reject a group diameter that reaches the first measured separation."""
    if (config.smooth_radius is None or config.smooth_radius == 0.0
            or "no-smooth-pivot" in options):
        return
    minimum, _ = config.native_limits("sphere")
    radius = 2.0*math.sin(
        0.5*math.radians(config.smooth_radius/60.0)
    )
    if 2.0*radius > minimum*(1.0 + 32.0*np.finfo(float).eps):
        raise ValueError(
            "spherical smooth-pivot requires 2*rsmooth <= min-sep; "
            f"got rsmooth={config.smooth_radius:g} arcmin and "
            f"min-sep={config.min_sep:g} {config.sep_units}"
        )


def _field_selector(value: str | int) -> str | int:
    if isinstance(value, int):
        return value
    text = str(value).strip()
    return int(text) if text.lstrip("+-").isdigit() else text


def _healpix_layout(path: Path) -> tuple[int, bool, list[str]]:
    try:
        from astropy.io import fits
    except ImportError as exc:
        raise RuntimeError("HEALPix FITS input requires astropy") from exc
    path = Path(path).expanduser().resolve()
    with fits.open(path, memmap=True, lazy_load_hdus=True) as hdus:
        tables = [hdu for hdu in hdus if hasattr(hdu, "columns") and hdu.data is not None]
        if not tables:
            raise ValueError(f"no binary table was found in {path}")
        table = tables[0]
        nside = int(table.header.get("NSIDE", 0))
        ordering = str(table.header.get("ORDERING", "RING")).strip().upper()
        names = list(table.columns.names)
    if nside <= 0:
        raise ValueError(f"missing or invalid NSIDE in {path}")
    if ordering not in {"RING", "NESTED", "NEST"}:
        raise ValueError(f"unsupported HEALPix ORDERING={ordering!r} in {path}")
    return nside, ordering != "RING", names


def _read_healpix_pixels(path: Path, field: str | int,
                         pixels: np.ndarray) -> np.ndarray:
    """Read selected implicit HEALPix pixels without materializing the full map."""
    try:
        from astropy.io import fits
    except ImportError as exc:
        raise RuntimeError("HEALPix FITS input requires astropy") from exc
    selector = _field_selector(field)
    with fits.open(Path(path).expanduser().resolve(), memmap=True,
                   lazy_load_hdus=True) as hdus:
        tables = [hdu for hdu in hdus if hasattr(hdu, "columns") and hdu.data is not None]
        if not tables:
            raise ValueError(f"no binary table was found in {path}")
        table = tables[0]
        try:
            column = table.data.field(selector)
        except (KeyError, IndexError, TypeError, ValueError) as exc:
            raise ValueError(
                f"field {field!r} is unavailable in {path}; fields are {table.columns.names}"
            ) from exc
        flat = np.asarray(column).reshape(-1)
        if pixels.size and int(pixels[-1]) >= flat.size:
            raise ValueError(f"HEALPix field {field!r} is shorter than NSIDE declares")
        return np.ascontiguousarray(flat[pixels], dtype=np.float64)


def _masked_healpix_pixels(path: Path, field: str | int, npix: int,
                           threshold: float, max_points: int,
                           chunk_size: int = 1_000_000) -> np.ndarray:
    """Select deterministic valid-mask ranks without an all-sky index array."""
    try:
        from astropy.io import fits
    except ImportError as exc:
        raise RuntimeError("HEALPix FITS input requires astropy") from exc
    selector = _field_selector(field)
    with fits.open(Path(path).expanduser().resolve(), memmap=True,
                   lazy_load_hdus=True) as hdus:
        tables = [hdu for hdu in hdus
                  if hasattr(hdu, "columns") and hdu.data is not None]
        if not tables:
            raise ValueError(f"no binary table was found in {path}")
        flat = np.asarray(tables[0].data.field(selector)).reshape(-1)
        if flat.size < npix:
            raise ValueError("the HEALPix mask is shorter than NSIDE declares")

        valid_count = 0
        for start in range(0, npix, chunk_size):
            values = np.asarray(flat[start:min(start + chunk_size, npix)])
            valid_count += int(np.count_nonzero(
                _valid_healpix(values) & (values > threshold)
            ))
        if valid_count == 0:
            return np.empty(0, dtype=np.int64)
        wanted = None
        if max_points and valid_count > max_points:
            wanted = np.linspace(0, valid_count, max_points,
                                 endpoint=False, dtype=np.int64)
        pieces: list[np.ndarray] = []
        rank_start = 0
        for start in range(0, npix, chunk_size):
            values = np.asarray(flat[start:min(start + chunk_size, npix)])
            local = np.flatnonzero(
                _valid_healpix(values) & (values > threshold)
            ).astype(np.int64, copy=False)
            if wanted is None:
                if local.size:
                    pieces.append(local + start)
            elif local.size:
                rank_stop = rank_start + local.size
                selected = wanted[(wanted >= rank_start) & (wanted < rank_stop)]
                if selected.size:
                    pieces.append(local[selected - rank_start] + start)
                rank_start = rank_stop
        return (np.concatenate(pieces) if pieces
                else np.empty(0, dtype=np.int64))


def _valid_healpix(values: np.ndarray) -> np.ndarray:
    try:
        import healpy as hp
        unseen = hp.UNSEEN
    except ImportError:
        unseen = -1.6375e30
    return np.isfinite(values) & (values != unseen) & (np.abs(values) < 1.0e29)


def stereographic_shear_patch(
    vectors: np.ndarray,
    gamma1: np.ndarray,
    gamma2: np.ndarray,
    *,
    center_ra_deg: float,
    center_dec_deg: float,
    rotate_from_local: bool = True,
    conjugate_input: bool = False,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Project unit vectors and rotate local east/north spin-2 components."""
    vectors = np.asarray(vectors, dtype=np.float64)
    gamma = np.asarray(gamma1, dtype=np.float64) + 1j*np.asarray(gamma2, dtype=np.float64)
    if conjugate_input:
        gamma = np.conjugate(gamma)
    ra0 = math.radians(center_ra_deg)
    dec0 = math.radians(center_dec_deg)
    center = np.array((math.cos(dec0)*math.cos(ra0),
                       math.cos(dec0)*math.sin(ra0), math.sin(dec0)))
    east0 = np.array((-math.sin(ra0), math.cos(ra0), 0.0))
    north0 = np.array((-math.sin(dec0)*math.cos(ra0),
                       -math.sin(dec0)*math.sin(ra0), math.cos(dec0)))
    denominator = 1.0 + vectors @ center
    if np.any(denominator <= np.finfo(float).eps):
        raise ValueError("the patch reaches the stereographic antipode")
    east_coordinate = vectors @ east0
    north_coordinate = vectors @ north0
    x = 2.0*east_coordinate/denominator
    y = 2.0*north_coordinate/denominator

    if rotate_from_local:
        longitude = np.arctan2(vectors[:, 1], vectors[:, 0])
        local_east = np.column_stack((
            -np.sin(longitude), np.cos(longitude), np.zeros(longitude.size)
        ))
        tangent_center = local_east @ center
        dx = 2.0*((local_east @ east0)*denominator
                  - east_coordinate*tangent_center)/(denominator*denominator)
        dy = 2.0*((local_east @ north0)*denominator
                  - north_coordinate*tangent_center)/(denominator*denominator)
        orientation = np.arctan2(dy, dx)
        gamma = gamma*np.exp(2j*orientation)

    positions = np.ascontiguousarray(
        np.column_stack((x, y, np.zeros(x.size))), dtype=np.float64
    )
    return positions, np.ascontiguousarray(gamma.real), np.ascontiguousarray(gamma.imag)


def catalog_from_healpix_patch(
    fits_path: Path,
    *,
    gamma1_field: str | int = "GAMMA1",
    gamma2_field: str | int = "GAMMA2",
    mask_path: Optional[Path] = None,
    mask_field: str | int = 0,
    mask_threshold: float = 0.0,
    weight_field: Optional[str | int] = None,
    center_ra_deg: float = 0.0,
    center_dec_deg: float = 0.0,
    patch_radius_deg: float = 5.0,
    max_points: int = 512,
    input_shear_frame: str = "local-east-north",
    conjugate_input: bool = False,
) -> ShearCatalog:
    try:
        import healpy as hp
    except ImportError as exc:
        raise RuntimeError("HEALPix FITS input requires healpy") from exc
    fits_path = Path(fits_path).expanduser().resolve()
    nside, nested, fields = _healpix_layout(fits_path)
    if not 0.0 < patch_radius_deg < 90.0:
        raise ValueError("patch_radius_deg must lie between 0 and 90")
    if max_points < 0 or 0 < max_points < 3:
        raise ValueError("max_points must be zero or at least three")
    center_vector = hp.ang2vec(center_ra_deg, center_dec_deg, lonlat=True)
    pixels = np.asarray(
        hp.query_disc(nside, center_vector, math.radians(patch_radius_deg),
                      inclusive=False, nest=nested),
        dtype=np.int64,
    )
    pixels.sort()
    if pixels.size < 3:
        raise ValueError("the requested patch contains fewer than three pixels")

    if mask_path is not None:
        mask_path = Path(mask_path).expanduser().resolve()
        mask_nside, mask_nested, _ = _healpix_layout(mask_path)
        if mask_nside != nside or mask_nested != nested:
            raise ValueError("the mask must have the same NSIDE and ORDERING as the shear map")
        mask_values = _read_healpix_pixels(mask_path, mask_field, pixels)
        keep = _valid_healpix(mask_values) & (mask_values > mask_threshold)
        pixels = pixels[keep]
        if pixels.size < 3:
            raise ValueError("the mask excludes all but fewer than three patch pixels")

    available = int(pixels.size)
    if max_points and pixels.size > max_points:
        indices = np.linspace(0, pixels.size, max_points, endpoint=False, dtype=np.int64)
        pixels = pixels[indices]

    gamma1 = _read_healpix_pixels(fits_path, gamma1_field, pixels)
    gamma2 = _read_healpix_pixels(fits_path, gamma2_field, pixels)
    if weight_field is None:
        weights = np.ones(pixels.size, dtype=np.float64)
    else:
        weights = _read_healpix_pixels(fits_path, weight_field, pixels)
    valid = _valid_healpix(gamma1) & _valid_healpix(gamma2)
    valid &= np.isfinite(weights) & (weights >= 0.0)
    pixels, gamma1, gamma2, weights = (
        values[valid] for values in (pixels, gamma1, gamma2, weights)
    )
    if pixels.size < 3:
        raise ValueError("the selected patch has fewer than three valid shear pixels")
    x, y, z = hp.pix2vec(nside, pixels, nest=nested)
    vectors = np.column_stack((x, y, z))
    positions, gamma1, gamma2 = stereographic_shear_patch(
        vectors, gamma1, gamma2,
        center_ra_deg=center_ra_deg,
        center_dec_deg=center_dec_deg,
        rotate_from_local=input_shear_frame == "local-east-north",
        conjugate_input=conjugate_input,
    )
    return ShearCatalog(
        positions=positions,
        gamma1=gamma1,
        gamma2=gamma2,
        weights=weights,
        metadata={
            "source": os.fspath(fits_path),
            "source_fields": fields,
            "gamma1_field": gamma1_field,
            "gamma2_field": gamma2_field,
            "weight_field": weight_field,
            "mask_source": os.fspath(mask_path) if mask_path is not None else None,
            "mask_field": mask_field if mask_path is not None else None,
            "nside": nside,
            "ordering": "NESTED" if nested else "RING",
            "center_ra_deg": center_ra_deg,
            "center_dec_deg": center_dec_deg,
            "patch_radius_deg": patch_radius_deg,
            "patch_points_before_thinning": available,
            "max_points": max_points,
            "input_shear_frame": input_shear_frame,
            "conjugate_input_shear": conjugate_input,
            "projection": "stereographic",
        },
    ).normalized()


def catalog_from_healpix_sphere(
    fits_path: Path,
    *,
    gamma1_field: str | int = "GAMMA1",
    gamma2_field: str | int = "GAMMA2",
    mask_path: Optional[Path] = None,
    mask_field: str | int = 0,
    mask_threshold: float = 0.0,
    weight_field: Optional[str | int] = None,
    max_points: int = 512,
    input_shear_frame: str = "local-east-north",
    conjugate_input: bool = False,
) -> ShearCatalog:
    """Load a HEALPix spin-2 map without projecting away the sky geometry."""
    try:
        import healpy as hp
    except ImportError as exc:
        raise RuntimeError("HEALPix FITS input requires healpy") from exc
    if input_shear_frame != "local-east-north":
        raise ValueError(
            "full-sky input requires gamma1+i*gamma2 in the local east/north frame"
        )
    if max_points < 0 or 0 < max_points < 3:
        raise ValueError("max_points must be zero or at least three")
    fits_path = Path(fits_path).expanduser().resolve()
    nside, nested, fields = _healpix_layout(fits_path)
    npix = hp.nside2npix(nside)

    if mask_path is not None:
        mask_path = Path(mask_path).expanduser().resolve()
        mask_nside, mask_nested, _ = _healpix_layout(mask_path)
        if mask_nside != nside or mask_nested != nested:
            raise ValueError("the mask must have the same NSIDE and ORDERING as the shear map")
        pixels = _masked_healpix_pixels(
            mask_path, mask_field, npix, mask_threshold, max_points,
        )
        available = None if max_points else int(pixels.size)
    else:
        available = int(npix)
        if max_points and npix > max_points:
            pixels = np.linspace(0, npix, max_points,
                                 endpoint=False, dtype=np.int64)
        else:
            pixels = np.arange(npix, dtype=np.int64)
    if pixels.size < 3:
        raise ValueError("the full-sky selection contains fewer than three pixels")

    gamma1 = _read_healpix_pixels(fits_path, gamma1_field, pixels)
    gamma2 = _read_healpix_pixels(fits_path, gamma2_field, pixels)
    weights = (np.ones(pixels.size, dtype=np.float64)
               if weight_field is None else
               _read_healpix_pixels(fits_path, weight_field, pixels))
    valid = _valid_healpix(gamma1) & _valid_healpix(gamma2)
    valid &= np.isfinite(weights) & (weights >= 0.0)
    pixels, gamma1, gamma2, weights = (
        values[valid] for values in (pixels, gamma1, gamma2, weights)
    )
    if pixels.size < 3:
        raise ValueError("the full-sky selection has fewer than three valid shear pixels")
    if conjugate_input:
        gamma2 = -gamma2
    x, y, z = hp.pix2vec(nside, pixels, nest=nested)
    return ShearCatalog(
        positions=np.column_stack((x, y, z)),
        gamma1=gamma1,
        gamma2=gamma2,
        weights=weights,
        geometry="sphere",
        metadata={
            "source": os.fspath(fits_path),
            "source_fields": fields,
            "gamma1_field": gamma1_field,
            "gamma2_field": gamma2_field,
            "weight_field": weight_field,
            "mask_source": os.fspath(mask_path) if mask_path is not None else None,
            "mask_field": mask_field if mask_path is not None else None,
            "nside": nside,
            "ordering": "NESTED" if nested else "RING",
            "full_sky_points_before_thinning": available,
            "max_points": max_points,
            "input_shear_frame": input_shear_frame,
            "conjugate_input_shear": conjugate_input,
            "projection": "none-full-sky",
        },
    ).normalized()


def catalog_from_npz(path: Path, geometry: Optional[str] = None) -> ShearCatalog:
    path = Path(path).expanduser().resolve()
    with np.load(path, allow_pickle=False) as archive:
        stored_geometry = (
            str(np.asarray(archive["geometry"]).reshape(()).item())
            if "geometry" in archive else "flat"
        )
        catalog = ShearCatalog(
            positions=archive["positions"],
            gamma1=archive["gamma1"],
            gamma2=archive["gamma2"],
            weights=archive["weights"] if "weights" in archive else None,
            geometry=geometry or stored_geometry,
            metadata={"source": os.fspath(path), "projection": "already-flat"},
        ).normalized()
        if "mask" in archive:
            mask = np.asarray(archive["mask"])
            if mask.ndim != 1 or mask.size != catalog.nbody:
                raise ValueError(f"mask must have shape ({catalog.nbody},)")
            if not np.all((mask == 0) | (mask == 1)):
                raise ValueError("mask values must be boolean or 0/1")
            keep = mask.astype(bool)
            catalog = ShearCatalog(
                positions=catalog.positions[keep], gamma1=catalog.gamma1[keep],
                gamma2=catalog.gamma2[keep], weights=catalog.weights[keep],
                geometry=catalog.geometry,
                metadata={**catalog.metadata, "mask_preselected": True},
            ).normalized()
    return catalog


def synthetic_shear_catalog(nbody: int = 128, seed: int = 8675309) -> ShearCatalog:
    if nbody < 16:
        raise ValueError("synthetic_nbody must be at least 16")
    rng = np.random.default_rng(seed)
    positions_2d = rng.uniform(-0.46, 0.46, size=(nbody, 2))
    x, y = positions_2d.T
    angle = np.arctan2(y, x)
    radius = np.hypot(x, y)
    gamma = 0.11*np.exp(2j*angle)*np.exp(-0.7*radius)
    gamma += (0.025*x - 0.015*y) + 1j*(-0.012*x + 0.021*y)
    gamma += 0.008*(rng.normal(size=nbody) + 1j*rng.normal(size=nbody))
    return ShearCatalog(
        positions=np.column_stack((x, y, np.zeros(nbody))),
        gamma1=gamma.real,
        gamma2=gamma.imag,
        weights=0.75 + 0.5*rng.random(nbody),
        geometry="flat",
        metadata={"source": "synthetic", "seed": seed, "projection": "already-flat"},
    ).normalized()


def synthetic_spherical_shear_catalog(
    nbody: int = 128, seed: int = 8675309,
) -> ShearCatalog:
    if nbody < 16:
        raise ValueError("synthetic_nbody must be at least 16")
    rng = np.random.default_rng(seed)
    positions = rng.normal(size=(nbody, 3))
    positions /= np.linalg.norm(positions, axis=1)[:, None]
    longitude = np.arctan2(positions[:, 1], positions[:, 0])
    latitude = np.arcsin(positions[:, 2])
    gamma = (0.08 + 0.025*np.cos(3.0*latitude)) \
        * np.exp(2j*(longitude + 0.2*np.sin(latitude)))
    gamma += 0.006*(rng.normal(size=nbody) + 1j*rng.normal(size=nbody))
    return ShearCatalog(
        positions=positions,
        gamma1=gamma.real,
        gamma2=gamma.imag,
        weights=0.75 + 0.5*rng.random(nbody),
        geometry="sphere",
        metadata={"source": "synthetic", "seed": seed,
                  "projection": "none-full-sky"},
    ).normalized()


def save_catalog_npz(path: Path, catalog: ShearCatalog) -> None:
    path = Path(path).expanduser().resolve()
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        path, positions=catalog.positions, gamma1=catalog.gamma1,
        gamma2=catalog.gamma2, weights=catalog.weights,
        geometry=np.asarray(catalog.geometry),
    )



def available_engines(cballs_executable: Path = DEFAULT_CBALLS) -> list[str]:
    result: list[str] = []
    try:
        from cyballs import search_method_id
        for engine in NATIVE_SHEAR_ENGINES:
            if search_method_id(engine) >= 0:
                result.append(engine)
    except (ImportError, OSError, AttributeError):
        pass
    executable = Path(cballs_executable).expanduser().resolve()
    if any(engine in result for engine in NATIVE_SHEAR_ENGINES) \
            and executable.is_file():
        probe = subprocess.run(
            [os.fspath(executable), "options=print-search-methods"],
            cwd=executable.parent, text=True, stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT, check=False,
        )
        for engine in NATIVE_SHEAR_ENGINES:
            if engine in result and f"- {engine} " not in probe.stdout:
                result.remove(engine)
    return result


def resolve_engines(tokens: Sequence[str], available: Sequence[str],
                    geometry: str = "sphere") -> list[str]:
    if geometry != "sphere":
        raise ValueError(
            "the active build exposes only full-sky shear dual-node engines; "
            "use --geometry sphere"
        )
    requested = split_options(tokens or ("all",))
    native_engines = list(SHEAR_SPHERE_ENGINES)
    if any(name in {"all", "all-shear"} for name in requested):
        requested = [name for name in native_engines if name in available]
    elif "all-omp" in requested:
        requested = [name for name in native_engines if name in available]
    unknown = [name for name in requested if name not in ENGINE_ORDER]
    if unknown:
        raise ValueError("unknown shear engine(s): " + ", ".join(unknown))
    missing = [name for name in requested if name not in available]
    if missing:
        raise RuntimeError("requested engine(s) are unavailable: " + ", ".join(missing))
    return [name for name in ENGINE_ORDER if name in requested]



def solve_mode_coupling(upsilon: np.ndarray, window: np.ndarray,
                        max_n: int) -> np.ndarray:
    """Apply C[ell,n] = N[ell-n]/N[0] with the addon's zero policy."""
    components, bins, bins_2, multipoles = upsilon.shape
    if components != 4 or bins != bins_2 or multipoles != 2*max_n + 1:
        raise ValueError("inconsistent shear multipole dimensions")
    if window.shape != (bins, bins, 4*max_n + 1):
        raise ValueError("inconsistent shear window dimensions")
    orders = np.arange(-max_n, max_n + 1)
    result = np.zeros_like(upsilon)
    for bin_1 in range(bins):
        for bin_2 in range(bins):
            n_zero = window[bin_1, bin_2, 2*max_n]
            if abs(n_zero) <= np.finfo(float).tiny:
                continue
            matrix = np.empty((multipoles, multipoles), dtype=np.complex128)
            matrix_scale = 0.0
            for row, ell in enumerate(orders):
                for column, order in enumerate(orders):
                    value = window[bin_1, bin_2, ell - order + 2*max_n]
                    matrix[row, column] = value/n_zero
                    matrix_scale = max(matrix_scale, abs(value)**2)
            rhs = (upsilon[:, bin_1, bin_2, :]/n_zero).T
            tolerance = 128.0*np.finfo(float).eps*(
                1.0 + math.sqrt(matrix_scale/abs(n_zero)**2)
            )
            singular = False
            for column in range(multipoles):
                pivot = column + int(np.argmax(np.abs(matrix[column:, column])**2))
                if abs(matrix[pivot, column])**2 <= tolerance**2:
                    singular = True
                    break
                if pivot != column:
                    matrix[[column, pivot], :] = matrix[[pivot, column], :]
                    rhs[[column, pivot], :] = rhs[[pivot, column], :]
                divisor = matrix[column, column]
                matrix[column, :] /= divisor
                rhs[column, :] /= divisor
                for row in range(multipoles):
                    if row == column:
                        continue
                    factor = matrix[row, column]
                    if factor != 0.0:
                        matrix[row, :] -= factor*matrix[column, :]
                        rhs[row, :] -= factor*rhs[column, :]
            if not singular:
                result[:, bin_1, bin_2, :] = rhs.T
    return result


def _timing_metadata(setup_wall: float, setup_cpu: float,
                     compute_wall: float, compute_cpu: float,
                     scope: str) -> dict[str, Any]:
    return {
        "setup_wall_time": float(setup_wall),
        "setup_cpu_time": float(setup_cpu),
        "compute_wall_time": float(compute_wall),
        "compute_cpu_time": float(compute_cpu),
        "total_wall_time": float(setup_wall + compute_wall),
        "total_cpu_time": float(setup_cpu + compute_cpu),
        "timing_scope": scope,
    }



def run_ctreeballs(catalog: ShearCatalog, config: RunConfig,
                   engine: Optional[str] = None) -> dict[str, Any]:
    try:
        from cyballs import cballs
    except (ImportError, OSError) as exc:
        raise RuntimeError("cannot import cyballs; rebuild this source profile") from exc
    engine = engine or SHEAR_SPHERE_TWO_BALLS_ENGINE
    minimum, maximum = config.native_limits(catalog.geometry)
    span = np.ptp(catalog.positions, axis=0)
    length_box = max(2.0*maximum, 1.1*float(np.max(span)), 1.0)
    statistic_option = {
        "2pcf": "only-2pcf",
        "3pcf": "only-3pcf",
        "both": "",
    }[config.statistics]
    options = split_options(
        ("GGGCorrelation", "no-out-Hist", statistic_option, *config.options)
    )
    if catalog.geometry == "sphere":
        validate_spherical_smooth_radius(config, options)
    with tempfile.TemporaryDirectory(prefix="ctreeballs-shear-") as output:
        setup_started = time.perf_counter()
        setup_cpu_started = time.process_time()
        model = cballs()
        try:
            parameters = {
                "searchMethod": engine,
                "iCatalogs": "1",
                "usePeriodic": "false",
                "useLogHist": str(config.use_log_bins).lower(),
                "rminHist": minimum,
                "rangeN": maximum,
                "sizeHistN": config.bins,
                "sizeHistPhi": config.phi_bins,
                "mChebyshev": config.multipoles,
                "lengthBox": length_box,
                "numberThreads": config.threads,
                "theta": config.tree_theta,
                "verbose": config.verbose,
                "verbose_log": config.verbose_log,
                "rootDir": output,
                "options": ",".join(options),
            }
            if config.smooth_radius is not None:
                parameters["rsmooth"] = config.smooth_radius
            model.set(parameters)
            model.set_catalog(
                catalog.positions, gamma1=catalog.gamma1,
                gamma2=catalog.gamma2, weights=catalog.weights,
            )
            model.Run(level=["SetNumberThreads"])
            setup_wall = time.perf_counter() - setup_started
            setup_cpu = time.process_time() - setup_cpu_started
            started = time.perf_counter()
            started_cpu = time.process_time()
            model.Run(level=["MainLoop"])
            elapsed = time.perf_counter() - started
            elapsed_cpu = time.process_time() - started_cpu
            result: dict[str, Any] = {
                "engine": engine,
                "geometry": catalog.geometry,
                f"elapsed_{config.statistics}": elapsed,
                f"cpu_{config.statistics}": elapsed_cpu,
                "native_reported_cpu_time": float(model.getCPUTime()),
                "radius": model.getrBins().copy(),
            }
            result.update(_timing_metadata(
                setup_wall, setup_cpu, elapsed, elapsed_cpu,
                "cTreeBalls object/catalog/thread setup plus one native MainLoop; "
                f"{engine} native {config.statistics} execution path",
            ))
            if config.statistics in {"2pcf", "both"}:
                result.update(
                    xi_plus=model.getShearXiPlus().copy(),
                    xi_minus=model.getShearXiMinus().copy(),
                    pair_weight=model.getShearXiWeight().copy(),
                )
            if config.statistics in {"3pcf", "both"}:
                result.update(
                    orders=model.getShearMultipoleOrders().copy(),
                    window_orders=np.arange(-2*config.multipoles,
                                            2*config.multipoles + 1,
                                            dtype=np.int32),
                    upsilon=np.ascontiguousarray(
                        model.getShearUpsilonXMultipoles().transpose(0, 2, 3, 1)
                    ),
                    gamma=np.ascontiguousarray(
                        model.getShearGammaXMultipoles().transpose(0, 2, 3, 1)
                    ),
                    window=np.ascontiguousarray(
                        model.getShearWindowMultipoles().transpose(1, 2, 0)
                    ),
                )
            return result
        finally:
            model.struct_cleanup()


def symmetric_relative_difference(candidate: np.ndarray, reference: np.ndarray,
                                  floor: float = 1.0e-10) -> np.ndarray:
    candidate = np.asarray(candidate)
    reference = np.asarray(reference)
    denominator = np.abs(candidate) + np.abs(reference)
    scale = max(float(np.max(np.abs(candidate), initial=0.0)),
                float(np.max(np.abs(reference), initial=0.0)),
                np.finfo(float).tiny)
    result = np.full(np.broadcast_shapes(candidate.shape, reference.shape), np.nan)
    valid = denominator > floor*scale
    result[valid] = 2.0*np.abs(candidate[valid] - reference[valid])/denominator[valid]
    return result


def relative_summary(values: np.ndarray) -> dict[str, float | int | None]:
    finite = np.abs(np.asarray(values, dtype=np.float64))
    finite = finite[np.isfinite(finite)]
    if not finite.size:
        return {"valid_bins": 0, "median": None, "p95": None, "max": None}
    return {
        "valid_bins": int(finite.size),
        "median": float(np.median(finite)),
        "p95": float(np.percentile(finite, 95.0)),
        "max": float(np.max(finite)),
    }


def compare_results(results: dict[str, dict[str, Any]]) -> dict[str, Any]:
    names = [name for name in NATIVE_SHEAR_ENGINES if name in results]
    if len(names) < 2:
        return {}
    reference_name = names[0]
    reference = results[reference_name]
    comparison: dict[str, Any] = {"reference": reference_name, "engines": {}}
    for candidate_name in names[1:]:
        candidate = results[candidate_name]
        if not np.allclose(candidate["radius"], reference["radius"],
                           rtol=5.0e-13, atol=5.0e-15):
            raise RuntimeError(
                f"{candidate_name} and {reference_name} returned different radial centers"
            )
        entry: dict[str, Any] = {}
        if "xi_plus" in candidate and "xi_plus" in reference:
            entry["2pcf"] = {
                name: relative_summary(
                    symmetric_relative_difference(candidate[name], reference[name])
                )
                for name in ("xi_plus", "xi_minus", "pair_weight")
            }
        if "gamma" in candidate and "gamma" in reference:
            if not np.array_equal(candidate["orders"], reference["orders"]):
                raise RuntimeError(
                    f"{candidate_name} and {reference_name} returned different orders"
                )
            entry["3pcf"] = {}
            for stage in ("upsilon", "gamma", "window"):
                relative = symmetric_relative_difference(
                    candidate[stage], reference[stage]
                )
                stage_summary: dict[str, Any] = {"all": relative_summary(relative)}
                if stage != "window":
                    stage_summary["components"] = {
                        label: relative_summary(relative[index])
                        for index, label in enumerate(COMPONENT_LABELS)
                    }
                entry["3pcf"][stage] = stage_summary
        comparison["engines"][candidate_name] = entry
    return comparison


def _save_result(path: Path, result: dict[str, Any]) -> None:
    arrays = {name: value for name, value in result.items()
              if isinstance(value, np.ndarray)}
    np.savez_compressed(path, **arrays)


def timing_summary(results: dict[str, dict[str, Any]]) -> dict[str, dict[str, Any]]:
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
    columns = (
        ("engine", 25), ("backend", 10), ("setup_wall_s", 14),
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
            engine, values["backend"],
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
        "The native shear addon computes 2PCF and 3PCF together in one MainLoop.",
    ))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _display_radius(radius: np.ndarray, config: RunConfig,
                    geometry: str = "flat") -> np.ndarray:
    if config.sep_units == "projected":
        return radius
    half_radius = np.asarray(radius)/2.0
    radians = (2.0*np.arcsin(np.clip(half_radius, -1.0, 1.0))
               if geometry == "sphere" else 2.0*np.arctan(half_radius))
    if config.sep_units == "degree":
        return np.degrees(radians)
    if config.sep_units == "arcmin":
        return 60.0*np.degrees(radians)
    return radians


def make_plots(results: dict[str, dict[str, Any]], config: RunConfig) -> list[str]:
    import matplotlib.pyplot as plt

    paths: list[str] = []
    colors = {
        SHEAR_SPHERE_TWO_BALLS_ENGINE: "C1",
        SHEAR_SPHERE_KDTREE_TWO_BALLS_ENGINE: "C4",
        SHEAR_SPHERE_BALLTREE_TWO_BALLS_ENGINE: "C5",
    }
    if any("xi_plus" in result for result in results.values()):
        figure, axes = plt.subplots(2, 2, figsize=(11.0, 7.8), sharex=True)
        for engine, result in results.items():
            if "xi_plus" not in result:
                continue
            radius = _display_radius(
                result["radius"], config, str(result.get("geometry", "flat")),
            )
            for column, key in enumerate(("xi_plus", "xi_minus")):
                axes[0, column].plot(radius, result[key].real, marker="o",
                                     color=colors[engine], label=engine)
                axes[1, column].plot(radius, result[key].imag, marker="o",
                                     color=colors[engine], label=engine)
                axes[0, column].set_title(r"$\xi_+$" if column == 0 else r"$\xi_-$")
        axes[0, 0].set_ylabel("Real")
        axes[1, 0].set_ylabel("Imaginary")
        label = "Projected separation" if config.sep_units == "projected" \
            else f"Tangent-plane equivalent separation [{config.sep_units}]"
        for axis in axes[1, :]:
            axis.set_xlabel(label)
        for axis in axes.flat:
            if config.use_log_bins:
                axis.set_xscale("log")
            axis.grid(True, which="both", ls=":", alpha=0.45)
        axes[0, 0].legend(frameon=False)
        figure.suptitle("Weak-lensing shear 2PCF")
        figure.tight_layout()
        path = config.output_dir / "shear_2pcf_comparison.png"
        figure.savefig(path, dpi=170)
        plt.close(figure)
        paths.append(path.name)

    if any("gamma" in result for result in results.values()):
        order_index = config.plot_order + config.multipoles
        figure, axes = plt.subplots(4, 2, figsize=(12.0, 12.0), sharex=True)
        for component, component_label in enumerate(COMPONENT_LABELS):
            for engine, result in results.items():
                if "gamma" not in result:
                    continue
                values = result["gamma"][component, :, :, order_index].reshape(-1)
                x = np.arange(values.size)
                axes[component, 0].plot(x, values.real, color=colors[engine],
                                        lw=1.2, label=engine)
                axes[component, 1].plot(x, values.imag, color=colors[engine],
                                        lw=1.2, label=engine)
            axes[component, 0].set_ylabel(component_label)
            axes[component, 0].grid(True, ls=":", alpha=0.4)
            axes[component, 1].grid(True, ls=":", alpha=0.4)
        axes[0, 0].set_title("Real")
        axes[0, 1].set_title("Imaginary")
        axes[-1, 0].set_xlabel("Flattened (r1, r2) bin")
        axes[-1, 1].set_xlabel("Flattened (r1, r2) bin")
        axes[0, 0].legend(frameon=False)
        figure.suptitle(
            f"Window-corrected Porth-x shear 3PCF, order n={config.plot_order}"
        )
        figure.tight_layout()
        path = config.output_dir / f"shear_3pcf_gamma_n{config.plot_order:+d}.png"
        figure.savefig(path, dpi=170)
        plt.close(figure)
        paths.append(path.name)

        comparable = [
            name for name in NATIVE_SHEAR_ENGINES
            if name in results and "gamma" in results[name]
        ]
        if len(comparable) > 1:
            reference_name, candidate_name = comparable[:2]
            relative = symmetric_relative_difference(
                results[candidate_name]["gamma"],
                results[reference_name]["gamma"],
            )[..., order_index]
            figure, axes = plt.subplots(1, 4, figsize=(15.0, 3.8), constrained_layout=True)
            for component, component_label in enumerate(COMPONENT_LABELS):
                image = axes[component].imshow(
                    100.0*relative[component], origin="lower", cmap="magma",
                    vmin=0.0,
                )
                axes[component].set_title(component_label)
                axes[component].set_xlabel("r2 bin")
                axes[component].set_ylabel("r1 bin")
                figure.colorbar(image, ax=axes[component], label="Relative difference [%]")
            path = config.output_dir / (
                f"shear_3pcf_relative_{candidate_name}_vs_{reference_name}"
                f"_n{config.plot_order:+d}.png"
            )
            figure.savefig(path, dpi=170)
            plt.close(figure)
            paths.append(path.name)
    return paths


def run_engine_suite(catalog: ShearCatalog, engines: Sequence[str],
                     config: RunConfig) -> dict[str, dict[str, Any]]:
    catalog = catalog.normalized()
    config = config.normalized()
    config.output_dir.mkdir(parents=True, exist_ok=True)
    results: dict[str, dict[str, Any]] = {}
    for engine in engines:
        print(f"[run] {engine}: N={catalog.nbody}, threads={config.threads}")
        try:
            result = run_ctreeballs(catalog, config, engine)
        except Exception as exc:
            if not config.continue_on_error:
                raise
            print(f"[failed] {engine}: {exc}", file=sys.stderr)
            continue
        results[engine] = result
        _save_result(config.output_dir / f"{engine}_histograms.npz", result)
        print(
            f"[done] {engine}: compute wall {result['compute_wall_time']:.6g} s, "
            f"process CPU {result['compute_cpu_time']:.6g} s"
        )
    timings = timing_summary(results)
    summary = {
        "catalog": {**catalog.metadata, "nbody": catalog.nbody,
                    "geometry": catalog.geometry},
        "config": {
            "statistics": config.statistics,
            "min_sep": config.min_sep,
            "max_sep": config.max_sep,
            "sep_units": config.sep_units,
            "native_limits": config.native_limits(catalog.geometry),
            "bins": config.bins,
            "multipoles": config.multipoles,
            "phi_bins": config.phi_bins,
            "threads": config.threads,
            "tree_theta": config.tree_theta,
            "smooth_radius": config.smooth_radius,
            "use_log_bins": config.use_log_bins,
            "options": list(config.options),
        },
        "engines": {
            name: {key: value for key, value in result.items()
                   if not isinstance(value, np.ndarray)}
            for name, result in results.items()
        },
        "comparison": compare_results(results),
        "timings": timings,
        "conventions": {
            "input": (
                "gamma=gamma1+i*gamma2 in each local east/north frame"
                if catalog.geometry == "sphere" else
                "gamma=gamma1+i*gamma2 in the common projected x/y frame"
            ),
            "projection": catalog.metadata.get("projection", "unknown"),
            "component_order": "[Gamma0,Gamma1,Gamma2,Gamma3]",
            "window_correction": "C[ell,n]=N[ell-n]/N[0]",
        },
    }
    summary["plots"] = make_plots(results, config) if config.plots else []
    (config.output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    write_timing_report(config.output_dir / "timing_report.txt", timings)
    return results


def parse_arguments(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    sources = parser.add_mutually_exclusive_group()
    sources.add_argument("--fits", type=Path, help="two-field HEALPix shear FITS map")
    sources.add_argument("--catalog-npz", type=Path,
                         help="positions/gamma1/gamma2 NPZ catalog")
    sources.add_argument("--synthetic-nbody", type=int, default=None)
    parser.add_argument(
        "--geometry", choices=("sphere",), default="sphere",
        help="retain full-sky unit vectors and use great-circle spin transport",
    )
    parser.add_argument("--gamma1-field", default="GAMMA1")
    parser.add_argument("--gamma2-field", default="GAMMA2")
    parser.add_argument("--weight-field", default=None,
                        help="optional weight column in the shear FITS file")
    parser.add_argument("--mask", type=Path, help="optional matching HEALPix mask")
    parser.add_argument("--mask-field", default="0")
    parser.add_argument("--mask-threshold", type=float, default=0.0)
    parser.add_argument("--center-ra-deg", type=float, default=0.0)
    parser.add_argument("--center-dec-deg", type=float, default=0.0)
    parser.add_argument("--patch-radius-deg", type=float, default=5.0)
    parser.add_argument(
        "--max-points", type=int, default=512,
        help="deterministically thin the selected sky to this size; 0 keeps every pixel",
    )
    parser.add_argument(
        "--input-shear-frame", choices=("local-east-north", "flat"),
        default="local-east-north",
        help="FITS shear frame; local fields are spin-2 rotated during projection",
    )
    parser.add_argument(
        "--conjugate-input-shear", action="store_true",
        help="flip the gamma2 sign before tangent-plane rotation",
    )
    parser.add_argument("--save-catalog-npz", type=Path)
    parser.add_argument(
        "--engine", "--engines", dest="engines", action="append", default=[],
        help="repeat or pass comma lists; all selects every enabled shear engine",
    )
    parser.add_argument("--list-engines", action="store_true")
    parser.add_argument(
        "--cballs", type=Path, default=DEFAULT_CBALLS,
        help=(
            "cballs executable; default lookup uses CTREEBALLS_CBALLS/CBALLS, "
            "the source tree, then PATH"
        ),
    )
    parser.add_argument("--statistics", choices=("2pcf", "3pcf", "both"),
                        default="both")
    parser.add_argument("--min-sep", type=float, default=0.05)
    parser.add_argument("--max-sep", type=float, default=2.0)
    parser.add_argument("--sep-units",
                        choices=("degree", "arcmin", "radian", "projected"),
                        default="degree")
    parser.add_argument("--nbins", type=int, default=6)
    parser.add_argument("--multipoles", type=int, default=3)
    parser.add_argument("--phi-bins", type=int, default=32)
    parser.add_argument("--threads", type=int,
                        default=max(1, (os.cpu_count() or 2) - 1))
    parser.add_argument(
        "--tree-theta", type=float, default=1.0,
        help="cTreeBalls accepted-node angular tolerance; use --exact-tree for body traversal",
    )
    parser.add_argument(
        "--exact-tree", action="store_true",
        help="add no-one-ball to native octree/KD-tree/ball-tree shear scans",
    )
    parser.add_argument(
        "--no-smooth-pivot", action="store_true",
        help="disable the default SMOOTHPIVOT behavior in every capable shear engine",
    )
    parser.add_argument(
        "--smooth-radius", type=float,
        help=(
            "set the cTreeBalls pivot-grouping radius; spherical values are "
            "arcmin and must satisfy 2*rsmooth <= min-sep"
        ),
    )
    parser.add_argument("--linear-bins", action="store_true")
    parser.add_argument(
        "--more-options", action="append", default=[],
        help="additional comma-separated cTreeBalls options; may be repeated",
    )
    parser.add_argument("--outdir", type=Path, default=Path("Output_shear_all_engines"))
    parser.add_argument("--plot-order", type=int, default=0)
    parser.add_argument("--no-plots", action="store_true")
    parser.add_argument("--verbose", type=int, default=1)
    parser.add_argument("--verbose-log", type=int, default=0)
    parser.add_argument("--continue-on-error", action="store_true")
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_arguments(argv)
    try:
        available = available_engines(args.cballs)
        if args.list_engines:
            print("Shear-compatible engines:")
            for name in ENGINE_ORDER:
                if name == SHEAR_SPHERE_TWO_BALLS_ENGINE:
                    note = "native full-sky spin-2 octree dual-node estimator"
                elif name == SHEAR_SPHERE_KDTREE_TWO_BALLS_ENGINE:
                    note = ("native full-sky spin-2 median KD tree with dual-node "
                            "2PCF and accepted-node Gamma 3PCF")
                else:
                    note = ("native full-sky spin-2 FCFC PCA ball tree with "
                            "dual-node 2PCF and accepted-node Gamma 3PCF")
                setting = SHEAR_ENGINE_SETTINGS[name]
                print(
                    f"- {name}: {'available' if name in available else 'unavailable'}; "
                    f"build={setting}; mask=preselect; edge=3PCF mode-coupling; "
                    f"smooth=supported; {note}"
                )
            return 0
        engines = resolve_engines(args.engines, available, args.geometry)
        native_options = list(args.more_options)
        if args.exact_tree:
            native_options.append("no-one-ball")
        if args.no_smooth_pivot:
            native_options.append("no-smooth-pivot")
        config = RunConfig(
            statistics=args.statistics,
            min_sep=args.min_sep,
            max_sep=args.max_sep,
            sep_units=args.sep_units,
            bins=args.nbins,
            multipoles=args.multipoles,
            phi_bins=args.phi_bins,
            threads=args.threads,
            tree_theta=args.tree_theta,
            smooth_radius=args.smooth_radius,
            use_log_bins=not args.linear_bins,
            options=tuple(native_options),
            output_dir=args.outdir,
            plots=not args.no_plots,
            plot_order=args.plot_order,
            verbose=args.verbose,
            verbose_log=args.verbose_log,
            continue_on_error=args.continue_on_error,
        ).normalized()
        if args.fits is not None:
            input_keywords = {
                "gamma1_field": _field_selector(args.gamma1_field),
                "gamma2_field": _field_selector(args.gamma2_field),
                "mask_path": args.mask,
                "mask_field": _field_selector(args.mask_field),
                "mask_threshold": args.mask_threshold,
                "weight_field": (
                    _field_selector(args.weight_field)
                    if args.weight_field is not None else None
                ),
                "max_points": args.max_points,
                "input_shear_frame": args.input_shear_frame,
                "conjugate_input": args.conjugate_input_shear,
            }
            if args.geometry == "sphere":
                catalog = catalog_from_healpix_sphere(args.fits, **input_keywords)
            else:
                catalog = catalog_from_healpix_patch(
                    args.fits, **input_keywords,
                    center_ra_deg=args.center_ra_deg,
                    center_dec_deg=args.center_dec_deg,
                    patch_radius_deg=args.patch_radius_deg,
                )
        elif args.catalog_npz is not None:
            if args.mask is not None:
                raise ValueError("--mask applies to HEALPix --fits input; put mask in the NPZ")
            catalog = catalog_from_npz(args.catalog_npz, geometry=args.geometry)
        else:
            if args.mask is not None:
                raise ValueError("--mask requires --fits")
            catalog = (
                synthetic_spherical_shear_catalog(args.synthetic_nbody or 128)
                if args.geometry == "sphere" else
                synthetic_shear_catalog(args.synthetic_nbody or 128)
            )
        if args.save_catalog_npz is not None:
            save_catalog_npz(args.save_catalog_npz, catalog)
        run_engine_suite(catalog, engines, config)
        print(f"Results written to {config.output_dir}")
        return 0
    except (OSError, RuntimeError, ValueError) as exc:
        print(f"ERROR: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
