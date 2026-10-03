"""Sparse DES/Takahashi shear tables and shared FITS-driver controls.

The DES workflow stores observer-frame x/y/z directions, not an implicit
HEALPix array. G2CONV describes the stored shear; WTSUM is provenance only.
Astropy is imported lazily, so synthetic and NPZ runs do not require it.
"""
from pathlib import Path
import re

import numpy as np


NAME = re.compile(r"DESY3_Takahashi_r(\d+)_bin([1-4])_region([1-4])\.fits$")
HEADER_KEYS = ("NSIDE", "ORDERING", "REALIZ", "REGION", "TOMOBIN", "NPIXCAT",
               "AREADEG", "KMEAN", "WTSUM", "NSRCPLN", "NZTYPE", "G2CONV",
               "G1MEAN", "G2MEAN")


def _fits():
    try:
        from astropy.io import fits
    except ImportError as exc:
        raise RuntimeError("FITS catalog input requires astropy") from exc
    return fits


def _table_candidates(hdus):
    for index, hdu in enumerate(hdus):
        columns = getattr(getattr(hdu, "columns", None), "names", None) or []
        names = {name.lower(): name for name in columns}
        if {"x", "y", "z"} <= names.keys():
            if len(names) != len(columns):
                raise ValueError("FITS table has ambiguous case-insensitive column names")
            yield index, hdu, names


def resolve_fits_format(path, requested="auto"):
    if requested not in ("auto", "healpix", "desy3"):
        raise ValueError("fits-format must be auto, healpix, or desy3")
    if requested != "auto":
        return requested
    with _fits().open(Path(path).expanduser(), memmap=True) as hdus:
        return "desy3" if next(_table_candidates(hdus), None) else "healpix"


def inspect_des_catalog(path):
    path = Path(path).expanduser().resolve()
    with _fits().open(path, memmap=True) as hdus:
        candidates = list(_table_candidates(hdus))
        if len(candidates) != 1:
            raise ValueError(f"{path}: expected exactly one x/y/z shear catalog table")
        index, hdu, names = candidates[0]
        required = {"x", "y", "z", "gamma1", "gamma2"}
        if not required <= names.keys():
            raise ValueError(f"{path}: DES shear table requires x/y/z/gamma1/gamma2")
        for name in required:
            column = hdu.columns[names[name]]
            if str(column.format) not in ("E", "D", "1E", "1D"):
                raise ValueError(f"{path}: {names[name]} must be a scalar floating-point column")
        header = {key: hdu.header[key] for key in HEADER_KEYS if key in hdu.header}
        rows = int(hdu.header["NAXIS2"])
    match = NAME.fullmatch(path.name)
    identity = dict(zip(("realization", "tomobin", "region"), map(int, match.groups()))) if match else {}
    for key, name in (("REALIZ", "realization"), ("TOMOBIN", "tomobin"), ("REGION", "region")):
        if key in header:
            value = int(header[key])
            if value != header[key] or value < 1 or (name != "realization" and value > 4):
                raise ValueError(f"{path}: invalid {key}")
            if name in identity and identity[name] != value:
                raise ValueError(f"{path}: filename/header {key} mismatch")
            identity[name] = value
    if header.get("NPIXCAT", rows) != rows:
        raise ValueError(f"{path}: NPIXCAT differs from row count")
    return dict(path=str(path), hdu=index, columns=names, header=header, rows=rows, **identity)


def shear_convention(header, requested="auto"):
    stored = str(header.get("G2CONV", "")).strip().upper()
    if requested == "auto":
        if stored == "T17RAW":
            requested = "takahashi"
        elif stored in ("LOCALEN", "EASTNORTH"):
            requested = "local-east-north"
        else:
            raise ValueError(f"unknown/missing G2CONV={stored!r}; specify --des-shear-convention takahashi or local-east-north")
    if requested not in ("takahashi", "local-east-north"):
        raise ValueError("unknown DES shear convention")
    return requested, requested == "takahashi"


def splitmix64(rows, seed):
    """Permutation of row IDs: bottom-k is independent of chunk boundaries."""
    with np.errstate(over="ignore"):
        value = rows.astype(np.uint64) + np.uint64(seed) + np.uint64(0x9E3779B97F4A7C15)
        value = (value ^ (value >> np.uint64(30))) * np.uint64(0xBF58476D1CE4E5B9)
        value = (value ^ (value >> np.uint64(27))) * np.uint64(0x94D049BB133111EB)
        return value ^ (value >> np.uint64(31))


def load_des_shear_catalog(path, *, max_points=0, seed=8675309,
                           chunk_rows=1 << 20, convention="auto"):
    """Load stored directions and shear with O(chunk_rows + selected rows) memory.

    Invalid directions/shears are removed before deterministic row selection.
    Unit pixel weights are retained; no centering, WTSUM division, dense-map
    reconstruction or additional footprint mask is applied.
    """
    if max_points < 0 or 0 < max_points < 3:
        raise ValueError("max_points must be zero or >=3")
    if not 0 <= seed < 2**64 or chunk_rows < 1:
        raise ValueError("sampling seed must be uint64; chunk_rows must be positive")
    info = inspect_des_catalog(path)
    resolved, flip = shear_convention(info["header"], convention)
    names, rows = info["columns"], info["rows"]
    keep = np.empty(0, np.int64)
    priority = np.empty(0, np.uint64)
    chunks, eligible = [], 0
    with _fits().open(info["path"], memmap=True) as hdus:
        table = hdus[info["hdu"]].data
        for start in range(0, rows, chunk_rows):
            stop = min(start + chunk_rows, rows)
            xyz = np.column_stack([table[names[a]][start:stop] for a in "xyz"]).astype(np.float64)
            norm = np.hypot(np.hypot(xyz[:, 0], xyz[:, 1]), xyz[:, 2])
            good = np.all(np.isfinite(xyz), axis=1) & np.isfinite(norm) & (norm > np.finfo(float).tiny)
            for field in ("gamma1", "gamma2"):
                good &= np.isfinite(table[names[field]][start:stop])
            selected = np.flatnonzero(good).astype(np.int64) + start
            eligible += len(selected)
            if max_points:
                combined = np.concatenate((keep, selected))
                hashed = np.concatenate((priority, splitmix64(selected, seed)))
                if len(combined) > max_points:
                    take = np.argpartition(hashed, max_points - 1)[:max_points]
                    combined, hashed = combined[take], hashed[take]
                keep, priority = combined, hashed
            else:
                chunks.append(selected)
        if not max_points:
            keep = np.concatenate(chunks) if chunks else np.empty(0, np.int64)
        keep.sort()
        if len(keep) < 3:
            raise ValueError("fewer than three valid selected shear rows")
        xyz = np.empty((len(keep), 3), np.float64)
        g1, g2 = np.empty(len(keep)), np.empty(len(keep))
        for start in range(0, len(keep), chunk_rows):
            take = keep[start:start + chunk_rows]
            sl = slice(start, start + len(take))
            for axis, name in enumerate("xyz"):
                xyz[sl, axis] = table[names[name]][take]
            g1[sl] = table[names["gamma1"]][take]
            g2[sl] = table[names["gamma2"]][take]
            block = xyz[sl]
            block /= np.hypot(np.hypot(block[:, 0], block[:, 1]), block[:, 2])[:, None]
    if flip:
        g2 *= -1
    metadata = dict(source=info["path"], loader="desy3-shear-table-memmap",
        fits_format="desy3", fits_hdu=info["hdu"], fits_header=info["header"],
        source_fields=list(names.values()), input_rows=rows, invalid_rows_removed=rows-eligible,
        eligible_pixels_before_thinning=eligible, footprint_preselected=True,
        statistical_weights="unit", source_plane_weights_renormalized=False,
        centered=False, projection="none-full-sky", input_shear_frame="local-east-north",
        stored_shear_convention=resolved, convention_request=convention,
        conjugate_input_shear=flip, shear_convention="gamma1+i*gamma2 in local east/north basis",
        max_points=max_points, sampling_seed=seed if max_points else None,
        sampling_method="splitmix64-bottom-k-row" if max_points else "none",
        fits_chunk_pixels=chunk_rows,
        **{k: info[k] for k in ("realization", "tomobin", "region") if k in info})
    return dict(positions=xyz, gamma1=g1, gamma2=g2, geometry="sphere", metadata=metadata)


def add_fits_table_arguments(parser):
    parser.add_argument("--fits-format", choices=("auto", "healpix", "desy3"), default="auto",
                        help="auto detects stored x/y/z tables; HEALPix uses implicit pixel directions")
    parser.add_argument("--des-shear-convention", choices=("auto", "takahashi", "local-east-north"), default="auto",
                        help="auto uses G2CONV; Takahashi gamma2 is negated exactly once")
    parser.add_argument("--sampling-seed", type=int, default=8675309,
                        help="uint64 seed for deterministic DES bottom-k row selection")
    parser.add_argument("--fits-chunk-rows", "--fits-chunk-pixels", dest="fits_chunk_rows", type=int, default=1 << 20)


def validate_des_controls(*, geometry="sphere", mask=None, weight_field=None,
                          input_shear_frame="local-east-north", conjugate_input=False,
                          gamma1_field="GAMMA1", gamma2_field="GAMMA2", nsides=()):
    if geometry != "sphere" or input_shear_frame != "local-east-north":
        raise ValueError("DES tables require spherical geometry and local-east-north shear")
    if mask is not None or weight_field is not None:
        raise ValueError("DES tables already encode the footprint and use unit weights; --mask/--weight-field are unsupported")
    if conjugate_input:
        raise ValueError("DES tables use --des-shear-convention; do not apply a second conjugation")
    if str(gamma1_field).lower() != "gamma1" or str(gamma2_field).lower() != "gamma2":
        raise ValueError("DES tables use named gamma1/gamma2 columns")
    if nsides:
        raise ValueError("--nsides requires a dense HEALPix map; DES tables use --sizes row samples")


def apply_binning_preset(name, minimum, maximum, bins, units, *, linear=False, geometry="sphere"):
    if name == "custom":
        return minimum, maximum, bins, units
    if name not in ("sofia-fig1", "paper-8-200-edges"):
        raise ValueError("unknown binning preset")
    if linear or geometry != "sphere":
        raise ValueError("named binning presets require spherical logarithmic bins")
    lo, hi = 2*np.sin(np.deg2rad(np.array([8., 200.])/60)/2)
    if name == "sofia-fig1":
        step = (hi/lo)**(1/19)
        lo, hi = lo/np.sqrt(step), hi*np.sqrt(step)
    limits = np.rad2deg(2*np.arcsin(np.array([lo, hi])/2))*60
    return float(limits[0]), float(limits[1]), 20, "arcmin"


def integer_selection(text, low=1, high=108):
    """Parse comma lists and inclusive ranges for DES batch identities."""
    values = set()
    for part in text.split(','):
        if ':' in part:
            a, b = map(int, part.split(':'))
            if a > b:
                raise ValueError('range endpoints must be increasing')
            if a < low or b > high:
                raise ValueError(f'selection must be within {low}..{high}')
            values.update(range(a, b+1))
        else:
            values.add(int(part))
    if not values or min(values) < low or max(values) > high:
        raise ValueError(f'selection must be within {low}..{high}')
    return sorted(values)


def discover_des_catalogs(root, realizations, tomobins, regions):
    found = {}
    for path in Path(root).expanduser().resolve().rglob('DESY3_Takahashi_r*_bin*_region*.fits'):
        match = NAME.fullmatch(path.name)
        if not match:
            continue
        key = tuple(map(int, match.groups()))
        if key[0] not in realizations or key[1] not in tomobins or key[2] not in regions:
            continue
        if key in found:
            raise ValueError(f'duplicate catalog identity {key}: {path}, {found[key]}')
        found[key] = path
    required = [(r,b,g) for r in realizations for b in tomobins for g in regions]
    missing = [key for key in required if key not in found]
    if missing:
        raise ValueError(f'missing requested (realization, bin, region): {missing}')
    return [inspect_des_catalog(found[key]) for key in required]
