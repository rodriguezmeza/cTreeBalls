"""Inspect catalogue headers only; never load survey arrays or unpickle input."""
import ast
import json
from pathlib import Path
import struct
import zipfile


def inspect_catalog(filename):
    path = Path(filename).expanduser().resolve()
    if not path.is_file():
        raise ValueError(f"Catalogue does not exist: {path}")
    columns, rows, details = [], None, {}
    if path.suffix.lower() == ".npz":
        with zipfile.ZipFile(path) as archive:
            entries = archive.infolist()
            if len(entries) > 256:
                raise ValueError("Too many NPZ members for header inspection (maximum 256).")
            for entry in entries:
                if not entry.filename.endswith(".npy"):
                    continue
                with archive.open(entry) as stream:
                    magic = stream.read(8)
                    if len(magic) != 8 or magic[:6] != b"\x93NUMPY" or magic[6] not in (1, 2, 3):
                        raise ValueError("Invalid NPY member header.")
                    size = 2 if magic[6] == 1 else 4
                    count = struct.unpack("<H" if size == 2 else "<I", stream.read(size))[0]
                    if count > 65536:
                        raise ValueError("NPY header exceeds 64 KiB.")
                    try:
                        info = ast.literal_eval(stream.read(count).decode("utf-8" if magic[6] == 3 else "latin1"))
                    except (SyntaxError, ValueError) as exc:
                        raise ValueError("Invalid NPY header information.") from exc
                    if not isinstance(info, dict):
                        raise ValueError("Invalid NPY header information.")
                    shape = info.get("shape")
                    if not isinstance(shape, tuple) or not all(isinstance(n, int) and n >= 0 for n in shape):
                        raise ValueError("Invalid NPY array shape.")
                    name = Path(entry.filename).stem.lower()
                    columns.append(name)
                    details[name] = dict(shape=list(shape), dtype=info.get("descr"))
                    if name in {"positions", "kappa", "gamma1", "delta"} and shape:
                        if rows is not None and rows != shape[0]:
                            raise ValueError("Catalogue columns have inconsistent row counts.")
                        rows = shape[0]
    elif any(str(path).lower().endswith(s) for s in (".fits", ".fit", ".fits.gz", ".fz")):
        try:
            from astropy.io import fits
        except ImportError as exc:
            raise ValueError("FITS header inspection needs astropy in the tutor's Python; offline documentation help does not.") from exc
        with fits.open(path, memmap=True, lazy_load_hdus=True) as hdus:
            for i, hdu in enumerate(hdus):
                names = list(getattr(getattr(hdu, "columns", None), "names", []) or [])
                if names:
                    lowered = [n.lower() for n in names]
                    columns.extend(lowered)
                    details[str(i)] = dict(name=hdu.name, rows=hdu.header.get("NAXIS2"), columns=names)
                    rows = max(rows or 0, int(hdu.header.get("NAXIS2", 0)))
    else:
        raise ValueError("Use a FITS or NPZ catalogue. /run reads saved JSON run information.")
    names = set(columns)
    candidates = []
    if "kappa" in names:
        candidates.append("kappa")
    if {"gamma1", "gamma2"} <= names or {"g1", "g2"} <= names:
        candidates.append("shear")
    if names & {"delta", "delta_blind"} and names & {"forest_ids", "forest_id", "los_id", "thing_id"}:
        candidates.append("lya")
    return dict(path=str(path), rows=rows, columns=sorted(names), fields=candidates,
                geometry_hint="sphere" if {"ra", "dec"} <= names else None,
                details=details, scope="Headers only; coordinates, units, masks, weights and data quality have not been validated.")


def read_run(filename):
    path = Path(filename).expanduser().resolve()
    if path.is_dir():
        choices = [path/"approximation-penalty.json", path/"qualification/approximation-penalty.json",
                   path/"summary.json", path/"run_metadata.json"]
        path = next((p for p in choices if p.is_file()), path)
    if not path.is_file() or path.stat().st_size > 8_000_000:
        raise ValueError("Choose a saved JSON run or approximation-penalty file no larger than 8 MB.")
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError("Run information must be a JSON object.")
    return path, data
