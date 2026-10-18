"""Bounded coordinate scratch and explicit native catalog registration.

Frontend preparation may materialize owning arrays from views; that cost belongs
to setup. The strict cyballs API itself never converts a transferred array.
"""
import argparse
import inspect

import numpy as np

DEFAULT_CHUNK_PIXELS = 1 << 20


def positive_chunk(value):
    value = int(value)
    if value < 1024:
        raise argparse.ArgumentTypeError("chunk size must be at least 1024 pixels")
    return value


def add_memory_arguments(parser, *, healpix=False):
    parser.add_argument(
        "--transfer-catalog-ownership", action="store_true",
        help="opt in to native registration without another NumPy catalog copy; "
             "buffers become read-only (views are materialized during setup); "
             "native particle storage is still allocated")
    if healpix:
        parser.add_argument(
            "--fits-chunk-pixels", type=positive_chunk, default=DEFAULT_CHUNK_PIXELS,
            help="HEALPix coordinate/gather chunk size (default 1048576); "
                 "bounds temporary memory without thinning or smoothing pixels")


def healpix_positions(hp, nside, pixels, *, nested=False,
                      chunk_pixels=DEFAULT_CHUNK_PIXELS):
    """Preserve pixel order and pix2vec arithmetic with bounded xyz scratch."""
    if chunk_pixels < 1:
        raise ValueError("chunk_pixels must be positive")
    positions = np.empty((len(pixels), 3), dtype=np.float64)
    for start in range(0, len(pixels), chunk_pixels):
        stop = min(start + chunk_pixels, len(pixels))
        x, y, z = hp.pix2vec(nside, pixels[start:stop], nest=nested)
        positions[start:stop, 0] = x
        positions[start:stop, 1] = y
        positions[start:stop, 2] = z
        del x, y, z
    return positions


def gather_field(values, pixels, chunk_pixels=DEFAULT_CHUNK_PIXELS):
    if chunk_pixels < 1:
        raise ValueError("chunk_pixels must be positive")
    result = np.empty(len(pixels), dtype=np.float64)
    for start in range(0, len(pixels), chunk_pixels):
        stop = min(start + chunk_pixels, len(pixels))
        result[start:stop] = values[pixels[start:stop]]
    return result


def require_transfer_support(method):
    if "take_ownership" not in inspect.signature(method).parameters:
        raise RuntimeError("--transfer-catalog-ownership requires rebuilding cyballs "
                           f"with {method.__name__}(..., take_ownership=True)")


def preflight_transfer(enabled, *, forest=False):
    if enabled:
        from cyballs import cballs
        require_transfer_support(cballs.set_forest_catalog if forest else cballs.set_catalog)


def _owned(values, dtype=np.float64):
    return None if values is None else np.require(values, dtype=dtype, requirements=["C", "A", "O"])


def register_catalog(model, transfer, positions, **fields):
    if not transfer:
        return model.set_catalog(positions, **fields)
    require_transfer_support(model.set_catalog)
    prepared = {key: (_owned(value, np.uint8 if key == "mask" else np.float64)
                      if key != "catalog" else value) for key, value in fields.items()}
    return model.set_catalog(_owned(positions), **prepared, take_ownership=True)


def register_forest_catalog(model, transfer, positions, delta, weights, forest_ids):
    if not transfer:
        return model.set_forest_catalog(positions, delta, weights, forest_ids)
    require_transfer_support(model.set_forest_catalog)
    return model.set_forest_catalog(_owned(positions), _owned(delta), _owned(weights),
                                    forest_ids, take_ownership=True)


def memory_metadata(transfer=False, chunk_pixels=None, *, native_input="numpy"):
    return {
        "native_catalog_registration": ("transfer" if transfer else "copy")
            if native_input == "numpy" else "file-input",
        "fits_chunk_pixels": chunk_pixels,
        "tree_node_allocation": "automatic; compact scalar octree/PCA/KD builders use required capacity",
        "native_particle_copy": native_input == "numpy",
        "scope": "Native backend only. Transfer preparation/registration is setup; "
                 "coordinate chunks do not change selection or smoothing. "
                 "Reserved node capacity is not measured peak RSS. "
                 "MPI retains a catalog/tree per rank. Peer backend ownership is unchanged.",
    }
