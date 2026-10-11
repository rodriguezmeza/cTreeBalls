"""Parallel KD construction with odd populations, masks and multiple catalogs."""
import os
from unittest.mock import patch

import numpy as np

import test_shear_sphere_octree_omp as reference

ENGINE = "kdtree-shear-sphere-2balls-omp"
OPTIONS = "no-out-Hist,no-smooth-pivot"


def compare_builds(catalogs, **kwargs):
    # Keep the search frontier fixed so this compares serial/parallel builds
    # and scheduling without changing the histogram partition.
    with patch.dict(os.environ, {"CBALLS_SHEAR_TASKS": "128",
                                "CBALLS_SHEAR_PIVOT_TOL": "3",
                                "CBALLS_SHEAR_BIN_THETA": ".5"}):
        serial = reference.run_native_catalogs(catalogs, 1, engine=ENGINE,
                                               use_log=True, **kwargs)
        parallel = reference.run_native_catalogs(catalogs, 4, engine=ENGINE,
                                                 use_log=True, **kwargs)
    for name in serial:
        np.testing.assert_array_equal(parallel[name], serial[name], err_msg=name)
        assert np.isfinite(parallel[name]).all(), name


def test_parallel_layout_and_moments():
    catalog = reference.fixture(32771)
    catalog[2][::17] = 0.0
    mask = (np.arange(len(catalog[0])) % 11 != 0).astype(np.int32)
    assert mask.sum() >= 16384  # Masked tree must still launch build tasks.
    with patch.object(reference, "RMIN", .003), patch.object(reference, "RMAX", .04):
        for capacity in (1, 8, 31):
            for masked in (False, True):
                compare_builds([catalog], nsmooth=capacity,
                    masks=[mask] if masked else None,
                    options=OPTIONS + (",read-mask" if masked else ""))


def test_parallel_reuse_cross_catalogs():
    first = reference.fixture(32771)
    second = reference.fixture(32773)
    second[0][:] = np.roll(second[0], 1, axis=1)
    second[1][:] *= np.exp(.37j)
    # The catalog selector has one entry per supplied catalog; the third
    # catalog is unused when the first catalog fills the third role again.
    catalogs = [first, second, reference.fixture(37)]
    masks = [(np.arange(len(c[0])) % 13 != 0).astype(np.int32)
             for c in catalogs]
    with patch.object(reference, "RMIN", .003), patch.object(reference, "RMAX", .04):
        compare_builds(catalogs, masks=masks, nsmooth=8,
            catalog_order="1,2,1", options=OPTIONS+",read-mask,shear-pivot-reuse")


if __name__ == "__main__":
    for test in (test_parallel_layout_and_moments, test_parallel_reuse_cross_catalogs):
        test()
        print("PASS", test.__name__, flush=True)
