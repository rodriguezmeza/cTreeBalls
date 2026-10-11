"""Ball-tree construction determinism and great-circle pair regressions."""
import os
from unittest.mock import patch

import numpy as np

import test_shear_sphere_octree_omp as reference

ENGINE = "balltree-shear-sphere-2balls-omp"
OPTIONS = "no-out-Hist,no-smooth-pivot"


def test_parallel_layout_and_moments():
    catalog = reference.fixture(32771)
    mask = (np.arange(len(catalog[0])) % 11 != 0).astype(np.int32)
    if not (mask.sum() >= 16384):
        raise AssertionError('required condition failed: mask.sum() >= 16384')
    with patch.object(reference, "RMIN", .003), patch.object(reference, "RMAX", .04):
        for capacity in (1, 8, 31):
            for masked in (False, True):
                options = OPTIONS + (",read-mask" if masked else "")
                kwargs = dict(engine=ENGINE, use_log=True, nsmooth=capacity,
                              masks=[mask] if masked else None)
                serial = reference.run_native_catalogs([catalog], 4,
                    options=options+",no-balltree-parallel-build", **kwargs)
                parallel = reference.run_native_catalogs([catalog], 4,
                    options=options, **kwargs)
                for name in serial:
                    np.testing.assert_array_equal(parallel[name], serial[name],
                                                   err_msg=name)
                    if not (np.isfinite(parallel[name]).all()):
                        raise AssertionError(name)
                if capacity == 8:
                    single = reference.run_native_catalogs([catalog], 1,
                        options=options, **kwargs)
                    for name in serial:
                        np.testing.assert_array_equal(single[name], serial[name],
                                                       err_msg=name)
                    for control in ("no-balltree-shear-member-cache",
                                    "no-balltree-shear-projection-cache"):
                        uncached = reference.run_native_catalogs([catalog], 4,
                            options=options+","+control, **kwargs)
                        for name in serial:
                            np.testing.assert_array_equal(uncached[name], serial[name],
                                                           err_msg=name)


def test_parallel_selection_cross_catalogs():
    first = reference.fixture(32773)
    second = reference.fixture(20003)
    second[0][:] = np.roll(second[0], 1, axis=1)
    first[2][::17] = 0.
    catalogs = [first, second, reference.fixture(35)]
    masks = [(np.arange(len(c[0])) % 11 != 0).astype(np.int32) for c in catalogs]
    # Empty source chunks exercise the prefix offsets in stable compaction.
    masks[0][4096:8192] = 0
    masks[1][:2048] = 0
    with patch.dict(os.environ, CBALLS_SHEAR_TASKS='128',
                    CBALLS_SHEAR_PIVOT_TOL='3', CBALLS_SHEAR_BIN_THETA='.5'), \
         patch.object(reference, 'RMIN', .003), patch.object(reference, 'RMAX', .04):
        for extra in ('', ',no-one-ball', ',shear-pivot-reuse'):
            kwargs = dict(engine=ENGINE, use_log=True, nsmooth=8, masks=masks,
                          catalog_order='1,2,1')
            serial = reference.run_native_catalogs(catalogs, 1,
                options=OPTIONS+',read-mask,no-balltree-shear-projection-cache'+extra,
                **kwargs)
            parallel = reference.run_native_catalogs(catalogs, 4,
                options=OPTIONS+',read-mask'+extra, **kwargs)
            for name in serial:
                np.testing.assert_array_equal(parallel[name], serial[name], err_msg=name)
                if not np.isfinite(parallel[name]).all():
                    raise AssertionError(name)


def test_projection_ties_and_small_selection():
    catalog = reference.fixture(16385)
    # A ring gives repeated projections along its principal axes, including
    # points on the split planes. The address tie break must stay unchanged.
    angles = np.linspace(0., 2.*np.pi, len(catalog[0]), endpoint=False)
    catalog[0][:] = np.column_stack((np.cos(angles), np.sin(angles),
                                    np.zeros(len(angles))))
    mask = np.zeros(len(angles), dtype=np.int32)
    mask[::127] = 1  # 130 selected: just above the projection-cache threshold.
    with patch.object(reference, 'RMIN', .003), patch.object(reference, 'RMAX', .015):
        for masks in (None, [mask]):
            options = OPTIONS + (',read-mask' if masks else '')
            kwargs = dict(engine=ENGINE, use_log=False, nsmooth=8, masks=masks)
            serial = reference.run_native_catalogs([catalog], 1,
                options=options+',no-balltree-shear-projection-cache', **kwargs)
            parallel = reference.run_native_catalogs([catalog], 4,
                options=options, **kwargs)
            for name in serial:
                np.testing.assert_array_equal(parallel[name], serial[name], err_msg=name)
                if not np.isfinite(parallel[name]).all():
                    raise AssertionError(name)


def test_geodesic_small_separation():
    # Exact near-coincident and near-polar pairs with independently transported
    # shears; use cross catalogs to test the imaginary xi+ component as well.
    rng = np.random.default_rng(3901)
    first, gamma, weight = reference.fixture(36)
    first[:2] = [[0, 0, 1], [0, 0, -1]]
    second = first + rng.normal(size=first.shape)*3.e-5
    second /= np.linalg.norm(second, axis=1)[:, None]
    other_gamma = .02*(rng.normal(size=len(first)) + 1j*rng.normal(size=len(first)))
    other_weight = rng.uniform(.3, 2., len(first))
    with patch.object(reference, "RMIN", 1.e-7), patch.object(reference, "RMAX", 2.e-4):
        expected = {key: np.zeros(reference.BINS, dtype=complex)
                    for key in ["xi_plus", "xi_minus", "xi_weight"]}
        for i, pivot in enumerate(first):
            for j, neighbor in enumerate(second):
                distance = np.linalg.norm(neighbor-pivot)
                bin_index = reference.radial_bin(distance)
                if bin_index is None:
                    continue
                _, phase, transported = reference.geometry_and_transport(
                    pivot, neighbor, other_gamma[j])
                w = weight[i]*other_weight[j]
                expected["xi_plus"][bin_index] += w*gamma[i]*np.conj(transported)
                expected["xi_minus"][bin_index] += w*gamma[i]*transported*phase**-4
                expected["xi_weight"][bin_index] += w
        for key in ["xi_plus", "xi_minus"]:
            np.divide(expected[key], expected["xi_weight"], out=expected[key],
                      where=expected["xi_weight"] != 0)
        catalogs = [(first, gamma, weight), (second, other_gamma, other_weight)]
        for threads in (1, 4):
            actual = reference.run_native_catalogs(catalogs, threads, engine=ENGINE,
                options=OPTIONS+",only-2pcf,no-one-ball", nsmooth=8)
            for name in expected:
                np.testing.assert_allclose(actual[name], expected[name],
                                           rtol=2.e-9, atol=3.e-12, err_msg=name)


if __name__ == "__main__":
    for test in [test_parallel_layout_and_moments, test_parallel_selection_cross_catalogs,
                 test_projection_ties_and_small_selection, test_geodesic_small_separation]:
        test()
        print("PASS", test.__name__, flush=True)
