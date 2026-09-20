"""Regression/calibration fixtures for opt-in native octree pivot reuse."""
import os
import re
from unittest.mock import patch

import numpy as np

import test_shear_sphere_octree_omp as reference

ENGINE = "octree-shear-sphere-2balls-omp"
BASE = "no-out-Hist,no-smooth-pivot,only-3pcf"
REUSE = BASE + ",shear-pivot-reuse"


def compare(a, b, rtol=2.e-11, atol=3.e-12):
    for key in a:
        np.testing.assert_allclose(b[key], a[key], rtol=rtol, atol=atol,
                                   err_msg=key)


def test_exact_limit_and_fallbacks():
    catalog = reference.fixture(80)
    with patch.dict(os.environ, {"CBALLS_SHEAR_PIVOT_TOL": "1e-12"}):
        for use_log in (False, True):
            exact = reference.run_native(*catalog, 1, options=BASE+",no-one-ball",
                                         use_log=use_log, engine=ENGINE)
            actual = reference.run_native(*catalog, 4, options=REUSE,
                                          use_log=use_log, engine=ENGINE)
            compare(exact, actual)
            direct = reference.oracle(*catalog, use_log=use_log)
            compare({key: direct[key] for key in actual}, actual,
                    rtol=1.e-10, atol=2.e-11)
    for token in ("no-one-ball", "no-two-balls", "only-2pcf"):
        options = "no-out-Hist,no-smooth-pivot,"+token
        a = reference.run_native(*catalog, 1, options=options, engine=ENGINE)
        b = reference.run_native(*catalog, 1, options=options+",shear-pivot-reuse",
                                 engine=ENGINE)
        for key in a:
            np.testing.assert_array_equal(a[key], b[key])
    a = reference.run_native(*catalog, 1, options=BASE, engine=ENGINE)
    with patch.dict(os.environ, {"CBALLS_SHEAR_PIVOT_TOL": "0"}):
        b = reference.run_native(*catalog, 1, options=REUSE, engine=ENGINE)
    compare(a, b, 0, 0)


def clustered_catalog():
    rng = np.random.default_rng(881)
    centers = rng.normal(size=(12, 3))
    centers /= np.linalg.norm(centers, axis=1)[:, None]
    position = np.repeat(centers, 32, axis=0) + rng.normal(size=(384, 3))*1.e-5
    position /= np.linalg.norm(position, axis=1)[:, None]
    gamma = .02*(rng.normal(size=384) + 1j*rng.normal(size=384))
    weight = rng.uniform(.1, 2., 384)
    weight[::31] = 0.0
    return position, gamma, weight


def test_aggregate_pivots_and_determinism():
    catalog = clustered_catalog()
    with patch.object(reference, "ENGINE", ENGINE), patch.dict(
            os.environ, {"CBALLS_SHEAR_PIVOT_TOL": "0.1"}):
        exact = reference.run_native(*catalog, 1, options=BASE+",no-one-ball")
        result, profile = reference.run_native_with_profile(*catalog, 1, REUSE)
        count = sum(map(int, re.findall(r"aggregate_pivots=(\d+)", profile)))
        assert count > 0, profile
        # This fixture has nearly coincident groups and deliberately weak modes.
        for key in result:
            relative_l2 = np.linalg.norm(result[key]-exact[key])/np.linalg.norm(exact[key])
            assert relative_l2 < .01, (key, relative_l2)
        threaded = reference.run_native(*catalog, 4, options=REUSE)
        compare(result, threaded, 0, 0)
        combined = reference.run_native(*catalog, 1,
            options="no-out-Hist,no-smooth-pivot,shear-pivot-reuse")
        compare(result, combined, 0, 0)
        pairs = reference.run_native(*catalog, 1,
            options="no-out-Hist,no-smooth-pivot,only-2pcf")
        compare(pairs, combined, 0, 0)


def test_cross_masks_and_polar_geometry():
    first = reference.fixture(42)
    position, gamma, weight = reference.fixture(51)
    position[:4] = [[0, 0, 1], [0, 0, -1], [1, 0, 0], [-1, 1.e-6, 0]]
    position /= np.linalg.norm(position, axis=1)[:, None]
    second = position, gamma, weight
    masks = [np.arange(len(c[0])) % 4 != 0 for c in (first, second)]
    with patch.dict(os.environ, {"CBALLS_SHEAR_PIVOT_TOL": "1e-10"}), \
            patch.object(reference, "RMAX", 2.0):
        for catalogs, order, selected_masks in (
                ([first, second], None, masks),
                ([first, second, first], "1,2,1", masks+[masks[0]])):
            for nmax in (2, 5):
                with patch.object(reference, "NMAX", nmax):
                    kwargs = dict(engine=ENGINE, masks=selected_masks,
                                  catalog_order=order, use_log=False)
                    a = reference.run_native_catalogs(catalogs, 1,
                        options=BASE+",read-mask,no-one-ball", **kwargs)
                    b = reference.run_native_catalogs(catalogs, 4,
                        options=REUSE+",read-mask", **kwargs)
                    compare(a, b, rtol=2.e-9, atol=1.e-10)


def test_inherited_rings_and_masked_aggregation():
    rng = np.random.default_rng(773)
    centers, _, _ = reference.fixture(64)
    position = np.repeat(centers, 32, axis=0) + rng.normal(size=(2048, 3))*.004
    position /= np.linalg.norm(position, axis=1)[:, None]
    gamma = .02*(rng.normal(size=2048) + 1j*rng.normal(size=2048))
    weight = rng.uniform(.5, 1.5, 2048)
    with patch.object(reference, "ENGINE", ENGINE), patch.dict(
            os.environ, {"CBALLS_SHEAR_PIVOT_TOL": "0.5"}):
        actual, profile = reference.run_native_with_profile(
            position, gamma, weight, 4, REUSE)
        assert sum(map(int, re.findall(r"ancestor_merges=(\d+)", profile))) > 0, profile
        exact = reference.run_native(position, gamma, weight, 1,
                                     options=BASE+",no-one-ball")
        for key in actual:
            error = np.linalg.norm(actual[key]-exact[key])/np.linalg.norm(exact[key])
            assert error < .05, (key, error)
    catalog = clustered_catalog()
    masks = [(np.arange(384) % 7 == 0).astype(np.int32)]
    with patch.dict(os.environ, {"CBALLS_SHEAR_PIVOT_TOL": "0.03"}):
        for catalogs, selected_masks in (([catalog], masks),
                                         ([catalog, catalog], masks*2)):
            exact = reference.run_native_catalogs(catalogs, 1,
                options=BASE+",no-one-ball,read-mask", masks=selected_masks, engine=ENGINE)
            actual = reference.run_native_catalogs(catalogs, 4,
                options=REUSE+",read-mask", masks=selected_masks, engine=ENGINE)
            for key in actual:
                error = np.linalg.norm(actual[key]-exact[key])/np.linalg.norm(exact[key])
                assert error < .02, (key, error)


def test_validation_and_smoothing_fallback():
    catalog = reference.smooth_fixture()
    options = "no-out-Hist,smooth-pivot,only-3pcf"
    kwargs = dict(engine=ENGINE, rsmooth=reference.SMOOTH_RMIN_ARCMIN)
    a = reference.run_native(*catalog, 1, options=options, **kwargs)
    b = reference.run_native(*catalog, 1, options=options+",shear-pivot-reuse", **kwargs)
    compare(a, b, 0, 0)
    try:
        reference.run_native(*catalog, 1, options=REUSE+",legacy-one-ball", engine=ENGINE)
    except Exception as error:
        assert "legacy-one-ball" in str(error), error
    else:
        raise AssertionError("reuse unexpectedly accepted in legacy mode")
    for text in ("nan", "inf", "-1", "3.1", "", "0.1junk"):
        with patch.dict(os.environ, {"CBALLS_SHEAR_PIVOT_TOL": text}):
            try:
                reference.run_native(*catalog, 1, options=REUSE, engine=ENGINE)
            except Exception as error:
                assert "CBALLS_SHEAR_PIVOT_TOL" in str(error), error
            else:
                raise AssertionError("invalid tolerance accepted: "+text)


if __name__ == "__main__":
    for test in (test_exact_limit_and_fallbacks, test_aggregate_pivots_and_determinism,
                 test_cross_masks_and_polar_geometry, test_inherited_rings_and_masked_aggregation,
                 test_validation_and_smoothing_fallback):
        test()
        print("PASS", test.__name__, flush=True)
