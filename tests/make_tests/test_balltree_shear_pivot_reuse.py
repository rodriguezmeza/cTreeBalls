"""Ball-tree coverage of the shared bounded-phase pivot reuse algorithm."""
from unittest.mock import patch
import os

import numpy as np

import test_shear_pivot_reuse as shared


def test_calibrated_complex_coefficients():
    # Deliberately clustered synthetic fixture, not a survey-wide guarantee.
    rng = np.random.default_rng(9202026)
    centers = rng.normal(size=(512, 3))
    centers /= np.linalg.norm(centers, axis=1)[:, None]
    position = np.repeat(centers, 128, axis=0) + rng.normal(size=(65536, 3))*1.e-6
    position /= np.linalg.norm(position, axis=1)[:, None]
    longitude = np.arctan2(position[:, 1], position[:, 0])
    latitude = np.arcsin(position[:, 2])
    gamma = (.018+.006*np.cos(3*latitude))*np.exp(2j*(longitude+.21*np.sin(latitude)))
    weight = rng.uniform(.65, 1.35, len(position))
    indices = np.concatenate([np.arange(i*128, i*128+32) for i in range(512)])
    catalog = position[indices], gamma[indices], weight[indices]
    reference = shared.reference
    settings = dict(RMIN=2*np.sin(np.deg2rad(5)/2),
                    RMAX=2*np.sin(np.deg2rad(60)/2), BINS=6, NMAX=3)
    kwargs = dict(engine="balltree-shear-sphere-2balls-omp", use_log=True, nsmooth=8)
    with patch.multiple(reference, **settings):
        exact = reference.run_native_catalogs([catalog], 4,
            options=shared.BASE+",no-one-ball,no-two-balls", **kwargs)
        for tolerance, expected_pass in [(0.3, True), (1., False)]:
            with patch.dict(os.environ, CBALLS_SHEAR_PIVOT_TOL=str(tolerance)):
                actual = reference.run_native_catalogs([catalog], 4,
                    options=shared.REUSE, **kwargs)
            passed = True
            for name in exact:
                delta = np.abs(actual[name] - exact[name])
                amp = np.abs(exact[name])
                passed &= bool(np.all(delta <= .05*amp))
            assert passed == expected_pass, (tolerance, passed)


def test_balltree_reuse():
    with patch.object(shared, "ENGINE", "balltree-shear-sphere-2balls-omp"):
        for test in (
            shared.test_exact_limit_and_fallbacks,
            shared.test_aggregate_pivots_and_determinism,
            shared.test_cross_masks_and_polar_geometry,
            shared.test_inherited_rings_and_masked_aggregation,
            shared.test_validation_and_smoothing_fallback,
        ):
            test()
            print("PASS balltree", test.__name__, flush=True)


if __name__ == "__main__":
    test_balltree_reuse()
    test_calibrated_complex_coefficients()
    print("PASS balltree test_calibrated_complex_coefficients", flush=True)
