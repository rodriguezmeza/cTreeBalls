"""DES table conventions, deterministic selection and public driver contracts."""
import contextlib
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'tests/python'))
import shear_fits_catalog as fits_catalog
import shear_corr_all_engines as driver
from shear_products import window_diagnostics


def write_catalog(path, *, convention='T17RAW', count=48, invalid=False,
                  names=('x', 'y', 'z', 'gamma1', 'gamma2'), vector=False):
    from astropy.io import fits
    rng = np.random.default_rng(1805)
    xyz = np.column_stack((np.ones(count), rng.uniform(-.04, .04, (count, 2))))
    xyz /= np.linalg.norm(xyz, axis=1)[:, None]
    g1, g2 = rng.normal(size=(2, count)) * .02
    if invalid:
        xyz[0] = 0
        xyz[1, 0] = np.nan
        g1[2] = np.inf
        g2[3] = np.nan
    arrays = [*xyz.T, g1, g2]
    columns = [fits.Column(name=name, format='2E' if vector and i == 4 else 'E',
                          array=np.column_stack((value, value)) if vector and i == 4 else value)
               for i, (name, value) in enumerate(zip(names, arrays))]
    # Optional kappa is intentionally absent; it is not required for shear.
    table = fits.BinTableHDU.from_columns(columns)
    for key, value in dict(REALIZ=2, TOMOBIN=1, REGION=1, NPIXCAT=count,
                           NSIDE=4096, WTSUM=73., G1MEAN=12., G2MEAN=15.).items():
        table.header[key] = value
    if convention is not None:
        table.header['G2CONV'] = convention
    fits.HDUList([fits.PrimaryHDU(), table]).writeto(path, overwrite=True)
    return xyz.astype('f4').astype('f8'), g1.astype('f4').astype('f8'), g2.astype('f4').astype('f8')


class ShearFitsCatalogTests(unittest.TestCase):
    def setUp(self):
        try:
            from astropy.io import fits
        except ImportError:
            self.skipTest('astropy is required for FITS fixtures')
        self.fits = fits
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.path = self.root/'DESY3_Takahashi_r002_bin1_region1.fits'

    def test_auto_detects_sparse_table_and_preserves_raw_science(self):
        xyz, g1, g2 = write_catalog(self.path, names=('X', 'Y', 'Z', 'GAMMA1', 'GAMMA2'))
        self.assertEqual(fits_catalog.resolve_fits_format(self.path), 'desy3')
        result = fits_catalog.load_des_shear_catalog(self.path, chunk_rows=7)
        np.testing.assert_allclose(result['positions'], xyz/np.linalg.norm(xyz, axis=1)[:, None], atol=2e-16)
        np.testing.assert_array_equal(result['gamma1'], g1)
        np.testing.assert_array_equal(result['gamma2'], -g2)
        catalog = driver.ShearCatalog(**result).normalized()
        np.testing.assert_array_equal(catalog.weights, np.ones(len(g1)))
        self.assertFalse(catalog.metadata['centered'])
        self.assertFalse(catalog.metadata['source_plane_weights_renormalized'])
        self.assertEqual(catalog.metadata['fits_header']['WTSUM'], 73.)
        self.assertEqual((catalog.metadata['realization'], catalog.metadata['tomobin'], catalog.metadata['region']), (2, 1, 1))

    def test_local_conventions_do_not_flip_again(self):
        for name in ('LOCALEN', 'EASTNORTH'):
            with self.subTest(convention=name):
                _, _, g2 = write_catalog(self.path, convention=name)
                np.testing.assert_array_equal(fits_catalog.load_des_shear_catalog(self.path)['gamma2'], g2)

    def test_unknown_convention_requires_explicit_choice(self):
        for convention in (None, 'unknown'):
            with self.subTest(convention=convention):
                _, _, g2 = write_catalog(self.path, convention=convention)
                with self.assertRaisesRegex(ValueError, 'G2CONV'):
                    fits_catalog.load_des_shear_catalog(self.path)
                np.testing.assert_array_equal(fits_catalog.load_des_shear_catalog(self.path, convention='takahashi')['gamma2'], -g2)
                np.testing.assert_array_equal(fits_catalog.load_des_shear_catalog(self.path, convention='local-east-north')['gamma2'], g2)

    def test_sampling_matches_global_bottom_k_after_filtering(self):
        xyz, g1, g2 = write_catalog(self.path, count=73, invalid=True)
        eligible = np.arange(4, 73, dtype=np.int64)
        take = np.sort(eligible[np.argsort(fits_catalog.splitmix64(eligible, 42))[:13]])
        for chunk in (1, 7, 100):
            with self.subTest(chunk=chunk):
                result = fits_catalog.load_des_shear_catalog(self.path, max_points=13, seed=42, chunk_rows=chunk)
                np.testing.assert_array_equal(result['gamma1'], g1[take])
                np.testing.assert_array_equal(result['gamma2'], -g2[take])
                self.assertEqual(result['metadata']['invalid_rows_removed'], 4)
                self.assertEqual(result['metadata']['eligible_pixels_before_thinning'], 69)

    def test_invalid_control_values_and_too_few_rows(self):
        write_catalog(self.path, count=6, invalid=True)
        with self.assertRaisesRegex(ValueError, 'fewer than three'):
            fits_catalog.load_des_shear_catalog(self.path)
        for kwargs in ({'max_points': 2}, {'max_points': -1}, {'seed': -1}, {'seed': 2**64}, {'chunk_rows': 0}):
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                fits_catalog.load_des_shear_catalog(self.path, **kwargs)

    def test_table_shape_and_identity_validation(self):
        for key, value, message in (('NPIXCAT', 99, 'row count'), ('TOMOBIN', 2, 'mismatch'), ('REGION', 7, 'invalid')):
            with self.subTest(key=key):
                write_catalog(self.path)
                with self.fits.open(self.path, mode='update') as hdus:
                    hdus[1].header[key] = value
                with self.assertRaisesRegex(ValueError, message):
                    fits_catalog.inspect_des_catalog(self.path)
        write_catalog(self.path, vector=True)
        with self.assertRaisesRegex(ValueError, 'scalar'):
            fits_catalog.inspect_des_catalog(self.path)

    def test_missing_fields_and_ambiguous_hdus_fail(self):
        write_catalog(self.path, names=('x', 'y', 'z', 'gamma1', 'not_gamma2'))
        self.assertEqual(fits_catalog.resolve_fits_format(self.path), 'desy3')
        with self.assertRaisesRegex(ValueError, 'gamma1/gamma2'):
            fits_catalog.load_des_shear_catalog(self.path)
        write_catalog(self.path)
        with self.fits.open(self.path, mode='update') as hdus:
            hdus.append(hdus[1].copy())
        with self.assertRaisesRegex(ValueError, 'exactly one'):
            fits_catalog.inspect_des_catalog(self.path)

    def test_dense_healpix_detection(self):
        table = self.fits.BinTableHDU.from_columns([
            self.fits.Column(name='GAMMA1', format='E', array=np.ones(12)),
            self.fits.Column(name='GAMMA2', format='E', array=np.ones(12)),
        ])
        self.fits.HDUList([self.fits.PrimaryHDU(), table]).writeto(self.path)
        self.assertEqual(fits_catalog.resolve_fits_format(self.path), 'healpix')

    def test_incompatible_des_controls_are_rejected(self):
        for kwargs in ({'geometry': 'flat'}, {'input_shear_frame': 'flat'}, {'mask': Path('mask')},
                       {'weight_field': 'w'}, {'conjugate_input': True}, {'gamma1_field': 'G1'},
                       {'gamma2_field': 2}, {'nsides': (128,)}):
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                fits_catalog.validate_des_controls(**kwargs)

    def test_binning_presets_match_centers_and_edges(self):
        for name in ('sofia-fig1', 'paper-8-200-edges'):
            lo, hi, count, units = fits_catalog.apply_binning_preset(name, 1, 2, 4, 'degree')
            config = driver.RunConfig(min_sep=lo, max_sep=hi, bins=count, sep_units=units)
            edges = np.geomspace(*config.native_limits('sphere'), count+1)
            chosen = np.sqrt(edges[:-1]*edges[1:])[[0, -1]] if name == 'sofia-fig1' else edges[[0, -1]]
            np.testing.assert_allclose(np.rad2deg(2*np.arcsin(chosen/2))*60, [8, 200], rtol=2e-14)
            with self.assertRaisesRegex(ValueError, 'spherical logarithmic'):
                fits_catalog.apply_binning_preset(name, 1, 2, 4, 'degree', linear=True)

    def test_driver_cli_routes_des_through_shared_loader(self):
        write_catalog(self.path)
        saved = self.root/'sample.npz'
        with patch.object(driver, 'inspect_cballs_runtime', return_value={}), \
             patch.object(driver, 'available_engines', return_value=driver.SHEAR_SPHERE_OMP_ENGINES), \
             patch.object(driver, 'run_engine_suite') as run, contextlib.redirect_stdout(io.StringIO()):
            status = driver.main(['--fits', str(self.path), '--fits-format', 'desy3', '--engine', 'all-omp',
                                  '--max-points', '13', '--sampling-seed', '42', '--fits-chunk-pixels', '7',
                                  '--binning', 'sofia-fig1', '--save-catalog-npz', str(saved)])
        self.assertEqual(status, 0)
        self.assertEqual(run.call_args.args[0].nbody, 13)
        self.assertEqual(run.call_args.args[2].bins, 20)
        self.assertTrue(saved.exists())

    def test_discovery_order_missing_and_duplicate_identities(self):
        write_catalog(self.path)
        self.assertEqual(len(fits_catalog.discover_des_catalogs(self.root, [2], [1], [1])), 1)
        with self.assertRaisesRegex(ValueError, 'missing requested'):
            fits_catalog.discover_des_catalogs(self.root, [1, 2], [1], [1])
        nested = self.root/'nested';nested.mkdir()
        write_catalog(nested/self.path.name)
        with self.assertRaisesRegex(ValueError, 'duplicate'):
            fits_catalog.discover_des_catalogs(self.root, [2], [1], [1])
        self.assertEqual(fits_catalog.integer_selection('2,1:3'), [1, 2, 3])
        with self.assertRaises(ValueError):
            fits_catalog.integer_selection('1:999999999')

    def test_relative_difference_keeps_valid_bins_beside_masked_bins(self):
        actual = driver.symmetric_relative_difference(np.array([np.nan, 1., 2., 0.]),
                                                       np.array([np.nan, 1., 3., 0.]))
        np.testing.assert_allclose(actual, [np.nan, 0., .4, np.nan], equal_nan=True)

    def test_window_validity_distinguishes_singular_and_bad_solves(self):
        n = 2
        window = np.zeros((2, 2, 4*n+1), complex)
        window[0, 0, 2*n] = 2
        window[0, 1, :] = 1  # Rank-one coupling.
        window[1, 0, 2*n] = 2
        upsilon = np.ones((4, 2, 2, 2*n+1), complex)
        gamma = upsilon/2
        gamma[:, 1, 0, :] = 0  # Deliberately incorrect solve.
        result = window_diagnostics(upsilon, window, gamma, n)
        np.testing.assert_array_equal(result['gamma_valid'], [[True, False], [False, False]])
        self.assertEqual(result['gamma_solve_relative_residual'][0, 0], 0)


if __name__ == '__main__':
    unittest.main()
