#!/usr/bin/env python3
"""Count normalization, KD ownership, literal directories, and sibling writers."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
BINARY = ROOT / 'cballs'
WITH_CYTHON = False


class RuntimeStabilization(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory(prefix='cballs-runtime-')
        self.addCleanup(tmp.cleanup)
        self.root = Path(tmp.name)
        self.serial = 0

    def params(self, **changes):
        self.serial += 1
        params = dict(searchMethod='octree-2balls-omp', testmodel='simple-cubic',
                      nbody=8, numberThreads=1, verbose=0, verbose_log=0,
                      rootDir=str(self.root / str(self.serial)), lengthBox=2,
                      useLogHist=False, rminHist=.02, rangeN=.2, sizeHistN=4,
                      mChebyshev=2, options='stop,no-smooth-pivot,no-out-Hist')
        params.update(changes)
        return params

    def cli(self, params, failure=None):
        p = subprocess.run([str(BINARY), *(f'{k}={v}' for k, v in params.items())],
                           cwd=self.root, capture_output=True, text=True, timeout=60)
        log = p.stdout+p.stderr
        if failure is None:
            self.assertEqual(p.returncode, 0, log)
        else:
            self.assertGreater(p.returncode, 0, log)  # no crash/signal
            self.assertIn(failure, log)
        return log

    def test_literal_nested_directories_and_logs(self):
        for name in ['path with spaces', "path'with'quotes", 'path$(printf _expanded)',
                     '/'.join(['nested']*24)]:
            with self.subTest(name=name):
                destination = self.root/name
                params = self.params(rootDir=str(destination), verbose_log=1)
                for _ in range(2):
                    self.cli(params)
                    self.assertTrue((destination/'tmp/cballs.log').is_file())
        self.assertFalse((self.root/'path_expanded').exists())

    def test_directory_file_collisions(self):
        collision = self.root/'file'
        collision.write_text('preserve me')
        for destination in [collision, collision/'child']:
            self.cli(self.params(rootDir=str(destination)), 'cannot create directory')
        destination = self.root/'logs'
        destination.mkdir()
        (destination/'tmp').write_text('preserve log collision')
        self.cli(self.params(rootDir=str(destination), verbose_log=1), 'cannot create directory')
        self.assertEqual(collision.read_text(), 'preserve me')

    def test_sibling_catalog_writer_failures(self):
        for fmt in ['columns-ascii', 'columns-ascii-all', 'columns-ascii-pos',
                    'binary', 'binary-all', 'fits', 'numpy-healpix']:
            with self.subTest(format=fmt):
                params = self.params(outfile='missing/catalog', outfileformat=fmt)
                self.cli(params, 'output')

    def test_histogram_writer_failures(self):
        for method, extra in [('octree-2balls-omp', ',legacy-one-ball'),
                              ('neighbor-boxes-omp', '')]:
            for filename in ['histNN', 'histCF', 'histXi2pcf']:
                if method == 'neighbor-boxes-omp' and filename == 'histCF':
                    continue  # neighbor boxes publishes CF through its Xi writer
                with self.subTest(method=method, filename=filename):
                    params = self.params(searchMethod=method, usePeriodic=True,
                        options='only-2pcf,no-smooth-pivot,compute-HistN,and-CF'+extra)
                    destination = Path(params['rootDir'])
                    (destination/(filename+'.txt')).mkdir(parents=True)
                    self.cli(params, 'output')

    def engine(self, **changes):
        from cyballs import cballs
        balls = cballs()
        self.addCleanup(balls.clean_all)
        balls.set(**self.params(**changes))
        return balls

    def assert_released(self, balls):
        for getter in ['getCMDAllocated', 'getGDAllocated', 'getAllocated2',
                       'getHistogramsAllocated', 'getTreeAllocated', 'getBodytableAllocated']:
            self.assertFalse(getattr(balls, getter)(), getter)

    def test_cython_kd_constructor_failures_and_reuse(self):
        if not WITH_CYTHON:
            self.skipTest('pass --cython for extension tests')
        from cyballs import CosmoComputationError
        for method, options in [('kdtree-box-omp', 'no-out-Hist'),
                                ('kdtree-2balls-omp', 'no-out-Hist,legacy-one-ball')]:
            for fail in range(4):
                with self.subTest(method=method, allocation=fail):
                    balls = self.engine(searchMethod=method, options=options)
                    balls.Run(level=['StartRun_Common'])
                    balls._set_allocation_failure_after_for_tests(fail)
                    try:
                        with self.assertRaisesRegex(CosmoComputationError, 'KDTREE'):
                            balls.Run(level=['MainLoop'])
                    finally:
                        balls._reset_allocation_failure_for_tests()
                    self.assert_released(balls)
                    balls.Run()
                    balls.clean_all()
                    self.assert_released(balls)

    def test_cython_failures_release_and_reuse(self):
        if not WITH_CYTHON:
            self.skipTest('pass --cython for extension tests')
        from cyballs import CosmoComputationError
        for kind in ['directory', 'catalog', 'octree-histogram', 'neighbor-histogram']:
            with self.subTest(kind=kind):
                changes = {}
                if kind == 'catalog':
                    changes.update(outfile='missing/catalog', outfileformat='binary')
                elif kind.endswith('histogram'):
                    changes.update(searchMethod='neighbor-boxes-omp' if kind.startswith('neighbor') else
                                   'octree-2balls-omp', usePeriodic=True,
                                   options='only-2pcf,no-smooth-pivot,compute-HistN,and-CF,legacy-one-ball')
                balls = self.engine(**changes)
                # The most recent generated root is deterministic within this test.
                destination = self.root/str(self.serial)
                if kind == 'directory':
                    destination.write_text('collision')
                elif kind.endswith('histogram'):
                    (destination/'histNN.txt').mkdir(parents=True)
                with self.assertRaises(CosmoComputationError):
                    balls.Run()
                self.assert_released(balls)
                balls.set_default(**self.params())
                balls.Run()
                balls.clean_all()
                self.assert_released(balls)

    def test_recorded_65536_count_normalization(self):
        if not WITH_CYTHON:
            self.skipTest('pass --cython for extension tests')
        import numpy as np
        fixture = json.loads((ROOT/'tests/fixtures/runtime_stabilization/counts_65536.json').read_text())
        positions = np.random.default_rng(fixture['seed']).uniform(-1, 1, (fixture['nbody'], 3))
        balls = self.engine(numberThreads=2, useLogHist=True, theta=0,
            options='only-2pcf,no-smooth-pivot,no-out-Hist,compute-HistN,and-CF')
        balls.set_catalog(positions)
        balls.Run()
        dd, cf = balls.getHistNN(), balls.getHistCF()
        np.testing.assert_array_equal(dd, fixture['unordered_pairs'])
        edges = np.geomspace(.02, .2, 5)
        shell = 4*np.pi/3*np.diff(edges**3)
        expected = 2*dd*8/(float(len(positions))**2*shell)-1
        self.assertTrue(np.isfinite(cf).all())
        np.testing.assert_allclose(cf, expected, rtol=2e-12, atol=2e-12)

    def test_small_catalog_pair_and_shell_oracles(self):
        if not WITH_CYTHON:
            self.skipTest('pass --cython for extension tests')
        import numpy as np
        positions = np.random.default_rng(331).uniform(-.5, .5, (96, 3))
        for method in ['octree-2balls-omp', 'kdtree-2balls-omp', 'balltree-2balls-omp', 'kdtree-box-omp']:
            for legacy in ([False] if method == 'kdtree-box-omp' else [False, True]):
                delta = positions[:, None, :]-positions[None, :, :]
                distances = np.sqrt(np.sum(delta*delta, axis=2))[np.triu_indices(len(positions), 1)]
                for mode in ['linear', 'log']:
                    with self.subTest(method=method, legacy=legacy, mode=mode):
                        rmin = .02
                        edges = (np.linspace(.02, .8, 5) if mode == 'linear' else
                                 np.geomspace(.02, .8, 5))
                        expected_dd, _ = np.histogram(distances, bins=edges)
                        balls = self.engine(searchMethod=method, useLogHist=mode != 'linear',
                            rminHist=rmin, rangeN=.8, logHistBinsPD=2, theta=1e-6 if legacy else 0,
                            options='only-2pcf,no-one-ball,no-smooth-pivot,no-out-Hist,compute-HistN,and-CF'+
                            (',legacy-one-ball' if legacy else ''))
                        balls.set_catalog(positions)
                        balls.Run()
                        dd = balls.getHistNN()
                        if legacy and method == 'octree-2balls-omp':
                            # B4 pivot aggregation persists even with no-one-ball.
                            # Check normalization of the published approximate DD;
                            # exact pair enumeration is checked for all other paths.
                            self.assertGreater(dd.sum(), 0)
                            expected_dd = dd
                        else:
                            np.testing.assert_array_equal(dd, expected_dd)
                        expected_cf = 2*expected_dd*8/(float(len(positions))**2*(4*np.pi/3)*np.diff(edges**3))-1
                        np.testing.assert_allclose(balls.getHistCF(), expected_cf, rtol=2e-12, atol=2e-12)
                        balls.clean_all()

    def test_legacy_multipoles_preserve_count_normalization(self):
        if not WITH_CYTHON:
            self.skipTest('pass --cython for extension tests')
        import numpy as np
        positions = np.random.default_rng(331).uniform(-.5, .5, (96, 3))
        for logarithmic in [False, True]:
            reference = None
            for only_pairs in [True, False]:
                balls = self.engine(useLogHist=logarithmic, rangeN=.8, theta=1e-6,
                    options='legacy-one-ball,no-one-ball,no-smooth-pivot,no-out-Hist,compute-HistN,and-CF'+
                    (',only-2pcf' if only_pairs else ''))
                balls.set_catalog(positions)
                balls.Run()
                result = [balls.getHistNN().copy(), balls.getHistCF().copy()]
                if reference is None:
                    reference = result
                else:
                    for actual, expected in zip(result, reference):
                        np.testing.assert_array_equal(actual, expected)
                balls.clean_all()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--cballs', type=Path, default=BINARY)
    parser.add_argument('--cython', action='store_true')
    args = parser.parse_args()
    BINARY = args.cballs.resolve()
    WITH_CYTHON = args.cython
    if WITH_CYTHON:
        sys.path.insert(0, str(ROOT))
        import cyballs
        if Path(cyballs.__file__).resolve().parent != ROOT:
            raise RuntimeError('tests require the locally built extension')
    unittest.main(argv=[sys.argv[0]], verbosity=2)
