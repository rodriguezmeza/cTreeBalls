"""Live frontier completion, output routing, and numerical invariance."""
import os
from pathlib import Path
import re
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

ROOT = Path(os.environ.get('CBALLS_PROGRESS_ROOT', Path(__file__).resolve().parents[2]))
sys.path.insert(0, str(ROOT))
from cyballs import cballs

ENGINES = ('balltree-shear-sphere-2balls-omp', 'kdtree-shear-sphere-2balls-omp',
           'octree-shear-sphere-2balls-omp')
CASES = ('body', 'exact', 'reuse', 'combined', 'masked', 'pair')
PATTERN = re.compile(r'3PCF progress: rank (\d+) completed tasks (\d+) / (\d+) '
                     r'\(([\d.]+)%\); elapsed ([\d.]+) s')


def _calculate(engine, case, threads, verbose=0, verbose_log=0):
    rng = np.random.default_rng(10926)
    centers = rng.normal(size=(32, 3))
    centers /= np.linalg.norm(centers, axis=1)[:, None]
    xyz = np.repeat(centers, 8, axis=0) + rng.normal(size=(256, 3))*1.e-5
    xyz /= np.linalg.norm(xyz, axis=1)[:, None]
    options = ['no-out-Hist', 'no-smooth-pivot']
    if case not in ('combined', 'pair'):
        options.append('only-3pcf')
    if case == 'pair':
        options.append('only-2pcf')
    if case == 'exact':
        options.append('no-one-ball')
    if case in ('reuse', 'combined'):
        options.append('shear-pivot-reuse')
    if case == 'masked':
        options.append('read-mask')
    with tempfile.TemporaryDirectory() as output:
        model = cballs()
        try:
            model.set(dict(searchMethod=engine, iCatalogs='1', usePeriodic='false',
                useLogHist='true', rminHist=.03, rangeN=1.75, sizeHistN=5,
                sizeHistPhi=16, mChebyshev=3, theta=1., lengthBox=2.2,
                nsmooth=1, numberThreads=threads, verbose=verbose,
                verbose_log=verbose_log, rootDir=output, options=','.join(options)))
            mask = {'mask': (np.arange(256) % 5 != 0).astype(np.int32)} if case == 'masked' else {}
            model.set_catalog(xyz, gamma1=rng.normal(size=256)*.02,
                gamma2=rng.normal(size=256)*.02, weights=rng.uniform(.5, 1.5, 256), **mask)
            model.Run(level=['MainLoop'])
            result = {}
            if case != 'pair':
                result.update(gamma=model.getShearGammaXMultipoles().copy(),
                    numerator=model.getShearUpsilonXMultipoles().copy(),
                    window=model.getShearWindowMultipoles().copy())
            if case in ('pair', 'combined'):
                result.update(xi=model.getShearXiPlus().copy(),
                    weight=model.getShearXiWeight().copy())
        finally:
            model.struct_cleanup()
        log = '\n'.join(p.read_text() for p in Path(output).rglob('cballs.log'))
    return result, log


def calculate(engine, case, threads, verbose=0, verbose_log=0):
    sys.stdout.flush()
    saved = os.dup(1)
    with tempfile.TemporaryFile() as output:
        try:
            os.dup2(output.fileno(), 1)
            result, log = _calculate(engine, case, threads, verbose, verbose_log)
        finally:
            sys.stdout.flush()
            os.dup2(saved, 1)
            os.close(saved)
        output.seek(0)
        terminal = output.read().decode('utf-8', errors='replace')
    return result, log, terminal


def check_progress(text):
    rows = PATTERN.findall(text)
    assert len(rows) >= 3, text
    ranks, done, totals, percents, seconds = zip(*rows)
    done, totals = list(map(int, done)), list(map(int, totals))
    assert set(ranks) == {'0'}
    assert len(set(totals)) == 1
    assert done == list(range(totals[0]+1)), done
    assert list(map(float, seconds)) == sorted(map(float, seconds))
    for n, pct in zip(done, percents):
        assert abs(float(pct)-100*n/totals[0]) <= .051
    assert percents[-1] == '100.0'
    assert text.index('tasks complete on rank 0; reducing histograms') > text.rindex('3PCF progress:')
    return rows


class ProgressTests(unittest.TestCase):
    def test_progress_and_results(self):
        with patch.dict(os.environ, {'CBALLS_SHEAR_PIVOT_TOL': '3',
                                    'CBALLS_SHEAR_BIN_THETA': '.1'}):
            for engine in ENGINES:
                for case in CASES:
                    for threads in (1, 4):
                        with self.subTest(engine=engine, case=case, threads=threads):
                            quiet, _, terminal = calculate(engine, case, threads)
                            self.assertFalse(PATTERN.search(terminal))
                            loud, log, terminal = calculate(engine, case, threads, 1, 1)
                            if case == 'pair':
                                self.assertFalse(PATTERN.search(terminal+log))
                            else:
                                self.assertEqual(check_progress(terminal), check_progress(log))
                            for name in quiet:
                                np.testing.assert_array_equal(quiet[name], loud[name])
                            reference = os.environ.get('CBALLS_PROGRESS_REFERENCE')
                            if reference:
                                with np.load(reference) as baseline:
                                    for name in loud:
                                        np.testing.assert_array_equal(loud[name], baseline[f'{engine}/{case}/{threads}/{name}'])

    def test_log_only(self):
        _, log, terminal = calculate(ENGINES[0], 'reuse', 4, 0, 1)
        self.assertFalse(PATTERN.search(terminal))
        check_progress(log)


if __name__ == '__main__':
    if len(sys.argv) != 3 or sys.argv[1] != '--snapshot':
        unittest.main()
    os.environ['CBALLS_SHEAR_PIVOT_TOL'] = '3'
    os.environ['CBALLS_SHEAR_BIN_THETA'] = '.1'
    values = {}
    for engine in ENGINES:
        for case in CASES:
            for threads in (1, 4):
                result, _, _ = calculate(engine, case, threads)
                values.update({f'{engine}/{case}/{threads}/{key}': value for key,value in result.items()})
    np.savez(sys.argv[2], **values)
    print('Saved original-extension reference:', len(values), 'arrays')
