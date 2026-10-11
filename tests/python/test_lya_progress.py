"""Native Ly-alpha progress routing, task completion, and numerical invariance.

Run with unittest. CBALLS_PROGRESS_BASELINE optionally points to a pre-change
binary; CBALLS_PROGRESS_MPI=1 includes two-rank coverage with mpiexec.
"""
import os
from pathlib import Path
import re
import shlex
import subprocess
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
BINARY = Path(os.environ.get('CBALLS', ROOT / 'cballs')).resolve()
BASELINE = Path(os.environ.get('CBALLS_PROGRESS_BASELINE', BINARY)).resolve()
METHODS = ('lya-2pcf', 'lya-3pcf', 'lya-2pcf-3pcf',
           'lya-los-tree-2pcf', 'lya-los-tree-3pcf', 'lya-los-tree-2pcf-3pcf',
           'lya-1d-2pcf', 'lya-1d-3pcf', 'lya-1d-2pcf-3pcf',
           'lya-1d-tree-2pcf', 'lya-1d-tree-3pcf')
OMP_METHODS = tuple(m+'-omp' for m in METHODS) + (
    'lya-1d-tree-same-los-2pcf-omp', 'lya-anisotropic-multipole-3pcf-omp')
PATTERN = re.compile(r'^(lya-\S+): (.+) progress: rank (\d+) completed tasks '
                     r'(\d+) / (\d+) \(([\d.]+)%\); elapsed ([\d.]+) s$', re.M)


def catalog(path, forests=32, pixels=20):
    rng = np.random.default_rng(91026)
    los = np.column_stack((np.ones(forests), rng.uniform(-.006, .006, (forests, 2))))
    los /= np.linalg.norm(los, axis=1)[:, None]
    radius = 100 + np.arange(forests)[:, None]*1.7 + np.arange(pixels)[None, :]*.2
    xyz = (radius[:, :, None]*los[:, None, :]).reshape(-1, 3)
    data = np.column_stack((xyz, rng.normal(0, .2, len(xyz)),
                            rng.uniform(.5, 2, len(xyz)), np.repeat(np.arange(forests), pixels)))
    data[::17, 4] = 0  # Filtered pixels must not leave progress incomplete.
    np.savetxt(path, data, fmt=['%.17g']*5+['%.0f'])


def command(binary, engine, cat, output, threads=4, verbose=1, verbose_log=0, extra=None):
    params = dict(searchMethod=engine, infile=cat, infileformat='lya-ascii', iCatalogs=1,
        rootDir=output, numberThreads=threads, usePeriodic='false', useLogHist='false',
        rangeN=6, rminHist=.01, sizeHistN=4, lya2RpMax=6, lya2RtMax=6,
        lya2RpBins=4, lya2RtBins=4, lya3RMax=6, lya3RBins=3, lya3ThetaBins=3,
        lya3MuBins=4, lya3LMax=3, verbose=verbose, verbose_log=verbose_log,
        options='no-smooth-pivot,lya-output-empty-bins')
    params.update(extra or {})
    prefix = shlex.split(os.environ.get('MPIEXEC', 'mpiexec')) + ['-n', '2'] if engine.endswith('-mpi') else []
    return prefix + [str(binary)] + [f'{k}={v}' for k, v in params.items()]


def run(binary, engine, cat, output, **kwargs):
    proc = subprocess.run(command(binary, engine, cat, output, **kwargs),
                          capture_output=True, text=True, timeout=120, cwd=cat.parent)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    products = {p.name: np.loadtxt(p, ndmin=2) for p in output.glob('hist*lya*.txt')}
    assert products, (engine, list(output.iterdir()))
    log = '\n'.join(p.read_text() for p in output.rglob('cballs.log'))
    return products, proc.stdout, log


def check_progress(text):
    phases = {}
    for engine, phase, rank, done, total, pct, elapsed in PATTERN.findall(text):
        phases.setdefault((engine, phase, rank), []).append(
            (int(done), int(total), float(pct), float(elapsed)))
    assert phases, text
    for (engine, phase, rank), rows in phases.items():
        assert rank == '0'  # Preserve MPI's rank-zero publication policy.
        totals = {r[1] for r in rows}
        assert len(totals) == 1, rows
        total = totals.pop()
        assert rows[0][0] == 0 and rows[-1][0] == total, rows
        assert rows[-1][2] == 100.0, rows
        done = [r[0] for r in rows]
        assert done == sorted(set(done)), rows
        elapsed = [r[3] for r in rows]
        assert elapsed == sorted(elapsed), rows
        if total > 1:
            assert len(rows) >= 3, rows
        for n, _, pct, _ in rows:
            assert abs(pct - (100*n/total if total else 100)) <= .051
        assert f'{engine}: {phase} counting complete on rank {rank};' in text
    return phases


class LyaProgressTests(unittest.TestCase):
    def compare(self, engine, *, threads=4, extra=None):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cat = root/'catalog.txt'
            catalog(cat)
            quiet, stdout, log = run(BASELINE, engine, cat, root/'quiet',
                threads=threads, verbose=0, extra=extra)
            self.assertFalse(PATTERN.search(stdout+log))
            loud, stdout, log = run(BINARY, engine, cat, root/'loud',
                threads=threads, verbose=1, verbose_log=1, extra=extra)
            terminal = check_progress(stdout)
            self.assertEqual(terminal, check_progress(log))
            self.assertEqual(quiet.keys(), loud.keys())
            for key in quiet:
                np.testing.assert_array_equal(quiet[key], loud[key], err_msg=f'{engine}: {key}')
            return terminal

    def test_registered_omp_engines(self):
        for engine in OMP_METHODS:
            for threads in (1, 4):
                with self.subTest(engine=engine, threads=threads):
                    self.compare(engine, threads=threads)

    def test_accelerated_and_frontier_paths(self):
        for prefix in ('lya', 'lya-los-tree'):
            for kernel in range(6):
                with self.subTest(prefix=prefix, kernel=kernel):
                    self.compare(prefix+'-2pcf-3pcf-omp', extra=dict(
                        lya2Kernel=1, lya3Kernel=kernel, lya3PivotCellMax=8))
        for extra in (dict(lyaScanLevel=3), dict(lyaPivotRadius=.3, lyaPivotMax=4),
                      dict(lya3Kernel=5, lya3PivotCellMax=1),
                      dict(lya3Kernel=4), dict(lya3Kernel=5),
                      dict(lya2Kernel=1, lya2RpSlop=.1, lya2RtSlop=.1)):
            with self.subTest(extra=extra):
                self.compare('lya-2pcf-3pcf-omp', extra=extra)
        for kernel in (0, 1, 2):
            with self.subTest(multipole_kernel=kernel):
                self.compare('lya-anisotropic-multipole-3pcf-omp',
                             extra=dict(lya3Kernel=kernel, lya3MuMode=1))

    def test_log_only_and_throttle(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); cat = root/'catalog.txt'; catalog(cat)
            _, stdout, log = run(BINARY, 'lya-3pcf-omp', cat, root/'log',
                verbose=0, verbose_log=1, extra=dict(lya3PivotBlock=1))
            self.assertFalse(PATTERN.search(stdout))
            phases = check_progress(log)
            for rows in phases.values():
                self.assertGreater(rows[-1][1], 100)
                self.assertLess(len(rows), rows[-1][1])
            _, stdout, log = run(BINARY, 'lya-3pcf-omp', cat, root/'quiet', verbose=0)
            self.assertFalse(PATTERN.search(stdout+log))

    def test_output_arrives_while_running(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); cat = root/'catalog.txt'; catalog(cat)
            with subprocess.Popen(command(BINARY, 'lya-3pcf-omp', cat, root/'out'),
                                  stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                  text=True, bufsize=1) as proc:
                lines = []
                live = False
                for line in proc.stdout:
                    lines.append(line)
                    match = PATTERN.match(line.strip())
                    if match and 0 < int(match[4]) < int(match[5]) and proc.poll() is None:
                        live = True
                self.assertEqual(proc.wait(timeout=120), 0, ''.join(lines))
                self.assertTrue(live, 'no intermediate progress observed while counting')
                check_progress(''.join(lines))

    @unittest.skipUnless(os.environ.get('CBALLS_PROGRESS_MPI') == '1', 'opt-in MPI tests')
    def test_mpi_engines_and_acceleration(self):
        for method in METHODS:
            with self.subTest(engine=method):
                self.compare(method+'-mpi', threads=2)
        for kernel in (3, 4, 5):
            with self.subTest(kernel=kernel):
                self.compare('lya-los-tree-2pcf-3pcf-mpi', threads=2,
                             extra=dict(lya2Kernel=1, lya3Kernel=kernel))


if __name__ == '__main__':
    unittest.main()
