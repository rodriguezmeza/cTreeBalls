#!/usr/bin/env python3
"""Fail-closed release gate for the resolved, unchanged active build profile.

No feature switches are enabled by this runner. Unknown registered methods,
missing required MPI execution, failed commands, and missing products fail it.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import signal
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
TESTS = ROOT / 'tests/make_tests'
from capabilities_generated import expected_registry, gate_plan


def parse_registry(text):
    rows = re.findall(r'^- (\S+) \(id=(\d+)\)', text, re.M)
    count = re.search(r'registered in this executable \((\d+)\)', text)
    if not count or len(rows) != int(count[1]) or len(dict(rows)) != len(rows):
        raise AssertionError('missing, duplicate or truncated registry')
    return {name: int(number) for name, number in rows}


def json_line(text):
    return next(json.loads(line) for line in text.splitlines() if line.startswith('{'))


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()


class Gate:
    def __init__(self, args):
        self.args = args
        self.out = args.output.resolve()
        self.out.mkdir(parents=True, exist_ok=False)
        (self.out/'logs').mkdir()
        self.scratch = self.out/'scratch'
        (self.scratch/'Output').mkdir(parents=True)
        self.env = dict(os.environ, PYTHONPATH=os.pathsep.join([str(ROOT)]+[str(Path(p).resolve()) for p in os.environ.get('PYTHONPATH', '').split(os.pathsep) if p]),
                        CBALLS=str(ROOT/'cballs'), PYTHON=sys.executable, OMP_DYNAMIC='FALSE')
        # Recursive make jobserver/overrides belong to the parent, not this fresh build.
        self.env.pop('MAKEFLAGS', None); self.env.pop('MFLAGS', None)
        self.env['MPLCONFIGDIR'] = str(self.scratch/'matplotlib')
        self.report = dict(schema_version=1, status='RUNNING', source_root=str(ROOT),
                           started_utc=datetime.now(timezone.utc).isoformat(),
                           commands=[], cases=[], comparisons=[], failures=[])
        self.save()

    def save(self):
        (self.out/'gate.json').write_text(json.dumps(self.report, indent=2, sort_keys=True)+'\n')

    def command(self, label, command, *, cwd=None, env=None, codes=(0,), timeout=1200):
        log = self.out/'logs'/f'{len(self.report["commands"]):03d}-{label}.log'
        record = dict(label=label, argv=[str(x) for x in command], log=str(log.relative_to(self.out)),
                      cwd=str(cwd or self.scratch), status='RUNNING')
        self.report['commands'].append(record); self.save()
        print(f'[{len(self.report["commands"])}] {label}', flush=True)
        start = time.monotonic()
        with log.open('w') as stream:
            process = subprocess.Popen(record['argv'], cwd=cwd or self.scratch,
                                       env=self.env | (env or {}), stdout=stream, stderr=subprocess.STDOUT,
                                       start_new_session=True)
            try:
                code = process.wait(timeout=timeout)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL); process.wait()
                code = None
        record.update(returncode=code, seconds=round(time.monotonic()-start, 3),
                      status='PASS' if code in codes else 'FAIL')
        self.save()
        if code not in codes:
            raise RuntimeError(f'{label}: exit={code}; see {log}')
        return log.read_text(errors='replace')

    def check(self, label, function):
        try:
            function()
        except Exception as exc:
            self.report['failures'].append(dict(check=label, error=str(exc)))
            print(f'FAIL {label}: {exc}', flush=True)
            self.save()

    def build(self):
        self.command('build-native', ['make', '-B', f'-j{self.args.jobs}', 'cballs', 'cyballs-static-lib',
                                     f'PYTHON={sys.executable}'], cwd=ROOT)
        self.command('build-cython', [sys.executable, 'setup.py', 'build_ext', '--inplace', '--force'],
                     cwd=ROOT, env={'CBALLS_STATIC_LIBRARY_READY': '1'})
        resolved = json_line(self.command('resolved-fingerprint', ['make', '--no-print-directory',
                             'print-build-fingerprint', f'PYTHON={sys.executable}'], cwd=ROOT))
        native = json_line(self.command('native-fingerprint', [ROOT/'cballs', 'options=build-fingerprint',
                           f'rootDir={self.scratch/"fingerprint"}', 'verbose=0', 'verbose_log=0']))
        wrapper = json_line(self.command('cython-fingerprint', [sys.executable, '-c',
                    'import cyballs,json; print(json.dumps(dict(build=cyballs.build_info(),file=cyballs.__file__))); cyballs.cballs()']))
        assert resolved == native == wrapper['build'], 'resolved/native/Cython build mismatch'
        assert Path(wrapper['file']).resolve().parent == ROOT, 'imported extension from another checkout'
        self.report['build'] = resolved
        self.report['artifacts'] = {p.name: sha(p) for p in (ROOT/'cballs', ROOT/'libcballs.a', Path(wrapper['file']))}
        registry = parse_registry(self.command('native-registry', [ROOT/'cballs', 'options=print-search-methods',
                                  f'rootDir={self.scratch/"registry"}', 'verbose=0', 'verbose_log=0'], codes=(1,)))
        assert registry == expected_registry(resolved['resolved_settings']), 'settings/registry mismatch or unhandled active engine'
        assert resolved['resolved_settings']['DEFDIMENSION'] == '3', 'this oracle matrix requires 3D'
        assert resolved['resolved_settings']['BUILD_PRECISION'] == 'double', 'this oracle matrix requires double precision'
        self.report['capability_gate_plan'] = gate_plan(registry)
        self.registry = registry
        self.report['registry'] = registry
        code = 'import cyballs,json; r=json.loads('+repr(json.dumps(registry))+'); assert all(cyballs.search_method_id(n)==i for n,i in r.items())'
        self.command('cython-registry', [sys.executable, '-c', code])
        candidates = [p for p in ROOT.glob('build/**/build-fingerprint.json') if json.loads(p.read_text())['id']==resolved['id']]
        assert len(candidates) == 1, 'ambiguous resolved build manifest'
        self.build_directory = candidates[0].parent
        (self.out/'build-fingerprint.json').write_bytes(candidates[0].read_bytes())
        self.save()

    def mpi_environment(self):
        """Retain the actual library used by two Python ranks, not just a launcher name."""
        self.command('python-environment', [sys.executable, '-m', 'pip', 'freeze'])
        if not any(e.endswith('-mpi') for e in self.registry):
            return
        code = ("from mpi4py import MPI; import json,mpi4py; "
                "c=MPI.COMM_WORLD; assert c.size==2; "
                "rows=c.gather(dict(rank=c.rank, library=MPI.Get_library_version(), "
                "vendor=MPI.get_vendor(), mpi4py=mpi4py.__version__), root=0); "
                "print(json.dumps(rows)) if c.rank==0 else None")
        output = self.command('mpi-environment', shlex.split(self.args.mpi_command)+
                              ['-n', '2', sys.executable, '-c', code])
        rows = next(json.loads(line) for line in output.splitlines() if line.startswith('[{'))
        assert {r['rank'] for r in rows} == {0, 1}
        assert rows[0]['vendor'] == rows[1]['vendor']
        wrapper = self.report['build']['toolchain'].get('mpi_version', '')
        version = re.search(r'Open MPI (\d+)\.(\d+)\.(\d+)', wrapper)
        if version:
            assert rows[0]['vendor'] == ['Open MPI', list(map(int, version.groups()))], \
                'mpi4py library does not match the Open MPI compiler wrapper'
        self.report['mpi_environment'] = rows
        self.save()

    def case(self, engine, ranks, threads):
        tag = f'{engine}-r{ranks}-t{threads}'
        output = self.out/'cases'/tag
        command = [sys.executable, ROOT/'scripts/release_gate_case.py', '--engine', engine,
                   '--output', output, '--threads', str(threads)]
        if engine.endswith('-mpi'):
            command = shlex.split(self.args.mpi_command)+['-n', str(ranks)]+command
        self.command(tag, command)
        metadata = json.loads((output/'result.json').read_text())
        assert metadata['build']['id'] == self.report['build']['id']
        assert metadata['engine'] == engine and metadata['engine_id'] == self.registry[engine]
        assert metadata['parallel']['estimator_ranks'] == ranks
        assert metadata['parallel']['openmp_probe_threads'] == threads
        for key in ('estimator', 'bin_edges', 'coordinate_convention', 'weights', 'masks',
                    'effective_smoothing', 'opening_tolerance', 'precision'):
            assert metadata[key], f'missing result provenance: {key}'
        record = dict(engine=engine, ranks=ranks, threads=threads, directory=str(output.relative_to(self.out)),
                      products_sha256=sha(output/'result.npz'), provenance_sha256=sha(output/'result.json'))
        self.report['cases'].append(record); self.save()
        return output

    def compare(self, left, right):
        import numpy as np
        a, b = np.load(left/'result.npz'), np.load(right/'result.npz')
        assert set(a.files) == set(b.files), 'product sets differ'
        maximum = 0.
        for key in a.files:
            np.testing.assert_array_equal(np.isfinite(a[key]), np.isfinite(b[key]), err_msg=key)
            # Text products are rounded by legacy writers. Raw full-precision arrays
            # use their existing oracle tolerances; forest/3D text compare at 3e-12.
            np.testing.assert_allclose(a[key], b[key], rtol=3e-12, atol=3e-12, equal_nan=True, err_msg=key)
            finite = np.isfinite(a[key]) & np.isfinite(b[key])
            if finite.any(): maximum = max(maximum, float(np.max(np.abs(a[key][finite]-b[key][finite]))))
        self.report['comparisons'].append(dict(left=str(left.relative_to(self.out)), right=str(right.relative_to(self.out)),
                                              rtol=3e-12, atol=3e-12, max_absolute_difference=maximum))
        self.save()

    def engine_cases(self):
        declarations=gate_plan(self.registry)
        for declaration in declarations:
            engine=declaration['engine']
            if declaration['parallel']=='omp':
                self.check(engine, lambda e=engine: self.compare(self.case(e, 1, 1), self.case(e, 1, 2)))
        for declaration in declarations:
            engine=declaration['engine']
            if declaration['parallel']=='mpi':
                def mpi_case(e=engine):
                    omp = self.out/'cases'/f'{e[:-3]}omp-r1-t1'
                    one = self.case(e, 1, 1)
                    many = self.case(e, 2, 2)
                    self.compare(omp, one); self.compare(one, many)
                self.check(engine, mpi_case)
        for engine in self.registry:
            if engine.startswith('lya-los-tree-'):
                original = engine.replace('lya-los-tree-', 'lya-')
                self.check(engine+'-baseline', lambda e=engine, o=original: self.compare(
                    self.out/'cases'/f'{e}-r1-t1', self.out/'cases'/f'{o}-r1-t1'))
        covered = {c['engine'] for c in self.report['cases']}
        assert covered == set(self.registry), f'missing successful engine cases: {set(self.registry)-covered}'

    def scripts(self):
        py = sys.executable
        self.check('affected-regressions', lambda: self.command('affected-regressions',
            [py, ROOT/'scripts/affected_regressions.py', '--changed', 'source/common_histogram.c',
             'source/smooth_pivots.c', 'source/mpi_runtime.c', '--output', self.out/'affected-regressions.json']))
        if any(e.endswith('-mpi') for e in self.registry):
            self.check('mpi-runtime-build', lambda: self.command('mpi-runtime-build',
                ['make','test-mpi-runtime-build',f'PYTHON={py}'],cwd=ROOT))
            executable=self.build_directory/'tests/test_mpi_runtime'
            for ownership in ('owned','borrowed'):
                self.check('mpi-runtime-'+ownership, lambda mode=ownership: self.command('mpi-runtime-'+mode,
                    shlex.split(self.args.mpi_command)+['-n','2',executable,mode]))
        self.check('phase-memory-accuracy-benchmark', lambda: self.command('phase-memory-accuracy-benchmark',
            [py,ROOT/'scripts/benchmark_contracts.py','--output',self.out/'benchmarks','--repeats','2']))
        binary = str(ROOT/'cballs')
        active = self.registry
        mpi = shlex.join(shlex.split(self.args.mpi_command)+['-n','2'])
        def script(name, *args, env=None):
            self.check(name + (env or {}).get('CBALLS_SHEAR_SPHERE_ENGINE',''),
                       lambda: self.command(name.replace('.py',''), [py, TESTS/name, *args], env=env))
        self.check('native-unit-contracts', lambda: self.command('native-unit-contracts',
            ['make', 'test-parameter-parser', 'test-option-cache', 'test-healpix-ordering',
             'test-runtime-stabilization-native', 'test-provenance-window-native', 'test-resource-contracts', 'test-runtime-context', f'PYTHON={py}'], cwd=ROOT))
        self.check('python-unit-contracts', lambda: self.command('python-unit-contracts', [py,'-m','pytest','-q',
            *[TESTS/name for name in ('test_p0_cython_instances.py','test_p2_cython.py',
              'test_provenance_window.py','test_two_ball_edge_cython.py','test_kappa_corr_all_engines.py',
              'test_shear_corr_all_engines.py','test_lya_corr_all_engines.py','test_release_gate.py',
              'test_cython_in_memory_catalog.py','test_p3_cython_startup.py',
              'test_scalar_numerical_contract.py',
              'test_release_packaging.py','test_resource_scientific_contracts.py','test_capability_contracts.py')], '-ra', '-k','not mpi']))
        if 'lya-los-tree-2pcf-omp' in active:
            self.check('los-tree-contracts', lambda: self.command('los-tree-contracts',
                [py, '-m', 'pytest', '-q', '-ra', TESTS/'test_lya_forest_los_tree.py']))
        scalar = [e for e in active if e in ('kdtree-2balls-omp','balltree-2balls-omp','octree-2balls-omp')]
        for engine in scalar:
            name = 'run_test_'+engine.replace('-','_')
            self.check(name, lambda n=name: self.command(n, ['bash',TESTS/n]))
        if scalar:
            script('test_two_ball_edge_corrections.py','--cballs',binary,
                   *[v for e in scalar for v in ('--engine',e)])
        mpi_scalar = [e for e in active if e.endswith('2balls-mpi')]
        if mpi_scalar:
            script('test_two_ball_edge_corrections.py','--cballs',binary,'--mpi-command',mpi,
                   *[v for e in mpi_scalar for v in ('--engine',e)])
        if 'octree-2balls-omp' in active:
            script('test_octree_2balls_mask.py','--cballs',binary,'--cython','--fits',
                   *(['--mpi'] if 'octree-2balls-mpi' in active else []),
                   env={'MPIEXEC':shlex.split(self.args.mpi_command)[0],
                        'MPIEXEC_ARGS':shlex.join(shlex.split(self.args.mpi_command)[1:])})
        for e in active:
            if 'shear' in e:
                script('test_shear_sphere_octree_omp.py', env={'CBALLS_SHEAR_SPHERE_ENGINE':e})
        if 'kdtree-box-omp' in active:
            self.check('box-frontier', lambda: self.command('box-frontier',['bash',TESTS/'run_test_kdtree_box_frontier']))
        if 'neighbor-boxes-omp' in active: script('test_neighbor_boxes_periodic.py')
        if 'lya-2pcf-omp' in active:
            script('test_lya_forest_omp.py'); script('test_lya_forest_1d_omp.py')
        if 'lya-2pcf-mpi' in active: script('test_lya_forest_mpi.py','--cballs',binary,'--mpi-command',mpi,'--cython')
        if 'octree-3pcf-3d-omp' in active:
            script('test_octree_3pcf_3d_omp.py','--cballs',binary)
            script('test_octree_3pcf_3d_survey.py','--cballs',binary)
            script('test_cython_octree_3pcf_3d_omp.py')
        if 'octree-3pcf-3d-mpi' in active:
            script('test_octree_3pcf_3d_mpi.py','--cballs',binary,'--mpi-command',mpi,'--cython','--fits')
        script('test_io_stabilization.py','--cballs',binary,'--cython')
        script('test_runtime_stabilization.py','--cballs',binary,'--cython')
        script('test_release_io_routes.py','--cballs',binary,'--output',str(self.out/'io'))
        self.check('accuracy-acceptance', lambda: self.command('accuracy-acceptance',
            [py, ROOT/'scripts/accuracy_acceptance.py', '--output', self.out/'accuracy']))
        self.drivers()

    def drivers(self):
        py = sys.executable
        directory = ROOT/'tests/python'
        common = ['--threads','2','--no-plots','--statistics','both']
        scalar = [e for e in self.registry if e.endswith('-omp') and e.count('2balls') and 'shear' not in e]
        if scalar:
            self.check('kappa-driver', lambda: self.command('kappa-driver', [py,directory/'kappa_corr_all_engines.py',
                '--synthetic-nside','2','--engines',','.join(scalar),'--outdir',self.out/'drivers/kappa',
                '--theta-min','2','--theta-max','100','--nbins','4','--multipoles','2',
                '--linear-bins','--no-smooth-pivot', *common]))
        shear = [e for e in self.registry if 'shear' in e]
        if shear:
            self.check('shear-driver', lambda: self.command('shear-driver', [py,directory/'shear_corr_all_engines.py',
                '--synthetic-nbody','45','--engines',','.join(shear),'--outdir',self.out/'drivers/shear',
                '--min-sep','2','--max-sep','100','--nbins','4','--multipoles','2',
                '--linear-bins','--no-smooth-pivot',*common]))
        forest = [e for e in self.registry if e.startswith('lya-') and e.endswith('-omp')]
        if forest:
            self.check('lya-driver', lambda: self.command('lya-driver', [py,directory/'lya_corr_all_engines.py',
                '--synthetic','--synthetic-forests','4','--synthetic-pixels','4','--engine',*forest,
                '--output',self.out/'drivers/lya',*common]))

    def finish(self):
        for name, expected in self.report.get('artifacts', {}).items():
            if name == 'libcballs.a':
                # Native unit targets recreate the archive; member timestamps can change.
                self.report['tested_archive_sha256'] = sha(ROOT/name)
            elif sha(ROOT/name) != expected:
                self.report['failures'].append(dict(check='final-artifact-identity', error=f'{name} changed during validation'))
        self.report['status'] = 'FAIL' if self.report['failures'] else 'PASS'
        self.report['finished_utc'] = datetime.now(timezone.utc).isoformat()
        self.save()
        print(f'{self.report["status"]}: {self.out/"gate.json"}', flush=True)
        return int(bool(self.report['failures']))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True, help='new directory; refuses stale results')
    parser.add_argument('--jobs', type=int, default=4)
    parser.add_argument('--mpi-command', default='mpiexec', help='launcher and host options, without -n')
    args = parser.parse_args()
    gate = Gate(args)
    gate.check('build-and-registry', gate.build)
    if not gate.report['failures']:
        gate.check('mpi-environment', gate.mpi_environment)
        gate.check('engine-coverage', gate.engine_cases)
        gate.check('script-matrix', gate.scripts)
        gate.check('final-source-identity', lambda: (
            None if json_line(gate.command('final-fingerprint',['make','--no-print-directory',
                  'print-build-fingerprint',f'PYTHON={sys.executable}'],cwd=ROOT)) == gate.report['build']
            else (_ for _ in ()).throw(AssertionError('source/profile changed during gate'))))
    return gate.finish()


if __name__ == '__main__':
    raise SystemExit(main())
