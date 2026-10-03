"""Profile equality and execution coverage must fail closed."""
from pathlib import Path
import sys
import pytest
sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'scripts'))
from check_installed_package import check_profile
from capabilities_generated import ENGINES, expected_registry
from active_release_gate import regression_coverage, ROOT


def test_installed_profile_is_exact_not_a_count():
    intended={'resolved_settings':{'CCFLAG':'-DADDONS -DLYAFORESTOMP'}}
    expected=expected_registry(intended['resolved_settings'])
    assert check_profile(intended,intended,lambda n:expected.get(n,-1))==expected
    with pytest.raises(AssertionError,match='profile mismatch'):
        check_profile({'resolved_settings':{'CCFLAG':''}},intended,lambda n:expected.get(n,-1))
    extra=next(n for n in ENGINES if n not in expected)
    with pytest.raises(AssertionError,match='registry mismatch'):
        check_profile(intended,intended,lambda n:55 if n==extra else expected.get(n,-1))


def test_declared_module_requires_successful_execution():
    registry={'lya-2pcf-omp':169}
    plan=regression_coverage(registry,[])
    assert 'tests/python/test_lya_pair_cells.py' in plan['missing']
    rows=[dict(status='PASS',argv=['python','-m','pytest',str(ROOT/p)]) for p in plan['required']]
    assert not regression_coverage(registry,rows)['missing']
    rows[-1]['status']='FAIL'
    assert regression_coverage(registry,rows)['missing']


def optimization_checks(directory):
    """Run real validators without pytest's assertion rewriting in this process."""
    import argparse
    import json
    import unittest
    import numpy as np
    from active_release_gate import Gate
    from check_sdist import inspect_archive
    from test_release_packaging import archive

    check = unittest.TestCase()
    intended = {'resolved_settings': {'CCFLAG': '-DADDONS -DLYAFORESTOMP'}}
    expected = expected_registry(intended['resolved_settings'])
    check.assertEqual(check_profile(intended, intended, lambda n: expected.get(n, -1)), expected)
    with check.assertRaisesRegex(AssertionError, 'profile mismatch'):
        check_profile({'resolved_settings': {'CCFLAG': ''}}, intended, lambda n: expected.get(n, -1))
    with check.assertRaisesRegex(AssertionError, 'registry mismatch'):
        check_profile(intended, intended, lambda n: -1)

    check.assertEqual(inspect_archive(archive(directory))['status'], 'PASS')
    for kwargs, message in [({'omit': 'python/ccyballs.pxd.in'}, 'missing public source'),
                            ({'extra': 'cyballs.so'}, 'build artifacts'),
                            ({'distribution': 'wrong'}, 'unexpected distribution')]:
        with check.assertRaisesRegex(AssertionError, message):
            inspect_archive(archive(directory, **kwargs))

    with check.assertRaisesRegex(AssertionError, 'external dependency profile'):
        inspect_archive(archive(directory, overrides={'Makefile_settings': 'GSLINTERNAL = 1\n'}))

    gate = Gate(argparse.Namespace(output=directory/'gate', jobs=1, mpi_command='mpiexec'))
    check.assertEqual(gate.env['PYTHONOPTIMIZE'], '0')
    output = gate.command('oracle-assertions', [sys.executable, '-c',
                          'import sys; print(sys.flags.optimize)'])
    check.assertEqual(output.strip(), '0')
    if sys.flags.optimize:
        from release_gate_case import run
        with check.assertRaisesRegex(RuntimeError, 'assertions enabled'):
            run('octree-sincos-omp', directory/'unsafe-worker', 1)
        check.assertFalse((directory/'unsafe-worker').exists())
    left, right = gate.out/'left', gate.out/'right'
    left.mkdir(); right.mkdir()
    np.savez(left/'result.npz', xi=[1.])
    np.savez(right/'result.npz', xi=[1.])
    gate.compare(left, right)
    np.savez(right/'result.npz', xi=[1.], unverified_extra=[2.])
    with check.assertRaisesRegex(AssertionError, 'product sets differ'):
        gate.compare(left, right)

    engine = 'octree-sincos-omp'
    gate.registry = {engine: 24}
    gate.report['build'] = {'id': 'expected-build'}
    metadata = dict(build={'id': 'expected-build'}, engine=engine, engine_id=24,
                    parallel=dict(estimator_ranks=1, openmp_probe_threads=1),
                    **{key: 'present' for key in ('estimator', 'bin_edges', 'coordinate_convention',
                       'weights', 'masks', 'effective_smoothing', 'opening_tolerance', 'precision')})
    case_output = gate.out/'cases'/f'{engine}-r1-t1'
    case_output.mkdir(parents=True)
    np.savez(case_output/'result.npz', xi=[1.])
    gate.command = lambda *args, **kwargs: ''
    (case_output/'result.json').write_text(json.dumps(metadata))
    gate.case(engine, 1, 1)
    for key, value in [('build', {'id': 'stale'}), ('engine', 'wrong'), ('engine_id', -1),
                       ('parallel', dict(estimator_ranks=2, openmp_probe_threads=1)),
                       ('parallel', dict(estimator_ranks=1, openmp_probe_threads=2)),
                       ('precision', '')]:
        (case_output/'result.json').write_text(json.dumps(metadata | {key: value}))
        with check.assertRaises(AssertionError):
            gate.case(engine, 1, 1)

    # The first build identity check must fail before reading any artifacts.
    def mismatched_build(label, *args, **kwargs):
        rows = {'resolved-fingerprint': {'id': 'new'}, 'native-fingerprint': {'id': 'old'},
                'cython-fingerprint': {'build': {'id': 'new'}, 'file': str(ROOT/'cyballs.so')}}
        return json.dumps(rows.get(label, {}))
    gate.command = mismatched_build
    with check.assertRaisesRegex(AssertionError, 'build mismatch'):
        gate.build()
    print('OPTIMIZATION_CHECKS_OK', sys.flags.optimize, flush=True)


@pytest.mark.parametrize('flags,optimize', [([], '0'), (['-O'], '0'), (['-OO'], '0'),
                                           ([], '1'), ([], '2')])
def test_release_validators_survive_optimization(tmp_path, flags, optimize):
    import os
    import subprocess
    completed = subprocess.run([sys.executable, *flags, str(Path(__file__).resolve()), str(tmp_path)],
                               cwd=tmp_path, env=dict(os.environ, PYTHONOPTIMIZE=optimize),
                               capture_output=True, text=True, timeout=30)
    assert completed.returncode == 0, completed.stdout + completed.stderr
    expected = max(int(optimize), 2 if '-OO' in flags else int('-O' in flags))
    assert f'OPTIMIZATION_CHECKS_OK {expected}' in completed.stdout


def test_release_scripts_have_no_removable_assert_statements():
    import ast
    # This also prevents future asserts in literal python -c command fragments.
    scripts = ('release_profile.py', 'verify_sdist_install.py', 'active_release_gate.py', 'check_installed_package.py', 'check_sdist.py',
               'release_gate_case.py', 'benchmark_contracts.py', 'benchmark_lya_exact.py',
               'workload_acceptance.py')
    for name in scripts:
        tree = ast.parse((ROOT/'scripts'/name).read_text())
        assert not any(isinstance(node, ast.Assert) for node in ast.walk(tree)), name
        for node in ast.walk(tree):
            if isinstance(node, ast.Constant) and isinstance(node.value, str):
                assert 'assert ' not in node.value, (name, node.value)


if __name__ == '__main__':
    optimization_checks(Path(sys.argv[1]))


@pytest.mark.parametrize('key', ['GSLINTERNAL', 'CFITSIOLIBON'])
def test_installed_dependency_profile_must_match_artifact(key):
    intended = {'resolved_settings': {'CCFLAG': '-DADDONS -DLYAFORESTOMP', key: '0'}}
    expected = expected_registry(intended['resolved_settings'])
    wrong = {'resolved_settings': dict(intended['resolved_settings'], **{key: '1'})}
    with pytest.raises(AssertionError, match='dependency profile mismatch'):
        check_profile(wrong, intended, lambda n: expected.get(n, -1))


def test_installed_source_identity_must_match_artifact():
    settings = {'CCFLAG': '-DADDONS -DLYAFORESTOMP'}
    with pytest.raises(AssertionError, match='source identity'):
        check_profile({'resolved_settings': settings, 'source_sha256': 'other'},
                      {'resolved_settings': settings, 'source_sha256': 'intended'}, lambda n: -1)
