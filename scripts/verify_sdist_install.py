#!/usr/bin/env python3
"""Build native + isolated wheel from a checked sdist, then test a clean install.

No dependency-profile overrides are injected. The archive's own Makefiles must
select its external libraries. Native toolchain/dependency discovery remains
host configuration; PEP 517 isolates Python build dependencies, not the OS.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import time

from check_sdist import inspect_archive


def extract(archive, directory):
    # inspect_archive has already rejected links, special files, absolute paths,
    # duplicate file entries and traversal. Extract only regular files ourselves.
    directory.mkdir()
    with tarfile.open(archive, 'r:gz') as source:
        for member in source.getmembers():
            target = directory / member.name
            if member.isdir():
                target.mkdir(parents=True, exist_ok=True)
            else:
                target.parent.mkdir(parents=True, exist_ok=True)
                with source.extractfile(member) as src, target.open('wb') as dst:
                    shutil.copyfileobj(src, dst)
                target.chmod(member.mode & 0o777)
    roots = list(directory.iterdir())
    if len(roots) != 1:
        raise AssertionError('archive extraction did not produce one root')
    return roots[0]


def json_line(text):
    return next(json.loads(line) for line in text.splitlines() if line.startswith('{'))


def verify(archive, output, python=sys.executable, jobs=4):
    archive = archive.resolve()
    output = output.resolve()
    if output.exists():
        raise RuntimeError(f'use a new evidence directory: {output}')
    output.mkdir(parents=True)
    report = {'status': 'RUNNING', 'steps': [], 'archive': inspect_archive(archive)}
    env = dict(os.environ)
    for key in list(env):
        if key in {'PYTHONPATH', 'PYTHONHOME', 'MAKEFLAGS', 'MFLAGS', 'MAKEOVERRIDES'} or key.startswith('CBALLS_'):
            env.pop(key)
    env.update(PYTHONNOUSERSITE='1', PYTHONOPTIMIZE='0', MAKEFLAGS=f'-j{jobs}')
    outside = output / 'outside'
    outside.mkdir()

    def save():
        (output/'verification.json').write_text(json.dumps(report, indent=2)+'\n')

    def run(label, argv, cwd=outside, timeout=1800):
        argv = [str(x) for x in argv]
        print(f'{label} ...', flush=True)
        started = time.monotonic()
        log = output / (label+'.log')
        with log.open('w') as stream:
            try:
                result = subprocess.run(argv, cwd=cwd, env=env, stdout=stream,
                                        stderr=subprocess.STDOUT, timeout=timeout)
                code = result.returncode
            except subprocess.TimeoutExpired:
                code = 124
        report['steps'].append(dict(name=label, argv=argv, cwd=str(cwd),
            returncode=code, seconds=round(time.monotonic()-started, 3), log=str(log)))
        save()
        if code:
            raise RuntimeError(f'{label} failed ({code}); see {log}')
        return log.read_text(errors='replace')

    try:
        native = extract(archive, output/'native-source')
        wheel_source = extract(archive, output/'wheel-source')
        run('create-native-env', [python, '-m', 'venv', output/'native-env'])
        native_python = output/'native-env/bin/python'
        run('native-python-dependencies', [native_python, '-m', 'pip', 'install',
                                            '-r', native/'requirements/build.txt'])
        run('native-default-build', ['make', 'cballs', f'PYTHON={native_python}'], cwd=native)
        intended = json_line(run('resolved-profile', ['make', '--no-print-directory',
            'print-build-fingerprint', f'PYTHON={native_python}'], cwd=native))
        for key in ('GSLINTERNAL', 'CFITSIOLIBON'):
            if intended['resolved_settings'].get(key) != '0':
                raise AssertionError(f'archive resolved {key} to a bundled dependency')
        native_build = json_line(run('native-profile', [native/'cballs',
            'options=build-fingerprint', f'rootDir={outside/"native-info"}',
            'verbose=0', 'verbose_log=0']))
        if native_build != intended:
            raise AssertionError('native binary differs from the resolved archive profile')
        (output/'expected-build.json').write_text(json.dumps(intended, indent=2)+'\n')
        report['native_build'] = intended
        run('native-window-oracle', [native_python,
            native/'tests/python/test_two_ball_edge_corrections.py',
            '--cballs', native/'cballs', '--engine', 'kdtree-2balls-omp',
            '--engine', 'balltree-2balls-omp', '--engine', 'octree-2balls-omp'])
        run('isolated-wheel', [python, '-m', 'build', '--wheel', '--outdir',
                                output/'wheels', wheel_source])
        wheels = list((output/'wheels').glob('*.whl'))
        if len(wheels) != 1:
            raise AssertionError(f'expected one built wheel: {wheels}')
        wheel = wheels[0]
        report['wheel'] = dict(path=str(wheel), sha256=hashlib.sha256(wheel.read_bytes()).hexdigest())
        run('create-install-env', [python, '-m', 'venv', output/'install-env'])
        installed_python = output/'install-env/bin/python'
        run('install-wheel', [installed_python, '-m', 'pip', 'install', wheel])
        run('installed-oracles', [installed_python, wheel_source/'scripts/check_installed_package.py',
            '--source-root', wheel_source, '--expected-build', output/'expected-build.json',
            '--output', output/'installation.json'])
        # Retain the pre-existing installed-package lifecycle regressions too.
        env['CBALLS_FAILURE_REPEATS'] = '60'
        for name in ('test_p0_cython_instances', 'test_p1_cython', 'test_p2_cython',
                     'test_cython_in_memory_catalog', 'test_p3_cython_startup'):
            run('installed-'+name, [installed_python, wheel_source/'tests/python'/(name+'.py')])
        env.pop('CBALLS_FAILURE_REPEATS')
        report['installation'] = json.loads((output/'installation.json').read_text())
        report['status'] = 'PASS'
    except Exception as exc:
        report['status'] = 'FAIL'
        report['error'] = str(exc)
        raise
    finally:
        save()
    print(f'PASS: checked sdist, native default build, isolated wheel and installed oracles: {output}', flush=True)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('archive', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--python', default=sys.executable)
    parser.add_argument('--jobs', type=int, default=4)
    args = parser.parse_args()
    if args.jobs < 1:
        parser.error('--jobs must be positive')
    verify(args.archive, args.output, args.python, args.jobs)


if __name__ == '__main__':
    main()
