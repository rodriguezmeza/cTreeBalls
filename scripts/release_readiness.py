#!/usr/bin/env python3
"""Run release components and combine their evidence without accepting missing checks.

A local active gate alone is not release approval. Collection requires active
configuration/incremental checks, a clean installed package, and ASan/UBSan for
one unchanged source identity. Components may run on separate CI hosts.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import shlex
import subprocess
import sys
import time


ROOT = Path(__file__).resolve().parents[1]
from acceptance_metadata import REQUIRED, source_identity, acceptance_policy


def json_line(text):
    return next(json.loads(line) for line in text.splitlines() if line.startswith('{'))


def collect(paths, source_id, support_path=None):
    report = dict(schema_version=1, source_id=source_id, status='PASS', components={}, acceptance_policy=acceptance_policy())
    states = []
    for component, required in REQUIRED.items():
        path = paths.get(component)
        entry = dict(path=str(path) if path else None, status='INCOMPLETE')
        try:
            if path is None or not path.is_file():
                raise FileNotFoundError('required evidence is missing')
            data = json.loads(path.read_text())
            entry['sha256'] = hashlib.sha256(path.read_bytes()).hexdigest()
            if data.get('schema_version') != 1 or data.get('component') != component:
                raise ValueError('wrong component or evidence schema')
            if data.get('source_id') != source_id:
                raise ValueError('evidence belongs to a different source identity')
            if data.get('acceptance_policy') != acceptance_policy():
                raise ValueError('evidence belongs to a different acceptance policy')
            if data.get('status') == 'FAIL':
                entry.update(status='FAIL', reason=data.get('error', 'component failed'))
            elif data.get('status') != 'PASS':
                entry['reason'] = 'component is not complete: ' + str(data.get('status'))
            elif any(data.get('checks', {}).get(name) != 'PASS' for name in required):
                entry['reason'] = 'required successful checks are missing'
            elif not data.get('commands') or any(c.get('returncode') != 0 for c in data['commands']):
                entry.update(status='FAIL', reason='component contains no successful command trail or a failed command')
            else:
                entry['status'] = 'PASS'
                if component=='active':
                    report['science_acceptance']=data.get('science_acceptance','UNMEASURED')
                    entry['science_acceptance']=report['science_acceptance']
        except FileNotFoundError as exc:
            entry['reason'] = str(exc)
        except (OSError, ValueError, TypeError, AttributeError, KeyError) as exc:
            entry.update(status='FAIL', reason=str(exc))
        report['components'][component] = entry
        states.append(entry['status'])
    report['platform_support'] = dict(status='NOT_RUN')
    if support_path is not None:
        entry=dict(status='INCOMPLETE',path=str(support_path))
        try:
            data=json.loads(support_path.read_text())
            from platform_support import matrix
            expected={r['id'] for r in matrix()}
            if data.get('source_id')!=source_id or data.get('acceptance_policy')!=acceptance_policy():
                raise ValueError('platform evidence uses another source or acceptance policy')
            if data.get('status')=='PASS' and set(data.get('rows',{}))==expected and all(v.get('status')=='PASS' for v in data['rows'].values()):
                entry['status']='PASS'
            else:
                entry['status']='FAIL' if data.get('status')=='FAIL' else 'INCOMPLETE'
        except FileNotFoundError:
            entry['reason']='required platform evidence is missing'
        except (OSError,ValueError,KeyError,TypeError,AttributeError) as exc:
            entry.update(status='FAIL',reason=str(exc))
        report['platform_support']=entry;states.append(entry['status'])
    report['status'] = 'FAIL' if 'FAIL' in states else ('PASS' if all(s == 'PASS' for s in states) else 'INCOMPLETE')
    report['finished_utc'] = datetime.now(timezone.utc).isoformat()
    return report


class Unavailable(RuntimeError):
    pass


def run_component(args):
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    env = dict(os.environ, PYTHONOPTIMIZE='0', PYTHONNOUSERSITE='1')
    report = dict(schema_version=1, component=args.component, status='RUNNING',
                  source_id=source_identity(), acceptance_policy=acceptance_policy(), checks={}, commands=[],
                  started_utc=datetime.now(timezone.utc).isoformat(),
                  host=platform.platform(), python=sys.version)

    def save():
        (output/'result.json').write_text(json.dumps(report, indent=2)+'\n')

    def command(label, argv, *, timeout=2400, unavailable=False):
        argv = [str(a) for a in argv]
        log = output/(label+'.log')
        record = dict(label=label, argv=argv, cwd=str(ROOT), log=log.name, returncode=None)
        report['commands'].append(record)
        save()
        print(label, flush=True)
        start = time.monotonic()
        with log.open('w') as stream:
            try:
                result = subprocess.run(argv, cwd=ROOT, env=env, stdout=stream,
                                        stderr=subprocess.STDOUT, timeout=timeout)
                code = result.returncode
            except (subprocess.TimeoutExpired, OSError) as exc:
                stream.write(str(exc)+'\n')
                code = 124 if isinstance(exc, subprocess.TimeoutExpired) else 127
        record.update(returncode=code, seconds=round(time.monotonic()-start, 3))
        report['checks'][label] = 'PASS' if code == 0 else 'FAIL'
        save()
        if code:
            error = f'{label} failed ({code}); see {log}'
            raise Unavailable(error) if unavailable else RuntimeError(error)
        return log.read_text(errors='replace')

    save()
    try:
        if args.component == 'active':
            command('active-gate', [sys.executable, ROOT/'scripts/active_release_gate.py',
                    '--output', output/'gate', '--jobs', args.jobs, '--mpi-command', args.mpi_command])
            gate = json.loads((output/'gate/gate.json').read_text())
            if gate.get('status') != 'PASS':
                raise RuntimeError('active gate did not report PASS')
            commands = {c['label']: c for c in gate['commands'] if c.get('status') == 'PASS'}
            native = commands.get('native-unit-contracts', {})
            if 'test-parameter-file-parser' not in native.get('argv', []):
                raise RuntimeError('active gate did not execute the parameter-file parser')
            for name in ('incremental-build', 'accuracy-documentation', 'survey-science-campaign'):
                if name not in commands:
                    raise RuntimeError(f'active gate did not execute {name}')
                report['checks'][name] = 'PASS'
            campaign=json.loads((output/'gate/survey-science/campaign.json').read_text())
            if campaign.get('status')!='PASS':raise RuntimeError('science campaign is incomplete')
            report['science_acceptance']=campaign['science_acceptance']
            report['checks']['parameter-file-parser'] = 'PASS'
            if 'test-configuration-readers' not in native.get('argv',[]):
                raise RuntimeError('active gate did not execute shared configuration readers')
            report['checks']['configuration-readers']='PASS'
            report['build'] = gate['build']['id']
        elif args.component == 'package':
            command('sdist-build', [sys.executable, '-m', 'build', '--sdist', '--outdir', output/'dist', ROOT])
            archives = list((output/'dist').glob('*.tar.gz'))
            if len(archives) != 1:
                raise RuntimeError('expected exactly one source archive')
            archive = archives[0]
            report['archive_sha256'] = hashlib.sha256(archive.read_bytes()).hexdigest()
            command('sdist-metadata', [sys.executable, '-m', 'twine', 'check', '--strict', archive])
            command('sdist-policy', [sys.executable, ROOT/'scripts/check_sdist.py', archive, '--output', output/'sdist.json'])
            command('clean-install', [sys.executable, ROOT/'scripts/verify_sdist_install.py', archive,
                                     '--output', output/'artifact-check', '--jobs', args.jobs])
            verification = json.loads((output/'artifact-check/verification.json').read_text())
            if verification.get('status') != 'PASS' or not verification.get('installation'):
                raise RuntimeError('installed scientific verification is incomplete')
            report['wheel_sha256'] = verification['wheel']['sha256']
        else:
            # The caller may select a compatible compiler/runtime; instrumentation
            # is mandatory and cannot be replaced by later command-line overrides.
            make = ['make', '--no-print-directory', *args.make_arg,
                    'DEBUGCOMPILINGON=1', 'SANITIZER=address,undefined', f'PYTHON={sys.executable}']
            settings = json_line(command('sanitizer-profile', make + ['print-build-fingerprint']))
            if settings['resolved_settings'].get('SANITIZER') != 'address,undefined':
                raise RuntimeError('required sanitizer instrumentation is missing')
            env.setdefault('ASAN_OPTIONS', ('detect_leaks=0' if sys.platform == 'darwin' else 'detect_leaks=1')
                           + ':halt_on_error=1:strict_string_checks=1')
            env.setdefault('UBSAN_OPTIONS', 'halt_on_error=1:print_stacktrace=1')
            report['sanitizer_environment'] = {k: env[k] for k in ('ASAN_OPTIONS', 'UBSAN_OPTIONS')}
            report['leak_detection_supported_by_default'] = sys.platform != 'darwin'
            probe = output/'sanitizer-probe.c'
            probe.write_text('int main(void) { return 0; }\n')
            compiler = shlex.split(settings['resolved_settings']['CC'])
            binary = output/'sanitizer-probe'
            command('sanitizer-probe-build', compiler + ['-fsanitize=address,undefined',
                    '-fno-omit-frame-pointer', str(probe), '-o', str(binary)], unavailable=True)
            command('sanitizer-probe-run', [binary], unavailable=True)
            command('sanitizer-build', make + [f'-j{args.jobs}', 'cballs'])
            for label, target in [('parameter-parser', 'test-parameter-parser'),
                                  ('parameter-file-parser', 'test-parameter-file-parser'),
                                  ('standalone-parser', 'test-standalone-parser'),
                                  ('sanitizer-smoke', 'test-sanitizer-smoke'),
                                  ('physical-oracle', 'test-octree-3pcf-3d-omp')]:
                command(label, make + [target])
            report['build'] = settings['id']
            command('exact-numerical-oracle', [sys.executable,
                ROOT/'tests/python/test_exact_numerical_cli.py','--cballs',ROOT/'cballs'])
            command('physical-survey-oracle', [sys.executable,
                ROOT/'tests/python/test_octree_3pcf_3d_survey.py','--cballs',ROOT/'cballs'])
        if source_identity() != report['source_id']:
            raise RuntimeError('source changed during component validation')
        if any(report['checks'].get(name) != 'PASS' for name in REQUIRED[args.component]):
            raise RuntimeError('required checks are missing')
        report['status'] = 'PASS'
    except Unavailable as exc:
        report.update(status='UNAVAILABLE', error=str(exc))
    except Exception as exc:
        report.update(status='FAIL', error=str(exc))
    finally:
        report['finished_utc'] = datetime.now(timezone.utc).isoformat()
        save()
    print(f"{report['status']}: {output/'result.json'}", flush=True)
    return 0 if report['status'] == 'PASS' else 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='action', required=True)
    run = sub.add_parser('run')
    run.add_argument('--component', choices=REQUIRED, required=True)
    run.add_argument('--output', type=Path, required=True)
    run.add_argument('--jobs', type=int, default=4)
    run.add_argument('--mpi-command', default='mpiexec')
    run.add_argument('--make-arg', action='append', default=[], help='sanitizer compiler/runtime Make override')
    combine = sub.add_parser('collect')
    for name in REQUIRED:
        combine.add_argument('--'+name, type=Path)
    combine.add_argument('--support',type=Path,help='required complete platform coverage receipt when supplied')
    combine.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.action == 'run':
        if args.jobs < 1:
            parser.error('--jobs must be positive')
        return run_component(args)
    report = collect({name: getattr(args, name) for name in REQUIRED}, source_identity(),args.support)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2)+'\n')
    print(f"{report['status']}: {args.output}")
    return 0 if report['status'] == 'PASS' else 1


if __name__ == '__main__':
    raise SystemExit(main())
