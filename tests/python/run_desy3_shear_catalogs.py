#!/usr/bin/env python3
"""Run one Slurm array item, or a serial sequence of DES-Y3 shear catalog jobs."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
from shear_fits_catalog import integer_selection, discover_des_catalogs as discover, shear_convention as convention
SOURCE_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(SOURCE_ROOT))


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--catalog-root', type=Path, required=True)
    p.add_argument('--realizations', default='2', help='comma list or inclusive ranges, e.g. 1:108')
    p.add_argument('--tomobins', default='1:4')
    p.add_argument('--regions', default='1:4')
    p.add_argument('--task-index', type=int, help='zero-based row in deterministic realization/bin/region order')
    p.add_argument('--outdir', type=Path, default=Path('results/desy3-shear'))
    p.add_argument('--list', action='store_true', help='validate headers and print array inventory without importing cyballs')
    p.add_argument('--resume', action='store_true', help='skip only matching completed inputs/configuration/code')
    args, extra = p.parse_known_args(argv)
    if extra and extra[0] == '--':
        extra.pop(0)
    if any(x.split('=')[0] in ('--fits', '--catalog-npz', '--synthetic-nbody', '--outdir', '--continue-on-error', '--fits-format', '--list-engines', '--save-catalog-npz') for x in extra):
        p.error('input/output and failure handling are owned by this batch driver')
    jobs = discover(args.catalog_root, integer_selection(args.realizations),
                    integer_selection(args.tomobins, 1, 4), integer_selection(args.regions, 1, 4))
    convention_parser = argparse.ArgumentParser(add_help=False)
    convention_parser.add_argument('--des-shear-convention', default='auto')
    convention_args, _ = convention_parser.parse_known_args(extra)
    for job in jobs:
        convention(job['header'], convention_args.des_shear_convention)
    if args.list:
        print('task\trealization\ttomobin\tregion\trows\tpath')
        for i, job in enumerate(jobs):
            print(f"{i}\t{job['realization']}\t{job['tomobin']}\t{job['region']}\t{job['rows']}\t{job['path']}")
        return 0
    if args.task_index is not None:
        if not 0 <= args.task_index < len(jobs):
            p.error(f'--task-index must be in 0..{len(jobs)-1}')
        jobs = [jobs[args.task_index]]
    here = Path(__file__).resolve().parent
    code = {name: hashlib.sha256((here/name).read_bytes()).hexdigest() for name in
            ('shear_corr_all_engines.py', 'shear_fits_catalog.py', 'shear_products.py', 'run_desy3_shear_catalogs.py')}
    defaults = ['--binning', 'sofia-fig1', '--statistics', 'both',
                '--no-plots', '--max-points', '0', '--multipoles', '7', '--phi-bins', '32']
    for job in jobs:
        target = args.outdir.resolve()/f"r{job['realization']:03d}"/f"bin{job['tomobin']}"/f"region{job['region']}"
        marker = target/'COMPLETE.json'
        stat = Path(job['path']).stat()
        config = dict(input=job, input_size=stat.st_size, input_mtime_ns=stat.st_mtime_ns,
                      arguments=defaults+extra, code=code,
                      native_source=str(SOURCE_ROOT),
                      native_executable=os.environ.get('CTREEBALLS_CBALLS'),
                      shear_pivot_budget=os.environ.get('CBALLS_SHEAR_PIVOT_TOL', '0.1'),
                      python=str(Path(sys.executable).resolve()))
        # Native extension content is included so a rebuilt library invalidates resume.
        import importlib.util
        spec = importlib.util.find_spec('cyballs')
        if spec is None or not spec.origin:
            raise RuntimeError('cyballs is unavailable; activate the matching cTreeBalls environment')
        config['extension'] = dict(path=spec.origin, sha256=hashlib.sha256(Path(spec.origin).read_bytes()).hexdigest())
        digest = hashlib.sha256(json.dumps(config, sort_keys=True).encode()).hexdigest()
        if marker.exists():
            old = json.loads(marker.read_text())
            intact = all(Path(name).is_file() and hashlib.sha256(Path(name).read_bytes()).hexdigest()==sha
                         for name, sha in old.get('products', {}).items())
            if args.resume and old.get('fingerprint') == digest and intact and old.get('products'):
                print(f'SKIP completed {target}', flush=True)
                continue
            raise ValueError(f'{target}: completed output exists; use matching --resume or a new --outdir')
        if (target/'request.json').exists():
            previous = json.loads((target/'request.json').read_text())
            previous.pop('command', None)
            if previous != config:
                raise ValueError(f'{target}: unfinished output has another configuration; use a new --outdir')
        target.mkdir(parents=True, exist_ok=True)
        command = [sys.executable, '-u', str(here/'shear_corr_all_engines.py'),
                   '--fits', job['path'], '--fits-format', 'desy3', '--outdir', str(target), *defaults, *extra]
        (target/'request.json').write_text(json.dumps(dict(command=command, **config), indent=2)+'\n')
        print(f"RUN r{job['realization']:03d} bin{job['tomobin']} region{job['region']}", flush=True)
        # Isolate native state and reclaim catalog/tree memory between jobs.
        with (target/'run.log').open('w') as stream:
            status = subprocess.call(command, stdout=stream, stderr=subprocess.STDOUT)
        if status:
            raise RuntimeError(f'{target}/run.log: driver failed with status {status}')
        summary = json.loads((target/'summary.json').read_text())
        if summary['failures'] or not summary['engines']:
            raise RuntimeError(f'{target}: incomplete engine suite')
        products = [target/'summary.json']
        for engine in summary['engines']:
            products += [target/f'{engine}_histograms.npz', target/f'{engine}_histograms.json']
        payload = dict(fingerprint=digest, products={str(f): hashlib.sha256(f.read_bytes()).hexdigest() for f in products})
        temporary = target/'COMPLETE.json.tmp'; temporary.write_text(json.dumps(payload, indent=2)+'\n')
        temporary.replace(marker)
        print(f'DONE {target}', flush=True)
    return 0


if __name__ == '__main__':
    try:
        raise SystemExit(main())
    except (ValueError, OSError, RuntimeError) as exc:
        print(f'ERROR: {exc}', file=sys.stderr); raise SystemExit(1)
