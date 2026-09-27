#!/usr/bin/env python3
"""Compare two compiled Ly-alpha checkouts on identical retained DESI pixels.

Each timed sample uses a fresh process, warmups force native recomputation,
and baseline/candidate order alternates. Loading and getters are outside the
Run timer. All six 3D OpenMP methods are supported. Requires NumPy and cyballs
in each checkout, and Astropy for FITS input.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import resource
import statistics
import subprocess
import sys
import time

METHODS = [f'{prefix}-{stat}-omp' for prefix in ('lya', 'lya-los-tree')
           for stat in ('2pcf', '3pcf', '2pcf-3pcf')]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def worker(case):
    import numpy as np
    root = Path(case['source']).resolve()
    sys.path.insert(0, str(root))
    import cyballs
    extension = Path(cyballs.__file__).resolve()
    if extension.parent != root:
        raise RuntimeError(f'cyballs imported from wrong checkout: {extension}')
    directory = Path(case['directory']); directory.mkdir()
    with np.load(case['fixture'], allow_pickle=False) as f:
        data = {k: np.asarray(f[k]) for k in ('positions', 'delta', 'weights', 'forest_ids')}
    model = cyballs.cballs()
    try:
        model.set_forest_catalog(**data)
        for index in range(case['warmups'] + 1):
            model.set(case['parameters'] | {'rootDir': str(directory / f'run-{index}')})
            start = time.perf_counter(); cpu = time.process_time()
            model.Run()
            sample = dict(wall_seconds=time.perf_counter()-start,
                          process_cpu_seconds=time.process_time()-cpu,
                          native_timings=model.getTimings())
        np.savez_compressed(directory / 'products.npz', **model.getForestResults()['arrays'])
        sample.update(peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
                      * (1 if sys.platform == 'darwin' else 1024),
                      extension=str(extension), extension_sha256=sha(extension),
                      build=cyballs.build_info(), metadata=model.getRunMetadata(),
                      cache_info=model.getCacheInfo(), parameters=case['parameters'])
        (directory / 'sample.json').write_text(json.dumps(sample, indent=2) + '\n')
    finally:
        model.struct_cleanup()


def compare(reference, candidate):
    import numpy as np
    if set(reference) != set(candidate):
        raise AssertionError('raw product names differ')
    metrics = {}
    def measure(key, a, b):
        if a.shape != b.shape: raise AssertionError(f'{key}: shape differs')
        finite = np.isfinite(a) & np.isfinite(b)
        difference = b[finite]-a[finite]
        metrics[key] = dict(passed=bool(np.allclose(a, b, rtol=3e-11, atol=1e-10, equal_nan=True)),
            same_finite_mask=bool(np.array_equal(np.isfinite(a), np.isfinite(b))),
            relative_l2=float(np.linalg.norm(difference)/max(np.linalg.norm(a[finite]), 1e-300)),
            max_absolute=float(np.max(np.abs(difference), initial=0)))
    for key in reference: measure(key, reference[key], candidate[key])
    for prefix in ('pair', 'triple'):
        nk, dk = prefix+'_numerator', prefix+'_denominator'
        if nk not in reference: continue
        n, d = reference[nk], reference[dk]; nn, dd = candidate[nk], candidate[dk]
        metrics[prefix+'_occupancy'] = dict(passed=bool(np.array_equal(d > 0, dd > 0)))
        measure(prefix+'_ratio', np.divide(n, d, out=np.zeros_like(d), where=d > 0),
                np.divide(nn, dd, out=np.zeros_like(dd), where=dd > 0))
    return dict(passed=all(v['passed'] and v.get('same_finite_mask', True) for v in metrics.values()),
                metrics=metrics)


def main(a):
    import numpy as np
    from benchmark_scaling import load, select
    a.output.mkdir(parents=True, exist_ok=False)
    data = load(a.catalog)
    roots = {'baseline': a.baseline_root.resolve(), 'candidate': a.candidate_root.resolve()}
    report = dict(status='RUNNING', schema_version=1, catalog=str(a.catalog), catalog_sha256=sha(a.catalog),
        host=dict(platform=platform.platform(), cpu_count=os.cpu_count(), python=sys.version),
        methodology='Fresh process per sample; alternating baseline/candidate order; identical retained rows. '
        'Run wall and process CPU exclude catalog loading/getters. Peak RSS includes interpreter, selected input, '
        'warmups and native allocations. Forest indexes rebuild each Run. Exact acceptance: per-element '
        'rtol=3e-11, atol=1e-10 plus matching occupancy and finite masks. No affinity pinning is imposed.',
        selection=dict(pixel_stride=a.pixel_stride, forest_order='ascending unique ID',
                       pixel_order='original order within each forest'), fixtures=[], samples=[], summary=[])
    def save():
        (a.output/'benchmark.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    save()
    try:
        for count in a.forest_counts:
            selected, rows = select(data, count, a.pixel_stride)
            fixture = a.output/f'forests-{count}.npz'
            np.savez_compressed(fixture, **selected)
            np.save(a.output/f'forests-{count}-source-rows.npy', rows)
            report['fixtures'].append(dict(path=str(fixture), sha256=sha(fixture), pixels=len(rows), forests=count))
            for method in a.methods:
                reference = None
                for threads in a.threads:
                    label = f'f{count}-{method}-t{threads}'
                    samples = {'baseline': [], 'candidate': []}
                    for repeat in range(a.repeats):
                        for name in (('baseline', 'candidate') if repeat % 2 == 0 else ('candidate', 'baseline')):
                            params = dict(searchMethod=method, numberThreads=threads, verbose=0, verbose_log=0,
                                usePeriodic=False, useLogHist=False, rangeN=a.radius, rminHist=.001, sizeHistN=6,
                                lya2RpMax=a.radius, lya2RtMax=a.radius, lya2RpBins=a.pair_bins, lya2RtBins=a.pair_bins,
                                lya3RMax=a.radius, lya3RBins=a.radial_bins, lya3ThetaBins=a.polar_bins, lya3MuBins=a.mu_bins,
                                options='no-out-Hist,no-smooth-pivot', lya2Kernel=0, lya3Kernel=0)
                            params.update(json.loads(getattr(a, name+'_parameters')))
                            directory = a.output/f'{label}-{name}-r{repeat}'
                            case = dict(source=str(roots[name]), directory=str(directory), fixture=str(fixture),
                                        warmups=a.warmups, parameters=params)
                            with directory.with_suffix('.log').open('w') as log:
                                subprocess.run([sys.executable, str(Path(__file__).resolve()), '--worker', json.dumps(case)],
                                    stdout=log, stderr=subprocess.STDOUT, cwd=roots[name], check=True, timeout=a.timeout,
                                    env=dict(os.environ, OMP_WAIT_POLICY='PASSIVE', OMP_DYNAMIC='FALSE', PYTHONDONTWRITEBYTECODE='1'))
                            row = json.loads((directory/'sample.json').read_text())
                            with np.load(directory/'products.npz') as f: arrays = {k: f[k] for k in f.files}
                            if reference is None: reference = arrays
                            row.update(label=label, variant=name, repeat=repeat, pixels=len(rows), forests=count,
                                       directory=str(directory), accuracy=compare(reference, arrays))
                            report['samples'].append(row); samples[name].append(row); save()
                            if not row['accuracy']['passed']: raise AssertionError(f'accuracy failed: {directory}')
                    summary = dict(method=method, threads=threads, pixels=len(rows), forests=count, status='PASS')
                    for name, records in samples.items():
                        summary[name] = {key: dict(median=statistics.median(r[key] for r in records),
                            minimum=min(r[key] for r in records), maximum=max(r[key] for r in records))
                            for key in ('wall_seconds', 'process_cpu_seconds', 'peak_rss_bytes')}
                    summary['wall_speedup'] = summary['baseline']['wall_seconds']['median']/summary['candidate']['wall_seconds']['median']
                    report['summary'].append(summary); save()
                    print(label, 'PASS', f"speedup={summary['wall_speedup']:.3f}", flush=True)
        report['status'] = 'PASS'
    except BaseException as exc:
        report.update(status='INTERRUPTED' if isinstance(exc, KeyboardInterrupt) else 'FAIL', error=str(exc)); raise
    finally: save()


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--worker', help=argparse.SUPPRESS)
    for name in ('baseline-root', 'candidate-root', 'catalog', 'output'): p.add_argument('--'+name, type=Path)
    p.add_argument('--methods', nargs='+', choices=METHODS, default=METHODS)
    p.add_argument('--forest-counts', type=int, nargs='+', default=[64, 256])
    p.add_argument('--pixel-stride', type=int, default=10)
    p.add_argument('--threads', type=int, nargs='+', default=[1, 4, 16])
    p.add_argument('--radius', type=float, default=200.)
    for flag, default in (('pair-bins', 8), ('radial-bins', 4), ('polar-bins', 4), ('mu-bins', 8), ('repeats', 3), ('warmups', 1)):
        p.add_argument('--'+flag, type=int, default=default)
    p.add_argument('--timeout', type=float, default=3600)
    p.add_argument('--baseline-parameters', default='{}'); p.add_argument('--candidate-parameters', default='{}')
    a = p.parse_args()
    if a.worker: worker(json.loads(a.worker))
    else:
        if any(getattr(a, k) is None for k in ('baseline_root', 'candidate_root', 'catalog', 'output')):
            p.error('--baseline-root, --candidate-root, --catalog and --output are required')
        if min(a.forest_counts+a.threads+[a.pixel_stride, a.repeats, a.pair_bins, a.radial_bins, a.polar_bins, a.mu_bins])<1 or a.warmups<0:
            p.error('counts/bins/threads/repeats must be positive; warmups must be nonnegative')
        for k in ('catalog', 'output'): setattr(a, k, getattr(a, k).resolve())
        main(a)
