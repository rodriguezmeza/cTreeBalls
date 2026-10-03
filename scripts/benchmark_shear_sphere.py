#!/usr/bin/env python3
"""Measure spherical shear with reproducible data and per-observable errors.

Each repeat uses a cold native context in this process. Run competing builds
in separate processes with --module; do not benchmark concurrent workers.
Approximation is never qualified solely by a timing or a phase parameter.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import resource
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def compare(reference, candidate):
    """Scaled norms avoid squaring large finite values before normalization."""
    if reference.shape != candidate.shape:
        raise ValueError('observable shapes differ')
    mask = np.isfinite(reference)
    same = bool(np.array_equal(mask, np.isfinite(candidate)))
    a, b = reference[mask], candidate[mask]
    if not same or not len(a):
        return dict(same_finite_mask=same, relative_l2=None, max_absolute=None,
                    passed=False)
    scale = max(float(np.max(np.abs(a), initial=0)),
                float(np.max(np.abs(b), initial=0)), 1e-300)
    delta = np.abs(b/scale - a/scale)
    error = float(np.linalg.norm(delta))
    norm = float(np.linalg.norm(a/scale))
    floor = max(1e-12/scale, 1e-6*float(np.max(np.abs(a/scale), initial=0)))
    strong = np.abs(a/scale) > floor
    return dict(same_finite_mask=True, nonfinite=int((~mask).sum()),
                relative_l2=error/max(norm, 1e-300),
                max_absolute=float(np.max(delta, initial=0)*scale),
                weak_max_absolute=float(np.max(delta[~strong], initial=0)*scale),
                max_bin_relative=float(np.max(delta[strong]/np.abs(a[strong]/scale), initial=0)),
                passed=bool(np.isfinite(error) and error <= .02*norm+1e-12/scale))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--module', type=Path, default=ROOT,
                        help='directory containing the extension to measure')
    parser.add_argument('--output', type=Path, required=True, help='new output directory')
    parser.add_argument('--engine', default='octree-shear-sphere-2balls-omp',
                        choices=('octree-shear-sphere-2balls-omp',
                                 'kdtree-shear-sphere-2balls-omp',
                                 'balltree-shear-sphere-2balls-omp'))
    parser.add_argument('--n', type=int, default=8192)
    parser.add_argument('--geometry', choices=('uniform', 'clustered'), default='uniform')
    parser.add_argument('--field', choices=('coherent', 'random'), default='coherent')
    parser.add_argument('--seed', type=int, default=84261001)
    parser.add_argument('--masked', action='store_true')
    parser.add_argument('--order', choices=('2', '3', 'both'), default='both')
    parser.add_argument('--threads', type=int, default=1)
    parser.add_argument('--theta', type=float, default=.1)
    parser.add_argument('--exact', action='store_true')
    parser.add_argument('--repeat', type=int, default=3)
    parser.add_argument('--bins', type=int, default=6)
    parser.add_argument('--modes', type=int, default=5)
    parser.add_argument('--rmin', type=float, default=.03)
    parser.add_argument('--rmax', type=float, default=1.8)
    parser.add_argument('--smooth', type=float, help='rsmooth in arcminutes; approximate pivots')
    parser.add_argument('--nsmooth', type=int)
    parser.add_argument('--pivot-reuse', action='store_true')
    parser.add_argument('--reference', type=Path, help='NPZ from a matching exact run')
    args = parser.parse_args()
    if min(args.n, args.repeat, args.threads, args.bins) < 1:
        parser.error('n, repeat, threads, and bins must be positive')
    if args.exact and (args.smooth is not None or args.pivot_reuse):
        parser.error('--exact excludes smoothing and pivot reuse')
    args.module = args.module.resolve()
    args.output = args.output.resolve()
    sys.path.insert(0, str(args.module))
    import cyballs
    if Path(cyballs.__file__).parent != args.module:
        raise RuntimeError('loaded extension differs from --module')
    args.output.mkdir(parents=True, exist_ok=False)
    rng = np.random.default_rng(args.seed)
    if args.geometry == 'clustered':
        centers = rng.normal(size=(64, 3))
        centers /= np.linalg.norm(centers, axis=1)[:, None]
        positions = centers[np.arange(args.n) % 64] + rng.normal(scale=.002, size=(args.n, 3))
    else:
        positions = rng.normal(size=(args.n, 3))
    positions /= np.linalg.norm(positions, axis=1)[:, None]
    lon = np.arctan2(positions[:, 1], positions[:, 0])
    lat = np.arcsin(positions[:, 2])
    gamma = .006*(rng.normal(size=args.n) + 1j*rng.normal(size=args.n))
    if args.field == 'coherent':
        gamma += .02*np.exp(2j*(lon+.21*np.sin(lat)))
    weights = rng.uniform(.25, 2., args.n)
    mask = np.arange(args.n) % 11 != 0 if args.masked else np.ones(args.n, dtype=bool)
    arrays = dict(positions=positions, gamma=gamma, weights=weights, mask=mask)
    np.savez_compressed(args.output/'fixture.npz', **arrays)
    options = ['no-out-Hist', 'no-smooth-pivot', 'no-normalize-HistZeta']
    if args.order != '2':
        options += ['edge-corrections']
    if args.exact:
        options += ['no-one-ball', 'no-two-balls']
    if args.order != 'both':
        options += ['only-'+args.order+'pcf']
    if args.masked:
        options += ['read-mask']
    if args.pivot_reuse:
        options += ['shear-pivot-reuse']
    params = dict(searchMethod=args.engine, usePeriodic=False,
                  useLogHist=True, rminHist=args.rmin, rangeN=args.rmax,
                  sizeHistN=args.bins, sizeHistPhi=32, mChebyshev=args.modes,
                  theta=args.theta, numberThreads=args.threads, verbose=0, verbose_log=0)
    if args.nsmooth is not None:
        params['nsmooth'] = args.nsmooth
    if args.smooth is not None:
        options.remove('no-smooth-pivot')
        options += ['smooth-pivot']
        params['rsmooth'] = str(args.smooth)
    params['options'] = ','.join(options)
    record = dict(schema_version=1, arguments={k: str(v) if isinstance(v, Path) else v
                                               for k, v in vars(args).items()},
                  parameters=params, build=cyballs.build_info(), host=platform.platform(),
                  environment={k: os.environ.get(k) for k in
                               ('CBALLS_SHEAR_PIVOT_TOL', 'CBALLS_SHEAR_BIN_THETA', 'CBALLS_SHEAR_PROFILE', 'OMP_PROC_BIND')},
                  fixture_sha256=hashlib.sha256(b''.join(x.tobytes() for x in arrays.values())).hexdigest(),
                  timed_scope='Run through MainLoop: initialization, tree, search, correction; cold context each repeat',
                  reference_criterion='each observable L2 <= .02*reference L2 + 1e-12, identical finite mask',
                  samples=[])
    reference = None
    if args.reference:
        with np.load(args.reference) as archive:
            reference = {k: archive[k].copy() for k in archive.files}
    for iteration in range(args.repeat):
        model = cyballs.cballs()
        model.set(dict(params, rootDir=str(args.output/f'native-{iteration}')))
        model.set_catalog(positions, gamma1=gamma.real, gamma2=gamma.imag,
                          weights=weights, **({'mask': mask} if args.masked else {}))
        wall, cpu = time.perf_counter(), time.process_time()
        try:
            model.Run(level=['MainLoop'])
            wall, cpu = time.perf_counter()-wall, time.process_time()-cpu
            values = {}
            if args.order != '3':
                values.update(xi_plus=model.getShearXiPlus(), xi_minus=model.getShearXiMinus(),
                              xi_weight=model.getShearXiWeight())
            if args.order != '2':
                values.update(upsilon=model.getShearUpsilonXMultipoles(),
                              window=model.getShearWindowMultipoles(),
                              corrected=model.getShearGammaXMultipoles())
            np.savez_compressed(args.output/f'products-{iteration}.npz', **values)
            sample = dict(wall_seconds=wall, process_cpu_seconds=cpu,
                          run_metadata=model.getRunMetadata())
            if reference is not None:
                if set(values) != set(reference):
                    raise ValueError('reference observable set differs')
                sample['accuracy'] = {k: compare(reference[k], values[k]) for k in values}
                sample['qualification'] = ('QUALIFIED' if all(v['passed'] for v in sample['accuracy'].values())
                                           else 'REJECTED')
            record['samples'].append(sample)
            peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
            record['process_peak_rss_bytes'] = int(peak if sys.platform == 'darwin' else peak*1024)
            (args.output/'timing.json').write_text(json.dumps(record, indent=2, allow_nan=False)+'\n')
            print(f'{iteration}: wall={wall:.6f}s cpu={cpu:.6f}s '
                  f'{sample.get("qualification", "UNQUALIFIED: no reference supplied")}', flush=True)
        finally:
            model.struct_cleanup()


if __name__ == '__main__':
    main()
