#!/usr/bin/env python3
"""Exercise the installed public profile from a fresh venv outside the checkout."""
import argparse
import hashlib
import importlib.metadata
import importlib.util
import json
from pathlib import Path
import sys
import subprocess


def check_profile(installed_build, intended_build, lookup):
    from capabilities_generated import expected_registry, ENGINES
    for key in ('GSLINTERNAL', 'CFITSIOLIBON'):
        if key in intended_build['resolved_settings'] and (
                installed_build['resolved_settings'].get(key) != intended_build['resolved_settings'][key]):
            raise AssertionError(f'installation dependency profile mismatch: {key}')
    if intended_build.get('source_sha256') and installed_build.get('source_sha256') != intended_build['source_sha256']:
        raise AssertionError('installed source identity differs from the tested archive')
    expected = expected_registry(intended_build['resolved_settings'])
    actual = expected_registry(installed_build['resolved_settings'])
    if not (actual == expected):
        raise AssertionError(f'installation profile mismatch: missing={set(expected) - set(actual)}, extra={set(actual) - set(expected)}')
    if not (expected):
        raise AssertionError('intended profile is empty')
    # Check absent canonical names too, so a stale extension cannot hide extras.
    for name in ENGINES:
        if not (lookup(name) == expected.get(name, -1)):
            raise AssertionError(f'installed registry mismatch: {name}')
    return expected


def check(source_root, output, expected_build=None):
    import numpy as np
    import cyballs
    from active_release_gate import expected_registry
    source_root = source_root.resolve()
    module = Path(cyballs.__file__).resolve()
    if not (not Path.cwd().resolve().is_relative_to(source_root)):
        raise AssertionError('run outside the checkout')
    if not (not module.is_relative_to(source_root)):
        raise AssertionError('imported cyballs from the checkout')
    if not (sys.prefix != sys.base_prefix):
        raise AssertionError('use a fresh virtual environment')
    if not (module.is_relative_to(Path(sys.prefix).resolve())):
        raise AssertionError('module is outside the install environment')
    build = cyballs.build_info()
    if expected_build is None:
        text = subprocess.check_output(['make', '--no-print-directory', 'print-build-fingerprint',
                                        f'PYTHON={sys.executable}'], cwd=source_root, text=True)
        intended = next(json.loads(line) for line in text.splitlines() if line.startswith('{'))
    else:
        intended = json.loads(Path(expected_build).read_text())
    expected = check_profile(build, intended, cyballs.search_method_id)
    distribution = importlib.metadata.distribution('cyballs')
    if not (distribution.metadata['Name'] == 'cyballs'):
        raise AssertionError('installed distribution name is not cyballs')
    products = output.resolve().with_suffix('.products')
    products.mkdir(parents=True, exist_ok=False)
    positions = np.array([[10., 0., 0.], [9., 2., 0.], [8., 0., 3.]])
    delta = np.array([2., 3., 5.])
    weights = np.array([1., 2., 4.])
    cases = []
    for order in ('2pcf', '3pcf', '2pcf-3pcf'):
        engine = f'lya-los-tree-{order}-omp'
        directory = products/order
        model = cyballs.cballs()
        try:
            model.set(dict(searchMethod=engine, rootDir=str(directory), numberThreads=2,
                           verbose=0, verbose_log=0, iCatalogs='1', usePeriodic=False,
                           useLogHist=False, rangeN=30., rminHist=.1, sizeHistN=4,
                           lya2RpMax=30., lya2RtMax=30., lya2RpBins=5, lya2RtBins=6,
                           lya3RMax=30., lya3RBins=4, lya3ThetaBins=5, lya3MuBins=6,
                           options=''))
            model.set_forest_catalog(positions, delta, weights, np.arange(3, dtype=np.int64))
            model.Run(level=['MainLoop'])
            for statistic, filename, columns, totals in (
                ('2pcf', 'histXi2pcf_lya.txt', (5, 6), (172., 14.)),
                ('3pcf', 'histZetaM_lya5d.txt', (11, 12), (1440., 48.))):
                if statistic in order:
                    data = np.atleast_2d(np.loadtxt(directory/filename))
                    np.testing.assert_allclose(data[:, columns].sum(axis=0), totals,
                                               rtol=2e-13, atol=2e-13)
            native = model.getForestResults()['arrays']
            for key, total in [('pair_numerator',172.), ('pair_denominator',14.),
                               ('triple_numerator',1440.), ('triple_denominator',48.)]:
                if key in native:
                    np.testing.assert_allclose(native[key].sum(), total, rtol=2e-13)
            plan = cyballs.resource_plan(engine, len(positions), 2)
            if not (plan['known_total_bytes'] > 0):
                raise AssertionError('resource plan has no known allocation bytes')
            if not (model.getAllocationInfo()['retained_result_bytes'] > 0):
                raise AssertionError('model retained no result bytes')
            cases.append(dict(engine=engine, status='PASS', arrays=sorted(native),
                              resource_plan=plan, metadata=model.getRunMetadata()))
        finally:
            model.struct_cleanup()
    # Independent ordered-triple oracle also covers the shared angular solver.
    specification = importlib.util.spec_from_file_location('installed_edge_oracle',
        source_root/'tests/python/test_two_ball_edge_corrections.py')
    oracle = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(oracle)
    data = oracle.catalog(count=32)
    signal, window = oracle.brute_force(data)
    reference = oracle.edge_solution(signal, window)
    for engine in ('kdtree-2balls-omp', 'balltree-2balls-omp', 'octree-2balls-omp'):
        model = cyballs.cballs()  # Construction verifies C/Cython ABI and build identity.
        try:
            model.set_catalog(data[0], kappa=data[1], weights=data[2])
            model.set(searchMethod=engine, rootDir=str(products/engine), numberThreads=2,
                verbose=0, verbose_log=0, rangeN=oracle.RMAX, rminHist=oracle.RMIN,
                sizeHistN=oracle.BINS, mChebyshev=oracle.MMAX, sizeHistPhi=8,
                usePeriodic=False, useLogHist=False, nsmooth=2, theta=0,
                options='KKKCorrelation,weights-norm,only-3pcf,no-smooth-pivot,'
                        'edge-corrections,no-normalize-HistZeta,no-out-Hist')
            model.Run()
            actual = np.array([model.getHistZetaM_EE_complex(m+1) for m in range(oracle.MMAX+1)])
            np.testing.assert_allclose(actual, reference, rtol=2e-10, atol=2e-10)
            cases.append(dict(engine=engine, status='PASS', oracle='independent ordered triples + NumPy solve',
                              abi_sizes=model.abi_sizes(), build_id=model.getRunMetadata()['build']['id']))
        finally:
            model.struct_cleanup()
    result = dict(status='PASS', distribution=distribution.metadata['Name'],
                  version=distribution.version, python=sys.executable, prefix=sys.prefix,
                  working_directory=str(Path.cwd()), module=str(module),
                  module_sha256=hashlib.sha256(module.read_bytes()).hexdigest(),
                  build=build, expected_registry=expected, cases=cases)
    output.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    print(f'PASS: installed cyballs {distribution.version}: {module}')
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--expected-build', type=Path, help='retained gate build-fingerprint.json; otherwise resolve saved Makefiles')
    args = parser.parse_args()
    check(args.source_root, args.output.resolve(), args.expected_build)


if __name__ == '__main__':
    main()
