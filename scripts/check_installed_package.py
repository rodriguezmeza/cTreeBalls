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


def check_exact_numerics(cyballs, products):
    import numpy as np
    points=np.vstack((np.eye(3),-np.ones(3)/np.sqrt(3)))
    cases=[]
    for engine in ('octree-2balls-omp','kdtree-2balls-omp','balltree-2balls-omp','octree-3pcf-3d-omp'):
        physical='3pcf-3d' in engine
        model=cyballs.cballs()
        options='no-out-Hist,no-smooth-pivot,weights-norm,no-normalize-HistZeta,no-one-ball,no-two-balls'
        options+=',compute-2pcf-3d,compute-3pcf-3d' if physical else ',edge-corrections'
        try:
            model.set(searchMethod=engine,rootDir=str(products/(engine+'-range')),
                sizeHistN=4,mChebyshev=2,rangeN=1.99,rminHist=.03,useLogHist=True,
                usePeriodic=False,theta=0,numberThreads=2,verbose=0,verbose_log=0,options=options)
            model.set_catalog(points,kappa=np.ones(4),weights=np.array([1e16,1,1,1]));model.Run()
            if physical:
                data=model.getPhysicalResults()['arrays']
                values=[data['triple_numerator'][0].sum(),data['triple_denominator'][0].sum()]
            else:
                values=[model.getHistZetaMsincos(1,1).sum(),model.getScalarWindowDiagnostics()['window_monopole'].sum()]
            np.testing.assert_allclose(values,[18e16+6]*2,rtol=3e-12,atol=0)
            model.set_catalog(points,kappa=np.ones(4),weights=np.full(4,1e155))
            try:
                model.Run()
            except Exception as error:
                if 'numerical range' not in str(error):raise
            else:raise AssertionError('installed raw overflow was accepted')
            if model.state:raise AssertionError('failed installed run has a valid state')
            model.set_catalog(points,kappa=np.ones(4),weights=np.ones(4));model.Run()
            if not model.state:raise AssertionError('installed numerical error did not recover')
            if not physical:
                model.set(options='only-2pcf,no-out-Hist,no-smooth-pivot,weights-norm,no-one-ball,no-two-balls')
                for weight in (1e-300,1e200,1e300):
                    model.set_catalog(points,kappa=np.ones(4),weights=np.full(4,weight));model.Run()
                    np.testing.assert_array_equal(model.getHistXi2pcf(),[0.,0.,0.,1.])
            cases.append(dict(engine=engine,status='PASS',oracle='N1 distinct products and N2 range/recovery'))
        finally:model.clean_all()
    return cases


def check_physical_survey(cyballs,products):
    import numpy as np
    solved=cyballs.solve_physical_survey([1.+.8/3,2.12],[1.,.8,.8],signal_lmax=1,padding=0)
    np.testing.assert_allclose(solved['zeta_corrected_legendre'],[1.,1.],rtol=2e-14)
    if solved['zeta_status']!=1:raise AssertionError('installed positive-window solve failed')
    rng=np.random.default_rng(61011);random=rng.uniform(0,1,(48,3));random[:,0]**=2;data=rng.uniform(.1,.9,(27,3))
    packets=[];model=cyballs.cballs()
    try:
        for enabled in (False,True):
            model.set(searchMethod='octree-3pcf-3d-omp',rootDir=str(products/('survey-'+str(enabled))),
                sizeHistN=2,mChebyshev=3,rangeN=1.2,rminHist=.05,useLogHist=False,usePeriodic=False,
                numberThreads=2,verbose=0,verbose_log=0,iCatalogs='1,2',surveySignalLmax=1,surveySolvePadding=1,
                options='survey-estimator-3d,compute-2pcf-3d,compute-3pcf-3d'+('' if enabled else ',no-out-Hist'))
            model.set_catalog(data,catalog=0);model.set_catalog(random,catalog=1);model.Run()
            packet=model.getResults();packets.append(packet)
            if not packet['arrays']['zeta_valid'].all():raise AssertionError('installed survey lost support')
            original=packet['arrays']['zeta_corrected_legendre'].copy()
            packet['arrays']['zeta_corrected_legendre'][0,0,0]+=1
            np.testing.assert_array_equal(model.getResults()['arrays']['zeta_corrected_legendre'],original)
            packet['arrays']['zeta_corrected_legendre']=original
            model.struct_cleanup()
        for name,value in packets[0]['arrays'].items():np.testing.assert_array_equal(value,packets[1]['arrays'][name])
        if cyballs.qualification_report(*packets,rtol=1e-12)['status']!='QUALIFIED':
            raise AssertionError('installed corrected survey qualification failed')
    finally:model.clean_all()
    return dict(engine='octree-3pcf-3d-omp',status='PASS',oracle='positive window, owned correction, output independence, reused-object identity and corrected qualification',arrays=sorted(packets[0]['arrays']))


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
    cases.extend(check_exact_numerics(cyballs,products))
    cases.append(check_physical_survey(cyballs,products))
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
