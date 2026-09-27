#!/usr/bin/env python3
"""Exercise the installed public profile from a fresh venv outside the checkout."""
import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import sys
import subprocess


def check_profile(installed_build, intended_build, lookup):
    from capabilities_generated import expected_registry, ENGINES
    expected = expected_registry(intended_build['resolved_settings'])
    actual = expected_registry(installed_build['resolved_settings'])
    assert actual == expected, f'installation profile mismatch: missing={set(expected)-set(actual)}, extra={set(actual)-set(expected)}'
    assert expected, 'intended profile is empty'
    # Check absent canonical names too, so a stale extension cannot hide extras.
    for name in ENGINES:
        assert lookup(name) == expected.get(name, -1), f'installed registry mismatch: {name}'
    return expected


def check(source_root, output, expected_build=None):
    import numpy as np
    import cyballs
    from active_release_gate import expected_registry
    source_root = source_root.resolve()
    module = Path(cyballs.__file__).resolve()
    assert not Path.cwd().resolve().is_relative_to(source_root), 'run outside the checkout'
    assert not module.is_relative_to(source_root), 'imported cyballs from the checkout'
    assert sys.prefix != sys.base_prefix, 'use a fresh virtual environment'
    assert module.is_relative_to(Path(sys.prefix).resolve()), 'module is outside the install environment'
    build = cyballs.build_info()
    if expected_build is None:
        text = subprocess.check_output(['make', '--no-print-directory', 'print-build-fingerprint',
                                        f'PYTHON={sys.executable}'], cwd=source_root, text=True)
        intended = next(json.loads(line) for line in text.splitlines() if line.startswith('{'))
    else:
        intended = json.loads(Path(expected_build).read_text())
    expected = check_profile(build, intended, cyballs.search_method_id)
    distribution = importlib.metadata.distribution('cyballs')
    assert distribution.metadata['Name'] == 'cyballs'
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
            assert plan['known_total_bytes'] > 0
            assert model.getAllocationInfo()['retained_result_bytes'] > 0
            cases.append(dict(engine=engine, status='PASS', arrays=sorted(native),
                              resource_plan=plan, metadata=model.getRunMetadata()))
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
