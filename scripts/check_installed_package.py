#!/usr/bin/env python3
"""Exercise the installed public profile from a fresh venv outside the checkout."""
import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import sys


def check(source_root, output):
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
    expected = expected_registry(build['resolved_settings'])
    assert len(expected) == 34, 'installation did not build the intended public profile'
    assert all(cyballs.search_method_id(name) == number for name, number in expected.items())
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
            cases.append(dict(engine=engine, status='PASS', metadata=model.getRunMetadata()))
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
    args = parser.parse_args()
    check(args.source_root, args.output.resolve())


if __name__ == '__main__':
    main()
