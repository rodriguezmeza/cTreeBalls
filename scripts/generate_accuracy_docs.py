#!/usr/bin/env python3
"""Generate bounded acceptance tables from the same cases the numerical gate runs."""
import argparse
from pathlib import Path
from accuracy_cases import SCALAR, SHEAR

ROOT = Path(__file__).resolve().parents[1]
START = '<!-- BEGIN GENERATED ACCURACY CASES -->'
END = '<!-- END GENERATED ACCURACY CASES -->'


def rows():
    result = []
    for case in SCALAR + SHEAR:
        label = case['engine']
        core = case.get('estimator_family') == 'core'
        if case.get('compatibility'):
            label += ' (legacy-one-ball)'
        if core:
            label += ' (core estimator)'
        smoothing = ('unsupported' if not case['smoothing'] else
                     '1 arcmin' if 'shear' in label else
                     '0.0003 arcmin' if core or case.get('compatibility') else
                     '0.0003 Cartesian units')
        result.append((label, f"{case['theta']:g}", f"{case['rtol']:g}", smoothing))
    return result


def generated(root=ROOT):
    table = ['| Engine / route | Theta | Relative L2 limit | Explicit smoothing case |',
             '| --- | ---: | ---: | --- |']
    table += ['| ' + ' | '.join(row) + ' |' for row in rows()]
    md = root / 'RESOURCE_AND_NUMERICAL_CHECKS.md'
    text = md.read_text()
    if text.count(START) != 1 or text.count(END) != 1:
        raise ValueError('accuracy table markers must each occur exactly once')
    before, rest = text.split(START)
    _, after = rest.split(END)
    text = before + START + '\n' + '\n'.join(table) + '\n' + END + after
    rst = ['Bounded numerical acceptance', '============================', '',
           'Generated from ``scripts/accuracy_cases.py`` by ``generate_accuracy_docs.py``.',
           'The same definitions drive the numerical acceptance runner. These are',
           'retained regression settings, not the fast runtime defaults or a universal',
           'error guarantee. Every listed route is a required acceptance case.', '',
           '.. list-table:: Retained acceptance cases', '   :header-rows: 1',
           '   :widths: 46 10 14 30', '']
    for row in [('Engine / route', 'Theta', 'Relative L2 limit', 'Explicit smoothing case')] + rows():
        rst += ['   * - ' + row[0]] + ['     - ' + v for v in row[1:]]
    rst += ['', 'The fixtures have 512 weighted points, six logarithmic chord bins from',
            '0.02 to 1.8, orders 0 through 5, and two OpenMP threads. Each pair and',
            '3PCF product must satisfy L2(error) <= 1e-10 + rtol * L2(reference).',
            'Smoothing is explicitly enabled only for the clustered fixture.', '',
            'Modern scalar products are raw distinct-neighbor Fourier moments;',
            'the core estimator retains its separate normalized neighbor-mean definition.',
            'Shear acceptance measures raw Upsilon, not corrected Gamma. Window',
            'conditioning, arbitrary data, different opening/smoothing settings and',
            'experimental reuse need separate qualification. The default theta is one;',
            'successful execution alone remains UNMEASURED.', '']
    return {md: text, root / 'docs/accuracy_acceptance.rst': '\n'.join(rst)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    stale = []
    for path, text in generated().items():
        if not path.exists() or path.read_text() != text:
            if args.check:
                stale.append(str(path.relative_to(ROOT)))
            else:
                path.write_text(text)
    if stale:
        raise SystemExit('stale accuracy documentation: ' + ', '.join(stale))
    print('PASS: bounded accuracy documentation')


if __name__ == '__main__':
    main()
