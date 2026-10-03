#!/usr/bin/env python3
"""Native ASCII and recorded Gadget/FITS regressions. CLI checks use only the standard library."""
import argparse
import math
from pathlib import Path
import struct
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
FIXTURES = ROOT / 'tests/fixtures/io_stabilization'
BINARY = ROOT / 'cballs'
WITH_CYTHON = False
POSITIONS = [(1, 2, 3), (2, 3, 4), (3, 4, 5), (4, 5, 6),
             (5, 6, 7), (6, 7, 8), (7, 8, 9), (8, 1, 2)]
OPTIONS = 'stop,no-smooth-pivot,no-out-Hist'


def snapshot(path, positions, box=10.0):
    """Make boundary variants of the recorded native Gadget format-1 data."""
    counts = [0, len(positions), 0, 0, 0, 0]
    header = struct.pack('<6i6d2d2i6i2i4d96s', *counts, *([1.0] * 6),
                         1.0, 0.0, 0, 0, *counts, 0, 1, box, .3, .7, .7, b'')
    coords = struct.pack('<' + 'f' * (3 * len(positions)),
                         *(x for point in positions for x in point))
    path.write_bytes(struct.pack('<i', 256) + header +
                     struct.pack('<ii', 256, len(coords)) + coords +
                     struct.pack('<i', len(coords)))


def read_text_catalog(path):
    lines = path.read_text().splitlines()
    header = [float(x) for x in lines[1].lstrip('#').split()]
    rows = [[float(x) for x in line.split()] for line in lines if line and not line.startswith('#')]
    return header, rows


def read_fits_export(path):
    """Decode this writer's primary HDU + four-double table independently."""
    data = path.read_bytes()
    offset = 0
    headers = []
    for _ in range(2):
        cards = {}
        while True:
            card = data[offset:offset + 80].decode('ascii')
            if len(card) != 80:
                raise AssertionError('truncated FITS header')
            offset += 80
            key = card[:8].strip()
            if key == 'END':
                break
            if card[8:10] == '= ':
                cards[key] = card[10:].split('/')[0].strip().strip("'").strip()
        offset = ((offset + 2879) // 2880) * 2880
        headers.append(cards)
    assert headers[0]['SIMPLE'] == 'T' and headers[0]['NAXIS'] == '0'
    h = headers[1]
    assert h['XTENSION'] == 'BINTABLE' and int(h['NAXIS1']) == 32
    assert int(h['TFIELDS']) == 4
    assert [h[f'TTYPE{i}'] for i in range(1, 5)] == ['X', 'Y', 'Z', 'KAPPA']
    assert all(h[f'TFORM{i}'] == '1D' for i in range(1, 5))
    return [struct.unpack_from('>4d', data, offset + 32 * i)
            for i in range(int(h['NAXIS2']))]


ASCII_ROWS = ['1 0 0 -1', '0 1 0 2', '0 0 1 3']
ASCII_ENGINES = ['kdtree-2balls-omp', 'balltree-2balls-omp', 'octree-2balls-omp']


def ascii_catalog(fmt):
    rows = [row + (' 1 1' if fmt == 'columns-ascii-all' else '') for row in ASCII_ROWS]
    return '# native ASCII regression\n# 3 3 4 4 4\n' + '\n'.join(rows) + '\n'


def ascii_invalid_catalogs(fmt):
    good = ascii_catalog(fmt)
    lines = good.splitlines()
    # Every header/body token boundary can be a truncated file, including a
    # missing scalar, weight or mask after otherwise complete coordinates.
    yield 'empty', '', 'header comment'
    tokens = '\n'.join(lines[1:]).split()
    for end in range(len(tokens)):
        yield f'truncated-{end}', lines[0] + '\n' + ' '.join(tokens[:end]), 'unexpected end of file'
    for index, field, bad in [(0, 'marker', '!'), (1, 'nbody', '0'), (1, 'nbody', '-3'),
                              (1, 'nbody', '9' * 80), (1, 'nbody', '3oops'),
                              (2, 'ndim', '2'), (3, 'box dimension', 'nan')]:
        header = lines[1].split()
        header[index] = bad
        yield f'header-{field}-{bad}', '\n'.join([lines[0], ' '.join(header), *lines[2:]]), 'header ' + field
    fields = [(0, 'position'), (1, 'position'), (2, 'position'), (3, 'kappa')]
    if fmt == 'columns-ascii-all':
        fields.append((4, 'weight'))
    for index, field in fields:
        for bad in ('nan', 'inf', '-inf', '1e309', '2oops', 'oops', '9' * 300):
            row = lines[-1].split()
            row[index] = bad
            yield f'{field}-{index}-{bad[:12]}', '\n'.join([*lines[:-1], ' '.join(row)]), 'data row 3, ' + field
    if fmt == 'columns-ascii-all':
        for bad in ('oops', '32768', '-32769'):
            row = lines[-1].split()
            row[-1] = bad
            yield 'mask-' + bad, '\n'.join([*lines[:-1], ' '.join(row)]), 'data row 3, mask'


def ascii_recovery_worker(directory):
    """Isolated so a regression to legacy exit(1) cannot kill the test runner."""
    sys.path.insert(0, str(ROOT))
    import numpy as np
    from cyballs import cballs, CosmoComputationError, search_method_id
    check = unittest.TestCase()
    checked = 0
    for engine in ASCII_ENGINES:
        if search_method_id(engine) < 0:
            continue
        for fmt in ('columns-ascii', 'columns-ascii-all'):
            good = directory / (fmt + '-good.txt')
            good.write_text(ascii_catalog(fmt).replace('\n# 3', '\n#3'))
            bad = directory / (fmt + '-bad.txt')
            model = cballs()
            params = dict(searchMethod=engine, infileformat=fmt,
                          rootDir=str(directory / (engine + '-' + fmt)),
                          numberThreads=1, sizeHistN=4, mChebyshev=3,
                          rangeN=2., rminHist=.01, lengthBox=4., useLogHist=False,
                          verbose=0, verbose_log=0,
                          options='only-2pcf,no-one-ball,no-two-balls,no-smooth-pivot,weights-norm')
            try:
                cases = list(ascii_invalid_catalogs(fmt))
                # A constant-field override must not hide invalid input.
                cases.append(('constant-nan', ascii_catalog(fmt).replace('0 0 1 3', '0 0 1 nan'), 'kappa'))
                for label, content, diagnostic in cases:
                    bad.write_text(content)
                    options = params['options'] + (',kappa-constant-one' if label == 'constant-nan' else '')
                    model.set(params | dict(infile=str(bad), options=options))
                    with check.assertRaises(CosmoComputationError, msg=f'{engine}/{fmt}/{label}') as caught:
                        model.Run()
                    check.assertIn(diagnostic, str(caught.exception))
                    check.assertIn(str(bad), str(caught.exception))
                    check.assertFalse(model.getBodytableAllocated())
                    check.assertFalse(model._runtime_bodytable_address())
                    # Reuse the exact same object after each rejected catalog.
                    model.set(params | dict(infile=str(good)))
                    model.Run()
                    np.testing.assert_allclose(model.getHistXi2pcf(), [0., 0., 1./3., 0.], rtol=1e-13, atol=1e-15)
                    model.struct_cleanup()
                    check.assertFalse(model.getBodytableAllocated())
                    checked += 1
            finally:
                model.struct_cleanup()
    check.assertGreater(checked, 0)
    print(f'ASCII_RECOVERY_OK {checked}', flush=True)


class IOStabilization(unittest.TestCase):
    def setUp(self):
        directory = tempfile.TemporaryDirectory(prefix='cballs-io-regression-')
        self.addCleanup(directory.cleanup)
        self.root = Path(directory.name)
        self.runs = 0

    def params(self, **changes):
        self.runs += 1
        result = dict(searchMethod='octree-2balls-omp', nbody=8,
                      testmodel='simple-cubic', sizeHistN=4, mChebyshev=3,
                      numberThreads=1, verbose=2, verbose_log=0,
                      rootDir=str(self.root / str(self.runs)),
                      rangeN=2, rminHist=.01, lengthBox=10,
                      outfile='catalog', outfileformat='columns-ascii-all', options=OPTIONS)
        result.update(changes)
        return result

    def run_case(self, failure=None, **changes):
        params = self.params(**changes)
        p = subprocess.run([str(BINARY), *(f'{k}={v}' for k, v in params.items())],
                           cwd=self.root, text=True, capture_output=True, timeout=30)
        output = Path(params['rootDir']) / (params['outfile'] + '.txt')
        log = p.stdout + p.stderr
        if failure is not None:
            self.assertGreater(p.returncode, 0, log)  # failure, not a crash/signal
            self.assertIn(failure, log)
            self.assertFalse(output.exists(), log)
        else:
            self.assertEqual(p.returncode, 0, log)
        return output, log

    def assert_rows_close(self, actual, expected, tol=2e-10):
        self.assertEqual(len(actual), len(expected))
        for row, reference in zip(actual, expected):
            self.assertGreaterEqual(len(row), len(reference))
            self.assertTrue(all(math.isfinite(x) for x in row))
            for x, y in zip(row, reference):
                self.assertAlmostEqual(x, y, delta=tol * max(1.0, abs(y)))

    def test_native_ascii_malformed_input_returns_failure(self):
        for fmt in ('columns-ascii', 'columns-ascii-all'):
            path = self.root / (fmt + '.txt')
            for label, content, diagnostic in ascii_invalid_catalogs(fmt):
                with self.subTest(format=fmt, case=label):
                    path.write_text(content)
                    _, log = self.run_case(failure=diagnostic, infile=path, infileformat=fmt)
                    self.assertIn(str(path), log)

    def test_native_ascii_valid_signed_values_and_long_comment(self):
        for fmt in ('columns-ascii', 'columns-ascii-all'):
            with self.subTest(format=fmt):
                path = self.root / (fmt + '.txt')
                path.write_text(ascii_catalog(fmt).replace('# native ASCII regression', '# ' + 'comment ' * 100).rstrip())
                out, _ = self.run_case(infile=path, infileformat=fmt, lengthBox=4)
                header, rows = read_text_catalog(out)
                self.assertEqual(header, [3, 3, 4, 4, 4])
                self.assert_rows_close(rows, [(*map(float, row.split()), 1, 1) for row in ASCII_ROWS])

    def test_native_ascii_compact_header_compatibility(self):
        catalogs = {fmt: ascii_catalog(fmt) for fmt in
                    ('columns-ascii', 'columns-ascii-all')}
        catalogs['columns-ascii-pos'] = '# positions\n# 3 3 4 4 4\n1 0 0\n0 1 0\n0 0 1\n'
        catalogs['columns-ascii-2d-to-3d'] = '# angles\n# 3 2 4 4\n1 0 -1\n1.2 1 2\n.8 2 3\n'
        for fmt, good in catalogs.items():
            path = self.root / (fmt + '.txt')
            path.write_text(good)
            out, _ = self.run_case(infile=path, infileformat=fmt)
            reference = read_text_catalog(out)
            for marker in ('#3', '#\t3', ' \t#3'):
                with self.subTest(format=fmt, marker=marker):
                    path.write_text(good.replace('\n# 3', '\n' + marker))
                    out, _ = self.run_case(infile=path, infileformat=fmt)
                    self.assertEqual(read_text_catalog(out), reference)
            for count in ('0', '-3', 'nan', '3oops', '9' * 300):
                with self.subTest(format=fmt, count=count[:16]):
                    path.write_text(good.replace('\n# 3', '\n#' + count))
                    self.run_case(failure='header nbody', infile=path, infileformat=fmt)
            for content, error in [('# comment\n#', 'unexpected end of file'),
                                   (good.replace('\n# 3', '\n!3'), 'expected #'),
                                   (good.replace('\n# 3', '\n#3').rsplit(' ', 1)[0], 'unexpected end of file')]:
                with self.subTest(format=fmt, error=error):
                    path.write_text(content)
                    self.run_case(failure=error, infile=path, infileformat=fmt)

    def test_native_ascii_header_inspection(self):
        for fmt in ('columns-ascii', 'columns-ascii-all'):
            path = self.root / (fmt + '.txt')
            path.write_text(ascii_catalog(fmt))
            self.run_case(infile=path, infileformat=fmt, options='header-info,stop', outfile='')
            path.write_text('')
            self.run_case(failure='unexpected end of file', infile=path, infileformat=fmt,
                          options='header-info,stop', outfile='')

    def test_native_ascii_cython_survival_cleanup_and_recovery(self):
        if not WITH_CYTHON:
            self.skipTest('enable with --cython')
        completed = subprocess.run([sys.executable, str(Path(__file__).resolve()),
                                    '--ascii-recovery-worker', str(self.root)],
                                   cwd=self.root, text=True, capture_output=True, timeout=120)
        log = completed.stdout + completed.stderr
        self.assertEqual(completed.returncode, 0, log)
        self.assertIn('ASCII_RECOVERY_OK', log)

    def test_recorded_single_box_and_default_field(self):
        out, log = self.run_case(infile=FIXTURES / 'single.snap', infileformat='gadget', lengthBox=23)
        header, rows = read_text_catalog(out)
        self.assertEqual(header, [8, 3, 10, 10, 10])
        self.assert_rows_close(rows, [(*p, 1, 1) for p in POSITIONS])
        self.assertIn('scalar field = constant 1', log)

    def test_recorded_multi_wraps_z_from_z(self):
        out, _ = self.run_case(infile=FIXTURES / 'multi.snap', infileformat='gadget')
        header, rows = read_text_catalog(out)
        self.assertEqual(header, [8, 3, 10, 10, 10])
        self.assert_rows_close(rows, [(*p, 1, 1) for p in [(1, 2, 1), *POSITIONS[1:]]])

    def test_gadget_wrapping_boundaries_and_one_particle(self):
        for points in [[(-1, 10, 31), (30, -22, -13), (1000001, -1000002, 11)], [(0, 0, 0)]]:
            with self.subTest(points=points):
                path = self.root / 'boundary.snap'
                snapshot(path, points)
                out, _ = self.run_case(infile=path, infileformat='gadget')
                header, rows = read_text_catalog(out)
                self.assertEqual(header[2:], [10, 10, 10])
                self.assert_rows_close(rows, [(*(x % 10 for x in p), 1, 1) for p in points])

    def test_explicit_gadget_scalar_sources(self):
        path = self.root / 'field.snap'
        points = [(1.125, 2.125, 3), (3.25, 4.125, 5)]
        snapshot(path, points)
        for option, scalar in [('kappa-constant', 2), ('kappa-constant-one', 1),
                               ('kappa-constant,kappa-constant-one', 1), ('gadget-kappa-synthetic', None)]:
            with self.subTest(option=option):
                out, log = self.run_case(infile=path, infileformat='gadget', options=OPTIONS + ',' + option)
                values = [scalar if scalar is not None else
                          1 + math.cos(40 * math.pi * x / 10) * math.sin(40 * math.pi * y / 10)
                          for x, y, _ in points]
                self.assert_rows_close(read_text_catalog(out)[1], [(*p, value, 1) for p, value in zip(points, values)])
                self.assertIn('scalar field =', log)
        for option in ['kappa-constant', 'kappa-constant-one']:
            self.run_case(failure='conflicts', infile=path, infileformat='gadget',
                          options=OPTIONS + ',gadget-kappa-synthetic,' + option)

    def test_invalid_gadget_boxes_and_positions(self):
        path = self.root / 'invalid.snap'
        for box in [0, -10, math.nan, math.inf]:
            with self.subTest(box=box):
                snapshot(path, POSITIONS, box=box)
                self.run_case(failure='BoxSize', infile=path, infileformat='gadget')
        for value in [math.nan, math.inf, -math.inf]:
            with self.subTest(coordinate=value):
                snapshot(path, [(1, 2, value), *POSITIONS[1:]])
                self.run_case(failure='non-finite coordinate', infile=path, infileformat='gadget')
        path.write_bytes((FIXTURES / 'single.snap').read_bytes()[:-8])
        self.run_case(failure='error reading binary', infile=path, infileformat='gadget')

    def test_inconsistent_multi_box(self):
        prefix = self.root / 'multi.snap'
        for i in range(2):
            data = bytearray((FIXTURES / f'multi.snap.{i}').read_bytes())
            if i == 1:
                struct.pack_into('<d', data, 4 + 128, 20.0)
            Path(str(prefix) + f'.{i}').write_bytes(data)
        self.run_case(failure='inconsistent BoxSize', infile=prefix, infileformat='gadget')

    def test_recorded_radec_and_optional_weight(self):
        for columns, weighted in [('1,2,3,4', True), ('1,2,3,99', False)]:
            out, _ = self.run_case(infile=FIXTURES / 'radec.fits', infileformat='fits-radec-field',
                                  columns=columns, options=OPTIONS + ',no-arfken' + (',with-weight' if weighted else ''))
            expected = []
            for i in range(8):
                ra, dec = .05 + .75 * i / 7, .3 + .5 * i / 7
                expected.append((math.cos(dec) * math.cos(ra), math.cos(dec) * math.sin(ra), math.sin(dec), 1, 1))
            self.assert_rows_close(read_text_catalog(out)[1], expected)

    def test_bad_columns_fail_in_each_table_reader(self):
        # Reuse recorded columns as XYZ/scalar and scalar RADecR catalogs.
        for fmt, columns in [('fits-radec-field', [1, 2, 3, 4]),
                             ('fits', [2, 3, 4, 1, 4, 1]), ('fits-radecr-field', [2, 3, 4, 1, 4])]:
            for index in range(len(columns)):
                with self.subTest(format=fmt, column_index=index):
                    bad = columns.copy()
                    bad[index] = 99
                    self.run_case(failure='column 99', infile=FIXTURES / 'radec.fits',
                                  infileformat=fmt, columns=','.join(map(str, bad)),
                                  options=OPTIONS + ',no-arfken,with-weight')

    def test_valid_xyz_and_radecr(self):
        for fmt, columns in [('fits', '2,3,4,1,4,1'), ('fits-radecr-field', '2,3,4,1,4')]:
            out, _ = self.run_case(infile=FIXTURES / 'radec.fits', infileformat=fmt,
                                  columns=columns, options=OPTIONS + ',no-arfken,with-weight')
            expected = []
            for i in range(8):
                a, b = .05 + .75 * i / 7, .3 + .5 * i / 7
                pos = ((a, b, 1) if fmt == 'fits' else
                       (math.cos(math.radians(b)) * math.cos(math.radians(a)),
                        math.cos(math.radians(b)) * math.sin(math.radians(a)), math.sin(math.radians(b))))
                expected.append((*pos, 1, 1))
            self.assert_rows_close(read_text_catalog(out)[1], expected)

    def test_fits_bad_hdu_truncation_and_missing_file(self):
        for fmt in ['fits', 'fits-radec-field', 'fits-radecr-field', 'fits-healpix']:
            with self.subTest(format=fmt):
                self.run_case(failure='status=', infile=FIXTURES / 'radec.fits',
                              infileformat=fmt, options=OPTIONS + ',fits-type-file')
                self.run_case(failure='status=', infile=self.root / 'missing.fits', infileformat=fmt)
        path = self.root / 'truncated.fits'
        path.write_bytes((FIXTURES / 'radec.fits').read_bytes()[:5760])
        self.run_case(failure='status=', infile=path, infileformat='fits-radec-field')

    def test_header_inspection_succeeds(self):
        for fmt in ['fits', 'fits-radec-field', 'fits-radecr-field', 'fits-healpix']:
            self.run_case(infile=FIXTURES / 'radec.fits', infileformat=fmt,
                          options='header-info,stop-fits', outfile='')

    def test_fits_export_failure_and_valid_roundtrip(self):
        for fmt in ['fits', 'numpy-healpix']:
            with self.subTest(format=fmt):
                self.run_case(failure='fits_create_file', outfile='absent-parent/catalog',
                              outfileformat=fmt, options=OPTIONS + ',kappa')
                self.run_case(failure='requires options=kappa', outfileformat=fmt)
                out, _ = self.run_case(infile=FIXTURES / 'single.snap', infileformat='gadget',
                                      outfileformat=fmt, options=OPTIONS + ',kappa')
                self.assert_rows_close(read_fits_export(out), [(*p, 1) for p in POSITIONS])
                back, _ = self.run_case(infile=out, infileformat='fits', columns='1,2,3,4')
                self.assert_rows_close(read_text_catalog(back)[1], [(*p, 1, 1) for p in POSITIONS])

    def test_cython_failure_cleanup_and_recovery(self):
        if not WITH_CYTHON:
            self.skipTest('enable with --cython')
        sys.path.insert(0, str(ROOT))
        from cyballs import CosmoComputationError, cballs
        errors = [dict(infile=FIXTURES / 'radec.fits', infileformat='fits-radec-field',
                       columns='1,2,3,99', options=OPTIONS + ',with-weight'),
                  dict(outfile='absent-parent/catalog', outfileformat='fits', options=OPTIONS + ',kappa')]
        for changes in errors:
            for _ in range(3):
                balls = cballs()
                self.addCleanup(balls.struct_cleanup)
                params = self.params(**changes)
                balls.set({k: str(v) if isinstance(v, Path) else v for k, v in params.items()})
                with self.assertRaises(CosmoComputationError):
                    balls.Run(level=['MainLoop'])
                self.assertFalse(balls.getBodytableAllocated())
                self.assertFalse(balls._runtime_bodytable_address())
                good = self.params(infile=str(FIXTURES / 'single.snap'), infileformat='gadget',
                                   options=OPTIONS + ',gadget-kappa-synthetic')
                balls.set(good)
                balls.Run(level=['MainLoop'])
                self.assertEqual(balls.getNBody(), 8)
                self.assertTrue(balls.getBodytableAllocated())
                balls.struct_cleanup()
                self.assertFalse(balls.getBodytableAllocated())


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cballs', type=Path, default=BINARY)
    parser.add_argument('--cython', action='store_true')
    parser.add_argument('--ascii-recovery-worker', type=Path, help=argparse.SUPPRESS)
    args, remaining = parser.parse_known_args()
    if args.ascii_recovery_worker is not None:
        ascii_recovery_worker(args.ascii_recovery_worker)
        raise SystemExit(0)
    BINARY = args.cballs.resolve()
    WITH_CYTHON = args.cython
    unittest.main(argv=[sys.argv[0], *remaining], verbosity=2)
