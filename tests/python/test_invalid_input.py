#!/usr/bin/env python3
"""Isolated malformed-catalog and two-rank failure/recovery contracts.

Use the active double/3D/LONGINT profile. All checks remain enabled under -O.
The parent owns timeouts: legacy exit(1) and MPI deadlocks fail this test.
"""
import argparse
import json
import os
from pathlib import Path
import shlex
import signal
import struct
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
ENGINES = ('kdtree-2balls', 'balltree-2balls', 'octree-2balls')
CHECK = unittest.TestCase()
BASE_OPTIONS = 'only-2pcf,no-one-ball,no-two-balls,no-smooth-pivot,weights-norm'


def parameters(directory, engine, path, fmt, extra=None):
    result = dict(searchMethod=engine, infile=str(path), infileformat=fmt,
                  rootDir=str(directory/'output'), numberThreads=1, sizeHistN=4,
                  mChebyshev=3, rangeN=2., rminHist=.01, lengthBox=4.,
                  useLogHist=False, verbose=0, verbose_log=0, options=BASE_OPTIONS)
    if extra:
        result.update({k: v for k, v in extra.items() if k != 'options'})
        result['options'] += ',' + extra.get('options', '')
    return result


def native_ascii(fmt):
    if fmt == 'columns-ascii-pos':
        return '# test\n# 3 3 4 4 4\n1 0 0\n0 1 0\n0 0 1\n'
    if fmt == 'columns-ascii-2d-to-3d':
        return '# test\n# 3 2 4 4\n1 0 -1\n1.2 1 2\n.8 2 3\n'
    return '# test\n# 3 3 4 4 4\n1 0 0 -1\n0 1 0 2\n0 0 1 3\n'


def binary(fmt):
    raw = struct.pack('=qi3d9d3d', 3, 3, 4, 4, 4,
                      1, 0, 0, 0, 1, 0, 0, 0, 1, -1, 2, 3)
    return raw + (struct.pack('=3d3h', 1, 1, 1, 1, 1, 1) if fmt == 'binary-all' else b'')


def fits_bytes(path, data):
    from astropy.io import fits
    columns = [fits.Column(name=f'C{i+1}', format='D', array=data[:, i])
               for i in range(data.shape[1])]
    fits.BinTableHDU.from_columns(columns).writeto(path, overwrite=True)
    return path.read_bytes()


def catalogs(directory):
    """Each group has a finite control and corruptions of that exact layout."""
    import numpy as np
    groups = []
    for fmt in ('binary', 'binary-all'):
        good = binary(fmt)
        bad = [(f'truncated-{i}', good[:i], {}) for i in range(len(good))]
        # Box dimensions, all coordinate components, scalar and optional weight.
        offsets = [12 + 8*i for i in range(3+9+3+(3 if fmt.endswith('all') else 0))]
        for offset in offsets:
            for invalid in (float('nan'), float('inf'), -float('inf')):
                bad.append((f'nonfinite-{offset}-{invalid}', good[:offset] + struct.pack('=d', invalid) + good[offset+8:], {}))
        bad += [('zero-count', struct.pack('=q', 0) + good[8:], {}),
                ('wrong-dimension', good[:8] + struct.pack('=i', 2) + good[12:], {}),
                ('oversize-count', struct.pack('=q', 2**63-1) + good[8:], {})]
        groups.append((fmt, {}, good, bad))
    for fmt in ('columns-ascii-pos', 'columns-ascii-2d-to-3d'):
        good = native_ascii(fmt).encode()
        tokens = good.decode().splitlines()[1:]
        tokens = '\n'.join(tokens).split()
        bad = [(f'truncated-{i}', ('# test\n'+' '.join(tokens[:i])).encode(), {}) for i in range(len(tokens))]
        for i in range(3, len(tokens)):
            for invalid in ('nan', 'inf', '-inf', '1e309', '2oops'):
                altered = tokens.copy(); altered[i] = invalid
                bad.append((f'invalid-{i}-{invalid}', ('# test\n'+' '.join(altered)).encode(), {}))
        groups.append((fmt, {}, good, bad))
    xyz = np.column_stack((np.eye(3), [-1., 2., 3.], [1., 1., 1.]))
    spherical = np.array([[0., 20., 1., 1., 1.], [50., 40., 1., 2., 1.], [100., 60., 1., 3., 1.]])
    layouts = [
        ('multi-columns-ascii', xyz, dict(columns='1,2,3,4,5', options='pos-and-convergence-weight')),
        ('multi-columns-ascii', xyz[:, :3], dict(columns='1,2,3', options='only-pos')),
        ('multi-columns-ascii', xyz, dict(columns='1,2,3,4,5', options='pos-and-shear')),
        ('ra-dec-ascii', spherical[:, [0,1,3]], dict(columns='1,2,3', options='in-degrees')),
        ('fits', xyz, dict(columns='1,2,3,4,5', options='with-weight')),
        ('fits-radec-field', spherical[:, [0,1,3,4]], dict(columns='1,2,3,4', options='with-weight,in-degrees')),
        ('fits-radecr-field', spherical, dict(columns='1,2,3,4,5', options='with-weight,in-degrees')),
    ]
    for fmt, data, extra in layouts:
        def encode(a):
            if fmt.startswith('fits'):
                return fits_bytes(directory/'fixture.fits', a)
            # A final newline is deliberately absent.
            return ('# comments\n% more comments\n\n'+'\n'.join(' '.join(map(str, row)) for row in a)).encode()
        good = encode(data)
        bad = []
        for column in range(data.shape[1]):
            for invalid in (np.nan, np.inf, -np.inf):
                altered = data.copy(); altered[-1, column] = invalid
                bad.append((f'nonfinite-{column}-{invalid}', encode(altered), {}))
        altered = data.copy(); altered[-1, -2 if data.shape[1] == 5 else -1] = np.nan
        bad.append(('constant-override', encode(altered), dict(options=extra['options']+',kappa-constant,kappa-constant-one')))
        if not fmt.startswith('fits'):
            bad.extend([('missing-column', good.rsplit(b' ', 1)[0], {}),
                        ('bad-token', good+b'junk', {}), ('embedded-nul', good+b'\0', {}),
                        ('column-range', good, dict(columns='1,2,99,4,5'))])
        groups.append((fmt, extra, good, bad))
    # LOS IDs cannot be null, nonfinite, fractional or vector-valued, and
    # integer identifiers beyond 2**53 must survive without float conversion.
    from astropy.io import fits
    def encode_los(values,format='K',null=None):
        cols=[fits.Column(name=f'C{i+1}',format='D',array=xyz[:,i]) for i in range(5)]
        cols.append(fits.Column(name='LOS',format=format,array=values,null=null))
        path=directory/'los.fits';fits.BinTableHDU.from_columns(cols).writeto(path,overwrite=True)
        return path.read_bytes()
    good=encode_los(np.array([2**53+1,2**53+2,2**53+3],dtype=np.int64))
    bad=[('null-los',encode_los([1,2,-999],null=-999),{}),
         ('vector-los',encode_los([[1,2],[3,4],[5,6]],format='2K'),{}),
         ('fractional-los',encode_los([1.,2.,3.5],format='D'),{})]
    bad += [(f'nonfinite-los-{x}',encode_los([1.,2.,x],format='D'),{}) for x in (np.nan,np.inf,-np.inf)]
    path=directory/'los-scaled.fits';path.write_bytes(good)
    with fits.open(path,mode='update') as hdus:hdus[1].header['TSCAL6']=.5
    bad.append(('fractional-scaling',path.read_bytes(),{}))
    groups.append(('fits',dict(columns='1,2,3,4,5,6',options='with-weight'),good,bad))
    # FITS null pixels must stay nonfinite until validation; never become zero.
    import healpy as hp
    map_data = np.arange(12, dtype=float) + 1
    for fmt in ('fits-healpix', 'numpy-healpix'):
        def encode_map(a):
            if fmt == 'numpy-healpix': return a.astype('=f8').tobytes()
            hp.write_map(str(directory/'map.fits'), a, dtype=np.float64, overwrite=True)
            return (directory/'map.fits').read_bytes()
        good = encode_map(map_data)
        bad = []
        for invalid in (np.nan, np.inf, -np.inf):
            altered = map_data.copy(); altered[-1] = invalid
            bad.append((f'nonfinite-{invalid}', encode_map(altered), {}))
        if fmt == 'numpy-healpix': bad.append(('truncated', good[:-1], {}))
        groups.append((fmt, dict(nbody=12), good, bad))
    # Native Takahashi layout: every byte boundary, every map, and both
    # geometry/selection failure branches. Overrides must not hide NaNs.
    values=np.arange(12,dtype=np.float32)+.5
    good=struct.pack('@iill',4,1,12,0)+values.tobytes()+(struct.pack('@l',0)+np.zeros(12,dtype=np.float32).tobytes())*3
    bad=[(f'truncated-{i}',good[:i],{}) for i in range(len(good))]
    for field in range(4):
        offset=struct.calcsize('@iill')+field*(48+8)
        for invalid in (np.nan,np.inf,-np.inf):
            altered=good[:offset]+struct.pack('@f',invalid)+good[offset+4:]
            bad.append((f'map-{field}-{invalid}',altered,{}))
            bad.append((f'override-{field}-{invalid}',altered,dict(options='kappa-constant-one')))
    bad += [('bad-nside',good[:4]+struct.pack('@i',0)+good[8:],{}),
            ('bad-npix',good[:8]+struct.pack('@l',11)+good[16:],{}),
            ('empty-patch',good,dict(options='patch',thetaL=.001,thetaR=.002,phiL=.001,phiR=.002))]
    groups.append(('takahashi',{},good,bad))
    return groups


def reader_worker(directory):
    import numpy as np
    sys.path.insert(0, str(ROOT))
    from cyballs import cballs, CosmoComputationError, CosmoSevereError, search_method_id
    groups = catalogs(directory)
    count = 0
    for stem in ENGINES:
        engine = stem+'-omp'
        if search_method_id(engine) < 0: continue
        model = cballs()
        for group, (fmt, extra, good, bads) in enumerate(groups):
            path = directory/f'catalog-{group}.dat'
            path.write_bytes(good)
            params = parameters(directory, engine, path, fmt, extra)
            model.set(params); model.Run()
            expected = model.getHistXi2pcf().copy()
            CHECK.assertTrue(np.isfinite(expected).all(), (engine, fmt, expected))
            if fmt in ('binary', 'binary-all', 'multi-columns-ascii', 'fits') and 'pos-and-shear' not in params['options'] and 'only-pos' not in params['options']:
                np.testing.assert_allclose(expected, [0., 0., 1./3., 0.], rtol=1e-13, atol=1e-15)
            model.struct_cleanup()
            for label, malformed, override in bads:
                path.write_bytes(malformed)
                updated = dict(extra, **override)
                model.set(parameters(directory, engine, path, fmt, updated))
                with CHECK.assertRaises(CosmoComputationError, msg=(engine, fmt, label)) as caught:
                    model.Run()
                CHECK.assertIn(str(path), str(caught.exception), (engine, fmt, label))
                CHECK.assertFalse(model.getBodytableAllocated())
                CHECK.assertFalse(model._runtime_bodytable_address())
                path.write_bytes(good)
                model.set(params); model.Run()
                np.testing.assert_allclose(model.getHistXi2pcf(), expected, rtol=1e-13, atol=1e-15)
                model.struct_cleanup()
                count += 1
        model.clean_all()
    CHECK.assertGreater(count, 0)
    print('INPUT_CONTRACTS_OK', json.dumps(dict(rejected_and_recovered=count)), flush=True)


def mpi_worker(directory, engine, ownership):
    # Owned mode lets cballs initialize MPI before importing MPI itself.
    if ownership == 'borrowed':
        from mpi4py import MPI
    else:
        import mpi4py
        mpi4py.rc.initialize = False
    import numpy as np
    sys.path.insert(0, str(ROOT))
    from cyballs import cballs, CosmoComputationError, CosmoSevereError
    rank = int(os.environ.get('OMPI_COMM_WORLD_RANK', os.environ.get('PMI_RANK', '0')))
    path = directory/f'rank-{rank}.dat'
    model = cballs()
    count = 0
    for failing_rank in (0, 1):
        for fmt, kind in [('columns-ascii', 'truncated'), ('columns-ascii', 'nonfinite'),
                          ('binary', 'truncated'), ('binary-all', 'nonfinite'),
                          ('columns-ascii-pos', 'truncated'), ('columns-ascii-2d-to-3d', 'truncated'),
                          ('columns-ascii', 'missing'), ('columns-ascii', 'parameter')]:
            good = binary(fmt) if fmt.startswith('binary') else native_ascii(fmt).encode()
            # Mixed legacy/canonical header spellings must parse identically.
            if not fmt.startswith('binary') and rank == 0:
                good = good.replace(b'\n# 3', b'\n#3')
            if rank == failing_rank:
                if kind == 'parameter': bad = good
                elif kind == 'truncated': bad = good[:-5]
                elif fmt.startswith('binary'): bad = good[:108] + struct.pack('=d', float('nan')) + good[116:]
                else: bad = good.replace(b'0 0 1 3', b'0 0 1 nan')
                path.write_bytes(bad)
                if kind == 'missing': path.unlink()
            else: path.write_bytes(good)
            params = parameters(directory, engine, path, fmt)
            model.set(params | (dict(sizeHistN=0) if kind == 'parameter' and rank == failing_rank else {}))
            with CHECK.assertRaises(CosmoComputationError) as caught: model.Run()
            from mpi4py import MPI
            comm = MPI.COMM_WORLD
            CHECK.assertEqual(comm.size, 2)
            messages = comm.allgather(str(caught.exception))
            CHECK.assertEqual(messages[0], messages[1])
            if kind != 'parameter':
                CHECK.assertIn(str(directory/f'rank-{failing_rank}.dat'), messages[rank])
            CHECK.assertFalse(model.getBodytableAllocated())
            CHECK.assertFalse(model._runtime_bodytable_address())
            CHECK.assertFalse(MPI.Is_finalized())
            path.write_bytes(good); comm.Barrier()
            model.set(params); model.Run()
            if rank == 0:
                result = model.getHistXi2pcf().copy()
                CHECK.assertTrue(np.isfinite(result).all())
                if fmt in ('columns-ascii', 'binary', 'binary-all'):
                    np.testing.assert_allclose(result, [0., 0., 1./3., 0.], rtol=1e-13, atol=1e-15)
            else:
                # Only the publisher exposes reduced scientific products.
                with CHECK.assertRaisesRegex(CosmoSevereError, 'not computed'):
                    model.getHistXi2pcf()
            model.struct_cleanup(); comm.Barrier()
            count += 1
    model.clean_all()
    if rank == 0: print('MPI_INPUT_CONTRACTS_OK', engine, ownership, count, flush=True)


def run_child(argv, directory, marker, timeout):
    with (directory/'worker.log').open('w') as log:
        process = subprocess.Popen(argv, cwd=directory, stdout=log, stderr=subprocess.STDOUT,
                                   start_new_session=True)
        try: process.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            os.killpg(process.pid, signal.SIGKILL); process.wait()
            raise AssertionError(f'worker timed out: {argv}; log: {directory}/worker.log')
    output = (directory/'worker.log').read_text(errors='replace')
    CHECK.assertEqual(process.returncode, 0, output[-12000:])
    CHECK.assertIn(marker, output, output[-12000:])
    print('\n'.join(line for line in output.splitlines() if marker in line), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--worker', choices=('readers', 'mpi'))
    parser.add_argument('--directory', type=Path)
    parser.add_argument('--engine')
    parser.add_argument('--ownership', choices=('owned', 'borrowed'))
    parser.add_argument('--mpi-command', help='two-rank MPI launcher, including -n 2')
    parser.add_argument('--output', type=Path)
    parser.add_argument('--skip-readers', action='store_true', help='focused MPI rerun')
    args = parser.parse_args()
    if args.worker == 'readers': return reader_worker(args.directory)
    if args.worker == 'mpi': return mpi_worker(args.directory, args.engine, args.ownership)
    def run(directory):
        directory.mkdir(parents=True, exist_ok=True)
        readers = directory/'readers'; readers.mkdir(exist_ok=True)
        if not args.skip_readers:
            run_child([sys.executable, __file__, '--worker', 'readers', '--directory', str(readers)], readers, 'INPUT_CONTRACTS_OK', 180)
        if args.mpi_command:
            sys.path.insert(0, str(ROOT))
            from cyballs import search_method_id
            for stem in ENGINES:
                engine = stem+'-mpi'
                if search_method_id(engine) < 0: continue
                for ownership in ('owned', 'borrowed'):
                    target = directory/(engine+'-'+ownership); target.mkdir(exist_ok=True)
                    run_child(shlex.split(args.mpi_command)+[sys.executable, __file__, '--worker', 'mpi',
                        '--engine', engine, '--ownership', ownership, '--directory', str(target)],
                        target, 'MPI_INPUT_CONTRACTS_OK', 45)
    if args.output: run(args.output.resolve())
    else:
        with tempfile.TemporaryDirectory(prefix='cballs-input-contract-') as temp: run(Path(temp))


if __name__ == '__main__': main()
