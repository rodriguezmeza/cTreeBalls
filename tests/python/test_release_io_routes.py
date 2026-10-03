#!/usr/bin/env python3
"""Retained successful round trips through active core/IOLIB/Gadget/CFITSIO routes."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import struct
import numpy as np
import healpy as hp
from test_io_stabilization import FIXTURES, POSITIONS, read_text_catalog, read_fits_export


def run(binary, root):
    root.mkdir(parents=True, exist_ok=False)
    records = []
    def case(label, fmt, path, *, output_format='columns-ascii-all', options='', **extra):
        output = root/label
        params = dict(search='octree-2balls-omp', infile=str(path), infileformat=fmt,
                      rootDir=str(output), outfile='catalog', outfileformat=output_format,
                      numberThreads=1, verbose=0, verbose_log=0, sizeHistN=4,
                      mChebyshev=2, rangeN=2, rminHist=.01, nbody=12, iCatalogs='1',
                      options='stop,no-smooth-pivot,no-out-Hist'+(','+options if options else ''))
        params.update(extra)
        command = [str(binary), *[f'{k}={v}' for k,v in params.items()]]
        process = subprocess.run(command, capture_output=True, text=True, timeout=60)
        output.mkdir(exist_ok=True)
        (output/'run.log').write_text(process.stdout+process.stderr)
        assert process.returncode == 0, process.stdout+process.stderr
        product = output/'catalog.txt'
        assert product.stat().st_size > 0
        metadata = json.loads((output/'run-metadata.json').read_text())
        records.append(dict(label=label, input_format=fmt, output_format=output_format,
                            product=str(product.relative_to(root)), sha256=hashlib.sha256(product.read_bytes()).hexdigest(),
                            provenance=metadata))
        return product
    base = case('gadget', 'gadget', FIXTURES/'single.snap')
    expected = np.array(read_text_catalog(base)[1])
    for fmt in ('columns-ascii', 'columns-ascii-all', 'binary', 'binary-all',
                'columns-ascii-pos', 'fits', 'numpy-healpix'):
        product = case('write-'+fmt, 'columns-ascii-all', base, output_format=fmt,
                       options='kappa' if fmt in ('fits','numpy-healpix') else '')
        if fmt in ('fits','numpy-healpix'):
            np.testing.assert_allclose(read_fits_export(product), expected[:, :4], rtol=2e-10, atol=2e-10)
        if fmt == 'columns-ascii-pos':
            np.testing.assert_allclose(np.loadtxt(product), expected[:, :3], rtol=2e-10, atol=2e-10)
            # This writer is headerless; the position reader requires a native header.
            headered = root/'positions-with-header.txt'
            headered.write_text('# position fixture\n# 8 3 10 10 10\n'+product.read_text())
            product = headered
        back = case('read-'+fmt, 'fits' if fmt=='numpy-healpix' else fmt, product,
                    columns='1,2,3,4', options='kappa-constant-one' if fmt=='columns-ascii-pos' else '')
        actual = np.array(read_text_catalog(back)[1])
        reference = expected[:, :5].copy()
        if fmt == 'columns-ascii-pos': reference[:, 3] = 2.0  # documented loader default
        np.testing.assert_allclose(actual[:, :5], reference, rtol=2e-10, atol=2e-10)
    points = np.array(POSITIONS, dtype=float)
    multi = root/'multi.txt'
    np.savetxt(multi, np.column_stack((points, np.ones(len(points)), np.ones(len(points)))))
    out = case('multi-columns','multi-columns-ascii',multi, columns='1,2,3,4,5',options='pos-and-convergence-weight')
    np.testing.assert_allclose(np.array(read_text_catalog(out)[1])[:, :5], expected[:, :5], rtol=2e-10, atol=2e-10)
    angles = np.column_stack((np.linspace(.3,.8,8), np.linspace(.1,1.4,8), np.ones(8)))
    angular = root/'angles.txt'
    np.savetxt(angular, angles, header='angular fixture\n8 2 2 2')
    out = case('2d-to-3d','columns-ascii-2d-to-3d',angular)
    coordinates = np.column_stack((np.sin(angles[:,0])*np.cos(angles[:,1]),
                                  np.sin(angles[:,0])*np.sin(angles[:,1]), np.cos(angles[:,0])))
    np.testing.assert_allclose(np.array(read_text_catalog(out)[1])[:, :3], coordinates, rtol=2e-10, atol=2e-10)
    ra_dec = root/'ra-dec.txt'
    np.savetxt(ra_dec, angles[:, [1,0,2]])
    out = case('ra-dec-ascii','ra-dec-ascii',ra_dec,columns='1,2,3')
    np.testing.assert_allclose(np.array(read_text_catalog(out)[1])[:, :3], coordinates, rtol=2e-10, atol=2e-10)
    for fmt, columns in (('fits-radec-field','1,2,3,4'), ('fits-radecr-field','2,3,4,1,4')):
        out = case(fmt,fmt,FIXTURES/'radec.fits',columns=columns, options='no-arfken,with-weight')
        assert np.array(read_text_catalog(out)[1]).shape[0] == 8
    values = np.arange(12, dtype=np.float64)+.5
    raw = root/'ring-native-f64.bin'; values.tofile(raw)
    fits = root/'ring.fits'; hp.write_map(fits, values, dtype=np.float64, nest=False)
    expected_xyz = np.array(hp.pix2vec(1, np.arange(12))).T
    for fmt, path in (('numpy-healpix',raw), ('fits-healpix',fits)):
        out = case('map-'+fmt,fmt,path)
        rows = np.array(read_text_catalog(out)[1])
        np.testing.assert_allclose(rows[:, :3], expected_xyz, rtol=2e-10, atol=2e-10)
        np.testing.assert_allclose(rows[:, 3], values, rtol=2e-10, atol=2e-10)
    takahashi = root/'takahashi-native.bin'
    assert struct.calcsize('@l') == 8, 'Takahashi fixture requires the active 64-bit-long ABI'
    takahashi.write_bytes(struct.pack('@iill', 4, 1, 12, 0) + values.astype(np.float32).tobytes()
                         + (struct.pack('@l',0)+np.zeros(12,dtype=np.float32).tobytes())*3)
    out = case('takahashi','takahashi',takahashi)
    rows = np.array(read_text_catalog(out)[1])
    np.testing.assert_allclose(rows[:, :3], expected_xyz, rtol=2e-10, atol=2e-10)
    np.testing.assert_array_equal(rows[:, 3], values)
    selected = np.arange(12)%3 != 0
    embedded = root/'embedded-mask.fits'
    hp.write_map(embedded, np.where(selected,values,0), dtype=np.float64)
    for fmt in ('fits-healpix','numpy-healpix'):
        # The legacy numpy-healpix mask-inside branch explicitly consumes FITS,
        # while its ordinary map branch consumes raw native-endian doubles.
        out = case('embedded-'+fmt,fmt,embedded,options='mask-inside')
        rows = np.array(read_text_catalog(out)[1])
        np.testing.assert_allclose(rows[:, :3], expected_xyz[selected], rtol=2e-10, atol=2e-10)
        np.testing.assert_array_equal(rows[:, 3], values[selected])
    mask = root/'companion-mask.fits'
    hp.write_map(mask, selected.astype(float), dtype=np.float64)
    for fmt, path in (('fits-healpix',fits),('numpy-healpix',raw)):
        out = case('companion-'+fmt,fmt+','+fmt,str(path)+','+str(mask),options='read-mask',iCatalogs='1,1')
        rows = np.array(read_text_catalog(out)[1])
        np.testing.assert_array_equal(rows[:, -1], selected)
    (root/'routes.json').write_text(json.dumps(records,indent=2,sort_keys=True)+'\n')
    print(f'PASS: {len(records)} retained I/O route cases')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cballs', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    run(args.cballs.resolve(), args.output.resolve())
