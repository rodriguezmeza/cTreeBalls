#!/usr/bin/env python3
"""One retained, independently checked active-engine case (also an MPI worker)."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / 'tests/make_tests')]


def run(engine, output, threads):
    import numpy as np
    from cyballs import cballs
    comm = None
    if engine.endswith('-mpi'):
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
    rank = comm.Get_rank() if comm else 0
    if rank == 0:
        output.mkdir(parents=True, exist_ok=False)
    if comm:
        comm.Barrier()
    model = cballs()
    p = dict(searchMethod=engine, rootDir=str(output), numberThreads=threads,
             verbose=0, verbose_log=0, iCatalogs='1', usePeriodic=False,
             useLogHist=False, sizeHistN=4, mChebyshev=2, nsmooth=2,
             sizeHistPhi=8, rangeN=1.5, rminHist=.02, theta=1.,
             options='no-smooth-pivot')
    kind = None
    if engine.startswith('lya-'):
        import test_lya_forest_mpi as forest
        same = 'same-los' in engine
        base = engine.rsplit('-', 1)[0].replace('lya-los-tree-', 'lya-')
        kind = 6 if same else forest.METHODS.index(base)
        data = forest.three.POINTS if kind < 3 else forest.radial.WIDE_ANGLE
        p.update(forest.params(kind, '', output, threads))
        p.pop('infile'); p.pop('infileformat')
        model.set(p)
        model.set_forest_catalog(data[:, :3], data[:, 3], data[:, 4], data[:, 5].astype(np.int64))
        fixture = dict(positions=data[:, :3], delta=data[:, 3], weights=data[:, 4], forest_ids=data[:, 5])
    elif 'shear' in engine:
        import test_shear_sphere_octree_omp as shear
        pos, gamma, weights = shear.fixture()
        p.update(rangeN=shear.RMAX, rminHist=shear.RMIN, sizeHistN=shear.BINS,
                 sizeHistPhi=shear.PHI_BINS, mChebyshev=shear.NMAX,
                 options='no-smooth-pivot,no-one-ball')
        model.set(p)
        model.set_catalog(pos, weights=weights, gamma1=gamma.real, gamma2=gamma.imag)
        fixture = dict(positions=pos, gamma=gamma, weights=weights)
    elif '3pcf-3d' in engine:
        import test_octree_3pcf_3d_omp as physical
        data = np.array(physical.CATALOG)
        p.update(rangeN=physical.RMAX, rminHist=physical.RMIN, sizeHistN=physical.NBINS,
                 mChebyshev=physical.LMAX, options='compute-2pcf-3d,compute-3pcf-3d')
        model.set(p)
        model.set_catalog(data[:, :3], kappa=data[:, 3], weights=data[:, 4])
        fixture = dict(positions=data[:, :3], kappa=data[:, 3], weights=data[:, 4])
    elif 'box' in engine:
        import test_neighbor_boxes_periodic as boxes
        pos = np.random.default_rng(83519).uniform(0., boxes.LBOX, (96, 3))
        p.update(lengthBox=boxes.LBOX, usePeriodic=True, rangeN=boxes.RANGE,
                 rminHist=0., sizeHistN=boxes.NBINS, options='compute-HistN,no-smooth-pivot')
        model.set(p)
        model.set_catalog(pos, kappa=np.ones(len(pos)))
        fixture = dict(positions=pos, kappa=np.ones(len(pos)))
    else:
        import test_two_ball_edge_corrections as scalar
        data = scalar.catalog()
        core = engine == 'octree-sincos-omp'
        options = 'KKKCorrelation,no-smooth-pivot,no-normalize-HistZeta,weights-norm'
        options += ',no-one-ball' if core else ',edge-corrections,no-one-ball,no-two-balls'
        p.update(rangeN=scalar.RMAX, rminHist=scalar.RMIN, sizeHistN=scalar.BINS,
                 mChebyshev=scalar.MMAX, options=options)
        model.set(p)
        model.set_catalog(data[0], kappa=data[1], weights=data[2], mask=data[3])
        fixture = dict(positions=data[0], kappa=data[1], weights=data[2], mask=data[3])
    try:
        model.Run(level=['MainLoop'])
        metadata = model.getRunMetadata()
        if rank != 0:
            return
        arrays = {}
        if engine.startswith('lya-'):
            if same:
                forest.radial.assert_histogram_close(
                    forest.radial.read_2pcf(output/'histXi2pcf_lya1d_same_los.txt'),
                    forest.radial.oracle_same_los_2pcf(), 'same LOS')
            else:
                forest.check_oracle(kind, output)
        elif 'shear' in engine:
            arrays = dict(xi_plus=model.getShearXiPlus(), xi_minus=model.getShearXiMinus(),
                          xi_weight=model.getShearXiWeight(), upsilon=model.getShearUpsilonXMultipoles(),
                          window=model.getShearWindowMultipoles(), multipoles=model.getShearGammaXMultipoles())
            expected = shear.oracle(pos, gamma, weights)
            for key in arrays:
                np.testing.assert_allclose(arrays[key], expected[key], rtol=5e-10, atol=5e-12, err_msg=key)
        elif '3pcf-3d' in engine:
            physical.check_oracle(output)
        elif 'box' in engine:
            arrays['pair_counts'] = model.getHistNN().copy()
            expected = boxes.expected_pair_counts(pos)
            if engine == 'kdtree-box-omp':
                expected = expected / 2
            np.testing.assert_array_equal(arrays['pair_counts'], expected)
        elif engine != 'octree-sincos-omp':
            actual = scalar.load_results(output)
            signal, window = scalar.brute_force(data)
            scalar.assert_results(actual, (signal, window, scalar.edge_solution(signal, window)))
            arrays.update(zip(('signal', 'window', 'corrected'), actual))
        else:
            arrays['pair_counts'] = model.getHistNN().copy()
            arrays['xi'] = model.getHistCF().copy()
            assert np.all(np.isfinite(arrays['xi']))
        for path in sorted(output.glob('hist*.txt')):
            # Shear and scalar raw arrays above retain full precision; retain text products too.
            values = np.loadtxt(path)
            if values.size:
                arrays[path.name] = values
        if not arrays:
            raise AssertionError('engine produced no numerical products')
        np.savez_compressed(output/'fixture.npz', **fixture)
        np.savez_compressed(output/'result.npz', **arrays)
        metadata['fixture'] = {'file': 'fixture.npz', 'sha256': hashlib.sha256((output/'fixture.npz').read_bytes()).hexdigest()}
        metadata['validation'] = 'independent oracle plus cross-execution comparison' if engine != 'octree-sincos-omp' else 'finite products and cross-execution comparison; core numerical scripts run separately'
        (output/'result.json').write_text(json.dumps(metadata, indent=2, sort_keys=True)+'\n')
    finally:
        model.struct_cleanup()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--engine', required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--threads', type=int, required=True)
    args = parser.parse_args()
    run(args.engine, args.output.resolve(), args.threads)


if __name__ == '__main__':
    main()
