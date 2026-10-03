#!/usr/bin/env python3
"""Retained MPI shear oracles, order/cross-catalog and root-publication contracts.

All checks survive Python -O. Each launcher invocation has a parent timeout.
MPI ranks hold identical catalogs; only rank zero publishes global results.
"""
import argparse
import json
import os
from pathlib import Path
import shlex
import sys
import traceback
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT), str(ROOT/'tests/python')]
CHECK = unittest.TestCase()


def worker(directory):
    import numpy as np
    from mpi4py import MPI
    from cyballs import cballs, search_method_id
    import test_shear_sphere_octree_omp as ref
    comm = MPI.COMM_WORLD
    rank = comm.rank
    records = []
    first = ref.fixture(45)
    second = ref.fixture(31)
    third = ref.fixture(23)
    masks = [np.ones(len(c[0]), dtype=np.uint8) for c in (first, second, third)]
    for m in masks: m[::4] = 0
    poles = ref.fixture(31)
    poles[0][:4] = [[0, 0, 1], [.1, 0, -np.sqrt(.99)], [1, 0, 0], [0, 1, 0]]
    scenarios = [
        dict(name='linear-oracle', catalogs=[first], oracle=True),
        dict(name='log-oracle', catalogs=[ref.fixture(60)], oracle=True, log=True, leaf=1),
        dict(name='poles-oracle', catalogs=[poles], oracle=True, leaf=32),
        dict(name='masked-cross-two', catalogs=[first, second], masks=masks[:2], log=True),
        dict(name='masked-cross-three', catalogs=[first, second, third], masks=masks, leaf=1),
        dict(name='repeated-role-121', catalogs=[first, second, third], masks=masks, roles='1,2,1'),
        dict(name='tiny-idle-ranks', catalogs=[ref.fixture(3)], leaf=1),
        dict(name='clustered-approximate', catalogs=[ref.octant_fixture(320)], theta=.05, approximate=True),
        dict(name='smooth-exact', catalogs=[ref.smooth_fixture()], smooth=True, oracle=True),
        dict(name='smooth-approximate', catalogs=[ref.smooth_fixture()], smooth=True, theta=.05, approximate=True),
    ]
    def configure(model, engine, case, order, threads):
        opt = 'no-out-Hist,' + ('smooth-pivot' if case.get('smooth') else 'no-smooth-pivot')
        if not case.get('approximate'): opt += ',no-one-ball,no-two-balls'
        if case.get('masks') is not None: opt += ',read-mask'
        if order != 'both': opt += ',only-'+order+'pcf'
        p = dict(searchMethod=engine, rootDir=str(directory/'native'), options=opt,
                 numberThreads=threads, iCatalogs=case.get('roles') or ','.join(str(i+1) for i in range(len(case['catalogs']))),
                 usePeriodic=False, useLogHist=case.get('log', False), rangeN=ref.RMAX,
                 rminHist=ref.RMIN, sizeHistN=ref.BINS, sizeHistPhi=ref.PHI_BINS,
                 mChebyshev=ref.NMAX, theta=case.get('theta', 0.), nsmooth=case.get('leaf', 8),
                 lengthBox=2.2, verbose=0, verbose_log=0)
        if case.get('smooth'): p['rsmooth']='40'
        model.set(p)
        for i, (pos, gamma, weight) in enumerate(case['catalogs']):
            kw = dict(catalog=i, weights=weight, gamma1=gamma.real, gamma2=gamma.imag)
            if case.get('masks') is not None: kw['mask']=case['masks'][i]
            model.set_catalog(pos, **kw)

    def extract(model, order):
        arrays={}
        if order!='3':
            arrays.update(xi_plus=model.getShearXiPlus(), xi_minus=model.getShearXiMinus(), xi_weight=model.getShearXiWeight())
        if order!='2':
            arrays.update(upsilon=model.getShearUpsilonXMultipoles(), window=model.getShearWindowMultipoles(), multipoles=model.getShearGammaXMultipoles())
        return {k:v.copy() for k,v in arrays.items()}

    def compare(actual, target):
        maximum=0.
        for key, value in actual.items():
            expected=target[key]
            np.testing.assert_array_equal(np.isfinite(value), np.isfinite(expected), err_msg=key)
            if key=='multipoles':
                # The retained oracle rejects ill-conditioned windows (cond>1e10).
                # Raw U/N sums still use the strict comparison above/below. A
                # forward-error assertion on an unstable solve is misleading;
                # require a small backward residual for each nonzero solution.
                for i in range(ref.BINS):
                    for j in range(ref.BINS):
                        indices=np.arange(2*ref.NMAX+1)
                        matrices=[v['window'][indices[:,None]-indices[None,:]+2*ref.NMAX,i,j]
                                  for v in (actual,target)]
                        if max(np.linalg.cond(m) for m in matrices) <= 1e8:
                            np.testing.assert_allclose(value[:,:,i,j],expected[:,:,i,j],rtol=5e-10,atol=5e-12)
                        else:
                            for v,m in zip((actual,target),matrices):
                                x=v[key][:,:,i,j].T; b=v['upsilon'][:,:,i,j].T
                                if not np.any(x):continue # Existing singular-window zero sentinel.
                                scale=np.linalg.norm(m)*np.linalg.norm(x)+np.linalg.norm(b)
                                CHECK.assertLessEqual(np.linalg.norm(m@x-b),1e-11*scale+1e-24)
                finite=np.isfinite(value)
                if finite.any():maximum=max(maximum,float(np.max(np.abs(value[finite]-expected[finite]))))
                continue
            np.testing.assert_allclose(value, expected, rtol=5e-10, atol=5e-12, equal_nan=True, err_msg=key)
            finite=np.isfinite(value)
            if finite.any(): maximum=max(maximum,float(np.max(np.abs(value[finite]-expected[finite]))))
        return maximum

    for tree in ('octree','kdtree','balltree'):
        engine=tree+'-shear-sphere-2balls-mpi'
        if search_method_id(engine)<0: continue
        omp=engine[:-3]+'omp'
        for case in scenarios:
            separate={}
            for order in ('2','3','both'):
                for threads in (1,2):
                    model=cballs(); configure(model,engine,case,order,threads)
                    comm.Barrier(); model.Run(level=['MainLoop'])
                    timing=comm.allreduce(model.getTimings()['wall_seconds'],op=MPI.MAX)
                    error=None
                    try:
                        metadata=model.getRunMetadata()
                        CHECK.assertEqual(metadata['parallel']['estimator_ranks'],comm.size)
                        CHECK.assertEqual(metadata['parallel']['openmp_probe_threads'],threads)
                        if rank:
                            with CHECK.assertRaisesRegex(Exception,'no published native products'):
                                model.getResults()
                        else:
                            arrays=extract(model,order)
                            reference=cballs(); configure(reference,omp,case,order,1)
                            try:
                                reference.Run(level=['MainLoop']); expected=extract(reference,order)
                            finally: reference.struct_cleanup()
                            maximum=compare(arrays,expected)
                            if case.get('oracle'):
                                compare(arrays,ref.oracle(*case['catalogs'][0],use_log=case.get('log',False),
                                    smooth_radius=2*np.sin(.5*np.deg2rad(40/60)) if case.get('smooth') else None))
                            if not case.get('approximate'):
                                if threads==1 and order!='both':separate.update(arrays)
                                if order=='both':compare(arrays,separate)
                            tag=f'{engine}-{case["name"]}-{order}-t{threads}'
                            np.savez_compressed(directory/(tag+'.npz'),**arrays)
                            ill_conditioned=[]
                            if 'window' in arrays:
                                indices=np.arange(2*ref.NMAX+1)
                                for i in range(ref.BINS):
                                    for j in range(ref.BINS):
                                        matrix=arrays['window'][indices[:,None]-indices[None,:]+2*ref.NMAX,i,j]
                                        condition=np.linalg.cond(matrix)
                                        if condition>1e8 or not np.isfinite(condition):
                                            ill_conditioned.append(dict(bins=[i,j],condition=float(condition) if np.isfinite(condition) else None))
                            records.append(dict(engine=engine,case=case['name'],order=order,threads=threads,ranks=comm.size,
                                                max_absolute_omp_difference=maximum,max_rank_mainloop_seconds=timing,
                                                independent_oracle=bool(case.get('oracle')),unstable_corrected_bins=ill_conditioned,build=metadata['build']['id']))
                    except Exception: error=traceback.format_exc()
                    finally: model.struct_cleanup()
                    errors=comm.allgather(error)
                    if any(errors):raise RuntimeError('\n'.join(e for e in errors if e))
            if rank==0: print('MPI_SHEAR_CASE_OK',engine,case['name'],flush=True)
        for unsupported in ('legacy-one-ball','shear-pivot-reuse'):
            model=cballs();configure(model,engine,scenarios[0],'both',1)
            model.set(options='no-out-Hist,no-smooth-pivot,'+unsupported)
            try:
                with CHECK.assertRaisesRegex(Exception,unsupported):model.Run(level=['MainLoop'])
                CHECK.assertFalse(model.getBodytableAllocated())
                model.clean();configure(model,engine,scenarios[0],'both',1)
                model.Run(level=['MainLoop'])
            finally:model.struct_cleanup()
    if rank==0:
        CHECK.assertTrue(records)
        (directory/'coverage.json').write_text(json.dumps(dict(status='PASS',cases=records),indent=2)+'\n')
        print('MPI_SHEAR_OK',len(records),flush=True)


def native_cases(directory, launcher):
    """Exercise the executable's MPI-owned startup and shared-file publication."""
    import numpy as np
    import subprocess
    import test_shear_sphere_octree_omp as ref
    from cyballs import search_method_id
    directory.mkdir(exist_ok=True)
    pos,gamma,_=ref.fixture(60)
    catalog=directory/'catalog.txt'
    np.savetxt(catalog,np.column_stack((pos,gamma.real,gamma.imag)),fmt='%.17g')
    expected=ref.oracle(pos,gamma,np.ones(len(pos)))
    rows=[]
    for tree in ('octree','kdtree','balltree'):
        engine=tree+'-shear-sphere-2balls-mpi'
        if search_method_id(engine)<0:continue
        for order in ('2','3','both'):
            products={}
            for suffix in ('omp','mpi'):
                method=engine[:-3]+suffix
                out=directory/f'{method}-{order}';out.mkdir(exist_ok=True)
                options='pos-and-shear,no-smooth-pivot,no-one-ball,no-two-balls'
                if order!='both':options+=',only-'+order+'pcf'
                argv=[str(ROOT/'cballs'),f'searchMethod={method}',f'infile={catalog}',
                      'infileformat=multi-columns-ascii','columns=1,2,3,4,5',f'rootDir={out}',
                      f'options={options}','numberThreads=2','iCatalogs=1','usePeriodic=false','useLogHist=false',
                      f'rangeN={ref.RMAX}',f'rminHist={ref.RMIN}',f'sizeHistN={ref.BINS}',
                      f'mChebyshev={ref.NMAX}',f'sizeHistPhi={ref.PHI_BINS}','theta=0','nsmooth=8',
                      'verbose=0','verbose_log=0']
                if suffix=='mpi':argv=shlex.split(launcher)+['-n','2']+argv
                with (out/'process.log').open('w') as log:
                    result=subprocess.run(argv,cwd=out,stdout=log,stderr=subprocess.STDOUT,timeout=45)
                CHECK.assertEqual(result.returncode,0,(out/'process.log').read_text()[-6000:])
                products[suffix]={p.name:np.loadtxt(p) for p in out.glob('histShear*.txt')}
                CHECK.assertEqual(len(products[suffix]),1 if order=='2' else (2 if order=='3' else 3))
                if order!='3':
                    path=out/'histShearXi.txt'
                    CHECK.assertIn('full-sky spin-2',path.read_text().splitlines()[0])
                    a=products[suffix][path.name]
                    for value,target in [(a[:,1]+1j*a[:,2],expected['xi_plus']),
                                         (a[:,3]+1j*a[:,4],expected['xi_minus']),
                                         (a[:,5],expected['xi_weight'])]:
                        np.testing.assert_allclose(value,target,rtol=5e-10,atol=5e-12)
                if order!='2':
                    a=products[suffix]['histShearGammaMultipoles.txt']
                    for col,key in [(4,'multipoles'),(6,'upsilon')]:
                        np.testing.assert_allclose(a[:,col]+1j*a[:,col+1],expected[key].ravel(),rtol=5e-10,atol=5e-12)
            CHECK.assertEqual(set(products['omp']),set(products['mpi']))
            for key in products['omp']:
                np.testing.assert_allclose(products['omp'][key],products['mpi'][key],rtol=5e-10,atol=5e-12)
            rows.append(dict(engine=engine,order=order,ranks=2,threads=2,oracle=True))
    (directory/'coverage.json').write_text(json.dumps(dict(status='PASS',cases=rows),indent=2)+'\n')
    return rows


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--worker',action='store_true');p.add_argument('--output',type=Path,required=True)
    p.add_argument('--mpi-command',default='mpiexec',help='launcher without rank count')
    a=p.parse_args()
    if a.worker:return worker(a.output)
    from test_input_failure_contracts import run_child
    a.output.mkdir(parents=True,exist_ok=True)
    records=[]
    for ranks in (1,2,4):
        directory=a.output/f'ranks-{ranks}';directory.mkdir(exist_ok=True)
        run_child(shlex.split(a.mpi_command)+['-n',str(ranks),sys.executable,__file__,'--worker','--output',str(directory)],
                  directory,'MPI_SHEAR_OK',300)
        records+=json.loads((directory/'coverage.json').read_text())['cases']
    CHECK.assertTrue(records)
    native=native_cases(a.output/'executable',a.mpi_command)
    (a.output/'coverage.json').write_text(json.dumps(dict(status='PASS',cases=records,native_cases=native),indent=2)+'\n')
    print('PASS: MPI shear numerical and publication contracts',len(records))

if __name__=='__main__':main()
