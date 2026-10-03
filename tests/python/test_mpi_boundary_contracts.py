#!/usr/bin/env python3
"""Enumerate and fault every reached application boundary of active MPI routes.

Each named boundary is failed on each rank, under owned and borrowed MPI,
followed by same-object recovery. This is finite application-error coverage;
process death, failed communicators and divergent calls made outside Run are not
recoverable MPI contracts. Assertions use unittest and survive Python -O.
"""
import argparse
import json
import os
from pathlib import Path
import shlex
import sys
import unittest
ROOT=Path(__file__).resolve().parents[2]
sys.path[:0]=[str(ROOT),str(ROOT/'scripts'),str(ROOT/'tests/python')]
from test_input_failure_contracts import run_child
CHECK=unittest.TestCase()


def configure(model,engine,directory,variant='ordinary',bad=None):
    import numpy as np
    from capabilities_generated import ENGINES
    model.clean_all()
    family=ENGINES[engine]['gate']['oracle']
    p=dict(searchMethod=engine,rootDir=str(directory/'out'),numberThreads=1,verbose=0,verbose_log=0,
           iCatalogs='1',usePeriodic=False,useLogHist=False,sizeHistN=4,mChebyshev=2,nsmooth=2,
           sizeHistPhi=8,rangeN=1.5,rminHist=.02,theta=0.,options='no-smooth-pivot,no-one-ball,no-two-balls,no-out-Hist')
    if family=='forest':
        import test_lya_forest_mpi as f
        same='same-los' in engine
        canonical=engine.replace('lya-los-tree-','lya-').rsplit('-',1)[0]
        kind=6 if same else f.METHODS.index(canonical)
        data=np.array(f.three.POINTS if kind<3 else f.radial.WIDE_ANGLE,copy=True)
        if bad=='catalog':data[0,3]+=.125
        p.update(f.params(kind,'',directory/'out',1));p.pop('infile');p.pop('infileformat')
        p['searchMethod']=engine
        if bad=='parameter':p['sizeHistN']=0
        model.set(p);model.set_forest_catalog(data[:,:3],data[:,3],data[:,4],data[:,5].astype(np.int64))
    elif family=='shear':
        from test_shear_sphere_octree_omp import fixture
        pos,gamma,weights=fixture()
        p['rminHist']=.03
        if bad=='catalog':gamma=gamma.copy();gamma[0]+=.125j
        if variant=='pairs':p['options']+=',only-2pcf'
        if variant=='triplets':p['options']+=',only-3pcf'
        if variant=='approximate':
            p['options']='no-smooth-pivot,no-out-Hist';p['theta']=.1
        if variant=='smooth':
            p['options']='smooth-pivot,no-one-ball,no-out-Hist';p['rsmooth']='40'
        if bad=='parameter':p['sizeHistN']=0
        if variant=='cross':p['iCatalogs']='1,2,3'
        model.set(p);model.set_catalog(pos,weights=weights,gamma1=gamma.real,gamma2=gamma.imag)
        if variant=='cross':
            for index in (1,2):
                q,g,w=fixture(31+index)
                model.set_catalog(q,weights=w,gamma1=g.real,gamma2=g.imag,catalog=index)
    elif family=='physical':
        import test_octree_3pcf_3d_omp as physical
        data=np.array(physical.CATALOG)
        if bad=='catalog':data[0,3]+=.125
        p.update(rangeN=physical.RMAX,rminHist=physical.RMIN,sizeHistN=physical.NBINS,mChebyshev=physical.LMAX,
                 options='compute-2pcf-3d,compute-3pcf-3d,no-smooth-pivot')
        if bad=='parameter':p['sizeHistN']=0
        model.set(p);model.set_catalog(data[:,:3],kappa=data[:,3],weights=data[:,4])
    else:
        from test_two_ball_edge_corrections import catalog
        pos,k,w,mask=catalog()
        if bad=='catalog':k=k.copy();k[0]+=.125
        p['options']+=',weights-norm,KKKCorrelation,no-normalize-HistZeta'
        if variant=='compatibility':p['options']=p['options'].replace(',no-two-balls','')+',legacy-one-ball'
        if variant=='edge':p['options']+=',edge-corrections'
        if bad=='parameter':p['sizeHistN']=0
        model.set(p);model.set_catalog(pos,kappa=k,weights=w,mask=mask)
    if bad=='controls':model.set(sizeHistN=5)
    if bad=='unread':model.set(unknown_parameter=1)
    if bad=='memory-conflict':model.set(infile='cannot-mix-file-and-memory')
    return p


def worker(directory,ownership):
    if ownership=='borrowed':
        from mpi4py import MPI
    else:
        import mpi4py
        mpi4py.rc.initialize=False
    import numpy as np
    from cyballs import cballs,search_method_id
    from capabilities_generated import ENGINES
    rank=int(os.environ.get('OMPI_COMM_WORLD_RANK','0'))
    model=cballs();records=[]
    for engine in ENGINES:
        if not engine.endswith('-mpi') or search_method_id(engine)<0:continue
        variants=['ordinary']
        if engine in ('kdtree-2balls-mpi','balltree-2balls-mpi','octree-2balls-mpi'):
            variants+=['compatibility','edge']
        if ENGINES[engine]['family']=='shear':
            variants+=['pairs','triplets','approximate','smooth','cross']
        for variant in variants:
            configure(model,engine,directory,variant)
            model._mpi_test_boundary();model.Run()
            from mpi4py import MPI
            comm=MPI.COMM_WORLD
            # Include the explicit finalization stage, without finalizing MPI.
            model.Run(level=['EndRun'])
            names=sorted(set(model._mpi_boundary_trace()))
            names=sorted(set().union(*map(set,comm.allgather(names))))
            CHECK.assertIn('MPI Python input preflight',names)
            CHECK.assertIn('MPI startup setup',names)
            CHECK.assertIn('MPI catalog loading',names)
            CHECK.assertIn('MPI catalog geometry and bins',names)
            CHECK.assertIn('MPI computation completion',names)
            CHECK.assertIn('MPI run finalization',names)
            model.struct_cleanup()
            for failing_rank in range(comm.size):
                for boundary in names:
                    configure(model,engine,directory,variant)
                    model._mpi_test_boundary(boundary if rank==failing_rank else '')
                    with CHECK.assertRaises(Exception,msg=(engine,variant,rank,boundary)) as caught:
                        model.Run(level=['EndRun'])
                    messages=comm.allgather(str(caught.exception))
                    CHECK.assertEqual(messages[0],messages[1],(engine,variant,boundary,messages))
                    CHECK.assertIn('injected application failure',messages[0],(engine,variant,boundary,messages))
                    CHECK.assertIn(boundary,messages[0])
                    CHECK.assertFalse(model.getBodytableAllocated())
                    CHECK.assertFalse(model._runtime_bodytable_address())
                    CHECK.assertFalse(MPI.Is_finalized())
                    configure(model,engine,directory,variant)
                    model._mpi_test_boundary();model.Run();model.struct_cleanup()
                    comm.Barrier()
                # Real failures before native startup, rather than injected flags.
                for bad in ('parameter','unread','memory-conflict','controls','catalog','level'):
                    configure(model,engine,directory,variant,bad if rank==failing_rank else None)
                    model._mpi_test_boundary()
                    with CHECK.assertRaises(Exception) as caught:
                        model.Run(level=['invalid-stage'] if bad=='level' and rank==failing_rank else ['MainLoop'])
                    messages=comm.allgather(str(caught.exception))
                    CHECK.assertEqual(messages[0],messages[1])
                    CHECK.assertFalse(model.getBodytableAllocated())
                    configure(model,engine,directory,variant);model.Run();model.struct_cleanup();comm.Barrier()
            # Cached calls must still detect a divergent rank's changed state.
            configure(model,engine,directory,variant);model.Run()
            if rank==1:model.set(rootDir=str(directory/'changed-output'))
            with CHECK.assertRaises(Exception) as caught:model.Run()
            messages=comm.allgather(str(caught.exception))
            CHECK.assertEqual(messages[0],messages[1]);CHECK.assertIn('MPI Python run state',messages[0])
            CHECK.assertFalse(model.getBodytableAllocated())
            configure(model,engine,directory,variant);model.Run();model.struct_cleanup();comm.Barrier()
            records.append(dict(engine=engine,variant=variant,ownership=ownership,boundaries=names,
                                failing_ranks=list(range(comm.size)),fault_cases=len(names)*comm.size,real_cases=6*comm.size+1))
            if rank==0:print('MPI_BOUNDARY_ROUTE_OK',engine,variant,len(names),flush=True)
    if rank==0:
        configure(model,'octree-2balls-omp',directory);model.Run();model.struct_cleanup()
    comm.Barrier()
    model.clean_all()
    if rank==0:
        (directory/'coverage.json').write_text(json.dumps(dict(status='PASS',routes=records),indent=2)+'\n')
        print('MPI_BOUNDARIES_OK',ownership,sum(x['fault_cases']+x['real_cases'] for x in records),flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--worker',choices=('owned','borrowed'));p.add_argument('--directory',type=Path)
    p.add_argument('--output',required=False,type=Path);p.add_argument('--mpi-command')
    a=p.parse_args()
    if a.worker:return worker(a.directory,a.worker)
    a.output.mkdir(parents=True,exist_ok=True)
    reports=[]
    for ownership in ('owned','borrowed'):
        directory=a.output/ownership;directory.mkdir(exist_ok=True)
        run_child(shlex.split(a.mpi_command)+[sys.executable,__file__,'--worker',ownership,'--directory',str(directory)],
                  directory,'MPI_BOUNDARIES_OK',600)
        reports+=json.loads((directory/'coverage.json').read_text())['routes']
    (a.output/'coverage.json').write_text(json.dumps(dict(status='PASS',routes=reports,
        scope='All reached named application boundaries for every active MPI method, native scalar, compatibility and edge routes, on either rank, owned/borrowed MPI, with same-object recovery.'),indent=2)+'\n')
if __name__=='__main__':main()
