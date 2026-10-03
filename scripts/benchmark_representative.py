#!/usr/bin/env python3
"""Representative scalar scaling with retained products and matching exact runs.

Sequential fresh workers, nested seeded full-sky/clustered catalogs, cold and
warm calls, per-rank RSS and max-wall/sum-CPU timing. Threads are per rank.
Approximation qualification is explicitly scoped to each measured catalog.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import shlex
import subprocess
import sys
import time
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT),str(ROOT/'scripts')]


def fixture(n,geometry):
    import numpy as np
    rng=np.random.default_rng(782104)
    if geometry=='clustered':
        centers=rng.normal(size=(128,3));centers/=np.linalg.norm(centers,axis=1)[:,None]
        p=centers[np.arange(n)%len(centers)]+np.random.default_rng(774).normal(scale=2e-4,size=(n,3))
    else:p=rng.normal(size=(n,3))
    p/=np.linalg.norm(p,axis=1)[:,None]
    k=.3+p[:,0]-.2*p[:,1]+np.random.default_rng(73).normal(scale=.02,size=n)
    w=np.random.default_rng(99).uniform(.5,1.5,n)
    mask=(np.arange(n)%11!=0).astype(np.uint8)
    return p,k,w,mask


def worker(a):
    import numpy as np
    from benchmark_contracts import peak_bytes
    from cyballs import cballs,save_result_packet,build_info
    comm=None
    if a.engine.endswith('-mpi'):
        from mpi4py import MPI
        comm=MPI.COMM_WORLD
    rank=comm.rank if comm else 0
    if rank==0:a.output.mkdir(parents=True,exist_ok=False)
    if comm:comm.Barrier()
    p,k,w,mask=fixture(a.count,a.geometry)
    if rank==0:np.savez_compressed(a.output/'catalog.npz',positions=p,kappa=k,weights=w,mask=mask)
    theta={'kdtree':.5,'balltree':.2,'octree':.025}[a.engine.split('-')[0]]
    options='KKKCorrelation,no-normalize-HistZeta,weights-norm,no-smooth-pivot,no-out-Hist,read-mask,dual-node-profile'
    if a.exact:options+=',no-one-ball,no-two-balls';theta=0.
    model=cballs();records=[]
    model.set(searchMethod=a.engine,sizeHistN=6,mChebyshev=5,useLogHist=True,usePeriodic=False,
              rminHist=.02,rangeN=1.8,nsmooth=8,theta=theta,numberThreads=a.threads,
              verbose=0,verbose_log=0,options=options)
    model.set_catalog(p,kappa=k,weights=w,mask=mask)
    try:
        for cache in ('cold','warm'):
            model.set(rootDir=str(a.output/cache))
            if comm:comm.Barrier()
            start=time.perf_counter();cpu=time.process_time()
            model.Run()
            timing=dict(rank=rank,wall=time.perf_counter()-start,cpu=time.process_time()-cpu,
                        peak_rss_bytes=peak_bytes(),native=model.getTimings(),cache=model.getCacheInfo())
            ranks=comm.allgather(timing) if comm else [timing]
            if rank==0:
                packet=model.getResults();save_result_packet(packet,a.output/cache/'packet')
                row=dict(cache=cache,max_rank_wall_seconds=max(r['wall'] for r in ranks),
                    sum_rank_cpu_seconds=sum(r['cpu'] for r in ranks),
                    max_rank_peak_rss_bytes=max(r['peak_rss_bytes'] for r in ranks),ranks=ranks)
                records.append(row)
        if rank==0:
            from cyballs import load_result_packet
            cold=load_result_packet(a.output/'cold/packet');warm=load_result_packet(a.output/'warm/packet')
            for name,target in cold['arrays'].items():
                np.testing.assert_allclose(warm['arrays'][name],target,rtol=2e-12,atol=1e-10,equal_nan=True)
            result=dict(status='PASS',engine=a.engine,count=a.count,geometry=a.geometry,threads_per_rank=a.threads,
                 mpi_ranks=comm.size if comm else 1,exact=a.exact,records=records,build=build_info(),
                 catalog_sha256=hashlib.sha256((a.output/'catalog.npz').read_bytes()).hexdigest(),
                 rss_semantics='Per-process high-water RSS including Python/imports; max across ranks, not node resident sum',
                 timing_scope='Run including native startup, tree/search and provenance; excludes registration, extraction and packet writes')
            (a.output/'sample.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    finally:model.struct_cleanup()


def main(a):
    from cyballs import qualification_report,load_result_packet
    a.output.mkdir(parents=True,exist_ok=False)
    report=dict(schema_version=1,status='RUNNING',samples=[],host=dict(platform=platform.platform(),machine=platform.machine(),logical_cpus=os.cpu_count()),
                methodology='One fresh worker at a time; nested deterministic catalogs; full-catalog exact reference for every N and geometry; cold and warm calls; no speedup acceptance threshold.')
    def save():(a.output/'scaling.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    save()
    for geometry in ('full-sky','clustered'):
        for n in a.counts:
            reference=None
            # Independent exact traversal of the full measured catalog is kept.
            matrix=[('kdtree-2balls-omp',1,1,True)]
            for stem in ('kdtree','balltree','octree'):
                matrix += [(stem+'-2balls-omp',t,1,False) for t in a.thread_counts]
                if a.mpi_command:matrix += [(stem+'-2balls-mpi',1,2,False)]
            for engine,threads,ranks,exact in matrix:
                label=f'{geometry}-{n}-{engine}-t{threads}-r{ranks}'+('-exact' if exact else '')
                directory=a.output/label
                command=[sys.executable,str(Path(__file__).resolve()),'--worker','--engine',engine,'--count',str(n),
                         '--geometry',geometry,'--threads',str(threads),'--output',str(directory)]
                if exact:command+=['--exact']
                if ranks>1:command=shlex.split(a.mpi_command)+['-n',str(ranks)]+command
                with (a.output/(label+'.log')).open('w') as stream:
                    done=subprocess.run(command,cwd=ROOT,stdout=stream,stderr=subprocess.STDOUT,timeout=600)
                if done.returncode:raise RuntimeError('scaling worker failed: '+label)
                sample=json.loads((directory/'sample.json').read_text())
                packet=load_result_packet(directory/'cold/packet')
                if exact:reference=packet
                sample['qualification']=qualification_report(packet,reference,rtol=.02,atol=1e-10)
                sample['directory']=label
                report['samples'].append(sample);save()
                print(label,sample['records'][0]['max_rank_wall_seconds'],sample['qualification']['status'],flush=True)
    report['status']='PASS';save()

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',required=True,type=Path);p.add_argument('--worker',action='store_true')
    p.add_argument('--engine');p.add_argument('--count',type=int);p.add_argument('--geometry')
    p.add_argument('--threads',type=int,default=1);p.add_argument('--exact',action='store_true')
    p.add_argument('--counts',type=int,nargs='+',default=[2048,8192,32768])
    p.add_argument('--thread-counts',type=int,nargs='+',default=[1,2,4]);p.add_argument('--mpi-command')
    a=p.parse_args();a.output=a.output.resolve()
    worker(a) if a.worker else main(a)
