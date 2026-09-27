#!/usr/bin/env python3
"""Retain exact/multipole accuracy, native process wall/CPU time and peak RSS.

Uses the same DESI/eBOSS/Cartesian FITS reader as lya_corr_all_engines.py.
No cyballs import is needed. Runs execute serially; each uses a fresh process.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import subprocess
import sys
import time
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
METHOD='lya-anisotropic-multipole-3pcf-omp'

def accuracy(reference,candidate,floor=1e-6,max_error=.05):
    np.testing.assert_array_equal(reference[:,:5],candidate[:,:5])
    n,d=reference[:,-2],reference[:,-1]
    cn,cd=candidate[:,-2],candidate[:,-1]
    eligible=(d>0)&(np.abs(reference[:,-3])>floor)
    valid=np.isfinite(candidate[:,-3])&(cd>0)
    good=eligible&valid
    errors=np.abs(candidate[good,-3]/reference[good,-3]-1)
    def relative_l2(a,b):
        norm=np.linalg.norm(a)
        return float(np.linalg.norm(b-a)/norm) if norm else (0. if np.array_equal(a,b) else None)
    maximum=float(errors.max()) if len(errors) else None
    return dict(numerator_relative_l2=relative_l2(n,cn),denominator_relative_l2=relative_l2(d,cd),
        eligible_bins=int(eligible.sum()),excluded_small_signal_bins=int(((d>0)&~eligible).sum()),
        invalid_eligible_bins=int((eligible&~valid).sum()),
        empty_reference_bins=int((d==0).sum()),empty_bin_reconstruction_l1=float(np.abs(cd[d==0]).sum()),
        max_relative_correlation_error=maximum,p95_relative_correlation_error=float(np.percentile(errors,95)) if len(errors) else None,
        relative_floor=floor,acceptance_limit=max_error,
        accepted=bool(eligible.any() and np.all(valid[eligible]) and maximum is not None and maximum<=max_error))

def run_process(command,log):
    start=time.perf_counter()
    with log.open('w') as stream:
        child=subprocess.Popen(command,stdout=stream,stderr=subprocess.STDOUT)
        if hasattr(os,'wait4'):
            _,status,usage=os.wait4(child.pid,0)
            child.returncode=os.waitstatus_to_exitcode(status)
            rss=usage.ru_maxrss*(1 if sys.platform=='darwin' else 1024)
            cpu=usage.ru_utime+usage.ru_stime
        else:
            child.wait();rss=cpu=None
    seconds=time.perf_counter()-start
    if child.returncode:raise RuntimeError(f'Native exit {child.returncode}; see {log}')
    text=log.read_text()
    match=re.search(r'cpusearch\s*=\s*([\d.eE+-]+)',text)
    result=dict(wall_seconds=seconds,process_cpu_seconds=cpu,peak_rss_bytes=rss,
                native_search_cpu_seconds=float(match[1]) if match else None)
    for key in ('aggregated_pairs','direct_pairs','segment_accepts'):
        found=re.search(rf'{key}=(\d+)',text)
        if found:result[key]=int(found[1])
    return result

def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    inputs=p.add_mutually_exclusive_group(required=True)
    inputs.add_argument('--fits',nargs='+');inputs.add_argument('--catalog',type=Path,help='NPZ saved by public forest reader')
    inputs.add_argument('--ascii',type=Path);inputs.add_argument('--synthetic',action='store_true')
    p.add_argument('--max-forests',type=int);p.add_argument('--pixel-stride',type=int,default=1)
    p.add_argument('--cballs',type=Path,default=ROOT/'cballs');p.add_argument('--outdir',type=Path,required=True)
    p.add_argument('--threads',type=int,default=1);p.add_argument('--warmups',type=int,default=1);p.add_argument('--repeats',type=int,default=3)
    p.add_argument('--r3-max',type=float,default=160);p.add_argument('--r3-bins',type=int,default=4)
    p.add_argument('--theta-bins',type=int,default=4);p.add_argument('--mu-bins',type=int,default=4)
    p.add_argument('--pivot-block',type=int,default=0,help='0 automatic, or 1..4096')
    p.add_argument('--lmax',type=int,nargs='*',default=[4,8,16],help='empty list runs exact kernels only')
    p.add_argument('--relative-floor',type=float,default=1e-6);p.add_argument('--max-relative-error',type=float,default=.05)
    p.add_argument('--require-accepted-multipole',action='store_true')
    args=p.parse_args(argv)
    if not 0<=args.pivot_block<=4096:p.error('pivot-block must be in 0..4096')
    if min(args.threads,args.repeats,args.r3_bins,args.theta_bins,args.mu_bins,args.pixel_stride)<1 or args.warmups<0:
        p.error('counts must be positive, warmups nonnegative')
    if not np.isfinite(args.r3_max) or args.r3_max<=0 or any(l<0 or l>32 for l in args.lmax):p.error('require r3-max>0 and lmax in 0..32')
    if args.max_forests is not None and args.max_forests<1:p.error('max-forests must be positive')
    if args.relative_floor<0 or not np.isfinite(args.relative_floor) or args.max_relative_error<0 or not np.isfinite(args.max_relative_error):p.error('accuracy limits must be finite and nonnegative')
    from lya_corr_all_engines import read_npz,read_ascii,synthetic_catalog,save_catalog
    if args.fits:
        from lya_fits import read_fits
        catalog=read_fits(args.fits,max_forests=args.max_forests,pixel_stride=args.pixel_stride)
    elif args.max_forests is not None or args.pixel_stride!=1:p.error('forest/stride selection requires FITS input')
    elif args.catalog:catalog=read_npz(args.catalog)
    elif args.ascii:catalog=read_ascii(args.ascii)
    else:catalog=synthetic_catalog()
    catalog=catalog.normalized()
    out=args.outdir.resolve();out.mkdir(parents=True,exist_ok=False)
    save_catalog(out/'catalog.npz',catalog)
    cat=out/'catalog.txt'
    # Integer IDs are formatted separately: never cast wide IDs through float.
    with cat.open('w') as stream:
        for xyz,delta,weight,fid in zip(catalog.positions,catalog.delta,catalog.weights,catalog.forest_ids):
            stream.write(' '.join(format(v,'.17g') for v in (*xyz,delta,weight))+f' {int(fid)}\n')
    binary=args.cballs.resolve()
    report=dict(platform=platform.platform(),binary=str(binary),binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),
        catalog_sha256=hashlib.sha256(cat.read_bytes()).hexdigest(),pixels=catalog.nbody,forests=len(np.unique(catalog.forest_ids)),
        selection=catalog.metadata,settings=vars(args).copy(),timing_scope='complete native process, including input/tree/search/output; preprocessing excluded',
        rss_scope='per-process high-water RSS from wait4; not Python reader memory',runs=[],comparisons=[])
    def save(): (out/'summary.json').write_text(json.dumps(report,indent=2,default=str,allow_nan=False)+'\n')
    cases=[('reference','lya-3pcf-omp',1,0),('tiled','lya-3pcf-omp',2,0),
           ('segments','lya-3pcf-omp',0,0),('los-segments','lya-los-tree-3pcf-omp',0,0)]
    cases.extend((f'multipole-L{l}',METHOD,0,l) for l in dict.fromkeys(args.lmax))
    reference=None
    for label,method,kernel,lmax in cases:
        timings=[];table=None
        for rep in range(-args.warmups,args.repeats):
            folder=out/f'{label}-{rep}';folder.mkdir()
            params=dict(search=method,infile=cat,infileformat='lya-ascii',iCatalogs=1,rootDir=folder,
                numberThreads=args.threads,usePeriodic='false',useLogHist='false',rangeN=args.r3_max,rminHist=.1,sizeHistN=4,
                lya3RMax=args.r3_max,lya3RBins=args.r3_bins,lya3ThetaBins=args.theta_bins,lya3MuBins=args.mu_bins,
                lya3Kernel=kernel,lya3LMax=lmax,lya3PivotBlock=args.pivot_block,verbose=2,verbose_log=1,options='no-smooth-pivot,lya-output-empty-bins')
            command=[str(binary),*(f'{k}={v}' for k,v in params.items())]
            measurement=run_process(command,folder/'process.log')
            filename='histZetaM_lya5d_multipole.txt' if method==METHOD else 'histZetaM_lya5d.txt'
            new=np.loadtxt(folder/filename,ndmin=2)
            if table is not None:np.testing.assert_allclose(new,table,rtol=0,atol=0,equal_nan=True)
            table=new
            report['runs'].append(dict(case=label,repeat=rep,warmup=rep<0,command=command,**measurement));save()
            if rep>=0:timings.append(measurement['wall_seconds'])
        if reference is None:reference=table.copy()
        metrics=accuracy(reference,table,args.relative_floor,args.max_relative_error)
        metrics.update(case=label,median_wall_seconds=float(np.median(timings)))
        if method!=METHOD:
            # Reassociation of hundreds of billions of contributions changes
            # rounding; raw sums must still agree tightly on the retained grid.
            np.testing.assert_allclose(table[:,-3:],reference[:,-3:],rtol=1e-9,atol=3e-11)
            metrics['exact_raw_check_passed']=True
        report['comparisons'].append(metrics);save();print(json.dumps(metrics),flush=True)
    accepted=[r['case'] for r in report['comparisons'] if r['case'].startswith('multipole') and r['accepted']]
    report['accepted_multipoles']=accepted;save()
    if args.require_accepted_multipole and not accepted:raise SystemExit('No multipole order met the requested accuracy; exact kernels remain available.')

if __name__=='__main__':main()
