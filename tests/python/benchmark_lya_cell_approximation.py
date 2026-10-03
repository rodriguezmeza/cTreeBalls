#!/usr/bin/env python3
"""Calibrate explicit cell geometry slop against the same catalog's exact 3PCF.

Retains raw sums, counts, per-bin errors, process wall/CPU time, peak RSS,
commands, hashes and run metadata. A slop value is NOT a correlation-error bound.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import platform
import re
import numpy as np
from benchmark_lya_triplet_kernels import accuracy,run_process,ROOT
from lya_corr_all_engines import read_npz,read_ascii,synthetic_catalog,save_catalog


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    source=p.add_mutually_exclusive_group(required=True)
    source.add_argument('--fits',nargs='+');source.add_argument('--catalog',type=Path)
    source.add_argument('--ascii',type=Path);source.add_argument('--synthetic',action='store_true')
    p.add_argument('--max-forests',type=int);p.add_argument('--pixel-stride',type=int,default=1)
    p.add_argument('--cballs',type=Path,default=ROOT/'cballs');p.add_argument('--outdir',type=Path,required=True)
    p.add_argument('--threads',type=int,default=1);p.add_argument('--warmups',type=int,default=1);p.add_argument('--repeats',type=int,default=3)
    p.add_argument('--r3-max',type=float,default=160);p.add_argument('--r3-bins',type=int,default=4)
    p.add_argument('--theta-bins',type=int,default=4);p.add_argument('--mu-bins',type=int,default=4)
    p.add_argument('--method',choices=['lya-3pcf-omp','lya-los-tree-3pcf-omp','lya-2pcf-3pcf-omp','lya-los-tree-2pcf-3pcf-omp'],default='lya-3pcf-omp')
    p.add_argument('--reference-kernel',type=int,choices=[0,1],default=1,help='1 independent direct; 0 previously validated exact segments for large catalogs')
    p.add_argument('--case',action='append',help='kernel:mu-slop:radial-slop:polar-slop[:pivot-cap], repeat for each case')
    p.add_argument('--relative-floor',type=float,default=1e-6);p.add_argument('--max-relative-error',type=float,default=.05)
    p.add_argument('--small-signal-atol',type=float,default=5e-8,help='absolute zeta error limit at or below relative-floor')
    p.add_argument('--require-accepted',action='store_true',help='fail unless a positive-slop case satisfies every accuracy check')
    args=p.parse_args(argv)
    if min(args.threads,args.repeats,args.r3_bins,args.theta_bins,args.mu_bins,args.pixel_stride)<1 or args.warmups<0:p.error('invalid count')
    if not np.isfinite(args.r3_max) or args.r3_max<=0:p.error('r3-max must be finite and positive')
    if args.max_forests is not None and args.max_forests<1:p.error('max-forests must be positive')
    for key in ['relative_floor','max_relative_error','small_signal_atol']:
        if not np.isfinite(getattr(args,key)) or getattr(args,key)<0:p.error(f'{key} must be finite and nonnegative')
    cases=[('exact',args.reference_kernel,0.,0.,0.,1)]
    for i,spec in enumerate(args.case or ['0:.01:0:0','0:.05:0:0','3:.01:0:0','4:.01:0:0','4:.01:.01:.01']):
        try:
            fields=spec.split(':');kernel=int(fields[0]);mu,radial,polar=map(float,fields[1:4]);cap=int(fields[4]) if len(fields)==5 else 8
            if len(fields) not in (4,5) or kernel not in (0,3,4,5) or not 1<=cap<=64 or not all(np.isfinite(v) and 0<=v<=1 for v in [mu,radial,polar]) or (kernel==0 and (radial or polar)):raise ValueError()
        except (ValueError,IndexError):p.error(f'invalid case {spec!r}')
        cases.append((f'case-{i}-k{kernel}-m{mu:g}-r{radial:g}-p{polar:g}-c{cap}',kernel,mu,radial,polar,cap))
    if args.fits:
        from lya_fits import read_fits
        catalog=read_fits(args.fits,max_forests=args.max_forests,pixel_stride=args.pixel_stride)
    elif args.max_forests is not None or args.pixel_stride!=1:p.error('selection requires FITS input')
    elif args.catalog:catalog=read_npz(args.catalog)
    elif args.ascii:catalog=read_ascii(args.ascii)
    else:catalog=synthetic_catalog()
    catalog=catalog.normalized()
    out=args.outdir.resolve();out.mkdir(parents=True,exist_ok=False)
    save_catalog(out/'catalog.npz',catalog);cat=out/'catalog.txt'
    with cat.open('w') as stream:
        for xyz,delta,weight,fid in zip(catalog.positions,catalog.delta,catalog.weights,catalog.forest_ids):
            stream.write(' '.join(format(v,'.17g') for v in (*xyz,delta,weight))+f' {int(fid)}\n')
    binary=args.cballs.resolve()
    report=dict(platform=platform.platform(),binary=str(binary),binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),
        catalog_sha256=hashlib.sha256(cat.read_bytes()).hexdigest(),pixels=catalog.nbody,forests=len(np.unique(catalog.forest_ids)),
        selection=catalog.metadata,settings=vars(args),timing_scope='complete native process, input/tree/search/output; reader preprocessing excluded',
        rss_scope='native process high-water RSS from wait4',runs=[],comparisons=[],accepted_cases=[])
    def save():(out/'summary.json').write_text(json.dumps(report,indent=2,default=str,allow_nan=False)+'\n')
    reference=None;reference_pairs=None;reference_count=None;baseline=None
    for label,kernel,mu,radial,polar,cap in cases:
        measurements=[];table=None;pair_table=None
        for rep in range(-args.warmups,args.repeats):
            folder=out/f'{label}-{rep}';folder.mkdir()
            params=dict(search=args.method,infile=cat,infileformat='lya-ascii',iCatalogs=1,rootDir=folder,
                numberThreads=args.threads,usePeriodic='false',useLogHist='false',rangeN=args.r3_max,rminHist=.1,sizeHistN=4,
                lya3RMax=args.r3_max,lya3RBins=args.r3_bins,lya3ThetaBins=args.theta_bins,lya3MuBins=args.mu_bins,
                lya3Kernel=kernel,lya3MuSlop=mu,lya3RadialSlop=radial,lya3PolarSlop=polar,lya3PivotCellMax=cap,
                verbose=2,verbose_log=1,options='no-smooth-pivot,lya-output-empty-bins')
            command=[str(binary),*(f'{key}={value}' for key,value in params.items())]
            print(f'Running {label}, repeat {rep}',flush=True)
            measurement=run_process(command,folder/'process.log')
            text=(folder/'process.log').read_text()
            for key in ['approximate_pairs','pivot_aggregates','nodes','pivot_tasks',
                        'evaluations','cache_hits','pruned_nodes','leaf_evaluations','pair_cache_hits']:
                match=re.search(rf'\b{key}=(\d+)',text)
                if match:measurement[key]=int(match[1])
            product=folder/'histZetaM_lya5d.txt';new=np.loadtxt(product,ndmin=2)
            if table is not None:np.testing.assert_array_equal(new,table)
            table=new
            count=int(re.search(r'ordered triplets: (\d+)',product.read_text())[1])
            if '2pcf-3pcf' in args.method:pair_table=np.loadtxt(folder/'histXi2pcf_lya.txt',ndmin=2)
            report['runs'].append(dict(case=label,repeat=rep,warmup=rep<0,command=command,**measurement));save()
            if rep>=0:measurements.append(measurement)
        wall=float(np.median([v['wall_seconds'] for v in measurements]))
        if reference is None:
            reference=table.copy();reference_count=count;reference_pairs=pair_table;baseline=wall
        if not np.isfinite(table).all():raise RuntimeError(f'Nonfinite output: {label}')
        metrics=accuracy(reference,table,args.relative_floor,args.max_relative_error)
        small=(reference[:,-1]>0)&(np.abs(reference[:,-3])<=args.relative_floor)
        absdiff=np.abs(table[:,-3]-reference[:,-3]);smallerr=float(absdiff[small].max()) if small.any() else 0.
        empty=reference[:,-1]==0;empty_ok=bool(np.all(table[empty,-1]==0))
        exact=(mu==radial==polar==0)
        if exact:np.testing.assert_allclose(table[:,-3:],reference[:,-3:],rtol=1e-9,atol=3e-11)
        if reference_pairs is not None:np.testing.assert_allclose(pair_table,reference_pairs,rtol=1e-11,atol=1e-8)
        metrics.update(case=label,kernel=kernel,mu_slop=mu,radial_slop=radial,polar_slop=polar,pivot_cap=cap,
            median_wall_seconds=wall,speedup_vs_exact=baseline/wall,median_process_cpu_seconds=float(np.median([v['process_cpu_seconds'] for v in measurements if v['process_cpu_seconds'] is not None])) if measurements[0]['process_cpu_seconds'] is not None else None,
            peak_rss_bytes=max(v['peak_rss_bytes'] or 0 for v in measurements),small_signal_max_absolute_error=smallerr,
            small_signal_atol=args.small_signal_atol,empty_bins_unchanged=empty_ok,ordered_triplets=count,counts_unchanged=count==reference_count,
            exact_raw_check_passed=exact)
        metrics['accepted']=bool(metrics['accepted'] and smallerr<=args.small_signal_atol and empty_ok and count==reference_count)
        report['comparisons'].append(metrics)
        if metrics['accepted'] and not exact:report['accepted_cases'].append(label)
        # Full per-bin evidence including small/empty bins; never hide cancellations.
        np.savetxt(out/f'{label}-errors.txt',np.column_stack([reference[:,:5],reference[:,-3:],table[:,-3:],absdiff]),
            header='b1 b2 t1 t2 mu ref_zeta ref_N ref_D cand_zeta cand_N cand_D abs_zeta_error')
        save();print(json.dumps(metrics),flush=True)
    if args.require_accepted and not report['accepted_cases']:raise SystemExit('No approximate case met all requested error checks; keep zero slop.')

if __name__=='__main__':main()
