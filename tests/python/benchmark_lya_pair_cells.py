#!/usr/bin/env python3
"""Exact/slop Ly-alpha pair-cell calibration with retained CPU, wall, RSS and raw sums."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import platform
import re
import numpy as np
from benchmark_lya_triplet_kernels import run_process,ROOT
from lya_corr_all_engines import read_npz,read_ascii,synthetic_catalog,save_catalog


def compare(reference,candidate,floor,limit,atol):
    np.testing.assert_array_equal(reference[:,:4],candidate[:,:4])
    x,y=reference[:,4],candidate[:,4]
    d,cd=reference[:,6],candidate[:,6]
    eligible=(d>0)&(np.abs(x)>floor);small=(d>0)&~eligible
    relative=np.abs(y[eligible]-x[eligible])/np.abs(x[eligible])
    maximum=float(relative.max()) if relative.size else None
    small_error=float(np.max(np.abs(x[small]-y[small]))) if small.any() else 0.
    finite=bool(np.isfinite(candidate).all())
    occupancy=bool(np.array_equal(d>0,cd>0))
    raw=bool(np.allclose(reference[:,-2:],candidate[:,-2:],rtol=3e-11,atol=1e-9))
    return dict(eligible_bins=int(eligible.sum()),small_signal_bins=int(small.sum()),
        empty_bins=int((d==0).sum()),finite=finite,same_occupancy=occupancy,
        max_relative_correlation_error=maximum,small_signal_max_absolute_error=small_error,
        raw_sums_match=raw,accepted=finite and occupancy and (maximum is None or maximum<=limit) and small_error<=atol)


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    source=p.add_mutually_exclusive_group(required=True)
    source.add_argument('--fits',nargs='+');source.add_argument('--catalog',type=Path)
    source.add_argument('--ascii',type=Path);source.add_argument('--synthetic',action='store_true')
    p.add_argument('--max-forests',type=int);p.add_argument('--pixel-stride',type=int,default=1)
    p.add_argument('--cballs',type=Path,default=ROOT/'cballs')
    p.add_argument('--baseline-cballs',type=Path,help='optional previous binary for legacy reference/timing')
    p.add_argument('--outdir',type=Path,required=True)
    p.add_argument('--threads',type=int,default=1);p.add_argument('--warmups',type=int,default=1);p.add_argument('--repeats',type=int,default=3)
    p.add_argument('--rp-max',type=float,default=160);p.add_argument('--rt-max',type=float,default=160)
    p.add_argument('--rp-bins',type=int,default=50);p.add_argument('--rt-bins',type=int,default=50)
    p.add_argument('--method',choices=['lya-2pcf-omp','lya-los-tree-2pcf-omp','lya-2pcf-3pcf-omp','lya-los-tree-2pcf-3pcf-omp'],default='lya-2pcf-omp')
    p.add_argument('--case',action='append',help='rp-slop:rt-slop; repeat to calibrate; defaults to exact 0:0 only')
    p.add_argument('--kernel3',type=int,choices=range(5),default=0)
    p.add_argument('--r3-max',type=float,default=160);p.add_argument('--r3-bins',type=int,default=4)
    p.add_argument('--theta-bins',type=int,default=4);p.add_argument('--mu-bins',type=int,default=4)
    p.add_argument('--relative-floor',type=float,default=1e-6);p.add_argument('--max-relative-error',type=float,default=.05)
    p.add_argument('--small-signal-atol',type=float,default=5e-8)
    p.add_argument('--require-accepted',action='store_true',help='fail unless a requested positive-slop case passes')
    args=p.parse_args(argv)
    if min(args.threads,args.repeats,args.rp_bins,args.rt_bins,args.r3_bins,args.theta_bins,args.mu_bins,args.pixel_stride)<1 or args.warmups<0:p.error('invalid count')
    if any(not np.isfinite(v) or v<=0 for v in [args.rp_max,args.rt_max,args.r3_max]):p.error('maxima must be positive and finite')
    if any(not np.isfinite(v) or v<0 for v in [args.relative_floor,args.max_relative_error,args.small_signal_atol]):p.error('invalid accuracy threshold')
    if args.max_forests is not None and args.max_forests<1:p.error('max-forests must be positive')
    cases=[('reference',0,0.,0.)]
    for i,spec in enumerate(args.case or ['0:0']):
        try:
            rp,rt=map(float,spec.split(':'))
            if not all(np.isfinite(x) and 0<=x<=1 for x in [rp,rt]):raise ValueError()
        except ValueError:p.error(f'invalid case {spec!r}')
        cases.append((f'case-{i}-rp{rp:g}-rt{rt:g}',1,rp,rt))
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
    binary=args.cballs.resolve();baseline=(args.baseline_cballs or binary).resolve()
    report=dict(platform=platform.platform(),binaries={str(b):hashlib.sha256(b.read_bytes()).hexdigest() for b in [binary,baseline]},
        catalog_sha256=hashlib.sha256(cat.read_bytes()).hexdigest(),pixels=catalog.nbody,forests=len(np.unique(catalog.forest_ids)),
        selection=catalog.metadata,settings=vars(args),timing_scope='complete native process; FITS preprocessing excluded',
        rss_scope='native process high-water RSS from wait4',runs=[],comparisons=[],accepted_cases=[])
    def save():(out/'summary.json').write_text(json.dumps(report,indent=2,default=str,allow_nan=False)+'\n')
    reference=None;reference3=None;reference_count=None;reference_wall=None;exact_failure=False
    results={label:dict(measurements=[],table=None,table3=None) for label,_,_,_ in cases}
    # Interleave reference/candidates within each repetition to reduce drift.
    for rep in range(-args.warmups,args.repeats):
        for label,kernel,rp,rt in cases:
            saved=results[label];table=saved['table'];table3=saved['table3']
            folder=out/f'{label}-{rep}';folder.mkdir()
            params=dict(search=args.method,infile=cat,infileformat='lya-ascii',iCatalogs=1,rootDir=folder,
                numberThreads=args.threads,usePeriodic='false',useLogHist='false',rangeN=args.r3_max,rminHist=.1,sizeHistN=4,
                lya2RpMax=args.rp_max,lya2RtMax=args.rt_max,lya2RpBins=args.rp_bins,lya2RtBins=args.rt_bins,
                lya3RMax=args.r3_max,lya3RBins=args.r3_bins,lya3ThetaBins=args.theta_bins,lya3MuBins=args.mu_bins,lya3Kernel=args.kernel3,
                verbose=2,verbose_log=1,options='no-smooth-pivot,lya-output-empty-bins')
            # Omit new controls for old binaries; their pair traversal is kernel 0.
            if kernel:params.update(lya2Kernel=kernel,lya2RpSlop=rp,lya2RtSlop=rt)
            command=[str(binary if kernel else baseline),*(f'{key}={value}' for key,value in params.items())]
            print(f'Running {label}, repeat {rep}',flush=True)
            measurement=run_process(command,folder/'process.log')
            text=(folder/'process.log').read_text()
            pair_line=next((line for line in text.splitlines() if line.startswith('Ly-alpha 2PCF cells:')),'')
            for key in ('aggregated_pairs','direct_pairs','segment_accepts'):measurement.pop(key,None)
            measurement['pair_work']={key:int(value) for key,value in re.findall(r'(\w+)=(\d+)(?: |$)',pair_line)}
            for label_rx,field in [('forest build:.*build_CPU','pair_build_cpu_seconds'),('cells:.*search_CPU','pair_search_cpu_seconds')]:
                match=re.search(rf'Ly-alpha 2PCF {label_rx}=([\d.eE+-]+)',text)
                measurement[field]=float(match[1]) if match else None
            product=folder/'histXi2pcf_lya.txt';new=np.loadtxt(product,ndmin=2)
            if table is not None:np.testing.assert_array_equal(new,table)
            table=new
            count=int(re.search(r'unordered pairs: (\d+)',product.read_text())[1])
            if '2pcf-3pcf' in args.method:
                new3=np.loadtxt(folder/'histZetaM_lya5d.txt',ndmin=2)
                if table3 is not None:np.testing.assert_array_equal(new3,table3)
                table3=new3
            report['runs'].append(dict(case=label,repeat=rep,warmup=rep<0,command=command,**measurement));save()
            saved.update(table=table,table3=table3,count=count)
            if rep>=0:saved['measurements'].append(measurement)
    for label,kernel,rp,rt in cases:
        saved=results[label];measurements=saved['measurements']
        table=saved['table'];table3=saved['table3'];count=saved['count']
        wall=float(np.median([v['wall_seconds'] for v in measurements]))
        if reference is None:
            reference=table.copy();reference3=table3;reference_count=count;reference_wall=wall
        result=compare(reference,table,args.relative_floor,args.max_relative_error,args.small_signal_atol)
        result['same_pair_count']=count==reference_count
        result['three_point_unchanged']=table3 is None or bool(np.allclose(reference3,table3,rtol=3e-11,atol=1e-9))
        result['accepted'] &= result['same_pair_count'] and result['three_point_unchanged']
        if not rp and not rt and not result['raw_sums_match']:result['accepted']=False
        if kernel and not rp and not rt and not result['accepted']:exact_failure=True
        error=np.column_stack((reference[:,:4],reference[:,4],table[:,4],np.abs(reference[:,4]-table[:,4]),reference[:,-2:],table[:,-2:]))
        np.savetxt(out/f'{label}-per-bin.txt',error,header='bp bt rp rt xi_ref xi_test abs_error N_ref D_ref N_test D_test')
        cpu=[v['process_cpu_seconds'] for v in measurements if v['process_cpu_seconds'] is not None]
        rss=[v['peak_rss_bytes'] for v in measurements if v['peak_rss_bytes'] is not None]
        result.update(case=label,kernel=kernel,rp_slop=rp,rt_slop=rt,pairs=count,median_wall_seconds=wall,
            median_cpu_seconds=float(np.median(cpu)) if cpu else None,peak_rss_bytes=max(rss) if rss else None,speedup_vs_reference=reference_wall/wall)
        report['comparisons'].append(result)
        if kernel and result['accepted']:report['accepted_cases'].append(label)
        save();print(json.dumps(result),flush=True)
    if exact_failure:return 1
    if args.require_accepted and not any(c['kernel'] and (c['rp_slop'] or c['rt_slop']) and c['accepted'] for c in report['comparisons']):return 2
    return 0

if __name__=='__main__':raise SystemExit(main())
