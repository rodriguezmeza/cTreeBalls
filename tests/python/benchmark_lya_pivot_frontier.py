#!/usr/bin/env python3
"""Calibrate forest-local pivot smoothing against exact pixels, with wall/CPU/RSS.

Each case retains its complete products and provenance. Positive radius is a
geometry approximation, not a promised relative-error tolerance. No Cython is
required; input readers and native process accounting are shared with the other
public Ly-alpha benchmarks.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import numpy as np
from benchmark_lya_triplet_kernels import ROOT,run_process
from lya_corr_all_engines import read_npz,read_ascii,synthetic_catalog,save_catalog

METHODS=tuple(f'{prefix}-{order}-omp' for prefix in ('lya','lya-los-tree')
              for order in ('2pcf','3pcf','2pcf-3pcf'))

def comparison(reference,candidate,floor,limit,atol,exact):
    # Both products are emitted densely for calibration, including empty bins.
    np.testing.assert_array_equal(reference[:,:-3],candidate[:,:-3])
    x,y=reference[:,-3],candidate[:,-3];d,cd=reference[:,-1],candidate[:,-1]
    eligible=(d>0)&(np.abs(x)>floor);small=(d>0)&~eligible
    errors=np.abs(y[eligible]-x[eligible])/np.abs(x[eligible])
    finite=bool(np.isfinite(candidate).all());occupancy=bool(np.array_equal(d>0,cd>0))
    def finite_value(value):
        return float(value) if np.isfinite(value) else None
    maximum=finite_value(errors.max()) if errors.size else None
    small_error=float(np.abs(y[small]-x[small]).max()) if small.any() else 0.
    raw=bool(np.allclose(reference[:,-2:],candidate[:,-2:],rtol=3e-11,atol=1e-9))
    relative_l2=[]
    for col in (-2,-1):
        scale=np.linalg.norm(reference[:,col]);error=np.linalg.norm(candidate[:,col]-reference[:,col])
        relative_l2.append(finite_value(error/scale) if scale else (0. if error==0 else None))
    accepted=(finite and occupancy and small_error<=atol and
              ((raw and (maximum is None or maximum<=limit)) if exact
               else (maximum is not None and maximum<=limit)))
    return dict(accepted=bool(accepted),finite=finite,same_occupancy=occupancy,raw_sums_match=raw,
        eligible_bins=int(eligible.sum()),small_signal_bins=int(small.sum()),empty_bins=int((d==0).sum()),
        max_relative_correlation_error=maximum,
        p95_relative_correlation_error=finite_value(np.percentile(errors,95)) if errors.size else None,
        small_signal_max_absolute_error=finite_value(small_error),numerator_relative_l2=relative_l2[0],denominator_relative_l2=relative_l2[1],
        newly_occupied_bins=int(((d==0)&(cd>0)).sum()),lost_occupied_bins=int(((d>0)&(cd==0)).sum()))

def plot_summary(report,out):
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        (out/'plot-notes.txt').write_text('Install matplotlib to generate the summary figures. Numeric results are in summary.json.\n');return
    for method in report['methods']:
        rows=[x for x in report['comparisons'] if x['method']==method];labels=[x['case'] for x in rows];x=np.arange(len(rows))
        fig,axes=plt.subplots(2,2,figsize=(max(10,len(rows)*1.6),7.5),layout='constrained')
        walls=[r['median_wall_seconds'] for r in rows]
        axes[0,0].bar(x,walls,color='#17698c');axes[0,0].set_ylabel('Whole-process wall seconds')
        axes[0,1].bar(x,[r['wall_speedup'] for r in rows],color='#23856d');axes[0,1].axhline(1,color='gray',ls='--');axes[0,1].set_ylabel('Wall speedup versus exact reference')
        for order in (2,3):
            errors=[r['products'].get(str(order),{}).get('max_relative_correlation_error') for r in rows]
            if any(v is not None for v in errors):axes[1,0].plot(x,[np.nan if v is None else 100*v for v in errors],marker='o',label=f'{order}PCF')
        axes[1,0].axhline(report['settings']['max_relative_error']*100,color='#ac3c3c',ls='--',label='Acceptance limit')
        axes[1,0].set_ylabel('Maximum relative error (%)');axes[1,0].set_yscale('symlog',linthresh=.01);axes[1,0].set_ylim(bottom=0);axes[1,0].legend()
        axes[1,1].bar(x,[np.nan if r['peak_rss_bytes'] is None else r['peak_rss_bytes']/2**20 for r in rows],color='#6b6498');axes[1,1].set_ylabel('Maximum process peak RSS (MiB)')
        for ax in axes.flat:ax.set_xticks(x,labels,rotation=30,ha='right');ax.grid(axis='y',alpha=.18)
        fig.suptitle(f'{method}: {report["pixels"]:,} pixels, {report["forests"]:,} forests, {report["settings"]["threads"]} threads')
        fig.savefig(out/f'{method}-calibration.png',dpi=180);plt.close(fig)

def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    source=p.add_mutually_exclusive_group(required=True)
    source.add_argument('--fits',nargs='+');source.add_argument('--catalog',type=Path)
    source.add_argument('--ascii',type=Path);source.add_argument('--synthetic',action='store_true')
    p.add_argument('--max-forests',type=int);p.add_argument('--pixel-stride',type=int,default=1)
    p.add_argument('--cballs',type=Path,default=ROOT/'cballs');p.add_argument('--baseline-cballs',type=Path)
    p.add_argument('--methods',choices=METHODS,nargs='+',default=list(METHODS[:3]))
    p.add_argument('--outdir',type=Path,required=True)
    p.add_argument('--threads',type=int,default=1);p.add_argument('--warmups',type=int,default=1);p.add_argument('--repeats',type=int,default=3)
    p.add_argument('--scan-level',type=int,default=2);p.add_argument('--pivot-max',type=int,default=8)
    p.add_argument('--radii',type=float,nargs='+',default=[0.],help='include positive catalog-unit radii to calibrate smoothing; zero tests exact scan tasks')
    p.add_argument('--kernel3',type=int,choices=(0,1,2),default=0)
    p.add_argument('--reference-kernel3',type=int,choices=(0,1,2),default=0)
    p.add_argument('--pivot-block',type=int,default=0)
    p.add_argument('--include-pair-cells',action='store_true',help='also measure the existing exact pair-cell backend when applicable')
    for name,default,typ in [('rp-max',160,float),('rt-max',160,float),('r3-max',160,float),('rp-bins',50,int),('rt-bins',50,int),('r3-bins',4,int),('theta-bins',4,int),('mu-bins',4,int)]:p.add_argument('--'+name,type=typ,default=default)
    p.add_argument('--relative-floor',type=float,default=1e-6);p.add_argument('--max-relative-error',type=float,default=.05)
    p.add_argument('--small-signal-atol',type=float,default=5e-8)
    p.add_argument('--require-accepted-smoothing',action='store_true',help='return 2 unless every requested method has an accepted positive-radius case that actually merges pixels')
    args=p.parse_args(argv)
    if not 1<=args.scan_level<=20 or not 1<=args.pivot_max<=1024 or not 0<=args.pivot_block<=4096:p.error('invalid frontier controls')
    if min(args.threads,args.repeats,args.pixel_stride,args.rp_bins,args.rt_bins,args.r3_bins,args.theta_bins,args.mu_bins)<1 or args.warmups<0:p.error('invalid count')
    if any(not np.isfinite(v) or v<=0 for v in (args.rp_max,args.rt_max,args.r3_max)):p.error('maxima must be positive and finite')
    if any(not np.isfinite(v) or v<0 for v in [*args.radii,args.relative_floor,args.max_relative_error,args.small_signal_atol]):p.error('radii and accuracy thresholds must be finite and nonnegative')
    if args.max_forests is not None and args.max_forests<1:p.error('max-forests must be positive')
    if args.fits:
        from lya_fits import read_fits
        catalog=read_fits(args.fits,max_forests=args.max_forests,pixel_stride=args.pixel_stride)
    elif args.max_forests is not None or args.pixel_stride!=1:p.error('selection requires FITS input')
    elif args.catalog:catalog=read_npz(args.catalog)
    elif args.ascii:catalog=read_ascii(args.ascii)
    else:catalog=synthetic_catalog()
    catalog=catalog.normalized();out=args.outdir.resolve();out.mkdir(parents=True,exist_ok=False)
    save_catalog(out/'catalog.npz',catalog);cat=out/'catalog.txt'
    with cat.open('w') as stream:
        for xyz,delta,weight,fid in zip(catalog.positions,catalog.delta,catalog.weights,catalog.forest_ids):
            stream.write(' '.join(format(v,'.17g') for v in (*xyz,delta,weight))+f' {int(fid)}\n')
    binary=args.cballs.resolve();baseline=(args.baseline_cballs or binary).resolve()
    methods=list(dict.fromkeys(args.methods))
    report=dict(platform=platform.platform(),methods=methods,pixels=catalog.nbody,forests=len(np.unique(catalog.forest_ids)),
        openmp_environment={key:os.environ.get(key) for key in ('OMP_NUM_THREADS','OMP_DYNAMIC','OMP_WAIT_POLICY','OMP_PROC_BIND','OMP_PLACES','GOMP_CPU_AFFINITY','KMP_AFFINITY')},
        catalog_sha256=hashlib.sha256(cat.read_bytes()).hexdigest(),selection=catalog.metadata,settings=vars(args),
        binaries={str(b):hashlib.sha256(b.read_bytes()).hexdigest() for b in (binary,baseline)},
        timing_scope='whole fresh native process: input, tree/frontier, search and output; Python preprocessing excluded',
        cpu_scope='total process user+system CPU, never divided by threads',rss_scope='per-child wait4 high-water RSS; unavailable on platforms without wait4',
        runs=[],comparisons=[],accepted_smoothing={m:[] for m in methods})
    def save():(out/'summary.json').write_text(json.dumps(report,indent=2,default=str,allow_nan=False)+'\n')
    exact_failure=False
    for method in methods:
        cases=[('reference',0.,0,0)]
        if baseline!=binary:cases.append(('catalog-order',0.,0,0))
        for radius in dict.fromkeys([0.,*args.radii]):cases.append((f'scan{args.scan_level}-r{radius:g}',radius,args.scan_level,0))
        if args.include_pair_cells and '2pcf' in method:cases.append(('exact-pair-cells',0.,0,1))
        results={label:dict(measurements=[],products={},counts={}) for label,_,_,_ in cases}
        for repeat in range(-args.warmups,args.repeats):
            for label,radius,level,pair_kernel in cases:
                folder=out/f'{method}-{label}-{repeat}';folder.mkdir()
                params=dict(search=method,infile=cat,infileformat='lya-ascii',iCatalogs=1,rootDir=folder,
                    numberThreads=args.threads,usePeriodic='false',useLogHist='false',rangeN=args.r3_max,rminHist=.1,sizeHistN=4,
                    lya2RpMax=args.rp_max,lya2RtMax=args.rt_max,lya2RpBins=args.rp_bins,lya2RtBins=args.rt_bins,
                    lya3RMax=args.r3_max,lya3RBins=args.r3_bins,lya3ThetaBins=args.theta_bins,lya3MuBins=args.mu_bins,
                    lya3Kernel=args.reference_kernel3 if label=='reference' else args.kernel3,lya3PivotBlock=args.pivot_block,
                    verbose=2,verbose_log=1,options='no-smooth-pivot,lya-output-empty-bins')
                if label!='reference':params.update(lyaScanLevel=level,lyaPivotRadius=radius,lyaPivotMax=args.pivot_max,lya2Kernel=pair_kernel)
                command=[str(baseline if label=='reference' else binary),*(f'{k}={v}' for k,v in params.items())]
                print(f'{method}: {label}, repeat {repeat}',flush=True)
                measurement=run_process(command,folder/'process.log');measurement.pop('native_search_cpu_seconds',None)
                text=(folder/'process.log').read_text();line=next((x for x in text.splitlines() if x.startswith('Ly-alpha pivot frontier:')),'')
                measurement['frontier']={k:float(v) for k,v in re.findall(r'(\w+)=([\d.eE+-]+)',line)}
                saved=results[label]
                for order,name in ((2,'histXi2pcf_lya.txt'),(3,'histZetaM_lya5d.txt')):
                    if f'{order}pcf' not in method:continue
                    path=folder/name;table=np.loadtxt(path,ndmin=2)
                    if order in saved['products']:np.testing.assert_array_equal(saved['products'][order],table)
                    saved['products'][order]=table
                    saved['counts'][order]=int(re.search(r'# distinct-forest .*: (\d+)',path.read_text())[1])
                if repeat>=0:saved['measurements'].append(measurement)
                report['runs'].append(dict(method=method,case=label,repeat=repeat,warmup=repeat<0,command=command,**measurement));save()
        baseline_wall=float(np.median([x['wall_seconds'] for x in results['reference']['measurements']]))
        for label,radius,level,pair_kernel in cases:
            saved=results[label];measurements=saved['measurements'];products={}
            for order,table in saved['products'].items():
                reference=results['reference']['products'][order]
                c=comparison(reference,table,args.relative_floor,args.max_relative_error,args.small_signal_atol,radius==0)
                c['count']=saved['counts'][order];c['reference_count']=results['reference']['counts'][order]
                c['same_count']=c['count']==c['reference_count']
                if radius==0:c['accepted'] &= c['same_count']
                products[str(order)]=c
                np.savetxt(out/f'{method}-{label}-{order}pcf-per-bin.txt',np.column_stack((reference[:,:-3],reference[:,-3:],table[:,-3:],np.abs(table[:,-3]-reference[:,-3]))),
                           header='reference axes; reference correlation N D; candidate correlation N D; absolute correlation error')
            walls=[x['wall_seconds'] for x in measurements];cpu=[x['process_cpu_seconds'] for x in measurements if x['process_cpu_seconds'] is not None];rss=[x['peak_rss_bytes'] for x in measurements if x['peak_rss_bytes'] is not None]
            accepted=all(x['accepted'] for x in products.values())
            frontier=measurements[0]['frontier']
            grouped=bool(radius>0 and frontier.get('representatives',0)<frontier.get('active',0))
            row=dict(method=method,case=label,radius=radius,scan_level=level,pair_kernel=pair_kernel,products=products,accepted=accepted,
                frontier=frontier,geometry_actually_grouped=grouped,
                median_wall_seconds=float(np.median(walls)),min_wall_seconds=min(walls),max_wall_seconds=max(walls),
                median_cpu_seconds=float(np.median(cpu)) if cpu else None,peak_rss_bytes=max(rss) if rss else None,wall_speedup=baseline_wall/float(np.median(walls)))
            report['comparisons'].append(row)
            if radius==0 and not accepted:exact_failure=True
            if grouped and accepted:report['accepted_smoothing'][method].append(label)
            save();print(json.dumps(row),flush=True)
    report['exact_validation_passed']=not exact_failure;save();plot_summary(report,out)
    if exact_failure:return 1
    if args.require_accepted_smoothing and not all(report['accepted_smoothing'].values()):return 2
    return 0

if __name__=='__main__':raise SystemExit(main())
