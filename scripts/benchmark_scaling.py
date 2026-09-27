#!/usr/bin/env python3
"""Fixed-selection Ly-alpha strong/size scaling with matched raw products.

Use --catalog on the existing XYZ FITS or ForestCatalog NPZ. Each configuration
runs in a fresh process; warmups and repeats force new computations in that
process. Outputs include exact same-input comparisons, errors, wall/CPU/RSS,
compiler/build identity and native cache measurements. No runtime is accepted
as a speed claim unless its scientific comparison passes.
"""
import argparse, hashlib, json, os, platform, statistics, subprocess, sys, time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from workload_acceptance import compare,peak


def load(path):
    import numpy as np
    if path.suffix=='.npz':
        with np.load(path,allow_pickle=False) as d:
            return {k:np.asarray(d[k]) for k in ('positions','delta','weights','forest_ids')}
    from astropy.io import fits
    with fits.open(path,memmap=True) as f:
        d=f[1].data;names={n.lower():n for n in d.names}
        return dict(positions=np.column_stack([d[names[n]] for n in ('x','y','z')]).astype('f8'),delta=np.asarray(d[names['delta']],dtype='f8'),weights=np.asarray(d[names.get('weight',names.get('weights'))],dtype='f8'),forest_ids=np.asarray(d[names.get('los_id',names.get('forest_ids'))],dtype='i8'))


def select(data,forests,stride):
    import numpy as np
    ids=data['forest_ids'];unique=np.unique(ids)
    if forests>len(unique):raise ValueError('requested more forests than catalog contains')
    # Stable input row order within each sorted forest ID; nested forest prefixes.
    order=np.argsort(ids,kind='stable');ordered=ids[order];edges=np.r_[0,np.flatnonzero(ordered[1:]!=ordered[:-1])+1,len(order)]
    rows=np.concatenate([order[edges[k]:edges[k+1]:stride] for k in range(forests)])
    return {k:v[rows] for k,v in data.items()},rows


def worker(a):
    import numpy as np
    import cyballs
    case=json.loads(a.worker);data=load(a.catalog);a.output.mkdir(parents=True,exist_ok=False)
    values=[];samples=[];m=cyballs.cballs()
    parameters=dict(searchMethod=case['engine'],numberThreads=case['threads'],verbose=0,verbose_log=0,
        usePeriodic=False,useLogHist=False,rangeN=case['radius'],rminHist=.001,sizeHistN=6,
        lya2RpMax=case['radius'],lya2RtMax=case['radius'],lya2RpBins=8,lya2RtBins=8,
        lya3RMax=case['radius'],lya3RBins=4,lya3ThetaBins=4,lya3MuBins=8,
        options='no-out-Hist,no-smooth-pivot',lya2Kernel=0,lya3Kernel=0)
    if case['mode']!='reference':parameters.update(lya2Kernel=1 if '2pcf' in case['engine'] else 0,lya3Kernel=case.get('triple_kernel',0) if '3pcf' in case['engine'] else 0)
    if case['mode']=='approx':parameters.update(lya2RpSlop=case['slop'] if '2pcf' in case['engine'] else 0.,lya2RtSlop=case['slop'] if '2pcf' in case['engine'] else 0.,lya3MuSlop=case['slop'] if '3pcf' in case['engine'] else 0.,lya3RadialSlop=case['slop'] if '3pcf' in case['engine'] else 0.,lya3PolarSlop=case['slop'] if '3pcf' in case['engine'] else 0.)
    forecast=cyballs.resource_plan(case['engine'],len(data['positions']),case['threads'],parameters)
    if forecast['known_exceeds_budget']:raise MemoryError('known plan exceeds configured budget')
    try:
        m.set_forest_catalog(data['positions'],data['delta'],data['weights'],data['forest_ids'])
        for index in range(case['warmups']+case['repeats']):
            params=parameters|dict(rootDir=str(a.output/f'run-{index}'))
            m.set(params) # Changed root forces recomputation; Run cache is not timed.
            cache_before=m.getCacheInfo();start=time.perf_counter();cpu=time.process_time();m.Run()
            sample=dict(iteration=index,warmup=index<case['warmups'],wall_seconds=time.perf_counter()-start,process_cpu_seconds=time.process_time()-cpu,
                native_timings=m.getTimings(),peak_rss_bytes=peak(),cache_before=cache_before,cache_after=m.getCacheInfo())
            result=m.getForestResults();arrays=result['arrays']
            if values:
                for k in arrays:np.testing.assert_allclose(arrays[k],values[0][k],rtol=3e-11,atol=1e-10)
            if not values:values=[arrays]
            if not sample['warmup']:samples.append(sample)
        np.savez_compressed(a.output/'products.npz',**values[0])
        record=dict(case=case,samples=samples,settings=parameters,pixels=len(data['positions']),forests=len(np.unique(data['forest_ids'])),
            fixture_sha256=hashlib.sha256(a.catalog.read_bytes()).hexdigest(),build=cyballs.build_info(),metadata=m.getRunMetadata(),forecast=forecast,
            cache_policy='fresh process for each configuration; first run cold; subsequent runs forced by changed rootDir on same object; forest trees currently rebuild per call',
            rss_semantics='per-process cumulative peak including interpreter/input/previous repetitions; no summed worker RSS')
        (a.output/'sample.json').write_text(json.dumps(record,indent=2)+'\n')
    finally:m.struct_cleanup()


def run(a):
    import numpy as np
    a.output.mkdir(parents=True,exist_ok=False);data=load(a.catalog)
    report=dict(schema_version=1,status='RUNNING',source_catalog=str(a.catalog),source_sha256=hashlib.sha256(a.catalog.read_bytes()).hexdigest(),
        host=dict(platform=platform.platform(),cpu_count=os.cpu_count(),python=sys.version),selection=dict(forest_order='ascending unique integer forest ID',pixel_order='original within-forest row order',pixel_stride=a.pixel_stride),
        methodology='one worker process at a time; full same-input exact reference; approximation evaluated separately; no noise filtering or tolerance tuning',samples=[],summary=[])
    def save():(a.output/'scaling.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    save()
    try:
        for count in a.forest_counts:
            selected,rows=select(data,count,a.pixel_stride);fixture=a.output/f'forests-{count}.npz';np.savez_compressed(fixture,**selected)
            np.save(a.output/f'forests-{count}-source-rows.npy',rows)
            for prefix in a.prefixes:
                engine=f'{prefix}-{a.statistic}-omp';reference=None;baseline_wall=None
                for mode in a.modes:
                    for threads in a.threads:
                        case=dict(engine=engine,mode=mode,threads=threads,radius=a.radius,slop=a.slop,warmups=a.warmups,repeats=a.repeats,triple_kernel=a.triple_kernel)
                        label=f'{count}-{engine}-{mode}-t{threads}';directory=a.output/label
                        with (a.output/(label+'.log')).open('w') as f:
                            result=subprocess.run([sys.executable,str(Path(__file__).resolve()),'--catalog',str(fixture),'--output',str(directory),'--worker',json.dumps(case)],cwd=ROOT,stdout=f,stderr=subprocess.STDOUT,timeout=a.timeout)
                        if result.returncode:raise RuntimeError('worker failed: '+label)
                        row=json.loads((directory/'sample.json').read_text());values=dict(np.load(directory/'products.npz'))
                        if reference is None:reference=values
                        metrics={key:compare(target,values[key],.05 if mode=='approx' else 3e-11,1e-10) for key,target in reference.items()}
                        for name in ('pair','triple'):
                            nk,dk=name+'_numerator',name+'_denominator'
                            if nk in reference:
                                nr,dr=reference[nk],reference[dk];nc,dc=values[nk],values[dk]
                                metrics[name+'_ratio']=compare(np.divide(nr,dr,out=np.zeros_like(dr),where=dr>0),np.divide(nc,dc,out=np.zeros_like(dc),where=dc>0),.05 if mode=='approx' else 3e-11,1e-12)
                                metrics[name+'_occupancy']=dict(passed=bool(np.array_equal(dr>0,dc>0)))
                        accepted=all(v['passed'] for v in metrics.values());measurements=row['samples'];wall=statistics.median(x['wall_seconds'] for x in measurements)
                        if baseline_wall is None:baseline_wall=wall
                        summary=dict(label=label,engine=engine,mode=mode,threads=threads,forests=count,pixels=row['pixels'],
                            status='QUALIFIED' if accepted else 'REJECTED',wall_median_seconds=wall,wall_min_seconds=min(x['wall_seconds'] for x in measurements),wall_max_seconds=max(x['wall_seconds'] for x in measurements),
                            cpu_median_seconds=statistics.median(x['process_cpu_seconds'] for x in measurements),peak_rss_bytes=max(x['peak_rss_bytes'] for x in measurements),
                            speedup_vs_reference_one_thread=baseline_wall/wall if accepted else None,metrics=metrics)
                        report['samples'].append(row|dict(directory=label));report['summary'].append(summary);save();print(label,summary['status'],wall,flush=True)
        report['status']='PASS' if all(r['status']=='QUALIFIED' for r in report['summary'] if r['mode']!='approx') else 'FAIL'
    except BaseException as exc:report.update(status='INTERRUPTED' if isinstance(exc,KeyboardInterrupt) else 'FAIL',error=str(exc));raise
    finally:save()
    plot(report,a.output)
    if report['status']!='PASS':raise SystemExit(1)


def plot(report,output):
    os.environ.setdefault('MPLCONFIGDIR',str(output/'matplotlib'))
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axs=plt.subplots(1,3,figsize=(14,4))
    labels=sorted({(r['engine'],r['mode'],r['threads']) for r in report['summary']})
    for engine,mode,threads in labels:
        rows=sorted([r for r in report['summary'] if (r['engine'],r['mode'],r['threads'])==(engine,mode,threads) and r['status']=='QUALIFIED'],key=lambda r:r['pixels'])
        if not rows:continue
        label=f'{engine.replace("-omp","")} {mode} t{threads}'
        for ax,key in zip(axs,['wall_median_seconds','cpu_median_seconds','peak_rss_bytes']):ax.loglog([r['pixels'] for r in rows],[r[key] for r in rows],'o-',label=label)
    for ax,title in zip(axs,['Wall seconds','Process CPU seconds','Peak RSS bytes']):ax.set(xlabel='Retained pixels',ylabel=title);ax.grid(True,alpha=.2)
    axs[-1].legend(fontsize=5);fig.tight_layout();fig.savefig(output/'size-scaling.png',dpi=180);plt.close(fig)
    fig,axs=plt.subplots(1,2,figsize=(10,4))
    for count in sorted({r['forests'] for r in report['summary']}):
        for engine,mode,_ in sorted({(r['engine'],r['mode'],0) for r in report['summary']}):
            rows=sorted([r for r in report['summary'] if r['forests']==count and r['engine']==engine and r['mode']==mode and r['status']=='QUALIFIED'],key=lambda r:r['threads'])
            if not rows:continue
            label=f'{engine} {mode} {rows[0]["pixels"]} pixels'
            axs[0].plot([r['threads'] for r in rows],[r['wall_median_seconds'] for r in rows],'o-',label=label)
            axs[1].plot([r['threads'] for r in rows],[rows[0]['wall_median_seconds']/r['wall_median_seconds'] for r in rows],'o-',label=label)
    axs[0].set(xlabel='Threads',ylabel='Wall seconds');axs[1].set(xlabel='Threads',ylabel='Speedup at fixed input');axs[1].legend(fontsize=5)
    for ax in axs:ax.grid(True,alpha=.2)
    fig.tight_layout();fig.savefig(output/'strong-scaling.png',dpi=180);plt.close(fig)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--catalog',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--forest-counts',type=int,nargs='+',default=[16,64,256]);p.add_argument('--pixel-stride',type=int,default=20)
    p.add_argument('--threads',type=int,nargs='+',default=[1,2,4]);p.add_argument('--prefixes',nargs='+',choices=['lya','lya-los-tree'],default=['lya','lya-los-tree'])
    p.add_argument('--statistic',choices=['2pcf','3pcf','2pcf-3pcf'],default='2pcf-3pcf');p.add_argument('--modes',nargs='+',choices=['reference','optimized','approx'],default=['reference','optimized'])
    p.add_argument('--triple-kernel',type=int,choices=range(5),default=0,help='candidate exact 3PCF kernel; 0 keeps segment traversal, 2 tiles direct pairs, 3 uses persistent nodes, 4 adds pivot cells');p.add_argument('--radius',type=float,default=80.);p.add_argument('--slop',type=float,default=.01);p.add_argument('--repeats',type=int,default=3);p.add_argument('--warmups',type=int,default=1);p.add_argument('--timeout',type=float,default=1800)
    p.add_argument('--worker',help=argparse.SUPPRESS);a=p.parse_args();a.catalog=a.catalog.resolve();a.output=a.output.resolve()
    if any(x<1 for x in a.threads+a.forest_counts+[a.pixel_stride,a.repeats]) or a.warmups<0:p.error('counts/stride/repeats/threads must be positive; warmups nonnegative')
    if not a.worker and 'approx' in a.modes and '3pcf' in a.statistic and a.triple_kernel not in (0,3,4):p.error('3PCF slop requires kernel 0, 3 or 4')
    if not a.worker and (a.modes[0]!='reference' or a.threads[0]!=1):p.error('reference mode and one thread must be first')
    (worker if a.worker else run)(a)
