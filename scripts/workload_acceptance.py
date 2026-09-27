#!/usr/bin/env python3
"""Held-out qualification of explicit settings; no tolerance is fitted to results.

PASS means every required exact/control case passed. Approximate candidates
are individually QUALIFIED or REJECTED for this workload, never silently
promoted. Retain raw/corrected products and weak-signal errors in either case.
"""
import argparse, hashlib, json, os, platform, resource, subprocess, sys, time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))


def peak():
    n=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if sys.platform=='darwin':return int(n)
    if sys.platform.startswith('linux'):return int(n*1024)
    raise RuntimeError('unknown RSS units')


def fixture(kind, seed=7260926):
    import numpy as np
    rng=np.random.default_rng(seed)
    if kind=='masked_weak_sky':
        p=rng.normal(size=(500,3));p/=np.linalg.norm(p,axis=1)[:,None];p=p[(p[:,2]>-.65)&(p[:,0]<.85)][:256]
    elif kind=='clustered_signed_sky':
        p=np.repeat(rng.normal(size=(64,3)),4,axis=0);p/=np.linalg.norm(p,axis=1)[:,None]
        p+=rng.normal(0,3e-5,p.shape);p/=np.linalg.norm(p,axis=1)[:,None]
    else:
        # Different radial lengths and directions; finite angular width + bent
        # rays exercise geometry bounds beyond a rectilinear training fixture.
        nf,npx=12,16
        directions=np.column_stack((np.ones(nf),rng.uniform(-.18,.18,(nf,2))))
        directions/=np.linalg.norm(directions,axis=1)[:,None]
        chi=np.linspace(95,145,npx)[None,:]+rng.uniform(-10,10,(nf,1))
        p=(directions[:,None,:]*chi[:,:,None]).reshape(-1,3)
        if kind=='bent_signed_forests':p[::7,1]+=rng.uniform(-1,1,len(p[::7]))
    w=rng.uniform(.25,2,len(p));field=rng.normal(0,.03,len(p));g=rng.normal(0,.001,len(p))+1j*rng.normal(0,.001,len(p))
    return p,field,w,g


def execute(case,root):
    import numpy as np
    import cyballs
    p,k,w,g=fixture(case['fixture']);m=cyballs.cballs();engine=case['engine'];family=case['family']
    root.mkdir(parents=True,exist_ok=True)
    fixture_arrays=dict(positions=p,kappa=k,weights=w,gamma=g)
    if family=='forest':fixture_arrays['forest_ids']=np.repeat(np.arange(12,dtype=np.int64),16)
    np.savez(root/'fixture.npz',**fixture_arrays)
    options=['no-out-Hist','no-smooth-pivot','weights-norm','KKKCorrelation','compute-HistN','no-normalize-HistZeta']
    params=dict(searchMethod=engine,rootDir=str(root/'native'),verbose=0,verbose_log=0,numberThreads=2,usePeriodic=False,
                sizeHistN=4,mChebyshev=2,rangeN=1.8,rminHist=.03,useLogHist=True,theta=0)
    if family!='forest':
        options+=['edge-corrections']
        if case['mode']=='exact':options+=['no-one-ball','no-two-balls']
        else:params['theta']=case['theta']
        if case['mode']=='smooth':
            options.remove('no-smooth-pivot');options+=['smooth-pivot'];params['rsmooth']='1' if family=='shear' else '.0003'
    else:
        params.update(useLogHist=False,rangeN=80.,rminHist=.1,lya2RpMax=80.,lya2RtMax=80.,lya2RpBins=5,lya2RtBins=6,
            lya3RMax=80.,lya3RBins=4,lya3ThetaBins=4,lya3MuBins=6,lya2Kernel=case.get('pair_kernel',0),lya3Kernel=case.get('triple_kernel',1),
            lya2RpSlop=case.get('slop',0),lya2RtSlop=case.get('slop',0),lya3MuSlop=case.get('slop',0),
            lya3RadialSlop=case.get('slop',0),lya3PolarSlop=case.get('slop',0))
    params['options']=','.join(options);m.set(params)
    if family=='forest':m.set_forest_catalog(p,k,w,fixture_arrays['forest_ids'])
    elif family=='shear':m.set_catalog(p,gamma1=g.real,gamma2=g.imag,weights=w)
    else:m.set_catalog(p,kappa=k,weights=w)
    start=time.perf_counter();cpu=time.process_time()
    try:
        m.Run();wall=time.perf_counter()-start;cpu=time.process_time()-cpu
        if family=='forest':values=m.getForestResults()['arrays']
        elif family=='shear':values=dict(pair=np.array([m.getShearXiPlus(),m.getShearXiMinus()]),raw=m.getShearUpsilonMultipoles(),corrected=m.getShearGammaMultipoles())
        else:values=dict(pair=m.getHistXi2pcf(),raw=np.array([m.getHistZetaMsincos(i,1)+m.getHistZetaMsincos(i,2)+1j*(m.getHistZetaMsincos(i,3)-m.getHistZetaMsincos(i,4)) for i in range(1,4)]),corrected=np.array([m.getHistZetaM_EE_complex(i) for i in range(1,4)]))
        if family=='forest':
            for prefix in ('pair','triple'):
                den=values[prefix+'_denominator'];values[prefix+'_ratio']=np.divide(values[prefix+'_numerator'],den,out=np.zeros_like(den),where=den>0)
        np.savez_compressed(root/'products.npz',**values)
        record=dict(case=case,settings=params,build=cyballs.build_info(),metadata=m.getRunMetadata(),
            wall_seconds=wall,process_cpu_seconds=cpu,peak_rss_bytes=peak(),cache_policy='fresh process; cold tree',
            fixture_sha256=hashlib.sha256((root/'fixture.npz').read_bytes()).hexdigest())
        (root/'sample.json').write_text(json.dumps(record,indent=2)+'\n')
    finally:m.struct_cleanup()


def compare(target,value,rtol=.05,atol=1e-12):
    import numpy as np
    finite=np.isfinite(target);same_finite=bool(np.array_equal(finite,np.isfinite(value)))
    ref=target[finite];candidate=value[finite];error=float(np.linalg.norm(candidate-ref));norm=float(np.linalg.norm(ref));delta=np.abs(candidate-ref)
    floor=max(atol,1e-6*float(np.max(np.abs(ref),initial=0)))
    eligible=np.abs(ref)>floor
    record = dict(absolute_l2=error,reference_l2=norm,relative_l2=error/max(norm,atol),rtol=rtol,atol=atol,
        same_finite_mask=same_finite,reference_nonfinite=int((~finite).sum()),
        signal_floor=floor,weak_bins=int((~eligible).sum()),weak_max_absolute=float(np.max(delta[~eligible],initial=0)),
        max_bin_relative=float(np.max(delta[eligible]/np.abs(ref[eligible]),initial=0)),
        passed=bool(same_finite and ref.size>0 and np.isfinite(error) and error<=rtol*norm+atol))
    return {k:(None if isinstance(v,float) and not np.isfinite(v) else v) for k,v in record.items()}


def run(output):
    import numpy as np
    import cyballs
    output.mkdir(parents=True,exist_ok=False)
    report=dict(schema_version=1,status='RUNNING',host=platform.platform(),criteria='per observable L2 <= 0.05*reference L2 + 1e-12; exact controls use 3e-11 + 1e-12; binwise errors also retained',
                seeds=[7260926],qualification_scope='held-out geometry/weak signed fields only; no survey-wide guarantee; corrected finite mask must match',samples=[],qualifications=[])
    families=[('scalar',[f'{t}-2balls-omp' for t in ('octree','kdtree','balltree')]),('shear',[f'{t}-shear-sphere-2balls-omp' for t in ('octree','kdtree','balltree')]),('forest',['lya-2pcf-3pcf-omp','lya-los-tree-2pcf-3pcf-omp'])]
    def save():(output/'acceptance.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    save();exact_fail=False
    try:
        for family,engines in families:
            fixtures=('masked_weak_sky','clustered_signed_sky') if family!='forest' else ('straight_signed_forests','bent_signed_forests')
            for engine in engines:
                if cyballs.search_method_id(engine)<0:continue
                for kind in fixtures:
                    modes=[dict(mode='exact')]
                    if family=='forest':modes += [dict(mode='aggregate_exact',pair_kernel=1,triple_kernel=4),dict(mode='approx',pair_kernel=1,triple_kernel=4,slop=.01)]
                    else:
                        modes += [dict(mode='approx',theta=.01)]
                        if engine!='octree-2balls-omp':modes += [dict(mode='smooth',theta=.01)]
                    exact=None;digest=None
                    for variant in modes:
                        case=dict(engine=engine,family=family,fixture=kind,**variant);label=f'{engine}-{kind}-{variant["mode"]}';directory=output/label
                        with (output/(label+'.log')).open('w') as log:
                            completed=subprocess.run([sys.executable,str(Path(__file__).resolve()),'--worker',json.dumps(case),'--output',str(directory)],cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,timeout=300)
                        if completed.returncode:raise RuntimeError('worker failed: '+label)
                        row=json.loads((directory/'sample.json').read_text());values=dict(np.load(directory/'products.npz'))
                        if exact is None:exact=values;digest=row['fixture_sha256']
                        assert digest==row['fixture_sha256'];assert set(values)==set(exact)
                        metrics={key:compare(target,values[key],3e-11 if variant['mode'] in ('exact','aggregate_exact') else .05) for key,target in exact.items()}
                        passed=all(v['passed'] for v in metrics.values())
                        if variant['mode'] in ('exact','aggregate_exact') and not passed:exact_fail=True
                        status='QUALIFIED' if passed else 'REJECTED'
                        report['samples'].append(dict(directory=label,**row,metrics=metrics,qualification=status))
                        report['qualifications'].append(dict(case=case,status=status,observables=list(metrics)))
                        save();print(label,status,flush=True)
        report['status']='FAIL' if exact_fail else 'PASS'
    except BaseException as exc:report.update(status='INTERRUPTED' if isinstance(exc,KeyboardInterrupt) else 'FAIL',error=str(exc));raise
    finally:save()
    if exact_fail:raise SystemExit(1)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--worker',help=argparse.SUPPRESS);a=p.parse_args();a.output=a.output.resolve()
    (execute(json.loads(a.worker),a.output) if a.worker else run(a.output))
