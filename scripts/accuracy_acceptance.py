#!/usr/bin/env python3
"""Retain exact-subsample evidence for a deliberately bounded accuracy envelope.

These are regression acceptance tolerances, not guarantees for arbitrary data.
No tolerance is inferred from the measured errors during this run.
"""
import argparse
import json
from pathlib import Path
import tempfile
from contextlib import nullcontext
import numpy as np
from cyballs import cballs, build_info, search_method_id

# Max complex Frobenius error, measured separately for 2PCF and 3PCF.
# The octree uses a much smaller opening angle than the binary trees.
SCALAR = [
    dict(engine='kdtree-2balls-omp',theta=.5,rtol=.02,smoothing=True),
    dict(engine='balltree-2balls-omp',theta=.2,rtol=.02,smoothing=True),
    dict(engine='octree-2balls-omp',theta=.025,rtol=.02,smoothing=False),
    dict(engine='octree-sincos-omp',theta=.025,rtol=.02,smoothing=True,diagnostic=True),
    dict(engine='octree-2balls-omp',theta=.025,rtol=.02,smoothing=True,compatibility=True,diagnostic=True),
]
SHEAR = [dict(engine=f'{tree}-shear-sphere-2balls-omp',theta=.05,rtol=.02,smoothing=True)
         for tree in ('octree','kdtree','balltree')]


def fixture(kind):
    rng=np.random.default_rng(20260921 if kind=='signed_sky' else 20260922)
    if kind=='signed_sky':
        p=rng.normal(size=(512,3)); p/=np.linalg.norm(p,axis=1)[:,None]
        k=rng.normal(size=len(p))
    else:
        # Small groups exercise smoothing and accepted cells, while the field
        # changes slowly across a group. Several wide angular scales remain.
        centers=rng.normal(size=(64,3));centers/=np.linalg.norm(centers,axis=1)[:,None]
        p=np.repeat(centers,8,axis=0)+rng.normal(scale=1e-5,size=(512,3))
        p/=np.linalg.norm(p,axis=1)[:,None]
        k=.3+p[:,0]+.2*p[:,1]-.1*p[:,2]**2
    w=rng.uniform(.5,1.5,len(p))
    gamma=(.2+.1*p[:,0]+.03*p[:,2])+1j*(.1*p[:,1]-.04*p[:,2])
    return p,k,w,gamma


def run(case,data,mode,root,*,phase=None,threads=2):
    observe=phase or (lambda name: nullcontext())
    p,k,w,g=data;engine=case['engine'];shear='shear' in engine
    options=['no-out-Hist','compute-HistN','KKKCorrelation','weights-norm','no-normalize-HistZeta']
    if case.get('compatibility'):options+=['legacy-one-ball']
    if mode=='exact':
        options+=['no-one-ball']
        if not case.get('compatibility'):options+=['no-two-balls']
    options+=['smooth-pivot' if mode=='smooth' else 'no-smooth-pivot']
    if phase:options+=['dual-node-profile']
    with observe('construct'):
        m=cballs()
    try:
        with observe('configure_and_register'):
            m.set(searchMethod=engine,rootDir=str(root),sizeHistN=6,mChebyshev=5,
                  rangeN=1.8,rminHist=.02,useLogHist=True,usePeriodic=False,
                  numberThreads=threads,theta=case['theta'] if mode!='exact' else 0.,
                  rsmooth='1' if shear else '0.0003',nsmooth=8,verbose=0,verbose_log=0,options=','.join(options))
            if shear:m.set_catalog(p,gamma1=g.real,gamma2=g.imag,weights=w)
            else:m.set_catalog(p,kappa=k,weights=w)
        if phase:
            with observe('initialize'):
                m.Run(level=['Initial'])
        with observe('main_loop'):
            m.Run()
        with observe('extract_products'):
            if shear:
                products={'pair':np.asarray([m.getShearXiPlus(),m.getShearXiMinus()]),
                          'triple':np.asarray(m.getShearUpsilonMultipoles())}
            else:
                products={'pair':np.asarray(m.getHistXi2pcf()),
                          'triple':np.asarray([m.getHistZetaMsincos(i,1)+m.getHistZetaMsincos(i,2)
                             +1j*(m.getHistZetaMsincos(i,3)-m.getHistZetaMsincos(i,4)) for i in range(1,7)])}
        return products,m.getRunMetadata(),m.getTimings()
    finally:
        with observe('cleanup'):
            m.struct_cleanup()


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',required=True,type=Path)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    report=dict(schema_version=1,status='RUNNING',build=build_info(),cases=[],
                metric='complex Frobenius norm(error)/norm(exact), separately for pair and triple products',
                absolute_l2_floor=1e-10,
                scope='512-point signed full-sky and clustered smooth-field weighted exact subsamples; 6 log bins .02..1.8; modes 0..5; 2 threads; no edge correction; smoothing radius .0003 Cartesian units for scalar, 1 arcmin for shear, only on clustered fixture; raw shear Upsilon multipoles',
                exclusions=['arbitrary larger theta or rsmooth','survey window conditioning','experimental shear-pivot-reuse','general accuracy guarantee outside these fixtures','octree-sincos and octree-GGG compatibility: diagnostic only; not accepted for clustered workload'])
    products={}; exact_families={}; report['exact_cross_engine']=[]
    try:
        for kind in ('signed_sky','clustered_smooth'):
            data=fixture(kind);np.savez(args.output/(kind+'-fixture.npz'),positions=data[0],kappa=data[1],weights=data[2],gamma=data[3])
            for index,case in enumerate(SCALAR+SHEAR):
                if search_method_id(case['engine'])<0:continue
                tag=f'{kind}-{index}-{case["engine"]}'
                with tempfile.TemporaryDirectory(dir=args.output) as tmp:
                    exact,metadata,timing=run(case,data,'exact',Path(tmp)/'exact')
                    for name,value in exact.items():products[tag+'-exact-'+name]=value
                    if not case.get('diagnostic'):
                        family=(kind,'shear' if 'shear' in case['engine'] else 'scalar')
                        if family not in exact_families: exact_families[family]=exact
                        for name,target in exact_families[family].items():
                            error=float(np.linalg.norm(exact[name]-target));norm=float(np.linalg.norm(target))
                            report['exact_cross_engine'].append(dict(engine=case['engine'],fixture=kind,product=name,
                                relative_l2=error/max(norm,1e-10),passed=bool(error<=1e-11*norm+1e-10)))
                    modes=['approx']+(['smooth'] if case['smoothing'] and kind=='clustered_smooth' else [])
                    for mode in modes:
                        value,metadata,timing=run(case,data,mode,Path(tmp)/mode)
                        record=dict(case,fixture=kind,mode=mode,metrics={},metadata=metadata,timings=timing)
                        for name,target in exact.items():
                            error=float(np.linalg.norm(value[name]-target));norm=float(np.linalg.norm(target))
                            record['metrics'][name]=dict(relative_l2=error/max(norm,1e-10),absolute_l2=error,reference_l2=norm,
                                passed=bool(np.isfinite(error) and error<=case['rtol']*norm+1e-10))
                            products[tag+'-'+mode+'-'+name]=value[name]
                        passed=all(v['passed'] for v in record['metrics'].values())
                        record['acceptance_status']=('DIAGNOSTIC_ONLY' if case.get('diagnostic') else ('PASS' if passed else 'FAIL'))
                        report['cases'].append(record)
                        print(tag,mode,{k:v['relative_l2'] for k,v in record['metrics'].items()},flush=True)
        report['status']='PASS' if report['cases'] and all(v['passed'] for v in report['exact_cross_engine']) and all(v['passed'] for c in report['cases'] if not c.get('diagnostic') for v in c['metrics'].values()) else 'FAIL'
    finally:
        np.savez_compressed(args.output/'products.npz',**products)
        # Strict JSON: nonfinite diagnostic errors are represented explicitly.
        def finite_json(value):
            if isinstance(value,float) and not np.isfinite(value):return None
            if isinstance(value,dict):return {k:finite_json(v) for k,v in value.items()}
            if isinstance(value,list):return [finite_json(v) for v in value]
            return value
        (args.output/'accuracy.json').write_text(json.dumps(finite_json(report),indent=2,allow_nan=False)+'\n')
    if report['status']!='PASS':raise SystemExit('accuracy acceptance failed; inspect retained metrics')

if __name__=='__main__':main()
