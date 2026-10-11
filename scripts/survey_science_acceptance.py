#!/usr/bin/env python3
"""Exact-reference campaign for fast defaults; rejection is a scientific result.

A passing execution receipt means all required comparisons were made and exact
controls passed. Each candidate independently remains QUALIFIED or REJECTED.
No default or tolerance is tuned to the outcome of this campaign.
"""
from acceptance_metadata import SURVEY_RTOL, SURVEY_ATOL, SURVEY_EXACT_RTOL, acceptance_policy
import argparse
import json
from pathlib import Path
import sys
import tempfile
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from cyballs import cballs, build_info, qualification_report, save_result_packet

ENGINES=[f'{t}-2balls-omp' for t in ('octree','kdtree','balltree')]+[
    f'{t}-shear-sphere-2balls-omp' for t in ('octree','kdtree','balltree')]+[
    'lya-2pcf-3pcf-omp','lya-los-tree-2pcf-3pcf-omp']
FIXTURES=('signed','clustered','masked','weak_signal')


def fixture(kind,forest=False):
    rng=np.random.default_rng(202610115)
    if forest:
        directions=np.column_stack((np.ones(12),rng.uniform(-.18,.18,(12,2))))
        directions/=np.linalg.norm(directions,axis=1)[:,None]
        chi=np.linspace(95,145,12)[None,:]+rng.uniform(-10,10,(12,1))
        p=(directions[:,None,:]*chi[:,:,None]).reshape(-1,3)
        if kind=='clustered':p[:,1:]*=.1
        ids=np.repeat(np.arange(12,dtype=np.int64),12)
    else:
        p=rng.normal(size=(192,3))
        if kind=='clustered':p=np.repeat(rng.normal(size=(24,3)),8,axis=0)+rng.normal(0,.002,(192,3))
        p/=np.linalg.norm(p,axis=1)[:,None];ids=None
    k=rng.normal(size=len(p));g=rng.normal(0,.03,len(p))+1j*rng.normal(0,.03,len(p))
    if kind=='weak_signal':k*=1e-5;g*=1e-3
    weights=rng.uniform(.25,2,len(p));mask=np.ones(len(p),dtype=bool)
    if kind=='masked':
        if forest:mask=(ids%3!=0)
        else:mask=(p[:,2]>-.4)&(p[:,0]<.75)
    return dict(positions=p,kappa=k,gamma=g,weights=weights,mask=mask,forest_ids=ids)


def execute(engine,data,root,exact,threads=1):
    forest=engine.startswith('lya-');shear='shear' in engine
    options=['no-out-Hist','no-smooth-pivot','weights-norm','KKKCorrelation','compute-HistN','no-normalize-HistZeta','read-mask']
    p=dict(searchMethod=engine,rootDir=str(root),verbose=0,verbose_log=0,numberThreads=threads,
        usePeriodic=False,useLogHist=True,sizeHistN=4,mChebyshev=3,rangeN=1.8,rminHist=.03)
    if not forest:
        options+=['edge-corrections']
        if exact:p['theta']=0.;options+=['no-one-ball','no-two-balls']
    else:
        p.update(useLogHist=False,rangeN=80.,rminHist=.1,lya2RpMax=80.,lya2RtMax=80.,lya2RpBins=4,lya2RtBins=4,
            lya3RMax=80.,lya3RBins=4,lya3ThetaBins=3,lya3MuBins=4)
        if exact:p.update(lya2Kernel=0,lya3Kernel=1,lya2RpSlop=0.,lya2RtSlop=0.,lya3MuSlop=0.,lya3RadialSlop=0.,lya3PolarSlop=0.)
    p['options']=','.join(options);m=cballs()
    try:
        m.set(p)
        if forest:
            # Forest loader selection is explicit before registration; preserve IDs.
            take=data['mask'];m.set_forest_catalog(data['positions'][take],data['kappa'][take],data['weights'][take],data['forest_ids'][take])
        elif shear:m.set_catalog(data['positions'],gamma1=data['gamma'].real,gamma2=data['gamma'].imag,weights=data['weights'],mask=data['mask'])
        else:m.set_catalog(data['positions'],kappa=data['kappa'],weights=data['weights'],mask=data['mask'])
        m.Run();return m.getResults()
    finally:m.struct_cleanup()


def selection(packet):
    """Predeclared interior radial bins 1 and 2; all angular modes retained."""
    result={};forest=packet['metadata']['engine'].startswith('lya-')
    for name,value in packet['arrays'].items():
        mask=np.ones(value.shape,dtype=bool)
        axes=(0,1) if forest and value.ndim>1 else tuple(range(max(0,value.ndim-2),value.ndim))
        for axis in axes:
            if value.shape[axis]!=4:continue
            shape=[1]*value.ndim;shape[axis]=4
            mask &= np.array([False,True,True,False]).reshape(shape)
        result[name]=mask
    return result


def run(output):
    output.mkdir(parents=True,exist_ok=False)
    report=dict(schema_version=1,acceptance_policy=acceptance_policy(),status='RUNNING',build=build_info(),cases=[],
        default_changes=False,thresholds=dict(rtol=SURVEY_RTOL,atol=SURVEY_ATOL),
        selection='interior radial bins 1 and 2 (zero based); every angular mode and shear component',
        scope='192-point sky / 144-pixel 12-forest fixtures; each returned raw, normalized and corrected observable; no transfer to another survey',
        timing_scope='single MainLoop sample per setting; diagnostic cost, not a scaling or speed ranking')
    try:
        for engine in ENGINES:
            for kind in FIXTURES:
                root=output/engine/kind;root.mkdir(parents=True)
                data=fixture(kind,engine.startswith('lya-'))
                np.savez_compressed(root/'catalog.npz',**{k:v for k,v in data.items() if v is not None})
                exact=execute(engine,data,root/'exact-native',True)
                control=execute(engine,data,root/'control-native',True,threads=2)
                candidate=execute(engine,data,root/'default-native',False)
                for name,packet in [('exact',exact),('control',control),('default',candidate)]:save_result_packet(packet,root/name)
                masks=selection(exact);np.savez(root/'selection.npz',**masks)
                control_report=qualification_report(control,exact,rtol=SURVEY_EXACT_RTOL,atol=SURVEY_ATOL)
                full=qualification_report(candidate,exact,rtol=SURVEY_RTOL,atol=SURVEY_ATOL)
                science=qualification_report(candidate,exact,rtol=SURVEY_RTOL,atol=SURVEY_ATOL,selection=masks,metric='pointwise')
                if set(candidate['arrays'])-set(full['observables']):raise AssertionError('an observable was omitted')
                record=dict(engine=engine,fixture=kind,exact_control=control_report,all_bins=full,science_bins=science,
                    exact_wall=exact['metadata']['timings']['wall_seconds'],default_wall=candidate['metadata']['timings']['wall_seconds'])
                report['cases'].append(record)
                print(engine,kind,'control='+control_report['status'],'all='+full['status'],'science='+science['status'],flush=True)
        report['status']='PASS' if len(report['cases'])==len(ENGINES)*len(FIXTURES) and all(c['exact_control']['status']=='QUALIFIED' for c in report['cases']) else 'FAIL'
        report['science_acceptance']='QUALIFIED' if report['status']=='PASS' and all(c['science_bins']['status']=='QUALIFIED' for c in report['cases']) else 'REJECTED'
    finally:
        (output/'campaign.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    return report

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    raise SystemExit(0 if run(a.output)['status']=='PASS' else 1)
