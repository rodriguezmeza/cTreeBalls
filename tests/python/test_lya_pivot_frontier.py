#!/usr/bin/env python3
"""Independent per-original-pixel oracle for forest-local pivot geometry.

The oracle does not aggregate weights: each original pivot contributes with its
original field, weight and row-ID ownership at its representative's coordinates.
This separately tests the C product-of-sums and multiplicity contracts.
"""
from __future__ import annotations
import itertools
import json
import math
import os
from pathlib import Path
import re
import subprocess
import numpy as np
import pytest

ROOT=Path(__file__).resolve().parents[2]
BINARY=Path(os.environ.get('CBALLS',ROOT/'cballs')).resolve()
METHODS=['lya-2pcf-omp','lya-3pcf-omp','lya-2pcf-3pcf-omp',
         'lya-los-tree-2pcf-omp','lya-los-tree-3pcf-omp','lya-los-tree-2pcf-3pcf-omp']

def fixture():
    rng=np.random.default_rng(4681)
    rows=[]
    for f,direction in enumerate(((1,.03,.04),(1,-.1,.03),(1,.01,-.1),(1,.12,.05))):
        u=np.asarray(direction)/np.linalg.norm(direction)
        for j in range(7):
            # Deliberately bent forests, signed fields and accepted zero weights.
            x=u*(90+1.2*j)+rng.normal(scale=.013,size=3)
            rows.append([*x,rng.normal(),0 if j==2 else rng.uniform(.2,1.6),f*37-119])
    rows=np.array(rows);rng.shuffle(rows);return rows

def representatives(rows,radius,cap):
    ids=list(range(len(rows)));result=np.arange(len(rows))
    if radius==0:return result
    ids.sort(key=lambda i:(rows[i,5],np.linalg.norm(rows[i,:3]),i))
    while ids:
        rep=ids.pop(0);group=[rep]
        while ids and len(group)<cap and rows[ids[0],5]==rows[rep,5] and np.linalg.norm(rows[ids[0],:3]-rows[rep,:3])<=radius:
            group.append(ids.pop(0))
        for i in group:result[i]=rep
    return result

def oracle(rows,radius=0,cap=8,rpmax=25.,rtmax=25.,rmax=25.):
    rep=representatives(rows,radius,cap);pair={};trip={};pc=tc=0
    def positive(x,limit,bins):
        return int(x/limit*bins) if 0<=x<limit else None
    def deposit(table,key,n,d):
        v=table.setdefault(key,np.zeros(2));v[0]+=n;v[1]+=d
    for i,row in enumerate(rows):
        x=rows[rep[i],:3];chi=np.linalg.norm(x);los=x/chi
        for j in range(i+1,len(rows)):
            q=rows[j]
            if row[5]==q[5] or np.linalg.norm(q[:3]-x)>=math.hypot(rpmax,rtmax):continue
            r=np.linalg.norm(q[:3]);c=np.clip(np.dot(los,q[:3]/r),-1,1)
            bp=positive(abs(chi-r)*math.sqrt((1+c)/2),rpmax,5)
            bt=positive((chi+r)*math.sqrt((1-c)/2),rtmax,6)
            if bp is None or bt is None:continue
            pc+=1;deposit(pair,(bp,bt),row[3]*row[4]*q[3]*q[4],row[4]*q[4])
        neighbors=[]
        for q in rows:
            if row[5]==q[5]:continue
            dr=q[:3]-x;r=np.linalg.norm(dr)
            if not 0<r<rmax:continue
            b=positive(r,rmax,4)
            theta=math.acos(np.clip(np.dot(dr,los)/r,-1,1))
            t=min(int(theta/math.pi*5),4)
            neighbors.append((q,dr,r,b,t))
        for (q,a,ra,ba,ta),(r,b,rb,bb,tb) in itertools.combinations(neighbors,2):
            if q[5]==r[5]:continue
            mu=np.clip(np.dot(a,b)/(ra*rb),-1,1);m=min(int((mu+1)*.5*6),5)
            n=row[3]*row[4]*q[3]*q[4]*r[3]*r[4];d=row[4]*q[4]*r[4]
            tc+=2
            deposit(trip,(ba,bb,ta,tb,m),n,d);deposit(trip,(bb,ba,tb,ta,m),n,d)
    return ({k:v for k,v in pair.items() if v[1]>0}, {k:v for k,v in trip.items() if v[1]>0},pc,tc)

def run(tmp_path,method,rows,threads=1,radius=0,level=0,cap=8,kernel=0,extra=()):
    tmp_path.mkdir(parents=True)
    path=tmp_path/'pixels.txt';np.savetxt(path,rows,fmt=('%.17g',)*5+('%.0f',))
    args=dict(search=method,infile=path,infileformat='lya-ascii',iCatalogs=1,rootDir=tmp_path/'out',numberThreads=threads,
              usePeriodic='false',useLogHist='false',rangeN=25,rminHist=.1,sizeHistN=4,
              lya2RpMax=25,lya2RtMax=25,lya2RpBins=5,lya2RtBins=6,
              lya3RMax=25,lya3RBins=4,lya3ThetaBins=5,lya3MuBins=6,
              lyaScanLevel=level,lyaPivotRadius=radius,lyaPivotMax=cap,lya3Kernel=kernel,
              verbose=2,options='no-smooth-pivot')
    args.update(dict(x.split('=',1) for x in extra))
    p=subprocess.run([str(BINARY),*(f'{k}={v}' for k,v in args.items())],text=True,capture_output=True,timeout=120)
    assert p.returncode==0,p.stdout+p.stderr
    return tmp_path/'out',p.stdout

def read(path,order):
    table={}
    for line in path.read_text().splitlines():
        if not line or line.startswith('#'):continue
        v=line.split();key=tuple(map(int,v[:2 if order==2 else 5]));value=np.array(list(map(float,v[-2:])))
        if value[1]>0:table[key]=value
    count=int(re.search(r'# distinct-forest .*: (\d+)',path.read_text())[1])
    return table,count

def compare(out,method,expected):
    for order,name in ((2,'histXi2pcf_lya.txt'),(3,'histZetaM_lya5d.txt')):
        if f'{order}pcf' not in method:continue
        table,count=read(out/name,order);want=expected[order-2]
        assert table.keys()==want.keys()
        for key in table:np.testing.assert_allclose(table[key],want[key],rtol=8e-12,atol=3e-11)
        assert count==expected[order]

@pytest.mark.parametrize('method',METHODS)
@pytest.mark.parametrize('radius',[0.,3.7])
def test_independent_oracle_and_threads(tmp_path,method,radius):
    rows=fixture();expected=oracle(rows,radius)
    one,_=run(tmp_path/'one',method,rows,radius=radius,level=2)
    many,_=run(tmp_path/'many',method,rows,radius=radius,level=2,threads=4)
    compare(one,method,expected);compare(many,method,expected)
    for p in one.glob('hist*.txt'):assert p.read_bytes()==(many/p.name).read_bytes()

@pytest.mark.parametrize('kernel',[0,1,2])
def test_smooth_triplet_kernels_and_combined_ownership(tmp_path,kernel):
    rows=fixture();out,log=run(tmp_path/'run','lya-2pcf-3pcf-omp',rows,radius=3.7,level=1,kernel=kernel,cap=3)
    compare(out,'lya-2pcf-3pcf-omp',oracle(rows,3.7,3))
    assert int(re.search(r'representatives=(\d+)',log)[1])<len(rows)

@pytest.mark.parametrize('case',['zero-weights','one-forest','coincident','boundary','antipodal','cancellation'])
def test_adversarial(tmp_path,case):
    rows=fixture()
    if case=='zero-weights':rows[:,4]=0
    if case=='one-forest':rows[:,5]=17
    if case=='coincident':
        rows[:,0:3]=[90,1,2]
        with pytest.raises(AssertionError,match='two bodies have same position'):
            run(tmp_path/'duplicate','lya-2pcf-3pcf-omp',rows,radius=4,level=3)
        return
    if case=='boundary':rows[0,:3]=rows[1,:3]+[25,0,0]
    if case=='antipodal':rows[::2,:3]*=-1
    if case=='cancellation':rows[:,4]=1;rows[:,3]=np.where(np.arange(len(rows))%2,1.,-1.)
    out,_=run(tmp_path/'run','lya-2pcf-3pcf-omp',rows,radius=4,level=3)
    compare(out,'lya-2pcf-3pcf-omp',oracle(rows,4))

@pytest.mark.parametrize('setting',[('lyaScanLevel=-1',),('lyaScanLevel=21',),('lyaPivotRadius=-1',),('lyaPivotRadius=nan',),('lyaPivotMax=0',),('lyaPivotMax=1025',),('lyaScanLevel=2','lya3Kernel=4')])
def test_invalid_controls(tmp_path,setting):
    with pytest.raises(AssertionError,match='require|invalid numeric value'):
        run(tmp_path/'bad','lya-2pcf-3pcf-omp',fixture(),extra=setting)

def test_single_member_cap_and_combined_exact_pair(tmp_path):
    rows=fixture();exact=oracle(rows)
    out,_=run(tmp_path/'one','lya-2pcf-3pcf-omp',rows,radius=30,cap=1,level=2)
    compare(out,'lya-2pcf-3pcf-omp',exact)
    mixed,_=run(tmp_path/'mixed','lya-2pcf-3pcf-omp',rows,radius=3.7,level=2,extra=('lya2Kernel=1',))
    approx=oracle(rows,3.7)
    compare(mixed,'lya-2pcf-3pcf-omp',(exact[0],approx[1],exact[2],approx[3]))

def test_frontier_budget_preflight(tmp_path,monkeypatch):
    rng=np.random.default_rng(44);rows=np.column_stack((rng.normal(size=(10000,3))+[100,0,0],np.ones((10000,2)),np.arange(10000)))
    # Select dimensions from the compiled body layout: extra numerical moments
    # may grow the catalog. Fail the frontier plan, not the earlier base plan.
    cyballs=pytest.importorskip('cyballs')
    monkeypatch.setenv('CBALLS_MEMORY_BUDGET_MB','7')
    budget=7*1024*1024
    frontier=(3*len(rows)+1)*cyballs.resource_policy()['type_bytes']['size_t']
    bins=next(b for b in range(4,1000) if
        budget-frontier < cyballs.resource_plan('lya-2pcf-omp',len(rows),parameters={
            'sizeHistN':4,'lya2RpBins':b,'lya2RtBins':b})['known_total_bytes'] <= budget)
    with pytest.raises(AssertionError,match='pivot frontier'):
        run(tmp_path/'budget','lya-2pcf-omp',rows,level=2,extra=(f'lya2RpBins={bins}',f'lya2RtBins={bins}'))

@pytest.mark.parametrize('method,kernel',[
    ('lya-2pcf-omp',1),('lya-2pcf-mpi',0),('lya-1d-2pcf-omp',0),
    ('lya-anisotropic-multipole-3pcf-omp',0)])
def test_unsupported_paths_are_explicit(tmp_path,method,kernel):
    with pytest.raises(AssertionError,match='pivot frontier/smoothing requires'):
        run(tmp_path/'bad',method,fixture(),level=1,extra=(f'lya2Kernel={kernel}',))

def test_cython_controls_oracle_and_reuse(tmp_path):
    cyballs=pytest.importorskip('cyballs');rows=fixture();balls=cyballs.cballs()
    try:
        for i,radius in enumerate((3.7,0.,3.7)):
            out=tmp_path/str(i)
            balls.set(dict(searchMethod='lya-2pcf-3pcf-omp',rootDir=str(out),numberThreads=3,
                usePeriodic=False,useLogHist=False,rangeN=25.,rminHist=.1,sizeHistN=4,
                lya2RpMax=25.,lya2RtMax=25.,lya2RpBins=5,lya2RtBins=6,
                lya3RMax=25.,lya3RBins=4,lya3ThetaBins=5,lya3MuBins=6,
                lyaScanLevel=0 if i==0 else 2,lyaPivotRadius=radius,lyaPivotMax=8,
                options='no-smooth-pivot',verbose=0,verbose_log=0))
            balls.set_forest_catalog(rows[:,:3],rows[:,3],rows[:,4],rows[:,5].astype(np.int64))
            balls.Run(level=['MainLoop'])
            compare(out,'lya-2pcf-3pcf-omp',oracle(rows,radius))
            meta=balls.getRunMetadata()
            assert meta['lya_pivot_frontier']['smoothing_radius']==radius
            assert meta['lya_pivot_frontier']['scan_level']==(0 if i==0 else 2)
            assert meta['lya_pivot_frontier']['geometry_approximate']==bool(radius)
            assert meta['lya_2pcf']['approximate']==bool(radius)
            balls.struct_cleanup()
            assert balls.getRunMetadata()==meta
    finally:balls.struct_cleanup()

def test_accuracy_acceptance_does_not_hide_bad_bins():
    import sys
    sys.path.insert(0,str(ROOT/'tests'/'python'))
    from benchmark_lya_pivot_frontier import comparison
    # axis, correlation, numerator, denominator; one normal, small and empty bin.
    ref=np.array([[0,.2,2,10],[1,1e-8,1e-7,10],[2,0,0,0]],float)
    def check(candidate,exact=False):
        result=comparison(ref,candidate,1e-6,.05,5e-8,exact)
        json.dumps(result,allow_nan=False)
        return result['accepted']
    assert check(ref,True)
    for index,value in [(0,.22),(1,1e-6),(2,.1)]:
        bad=ref.copy();bad[index,1]=value;bad[index,-1]=10
        assert not check(bad)
    bad=ref.copy();bad[0,1]=np.nan
    assert not check(bad)
    bad=ref.copy();bad[0,-1]=0
    assert not check(bad)

def test_signed_64_bit_forest_ids(tmp_path):
    rows=fixture();cat=tmp_path/'wide.txt'
    wide=[2**62+1,2**62+2,-2**62-1,-2**62-2]
    mapping=dict(zip(np.unique(rows[:,5]),wide))
    cat.write_text(''.join(' '.join(format(v,'.17g') for v in row[:5])+f' {mapping[row[5]]}\n' for row in rows))
    out,_=run(tmp_path/'run','lya-2pcf-3pcf-omp',rows,radius=3.7,level=2,extra=(f'infile={cat}',))
    compare(out,'lya-2pcf-3pcf-omp',oracle(rows,3.7))
