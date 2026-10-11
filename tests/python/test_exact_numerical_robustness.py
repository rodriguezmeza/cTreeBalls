"""Independent Decimal tuple sums and extreme-range/recovery contracts.

No project numerical reference helper is imported. Geometry is evaluated in
binary64; weight/field products and summation use Decimal at 100 digits.
"""
from decimal import Decimal, localcontext
from pathlib import Path
import math
import os
import sys
import tempfile
import numpy as np
import pytest
from cyballs import cballs

SCALAR=['octree-2balls-omp','kdtree-2balls-omp','balltree-2balls-omp']
PHYSICAL='octree-3pcf-3d-omp'
ENGINES=SCALAR+[PHYSICAL]
POINTS=np.vstack((np.eye(3),-np.ones(3)/np.sqrt(3)))
EDGES=np.geomspace(.03,1.99,5)


def settings(engine, path, order=4, threads=1, pair=False):
    physical='3pcf-3d' in engine
    options='no-out-Hist,no-smooth-pivot,weights-norm,no-normalize-HistZeta,no-one-ball,no-two-balls'
    options+=(',only-2pcf-3d' if pair else ',compute-2pcf-3d,compute-3pcf-3d') if physical else (',only-2pcf' if pair else ',edge-corrections')
    return dict(searchMethod=engine,rootDir=str(path),sizeHistN=4,mChebyshev=max(2,order),
                rangeN=1.99,rminHist=.03,useLogHist=True,usePeriodic=False,theta=0,
                numberThreads=threads,verbose=0,verbose_log=0,options=options)


def extract(model,engine,order,pair=False):
    if '-mpi' in engine and __import__('mpi4py.MPI',fromlist=['COMM_WORLD']).COMM_WORLD.rank!=0:return {}
    if '3pcf-3d' in engine:return model.getPhysicalResults()['arrays']
    result={'pair':model.getHistXi2pcf(),'count':model.getHistNN()}
    if not pair:
        result['triple']=np.array([model.getHistZetaMsincos(m+1,1)+model.getHistZetaMsincos(m+1,2)
          +1j*(model.getHistZetaMsincos(m+1,3)-model.getHistZetaMsincos(m+1,4)) for m in range(order+1)])
        result['denominator']=model.getScalarWindowDiagnostics()['window_monopole']
    return result


def run(engine,positions,field,weights,path,order=4,threads=1,pair=False):
    model=cballs()
    try:
        model.set(settings(engine,path,order,threads,pair))
        model.set_catalog(positions,kappa=field,weights=weights);model.Run()
        return extract(model,engine,order,pair)
    finally:model.clean_all()


def reference(positions,field,weights,order,physical,los_ids=None):
    def D(x):return Decimal.from_float(float(x))
    with localcontext() as ctx:
        ctx.prec=100
        n=len(weights);w=list(map(D,weights));f=[w[i]*D(field[i]) for i in range(n)]
        re=[Decimal(0)]*((order+1)*16);im=list(re);den=[Decimal(0)]*16
        for i in range(n):
            dr=positions-positions[i];r=np.linalg.norm(dr,axis=1)
            bins=np.searchsorted(EDGES,r,side='right')-1
            selected=[j for j in range(n) if j!=i and EDGES[0]<r[j]<EDGES[-1]
                      and (los_ids is None or los_ids[j]!=los_ids[i])]
            normal=positions[i]/np.linalg.norm(positions[i])
            for j in selected:
                for k in selected:
                    if j==k or (los_ids is not None and los_ids[j]==los_ids[k]):continue
                    cell=int(bins[j]*4+bins[k]);den[cell]+=w[i]*w[j]*w[k]
                    if physical:
                        mu=np.clip(np.dot(dr[j],dr[k])/(r[j]*r[k]),-1,1)
                        poly=[1.,mu]
                        for m in range(1,order):poly.append(((2*m+1)*mu*poly[m]-m*poly[m-1])/(m+1))
                        angular=[complex((2*m+1)*poly[m]) for m in range(order+1)]
                    else:
                        a=dr[j]-normal*np.dot(dr[j],normal);b=dr[k]-normal*np.dot(dr[k],normal)
                        phi=np.arctan2(np.dot(normal,np.cross(a,b)),np.dot(a,b))
                        angular=np.exp(1j*np.arange(order+1)*phi)
                    product=f[i]*f[j]*f[k]
                    for m,z in enumerate(angular):
                        re[m*16+cell]+=product*D(z.real);im[m*16+cell]+=product*D(z.imag)
        num=np.array([complex(float(a),float(b)) for a,b in zip(re,im)]).reshape(order+1,4,4)
        return num.real if physical else num,np.array(list(map(float,den))).reshape(4,4)


def compare(actual,expected):
    # Same strict ordinary-exact criterion; no tolerance relaxation for H=1e16.
    assert np.isfinite(actual).all()
    error=np.max(np.abs(actual-expected)); scale=np.max(np.abs(expected))
    assert error <= 3e-12*scale+1e-300,(error,scale)


@pytest.mark.parametrize('engine',ENGINES)
@pytest.mark.parametrize('heavy',[1.,1e4,1e8,1e12,1e16])
def test_four_point_positive_monopole(tmp_path,engine,heavy):
    a=run(engine,POINTS,np.ones(4),np.array([heavy,1,1,1]),tmp_path,order=0)
    target=float(Decimal(18)*Decimal.from_float(heavy)+Decimal(6))
    if '3pcf-3d' in engine:
        compare(a['triple_numerator'][0].sum(),target)
        compare(a['triple_denominator'][0].sum(),target)
        occupied=a['triple_denominator'][0]>0
        compare(a['triple_numerator'][0][occupied]/a['triple_denominator'][0][occupied],np.ones(occupied.sum()))
    else:
        compare(a['triple'][0].real.sum(),target)
        compare(a['denominator'].sum(),target)


@pytest.mark.parametrize('engine',ENGINES)
@pytest.mark.parametrize('threads',[1,2,4])
def test_signed_multibin_permutations(tmp_path,engine,threads):
    rng=np.random.default_rng(552026);p=rng.normal(size=(9,3));p/=np.linalg.norm(p,axis=1)[:,None]
    f=rng.uniform(-1,1,9);w=rng.uniform(.5,2,9);w[0]=1e16
    target,den=reference(p,f,w,4,'3pcf-3d' in engine)
    for index,perm in enumerate([np.arange(9),np.arange(8,-1,-1),rng.permutation(9)]):
        a=run(engine,p[perm],f[perm],w[perm],tmp_path/str(index),threads=threads)
        compare(a['triple_numerator'] if '3pcf-3d' in engine else a['triple'],target)
        compare(a['triple_denominator'][0] if '3pcf-3d' in engine else a['denominator'],den)


@pytest.mark.parametrize('engine',SCALAR)
@pytest.mark.parametrize('weight',[0.,1.,1e-300,1e-200,1e155,1e200,1e300])
def test_normalized_pair_scale_and_empty_bins(tmp_path,engine,weight):
    a=run(engine,POINTS,np.ones(4),np.full(4,weight),tmp_path,pair=True)
    target=np.zeros(4);target[-1]=1. if weight else 0.
    np.testing.assert_array_equal(a['pair'],target)
    np.testing.assert_array_equal(a['count'],[0.,0.,0.,6.])


@pytest.mark.parametrize('engine',ENGINES)
@pytest.mark.parametrize('weight',[1e-100,1e-105])
def test_representable_small_raw_products(tmp_path,engine,weight):
    target=float(Decimal(24)*Decimal.from_float(weight)**3)
    p=np.random.default_rng(1922).normal(size=(4,3));p/=np.linalg.norm(p,axis=1)[:,None]
    a=run(engine,p,np.ones(4),np.full(4,weight),tmp_path,order=0)
    num=a['triple_numerator'] if '3pcf-3d' in engine else a['triple'].real
    den=a['triple_denominator'][0] if '3pcf-3d' in engine else a['denominator']
    # Relative quantization error of a subnormal cannot satisfy normal ulp bounds.
    np.testing.assert_allclose(num[0].sum()/target,1,rtol=3e-12,atol=64*np.nextafter(0.,1.)/target)
    np.testing.assert_allclose(den.sum()/target,1,rtol=3e-12,atol=64*np.nextafter(0.,1.)/target)


@pytest.mark.parametrize('engine',ENGINES)
@pytest.mark.parametrize('weight',[1e-120,1e110,1e155,1e200])
def test_raw_range_failure_and_same_object_recovery(tmp_path,engine,weight):
    model=cballs()
    try:
        model.set(settings(engine,tmp_path,order=0))
        model.set_catalog(POINTS,kappa=np.ones(4),weights=np.full(4,weight))
        with pytest.raises(Exception,match='numerical range'):
            model.Run()
        assert not model.state
        model.set_catalog(POINTS,kappa=np.ones(4),weights=np.ones(4));model.Run()
        assert model.state
        a=extract(model,engine,0)
        if '-mpi' not in engine or __import__('mpi4py.MPI',fromlist=['COMM_WORLD']).COMM_WORLD.rank==0:
            compare((a['triple_numerator'] if '3pcf-3d' in engine else a['triple'])[0].real.sum(),24.)
    finally:model.clean_all()



@pytest.mark.parametrize('engine',ENGINES)
def test_raw_overflow_boundary(tmp_path,engine):
    a=run(engine,POINTS,np.ones(4),np.full(4,1e102),tmp_path/'safe',order=0)
    num=a['triple_numerator'] if '3pcf-3d' in engine else a['triple'].real
    compare(num[0].sum(),2.4e307)
    test_raw_range_failure_and_same_object_recovery(tmp_path/'overflow',engine,2e102)


@pytest.mark.parametrize('engine',SCALAR)
def test_heterogeneous_scaled_pairs(tmp_path,engine):
    w=np.array([1e300,1e-300,2e-300,3e-300]);f=np.array([1.,-2.,3.,-4.])
    with localcontext() as ctx:
        ctx.prec=100
        D=Decimal.from_float
        den=sum(D(w[i])*D(w[j]) for i in range(4) for j in range(i+1,4))
        num=sum(D(w[i])*D(w[j])*D(f[i])*D(f[j]) for i in range(4) for j in range(i+1,4))
        target=float(num/den)
    a=run(engine,POINTS,f,w,tmp_path,pair=True)
    compare(a['pair'][-1],target)


@pytest.mark.parametrize('engine',SCALAR)
@pytest.mark.parametrize('field',[1e-200,1e200])
def test_pair_result_outside_range_recovers(tmp_path,engine,field):
    m=cballs()
    try:
        m.set(settings(engine,tmp_path,pair=True))
        m.set_catalog(POINTS,kappa=np.full(4,field),weights=np.ones(4))
        with pytest.raises(Exception,match='numerical range'):m.Run()
        assert not m.state
        m.set_catalog(POINTS,kappa=np.ones(4),weights=np.ones(4));m.Run()
        np.testing.assert_array_equal(m.getHistXi2pcf(),[0,0,0,1])
    finally:m.clean_all()


@pytest.mark.parametrize('threads',[1,2,4])
def test_distinct_los_high_dynamic_range(tmp_path,threads):
    rng=np.random.default_rng(121)
    p=rng.normal(size=(9,3));p/=np.linalg.norm(p,axis=1)[:,None]
    f=rng.uniform(-1.,1.,9);w=np.ones(9);w[0]=1e16
    ids=np.repeat(np.arange(3,dtype=np.int64),3)
    num,den=reference(p,f,w,4,True,ids)
    m=cballs()
    try:
        params=settings(PHYSICAL,tmp_path,threads=threads)
        params['options']+=',exclude-all-same-los'
        m.set(params);m.set_forest_catalog(p,f,w,ids);m.Run()
        a=extract(m,PHYSICAL,4)
        compare(a['triple_numerator'],num);compare(a['triple_denominator'][0],den)
    finally:m.clean_all()


if __name__=='__main__':
    import argparse,json
    parser=argparse.ArgumentParser();parser.add_argument('--mpi',action='store_true');args=parser.parse_args()
    if not args.mpi:raise SystemExit('use pytest, or --mpi under mpiexec')
    from mpi4py import MPI
    rank=MPI.COMM_WORLD.rank
    with tempfile.TemporaryDirectory(prefix='exact-robustness-mpi-') as directory:
        base=Path(directory)
        for engine0 in ENGINES:
            engine=engine0.replace('-omp','-mpi')
            for threads in [1,2]:
                a=run(engine,POINTS,np.ones(4),np.array([1e16,1,1,1]),base/engine,order=0,threads=threads)
                if rank==0:
                    compare((a['triple_numerator'] if engine0==PHYSICAL else a['triple'])[0].real.sum(),18e16+6)
                    compare((a['triple_denominator'][0] if engine0==PHYSICAL else a['denominator']).sum(),18e16+6)
            rng=np.random.default_rng(552026);p=rng.normal(size=(9,3));p/=np.linalg.norm(p,axis=1)[:,None]
            f=rng.uniform(-1,1,9);w=rng.uniform(.5,2,9);w[0]=1e16
            target,den=reference(p,f,w,4,engine0==PHYSICAL)
            a=run(engine,p[::-1],f[::-1],w[::-1],base/engine/'signed',threads=2)
            if rank==0:
                compare(a['triple_numerator'] if engine0==PHYSICAL else a['triple'],target)
                compare(a['triple_denominator'][0] if engine0==PHYSICAL else a['denominator'],den)
            if engine0!=PHYSICAL:
                a=run(engine,POINTS,np.ones(4),np.full(4,1e200),base/engine/'pairs',pair=True)
                if rank==0:np.testing.assert_array_equal(a['pair'],[0,0,0,1])
            # All ranks must reject the same arithmetic phase, then recover.
            test_raw_range_failure_and_same_object_recovery(base,engine,1e155)
        if rank==0:print('PASS: MPI exact numerical robustness and recovery')
