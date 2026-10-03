"""Independent product totals, recovery, cache isolation and resource plans."""
from pathlib import Path
import numpy as np
import pytest
import cyballs


def forest(root,method='lya-los-tree-2pcf-3pcf-omp',no_out=True):
    m=cyballs.cballs();m.set(searchMethod=method,rootDir=str(root),numberThreads=2,
        verbose=0,verbose_log=0,usePeriodic=False,useLogHist=False,rangeN=30.,rminHist=.1,sizeHistN=4,
        lya2RpMax=30.,lya2RtMax=30.,lya2RpBins=5,lya2RtBins=6,
        lya3RMax=30.,lya3RBins=4,lya3ThetaBins=5,lya3MuBins=6,
        options='no-out-Hist' if no_out else '')
    m.set_forest_catalog(np.array([[10.,0,0],[9.,2,0],[8.,0,3]]),np.array([2.,3.,5.]),np.array([1.,2.,4.]),np.arange(3,dtype=np.int64))
    m.Run();return m


@pytest.mark.parametrize('method',['lya-2pcf-3pcf-omp','lya-los-tree-2pcf-3pcf-omp','lya-1d-2pcf-3pcf-omp','lya-1d-tree-2pcf-omp','lya-1d-tree-3pcf-omp'])
def test_native_forest_products_and_lifetime(tmp_path,method):
    m=forest(tmp_path,method);arrays=m.getForestResults()['arrays']
    for key,total in [('pair_numerator',172.),('pair_denominator',14.),('triple_numerator',1440.),('triple_denominator',48.)]:
        if key in arrays:np.testing.assert_allclose(arrays[key].sum(),total,rtol=3e-13)
    assert m.getAllocationInfo()['retained_result_bytes']>0
    snapshot={k:v.copy() for k,v in arrays.items()};arrays[next(iter(arrays))].flat[0]+=999
    assert any(not np.array_equal(arrays[k],v) for k,v in snapshot.items())
    for k,v in snapshot.items():np.testing.assert_array_equal(m.getForestResults()['arrays'][k],v)
    with pytest.raises(Exception,match='does not match'):m.getPhysicalResults()
    m.struct_cleanup();assert m.getAllocationInfo()['retained_result_bytes']==0
    with pytest.raises(Exception,match='unavailable'):m.getForestResults()
    for k,v in snapshot.items():assert np.isfinite(v).all()


def test_physical_native_raw_arrays_match_export(tmp_path):
    m=forest(tmp_path,'octree-3pcf-3d-omp',no_out=False)
    arrays=m.getPhysicalResults()['arrays'];table=np.loadtxt(tmp_path/'histZetaM_3d.txt')
    np.testing.assert_allclose(arrays['triple_numerator'].ravel(),table[:,-2],rtol=2e-11,atol=1e-11)
    np.testing.assert_allclose(arrays['triple_denominator'].ravel(),table[:,-1],rtol=2e-11,atol=1e-11)
    m.struct_cleanup()


@pytest.mark.parametrize('engine',['octree-2balls-omp','balltree-2balls-omp'])
def test_cache_isolation_and_release(tmp_path,engine):
    rng=np.random.default_rng(41);p=rng.normal(size=(120,3));p/=np.linalg.norm(p,axis=1)[:,None]
    a,b=cyballs.cballs(),cyballs.cballs()
    for index,m in enumerate((a,b)):
        m.set(searchMethod=engine,rootDir=str(tmp_path/str(index)),numberThreads=2,
            options='only-2pcf,no-out-Hist,no-smooth-pivot',theta=0,verbose=0,verbose_log=0)
        m.set_catalog(p,kappa=np.ones(120));m.Run()
    ab=a.getCacheInfo()['bytes'];bb=b.getCacheInfo()['bytes'];assert ab>0 and bb>0
    ref=b.getHistXi2pcf().copy();a.struct_cleanup();assert a.getCacheInfo()['bytes']==0
    assert b.getCacheInfo()['bytes']==bb;np.testing.assert_array_equal(b.getHistXi2pcf(),ref)
    b.set(rootDir=str(tmp_path/'warm'));b.Run()
    assert b.getCacheInfo()['bytes']==bb
    np.testing.assert_allclose(b.getHistXi2pcf(),ref,rtol=3e-13,atol=1e-13)
    b.clearCaches();assert b.getCacheInfo()['bytes']==0;np.testing.assert_array_equal(b.getHistXi2pcf(),ref)
    b.struct_cleanup()


def test_forecast_default_overflow_and_budget(monkeypatch,tmp_path):
    monkeypatch.delenv('CBALLS_MEMORY_BUDGET_MB',raising=False)
    assert cyballs.resource_policy()['default_budget_mib']==65536
    assert cyballs.resource_policy()['budget_bytes']==65536*1048576
    plan=cyballs.resource_plan('lya-3pcf-omp',1000,4)
    assert plan['known_components_bytes']['estimator_global_and_worker_histograms']==800000*(16+4*24)
    with pytest.raises(OverflowError):cyballs.resource_plan('lya-3pcf-omp',2**63,4)
    monkeypatch.setenv('CBALLS_MEMORY_BUDGET_MB','1')
    assert cyballs.resource_plan('lya-3pcf-omp',1000,4)['known_exceeds_budget']
    monkeypatch.setenv('CBALLS_MEMORY_BUDGET_MB','bad')
    with pytest.raises(ValueError):cyballs.resource_plan('lya-3pcf-omp',1000,4)
    monkeypatch.setenv('CBALLS_MEMORY_BUDGET_MB','64')
    m=forest(tmp_path);m.struct_cleanup()


def test_same_los_count_dtype_and_equal_forest_mean(tmp_path):
    m=cyballs.cballs();m.set(searchMethod='lya-1d-tree-same-los-2pcf-omp',rootDir=str(tmp_path),
        usePeriodic=False,numberThreads=2,lya2RpMax=10.,lya2RpBins=2,verbose=0,verbose_log=0,options='no-out-Hist')
    p=np.zeros((4,3));p[:,0]=[100,101,103,104]
    m.set_forest_catalog(p,np.array([2.,3.,5.,7.]),np.array([1.,2.,4.,8.]),np.array([1,1,2,2],dtype=np.int64))
    m.Run();a=m.getForestResults()['arrays']
    assert a['contributing_forests'].dtype==np.uint64
    assert a['contributing_forests'][0]==2
    np.testing.assert_allclose(a['correlation_sum'][0],41.)
    m.struct_cleanup()


@pytest.mark.parametrize('engine,option',[('octree-2balls-omp','no-native-tree-cache'),('balltree-2balls-omp','no-balltree-tree-cache')])
def test_uncached_tree_release_does_not_create_owner(tmp_path,engine,option):
    rng=np.random.default_rng(413);p=rng.normal(size=(80,3));p/=np.linalg.norm(p,axis=1)[:,None]
    m=cyballs.cballs()
    m.set(searchMethod=engine,rootDir=str(tmp_path),numberThreads=2,
          options='only-2pcf,no-out-Hist,no-smooth-pivot,'+option,theta=0,verbose=0,verbose_log=0)
    m.set_catalog(p,kappa=np.ones(80));m.Run()
    assert m.getCacheInfo()['bytes']==0
    m.struct_cleanup();assert m.getCacheInfo()['bytes']==0
