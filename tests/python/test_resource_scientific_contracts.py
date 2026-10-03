"""Resource, lifecycle and timing regressions; huge requests never need huge inputs."""
import numpy as np
import pytest
from cyballs import cballs, search_method_id


def model(tmp_path, engine='octree-2balls-omp', **changes):
    m=cballs()
    m.set(searchMethod=engine, rootDir=str(tmp_path), sizeHistN=4, mChebyshev=3,
          rangeN=1.8,rminHist=.02,useLogHist=True,numberThreads=2,verbose=0,verbose_log=0,
          options='no-out-Hist,no-smooth-pivot,no-one-ball,no-two-balls,no-normalize-HistZeta')
    p=np.array([[1,0,0],[.8,.6,0],[.6,0,.8],[0,1,0]],float)
    m.set_catalog(p,kappa=np.array([1.,-.2,.3,.8]))
    if changes: m.set(**changes)
    return m


@pytest.mark.parametrize('engine',['octree-2balls-omp','kdtree-2balls-omp','balltree-2balls-omp',
                                  'octree-shear-sphere-2balls-omp','kdtree-shear-sphere-2balls-omp',
                                  'balltree-shear-sphere-2balls-omp'])
def test_pair_plan_has_no_scalar_tensors(tmp_path,engine):
    if search_method_id(engine)<0: pytest.skip('engine absent from this profile')
    m=model(tmp_path,engine,sizeHistN=512,mChebyshev=31,
            options='only-2pcf,no-out-Hist,no-smooth-pivot')
    try:
        # Startup is enough to inspect allocations, avoiding massive scientific outputs.
        m.Run(level=['StartRun_Common'])
        p=m.getAllocationInfo()
        assert p['live'] and not p['scalar_3pcf_planned']
        assert not p['scalar_tensor_allocated'] and not p['square_export_allocated']
        assert p['common_histogram_bytes']<100000
        assert not m.state
    finally: m.struct_cleanup()
    assert not m.getAllocationInfo()['live']


def test_lifecycle_timings_and_missing_product(tmp_path):
    m=cballs(); assert not m.state
    with pytest.raises(Exception,match='timings unavailable'): m.getTimings()
    m=model(tmp_path)
    try:
        m.Run(level=['StartRun_Common']); assert not m.state
        returned=m.Run(); assert m.state
        t=m.getTimings()
        assert t['wall_seconds']>0 and t['process_cpu_seconds']>=0
        assert returned==t['legacy_cpu_per_requested_thread_seconds']
        assert t['process_cpu_seconds']/t['requested_threads']==returned
        assert m.Run()==returned and m.getTimings()==t
        assert m.getAllocationInfo()['scalar_tensor_allocated']
        metadata=m.getRunMetadata()
        assert metadata['geometry']['scalar_observer_frame']
        assert metadata['resources']['common_scalar_3pcf']
        m.struct_cleanup(); assert not m.state and m.getTimings()==t
        m.set(options='only-2pcf,no-out-Hist,no-smooth-pivot')
        assert not m.state
        with pytest.raises(Exception,match='timings unavailable'):m.getTimings()
        m.Run(); assert m.state
        for f in (lambda:m.getHistZetaMsincos(1,1),lambda:m.getHistZetaM_EE(1),lambda:m.getHistZetaM_EE_Im(1)):
            with pytest.raises(Exception,match='not computed'):f()
        m.set(sizeHistN=2147483647)
        with pytest.raises(Exception,match='dimensions'): m.Run()
        assert not m.state
        m.set(sizeHistN=4);m.Run();assert m.state
    finally: m.struct_cleanup()


@pytest.mark.parametrize('changes',[dict(sizeHistN=46341),dict(mChebyshev=2147483647),
                                     dict(sizeHistN=40000,mChebyshev=20)])
def test_dimension_errors_precede_allocation(tmp_path,changes):
    m=model(tmp_path,**changes)
    try:
        with pytest.raises(Exception,match='dimensions'):m.Run()
        assert not m.state and not m.getAllocationInfo()['live']
    finally:m.struct_cleanup()


def test_common_plan_budget_and_recovery(tmp_path,monkeypatch):
    m=model(tmp_path,sizeHistN=128,mChebyshev=15)
    monkeypatch.setenv('CBALLS_MEMORY_BUDGET_MB','1')
    try:
        with pytest.raises(Exception,match='complete common histogram plan'):m.Run()
        assert not m.state and not m.getAllocationInfo()['live']
        monkeypatch.setenv('CBALLS_MEMORY_BUDGET_MB','32')
        m.Run(level=['StartRun_Common'])
        assert m.getAllocationInfo()['scalar_tensor_allocated']
    finally:m.struct_cleanup()


def test_sincos_does_not_silently_ignore_pair_only(tmp_path):
    m=model(tmp_path,'octree-sincos-omp',options='only-2pcf,no-out-Hist')
    try:
        with pytest.raises(Exception,match='does not implement only-2pcf'):m.Run()
    finally:m.struct_cleanup()

@pytest.mark.parametrize('engine',['kdtree-2balls-omp','balltree-2balls-omp'])
def test_weighted_smoothing_keeps_constant_scalar(tmp_path,engine):
    m=model(tmp_path,engine,rsmooth='0.0003',options='smooth-pivot,no-out-Hist,weights-norm,no-one-ball,no-two-balls')
    rng=np.random.default_rng(122)
    centers=rng.normal(size=(12,3));centers/=np.linalg.norm(centers,axis=1)[:,None]
    p=np.repeat(centers,2,axis=0)+rng.normal(scale=1e-5,size=(24,3));p/=np.linalg.norm(p,axis=1)[:,None]
    m.set_catalog(p,kappa=np.ones(24),weights=np.linspace(.5,2.,24))
    try:
        m.Run()
        xi=m.getHistXi2pcf(); nonempty=m.getHistNN()>0
        assert np.count_nonzero(nonempty)>0
        np.testing.assert_allclose(xi[nonempty],1.,atol=2e-12,rtol=2e-12)
    finally:m.struct_cleanup()


def test_unused_forest_axes_do_not_overflow_metadata(tmp_path):
    m=cballs()
    m.set(searchMethod='lya-1d-2pcf-omp',rootDir=str(tmp_path),numberThreads=2,
          lya2RpBins=4,lya2RpMax=20.,lya3RBins=2147483647,
          options='no-out-Hist',verbose=0,verbose_log=0)
    m.set_forest_catalog(np.array([[0.,0.,10.],[1.,0.,11.],[0.,1.,12.]]),
                         delta=np.array([.1,.2,-.1]),weights=np.ones(3),
                         forest_ids=np.array([1,2,3],dtype=np.int64))
    try:
        m.Run()
        edges=m.getRunMetadata()['bin_edges']
        assert set(edges)=={'two_point_abs_parallel'}
        assert not m.getAllocationInfo()['scalar_tensor_allocated']
    finally:m.struct_cleanup()


def test_legacy_balltree_pair_storage(tmp_path):
    m=model(tmp_path,'balltree-2balls-omp',options='legacy-one-ball,only-2pcf,no-out-Hist,no-smooth-pivot,no-one-ball')
    try:
        m.Run()
        assert m.state and np.isfinite(m.getHistXi2pcf()).all()
        assert not m.getAllocationInfo()['scalar_tensor_allocated']
    finally:m.struct_cleanup()
