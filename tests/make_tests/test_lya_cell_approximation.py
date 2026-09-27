"""Persistent forest/pivot cells against independent direct and Python oracles."""
import json
import re
import subprocess
import numpy as np
import pytest
import test_lya_forest_los_tree as los

METHODS=['lya-3pcf-omp','lya-los-tree-3pcf-omp','lya-2pcf-3pcf-omp','lya-los-tree-2pcf-3pcf-omp']

@pytest.mark.parametrize('kernel',[3,4])
@pytest.mark.parametrize('method',METHODS)
def test_cells_exact_and_threads(tmp_path,kernel,method):
    points=los.make_forests(6,20)
    points[::11,4]=0;points[::9,3]*=100
    catalog=los.catalog_file(tmp_path,points)
    ref=tmp_path/'ref';a=tmp_path/'a';b=tmp_path/'b'
    los.run(catalog,ref,method,radius=160,extra=('lya3Kernel=1',))
    log=los.run(catalog,a,method,1,radius=160,extra=(f'lya3Kernel={kernel}',))
    los.run(catalog,b,method,4,radius=160,extra=(f'lya3Kernel={kernel}',))
    orders=(2,3) if '2pcf-3pcf' in method else (3,)
    los.compare_products(ref,a,orders)
    for order in orders:
        assert np.array_equal(np.loadtxt(a/los.PRODUCTS[order]),np.loadtxt(b/los.PRODUCTS[order]))
    assert 'Persistent forest cells:' in log
    assert not json.loads((a/'run-metadata.json').read_text())['lya_geometry']['three_point_approximate']

@pytest.mark.parametrize('kernel',[3,4])
@pytest.mark.parametrize('case',['oracle','skew','boundary','tiny','large','two-forests'])
def test_cells_adversarial(tmp_path,kernel,case):
    points=los.reference.POINTS.copy() if case=='oracle' else los.make_forests(4,8)
    if case=='skew':points[:,:3]+=np.random.default_rng(17).normal(size=(len(points),3))*5
    if case=='boundary':
        points[:,:3]=np.array([100,0,0])+np.arange(len(points))[:,None]*np.array([10.,0.,0.])
    if case=='tiny':points[:,4]=1e-110
    if case=='large':points[:,4]=1e100
    if case=='two-forests':points[:,5]=np.arange(len(points))%2
    catalog=los.catalog_file(tmp_path,points)
    for k in [1,kernel]:
        los.run(catalog,tmp_path/str(k),'lya-3pcf-omp',radius=30,extra=(f'lya3Kernel={k}','options=lya-output-empty-bins'))
    los.compare_products(tmp_path/'1',tmp_path/str(kernel),orders=(3,))
    if case=='oracle':
        actual=los.reference.read_3pcf(tmp_path/str(kernel)/los.PRODUCTS[3])
        los.reference.assert_histogram_close({key:value for key,value in actual.items() if value[1]>0},los.reference.oracle_3pcf(),'cell oracle')

@pytest.mark.parametrize('kernel',[0,3,4])
def test_mu_only_conserves_leg_marginals(tmp_path,kernel):
    catalog=los.catalog_file(tmp_path,los.make_forests(8,24))
    tables=[]
    for label,slop in [('exact',0),('slop',.25)]:
        path=tmp_path/label
        log=los.run(catalog,path,'lya-2pcf-3pcf-omp',4,radius=160,extra=(f'lya3Kernel={kernel}',f'lya3MuSlop={slop}','options=lya-output-empty-bins'))
        tables.append(np.loadtxt(path/los.PRODUCTS[3]))
    assert int(re.search(r'approximate_pairs=(\d+)',log)[1])>0
    a,b=tables
    assert np.array_equal(a[:,:5],b[:,:5])
    np.testing.assert_allclose(a[:,-2:].reshape(-1,los.reference.MU_BINS,2).sum(axis=1),b[:,-2:].reshape(-1,los.reference.MU_BINS,2).sum(axis=1),rtol=1e-11,atol=1e-8)
    los.compare_products(tmp_path/'exact',tmp_path/'slop',orders=(2,))
    meta=json.loads((tmp_path/'slop'/'run-metadata.json').read_text())['lya_geometry']
    assert meta['three_point_approximate'] and meta['pair_geometry_exact']

@pytest.mark.parametrize('radial,polar',[(0,.3),(.3,0),(.3,.3)])
def test_independent_slops_and_cutoffs(tmp_path,radial,polar):
    catalog=los.catalog_file(tmp_path,los.make_forests(5,16))
    tables=[];logs=[]
    for label,kernel,threads in [('ref',1,1),('one',4,1),('many',4,4)]:
        extra=(f'lya3Kernel={kernel}','options=lya-output-empty-bins')
        if kernel==4:extra+=(f'lya3MuSlop=.2',f'lya3RadialSlop={radial}',f'lya3PolarSlop={polar}')
        log=los.run(catalog,tmp_path/label,'lya-3pcf-omp',threads,radius=73,extra=extra)
        logs.append((tmp_path/label/los.PRODUCTS[3]).read_text())
        tables.append(np.loadtxt(tmp_path/label/los.PRODUCTS[3]))
    # Geometry slop may move bins, but no triplet may enter/leave the hard domain.
    counts=[re.search(r'ordered triplets: (\d+)',s)[1] for s in logs]
    assert len(set(counts))==1
    np.testing.assert_allclose(tables[0][:,-2:].sum(axis=0),tables[1][:,-2:].sum(axis=0),rtol=1e-11,atol=1e-9)
    assert np.array_equal(tables[1],tables[2])

@pytest.mark.parametrize('extra',[
    ('lya3MuSlop=-.1',),('lya3MuSlop=nan',),('lya3RadialSlop=1.1',),
    ('lya3PolarSlop=.1','lya3Kernel=0'),('lya3MuSlop=.1','lya3Kernel=1'),
    ('lya3Kernel=4','lya3PivotCellMax=0'),
    ('lya3Kernel=4','search=lya-3pcf-mpi'),('lya3Kernel=4','search=lya-2pcf-omp')])
def test_invalid_controls(tmp_path,extra):
    catalog=los.catalog_file(tmp_path,los.reference.POINTS)
    with pytest.raises(AssertionError):los.run(catalog,tmp_path/'bad','lya-3pcf-omp',extra=extra)


def test_pivot_aggregation_is_used(tmp_path):
    # Tight clusters inside three forests guarantee certified multi-pivot work.
    points=[]
    rng=np.random.default_rng(77)
    for forest,center in enumerate([[100,10,0],[130,30,4],[150,2,8]]):
        for i in range(16):points.append([*(np.asarray(center)+rng.normal(size=3)*.001),.5,1,forest])
    catalog=los.catalog_file(tmp_path,np.asarray(points))
    for kernel in [1,4]:
        log=los.run(catalog,tmp_path/str(kernel),'lya-3pcf-omp',radius=160,extra=(f'lya3Kernel={kernel}',))
    assert int(re.search(r'pivot_aggregates=(\d+)',log)[1])>0
    los.compare_products(tmp_path/'1',tmp_path/'4',orders=(3,))


def test_persistent_memory_preflight(tmp_path,monkeypatch):
    monkeypatch.setenv('CBALLS_MEMORY_BUDGET_MB','1')
    catalog=los.catalog_file(tmp_path,los.make_forests(100,30))
    with pytest.raises(AssertionError,match='persistent forest nodes'):
        los.run(catalog,tmp_path/'budget','lya-3pcf-omp',extra=('lya3Kernel=4','lya3RBins=1','lya3ThetaBins=1','lya3MuBins=1'))

@pytest.mark.parametrize('kernel',[3,4])
def test_cython_geometry_metadata(tmp_path,kernel):
    cyballs=pytest.importorskip('cyballs')
    points=los.make_forests(4,10)
    params=dict(searchMethod='lya-2pcf-3pcf-omp',rootDir=str(tmp_path/'memory'),
        numberThreads=2,usePeriodic=False,useLogHist=False,rangeN=160.,rminHist=.1,sizeHistN=4,
        lya3RMax=160.,lya3RBins=4,lya3ThetaBins=5,lya3MuBins=6,
        lya3Kernel=kernel,lya3MuSlop=.02,lya3RadialSlop=.01,lya3PolarSlop=.03,lya3PivotCellMax=4,
        options='no-smooth-pivot',verbose=0,verbose_log=0)
    balls=cyballs.cballs()
    try:
        balls.set(params)
        balls.set_forest_catalog(points[:,:3],points[:,3],points[:,4],points[:,5].astype(np.int64))
        balls.Run(level=['MainLoop'])
        meta=balls.getRunMetadata()
        assert meta['lya_3pcf']['kernel']==kernel
        assert meta['lya_geometry']==dict(mu_slop=.02,radial_slop=.01,polar_slop=.03,
            pivot_cell_max=4,three_point_approximate=True,pair_geometry_exact=True)
        balls.struct_cleanup()
        assert balls.getRunMetadata()==meta
    finally:balls.struct_cleanup()


def test_pruning_precedes_pixel_cache(tmp_path):
    # Well separated forest groups should never populate a full all-pairs cache.
    points=los.make_forests(8,24)
    ids=np.unique(points[:,5])
    for j,fid in enumerate(ids):
        points[points[:,5]==fid,:3]+=np.array([10000.*j,7000.*j,2000.*j])
    catalog=los.catalog_file(tmp_path,points)
    for kernel in [1,3]:
        log=los.run(catalog,tmp_path/str(kernel),'lya-3pcf-omp',radius=160,
                    extra=(f'lya3Kernel={kernel}',))
    los.compare_products(tmp_path/'1',tmp_path/'3',orders=(3,))
    assert int(re.search(r'pruned_nodes=(\d+)',log)[1])>0
    assert int(re.search(r'leaf_evaluations=(\d+)',log)[1])<len(points)


@pytest.mark.parametrize('cap',[2,8,64])
def test_pivot_child_cache_and_changing_los(tmp_path,cap):
    rng=np.random.default_rng(558)
    # Nearby multi-pixel pivots still span polar and mu boundaries. The forest
    # identity is not an assumption of collinearity with the observer.
    points=[]
    for fid,center in enumerate([[11.,3,4],[16,8,2],[21,2,7],[6,7,8],[18,3,9]]):
        for i in range(20):
            xyz=np.asarray(center)+rng.normal(size=3)*.7
            points.append([*xyz,rng.normal(),rng.uniform(.1,2),fid])
    catalog=los.catalog_file(tmp_path,np.asarray(points))
    for kernel in [1,4]:
        log=los.run(catalog,tmp_path/str(kernel),'lya-2pcf-3pcf-omp',radius=160,
            extra=(f'lya3Kernel={kernel}',f'lya3PivotCellMax={cap}',
                   'lya3ThetaBins=11','lya3MuBins=13'))
    los.compare_products(tmp_path/'1',tmp_path/'4')
    if cap>2:
        assert int(re.search(r'pair_cache_hits=(\d+)',log)[1])>0
    assert int(re.search(r'cache_hits=(\d+)',log)[1])>0
    assert int(re.search(r'leaf_evaluations=(\d+)',log)[1])<=len(points)**2


@pytest.mark.parametrize('geometry',['oblique','bent','closed','boundary','translated'])
@pytest.mark.parametrize('kernel',[3,4])
def test_forest_capsule_enclosures(tmp_path,geometry,kernel):
    rng=np.random.default_rng(129)
    points=[]
    for fid in range(5):
        t=np.linspace(-1,1,19)
        direction=np.array([1.,.8,.7]);center=np.array([100.,30.+fid*7,10.-fid*3])
        xyz=center+t[:,None]*direction*15
        if geometry=='bent':xyz[:,1]+=15*np.sin(t*3)
        if geometry=='closed':xyz=center+np.column_stack([np.cos(t*np.pi),np.sin(t*np.pi),t*0])*10
        if geometry=='boundary':xyz[:,0]=100+fid*10+np.nextafter(t*20,np.inf)
        if geometry=='translated':xyz+=np.array([1e8,-2e8,3e8])
        for v in xyz:points.append([*v,rng.normal(),rng.uniform(.1,2),fid])
    catalog=los.catalog_file(tmp_path,np.asarray(points))
    for k in [1,kernel]:
        los.run(catalog,tmp_path/str(k),'lya-3pcf-omp',radius=40,
            extra=(f'lya3Kernel={k}','lya3PivotCellMax=64','lya3RBins=7',
                   'lya3ThetaBins=9','lya3MuBins=11','options=no-check-two-bodies-eq-pos'))
    los.compare_products(tmp_path/'1',tmp_path/str(kernel),orders=(3,))


def test_oblique_capsule_prunes_overlapping_box(tmp_path):
    # The cluster is inside the line's Cartesian box, yet outside its R=20
    # tube. A box-only root rejection cannot establish this separation.
    line=np.array([[100.+t,t,t,.5,1,0] for t in np.linspace(-30,30,25)])
    cluster=np.array([[100.+t,25.,-25.,.2,1,1] for t in np.linspace(-.1,.1,8)])
    catalog=los.catalog_file(tmp_path,np.concatenate([line,cluster]))
    for kernel in [1,3]:
        log=los.run(catalog,tmp_path/str(kernel),'lya-3pcf-omp',radius=20,
                    extra=(f'lya3Kernel={kernel}',))
    los.compare_products(tmp_path/'1',tmp_path/'3',orders=(3,))
    assert int(re.search(r'pruned_nodes=(\d+)',log)[1])>0
    assert int(re.search(r'leaf_evaluations=(\d+)',log)[1])==0


def test_active_pivot_cache_memory_preflight(tmp_path,monkeypatch):
    monkeypatch.setenv('CBALLS_MEMORY_BUDGET_MB','12')
    # Large, tight pivot tasks make pixel-cache storage dominate the plan.
    points=np.array([[100.+20*fid+i*.001,5.+fid,2.,.2,1.,fid]
                     for fid in range(4) for i in range(250)])
    catalog=los.catalog_file(tmp_path,points)
    with pytest.raises(AssertionError,match='persistent nodes, histograms and geometry caches'):
        los.run(catalog,tmp_path/'budget','lya-3pcf-omp',threads=4,radius=160,
            extra=('lya3Kernel=4','lya3PivotCellMax=64','lya3RBins=1',
                   'lya3ThetaBins=1','lya3MuBins=1'))
    assert not (tmp_path/'budget'/los.PRODUCTS[3]).exists()


def test_large_pivot_extent_near_observer(tmp_path):
    # A small bounding-box center must not set the rounding scale for a pivot
    # cell whose actual extent is many orders of magnitude larger.
    rng=np.random.default_rng(471)
    points=[]
    for fid,scale in enumerate([1e8,1e3,2e3,3e3]):
        xyz=rng.normal(size=(16,3))*scale
        xyz=np.concatenate([xyz,-xyz])
        for v in xyz:points.append([*v,rng.normal(),rng.uniform(.2,2),fid])
    catalog=los.catalog_file(tmp_path,np.asarray(points))
    for kernel in [1,4]:
        los.run(catalog,tmp_path/str(kernel),'lya-3pcf-omp',radius=1e11,
                extra=(f'lya3Kernel={kernel}','lya3PivotCellMax=64'))
    los.compare_products(tmp_path/'1',tmp_path/'4',orders=(3,))
