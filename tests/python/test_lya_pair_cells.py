"""Pair-cell geometry, independent sums, reference fallback and shared 3PCF."""
import json
import re
import numpy as np
import pytest
import test_lya_forest_omp as oracle
import test_lya_forest_los_tree as los

METHODS=['lya-2pcf-omp','lya-los-tree-2pcf-omp']

@pytest.mark.parametrize('method',METHODS)
@pytest.mark.parametrize('case',['random','bent','antipodal','zero','one-forest','cancellation','boundaries'])
def test_pair_cells_against_pixel_reference(tmp_path,method,case):
    points=los.make_forests(9,31,seed=71)
    if case=='bent':points[::3,1:3]+=17
    elif case=='antipodal':points[::2,:3]*=-1
    elif case=='zero':points[::3,4]=0
    elif case=='one-forest':points[:,5]=101
    elif case=='cancellation':points[:,3]=np.where(np.arange(len(points))%2,1.,-1.)
    elif case=='boundaries':
        points=np.array([(100,0,0,.2,1,0),(106,0,0,-.1,2,1),
            (130,0,0,.3,0,2),(130-1e-10,0,1e-5,.4,1,3),
            (100,1e-5,0,.1,1,4),(-100,0,0,.2,1,5),
            (105,7,0,.6,1,6),(95,-7,0,.1,1,7)])
    cat=los.catalog_file(tmp_path,points)
    for name,threads,kernel in [('reference',1,0),('one',1,1),('many',4,1)]:
        los.run(cat,tmp_path/name,method,threads,extra=(f'lya2Kernel={kernel}',))
    los.compare_products(tmp_path/'reference',tmp_path/'one',(2,))
    assert (tmp_path/'one'/los.PRODUCTS[2]).read_bytes()==(tmp_path/'many'/los.PRODUCTS[2]).read_bytes()
    if case not in ('boundaries',):
        oracle.assert_histogram_close(oracle.read_2pcf(tmp_path/'one'/los.PRODUCTS[2]),
                                     oracle.oracle_2pcf(points),'independent pairs')

@pytest.mark.parametrize('kernel3',[0,1,2,3,4])
@pytest.mark.parametrize('prefix',['lya','lya-los-tree'])
def test_combined_shares_pair_and_preserves_triplets(tmp_path,kernel3,prefix):
    cat=los.catalog_file(tmp_path,los.make_forests(5,15))
    method=f'{prefix}-2pcf-3pcf-omp'
    los.run(cat,tmp_path/'ref',method,radius=60,extra=(f'lya3Kernel={kernel3}',))
    log=los.run(cat,tmp_path/'both',method,3,radius=60,extra=(f'lya3Kernel={kernel3}','lya2Kernel=1'))
    los.run(cat,tmp_path/'pair',f'{prefix}-2pcf-omp',2,extra=('lya2Kernel=1',))
    los.compare_products(tmp_path/'both',tmp_path/'ref')
    assert (tmp_path/'pair'/los.PRODUCTS[2]).read_bytes()==(tmp_path/'both'/los.PRODUCTS[2]).read_bytes()
    assert f'shared_with_3pcf={int(kernel3>=3)}' in log


def test_aggregation_reuse_pruning_and_slop_counts(tmp_path):
    points=los.make_forests(5,64)
    # Identical rays, narrow per-forest radial intervals, and a distant forest.
    for i,fid in enumerate(np.unique(points[:,5])):
        select=points[:,5]==fid
        points[select,:3]=np.column_stack((np.linspace(100+6*i,100.1+6*i,select.sum()),np.zeros((select.sum(),2))))
    points[points[:,5]==points[-1,5],:3]*=100
    cat=los.catalog_file(tmp_path,points)
    los.run(cat,tmp_path/'ref',METHODS[0])
    logs=[]
    for label,slop in [('exact',0),('approx',.1)]:
        logs.append(los.run(cat,tmp_path/label,METHODS[0],4,extra=('lya2Kernel=1',f'lya2RpSlop={slop}',f'lya2RtSlop={slop}')))
    los.compare_products(tmp_path/'ref',tmp_path/'exact',(2,))
    for key in ['aggregates','pruned','angle_reuses']:
        assert int(re.search(rf'\b{key}=(\d+)',logs[0])[1])>0
    for label in ['exact','approx']:
        p=tmp_path/label/los.PRODUCTS[2];r=tmp_path/'ref'/los.PRODUCTS[2]
        assert re.search(r'pairs: (\d+)',p.read_text())[1]==re.search(r'pairs: (\d+)',r.read_text())[1]
        data=np.loadtxt(p);ref=np.loadtxt(r)
        np.testing.assert_allclose(data[:,-2:].sum(axis=0),ref[:,-2:].sum(axis=0),rtol=2e-12,atol=2e-12)
    assert int(re.search(r'approximate_pairs=(\d+)',logs[1])[1])>0
    assert 'approximate geometry; exact forest exclusions/cutoffs' in logs[1]

@pytest.mark.parametrize('extra,reason',[
    (('lya2Kernel=2',),'lya2Kernel'),(('lya2Kernel=1','lya2RpSlop=nan'),'lya2RpSlop'),
    (('lya2Kernel=1','lya2RtSlop=-.1'),'pair slops'),(('lya2RpSlop=.1',),'requires lya2Kernel'),
    (('lya2Kernel=1','search=lya-3pcf-omp'),'pair cells require'),
    (('lya2Kernel=1','search=lya-1d-2pcf-omp'),'pair cells require')])
def test_validation(tmp_path,extra,reason):
    cat=los.catalog_file(tmp_path,oracle.POINTS)
    with pytest.raises(AssertionError,match=reason):los.run(cat,tmp_path/'bad',METHODS[0],extra=extra)


def test_pair_memory_budget(tmp_path,monkeypatch):
    cat=los.catalog_file(tmp_path,los.make_forests(8,250))
    monkeypatch.setenv('CBALLS_MEMORY_BUDGET_MB','1')
    with pytest.raises(AssertionError,match='persistent forest nodes and histograms'):
        los.run(cat,tmp_path/'budget',METHODS[0],extra=('lya2Kernel=1',))


def test_cython_pair_controls_repeated(tmp_path):
    cyballs=pytest.importorskip('cyballs');points=los.make_forests(4,15)
    balls=cyballs.cballs()
    try:
        for i,kernel in enumerate((1,0,1)):
            balls.set(dict(searchMethod='lya-2pcf-omp',rootDir=str(tmp_path/str(i)),numberThreads=3,
                usePeriodic=False,useLogHist=False,rangeN=30.,rminHist=.1,sizeHistN=4,
                lya2RpMax=30.,lya2RtMax=30.,lya2RpBins=5,lya2RtBins=6,lya2Kernel=kernel,
                lya2RpSlop=.01 if kernel else 0.,lya2RtSlop=0.,options='no-smooth-pivot',verbose=0,verbose_log=0))
            balls.set_forest_catalog(points[:,:3],points[:,3],points[:,4],points[:,5].astype(np.int64))
            balls.Run(level=['MainLoop'])
            meta=balls.getRunMetadata()
            assert meta['lya_2pcf']==dict(kernel=kernel,rp_slop=.01 if kernel else 0.,rt_slop=0.,approximate=bool(kernel))
            assert meta['lya_geometry']['pair_geometry_exact']==(kernel==0)
            balls.struct_cleanup()
    finally:balls.struct_cleanup()


@pytest.mark.parametrize('seed',range(6))
def test_variable_windows_and_slop_conserve_products(tmp_path,seed):
    points=los.make_forests(7,35,seed=seed)
    points[::4,1:3]+=seed*2
    cat=los.catalog_file(tmp_path,points)
    domain=dict(rp=10+seed*13,rt=8+seed*5)
    bins=(f'lya2RpBins={1+seed}',f'lya2RtBins={2+seed}')
    los.run(cat,tmp_path/'ref',METHODS[0],**domain,extra=bins)
    los.run(cat,tmp_path/'exact',METHODS[0],3,**domain,extra=bins+('lya2Kernel=1',))
    los.compare_products(tmp_path/'ref',tmp_path/'exact',(2,))
    los.run(cat,tmp_path/'approx',METHODS[0],3,**domain,extra=bins+('lya2Kernel=1','lya2RpSlop=.5','lya2RtSlop=.2'))
    paths=[tmp_path/name/los.PRODUCTS[2] for name in ['ref','approx']]
    assert re.search(r'pairs: (\d+)',paths[0].read_text())[1]==re.search(r'pairs: (\d+)',paths[1].read_text())[1]
    np.testing.assert_allclose(np.loadtxt(paths[0])[:,-2:].sum(axis=0),np.loadtxt(paths[1])[:,-2:].sum(axis=0),rtol=3e-12,atol=3e-10)


def test_wide_forest_ids(tmp_path):
    points=los.make_forests(4,23);ids=[2**62+1,2**62+2,-2**62-1,-2**62-2]
    mapping=dict(zip(np.unique(points[:,5]),ids));cat=tmp_path/'wide.txt'
    cat.write_text(''.join(' '.join(format(v,'.17g') for v in row[:5])+f' {mapping[row[5]]}\n' for row in points))
    los.run(cat,tmp_path/'ref',METHODS[0])
    los.run(cat,tmp_path/'cell',METHODS[0],extra=('lya2Kernel=1',))
    los.compare_products(tmp_path/'ref',tmp_path/'cell',(2,))
