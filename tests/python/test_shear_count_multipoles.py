"""N0-only shear: independent counts, full-window precedence, and owned output."""

import sys
from pathlib import Path

# Locate sibling analysis and benchmark helpers for direct script execution.
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python"))

import tempfile
from unittest.mock import patch

import numpy as np
import pytest
from cyballs import cballs

import test_shear_sphere_octree_omp as ref
import test_shear_pivot_reuse as reuse

TREES = ('octree', 'kdtree', 'balltree')
BASE = 'no-out-Hist,no-smooth-pivot,only-3pcf'
OFF = 'no-count-multipoles'


def measure(tree, catalog, options, nmax=7, threads=2, suffix='omp'):
    model = cballs()
    with tempfile.TemporaryDirectory(prefix='shear-n0-test-') as directory:
        model.set(dict(searchMethod=f'{tree}-shear-sphere-2balls-{suffix}',
            options=options, iCatalogs='1', usePeriodic=False, useLogHist=False,
            rminHist=ref.RMIN, rangeN=ref.RMAX, sizeHistN=ref.BINS,
            sizeHistPhi=ref.PHI_BINS, mChebyshev=nmax, theta=1., nsmooth=1,
            lengthBox=2.2, numberThreads=threads, verbose=0, verbose_log=0,
            rootDir=directory))
        p,g,w=catalog
        model.set_catalog(p,gamma1=g.real,gamma2=g.imag,weights=w)
        try:
            model.Run(level=['MainLoop'])
            result=model.getResults()
            result['orders']=model.getShearWindowMultipoleOrders()
            result['gamma']=model.getShearGammaXMultipoles()
            return result
        finally:
            model.struct_cleanup()


@pytest.mark.parametrize('tree', TREES)
@pytest.mark.parametrize('nmax', (2,7))
def test_monopole_matches_independent_triples(tree,nmax):
    catalog=ref.fixture(36)
    # Force exact traversal so changing approximation controls cannot hide an
    # incorrect count, diagonal subtraction, or normalization.
    actual=measure(tree,catalog,BASE+',no-one-ball,'+OFF,nmax=nmax)
    with patch.object(ref,'NMAX',nmax):
        oracle=ref.oracle(*catalog)
    arrays=actual['arrays']
    assert 'gamma_normalized' in arrays and 'gamma_corrected' not in arrays
    np.testing.assert_array_equal(actual['orders'],[0])
    np.testing.assert_allclose(arrays['window'],oracle['window'][2*nmax:2*nmax+1],rtol=2e-12,atol=1e-10)
    np.testing.assert_allclose(arrays['upsilon_raw'],oracle['upsilon'],rtol=2e-12,atol=2e-12)
    n0=oracle['window'][2*nmax]
    expected=np.zeros_like(oracle['upsilon'])
    np.divide(oracle['upsilon'],n0[None,None,:,:],out=expected,where=n0[None,None,:,:]>0)
    np.testing.assert_allclose(arrays['gamma_normalized'],expected,rtol=2e-12,atol=2e-12)
    assert actual['metadata']['shear_window']['max_order']==0
    assert actual['metadata']['shear_window']['edge_correction_applied'] is False


@pytest.mark.parametrize('tree',TREES)
def test_correction_overrides_opt_out(tree):
    catalog=ref.fixture(70)
    baseline=measure(tree,catalog,BASE,nmax=3)
    forced=measure(tree,catalog,BASE+','+OFF+',edge-corrections',nmax=3)
    for key,value in baseline['arrays'].items():
        np.testing.assert_array_equal(value,forced['arrays'][key])
    assert forced['metadata']['shear_window']['edge_correction_applied'] is True
    np.testing.assert_array_equal(forced['orders'],np.arange(-6,7))


@pytest.mark.parametrize('tree',TREES)
def test_reuse_monopole_and_thread_determinism(tree):
    catalog=reuse.clustered_catalog()
    with patch.dict('os.environ',CBALLS_SHEAR_PIVOT_TOL='0.1',CBALLS_SHEAR_BIN_THETA='0.1'):
        one=measure(tree,catalog,BASE+','+OFF+',shear-pivot-reuse',threads=1)
        many=measure(tree,catalog,BASE+','+OFF+',shear-pivot-reuse',threads=4)
    for key,value in one['arrays'].items():
        np.testing.assert_array_equal(value,many['arrays'][key])
    arrays=one['arrays'];n0=arrays['window'][0]
    expected=np.zeros_like(arrays['upsilon_raw'])
    np.divide(arrays['upsilon_raw'],n0[None,None,:,:],out=expected,where=n0[None,None,:,:]>0)
    np.testing.assert_allclose(arrays['gamma_normalized'],expected,rtol=2e-15,atol=1e-18)


@pytest.mark.parametrize('tree',TREES)
def test_pairs_unchanged(tree):
    catalog=ref.fixture(40)
    base='no-out-Hist,no-smooth-pivot,only-2pcf'
    a=ref.run_native(*catalog,1,options=base,engine=f'{tree}-shear-sphere-2balls-omp')
    b=ref.run_native(*catalog,1,options=base+','+OFF,engine=f'{tree}-shear-sphere-2balls-omp')
    for key in a:np.testing.assert_array_equal(a[key],b[key])


@pytest.mark.parametrize('tree',TREES)
def test_cross_catalog_monopole(tree):
    catalogs=[ref.fixture(25),ref.fixture(30),ref.fixture(32)]
    masks=[(np.arange(len(c[0]))%3!=0).astype(np.uint8) for c in catalogs]
    with patch.object(ref,'NMAX',7):
        for roles in ('1,2,3','1,2,1'):
            opts=BASE+',no-one-ball,read-mask'
            kw=dict(engine=f'{tree}-shear-sphere-2balls-omp',masks=masks,catalog_order=roles)
            full=ref.run_native_catalogs(catalogs,1,options=opts,**kw)
            off=ref.run_native_catalogs(catalogs,3,options=opts+','+OFF,**kw)
            np.testing.assert_allclose(off['window'][0],full['window'][14],rtol=2e-13,atol=1e-10)
            np.testing.assert_allclose(off['upsilon'],full['upsilon'],rtol=2e-13,atol=1e-12)


def test_driver_marks_monopole_output(tmp_path):
    from shear_corr_all_engines import RunConfig, ShearCatalog, run_ctreeballs
    from shear_products import attach_products
    p,g,w=ref.fixture(70)
    catalog=ShearCatalog(p,g.real,g.imag,w,geometry='sphere').normalized()
    config=RunConfig(output_dir=tmp_path,statistics='3pcf',bins=3,multipoles=7,
                     options=(OFF,'no-smooth-pivot','no-one-ball'))
    result=attach_products(run_ctreeballs(catalog,config,'balltree-shear-sphere-2balls-omp'),config,catalog)
    assert result['window'].shape[-1]==1
    assert result['window_correction_applied'] is False
    assert 'window_condition' not in result
    assert 'gamma_solve_relative_residual' not in result
    assert 'no edge correction' in result['normalizations']['gamma']


@pytest.mark.parametrize('tree',TREES)
def test_scalar_edge_correction_still_required(tree):
    # The generic opt-out cannot remove window modes from a scalar solve.
    from test_two_ball_edge_corrections import catalog
    p,k,w,mask=catalog(count=45)
    results=[]
    for extra in ('',','+OFF):
        model=cballs()
        with tempfile.TemporaryDirectory() as directory:
            model.set(dict(searchMethod=f'{tree}-2balls-omp',options='no-out-Hist,no-smooth-pivot,only-3pcf,no-two-balls,no-normalize-HistZeta,weights-norm,edge-corrections'+extra,
                iCatalogs='1',usePeriodic=False,useLogHist=False,rminHist=.02,rangeN=1.5,
                sizeHistN=4,mChebyshev=2,lengthBox=2.,numberThreads=2,verbose=0,verbose_log=0,rootDir=directory))
            model.set_catalog(p,kappa=k,weights=w)
            try:
                model.Run(level=['MainLoop'])
                results.append(model.getResults()['arrays']['zeta_corrected'])
            finally:model.struct_cleanup()
    np.testing.assert_array_equal(*results)


@pytest.mark.parametrize('tree',TREES)
def test_zero_weight_monopole_is_finite_and_empty(tree):
    p,g,w=ref.fixture(20)
    result=measure(tree,(p,g,np.zeros_like(w)),BASE+','+OFF)
    for key,value in result['arrays'].items():
        np.testing.assert_array_equal(value,np.zeros_like(value),err_msg=key)


@pytest.mark.parametrize('tree',TREES)
def test_smoothed_monopole_matches_full_window(tree):
    catalog=ref.smooth_fixture()
    opts='no-out-Hist,smooth-pivot,no-one-ball,only-3pcf'
    kwargs=dict(engine=f'{tree}-shear-sphere-2balls-omp',rsmooth=ref.SMOOTH_RMIN_ARCMIN)
    full=ref.run_native(*catalog,1,options=opts,**kwargs)
    n0=ref.run_native(*catalog,2,options=opts+','+OFF,**kwargs)
    np.testing.assert_allclose(n0['window'][0],full['window'][2*ref.NMAX],rtol=2e-13,atol=1e-10)
    np.testing.assert_allclose(n0['upsilon'],full['upsilon'],rtol=2e-13,atol=1e-12)
