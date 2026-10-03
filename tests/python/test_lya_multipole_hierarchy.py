"""Forest-disjoint Legendre hierarchy versus independent triangle recurrence."""
import json
import re
import numpy as np
import pytest
import test_lya_forest_los_tree as los
from test_lya_triplet_acceleration import MULTIPOLE, moment_oracle


def run_case(tmp_path,points,lmax=4,rbins=4,tbins=5,radius=30):
    catalog=los.catalog_file(tmp_path,points)
    expected=moment_oracle(points,lmax,radius,rbins,tbins)
    products=[]
    expected_count=int(round(moment_oracle(np.column_stack((points[:,:3],np.ones((len(points),2)),points[:,5])),0,radius,rbins,tbins)[...,1].sum()))
    for kernel,threads in [(1,1),(0,1),(2,1),(0,4),(2,3)]:
        out=tmp_path/f'{kernel}-{threads}'
        los.run(catalog,out,MULTIPOLE,threads,radius=radius,extra=(
            f'lya3Kernel={kernel}',f'lya3LMax={lmax}',f'lya3RBins={rbins}',f'lya3ThetaBins={tbins}'))
        text=(out/'histZetaM_lya_multipoles.txt').read_text()
        if int(re.search(r'ordered triplets: (\d+)',text)[1])!=expected_count:
            raise AssertionError('geometric triplet count changed')
        a=np.loadtxt(out/'histZetaM_lya_multipoles.txt')[:,-2:].reshape(expected.shape)
        # A global scale floor is used only for signed moments near zero;
        # populated, noncancelled bins also receive a relative comparison.
        for component in (0,1):
            scale=max(float(np.max(np.abs(expected[...,component]))),1e-300)
            np.testing.assert_allclose(a[...,component]/scale,expected[...,component]/scale,
                                       rtol=3e-11,atol=2e-13)
        # Monopole denominators are positive sums: require relative accuracy
        # even in bins tiny compared with a dominant forest elsewhere.
        np.testing.assert_allclose(a[...,0,1],expected[...,0,1],rtol=3e-12,atol=1e-300)
        np.testing.assert_array_equal(a[...,0,1]>0,expected[...,0,1]>0)
        reconstructed=np.loadtxt(out/'histZetaM_lya5d_multipole.txt')
        products.append((kernel,a,reconstructed))
        meta=json.loads((out/'run-metadata.json').read_text())
        if not meta['lya_3pcf']['mu_reconstruction_approximate']:
            raise AssertionError('lost reconstruction qualification')
    for kernel in (0,2):
        matches=[(a,b) for k,a,b in products if k==kernel]
        np.testing.assert_array_equal(matches[0][0],matches[1][0])
        np.testing.assert_array_equal(matches[0][1],matches[1][1])


@pytest.mark.parametrize('case',['ordinary','bent','clustered','zero','one-forest','two-forests','dominant','signed','tiny','boundary'])
def test_triangle_contracts(tmp_path,case):
    points=los.make_forests(6,4,seed=112358)
    points[:,:3]=[100,0,0]+(points[:,:3]-[150,0,0])*.1
    if case=='bent':points[::3,1:3]+=[.5,-1.7]
    if case=='clustered':points[:,:3]=[100,0,0]+(points[:,:3]-[100,0,0])*.001
    if case=='zero':points[::2,4]=0
    if case=='one-forest':points[:,5]=77
    if case=='two-forests':points[:,5]=np.arange(len(points))%2
    if case=='dominant':points[points[:,5]==points[0,5],4]=1e24
    if case=='signed':points[:,3]=np.resize([1e12,-1e12,1e-8,-1e-8],len(points))
    if case=='tiny':points[:,4]*=1e-90
    if case=='boundary':
        points=np.array([(100,0,0,.2,1,0),(107.5,0,0,-.4,1,1),
            (100,7.5,0,.7,0,2),(100,0,7.5,1.2,2,3),
            (130,0,0,.8,1,4),(130-1e-10,0,0,.3,1,5),
            (100,-15,0,-.2,1,6),(100,0,-22.5,.4,1,7)],dtype=float)
    run_case(tmp_path,points)


@pytest.mark.parametrize('lmax',[0,1,12,32])
def test_all_supported_orders(tmp_path,lmax):
    run_case(tmp_path,los.reference.POINTS.copy(),lmax=lmax)


def test_hierarchy_selected_and_rotation(tmp_path):
    points=los.make_forests(64,8,seed=919)
    cat=los.catalog_file(tmp_path,points)
    arrays=[]
    for k in (0,1,2):
        out=tmp_path/str(k)
        los.run(cat,out,MULTIPOLE,4,radius=160,extra=(f'lya3Kernel={k}','lya3RBins=8','lya3ThetaBins=8'))
        arrays.append(np.loadtxt(out/'histZetaM_lya_multipoles.txt')[:,-2:])
        if k==0:
            m=json.loads((out/'run-metadata.json').read_text())['lya_multipole_reuse']
            if not (m['hierarchical_pivots']>0 and m['products']<m['prefix_products']):
                raise AssertionError(m)
    for a in arrays[1:]:np.testing.assert_allclose(a,arrays[0],rtol=3e-11,atol=1e-10)
    # Use a small independent triangle oracle after a proper rigid rotation.
    p=los.reference.POINTS.copy()
    q,_=np.linalg.qr(np.random.default_rng(42).normal(size=(3,3)))
    p[:,:3]=p[:,:3]@q
    sub=tmp_path/'rotated';sub.mkdir();run_case(sub,p,lmax=12)


def test_owned_results_and_failure_recovery(tmp_path,monkeypatch):
    from cyballs import cballs
    points=los.reference.POINTS.copy()
    expected=moment_oracle(points,4)
    model=cballs()
    saved=[]
    try:
        for i,kernel in enumerate((2,0,1,0)):
            model.set(dict(searchMethod=MULTIPOLE,rootDir=str(tmp_path/str(i)),
                numberThreads=3,rangeN=30.,rminHist=.1,sizeHistN=4,
                lya3RMax=30.,lya3RBins=4,lya3ThetaBins=5,lya3MuBins=4,
                lya3LMax=4,lya3Kernel=kernel,verbose=0,verbose_log=0,
                options='no-smooth-pivot,no-out-Hist'))
            model.set_forest_catalog(points[:,:3],points[:,3],points[:,4],points[:,5].astype(np.int64))
            if i==0:
                monkeypatch.setenv('CBALLS_MEMORY_BUDGET_MB','.001')
                with pytest.raises(Exception,match='budget|memory|resource'):
                    model.Run(level=['MainLoop'])
                monkeypatch.delenv('CBALLS_MEMORY_BUDGET_MB')
                model.struct_cleanup()
                continue
            model.Run(level=['MainLoop'])
            arrays=model.getForestResults()['arrays']
            saved.append(arrays)
            np.testing.assert_allclose(arrays['moments_numerator'],expected[...,0],rtol=3e-11,atol=1e-10)
            np.testing.assert_allclose(arrays['moments_denominator'],expected[...,1],rtol=3e-11,atol=1e-10)
            model.struct_cleanup()
        for arrays in saved:
            np.testing.assert_allclose(arrays['moments_numerator'],expected[...,0],rtol=3e-11,atol=1e-10)
    finally:
        model.struct_cleanup()


def test_calibration_rejects_unknown_multipole_kernel(tmp_path):
    import importlib.util
    path=los.ROOT/'tests/python/benchmark_lya_triplet_kernels.py'
    spec=importlib.util.spec_from_file_location('moment_calibration',path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    with pytest.raises(SystemExit) as error:
        module.main(['--synthetic','--outdir',str(tmp_path),'--multipole-kernel','3'])
    if error.value.code!=2:raise AssertionError(error.value)
