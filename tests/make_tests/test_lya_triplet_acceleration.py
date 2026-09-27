"""Reference-versus-accelerated raw sums, counts, and anisotropic moments."""
import re
import numpy as np
import pytest
import test_lya_forest_los_tree as los


def test_typed_neighbor_sort_total_order(tmp_path):
    """Tied keys, negative IDs, reverse order, empty arrays and heap fallback."""
    import os, shlex, subprocess
    executable=tmp_path/'sort-check'
    subprocess.run(shlex.split(os.environ.get('CC','cc'))+['-std=c99','-O2',
        str(los.ROOT/'tests/test_lya_neighbor_sort.c'),'-o',str(executable)],check=True)
    subprocess.run([str(executable)],check=True,timeout=30)


def test_cosine_polar_lookup_keeps_reference_boundary_rounding(tmp_path):
    import os, shlex, subprocess
    executable=tmp_path/'geometry-check'
    subprocess.run(shlex.split(os.environ.get('CC','cc'))+['-std=c99','-O3','-fno-fast-math',
        str(los.ROOT/'tests/test_lya_neighbor_geometry.c'),'-o',str(executable),'-lm'],check=True)
    subprocess.run([str(executable)],check=True,timeout=30)


@pytest.mark.parametrize('prefix',['lya','lya-los-tree'])
def test_pair_owner_pruning_preserves_shuffled_zero_weight_products(tmp_path,prefix):
    """Independent pair oracle; the combined path cannot prune reverse legs."""
    points=los.make_forests(5,29,seed=92627)
    points[:,:3]*=.1
    points[::7,4]=0
    catalog=los.catalog_file(tmp_path,points)
    expected=los.reference.oracle_2pcf(points)
    for statistic in ['2pcf','2pcf-3pcf']:
        path=tmp_path/statistic
        los.run(catalog,path,f'{prefix}-{statistic}-omp',4)
        los.reference.assert_histogram_close(los.reference.read_2pcf(path/los.PRODUCTS[2]),expected,statistic)

@pytest.mark.parametrize('method', ['lya-3pcf-omp','lya-los-tree-3pcf-omp','lya-2pcf-3pcf-omp'])
@pytest.mark.parametrize('threads',[1,4])
def test_exact_kernels(tmp_path, method, threads):
    points=los.make_forests(9,32)
    points[::13,4]=0
    points[::7,3]*=100
    catalog=los.catalog_file(tmp_path,points)
    paths=[]
    logs=[]
    for kernel in (1,0,2):
        path=tmp_path/f'kernel-{kernel}'; paths.append(path)
        logs.append(los.run(catalog,path,method,threads,radius=160,extra=(f'lya3Kernel={kernel}',)))
    for path in paths[1:]:los.compare_products(paths[0],path,orders=(3,))
    assert int(re.search(r'aggregated_pairs=(\d+)',logs[1])[1])>0
    # A fixed block reduction must not depend on thread assignment.
    path=tmp_path/'other-threads'
    los.run(catalog,path,method,3,radius=160,extra=('lya3Kernel=0',))
    assert np.array_equal(np.loadtxt(paths[1]/los.PRODUCTS[3]),np.loadtxt(path/los.PRODUCTS[3]))

@pytest.mark.parametrize('value',[1e-110,1e100])
def test_extreme_weights_fallback(tmp_path,value):
    points=los.make_forests(4,12)
    points[:,4]=value
    catalog=los.catalog_file(tmp_path,points)
    # Only use finite triple products; underflow is permitted in reference.
    for kernel in (0,1):los.run(catalog,tmp_path/str(kernel),'lya-3pcf-omp',radius=160,extra=(f'lya3Kernel={kernel}','options=lya-output-empty-bins'))
    a=np.loadtxt(tmp_path/'0'/los.PRODUCTS[3],ndmin=2)
    b=np.loadtxt(tmp_path/'1'/los.PRODUCTS[3],ndmin=2)
    assert np.allclose(a,b,rtol=3e-12,atol=3e-12,equal_nan=True)

MULTIPOLE='lya-anisotropic-multipole-3pcf-omp'

def moment_oracle(points, lmax, radius=30., rbins=4, tbins=5):
    """Independent O(N^3) Legendre recurrence on individual ordered triplets."""
    out=np.zeros((rbins,rbins,tbins,tbins,lmax+1,2))
    for pi,p in enumerate(points):
        neighbors=[]
        los=p[:3]/np.linalg.norm(p[:3])
        for qi,q in enumerate(points):
            v=q[:3]-p[:3];r=np.linalg.norm(v)
            if pi==qi or q[5]==p[5] or not 0<r<radius:continue
            b=min(int(r/radius*rbins),rbins-1)
            t=min(int(np.arccos(np.clip(v@los/r,-1,1))/np.pi*tbins),tbins-1)
            neighbors.append((q,v/r,b,t))
        for qi,(q,u,b,t) in enumerate(neighbors):
            for r,v,c,s in neighbors[qi+1:]:
                if q[5]==r[5]:continue
                mu=np.clip(u@v,-1,1)
                leg=np.polynomial.legendre.legvander(mu,lmax).ravel()
                d=p[4]*q[4]*r[4];n=d*p[3]*q[3]*r[3]
                out[b,c,t,s,:,0]+=n*leg;out[b,c,t,s,:,1]+=d*leg
                out[c,b,s,t,:,0]+=n*leg;out[c,b,s,t,:,1]+=d*leg
    return out

def check_moment_output(folder,points,lmax,radius=30.,rbins=4,tbins=5):
    expected=moment_oracle(points,lmax,radius,rbins,tbins)
    table=np.loadtxt(folder/'histZetaM_lya_multipoles.txt',ndmin=2)
    actual=np.zeros_like(expected)
    for row in table:actual[tuple(row[:5].astype(int))]=row[-2:]
    np.testing.assert_allclose(actual,expected,rtol=2e-11,atol=3e-11)
    return actual

@pytest.mark.parametrize('lmax',[0,1,4,12,32])
def test_moments_and_window(tmp_path,lmax):
    points=los.reference.POINTS.copy()
    # Zero weight pixels retain geometric counts and cannot mark active bins.
    points[0,4]=0
    catalog=los.catalog_file(tmp_path,points)
    output=tmp_path/'moments'
    los.run(catalog,output,MULTIPOLE,4,extra=(f'lya3LMax={lmax}',))
    import json
    meta=json.loads((output/'run-metadata.json').read_text())
    assert meta['lya_3pcf']['mu_reconstruction_approximate']
    assert meta['lya_3pcf']['lmax']==meta['multipole_max']==lmax
    moment=check_moment_output(output,points,lmax)
    reconstructed=np.loadtxt(output/'histZetaM_lya5d_multipole.txt')
    for row in reconstructed:
        b,c,t,s,m=row[:5].astype(int);edges=np.linspace(-1,1,los.reference.MU_BINS+1)
        coefficients=[]
        for l in range(lmax+1):
            poly=np.polynomial.legendre.Legendre.basis(l).integ()
            coefficients.append((2*l+1)/2*(poly(edges[m+1])-poly(edges[m])))
        expected=np.asarray(coefficients)@moment[b,c,t,s]
        np.testing.assert_allclose(row[-2:],expected,rtol=2e-11,atol=3e-11)
        if row[-1]<=0:assert np.isnan(row[-3])
    # Window bins partition [-1,1], so their raw sum must equal the monopole.
    summed=reconstructed[:,-2:].reshape(4,4,5,5,los.reference.MU_BINS,2).sum(axis=-2)
    np.testing.assert_allclose(summed,moment[:,:,:,:,0,:],rtol=2e-11,atol=3e-11)


def test_multipole_two_forests_zero(tmp_path):
    points=los.make_forests(2,12)
    catalog=los.catalog_file(tmp_path,points)
    output=tmp_path/'two'
    los.run(catalog,output,MULTIPOLE,radius=160)
    assert np.all(np.loadtxt(output/'histZetaM_lya_multipoles.txt')[:,-2:]==0)


def test_multipole_threads(tmp_path):
    catalog=los.catalog_file(tmp_path,los.make_forests(9,16))
    arrays=[]
    for threads in (1,4):
        out=tmp_path/str(threads);los.run(catalog,out,MULTIPOLE,threads,radius=160)
        arrays.append(np.loadtxt(out/'histZetaM_lya_multipoles.txt'))
    assert np.array_equal(*arrays)


def test_multipole_budget(tmp_path, monkeypatch):
    import subprocess
    catalog=los.catalog_file(tmp_path,los.make_forests(5,10))
    # Catalog + common/global/worker histograms fit; multipole scratch does not.
    monkeypatch.setenv('CBALLS_MEMORY_BUDGET_MB','2')
    run=subprocess.run([str(los.BINARY),f'search={MULTIPOLE}',f'infile={catalog}',
        'infileformat=lya-ascii','iCatalogs=1',f'rootDir={tmp_path/"budget"}',
        'numberThreads=4','lya3RBins=4','lya3ThetaBins=4','lya3MuBins=4',
        'lya3LMax=32','lya3RMax=160','sizeHistN=4','options=no-smooth-pivot'],
        capture_output=True,text=True,timeout=30)
    assert run.returncode!=0
    assert 'multipole scratch' in run.stdout+run.stderr
    assert not (tmp_path/'budget/histZetaM_lya_multipoles.txt').exists()


def test_calibration_flags_invalid_bins_and_floor():
    import importlib.util
    from pathlib import Path
    path=Path(__file__).resolve().parents[1]/'python/benchmark_lya_triplet_kernels.py'
    spec=importlib.util.spec_from_file_location('triplet_calibration',path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    ref=np.zeros((3,13));ref[:,0]=np.arange(3)
    ref[:,-3:]=[[1,2,2],[1e-8,2e-8,2],[0,0,0]]
    cand=ref.copy();cand[0,-3:]=[np.nan,-1,-1]
    metrics=module.accuracy(ref,cand)
    assert metrics['eligible_bins']==1
    assert metrics['excluded_small_signal_bins']==1
    assert metrics['invalid_eligible_bins']==1
    assert not metrics['accepted']

@pytest.mark.parametrize('block',[1,8,64])
def test_explicit_scheduling_blocks(tmp_path,block):
    catalog=los.catalog_file(tmp_path,los.make_forests(5,20))
    a,b=tmp_path/'a',tmp_path/'b'
    los.run(catalog,a,'lya-los-tree-3pcf-omp',4,radius=160,extra=('lya3PivotBlock=0',))
    text=los.run(catalog,b,'lya-los-tree-3pcf-omp',4,radius=160,extra=(f'lya3PivotBlock={block}',))
    assert f'pivot_block={block}' in text
    los.compare_products(a,b,orders=(3,))


def test_large_mu_grid_fallback(tmp_path):
    catalog=los.catalog_file(tmp_path,los.make_forests(4,12))
    for kernel in (0,1):los.run(catalog,tmp_path/str(kernel),'lya-3pcf-omp',2,radius=160,
        extra=(f'lya3Kernel={kernel}','lya3MuBins=129'))
    los.compare_products(tmp_path/'0',tmp_path/'1',orders=(3,))


def test_multipole_rotation(tmp_path):
    points=los.reference.POINTS.copy()
    rotation,_=np.linalg.qr(np.random.default_rng(914).normal(size=(3,3)))
    arrays=[]
    for i in range(2):
        folder=tmp_path/str(i);folder.mkdir()
        catalog=los.catalog_file(folder,points)
        out=folder/'result';los.run(catalog,out,MULTIPOLE,2,extra=('lya3LMax=12',))
        arrays.append(np.loadtxt(out/'histZetaM_lya_multipoles.txt'))
        points[:,:3]=points[:,:3]@rotation
    np.testing.assert_allclose(*arrays,rtol=2e-11,atol=3e-11)
