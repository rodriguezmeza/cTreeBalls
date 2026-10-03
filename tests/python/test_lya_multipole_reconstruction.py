"""Exact angular companions must not inherit finite-L reconstruction error."""
import json
import os
import re
import subprocess
import numpy as np
import pytest
import test_lya_forest_los_tree as los
from test_lya_triplet_acceleration import MULTIPOLE, moment_oracle


def triangle_bins(points,radius=30.,rbins=4,tbins=5,mbins=4):
    """Independent individual ordered-triangle oracle, with no harmonic sums."""
    sums=np.zeros((rbins,rbins,tbins,tbins,mbins,2));count=0
    for p in points:
        direction=p[:3]/np.linalg.norm(p[:3]);neighbors=[]
        for q in points:
            v=q[:3]-p[:3];r=np.linalg.norm(v)
            if p[5]==q[5] or not 0<r<radius:continue
            b=int(r/radius*rbins)
            t=min(int(np.arccos(np.clip(v@direction/r,-1,1))/np.pi*tbins),tbins-1)
            neighbors.append((q,v,r,b,t))
        for q,u,ru,b,t in neighbors:
            for s,v,rv,c,h in neighbors:
                if q[5]==s[5]:continue
                m=min(int((np.clip(u@v/(ru*rv),-1,1)+1)*.5*mbins),mbins-1)
                d=p[4]*q[4]*s[4];n=d*p[3]*q[3]*s[3]
                if d>0:sums[b,c,t,h,m]+=(n,d)
                count+=1
    return sums,count


def run_native(tmp,points,*,lmax=4,threads=1,kernel=0,mode=1,mbins=4,extra=None,env=None):
    tmp.mkdir();catalog=los.catalog_file(tmp,points)
    params=dict(search=MULTIPOLE,infile=catalog,infileformat='lya-ascii',iCatalogs=1,
        rootDir=tmp/'out',numberThreads=threads,usePeriodic='false',useLogHist='false',
        rangeN=30,rminHist=.1,sizeHistN=4,lya3RMax=30,lya3RBins=4,lya3ThetaBins=5,
        lya3MuBins=mbins,lya3LMax=lmax,lya3Kernel=kernel,lya3MuMode=mode,
        verbose=2,verbose_log=0,options='no-smooth-pivot,lya-output-empty-bins')
    params.update(extra or {})
    proc=subprocess.run([str(los.BINARY),*(f'{k}={v}' for k,v in params.items())],
        capture_output=True,text=True,timeout=120,env=env)
    return proc,tmp/'out'


def require_success(proc):
    if proc.returncode:raise AssertionError(proc.stdout+proc.stderr)


def check_case(tmp,points,lmax=4,mbins=4):
    expected,count=triangle_bins(points,mbins=mbins)
    moments=moment_oracle(points,lmax)
    results=[]
    for kernel,threads in [(0,1),(0,4),(1,2),(2,3)]:
        proc,out=run_native(tmp/f'{kernel}-{threads}',points,lmax=lmax,threads=threads,kernel=kernel,mbins=mbins)
        require_success(proc)
        actual=np.loadtxt(out/'histZetaM_lya5d.txt')[:,-2:].reshape(expected.shape)
        raw=np.loadtxt(out/'histZetaM_lya_multipoles.txt')[:,-2:].reshape(moments.shape)
        for component in (0,1):
            scale=max(np.max(np.abs(expected[...,component])),1e-300)
            np.testing.assert_allclose(actual[...,component]/scale,expected[...,component]/scale,rtol=3e-11,atol=2e-13)
            scale=max(np.max(np.abs(moments[...,component])),1e-300)
            np.testing.assert_allclose(raw[...,component]/scale,moments[...,component]/scale,rtol=3e-11,atol=2e-13)
        np.testing.assert_allclose(actual[...,1],expected[...,1],rtol=3e-12,atol=1e-300)
        np.testing.assert_array_equal(actual[...,1]>0,expected[...,1]>0)
        np.testing.assert_allclose(actual.sum(axis=-2),raw[...,0,:],rtol=3e-11,atol=2e-13*max(1.,np.max(np.abs(expected))))
        header=(out/'histZetaM_lya5d.txt').read_text()
        if 'EXACT hard mu bins' not in header:raise AssertionError(header[:1000])
        if int(re.search(r'ordered triplets: (\d+)',header)[1])!=count:raise AssertionError('double-counted or lost triangles')
        metadata=json.loads((out/'run-metadata.json').read_text())
        if metadata['lya_mu_output']!={'mode':1,'exact_bins_available':True,'exact_product':'_lya5d','finite_l_product':'_lya5d_multipole','discovery_shared':True}:
            raise AssertionError(metadata['lya_mu_output'])
        if not metadata['lya_3pcf']['mu_reconstruction_approximate']:raise AssertionError('finite-L diagnostic mislabeled')
        if 'APPROXIMATE' not in (out/'histZetaM_lya5d_multipole.txt').read_text():raise AssertionError('lost diagnostic qualification')
        results.append((actual,raw))
    np.testing.assert_array_equal(results[0][0],results[1][0])
    np.testing.assert_array_equal(results[0][1],results[1][1])


@pytest.mark.parametrize('case',['ordinary','bent','clustered','zero','one-forest','two-forests','dominant','signed','tiny','boundary'])
def test_exact_mu_triangle_contracts(tmp_path,case):
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
    check_case(tmp_path,points)


@pytest.mark.parametrize('lmax,mbins',[(0,1),(1,3),(12,20),(32,65)])
def test_exact_mu_orders_and_bin_counts(tmp_path,lmax,mbins):
    check_case(tmp_path,los.reference.POINTS.copy(),lmax,mbins)


def test_shared_discovery_and_exact_L_independence(tmp_path):
    points=los.make_forests(48,2,seed=222)
    points[:,:3]=[100,0,0]+(points[:,:3]-[150,0,0])*.08
    products=[];metas=[]
    for lmax,mode in [(0,0),(0,1),(8,1),(32,1)]:
        proc,out=run_native(tmp_path/f'{lmax}-{mode}',points,lmax=lmax,mode=mode,threads=3)
        require_success(proc);metas.append(json.loads((out/'run-metadata.json').read_text()))
        if mode:products.append(np.loadtxt(out/'histZetaM_lya5d.txt'))
        elif (out/'histZetaM_lya5d.txt').exists():raise AssertionError('default acquired extra product')
    for p in products[1:]:np.testing.assert_array_equal(p,products[0])
    if metas[0]['lya_los_tree']!=metas[1]['lya_los_tree']:raise AssertionError('duplicate discovery')
    if not metas[1]['lya_hierarchy']['certificates']:raise AssertionError('mixed-forest hierarchy untested')


def test_python_owned_products_mode_switch_and_recovery(tmp_path,monkeypatch):
    from cyballs import cballs
    points=los.reference.POINTS.copy();expected,_=triangle_bins(points)
    model=cballs();saved=[]
    try:
        for i,mode in enumerate((1,0,1,0,1)):
            model.set(dict(searchMethod=MULTIPOLE,rootDir=str(tmp_path/str(i)),
                numberThreads=3,rangeN=30.,rminHist=.1,sizeHistN=4,
                lya3RMax=30.,lya3RBins=4,lya3ThetaBins=5,lya3MuBins=4,
                lya3LMax=4,lya3MuMode=mode,lya3Kernel=0,verbose=0,verbose_log=0,
                options='no-smooth-pivot,no-out-Hist'))
            model.set_forest_catalog(points[:,:3],points[:,3],points[:,4],points[:,5].astype(np.int64))
            if i==0:
                monkeypatch.setenv('CBALLS_MEMORY_BUDGET_MB','.001')
                with pytest.raises(Exception,match='budget|memory|resource'):model.Run(level=['MainLoop'])
                monkeypatch.delenv('CBALLS_MEMORY_BUDGET_MB');model.struct_cleanup();continue
            model.Run(level=['MainLoop']);products=model.getForestResults();arrays=products['arrays']
            if ('triple_numerator' in arrays)!=bool(mode):raise AssertionError('stale or missing exact product')
            if mode:
                saved.append(arrays)
                np.testing.assert_allclose(arrays['triple_numerator'],expected[...,0],rtol=3e-12,atol=1e-12)
                np.testing.assert_allclose(arrays['triple_denominator'],expected[...,1],rtol=3e-12,atol=1e-12)
            model.struct_cleanup()
        for arrays in saved:np.testing.assert_allclose(arrays['triple_numerator'],expected[...,0],rtol=3e-12,atol=1e-12)
    finally:model.struct_cleanup()


@pytest.mark.parametrize('extra,reason',[
    ({'lya3MuMode':-1},'lya3MuMode=0 or 1'),({'lya3MuMode':2},'lya3MuMode=0 or 1'),
    ({'search':'lya-3pcf-omp'},'requires lya-anisotropic'),
    ({'lya3MuSlop':.1},'hard-bin 3PCF method'),({'lyaPivotRadius':.1},'pivot frontier/smoothing')])
def test_exact_mode_rejects_unsupported_controls(tmp_path,extra,reason):
    proc,_=run_native(tmp_path/'run',los.reference.POINTS.copy(),extra=extra)
    if proc.returncode==0 or reason not in proc.stdout+proc.stderr:raise AssertionError(proc.stdout+proc.stderr)


def test_additional_grid_budget_is_enforced(tmp_path):
    env=dict(os.environ,CBALLS_MEMORY_BUDGET_MB='2')
    for mode in (0,1):
        proc,_=run_native(tmp_path/str(mode),los.reference.POINTS.copy(),lmax=0,mode=mode,
            threads=4,mbins=65,env=env)
        if mode==0:require_success(proc)
        elif proc.returncode==0 or 'global and worker histogram plan' not in proc.stdout+proc.stderr:
            raise AssertionError(proc.stdout+proc.stderr)


def test_calibration_exact_mode_passes_while_finite_L_fails(tmp_path):
    # This is the same synthetic calibration that exposed L=4 reconstruction
    # errors; do not turn an exact companion pass into a finite-L claim.
    path=los.ROOT/'tests/python/benchmark_lya_triplet_kernels.py'
    proc=subprocess.run([os.sys.executable,str(path),'--synthetic','--outdir',str(tmp_path/'calibration'),
        '--mu-mode','exact','--lmax','4','--threads','2','--repeats','1','--warmups','0',
        '--require-accepted-multipole'],capture_output=True,text=True,timeout=120)
    require_success(proc)
    report=json.loads((tmp_path/'calibration/summary.json').read_text())
    exact=report['comparisons'][-1];finite=report['reconstruction_comparisons'][-1]
    if not exact['accepted'] or not exact['exact_raw_check_passed'] or finite['accepted']:
        raise AssertionError((exact,finite))
    if exact['invalid_eligible_bins'] or exact['empty_bin_reconstruction_l1']!=0:
        raise AssertionError(exact)


def test_finite_moments_do_not_determine_hard_bins():
    # Two positive discrete measures can share all moments through L, while
    # retaining different exact masses in an angular bin. No postprocessor
    # using only these moments can identify both hard-bin answers correctly.
    L=8;moments=[];masses=[]
    for n in (5,6):
        mu,weight=np.polynomial.legendre.leggauss(n)
        moments.append(weight@np.polynomial.legendre.legvander(mu,L))
        masses.append(weight[(mu>=0)&(mu<.5)].sum())
    np.testing.assert_allclose(moments[0],moments[1],rtol=0,atol=4e-15)
    if abs(masses[0]-masses[1])<.1:raise AssertionError(masses)


def test_python_resource_plan_includes_exact_grid():
    from cyballs import resource_plan
    params=dict(lya3RBins=4,lya3ThetaBins=5,lya3MuBins=65,lya3LMax=0)
    a=resource_plan(MULTIPOLE,24,threads=3,parameters=params)
    b=resource_plan(MULTIPOLE,24,threads=3,parameters=dict(params,lya3MuMode=1))
    sizes=a['policy']['type_bytes']
    expected=4**2*5**2*65*(2*sizes['real']+3*(2*sizes['real']+sizes['size_t']))
    if b['known_total_bytes']-a['known_total_bytes']!=expected:raise AssertionError((a,b))
    with pytest.raises(ValueError,match='lya3MuMode'):
        resource_plan(MULTIPOLE,24,parameters=dict(lya3MuMode=2))


def test_parameter_file_uses_exact_mode(tmp_path):
    cat=los.catalog_file(tmp_path,los.reference.POINTS.copy())
    par=tmp_path/'params.ini'
    par.write_text('\n'.join(f'{k}={v}' for k,v in dict(searchMethod=MULTIPOLE,
        infile=cat,infileformat='lya-ascii',iCatalogs=1,rootDir=tmp_path/'out',
        numberThreads=1,lya3MuMode=1,lya3LMax=0,lya3RBins=4,lya3ThetaBins=5,
        lya3MuBins=4,lya3RMax=30,sizeHistN=4,verbose=0,verbose_log=0,
        options='no-smooth-pivot,lya-output-empty-bins').items())+'\n')
    proc=subprocess.run([str(los.BINARY),str(par)],capture_output=True,text=True,timeout=30)
    require_success(proc)
    meta=json.loads((tmp_path/'out/run-metadata.json').read_text())
    if not meta['lya_mu_output']['exact_bins_available']:raise AssertionError(meta)
