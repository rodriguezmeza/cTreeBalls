"""Independent physical survey product, bandwidth, support and ownership oracles."""
from pathlib import Path
import copy
import json
import sys
import tempfile
import numpy as np
import pytest
from numpy.polynomial.legendre import leggauss, legvander, legval
from cyballs import cballs, solve_physical_survey, qualification_report, save_result_packet, load_result_packet


def multiply_window(signal, window, solve_lmax):
    """Gauss quadrature integrates these polynomial products exactly."""
    x,w=leggauss(96)
    product=legval(x,signal)*legval(x,window)
    return ((legvander(x,solve_lmax).T @ (w*product))/2)*(2*np.arange(solve_lmax+1)+1)


def test_positive_window_counterexample():
    window=np.array([1.,.8,.8])
    x=np.linspace(-1,1,1001)
    assert legval(x,window).min()>.4666
    n=multiply_window([1.,1.],window,1)
    np.testing.assert_allclose(n,[1.2666666666666666,2.12],rtol=1e-14)
    old_matrix=np.array([[1.,.8/3],[.8,1.]])
    assert abs(np.linalg.solve(old_matrix,n)[0]-1)>.1
    result=solve_physical_survey(n,window,signal_lmax=1,padding=0)
    np.testing.assert_allclose(result['zeta_corrected_legendre'],[1.,1.],rtol=1e-14)
    assert result['zeta_status']==1 and result['zeta_relative_residual']<1e-14
    expected=1/np.linalg.cond(np.array([[1.,.8/3],[.8,1.32]]),np.inf)
    np.testing.assert_allclose(result['zeta_rcond'],expected,rtol=1e-14)
    with pytest.raises(ValueError,match='bandwidth'):
        solve_physical_survey(n,window[:2],signal_lmax=1,padding=0)


@pytest.mark.parametrize('scale',[1e-250,1.,1e250])
@pytest.mark.parametrize('padding',[0,1,3])
def test_padded_polynomial_recovery(scale,padding):
    signal=np.array([.25,-.12,.08]);degree=2+padding
    window=np.zeros(2*degree+1);window[:3]=[1.,.3,.4]
    n=multiply_window(signal,window,degree)
    out=solve_physical_survey(n*scale,window*scale,signal_lmax=2,padding=padding)
    assert out['zeta_status']==1
    np.testing.assert_allclose(out['zeta_corrected_legendre'],signal,rtol=2e-12,atol=1e-14)


def test_padding_recovers_leakage_from_omitted_signal():
    window=np.array([1.,.6,.3,0.,0.,0.,0.]);signal=[.2,-.1,.5,.4]
    n=multiply_window(signal,window,3)
    low=solve_physical_survey(n[:2],window,signal_lmax=1,padding=0)
    full=solve_physical_survey(n,window,signal_lmax=1,padding=2)
    assert np.linalg.norm(low['zeta_corrected_legendre']-signal[:2])>.1
    np.testing.assert_allclose(full['zeta_corrected_legendre'],signal[:2],atol=2e-14)


@pytest.mark.parametrize('window,status',[(np.zeros(3),0),(np.array([1.,3.,5.]),2),
    (np.array([1.,0.,(-1+1e-8)*2.5]),3),(np.array([1.,np.nan,0.]),5)])
def test_invalid_support(window,status):
    out=solve_physical_survey([1.,1.],window,signal_lmax=1,padding=0,rcond_min=1e-6)
    assert out['zeta_status']==status and not out['zeta_valid']
    assert np.isnan(out['zeta_corrected_legendre']).all()


def catalogs(n=48):
    rng=np.random.default_rng(61011)
    random=rng.uniform(0,1,(n,3));random[:,0]**=2  # anisotropic selection/window
    data=rng.uniform(.1,.9,(27,3))
    return data,random


def run_survey(root,*,output=False,threads=1,engine='octree-3pcf-3d-omp',n=48,model=None,**controls):
    m=model or cballs();data,random=catalogs(n)
    p=dict(searchMethod=engine,rootDir=str(root),sizeHistN=2,mChebyshev=3,
        rangeN=1.2,rminHist=.05,useLogHist=False,usePeriodic=False,numberThreads=threads,
        verbose=0,verbose_log=0,iCatalogs='1,2',surveySignalLmax=1,surveySolvePadding=1,
        options='survey-estimator-3d,compute-2pcf-3d,compute-3pcf-3d'+('' if output else ',no-out-Hist'))
    p.update(controls);m.set(p)
    m.set_catalog(data,weights=np.ones(len(data)),catalog=0)
    m.set_catalog(random,weights=np.ones(len(random)),catalog=1)
    try:
        m.Run()
        if '-mpi' in engine and m.getRunMetadata()['parallel']['rank']!=0:return None
        return m.getResults()
    finally:
        m.struct_cleanup()


def direct_injected(random,alpha,degree,signal):
    """Catalog tuples inject zeta(mu) into the measured random-window measure.

    This isolates correction recovery from stochastic D-R sampling noise. It is
    a finite-catalog known-signal injection, not a claim of unbiased cosmic data.
    No Gaunt coefficients or native moments are used to form the numerator.
    """
    out=np.zeros((degree+1,2,2))
    for i,p in enumerate(random):
        dr=random-p;radius=np.linalg.norm(dr,axis=1)
        keep=(radius>.05)&(radius<1.2)
        ids=np.flatnonzero(keep);dirs=dr[keep]/radius[keep,None]
        bins=((radius[keep]-.05)*2/(1.2-.05)).astype(int)
        mu=np.clip(dirs@dirs.T,-1,1)
        valid=ids[:,None]!=ids[None,:]
        grid=legvander(mu,degree);injected=legval(mu,signal)
        for a in range(2):
            for b in range(2):
                mask=valid & (bins[:,None]==a)&(bins[None,:]==b)
                out[:,a,b]+=np.sum(grid[mask]*injected[mask,None],axis=0)*(2*np.arange(degree+1)+1)*alpha**3
    return out


@pytest.mark.parametrize('density',[24,48,96])
@pytest.mark.parametrize('window_lmax',[4,6])
def test_catalog_injection_independent_density_and_bandwidth(tmp_path,density,window_lmax):
    packet=run_survey(tmp_path,n=density,surveyWindowLmax=window_lmax)
    data,random=catalogs(density);signal=np.array([.2,-.07,.04])
    # Full solve includes ell=2 padding; retain the first two coefficients.
    injected=direct_injected(random,len(data)/len(random),2,signal)
    window=packet['arrays']['random_triple_numerator']
    out=solve_physical_survey(injected,window,signal_lmax=1,padding=1)
    assert out['zeta_valid'].all()
    np.testing.assert_allclose(out['zeta_corrected_legendre'],np.broadcast_to(signal[:2,None,None],(2,2,2)),atol=3e-12)
    assert out['zeta_relative_residual'].max()<1e-14


def test_native_owned_output_and_packet_semantics(tmp_path):
    off=run_survey(tmp_path/'off')
    on=run_survey(tmp_path/'on',output=True)
    threaded=run_survey(tmp_path/'threaded',threads=3)
    assert len(off['arrays'])==17
    assert off['arrays']['data_triple_numerator'].shape==(3,2,2)
    assert off['arrays']['random_triple_numerator'].shape==(5,2,2)
    for name,value in off['arrays'].items():
        np.testing.assert_array_equal(value,on['arrays'][name])
        np.testing.assert_allclose(value,threaded['arrays'][name],rtol=3e-12,atol=1e-12,equal_nan=True)
    assert not list((tmp_path/'off').glob('hist*.txt'))
    rows=np.loadtxt(tmp_path/'on/histZetaM_3d_survey.txt')
    for row in rows:
        ell,a,b=map(int,row[:3]);a-=1;b-=1
        assert row[5]==on['arrays']['zeta_corrected_legendre'][ell,a,b]
        assert row[6]==on['arrays']['zeta_corrected_encore'][ell,a,b]
        assert row[14]==on['arrays']['zeta_rcond'][a,b]
        assert row[15]==on['arrays']['zeta_relative_residual'][a,b]
    independent=solve_physical_survey(on['arrays']['data_triple_numerator'],on['arrays']['random_triple_numerator'],signal_lmax=1,padding=1)
    np.testing.assert_allclose(independent['zeta_corrected_legendre'],on['arrays']['zeta_corrected_legendre'],rtol=1e-14)
    save_result_packet(off,tmp_path/'packet');saved=load_result_packet(tmp_path/'packet')
    assert qualification_report(saved,on,rtol=1e-11)['status']=='QUALIFIED'
    changed=copy.deepcopy(saved);changed['arrays']['zeta_corrected_legendre'][0,0,1]+=.1
    q=qualification_report(changed,on,rtol=1e-11)
    assert q['observables']['zeta_corrected_legendre']['status']=='REJECTED'
    assert q['status']=='REJECTED'
    del changed['arrays']['zeta_corrected_legendre']
    with pytest.raises(ValueError):qualification_report(changed,changed)


@pytest.mark.parametrize('controls',[
    dict(surveyWindowLmax=3),dict(surveySignalLmax=-2),dict(surveySolvePadding=-1),
    dict(surveySignalLmax=32,surveySolvePadding=1),dict(surveyWindowLmax=65),
    dict(surveyRcondMin=0),dict(surveyRcondMin=1.1),dict(surveyRcondMin=float('nan'))])
def test_bad_controls_and_recovery(tmp_path,controls):
    model=cballs()
    with pytest.raises(Exception):run_survey(tmp_path/'bad',model=model,**controls)
    good=run_survey(tmp_path/'good',model=model,surveySignalLmax=1,surveySolvePadding=1,surveyWindowLmax=-1,surveyRcondMin=1e-10)
    assert good['arrays']['zeta_valid'].all()


def test_condition_policy_is_qualification_scope(tmp_path):
    a=run_survey(tmp_path/'a');b=run_survey(tmp_path/'b',surveyRcondMin=.9)
    assert not b['arrays']['zeta_valid'].any()
    assert np.isnan(b['arrays']['zeta_corrected_legendre']).all()
    with pytest.raises(ValueError,match='matching'):qualification_report(a,b)
    assert qualification_report(b,b)['status']=='REJECTED'


def test_selected_pointwise_science_metric(tmp_path):
    a=run_survey(tmp_path);b=copy.deepcopy(a)
    name='zeta_corrected_legendre';b['arrays'][name][0,0,0]+=.01
    select={name:np.ones(a['arrays'][name].shape,dtype=bool)};select[name][0,0,0]=False
    assert qualification_report(b,a,rtol=1e-12,metric='pointwise')['status']=='REJECTED'
    report=qualification_report(b,a,rtol=1e-12,metric='pointwise',selection=select)
    assert report['status']=='QUALIFIED' and report['selection'][name]['selected_bins']==7
    with pytest.raises(ValueError):qualification_report(a,a,selection={name:np.zeros((2,2,2),dtype=bool)})



def test_reused_object_preserves_catalog_identity(tmp_path):
    model=cballs()
    a=run_survey(tmp_path/'first',model=model)
    b=run_survey(tmp_path/'second',model=model,output=True)
    assert qualification_report(a,b,rtol=1e-12)['status']=='QUALIFIED'


def test_default_mapping_and_high_window_moments(tmp_path):
    from numpy.polynomial.legendre import legvander
    out=run_survey(tmp_path,surveySignalLmax=-1,surveySolvePadding=1,surveyWindowLmax=64)
    config=out['metadata']['physical_survey']
    assert config['signal_lmax']==2 and config['solve_lmax']==3 and config['window_lmax']==64
    data,random=catalogs()
    direct=direct_injected(random,len(data)/len(random),64,[1.])
    np.testing.assert_allclose(out['arrays']['random_triple_numerator'],direct,rtol=2e-10,atol=3e-9)


def test_monopole_only_and_empty_bins(tmp_path):
    mono=run_survey(tmp_path/'mono',surveySignalLmax=0,surveySolvePadding=0,surveyWindowLmax=0)
    arrays=mono['arrays'];assert arrays['zeta_valid'].all()
    np.testing.assert_allclose(arrays['zeta_corrected_legendre'],arrays['data_triple_numerator']/arrays['random_triple_numerator'],rtol=2e-14)
    empty=run_survey(tmp_path/'empty',rangeN=10.,rminHist=5.)['arrays']
    assert not empty['zeta_valid'].any() and not empty['xi_valid'].any()
    assert np.isnan(empty['zeta_corrected_legendre']).all() and np.isnan(empty['xi_corrected']).all()
    assert (empty['zeta_status']==0).all() and (empty['xi_status']==0).all()


def test_legacy_raw_only_survey_packets_require_rerun(tmp_path):
    packet=run_survey(tmp_path)
    packet['metadata'].pop('physical_survey')
    packet['arrays']={k:v for k,v in packet['arrays'].items() if k.startswith(('data_','random_'))}
    with pytest.raises(ValueError,match='legacy physical survey'):
        qualification_report(packet,packet)


def test_legendre_solver_rejects_complex_inputs():
    with pytest.raises(ValueError,match='real'):
        solve_physical_survey([1+1j,0],[1,0,0],signal_lmax=1,padding=0)


def mpi_main():
    from mpi4py import MPI
    comm=MPI.COMM_WORLD
    root=comm.bcast(tempfile.mkdtemp(prefix='physical-survey-mpi-') if comm.rank==0 else None,root=0)
    model=cballs()
    off=run_survey(Path(root)/'off',engine='octree-3pcf-3d-mpi',model=model)
    on=run_survey(Path(root)/'on',engine='octree-3pcf-3d-mpi',output=True,threads=2,model=model)
    serial=run_survey(Path(root)/('serial-'+str(comm.rank)))
    if comm.rank==0:
        for name in off['arrays']:np.testing.assert_allclose(off['arrays'][name],on['arrays'][name],rtol=1e-12,atol=1e-12,equal_nan=True)
        for name in off['arrays']:np.testing.assert_allclose(off['arrays'][name],serial['arrays'][name],rtol=2e-11,atol=1e-12,equal_nan=True)
        assert qualification_report(off,on,rtol=1e-11)['status']=='QUALIFIED'
    # Different dimensions must fail collectively before any mismatched reduction.
    try:
        run_survey(Path(root)/'mismatch',engine='octree-3pcf-3d-mpi',model=model,surveyWindowLmax=4+2*(comm.rank%2))
        if comm.size>1:raise AssertionError('mismatched survey controls accepted')
    except Exception as exc:
        if comm.size==1:raise
        assert 'controls' in str(exc)
    run_survey(Path(root)/'retry',engine='octree-3pcf-3d-mpi',model=model,surveyWindowLmax=4)
    comm.Barrier()
    if comm.rank==0:print('PASS: physical survey ownership, corrected products, controls and recovery at',comm.size,'ranks')

if __name__=='__main__':
    if '--mpi' in sys.argv:mpi_main()
    else:raise SystemExit(pytest.main([__file__,'-q']))
