"""Continuum distortion checked against independent four-index projection."""
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT), str(ROOT / 'tests' / 'python')]
from cyballs import lya_distortion_matrix


def fixture():
    rng = np.random.default_rng(913)
    ids = np.repeat(np.array([-2**59+1, 7, 2**58+3], dtype=np.int64), [4, 3, 5])
    radius = np.concatenate([100 + np.arange(n)*3.7 + k for k,n in enumerate([4,3,5])])
    angle = np.repeat([0., .075, .14], [4,3,5])
    pos = radius[:,None] * np.column_stack([np.cos(angle), np.sin(angle), angle*0])
    lam = 3.5 + (radius-100)*.0007
    w = rng.uniform(.3, 2., len(ids))
    p = rng.uniform(.4, 3., len(ids))
    return pos, w, ids, lam, p


def bins(pos, rpmax, rtmax, nrp, nrt):
    radius = np.sqrt(np.sum(pos*pos, axis=1))
    los = pos/radius[:,None]
    mu = np.clip(los@los.T, -1., 1.)
    rp = abs(radius[:,None]-radius[None,:])*np.sqrt((1+mu)/2)
    rt = (radius[:,None]+radius[None,:])*np.sqrt((1-mu)/2)
    result = (rp/rpmax*nrp).astype(int)*nrt+(rt/rtmax*nrt).astype(int)
    result[(rp >= rpmax) | (rt >= rtmax)] = -1
    return result


def oracle(pos, w, ids, lam, p, order, grid=(20.,40.,4,3), model=None):
    """Deliberately literal O(m²n²) reference, no low-rank expansion."""
    model = grid if model is None else model
    data_bin, model_bin = bins(pos,*grid), bins(pos,*model)
    na, nb = grid[2]*grid[3], model[2]*model[3]
    raw, den = np.zeros((na,nb)), np.zeros(na)
    forests = [np.flatnonzero(ids == f) for f in np.unique(ids)]
    projectors = []
    for indices in forests:
        q, l = p[indices], lam[indices]
        proj = np.eye(len(indices))
        if order >= 0:
            proj -= np.ones((len(indices),1))*q[None,:]/q.sum()
            x = l - np.average(l,weights=q)
            var = np.dot(q,x*x)
            if order == 1 and np.ptp(l[q>0]) > 0:
                proj -= x[:,None]*(q*x)[None,:]/var
        projectors.append(proj)
    for f, fi in enumerate(forests):
        for g in range(f+1,len(forests)):
            gj = forests[g]
            for ii,i in enumerate(fi):
                for jj,j in enumerate(gj):
                    a = data_bin[i,j]
                    if a < 0:
                        continue
                    den[a] += w[i]*w[j]
                    for kk,k in enumerate(fi):
                        for ll,l in enumerate(gj):
                            b = model_bin[k,l]
                            if b >= 0:
                                raw[a,b] += w[i]*w[j]*projectors[f][ii,kk]*projectors[g][jj,ll]
    matrix = np.full_like(raw,np.nan)
    np.divide(raw,den[:,None],out=matrix,where=den[:,None]>0)
    return raw,den,matrix


@pytest.mark.parametrize('order',[-1,0,1])
@pytest.mark.parametrize('wide_model',[False,True])
def test_four_index_oracle(order,wide_model):
    pos,w,ids,lam,p = fixture()
    w[[1,6,9]] = 0  # full projection domain includes pixels absent in data
    grid = (9.,25.,3,2)
    model = (30.,40.,5,3) if wide_model else grid
    result = lya_distortion_matrix(pos,w,ids,lam,projection_weights=p,projection_order=order,
        rp_max=grid[0],rt_max=grid[1],rp_bins=grid[2],rt_bins=grid[3],
        model_rp_max=model[0],model_rt_max=model[1],model_rp_bins=model[2],model_rt_bins=model[3])
    raw,den,matrix = oracle(pos,w,ids,lam,p,order,grid,model)
    np.testing.assert_allclose(result['weights'],den,rtol=2e-14,atol=2e-14)
    np.testing.assert_allclose(result['numerator'],raw,rtol=2e-11,atol=2e-11)
    np.testing.assert_allclose(result['matrix'],matrix,rtol=2e-11,atol=2e-12)


def test_full_support_null_mode_threads_and_affine_log_wavelength():
    pos,w,ids,lam,p = fixture()
    kw = dict(projection_weights=p,rp_max=40,rt_max=60,rp_bins=5,rt_bins=5)
    one = lya_distortion_matrix(pos,w,ids,lam,**kw)
    four = lya_distortion_matrix(pos,w,ids,lam,number_threads=4,**kw)
    np.testing.assert_array_equal(one['matrix'],four['matrix'])
    valid = one['valid_rows']
    np.testing.assert_allclose(one['matrix'][valid].sum(axis=1),0,atol=2e-14)
    assert np.all(np.isnan(one['matrix'][~valid]))
    assert np.any(one['matrix'][valid]<0)
    shift = lya_distortion_matrix(pos,w,ids,lam*np.log(10)+7,**kw)
    np.testing.assert_allclose(one['matrix'],shift['matrix'],atol=3e-12)
    scaled = lya_distortion_matrix(pos,w*3,ids,lam,**{**kw,'projection_weights':p*8})
    np.testing.assert_allclose(one['matrix'],scaled['matrix'],atol=2e-14)
    perm = np.random.default_rng(9).permutation(len(w))
    shuffled = lya_distortion_matrix(pos[perm],w[perm],ids[perm],lam[perm],
        **{**kw,'projection_weights':p[perm]})
    np.testing.assert_allclose(one['matrix'],shuffled['matrix'],atol=2e-14)


def test_degenerate_slope_and_identity():
    pos,w,ids,lam,p = fixture()
    lam[:] = 3.7
    kw = dict(rp_bins=4,rt_bins=3,rp_max=40,rt_max=60)
    mean = lya_distortion_matrix(pos,w,ids,projection_order=0,**kw)
    slope = lya_distortion_matrix(pos,w,ids,lam,**kw)
    np.testing.assert_array_equal(mean['matrix'],slope['matrix'])
    assert slope['stats']['degenerate_slope_forests'] == 3
    identity = lya_distortion_matrix(pos,w,ids,projection_order=-1,**kw)
    valid = identity['valid_rows']
    np.testing.assert_allclose(identity['matrix'][valid],np.eye(12)[valid],atol=1e-15)
    # A singleton projector vanishes exactly.
    single = lya_distortion_matrix(pos[:3],w[:3],np.arange(3),lam[:3],**kw)
    np.testing.assert_allclose(single['matrix'][single['valid_rows']],0,atol=2e-15)


def test_input_and_workspace_validation(monkeypatch):
    pos,w,ids,lam,p = fixture()
    with pytest.raises(Exception,match='true log_lambda'):
        lya_distortion_matrix(pos,w,ids)
    for kwargs,match in [({'max_memory_mb':.001},'workspace'),
                         ({'rp_bins':0},'invalid'),({'number_threads':0},'invalid'),
                         ({'projection_weights':p*0},'projection-weight')]:
        with pytest.raises(Exception,match=match):
            lya_distortion_matrix(pos,w,ids,lam,**kwargs)
    with pytest.raises(Exception,match='integ'):
        lya_distortion_matrix(pos,w,ids.astype(float),lam)
    with pytest.raises(Exception,match='int64'):
        lya_distortion_matrix(pos,w,np.full(len(ids),2**63,dtype='uint64'),lam)
    for bad in [np.full_like(w,-1),np.full_like(w,np.nan)]:
        with pytest.raises(Exception):
            lya_distortion_matrix(pos,bad,ids,lam)
    with pytest.raises(Exception,match='observer distance'):
        lya_distortion_matrix(pos*0,w,ids,lam)
    monkeypatch.setenv('CBALLS_MEMORY_BUDGET_MB','.01')
    with pytest.raises(Exception,match='budget|BUDGET'):
        lya_distortion_matrix(pos,w,ids,lam)


def test_fitting_domain_cache_selection_and_reprojection(tmp_path):
    import lya_corr_all_engines as driver
    import lya_fits
    from test_lya_analysis import eboss_file
    from lya_distortion import compute_catalog_distortion, DOMAIN_FIELDS
    path = eboss_file(tmp_path/'forests.fits')
    cat = lya_fits.read_eboss([str(path)],pixel_stride=2,z_max=1.97,
                            retain_projection_domain=True,redshift_weight_exponent=1.9)
    domain = cat.distortion_domain
    assert len(domain['positions']) == 9 and cat.nbody == 3
    np.testing.assert_array_equal(domain['projection_weights'],np.tile([1,2,3],3))
    cfg = driver.RunConfig(('lya-2pcf-omp',),tmp_path/'run',compute_distortion=True,
        rp_bins=3,rt_bins=3,rp_max=40,rt_max=150)
    computed = compute_catalog_distortion(cat,cfg)
    weights = np.zeros(9)
    weights[domain['data_indices']] = cat.weights
    raw,den,dm = oracle(domain['positions'],weights,domain['forest_ids'],
        domain['log_lambda'],domain['projection_weights'],1,grid=(40,150,3,3))
    np.testing.assert_allclose(computed['matrix'],dm,atol=2e-12)
    driver.save_catalog(tmp_path/'cache.npz',cat)
    cached = driver.read_npz(tmp_path/'cache.npz')
    np.testing.assert_array_equal(cached.log_lambda,cat.log_lambda)
    for key in DOMAIN_FIELDS:
        np.testing.assert_array_equal(cached.distortion_domain[key],domain[key])
    np.testing.assert_array_equal(compute_catalog_distortion(cached,cfg)['matrix'],computed['matrix'])
    projected = lya_fits.read_eboss([str(path)],pixel_stride=2,project_delta=True,
        retain_projection_domain=True,redshift_weight_exponent=1.9)
    assert projected.distortion_domain is None
    result = compute_catalog_distortion(projected,cfg)
    direct = lya_distortion_matrix(projected.positions,projected.weights,projected.forest_ids,
        projected.log_lambda,rp_max=40,rt_max=150,rp_bins=3,rt_bins=3)
    np.testing.assert_array_equal(result['matrix'],direct['matrix'])
    cfg.distortion_projection_order = 0
    with pytest.raises(ValueError,match='mean.slope'):
        compute_catalog_distortion(projected,cfg)
    incomplete = lya_fits.read_eboss([str(path)],pixel_stride=2)
    with pytest.raises(ValueError,match='full continuum-fitting domain'):
        compute_catalog_distortion(incomplete,cfg)


def test_desi_wavelength_domain_and_export(tmp_path):
    import lya_corr_all_engines as driver
    import lya_analysis
    from lya_distortion import compute_catalog_distortion,save_distortion
    from astropy.io import fits
    wave = np.array([3700.,3800.,3900.,4000.])
    meta = fits.BinTableHDU.from_columns([
        fits.Column(name='LOS_ID',format='K',array=np.array([2**55+1,2**55+2,2**55+3])),
        fits.Column(name='RA',format='D',unit='rad',array=[0,.01,.02]),
        fits.Column(name='DEC',format='D',unit='rad',array=[0,0,0])],name='METADATA')
    path = tmp_path/'desi.fits'
    fits.HDUList([fits.PrimaryHDU(),fits.ImageHDU(wave,name='LAMBDA'),meta,
        fits.ImageHDU(np.arange(12.).reshape(3,4)/10,name='DELTA'),
        fits.ImageHDU(np.ones((3,4)),name='WEIGHT')]).writeto(path)
    cat = driver.read_desi([str(path)],pixel_stride=2,retain_projection_domain=True)
    assert cat.nbody == 6 and len(cat.distortion_domain['positions']) == 12
    np.testing.assert_array_equal(cat.log_lambda,np.tile(np.log10(wave[::2]),3))
    cfg = driver.RunConfig(('lya-2pcf-omp',),tmp_path/'run',compute_distortion=True,
        rp_bins=3,rt_bins=3,rp_max=300,rt_max=150)
    response = compute_catalog_distortion(cat,cfg)
    report = save_distortion(tmp_path/'dm',response)
    for kind in ('fits','npz'):
        loaded = lya_analysis.load_analysis_array(report[kind],'distortion',cfg)
        np.testing.assert_array_equal(loaded,response['matrix'])
        model = np.arange(9.)/100
        np.testing.assert_allclose(lya_analysis.apply_distortion(loaded,model),response['matrix']@model)
    with fits.open(report['fits']) as hdus:
        np.testing.assert_array_equal(hdus['DISTORTION'].data['WDM'],response['weights'])
        assert hdus['DISTORTION'].header['PROJORD'] == 1
    with pytest.raises(FileExistsError):
        save_distortion(tmp_path/'dm',response)
    broken = dict(cat.distortion_domain)
    broken['data_indices'] = np.zeros(cat.nbody,dtype=int)
    cat.distortion_domain = broken
    with pytest.raises(ValueError,match='data_indices'):
        cat.normalized()


def test_configuration_and_projected_covariance_response(tmp_path):
    import lya_corr_all_engines as driver
    cfg = driver.RunConfig(('lya-2pcf-omp',),tmp_path,compute_distortion=True,
        rp_bins=4,rt_bins=3,rp_max=20,rt_max=40)
    assert cfg.validate() is cfg
    bad = driver.RunConfig(('lya-3pcf-omp',),tmp_path,compute_distortion=True)
    with pytest.raises(ValueError,match='3D 2PCF'):
        bad.validate()
    # Propagate an arbitrary cross-forest covariance directly through P C P^T,
    # then bin it. This checks matrix orientation independently of the expansion.
    pos,w,ids,lam,p = fixture()
    model = np.random.default_rng(16).normal(size=12)
    result = lya_distortion_matrix(pos,w,ids,lam,projection_weights=p,
        rp_max=20,rt_max=40,rp_bins=4,rt_bins=3)
    pair_bins = bins(pos,20,40,4,3)
    predicted = np.zeros(12)
    labels = np.unique(ids)
    for f,label in enumerate(labels):
        a = np.flatnonzero(ids == label)
        for other in labels[f+1:]:
            b = np.flatnonzero(ids == other)
            pp = []
            for k in (a,b):
                x = lam[k]-np.average(lam[k],weights=p[k])
                pp.append(np.eye(len(k))-np.ones((len(k),1))*p[k]/p[k].sum()
                          -x[:,None]*(p[k]*x)/np.dot(p[k],x*x))
            index = pair_bins[np.ix_(a,b)]
            covariance = np.where(index >= 0,model[np.maximum(index,0)],0)
            covariance = pp[0]@covariance@pp[1].T
            for ii,i in enumerate(a):
                for jj,j in enumerate(b):
                    if index[ii,jj] >= 0:
                        predicted[index[ii,jj]] += w[i]*w[j]*covariance[ii,jj]
    valid = result['valid_rows']
    np.testing.assert_allclose((result['matrix']@model)[valid],
        predicted[valid]/result['weights'][valid],atol=3e-12)


def test_driver_computes_once_and_matches_native_pair_weights(tmp_path,capsys):
    import json
    import lya_corr_all_engines as driver
    from test_lya_analysis import eboss_file
    path = eboss_file(tmp_path/'forests.fits')
    output = tmp_path/'run'
    model = tmp_path/'model.npy'
    np.save(model,np.linspace(-.01,.02,9))
    command = ['--fits',str(path),'--pixel-stride','2','--engine','lya-2pcf-omp','lya-los-tree-2pcf-omp',
        '--threads','2','--compute-distortion','--model-correlation',str(model),
        '--rp-max','40','--rt-max','150','--rp-bins','3','--rt-bins','3',
        '--wedge-bins','3','--wedge-subsamples','2','--no-plots','--verbose','0',
        '--output',str(output)]
    assert driver.main(command) == 0
    assert capsys.readouterr().out.count('Computing continuum distortion') == 1
    summary = json.loads((output/'summary.json').read_text())
    assert summary['computed_distortion']['shape'] == [9,9]
    assert summary['validation']['passed']
    with np.load(output/'distortion/distortion.npz',allow_pickle=False) as data:
        weights = data['weights']
    for engine in ('lya-2pcf-omp','lya-los-tree-2pcf-omp'):
        table = np.loadtxt(output/engine/'histXi2pcf_lya.txt',ndmin=2)
        # Histograms publish zero-based (rp,rt) bin indices.
        indices = table[:,0].astype(int)*3+table[:,1].astype(int)
        populated = np.zeros(9)
        populated[indices] = table[:,-1]
        np.testing.assert_allclose(weights,populated,atol=2e-13)


def test_cartesian_mean_only_full_domain_and_cache(tmp_path):
    from astropy.io import fits
    import lya_corr_all_engines as driver
    from lya_fits import read_fits
    from lya_distortion import compute_catalog_distortion
    pos,w,ids,lam,p = fixture()
    path = tmp_path/'cartesian.fits'
    hdu = fits.BinTableHDU.from_columns([
        *[fits.Column(name=name,format='D',array=pos[:,i]) for i,name in enumerate(('x','y','z'))],
        fits.Column(name='delta',format='D',array=w*0),
        fits.Column(name='weight',format='D',array=w),
        fits.Column(name='los_id',format='K',array=ids)])
    fits.HDUList([fits.PrimaryHDU(),hdu]).writeto(path)
    cat = read_fits([str(path)],pixel_stride=2,retain_projection_domain=True)
    assert cat.log_lambda is None and 'log_lambda' not in cat.distortion_domain
    cfg = driver.RunConfig(('lya-2pcf-omp',),tmp_path/'run',compute_distortion=True,
        distortion_projection_order=0,rp_max=20,rt_max=40,rp_bins=4,rt_bins=3)
    result = compute_catalog_distortion(cat,cfg)
    selected = np.zeros(len(w))
    selected[cat.distortion_domain['data_indices']] = cat.weights
    raw,den,matrix = oracle(pos,selected,ids,lam,w,0)
    np.testing.assert_allclose(result['matrix'],matrix,atol=1e-14)
    driver.save_catalog(tmp_path/'cache.npz',cat)
    cached = driver.read_npz(tmp_path/'cache.npz')
    np.testing.assert_array_equal(compute_catalog_distortion(cached,cfg)['matrix'],result['matrix'])
    cfg.distortion_projection_order = 1
    with pytest.raises(Exception,match='true log_lambda'):
        compute_catalog_distortion(cached,cfg)


def test_multiple_worker_task_reuse():
    rng = np.random.default_rng(225)
    ids = np.repeat(np.arange(19),7)
    angles = rng.uniform(-.1,.1,(19,2))[ids]
    los = np.column_stack([np.cos(angles[:,0])*np.cos(angles[:,1]),
        np.sin(angles[:,0])*np.cos(angles[:,1]),np.sin(angles[:,1])])
    radius = rng.uniform(80,120,len(ids))
    pos,lam = radius[:,None]*los,3.5+radius*.0001
    w,p = rng.uniform(.5,3.,(2,len(ids)))
    kw = dict(projection_weights=p,rp_max=40,rt_max=30,rp_bins=5,rt_bins=6)
    serial = lya_distortion_matrix(pos,w,ids,lam,**kw)
    parallel = lya_distortion_matrix(pos,w,ids,lam,number_threads=4,**kw)
    for name in ('matrix','numerator','weights'):
        np.testing.assert_array_equal(serial[name],parallel[name])
