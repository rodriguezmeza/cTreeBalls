"""Production frontend contracts for signed anisotropic moments and exact bins."""
import json
from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lya_corr_all_engines as driver
import lya_anisotropic_multipoles as products

METHOD = driver.ANISOTROPIC_MULTIPOLE


def config(path, **kwargs):
    return driver.RunConfig((METHOD,), path, r3_max=10., r3_bins=1,
                            theta_bins=1, mu_bins=2, lya3_lmax=1,
                            plots=False, **kwargs)


def files(path):
    path.mkdir(exist_ok=True)
    # A positive angular measure at mu=-1 has D0=4, D1=-4. Finite L=1
    # reconstruction has a negative denominator in the other angular bin.
    raw = np.array([[0, 0, 0, 0, 0, 2., 4.], [0, 0, 0, 0, 1, -2., -4.]])
    finite = np.array([
        [0, 0, 0, 0, 0, 5, 5, np.pi/2, np.pi/2, -.5, .5, 2.5, 5.],
        [0, 0, 0, 0, 1, 5, 5, np.pi/2, np.pi/2,  .5, np.nan, -.5, -1.]])
    exact = np.array([[0, 0, 0, 0, 0, 5, 5, np.pi/2, np.pi/2, -.5, .5, 2., 4.]])
    np.savetxt(path/'histZetaM_lya_multipoles.txt', raw, fmt='%.17g')
    np.savetxt(path/'histZetaM_lya5d_multipole.txt', finite, fmt='%.17g')
    np.savetxt(path/'histZetaM_lya5d.txt', exact, fmt='%.17g')
    return raw, finite, exact


def test_explicit_registration_preserves_bulk_engine_selection():
    available = list(driver.LYA_ENGINES)
    assert driver.resolve_engines([METHOD], available, '3pcf') == (METHOD,)
    assert METHOD not in driver.resolve_engines(['all'], available, 'both')
    assert driver.LYA_ENGINES[METHOD].explicit_only
    assert METHOD not in driver.INCOMPATIBLE_ENGINE_REASONS


def test_kernel_dispatch_and_parameter_reset(tmp_path):
    cfg = driver.RunConfig((METHOD, 'lya-los-tree-3pcf-omp'), tmp_path,
                            lya3_kernel=5, lya3_lmax=12, lya3_multipole_kernel=2)
    cfg.validate()
    p = driver.engine_parameters(cfg, METHOD, tmp_path)
    assert p['lya3Kernel'] == 2 and p['lya3LMax'] == 12 and p['lya3MuMode'] == 1
    assert p['rangeN'] == cfg.r3_max
    for method in ('lya-los-tree-3pcf-omp', 'octree-3pcf-3d-omp', 'lya-2pcf-omp'):
        p = driver.engine_parameters(cfg, method, tmp_path)
        assert p['lya3MuMode'] == 0


@pytest.mark.parametrize('kwargs,reason', [
    ({'lya3_kernel':5},'hard-bin'), ({'lya3_lmax':33},'lya3_lmax'),
    ({'lya3_lmax':-1},'lya3_lmax'), ({'lya3_mu_mode':'invalid'},'lya3_mu_mode'),
    ({'lya3_multipole_kernel':3},'lya3_multipole_kernel'),
    ({'qualification_reference':Path('reference')},'moment packets')])
def test_unsupported_controls_fail_before_loading_catalog(tmp_path, kwargs, reason):
    with pytest.raises(ValueError, match=reason):
        driver.RunConfig((METHOD,), tmp_path, **kwargs).validate()


def test_memory_estimate_counts_both_moment_and_exact_grids(tmp_path):
    args = dict(engines=(METHOD,),output_dir=tmp_path,r3_bins=10,theta_bins=10,
                mu_bins=64,lya3_lmax=0,threads=4,max_hist_mib=10)
    driver.RunConfig(**args,lya3_mu_mode='finite-multipole').validate()
    with pytest.raises(ValueError, match='histogram estimate'):
        driver.RunConfig(**args,lya3_mu_mode='exact').validate()


def test_signed_moments_and_invalid_reconstruction_remain_distinct(tmp_path):
    raw, finite, exact = files(tmp_path)
    got_raw, got_finite = products.read_outputs(tmp_path, config(tmp_path))
    np.testing.assert_array_equal(got_raw, raw)
    np.testing.assert_allclose(got_finite, finite, equal_nan=True)
    report = products.inspect_outputs(tmp_path, config(tmp_path), exact)
    assert report['finite_reconstruction']['negative_denominator_bins'] == 1
    assert report['finite_reconstruction']['undefined_ratio_bins'] == 1
    diagnostic = report['finite_reconstruction_vs_exact']
    assert diagnostic['empty_bin_denominator_leakage_l1'] == 1
    assert diagnostic['numerator_relative_l2'] == pytest.approx(np.sqrt(.5)/2)
    assert not report['finite_reconstruction']['included_in_exact_comparisons']
    assert report['exact_bins']['available']
    np.testing.assert_array_equal(driver.read_products(tmp_path,METHOD,config(tmp_path))['3d_3pcf'], exact)
    finite_cfg = config(tmp_path,lya3_mu_mode='finite-multipole')
    assert driver.read_products(tmp_path,METHOD,finite_cfg) == {}  # Ignore even a stale exact filename.
    report = products.inspect_outputs(tmp_path,finite_cfg)
    assert not report['exact_bins']['available']
    assert 'finite_reconstruction_vs_exact' not in report
    json.dumps(report,allow_nan=False)


@pytest.mark.parametrize('change,reason', [
    ('raw_nan','non-finite raw'), ('duplicate','duplicate'), ('ratio','ratio'),
    ('nan_convention','NaN convention'), ('centers','bin coordinates')])
def test_malformed_auxiliary_outputs_are_rejected(tmp_path, change, reason):
    raw, finite, _ = files(tmp_path)
    if change == 'raw_nan':raw[0,-2] = np.nan
    if change == 'duplicate':raw[1,:5] = raw[0,:5]
    if change == 'ratio':finite[0,-3] = 123
    if change == 'nan_convention':finite[1,-3] = 0
    if change == 'centers':finite[:,5] *= 2
    np.savetxt(tmp_path/'histZetaM_lya_multipoles.txt',raw)
    np.savetxt(tmp_path/'histZetaM_lya5d_multipole.txt',finite)
    with pytest.raises(ValueError,match=reason):
        products.read_outputs(tmp_path,config(tmp_path))


@pytest.mark.parametrize('mode',['exact','finite-multipole'])
def test_public_driver_matches_independent_ordered_triangle_moments(tmp_path,mode):
    cy = pytest.importorskip('cyballs')
    if cy.search_method_id(METHOD)<0:pytest.skip('anisotropic addon not built')
    from test_lya_triplet_acceleration import moment_oracle
    cat = driver.synthetic_catalog(4,4,seed=253)
    engines = (METHOD,'octree-3pcf-3d-omp','lya-los-tree-3pcf-omp')
    cfg = driver.RunConfig(engines,tmp_path/mode,threads=2,plots=False,analysis=False,
            r3_max=80.,r3_bins=2,theta_bins=3,mu_bins=4,lya3_lmax=4,lya3_mu_mode=mode)
    result = driver.run_engine_suite(cat,cfg)
    points = np.column_stack([cat.positions,cat.delta,cat.weights,cat.forest_ids])
    oracle = moment_oracle(points,4,80.,2,3)
    raw = np.loadtxt(cfg.output_dir/METHOD/'histZetaM_lya_multipoles.txt')
    np.testing.assert_allclose(raw[:,-2:].reshape(oracle.shape),oracle,rtol=2e-11,atol=3e-11)
    summary = json.loads((cfg.output_dir/'summary.json').read_text())
    assert summary['validation']['passed']
    assert summary['engines'][METHOD]['anisotropic_multipole']['mode'] == mode
    if mode == 'exact':
        assert result[METHOD]['products'].keys() == {'3d_3pcf'}
        assert summary['validation']['comparisons'] == 1
    else:
        assert result[METHOD]['products'] == {}
        assert summary['validation']['comparisons'] == 0
        assert 'not an accuracy certificate' in summary['validation']['scope']


def test_public_cli_defaults_and_auxiliary_plots(tmp_path):
    args = driver.parse_arguments(['--engine',METHOD,'--statistics','3pcf'])
    assert (args.lya3_lmax,args.lya3_mu_mode,args.lya3_multipole_kernel) == (8,'exact',0)
    _,_,exact = files(tmp_path)
    cfg = config(tmp_path)
    cfg.plots = True
    report = products.inspect_outputs(tmp_path,cfg,exact)
    assert len(report['plots']) == 8
    assert all(Path(p).stat().st_size>1000 for p in report['plots'])
