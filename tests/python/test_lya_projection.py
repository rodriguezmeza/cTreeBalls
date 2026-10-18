"""Match fill_dw and apply the same corrected field to 2PCF/3PCF engines."""
import ast
import json
from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lya_projection import project_lya_delta, project_forest_delta, forest_projection_diagnostics
import lya_corr_all_engines as driver


REFERENCE = Path('/Users/mar/Documents/Codex/lyman_alpha/lya2pcf_2026-03-14/forest_class.py')


def upstream_fill_dw():
    if not REFERENCE.exists():
        pytest.skip('user lya2pcf checkout not present')
    tree = ast.parse(REFERENCE.read_text())
    cls = next(node for node in tree.body if isinstance(node, ast.ClassDef) and node.name == 'quasar')
    method = next(node for node in cls.body if isinstance(node, ast.FunctionDef) and node.name == 'fill_dw')
    namespace = dict(np=np)
    exec(compile(ast.Module(body=[method], type_ignores=[]), str(REFERENCE), 'exec'), namespace)
    return namespace['fill_dw']


def sample():
    rng = np.random.default_rng(1782)
    ids = np.repeat(np.array([2**60+17, -2**59+1, 5], dtype=np.int64), [5, 7, 4])
    lam = np.concatenate([3.5+np.arange(n)*.003 for n in [5, 7, 4]])
    w = rng.uniform(.15, 5, len(ids))
    delta = .3+4*(lam-3.5)+rng.normal(0, .2, len(ids))
    return delta, w, ids, lam


def matrix_reference(delta, w, ids, lam):
    out = np.zeros_like(delta)
    for label in np.unique(ids):
        index = np.flatnonzero(ids == label)
        p = w[index]
        x = lam[index]-np.average(lam[index], weights=p)
        projector = np.eye(len(index))-np.ones((len(index), 1))*p[None, :]/p.sum()
        if np.ptp(lam[index][p > 0]) > 0:
            projector -= x[:, None]*(p*x)[None, :]/np.dot(p, x*x)
        out[index] = projector@delta[index]
    return out


def sightline_positions(ids, lam):
    labels = np.unique(ids, return_inverse=True)[1]
    angle = labels*.05
    radius = 100+(lam-3.5)*3000
    return radius[:, None]*np.column_stack((np.cos(angle), np.sin(angle), angle*0))


def test_matches_actual_reference_fill_dw_and_dense_projector():
    delta, w, ids, lam = sample()
    saved = [a.copy() for a in (delta, w, ids, lam)]
    got, report = project_lya_delta(delta, w, ids, lam, return_diagnostics=True)
    fill = upstream_fill_dw()
    for label in np.unique(ids):
        indices = ids == label
        qso = SimpleNamespace(we=w[indices])
        fill(qso, delta[indices], lam[indices], True)
        np.testing.assert_allclose(got[indices]*w[indices], qso.dw, rtol=2e-11, atol=2e-13)
    np.testing.assert_allclose(got, matrix_reference(delta, w, ids, lam), atol=3e-13)
    assert report['forests'] == 3
    assert report['max_abs_weighted_mean_after'] < 2e-16
    assert report['max_abs_normalized_first_moment_after'] < 2e-16
    for original, retained in zip(saved, (delta, w, ids, lam)):
        np.testing.assert_array_equal(original, retained)
    json.dumps(report, allow_nan=False)


def test_grouping_permutation_affine_log_coordinate_weight_scale_and_idempotence():
    delta, w, ids, lam = sample()
    got = project_lya_delta(delta, w, ids, lam)
    perm = np.random.default_rng(44).permutation(len(ids))
    shuffled = project_lya_delta(delta[perm], w[perm], ids[perm], lam[perm])
    np.testing.assert_allclose(shuffled, got[perm], atol=1e-15)
    np.testing.assert_allclose(project_lya_delta(delta, w*1e307, ids, lam), got, atol=1e-15)
    np.testing.assert_allclose(project_lya_delta(delta, w, ids, 7+np.log(10)*lam), got, atol=2e-13)
    np.testing.assert_allclose(project_lya_delta(got, w, ids, lam), got, atol=2e-16)


def test_constant_or_single_pixel_forests_and_zero_weights_are_finite():
    delta = np.array([1., 4., 9., 3., 7., 8.])
    weights = np.array([1., 3., 0., 1., 0., 0.])
    ids = np.array([4, 4, 4, 8, 9, 9])
    lam = np.array([3.6, 3.6, 1e308, 3.8, 4., 5.])
    got, report = project_lya_delta(delta, weights, ids, lam, return_diagnostics=True)
    np.testing.assert_array_equal(got, [-2.25, .75, 9., 0., 7., 8.])
    assert report['mean_only_forests'] == 2 and report['zero_weight_forests'] == 1


@pytest.mark.parametrize('field,value,message', [
    ('delta', [1., np.nan], 'finite'), ('weights', [1., -1.], 'non-negative'),
    ('forest_ids', [1., 1.], 'integer'), ('log_lambda', [3.], 'same length'),
    ('forest_ids', np.array([2**63, 2**63], dtype='u8'), 'signed int64'),
    ('delta', [1+1j, 2+0j], 'real array')])
def test_invalid_projection_inputs_are_rejected(field, value, message):
    values = dict(delta=[1., 2.], weights=[1., 2.], forest_ids=np.array([1, 1]), log_lambda=[3., 4.])
    values[field] = value
    with pytest.raises(ValueError, match=message):
        project_lya_delta(**values)


def test_fits_preprocessing_uses_final_weights_before_projection():
    from lya_fits import preprocess_forest
    delta, w, ids, lam = sample()
    chosen = ids == ids[0]
    delta, w, lam = delta[chosen], w[chosen], lam[chosen]
    original_delta, original_w = delta.copy(), w.copy()
    got, adjusted_w = preprocess_forest(delta, w, lam, project_delta=True,
        redshift_weight_exponent=1.9, weight_z_ref=2.25)
    expected_w = w*(10**lam/1215.67/3.25)**1.9
    np.testing.assert_allclose(adjusted_w, expected_w)
    qso = SimpleNamespace(we=expected_w)
    upstream_fill_dw()(qso, delta, lam, True)
    np.testing.assert_allclose(got*adjusted_w, qso.dw, atol=2e-13)
    np.testing.assert_array_equal(delta, original_delta)
    np.testing.assert_array_equal(w, original_w)


def test_npz_projection_matches_manual_field_and_records_domain(tmp_path):
    delta, w, ids, lam = sample()
    cat = driver.synthetic_catalog(4, 4)
    cat.delta, cat.weights, cat.forest_ids, cat.log_lambda = delta, w, ids, lam
    cat.positions = sightline_positions(ids, lam)
    cache = tmp_path/'unprojected.npz'
    driver.save_catalog(cache, cat)
    got = driver.load_catalog(driver.parse_arguments(['--catalog', str(cache), '--project-delta']))
    np.testing.assert_allclose(got.delta, project_lya_delta(delta, w, ids, lam), atol=1e-15)
    np.testing.assert_array_equal(got.weights, w)
    assert got.metadata['project_delta'] and got.distortion_domain is None
    assert got.metadata['projection_constraints']['forests'] == 3
    assert got.metadata['projection_constraints']['max_abs_weighted_mean_after'] < 1e-15
    unprojected = driver.load_catalog(driver.parse_arguments(['--catalog', str(cache)]))
    np.testing.assert_array_equal(unprojected.delta, delta)
    cat.log_lambda = None
    driver.save_catalog(cache, cat)
    with pytest.raises(ValueError, match='true log_lambda'):
        driver.load_catalog(driver.parse_arguments(['--catalog', str(cache), '--project-delta']))


def test_compiled_projector_matches_shared_source():
    import cyballs
    delta, w, ids, lam = sample()
    np.testing.assert_allclose(cyballs.project_lya_delta(delta, w, ids, lam),
                               project_lya_delta(delta, w, ids, lam), atol=1e-15)
    np.testing.assert_allclose(cyballs.project_forest_delta(delta, w, lam),
                               project_forest_delta(delta, w, lam), atol=1e-15)


@pytest.mark.parametrize('method', ['lya-2pcf-3pcf-omp', 'lya-los-tree-2pcf-3pcf-omp',
    'lya-rprt-quadtree-2pcf-omp', 'lya-1d-2pcf-3pcf-omp', 'lya-1d-tree-3pcf-omp',
    'lya-1d-tree-same-los-2pcf-omp', 'lya-anisotropic-multipole-3pcf-omp'])
def test_native_registration_projection_matches_manual_for_engines(tmp_path, method):
    from cyballs import cballs, search_method_id
    if search_method_id(method) < 0:
        pytest.skip('engine not built')
    delta, w, ids, lam = sample()
    cat = driver.synthetic_catalog(4, 4, seed=253)
    cat.positions = sightline_positions(ids, lam)
    cfg = driver.RunConfig((method,), tmp_path, r3_max=80, r3_bins=2,
        theta_bins=2, mu_bins=4, lya3_lmax=2, plots=False, analysis=False)
    projected = project_lya_delta(delta, w, ids, lam)
    outputs = []
    for label, field, automatic in [('automatic', delta, True), ('manual', projected, False)]:
        balls = cballs()
        try:
            balls.set_forest_catalog(cat.positions, field, w, ids,
                project_delta=automatic, log_lambda=lam if automatic else None)
            root = tmp_path/label
            root.mkdir()
            balls.set(driver.engine_parameters(cfg, method, root))
            balls.Run()
            if method == driver.ANISOTROPIC_MULTIPOLE:
                outputs.append(np.loadtxt(root/'histZetaM_lya_multipoles.txt'))
            else:
                outputs.append(driver.read_products(root, method, cfg))
        finally:
            balls.struct_cleanup()
    if method == driver.ANISOTROPIC_MULTIPOLE:
        np.testing.assert_allclose(outputs[0], outputs[1], rtol=2e-12, atol=2e-13)
    else:
        for key in outputs[0]:
            np.testing.assert_allclose(outputs[0][key], outputs[1][key], rtol=2e-12, atol=2e-13)


def test_native_registration_requires_wavelength_for_requested_projection():
    from cyballs import cballs
    cat = driver.synthetic_catalog(3, 3)
    balls = cballs()
    try:
        with pytest.raises(Exception, match='true log_lambda'):
            balls.set_forest_catalog(cat.positions, cat.delta, cat.weights, cat.forest_ids, project_delta=True)
        assert balls.catalog_count == 0
    finally:
        balls.struct_cleanup()
