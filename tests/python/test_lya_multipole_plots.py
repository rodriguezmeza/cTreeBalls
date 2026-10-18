"""Scientific normalization and saved-output plotting contracts."""
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lya_corr_all_engines as driver
import lya_multipole_plots as plots


def config(path, **kwargs):
    return driver.RunConfig((driver.ANISOTROPIC_MULTIPOLE,), path, r3_max=10.,
        r3_bins=1, theta_bins=2, mu_bins=4, lya3_lmax=4, **kwargs)


def test_polar_sums_precede_monopole_normalization(tmp_path):
    cfg = config(tmp_path)
    index = np.indices((1, 1, 2, 2, 5)).reshape(5, -1).T
    raw = np.column_stack((index, np.zeros((len(index), 2))))
    field = raw[:, -2:].reshape(1, 1, 2, 2, 5, 2)
    field[0, 0, 0, 0, :, 0] = [1, 2, -3, 0, 0]
    field[0, 0, 1, 1, :, 0] = [9, -1, 5, 0, 0]
    field[0, 0, 0, 0, :, 1] = [1, -2, 0, 0, 0]
    field[0, 0, 1, 1, :, 1] = [3, 2, -4, 0, 0]
    num, den, q, window = plots.anisotropic_radial_moments(raw, cfg)
    np.testing.assert_array_equal(num[0, 0, :3], [10, 1, 2])
    np.testing.assert_array_equal(den[0, 0, :3], [4, 0, -4])
    np.testing.assert_allclose(q[0, 0, :3], [2.5, .75, 2.5])
    np.testing.assert_allclose(window[0, 0, :3], [1, 0, -5])
    raw[:, -2:] = 0
    assert np.isnan(plots.anisotropic_radial_moments(raw, cfg)[2]).all()


def linear_binned_table(cfg):
    rows = []
    # Two polar cells with unequal weights: pool N/D, never average ratios.
    for mu_bin in range(cfg.mu_bins):
        mu = -1+(mu_bin+.5)*2/cfg.mu_bins
        value = .2-.04*mu
        for polar, weight, multiplier in [(0, 1., 2.), (1, 3., 2/3)]:
            zeta = value*multiplier
            rows.append([0, 0, polar, polar, mu_bin, 5, 5, .5, .5, mu,
                         zeta, zeta*weight, weight])
    return np.asarray(rows)


def test_binned_projection_integrates_legendre_windows_and_masks_missing_mu(tmp_path):
    cfg = config(tmp_path)
    exact = linear_binned_table(cfg)
    coefficients, coverage = plots.binned_correlation_multipoles(exact, cfg)
    # A linear function sampled at bin midpoints then held constant has this
    # exact l=1 attenuation. This intentionally exposes finite-bin resolution.
    np.testing.assert_allclose(coefficients[0, 0, :3], [.2, -.04*(1-1/16), 0], atol=1e-16)
    assert coefficients.shape == (1, 1, 4)  # no l>=N_mu extrapolation
    assert coverage[0, 0] == 1
    missing = exact[exact[:, 4] != 2]
    coefficients, coverage = plots.binned_correlation_multipoles(missing, cfg)
    assert np.isnan(coefficients).all()
    assert coverage[0, 0] == .75


@pytest.mark.parametrize('orders', [(0, 0), (-1,), (33,), ()])
def test_bad_order_selection_is_rejected(tmp_path, orders):
    with pytest.raises(ValueError, match='distinct integers'):
        plots.selected_orders(config(tmp_path, plot_multipole_orders=orders), range(5))


def test_native_scalar_normalization_and_empty_shells_are_preserved(tmp_path, monkeypatch):
    captured = {}
    def capture(root, stem, values, centers, orders, *args, **kwargs):
        captured.update(values=values, centers=centers, orders=orders, bounds=kwargs['bounds'])
        return ['a.png']
    monkeypatch.setattr(plots, 'radial_views', capture)
    cfg = driver.RunConfig(('octree-3pcf-3d-omp',), tmp_path,
        r3_max=10., r3_bins=2, multipole_lmax=2, multipole_rmin=2.)
    table = np.array([[0, 1, 2, 4, 8, .1, 1, 10],
                      [1, 1, 2, 4, 8, -.2, -2, 10],
                      [2, 1, 2, 4, 8, 0., 0, 0]])
    assert plots.write_scalar_plots(tmp_path, table, cfg) == ['a.png']
    np.testing.assert_array_equal(captured['values'][0, 1, :2], [.1, -.2])
    assert np.isnan(captured['values'][0, 1, 2])
    assert np.isnan(captured['values'][1, 0]).all()
    np.testing.assert_array_equal(captured['centers'], [4, 8])
    assert captured['bounds'] == (2., 10.)


def test_plot_existing_never_loads_catalog_or_calls_discovery(tmp_path, monkeypatch):
    from test_lya_anisotropic_driver import files, config as small_config
    source, destination = tmp_path/'saved', tmp_path/'plots'
    engine = source/driver.ANISOTROPIC_MULTIPOLE
    engine.mkdir(parents=True)
    files(engine)
    cfg = small_config(source)
    settings = dict(engines=list(cfg.engines), output_dir=str(source),
        r3_max=10., r3_bins=1, theta_bins=1, mu_bins=2, lya3_lmax=1, lya3_mu_mode='exact')
    summary = source/'summary.json'
    summary.write_text(json.dumps(dict(config=settings, engines={cfg.engines[0]:{}})))
    original = {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in source.rglob('*') if p.is_file()}
    def forbidden(*args, **kwargs):
        pytest.fail('plot-existing invoked catalog/native search setup')
    monkeypatch.setattr(driver, 'load_catalog', forbidden)
    monkeypatch.setattr(driver, 'discover_cython_methods', forbidden)
    monkeypatch.setattr(driver, 'get_mpi_comm', forbidden)
    assert driver.main(['--plot-existing', str(source), '--output', str(destination),
                        '--plot-multipole-orders', '0', '1']) == 0
    report = json.loads((destination/'plot_summary.json').read_text())
    assert report['native_search_run'] is False
    assert len(report['plots']) == 8
    assert all(Path(p).stat().st_size > 1000 for p in report['plots'])
    assert all(hashlib.sha256(p.read_bytes()).hexdigest() == digest for p, digest in original.items())
    csv = np.genfromtxt(destination/cfg.engines[0]/'anisotropic_multipoles.csv', delimiter=',', names=True)
    np.testing.assert_allclose(csv['Q_ell'], [.5, -1.5])
    details = report['engines'][cfg.engines[0]]['multipole_plotting']
    assert details['window_corrected'] is False
    assert details['binned_correlation']['complete_radial_cells'] == 0


def test_order_pagination_and_unsupported_only_selection(tmp_path):
    cfg = config(tmp_path, plot_multipole_orders=(8,))
    raw = np.column_stack((np.indices((1,1,2,2,5)).reshape(5,-1).T, np.zeros((20,2))))
    with pytest.raises(ValueError, match='none of'):
        plots.write_anisotropic_plots(tmp_path, raw, cfg)
    paths = plots.radial_views(tmp_path, 'empty', np.full((1,1,10), np.nan),
        np.array([5.]), list(range(10)), 'Pagination with undefined bins', 'test')
    assert len(paths) == 4
    assert all(Path(p).exists() for p in paths)


def test_regular_suite_writes_multipole_views_and_csv(tmp_path):
    cy = pytest.importorskip('cyballs')
    method = driver.ANISOTROPIC_MULTIPOLE
    if cy.search_method_id(method) < 0:
        pytest.skip('anisotropic addon not built')
    cfg = driver.RunConfig((method,), tmp_path, threads=2, analysis=False,
        r3_max=80., r3_bins=2, theta_bins=2, mu_bins=4, lya3_lmax=2,
        plot_multipole_orders=(0, 2))
    driver.run_engine_suite(driver.synthetic_catalog(4, 4, seed=253), cfg)
    saved = json.loads((tmp_path/'summary.json').read_text())
    detail = saved['engines'][method]['anisotropic_multipole']['multipole_plotting']
    assert detail['plotted_orders'] == [0, 2]
    assert Path(detail['csv']).exists()
    assert any('anisotropic_multipoles_maps' in p for p in saved['plots'])
