"""Input, estimator, reference-kernel, and forward-model regression tests."""
from __future__ import annotations

from itertools import combinations
import importlib.util
import json
import os
from pathlib import Path
import sys

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/"tests/python"))
import lya_corr_all_engines as driver
import lya_analysis as analysis
import lya_fits
from lya_reference import run_references, product_table


def config(tmp_path, **kwargs):
    return driver.RunConfig(("lya-2pcf-omp",), tmp_path, rp_bins=2, rt_bins=3,
                            rp_max=20., rt_max=30., wedge_bins=3, **kwargs)


def eboss_file(path, *, duplicate=False, floating_id=False):
    from astropy.io import fits
    hdus = [fits.PrimaryHDU()]
    for i in range(3):
        hdu = fits.BinTableHDU.from_columns([
            fits.Column(name="LOGLAM", format="D", array=np.log10([3600., 3610., 3620., 3630.])),
            fits.Column(name="DELTA", format="D", array=[.1+.02*i, -.2-.01*i, .3+.025*i, .4]),
            fits.Column(name="WEIGHT", format="D", array=[1., 2., 3., 0.]),
        ])
        identifier = 2**55+1+(0 if duplicate else i)
        hdu.header.update(THING_ID=float(identifier) if floating_id else identifier,
                          RA=.01*i, DEC=.003*i)
        hdus.append(hdu)
    fits.HDUList(hdus).writeto(path)
    return path


def test_eboss_auto_input_int64_cosmology_and_cache(tmp_path):
    from astropy.cosmology import FlatLambdaCDM
    path = eboss_file(tmp_path/"eboss.fits.gz")
    cat = lya_fits.read_fits([str(path)])
    assert cat.metadata["input_format"] == "eBOSS/PICCA forest HDUs"
    assert cat.nbody == 9
    np.testing.assert_array_equal(np.unique(cat.forest_ids), np.arange(3, dtype=np.int64)+2**55+1)
    np.testing.assert_allclose(cat.delta[:3], [.1, -.2, .3])
    np.testing.assert_allclose(cat.weights[:3], [1., 2., 3.])
    expected = FlatLambdaCDM(H0=67.4, Om0=.315, Tcmb0=0).comoving_distance(
        np.array([3600., 3610., 3620.])/1215.67-1).value*.674
    np.testing.assert_allclose(np.linalg.norm(cat.positions[:3], axis=1), expected)
    driver.save_catalog(tmp_path/"cache.npz", cat)
    cached = driver.read_npz(tmp_path/"cache.npz")
    np.testing.assert_array_equal(cached.forest_ids, cat.forest_ids)
    assert cached.metadata == cat.metadata


@pytest.mark.parametrize("option,match", [("duplicate", "duplicate forest ID"), ("floating_id", "integer THING_ID")])
def test_eboss_rejects_ambiguous_ids(tmp_path, option, match):
    path = eboss_file(tmp_path/"bad.fits", **{option: True})
    with pytest.raises(ValueError, match=match):
        lya_fits.read_eboss([str(path)])


def test_explicit_projection_and_redshift_weighting(tmp_path):
    path = eboss_file(tmp_path/"eboss.fits")
    original = lya_fits.read_eboss([str(path)], max_forests=1)
    projected = lya_fits.read_eboss([str(path)], max_forests=1,
                                  project_delta=True, redshift_weight_exponent=1.9)
    loglam = np.log10([3600., 3610., 3620.])
    factor = (10**loglam/1215.67/3.25)**1.9
    np.testing.assert_allclose(projected.weights, original.weights*factor)
    assert abs(np.dot(projected.delta, projected.weights)) < 1e-13
    assert abs(np.dot(projected.delta*projected.weights, loglam)) < 1e-13
    assert projected.metadata["project_delta"]
    with pytest.raises(ValueError, match="no valid forest pixels"):
        lya_fits.read_eboss([str(path)], z_min=3., z_max=4.)


def test_projection_degenerate_forest():
    d, w = lya_fits.preprocess_forest([.2], [1.], [3.6], project_delta=True)
    np.testing.assert_allclose(d, 0.)
    np.testing.assert_allclose(w, 1.)


def test_all_pair_comparisons_and_near_zero_policy(tmp_path):
    cfg = config(tmp_path, plots=False, relative_floor=1e-6)
    table = np.array([[0, 1., 1e-9, 2e-9, 2.]])
    results = {name: dict(products={"1d_2pcf": table.copy()}) for name in ("A", "B", "C")}
    metrics, paths = analysis.write_comparisons(results, cfg)
    assert len(metrics) == 3 and not paths
    assert all(m["passed"] and m["max_abs_relative"] is None for m in metrics.values())
    bad = table.copy()
    bad[0, -1] = 0
    m, _ = driver.compare_products(table, bad, "1d_2pcf")
    assert not m["passed"] and m["occupancy_mismatches"] == 1
    with pytest.raises(ValueError, match="duplicate histogram"):
        driver.compare_products(np.vstack((table, table)), table, "1d_2pcf")


def test_distortion_orientation_metadata_and_validation(tmp_path):
    cfg = config(tmp_path)
    matrix = np.eye(6)
    matrix[0, 1] = .25
    model = np.arange(6.)
    np.testing.assert_allclose(analysis.apply_distortion(matrix, model), [0.25, 1, 2, 3, 4, 5])
    np.savez(tmp_path/"matrix.npz", distortion=matrix, rp_bins=2, rt_bins=3,
             rp_max=20., rt_max=30.)
    np.testing.assert_array_equal(analysis.load_analysis_array(tmp_path/"matrix.npz", "distortion", cfg), matrix)
    np.savez(tmp_path/"wrong.npz", distortion=matrix, rp_max=25.)
    with pytest.raises(ValueError, match="does not match"):
        analysis.load_analysis_array(tmp_path/"wrong.npz", "distortion", cfg)
    np.save(tmp_path/"wrong.npy", matrix[:3])
    with pytest.raises(ValueError, match="shape"):
        analysis.load_analysis_array(tmp_path/"wrong.npy", "distortion", cfg)
    matrix[0, :] = np.nan
    np.save(tmp_path/"empty.npy", matrix)
    assert np.isnan(analysis.load_analysis_array(tmp_path/"empty.npy", "distortion", cfg)[0]).all()
    matrix[0, 0] = 1
    np.save(tmp_path/"partial.npy", matrix)
    with pytest.raises(ValueError, match="partly non-finite"):
        analysis.load_analysis_array(tmp_path/"partial.npy", "distortion", cfg)


def test_fits_analysis_inputs(tmp_path):
    from astropy.io import fits
    cfg = config(tmp_path)
    hdu = fits.BinTableHDU.from_columns([
        fits.Column(name="DA", format="D", array=np.arange(6.)),
        fits.Column(name="DM", format="6D", array=np.eye(6)),
        fits.Column(name="CO", format="6D", array=np.eye(6)*.1),
    ], name="COR")
    hdu.header.update(NP=2, NT=3, RPMAX=20., RTMAX=30., RPMIN=0.)
    hdu.writeto(tmp_path/"cor.fits")
    for kind, expected in (("model", np.arange(6.)), ("distortion", np.eye(6)), ("covariance", np.eye(6)*.1)):
        np.testing.assert_array_equal(analysis.load_analysis_array(tmp_path/"cor.fits", kind, cfg), expected)


def test_wedge_covariance_and_empty_bins():
    xi, den = np.array([1., 3., np.nan]), np.array([1., 2., 0.])
    covariance = np.diag([1., 4., 0.])
    geometry = np.array([[1., 1., 1.], [0., 0., 1.]])
    value, error, projection, wedge_cov = analysis.project_wedges(xi, den, geometry, covariance)
    np.testing.assert_allclose(projection[0], [.8, .2, 0.])
    np.testing.assert_allclose(value[0], 1.4)
    np.testing.assert_allclose(error[0]**2, .8)
    np.testing.assert_allclose(wedge_cov, projection @ covariance @ projection.T)
    assert np.isnan(value[1]) and np.isnan(error[1])
    value, error, _, _ = analysis.project_wedges(xi, den, geometry)
    np.testing.assert_allclose(value[0], 7/3)
    assert np.isnan(error).all()
    with pytest.raises(ValueError, match="negative wedge variance"):
        analysis.project_wedges([1., 2.], [1., 1.], np.ones((1, 2)), np.array([[1., -2.], [-2., 1.]]))


def test_analysis_exports_and_forward_model(tmp_path):
    cfg = config(tmp_path, plots=True)
    model = np.arange(6.)*.001
    np.save(tmp_path/"model.npy", model)
    np.save(tmp_path/"distortion.npy", np.eye(6))
    cfg.model_correlation, cfg.distortion_matrix = tmp_path/"model.npy", tmp_path/"distortion.npy"
    table = product_table(model.reshape(2, 3), np.ones((2, 3)), vars(cfg), "3d_2pcf")
    report, paths = analysis.analyse_2pcf({"test": dict(products={"3d_2pcf": table})}, cfg)
    assert len(paths) == 3 and all(Path(p).is_file() for p in paths)
    assert report["model_residuals"]["test"]["max_abs"] == 0
    with np.load(report["archives"]["test"]) as data:
        np.testing.assert_allclose(data["model_residual"], 0.)
        assert np.isnan(data["wedge_error"]).all()
    cfg.max_hist_mib = .00001
    with pytest.raises(ValueError, match="workspace exceeds"):
        analysis.analyse_2pcf({"test": dict(products={"3d_2pcf": table})}, cfg)


def test_invalid_analysis_options(tmp_path):
    with pytest.raises(ValueError, match="together"):
        config(tmp_path, distortion_matrix=tmp_path/"missing.npy").validate()
    with pytest.raises(ValueError, match="covariance requires"):
        config(tmp_path, reference_covariance=True).validate()
    with pytest.raises(ValueError, match="increase strictly"):
        config(tmp_path, wedge_mu_edges=(0., 1., .5)).validate()


def test_multipole_map_is_monopole_not_sum(tmp_path):
    cfg = config(tmp_path)
    table = np.array([[0, 1, 2, 1., 2., .5, 1., 2.],
                      [1, 1, 2, 1., 2., 5., 10., 2.]])
    assert driver.projected_plot_data(table, "multipole_3pcf", cfg)[0, 1] == .5


REFERENCE = os.environ.get("LYA2PCF_SOURCE")


@pytest.mark.skipif(not REFERENCE, reason="set LYA2PCF_SOURCE for external kernel integration")
def test_original_kernel_all_estimator_families_and_threads(tmp_path):
    # Equal-RA distinct forests exercise the upstream neighborhood tie omission.
    ra, dec = np.array([0., 0., .025, .04]), np.array([0., .03, .01, -.02])
    los = np.column_stack((np.cos(dec)*np.cos(ra), np.cos(dec)*np.sin(ra), np.sin(dec)))
    radii = np.array([[100., 100., 107.3], [103.1, 104.7, 108.],
                      [99.3, 103.8, 109.2], [101.2, 108.3, 110.1]])
    pos = (radii[..., None]*los[:, None]).reshape(-1, 3)
    rng = np.random.default_rng(123)
    weights = rng.uniform(.5, 3, len(pos))
    weights[4] = 0
    cat = driver.ForestCatalog(pos, rng.normal(size=len(pos)), weights,
                              np.repeat(np.arange(4, dtype=np.int64)+2**55+1, 3)).normalized()
    families = ("3d_2pcf", "1d_2pcf", "1d_same_los_2pcf")
    expected = {}
    for family in families:
        shape = (2, 3) if family == "3d_2pcf" else (2,)
        num, den, per_forest = np.zeros(shape), np.zeros(shape), {}
        for i, j in combinations(range(len(pos)), 2):
            same = cat.forest_ids[i] == cat.forest_ids[j]
            if same != (family == "1d_same_los_2pcf"):
                continue
            r1, r2 = radii.ravel()[[i, j]]
            lag = abs(r1-r2)
            if family == "3d_2pcf":
                cosine = np.dot(pos[i]/r1, pos[j]/r2)
                rp = lag*np.sqrt((1+cosine)/2)
                rt = (r1+r2)*np.sqrt(max(0., (1-cosine)/2))
                bins = (int(rp/10), int(rt/10))
                if rp >= 20 or rt >= 30:
                    continue
            else:
                bins = (int(lag/10),)
                if lag >= 20:
                    continue
            w = weights[i]*weights[j]
            if same:
                n, d = per_forest.setdefault(int(cat.forest_ids[i]), (np.zeros(shape), np.zeros(shape)))
                n[bins] += w*cat.delta[i]*cat.delta[j]
                d[bins] += w
            else:
                num[bins] += w*cat.delta[i]*cat.delta[j]
                den[bins] += w
        for n, d in per_forest.values():
            num += np.divide(n, d, out=np.zeros_like(n), where=d > 0)
            den += d > 0
        expected[family] = product_table(num, den, vars(config(tmp_path)), family)
    previous = None
    for threads in (1, 2):
        out = tmp_path/str(threads)
        out.mkdir()
        cfg = config(out, lya2pcf_source=Path(REFERENCE), threads=threads,
                     reference_covariance=True, reference_nside=128, plots=False)
        result = run_references(cat, cfg, families)
        for row in result.values():
            for family, table in row["products"].items():
                np.testing.assert_allclose(table, expected[family], rtol=1e-12, atol=1e-12)
                if previous is not None:
                    np.testing.assert_array_equal(table, previous[family])
        previous = {key: table for row in result.values() for key, table in row["products"].items()}
        with np.load(result["lya2pcf-cpu"]["reference_archive"]) as data:
            w = data["partial_denominator"].reshape(-1, 6)
            n = data["partial_numerator"].reshape(-1, 6)
            total = w.sum(axis=0)
            means = np.divide(n.sum(axis=0), total, out=np.zeros(6), where=total > 0)
            residual = n-w*means
            covariance = np.divide(residual.T@residual, total[:, None]*total[None, :],
                                   out=np.zeros((6, 6)), where=total[:, None]*total[None, :] > 0)
            np.testing.assert_allclose(data["covariance"], covariance, rtol=1e-12, atol=1e-12)


@pytest.mark.skipif(not REFERENCE, reason="set LYA2PCF_SOURCE for notebook formula comparison")
def test_wedges_match_reference_notebook(tmp_path):
    spec = importlib.util.spec_from_file_location("upstream_wedges", Path(REFERENCE)/"plot_auxiliars.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    # The old helper uses scipy.diagonal, removed from current SciPy.
    module.sp = np
    cfg = config(tmp_path, wedge_subsamples=100, wedge_mu_edges=(0., .5, 1.))
    geometry, radii = analysis.wedge_geometry(cfg)
    xi, cov = np.arange(6.)*.01, np.diag(np.arange(6.)+1)
    value, _, _, projected_cov = analysis.project_wedges(xi, np.ones(6), geometry, cov)
    reference_geometry = []
    for lo, hi in zip(cfg.wedge_mu_edges[:-1], cfg.wedge_mu_edges[1:]):
        for r in range(cfg.wedge_bins):
            with np.errstate(invalid="ignore"):
                g = module.weight_2d_wedge(lo, hi, r*20/3, (r+1)*20/3,
                                          (2, 3), 20., 30.)
            reference_geometry.append(np.nan_to_num(g).ravel())
    ref_value, ref_cov = module.wedge(np.array(reference_geometry), xi, cov)
    active = np.isfinite(value)
    np.testing.assert_allclose(value[active], ref_value[active], rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(projected_cov, ref_cov, rtol=1e-12, atol=1e-12)
