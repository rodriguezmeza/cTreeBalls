"""RA-DEC-R vector forests: geometry, identity, chunk joins and selection."""
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest

fits = pytest.importorskip("astropy.io.fits")
import lya_corr_all_engines as driver
from lya_fits import read_fits, read_radec_r


def catalog_file(path, *, angle_unit=None, radial_unit=None, ids=None, id_format="D",
                 blind=False, inconsistent=False, duplicate=False):
    ids = np.array([11, 12, 13] if ids is None else ids)
    hdus = [fits.PrimaryHDU()]
    for chunk, rows in enumerate(([0, 1, 2], [0, 2])):
        z = np.array([2.1, 2.15, 2.2, 2.25] if chunk == 0 else [2.3, 2.35, 2.4])
        if duplicate and chunk: z[0] = 2.1
        n, m = len(rows), len(z)
        ra, dec = np.array([0., 90., 0.])[rows], np.array([0., 0., 90.])[rows]
        if inconsistent and chunk: ra[0] = 1.
        if angle_unit == "rad": ra, dec = np.deg2rad(ra), np.deg2rad(dec)
        weight = np.ones((n, m))
        if not chunk: weight[:, 1] = 0
        delta = np.tile(np.arange(1, m+1)*.1, (n, 1))+np.array(rows)[:, None]
        columns = [
            fits.Column(name="ra", format="D", unit=angle_unit, array=ra),
            fits.Column(name="dec", format="D", unit=angle_unit, array=dec),
            fits.Column(name="tid", format=id_format, array=ids[rows]),
            fits.Column(name="cid", format="K", array=np.full(n, chunk+1)),
            fits.Column(name="rcdis", format=f"{m}D", unit=radial_unit,
                        array=np.tile(np.arange(1, m+1)*10+chunk*40, (n, 1))),
            fits.Column(name="delta", format=f"{m}D", array=delta),
            fits.Column(name="wdelta", format=f"{m}D", array=weight),
            fits.Column(name="wobs", format=f"{m}D", array=np.tile(1215.67*(1+z), (n, 1))),
            fits.Column(name="zobs", format=f"{m}D", array=np.tile(z, (n, 1))),
        ]
        if blind:
            columns.append(fits.Column(name="DELTA_BLIND", format=f"{m}D", array=delta+10))
        hdus.append(fits.BinTableHDU.from_columns(columns, name=f"SPECTRA_{2.2+chunk*.2:.1f}"))
    fits.HDUList(hdus).writeto(path)
    return path


def test_auto_geometry_chunk_identity_weights_and_cache(tmp_path):
    path = catalog_file(tmp_path/"spectra.fits")
    cat = read_fits([str(path)])
    assert cat.metadata["input_format"] == "RA-DEC-R forest FITS"
    assert cat.nbody == 15
    assert cat.forest_ids.dtype == np.int64
    np.testing.assert_array_equal(cat.forest_ids, np.repeat([11, 12, 13], [6, 3, 6]))
    radius = np.array([10, 30, 40, 50, 60, 70])
    np.testing.assert_allclose(cat.positions[:6], np.column_stack([radius, radius*0, radius*0]))
    np.testing.assert_allclose(cat.positions[6:9], np.column_stack([radius[:3]*0, radius[:3], radius[:3]*0]), atol=1e-13)
    np.testing.assert_allclose(cat.positions[9:], np.column_stack([radius*0, radius*0, radius]), atol=1e-13)
    np.testing.assert_array_equal(cat.weights, 1)
    np.testing.assert_allclose(cat.delta[:6], [.1, .3, .4, .1, .2, .3])
    np.testing.assert_allclose(10**cat.log_lambda[:6], 1215.67*np.array([3.1, 3.2, 3.25, 3.3, 3.35, 3.4]))
    driver.save_catalog(tmp_path/"cache.npz", cat)
    loaded = driver.read_npz(tmp_path/"cache.npz")
    np.testing.assert_array_equal(loaded.positions, cat.positions)
    np.testing.assert_array_equal(loaded.forest_ids, cat.forest_ids)
    assert loaded.metadata == cat.metadata


def test_max_forests_and_stride_apply_to_joined_forest(tmp_path):
    path = catalog_file(tmp_path/"spectra.fits.gz")
    cat = read_fits([str(path)], max_forests=1, pixel_stride=2, retain_projection_domain=True)
    np.testing.assert_array_equal(cat.forest_ids, 11)
    np.testing.assert_array_equal(cat.positions[:, 0], [10, 40, 60])
    np.testing.assert_array_equal(cat.distortion_domain["positions"][:, 0], [10, 30, 40, 50, 60, 70])
    np.testing.assert_array_equal(cat.distortion_domain["data_indices"], [0, 2, 4])
    # Redshift selection starts in the SECOND HDU; keep the earlier fit domain.
    later = read_fits([str(path)], max_forests=2, z_min=2.3, pixel_stride=2,
                      retain_projection_domain=True)
    np.testing.assert_array_equal(later.forest_ids, [11, 11, 13, 13])
    np.testing.assert_allclose(np.linalg.norm(later.positions, axis=1), [50, 70, 50, 70])
    assert len(later.distortion_domain["positions"]) == 12
    np.testing.assert_array_equal(later.distortion_domain["data_indices"], [3, 5, 9, 11])


@pytest.mark.parametrize("declared", [False, True])
def test_mpc_conversion_and_radian_angles(tmp_path, declared):
    path = catalog_file(tmp_path/"spectra.fits", angle_unit="rad" if declared else None,
                        radial_unit="Mpc" if declared else None)
    if not declared:
        with fits.open(path, mode="update") as hdus:
            for hdu in hdus[1:]:
                hdu.data["ra"][:] = np.deg2rad(hdu.data["ra"])
                hdu.data["dec"][:] = np.deg2rad(hdu.data["dec"])
    options = {} if declared else dict(radec_r_angle_unit="rad", radec_r_distance_unit="Mpc")
    cat = read_fits([str(path)], h=.7, **options)
    np.testing.assert_allclose(np.linalg.norm(cat.positions[:6], axis=1), .7*np.array([10, 30, 40, 50, 60, 70]))
    np.testing.assert_allclose(cat.positions[6:9, 1], .7*np.array([10, 30, 40]))


def test_joint_padding_mask(tmp_path):
    path = catalog_file(tmp_path/"spectra.fits")
    with fits.open(path, mode="update") as hdus:
        hdus[1].data["delta"][0, 0] = np.nan
        hdus[1].data["rcdis"][0, 2] = np.nan
        hdus[1].data["wdelta"][0, 3] = -1
    cat = read_fits([str(path)])
    # Only the valid later chunk remains; independent column NaN removal
    # would instead misalign delta, radius, wavelength and weight.
    np.testing.assert_array_equal(cat.positions[cat.forest_ids == 11, 0], [50, 60, 70])
    np.testing.assert_allclose(cat.delta[cat.forest_ids == 11], [.1, .2, .3])


def test_projection_and_redshift_weights_use_joined_wavelengths(tmp_path):
    path = catalog_file(tmp_path/"spectra.fits")
    raw = read_fits([str(path)], max_forests=1)
    cat = read_fits([str(path)], max_forests=1, project_delta=True,
                    redshift_weight_exponent=1.5, retain_projection_domain=True)
    expected_weight = raw.weights*(10**raw.log_lambda/1215.67/3.25)**1.5
    design = np.column_stack([np.ones(raw.nbody), raw.log_lambda-raw.log_lambda[0]])
    coeff = np.linalg.lstsq(design*np.sqrt(expected_weight[:, None]),
                           raw.delta*np.sqrt(expected_weight), rcond=None)[0]
    np.testing.assert_allclose(cat.weights, expected_weight)
    np.testing.assert_allclose(cat.delta, raw.delta-design@coeff, atol=1e-13)
    assert cat.distortion_domain is None  # fitting is on the retained pixels


def test_blind_field_and_int64_identifiers(tmp_path):
    ids = np.arange(3, dtype=np.int64)+2**55+1
    path = catalog_file(tmp_path/"spectra.fits", ids=ids, id_format="K", blind=True)
    cat = read_fits([str(path)])
    np.testing.assert_array_equal(np.unique(cat.forest_ids), ids)
    assert cat.delta.min() > 10
    raw = read_fits([str(path)], delta_field="DELTA")
    np.testing.assert_allclose(cat.delta, raw.delta+10)


@pytest.mark.parametrize("identifier", [11.5, np.nan, float(2**53)])
def test_unsafe_float_ids_fail(tmp_path, identifier):
    path = catalog_file(tmp_path/"spectra.fits", ids=[identifier, 12, 13])
    with pytest.raises(ValueError, match="tid.*exact integers"):
        read_fits([str(path)])


@pytest.mark.parametrize("fault", ["inconsistent", "duplicate"])
def test_conflicting_chunks_fail(tmp_path, fault):
    path = catalog_file(tmp_path/"spectra.fits", **{fault: True})
    with pytest.raises(ValueError, match="inconsistent sky|duplicate wavelength"):
        read_fits([str(path)])


def test_options_validation_and_explicit_layout(tmp_path):
    path = catalog_file(tmp_path/"spectra.fits")
    args = driver.parse_arguments(["--fits", str(path), "--fits-layout", "radec-r",
        "--radec-r-angle-unit", "deg", "--radec-r-distance-unit", "Mpc", "--h", ".7"])
    cat = driver.load_catalog(args)
    np.testing.assert_allclose(cat.positions[0], [7, 0, 0])
    for options, message in [(dict(pixel_stride=0), "positive"), (dict(max_forests=0), "positive"),
                              (dict(z_min=5, z_max=2), "z_min"), (dict(z_min=5), "no valid"),
                              (dict(omega_m=.3), "supplies rcdis")]:
        with pytest.raises(ValueError, match=message):
            read_radec_r([str(path)], **options)
    with pytest.raises(ValueError, match="FITS.*require --fits"):
        driver.load_catalog(driver.parse_arguments(["--synthetic", "--radec-r-distance-unit", "Mpc"]))


def test_driver_radec_r_and_cached_cartesian_results_agree(tmp_path):
    pytest.importorskip("cyballs")
    path = catalog_file(tmp_path/"spectra.fits")
    cat = read_fits([str(path)])
    cache = tmp_path/"cartesian.npz"
    driver.save_catalog(cache, cat)
    products = []
    for name, source in [("fits", ["--fits", str(path)]), ("cache", ["--catalog", str(cache)])]:
        output = tmp_path/name
        result = subprocess.run([sys.executable, driver.__file__, *source,
            "--engine", "lya-2pcf-omp", "lya-los-tree-2pcf-omp", "--threads", "2",
            "--lya2-kernel", "1", "--no-analysis", "--no-plots", "--fail-on-mismatch",
            "--output", str(output)], text=True, capture_output=True, timeout=90)
        assert result.returncode == 0, result.stdout+result.stderr
        summary = json.loads((output/"summary.json").read_text())
        assert summary["validation"]["passed"]
        products.append(np.loadtxt(output/"lya-2pcf-omp"/"histXi2pcf_lya.txt"))
    np.testing.assert_array_equal(products[0], products[1])


def test_redshift_mode_uses_observer_distances_and_shared_domain(tmp_path):
    from astropy.cosmology import FlatLambdaCDM
    path = catalog_file(tmp_path/"spectra.fits")
    cat = read_fits([str(path)], radec_r_distance_source="redshift", omega_m=.3, h=.7,
                    max_forests=2, z_min=2.3, pixel_stride=2, retain_projection_domain=True)
    cosmo = FlatLambdaCDM(H0=70, Om0=.3, Tcmb0=0)
    z = np.array([2.1, 2.2, 2.25, 2.3, 2.35, 2.4])
    chi = cosmo.comoving_distance(z).value*.7
    np.testing.assert_allclose(np.linalg.norm(cat.positions,axis=1), np.tile(chi[[3, 5]],2), rtol=2e-14)
    np.testing.assert_allclose(np.linalg.norm(cat.distortion_domain["positions"],axis=1), np.tile(chi,2), rtol=2e-14)
    np.testing.assert_array_equal(cat.distortion_domain["data_indices"], [3, 5, 9, 11])
    assert cat.metadata["distance_source_mode"] == "redshift"
    assert cat.metadata["cosmology"] == dict(model="FlatLambdaCDM", omega_m=.3, h=.7, Tcmb0=0)
    assert cat.metadata["files"][0]["units"][0]["rcdis_unit_is_assumed"]
    driver.save_catalog(tmp_path/"redshift.npz", cat)
    cached = driver.read_npz(tmp_path/"redshift.npz")
    np.testing.assert_array_equal(cached.positions, cat.positions)
    np.testing.assert_array_equal(cached.distortion_domain["positions"], cat.distortion_domain["positions"])


def test_redshift_mode_ignores_unused_rcdis_but_checks_absorber_redshift(tmp_path):
    path = catalog_file(tmp_path/"spectra.fits", radial_unit="km")
    with fits.open(path, mode="update") as hdus:
        hdus[1].data["rcdis"][:] = np.nan
    cat = read_radec_r([str(path)], distance_source="redshift", max_forests=1)
    assert cat.nbody == 6  # All otherwise valid pixels survive despite meaningless rcdis.
    with fits.open(path, mode="update") as hdus:
        hdus[1].data["zobs"][0, 0] += .1
    with pytest.raises(ValueError, match="zobs disagrees"):
        read_radec_r([str(path)], distance_source="redshift")


def test_redshift_projection_uses_same_final_weights_as_estimator(tmp_path):
    path = catalog_file(tmp_path/"spectra.fits")
    cat = read_radec_r([str(path)], distance_source="redshift", project_delta=True,
                       redshift_weight_exponent=1.9, pixel_stride=2)
    for identifier in np.unique(cat.forest_ids):
        m = cat.forest_ids == identifier
        np.testing.assert_allclose(np.dot(cat.delta[m],cat.weights[m]),0,atol=1e-13)
        np.testing.assert_allclose(np.dot(cat.delta[m]*cat.weights[m],cat.log_lambda[m]),0,atol=1e-13)
    assert cat.distortion_domain is None


def test_new_cli_option_alias_and_field_diagnostics(tmp_path):
    path = catalog_file(tmp_path/"spectra.fits")
    args = driver.parse_arguments(["--fits", str(path), "--fits-layout", "radec-r",
        "--radec-r-distance-source", "redshift", "--omega-m", ".3", "--project-delta",
        "--outdir", str(tmp_path/"out")])
    assert args.output == tmp_path/"out"
    cat = driver.load_catalog(args)
    assert cat.metadata["distance_source_mode"] == "redshift"
    assert abs(cat.metadata["field_diagnostics"]["weighted_mean_delta"]) < 1e-13
    with pytest.raises(ValueError, match="FITS.*require --fits"):
        driver.load_catalog(driver.parse_arguments(["--synthetic", "--radec-r-distance-source", "redshift"]))
    with pytest.raises(ValueError, match="requires RA-DEC-R"):
        read_fits([str(path)], fits_layout="eboss", radec_r_distance_source="redshift")
    for value in [np.nan, -1., 1.]:
        with pytest.raises(ValueError, match="omega_m"):
            read_radec_r([str(path)], distance_source="redshift", omega_m=value)


def test_diagnostics_match_weighted_field_across_large_catalog():
    # Cross more than one diagnostic chunk with nonconstant signed data and
    # a broad weight range. This requires no native histogram workspace.
    cat = driver.synthetic_catalog(forests=3, pixels=100001)
    result = driver.catalog_field_diagnostics(cat)
    assert result["weighted_mean_delta"] == pytest.approx(np.average(cat.delta,weights=cat.weights),abs=1e-16)
    assert result["unweighted_mean_delta"] == pytest.approx(np.mean(cat.delta),abs=1e-16)
    radius = np.linalg.norm(cat.positions,axis=1)
    assert result["min_observer_distance_mpc_h"] == pytest.approx(radius.min())
    assert result["max_observer_distance_mpc_h"] == pytest.approx(radius.max())
