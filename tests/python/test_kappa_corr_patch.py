#!/usr/bin/env python3
"""Angular patch selection shared by all convergence engines."""

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tests" / "python"))

import kappa_corr_all_engines as driver


def test_patch_bounds_and_native_open_boundary_convention():
    patch = driver.AngularPatch(20, 80, 30, 100)
    theta = np.radians([30, 100, 60, 60, 60, 60, 60])
    phi = np.radians([50, 50, 20, 80, 50, 10, 90])
    np.testing.assert_array_equal(
        patch.select(theta, phi), [False, False, False, False, True, False, False],
    )
    full = driver.AngularPatch(0, 360, 0, 180)
    assert not full.select(np.array([np.pi / 2]), np.array([0]))[0]


@pytest.mark.parametrize("bounds", [
    (-1, 90, 0, 90), (0, 361, 0, 90), (350, 10, 0, 90),
    (20, 20, 0, 90), (0, 90, -1, 90), (0, 90, 0, 181),
    (0, 90, 90, 90), (0, 90, 100, 10),
    (float("nan"), 90, 0, 90), (0, float("inf"), 0, 90),
])
def test_patch_rejects_invalid_bounds(bounds):
    with pytest.raises(ValueError, match="patch"):
        driver.AngularPatch(*bounds)


def test_patch_is_opt_in_and_consumed_before_native_engine(tmp_path):
    engine = "octree-2balls-omp"
    config = driver.RunConfig((engine,), tmp_path)
    assert config.normalized().angular_patch is None
    assert not driver.parse_arguments([]).patch
    assert driver.parse_arguments(["--patch"]).patch
    for updates in ({"patch": True}, {"options": ("patch,no-smooth-pivot",)}):
        config = driver.RunConfig((engine,), tmp_path, **updates).normalized()
        assert config.angular_patch == driver.AngularPatch()
        params = driver.engine_parameters(config, engine, False)
        assert "patch" not in params["options"].split(",")
    with pytest.raises(ValueError, match="patch-with-all"):
        driver.RunConfig((engine,), tmp_path, patch=True,
                         options=("patch-with-all",)).normalized()


@pytest.mark.parametrize("nested", [False, True])
@pytest.mark.parametrize("masked", [False, True])
@pytest.mark.parametrize("max_points", [0, 25])
def test_native_fits_patch_precedes_sampling_and_centering(
    tmp_path, nested, masked, max_points,
):
    hp = pytest.importorskip("healpy")
    nside = 16
    pixels = np.arange(hp.nside2npix(nside))
    values = pixels.astype(float) + 0.125
    values[[100, 900]] = hp.UNSEEN
    values[1300] = np.nan
    mask = (pixels % 3 != 0).astype(float)
    mask[500] = hp.UNSEEN
    path = tmp_path / "map.fits"
    mask_path = tmp_path / "mask.fits" if masked else None
    hp.write_map(path, values, nest=nested, dtype=np.float64, overwrite=True)
    if masked:
        hp.write_map(mask_path, mask, nest=nested, dtype=np.float64, overwrite=True)
    patch = driver.AngularPatch(25, 290, 35, 130)
    theta, phi = hp.pix2ang(nside, pixels, nest=nested)
    eligible = driver._healpix_valid(values, hp)
    if masked:
        eligible &= driver._healpix_valid(mask, hp) & (mask > 0)
    before_patch = int(eligible.sum())
    expected = pixels[eligible & patch.select(theta, phi)]
    before_thinning = expected.size
    if max_points:
        order = np.argsort(driver._splitmix64(expected, 123))[:max_points]
        expected = np.sort(expected[order])
    reference = None
    for chunk in (1024, 4096):
        catalog = driver.catalog_from_healpix(
            path, mask_path=mask_path, patch=patch, max_points=max_points,
            sampling_seed=123, chunk_pixels=chunk, center_field=False,
        )
        np.testing.assert_array_equal(catalog.kappa, values[expected])
        np.testing.assert_allclose(
            catalog.positions, np.column_stack(hp.pix2vec(nside, expected, nest=nested)),
        )
        assert catalog.metadata["eligible_pixels_before_patch"] == before_patch
        assert catalog.metadata["eligible_pixels_before_thinning"] == before_thinning
        assert catalog.metadata["patch"] == patch.metadata()
        assert catalog.metadata["loader"] == "memory-mapped-native-resolution"
        if masked:
            np.testing.assert_array_equal(catalog.mask, 1)
        else:
            assert catalog.mask is None
        if reference is not None:
            np.testing.assert_array_equal(catalog.kappa, reference.kappa)
        reference = catalog
    centered = driver.catalog_from_healpix(
        path, mask_path=mask_path, patch=patch, max_points=max_points,
        sampling_seed=123, center_field=True,
    )
    np.testing.assert_allclose(centered.kappa, values[expected] - values[expected].mean())


@pytest.mark.parametrize("fallback", ["downgrade", "mask-ordering", "explicit"])
def test_healpy_fallback_patch(tmp_path, fallback):
    hp = pytest.importorskip("healpy")
    nside = 16
    pixels = np.arange(hp.nside2npix(nside))
    values = 0.1 + np.sin(pixels)
    mask = (pixels % 3 != 0).astype(float)
    path, mask_path = tmp_path / "map.fits", tmp_path / "mask.fits"
    if fallback == "explicit":
        values[::7] = hp.UNSEEN
    hp.write_map(path, values, partial=fallback == "explicit", dtype=np.float64)
    if fallback == "mask-ordering":
        hp.write_map(mask_path, hp.reorder(mask, r2n=True), nest=True, dtype=np.float64)
    else:
        hp.write_map(mask_path, mask, dtype=np.float64)
    target = 8 if fallback == "downgrade" else nside
    expected_values = hp.read_map(path, dtype=np.float64)
    expected_mask = hp.read_map(mask_path, dtype=np.float64)
    if target != nside:
        expected_values = hp.ud_grade(expected_values, target)
        expected_mask = hp.ud_grade(expected_mask, target)
    patch = driver.AngularPatch(45, 220, 20, 150)
    theta, phi = hp.pix2ang(target, np.arange(expected_values.size))
    selected = (driver._healpix_valid(expected_values, hp)
                & driver._healpix_valid(expected_mask, hp) & (expected_mask > 0.5)
                & patch.select(theta, phi))
    catalog = driver.catalog_from_healpix(
        path, nside_down=target, mask_path=mask_path, mask_threshold=0.5,
        center_field=False, patch=patch,
    )
    np.testing.assert_array_equal(catalog.kappa, expected_values[selected])
    assert catalog.metadata["loader"].startswith("healpy-")
    assert catalog.metadata["patch"] == patch.metadata()


def test_npz_patch_preserves_rows_weights_masks_and_centers_selected_field(tmp_path):
    hp = pytest.importorskip("healpy")
    positions = np.column_stack(hp.pix2vec(4, np.arange(hp.nside2npix(4))))
    values = np.arange(len(positions), dtype=float)
    weights = values + 1
    mask = (values % 3 != 0).astype(np.uint8)
    path = tmp_path / "catalog.npz"
    np.savez(path, positions=positions, kappa=values, weights=weights, mask=mask)
    patch = driver.AngularPatch(20, 140, 15, 120)
    theta, phi = hp.vec2ang(positions)
    selected = patch.select(theta, phi)
    catalog = driver.catalog_from_npz(path, patch=patch, center_field=True)
    np.testing.assert_array_equal(catalog.positions, positions[selected])
    np.testing.assert_array_equal(catalog.weights, weights[selected])
    np.testing.assert_array_equal(catalog.mask, mask[selected])
    np.testing.assert_allclose(
        catalog.kappa, values[selected] - values[selected & mask.astype(bool)].mean(),
    )
    assert driver.filter_catalog_patch(catalog, patch) is catalog
    original = driver.KappaCatalog(positions, values - values.mean(), weights, mask,
                                   {"centered": True})
    filtered = driver.filter_catalog_patch(original, patch)
    np.testing.assert_allclose(filtered.kappa, catalog.kappa)
    np.testing.assert_array_equal(original.kappa, values - values.mean())


def test_empty_patch_reports_clear_error(tmp_path):
    hp = pytest.importorskip("healpy")
    patch = driver.AngularPatch(1, 1.001, 1, 1.001)
    with pytest.raises(ValueError, match="fewer than three"):
        driver.synthetic_healpix_catalog(4, patch=patch)
    path = tmp_path / "map.fits"
    hp.write_map(path, np.ones(hp.nside2npix(4)), dtype=np.float64)
    with pytest.raises(ValueError, match="fewer than three"):
        driver.catalog_from_healpix(path, patch=patch)
    catalog = driver.synthetic_healpix_catalog(4)
    with pytest.raises(ValueError, match="fewer than three"):
        driver.filter_catalog_patch(catalog, patch)
    catalog.positions[0] = 0
    with pytest.raises(ValueError, match="nonzero"):
        driver.filter_catalog_patch(catalog, driver.AngularPatch())


@pytest.mark.parametrize("flag", [["--patch"], ["--more-options", "patch"]])
def test_main_passes_one_filtered_catalog_to_suite_and_npz(tmp_path, monkeypatch, flag):
    pytest.importorskip("healpy")
    engine = "octree-2balls-omp"
    monkeypatch.setattr(driver, "discover_search_methods", lambda path: [engine])
    monkeypatch.setattr(driver, "discover_make_settings", lambda path: {})
    monkeypatch.setattr(driver, "discover_cython_methods", lambda names: [engine])
    if hasattr(driver, "inspect_cballs_runtime"):
        monkeypatch.setattr(driver, "inspect_cballs_runtime", lambda path: {
            "search_methods": [engine], "make_settings": {},
        })
    monkeypatch.setattr(driver, "get_mpi_comm", lambda required: driver.SerialComm())
    calls = []
    monkeypatch.setattr(driver, "run_engine_suite",
                        lambda catalog, config, **kwargs: calls.append((catalog, config)))
    output = tmp_path / "patch.npz"
    monkeypatch.setattr(sys, "argv", [
        "driver", "--synthetic-nside", "8", "--engine", engine,
        "--phiL", "25", "--phiR", "130", "--thetaL", "30", "--thetaR", "110",
        "--save-catalog-npz", str(output), "--outdir", str(tmp_path / "out"),
        "--no-plots", *flag,
    ])
    assert driver.main() == 0
    assert len(calls) == 1
    catalog, config = calls[0]
    expected = driver.synthetic_healpix_catalog(8, patch=driver.AngularPatch(25, 130, 30, 110))
    np.testing.assert_array_equal(catalog.positions, expected.positions)
    np.testing.assert_array_equal(catalog.kappa, expected.kappa)
    assert config.patch
    with np.load(output) as saved:
        np.testing.assert_array_equal(saved["positions"], catalog.positions)
        np.testing.assert_array_equal(saved["kappa"], catalog.kappa)


def test_no_patch_keeps_full_catalog(tmp_path):
    hp = pytest.importorskip("healpy")
    values = np.arange(hp.nside2npix(4), dtype=float)
    path = tmp_path / "map.fits"
    hp.write_map(path, values, dtype=np.float64)
    catalog = driver.catalog_from_healpix(path, center_field=False)
    np.testing.assert_array_equal(catalog.kappa, values)
    assert catalog.metadata["patch"] is None
    synthetic = driver.synthetic_healpix_catalog(4)
    assert synthetic.nbody == len(values)
