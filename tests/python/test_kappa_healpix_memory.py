"""HEALPix loader equivalence and bounded coordinate scratch regressions."""

import weakref

import numpy as np
import pytest

import kappa_corr_all_engines as driver

hp = pytest.importorskip("healpy")


@pytest.mark.parametrize("nested", [False, True])
@pytest.mark.parametrize("masked", [False, True])
@pytest.mark.parametrize("patch_enabled", [False, True])
@pytest.mark.parametrize("max_points", [0, 1500, 10000])
@pytest.mark.parametrize("centered", [False, True])
def test_native_catalog_matches_full_array_reference(
    tmp_path, nested, masked, patch_enabled, max_points, centered,
):
    nside = 16
    pixels = np.arange(hp.nside2npix(nside))
    values = np.sin(pixels * 0.13) + pixels * 0.0001
    values[[3, 1200]] = hp.UNSEEN
    values[2000] = np.nan
    values[2300] = np.inf
    mask = (pixels % 3 != 0).astype(float)
    mask[100] = np.nan
    mask[1500] = hp.UNSEEN
    path = tmp_path / "map.fits"
    mask_path = tmp_path / "mask.fits" if masked else None
    hp.write_map(path, values, nest=nested, dtype=np.float64)
    if masked:
        hp.write_map(mask_path, mask, nest=nested, dtype=np.float64)
    patch = driver.AngularPatch(25, 290, 35, 130) if patch_enabled else None
    valid = np.isfinite(values) & ~hp.mask_bad(values)
    if masked:
        valid &= np.isfinite(mask) & ~hp.mask_bad(mask) & (mask > 0)
    before_patch = int(valid.sum())
    if patch is not None:
        theta, phi = hp.pix2ang(nside, pixels, nest=nested)
        valid &= patch.select(theta, phi)
    expected_pixels = pixels[valid]
    before_thinning = expected_pixels.size
    if max_points and expected_pixels.size > max_points:
        order = np.argsort(driver._splitmix64(expected_pixels, 123))[:max_points]
        expected_pixels = np.sort(expected_pixels[order])
    expected_positions = np.column_stack(hp.pix2vec(nside, expected_pixels, nest=nested))
    expected_kappa = values[expected_pixels]
    if centered:
        expected_kappa -= expected_kappa.mean(dtype=np.float64)
    for chunk in (1024, 4096):
        catalog = driver.catalog_from_healpix(
            path, mask_path=mask_path, patch=patch, max_points=max_points,
            sampling_seed=123, chunk_pixels=chunk, center_field=centered,
        )
        np.testing.assert_array_equal(catalog.positions, expected_positions)
        np.testing.assert_array_equal(catalog.kappa, expected_kappa)
        assert catalog.positions.dtype == np.float64
        assert catalog.positions.flags.c_contiguous
        assert catalog.kappa.flags.owndata
        assert catalog.metadata["eligible_pixels_before_patch"] == before_patch
        assert catalog.metadata["eligible_pixels_before_thinning"] == before_thinning
        assert catalog.metadata["loader"] == "memory-mapped-native-resolution"
        assert catalog.metadata["mask_preselected"] is masked
        if masked:
            np.testing.assert_array_equal(catalog.mask, np.ones(len(expected_pixels), dtype=np.uint8))
        else:
            assert catalog.mask is None
    # In-place centering must only touch the owned catalog, never the source.
    np.testing.assert_array_equal(hp.read_map(path, nest=nested, dtype=np.float64), values)


@pytest.mark.parametrize("mode", ["downgrade", "mask-ordering", "mask-resolution", "explicit"])
@pytest.mark.parametrize("masked", [False, True])
@pytest.mark.parametrize("centered", [False, True])
@pytest.mark.parametrize("max_points", [0, 1500])
def test_fallback_catalog_matches_full_array_reference(
    tmp_path, mode, masked, centered, max_points,
):
    if mode.startswith("mask-") and not masked:
        pytest.skip("mask fallback requires a mask")
    nside = 32
    pixels = np.arange(hp.nside2npix(nside))
    values = np.sin(pixels * 0.13).astype(np.float32)
    values[::17] = hp.UNSEEN
    path = tmp_path / "map.fits"
    hp.write_map(path, values, partial=mode == "explicit", dtype=np.float32)
    mask_path = tmp_path / "mask.fits" if masked else None
    if masked:
        mask_nside = 16 if mode == "mask-resolution" else nside
        mask = (np.arange(hp.nside2npix(mask_nside)) % 3 != 0).astype(float)
        nested_mask = mode == "mask-ordering"
        hp.write_map(mask_path, hp.reorder(mask, r2n=True) if nested_mask else mask,
                     nest=nested_mask, dtype=np.float64)
    target = 16 if mode == "downgrade" else nside
    expected_values = hp.read_map(path, dtype=np.float64)
    if target != nside:
        expected_values = hp.ud_grade(expected_values, target, power=0.0)
    valid = np.isfinite(expected_values) & ~hp.mask_bad(expected_values)
    if masked:
        expected_mask = hp.read_map(mask_path, dtype=np.float64)
        if hp.get_nside(expected_mask) != target:
            expected_mask = hp.ud_grade(expected_mask, target, power=0.0)
        valid &= np.isfinite(expected_mask) & ~hp.mask_bad(expected_mask) & (expected_mask > 0.5)
    patch = driver.AngularPatch(15, 325, 20, 160)
    theta, phi = hp.pix2ang(target, np.arange(expected_values.size))
    valid &= patch.select(theta, phi)
    expected_pixels = np.flatnonzero(valid)
    if max_points and expected_pixels.size > max_points:
        order = np.argsort(driver._splitmix64(expected_pixels, 123))[:max_points]
        expected_pixels = np.sort(expected_pixels[order])
    expected_kappa = expected_values[expected_pixels]
    if centered:
        expected_kappa -= expected_kappa.mean(dtype=np.float64)
    catalog = driver.catalog_from_healpix(
        path, nside_down=target, mask_path=mask_path, mask_threshold=0.5,
        center_field=centered, patch=patch, max_points=max_points,
        sampling_seed=123, chunk_pixels=1024,
    )
    np.testing.assert_array_equal(catalog.kappa, expected_kappa)
    np.testing.assert_array_equal(catalog.positions, np.column_stack(hp.pix2vec(target, expected_pixels)))
    assert catalog.metadata["loader"].startswith("healpy-")
    assert catalog.metadata["mask_preselected"] is masked
    if masked:
        np.testing.assert_array_equal(catalog.mask, np.ones(len(expected_pixels), dtype=np.uint8))
    else:
        assert catalog.mask is None


@pytest.mark.parametrize("nside_down", [0, 16])
def test_coordinate_scratch_is_bounded_and_released_between_chunks(tmp_path, monkeypatch, nside_down):
    nside = 32
    path = tmp_path / "map.fits"
    hp.write_map(path, np.ones(hp.nside2npix(nside)), dtype=np.float64)
    original = hp.pix2vec
    previous = []
    calls = []

    def checked_pix2vec(nside, pixels, **kwargs):
        assert len(pixels) <= 1024
        assert all(ref() is None for ref in previous)
        result = original(nside, pixels, **kwargs)
        previous[:] = [weakref.ref(array) for array in result]
        calls.append(len(pixels))
        return result

    monkeypatch.setattr(hp, "pix2vec", checked_pix2vec)
    catalog = driver.catalog_from_healpix(path, nside_down=nside_down, chunk_pixels=1024)
    assert len(calls) > 1
    assert sum(calls) == catalog.nbody
    assert all(ref() is None for ref in previous)


@pytest.mark.parametrize("bad_selection", ["invalid", "empty-mask", "empty-patch"])
def test_empty_selection_still_reports_clear_error(tmp_path, bad_selection):
    path = tmp_path / "map.fits"
    values = np.ones(hp.nside2npix(16))
    if bad_selection == "invalid":
        values[:] = hp.UNSEEN
    hp.write_map(path, values, dtype=np.float64)
    kwargs = {}
    if bad_selection == "empty-mask":
        mask_path = tmp_path / "mask.fits"
        hp.write_map(mask_path, np.zeros_like(values), dtype=np.float64)
        kwargs["mask_path"] = mask_path
    if bad_selection == "empty-patch":
        kwargs["patch"] = driver.AngularPatch(1, 1.001, 1, 1.001)
    with pytest.raises(ValueError, match="fewer than three"):
        driver.catalog_from_healpix(path, chunk_pixels=1024, **kwargs)
