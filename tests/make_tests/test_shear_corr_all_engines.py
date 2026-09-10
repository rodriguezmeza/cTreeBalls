#!/usr/bin/env python3
"""Contracts for the active full-sky shear all-engines driver."""

import os
from pathlib import Path
import sys

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, os.fspath(ROOT / "tests" / "python"))

from shear_corr_all_engines import (  # noqa: E402
    ENGINE_ORDER,
    RunConfig,
    SHEAR_SPHERE_BALLTREE_TWO_BALLS_ENGINE,
    SHEAR_SPHERE_KDTREE_TWO_BALLS_ENGINE,
    SHEAR_SPHERE_TWO_BALLS_ENGINE,
    catalog_from_healpix_sphere,
    compare_results,
    parse_arguments,
    resolve_engines,
    stereographic_shear_patch,
    synthetic_spherical_shear_catalog,
    validate_spherical_smooth_radius,
)


ACTIVE_ENGINES = (
    SHEAR_SPHERE_TWO_BALLS_ENGINE,
    SHEAR_SPHERE_KDTREE_TWO_BALLS_ENGINE,
    SHEAR_SPHERE_BALLTREE_TWO_BALLS_ENGINE,
)


def test_registry_is_limited_to_active_full_sky_addons():
    assert ENGINE_ORDER == ACTIVE_ENGINES
    assert resolve_engines(("all",), ACTIVE_ENGINES, "sphere") == list(ACTIVE_ENGINES)
    assert resolve_engines((ACTIVE_ENGINES[1],), ACTIVE_ENGINES, "sphere") == [
        ACTIVE_ENGINES[1]
    ]
    with pytest.raises(ValueError, match="sphere"):
        resolve_engines(("all",), ACTIVE_ENGINES, "flat")


def test_projection_at_patch_center():
    epsilon = 0.03
    vectors = np.array((
        (1.0, 0.0, 0.0),
        (np.cos(epsilon), np.sin(epsilon), 0.0),
        (np.cos(epsilon), -np.sin(epsilon), 0.0),
    ))
    positions, gamma1, gamma2 = stereographic_shear_patch(
        vectors, np.ones(3), np.zeros(3),
        center_ra_deg=0.0, center_dec_deg=0.0,
    )
    np.testing.assert_allclose(positions[0], 0.0, atol=1.0e-15)
    np.testing.assert_allclose(
        positions[1:, 0],
        (2.0*np.tan(epsilon/2.0), -2.0*np.tan(epsilon/2.0)),
        rtol=2.0e-14, atol=1.0e-15,
    )
    np.testing.assert_allclose(positions[:, 1:], 0.0, atol=1.0e-15)
    np.testing.assert_allclose(gamma1, 1.0, atol=2.0e-15)
    np.testing.assert_allclose(gamma2, 0.0, atol=2.0e-15)


def test_spherical_catalog_normalizes_positions():
    catalog = synthetic_spherical_shear_catalog(32).normalized()
    assert catalog.geometry == "sphere"
    np.testing.assert_allclose(np.linalg.norm(catalog.positions, axis=1), 1.0,
                               rtol=0.0, atol=3.0e-15)


def test_spherical_smooth_radius_guard():
    safe = RunConfig(
        min_sep=3.0, max_sep=120.0, sep_units="arcmin",
        smooth_radius=1.0, bins=4, multipoles=2, plots=False,
    ).normalized()
    validate_spherical_smooth_radius(safe, ())
    unsafe = RunConfig(
        min_sep=3.0, max_sep=120.0, sep_units="arcmin",
        smooth_radius=2.0, bins=4, multipoles=2, plots=False,
    ).normalized()
    with pytest.raises(ValueError, match=r"2\*rsmooth <= min-sep"):
        validate_spherical_smooth_radius(unsafe, ())
    validate_spherical_smooth_radius(unsafe, ("no-smooth-pivot",))


def test_arcmin_units_match_degrees():
    arcmin = RunConfig(min_sep=3, max_sep=120, sep_units="arcmin",
                       bins=4, multipoles=2).normalized()
    degree = RunConfig(min_sep=.05, max_sep=2, sep_units="degree",
                       bins=4, multipoles=2).normalized()
    np.testing.assert_allclose(arcmin.native_limits("sphere"),
                               degree.native_limits("sphere"), rtol=1e-15)


def test_full_sky_healpix_reader(tmp_path):
    hp = pytest.importorskip("healpy")
    nside = 4
    pixels = np.arange(hp.nside2npix(nside))
    x, y, z = hp.pix2vec(nside, pixels)
    path = tmp_path / "shear.fits"
    hp.write_map(path, (0.02*x + 0.01*z, 0.03*y), dtype=np.float64,
                 column_names=("GAMMA1", "GAMMA2"), overwrite=True)
    catalog = catalog_from_healpix_sphere(path, max_points=17)
    assert catalog.geometry == "sphere" and catalog.nbody == 17
    np.testing.assert_allclose(np.linalg.norm(catalog.positions, axis=1), 1.0,
                               rtol=0.0, atol=3.0e-15)


def test_comparison_uses_first_active_engine_as_reference():
    base = {
        "radius": np.array([0.1, 0.2]),
        "xi_plus": np.array([1+0j, 2+0j]),
        "xi_minus": np.array([.5+0j, 1+0j]),
        "pair_weight": np.array([2.0, 3.0]),
        "orders": np.array([-1, 0, 1]),
        "upsilon": np.ones((4, 2, 2, 3), dtype=complex),
        "gamma": np.ones((4, 2, 2, 3), dtype=complex),
        "window": np.ones((2, 2, 3), dtype=complex),
    }
    candidate = {key: value.copy() for key, value in base.items()}
    comparison = compare_results({ACTIVE_ENGINES[0]: base,
                                  ACTIVE_ENGINES[1]: candidate})
    assert comparison["reference"] == ACTIVE_ENGINES[0]
    assert comparison["engines"][ACTIVE_ENGINES[1]]["2pcf"]["xi_plus"][
        "max"
    ] == 0.0


def test_cli_defaults_to_spherical_active_engines():
    args = parse_arguments(["--synthetic-nbody", "32", "--list-engines"])
    assert args.geometry == "sphere"
    assert args.engines == []


if __name__ == "__main__":
    raise SystemExit(pytest.main([os.fspath(Path(__file__))]))
