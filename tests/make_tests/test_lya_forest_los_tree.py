#!/usr/bin/env python3
"""Exact 3D equivalence of octree discovery + per-forest radial trees."""
from __future__ import annotations

import os
from pathlib import Path
import re
import subprocess

import numpy as np
import pytest

import test_lya_forest_omp as reference

ROOT = Path(__file__).resolve().parents[2]
BINARY = Path(os.environ.get("CBALLS", ROOT / "cballs")).resolve()
PRODUCTS = {2: "histXi2pcf_lya.txt", 3: "histZetaM_lya5d.txt"}


def run(catalog, output, method, threads=1, rp=30, rt=30, radius=30, extra=()):
    output.mkdir()
    args = [
        str(BINARY), f"search={method}", f"infile={catalog}",
        "infileformat=lya-ascii", "iCatalogs=1", f"rootDir={output}",
        f"numberThreads={threads}", "usePeriodic=false", "useLogHist=false",
        "rangeN=30", "rminHist=0.1", "sizeHistN=4",
        f"lya2RpMax={rp}", f"lya2RtMax={rt}",
        f"lya2RpBins={reference.RP_BINS}", f"lya2RtBins={reference.RT_BINS}",
        f"lya3RMax={radius}", f"lya3RBins={reference.R3_BINS}",
        f"lya3ThetaBins={reference.THETA_BINS}", f"lya3MuBins={reference.MU_BINS}",
        "verbose=2", "verbose_log=2", "options=",
    ]
    parameters = dict(arg.split("=", 1) for arg in args[1:])
    parameters.update(arg.split("=", 1) for arg in extra)
    args = [args[0], *(f"{key}={value}" for key, value in parameters.items())]
    result = subprocess.run(args, text=True, capture_output=True, timeout=120)
    assert result.returncode == 0, result.stdout + result.stderr
    return result.stdout


def catalog_file(tmp_path, points):
    path = tmp_path / "forest.txt"
    np.savetxt(path, points, fmt=("%.17g",) * 5 + ("%.0f",))
    return path


def compare_products(left, right, orders=(2, 3)):
    for order in orders:
        name = PRODUCTS[order]
        reader = reference.read_2pcf if order == 2 else reference.read_3pcf
        reference.assert_histogram_close(reader(left / name), reader(right / name), name)
        # Counts include accepted zero-weight contributions, not just occupied bins.
        count = re.compile(r"# distinct-forest .*: (\d+)")
        assert count.search((left / name).read_text())[1] == count.search(
            (right / name).read_text())[1]


def make_forests(forests=7, pixels=25, seed=919):
    rng = np.random.default_rng(seed)
    directions = np.column_stack((np.ones(forests), rng.uniform(-0.13, 0.13, (forests, 2))))
    directions /= np.linalg.norm(directions, axis=1)[:, None]
    chi = rng.uniform(80, 260, (forests, pixels))
    xyz = (directions[:, None, :] * chi[:, :, None]).reshape(-1, 3)
    points = np.column_stack((xyz, rng.normal(size=forests * pixels),
                             rng.uniform(0.1, 2, forests * pixels),
                             np.repeat(np.arange(forests) * 37 - 2**40, pixels)))
    rng.shuffle(points)
    return points


@pytest.mark.parametrize("order", ["2pcf", "3pcf", "2pcf-3pcf"])
def test_oracle_modes_and_thread_determinism(tmp_path, order):
    catalog = catalog_file(tmp_path, reference.POINTS)
    old, one, many = (tmp_path / name for name in ("old", "one", "many"))
    run(catalog, old, f"lya-{order}-omp")
    run(catalog, one, f"lya-los-tree-{order}-omp")
    run(catalog, many, f"lya-los-tree-{order}-omp", threads=4)
    orders = (2, 3) if order == "2pcf-3pcf" else (int(order[0]),)
    compare_products(one, old, orders)
    for statistic in orders:
        reader = reference.read_2pcf if statistic == 2 else reference.read_3pcf
        oracle = reference.oracle_2pcf if statistic == 2 else reference.oracle_3pcf
        reference.assert_histogram_close(reader(one / PRODUCTS[statistic]), oracle(), order)
        assert (one / PRODUCTS[statistic]).read_bytes() == (many / PRODUCTS[statistic]).read_bytes()


@pytest.mark.parametrize("case", ["long", "skew", "zero-weights", "one-forest", "outside", "boundary"])
def test_adversarial_geometry(tmp_path, case):
    points = make_forests()
    if case == "skew":
        # The input contract permits a forest ID with non-collinear pixels.
        points[::4, 1] += 19
    elif case == "zero-weights":
        points[::3, 4] = 0
    elif case == "one-forest":
        points[:, 5] = 2**42
    elif case == "outside":
        points = np.array([(100, 0, 0, 1, 1, 0), (0, 100, 0, 2, 1, 1),
                           (-100, 0, 0, 3, 1, 2)], dtype=float)
    elif case == "boundary":
        # Tangency, both sphere boundaries, near-collinear LOS, and points in
        # the 2PCF rectangle whose distance exceeds either individual axis max.
        points = np.array([
            (100, 0, 0, .2, 1, 0), (80, 0, 0, .3, 1, 0),
            (100, 7, 0, .4, 1, 1), (100, 7 - 1e-10, .01, .5, 1, 2),
            (97, .000001, 0, -.1, .8, 3), (103, .000001, 0, .1, 1, 3),
            (102.7, 6, 0, .9, 1, 4), (97.3, -6, 0, .8, .6, 5),
            # A transverse offset keeps these distinct under float storage too.
            (107, 0, 0, -.4, 1, 6), (107 - 1e-10, 0, 1e-6, .6, .7, 7),
        ])
    catalog = catalog_file(tmp_path, points)
    old, new, many = (tmp_path / name for name in ("old", "new", "many"))
    run(catalog, old, "lya-2pcf-3pcf-omp", rp=3, rt=7, radius=7)
    stdout = run(catalog, new, "lya-los-tree-2pcf-3pcf-omp", rp=3, rt=7, radius=7)
    run(catalog, many, "lya-los-tree-2pcf-3pcf-omp", threads=4, rp=3, rt=7, radius=7)
    compare_products(new, old)
    for name in PRODUCTS.values():
        assert (new / name).read_bytes() == (many / name).read_bytes()
    assert "LOS-tree: forests=" in stdout
    if case in ("long", "one-forest"):
        assert int(re.search(r"forest_skips=(\d+)", stdout)[1]) > 0


def test_combined_matches_independent_domains(tmp_path):
    catalog = catalog_file(tmp_path, make_forests(forests=5, pixels=17))
    both = tmp_path / "both"
    run(catalog, both, "lya-los-tree-2pcf-3pcf-omp", rp=22, rt=19, radius=11)
    for order in (2, 3):
        separate = tmp_path / str(order)
        run(catalog, separate, f"lya-los-tree-{order}pcf-omp", rp=22, rt=19, radius=11)
        compare_products(both, separate, (order,))


def test_new_methods_are_advertised():
    result = subprocess.run([str(BINARY), "options=print-search-methods"],
                            text=True, capture_output=True, timeout=30)
    # Informational early exits use status 1 in the native command line.
    assert result.returncode in (0, 1)
    for order in ("2pcf", "3pcf", "2pcf-3pcf"):
        assert f"lya-los-tree-{order}-omp" in result.stdout


@pytest.mark.parametrize("extra,reason", [("usePeriodic=true", "periodic"),
                                          ("lya3ThetaBins=0", "invalid Lyman-alpha 3PCF"),
                                          ("lya2RtMax=0", "invalid Lyman-alpha 2PCF")])
def test_3d_validation_is_not_radial_only(tmp_path, extra, reason):
    catalog = catalog_file(tmp_path, reference.POINTS)
    with pytest.raises(AssertionError, match=reason):
        run(catalog, tmp_path / "bad", "lya-los-tree-2pcf-3pcf-omp", extra=(extra,))


def test_forest_ids_are_not_narrowed(tmp_path):
    ids = [2**62 + 1, 2**62 + 2, -(2**62) - 1, -(2**62) - 2]
    path = tmp_path / "large-ids.txt"
    path.write_text("".join(" ".join(format(x, ".17g") for x in row[:5])
                            + f" {ids[int(row[5])]}\n" for row in reference.POINTS))
    old, new = tmp_path / "old", tmp_path / "new"
    run(path, old, "lya-2pcf-3pcf-omp")
    run(path, new, "lya-los-tree-2pcf-3pcf-omp", threads=3)
    compare_products(new, old)
    reference.assert_histogram_close(reference.read_2pcf(new / PRODUCTS[2]),
                                      reference.oracle_2pcf(), "wide IDs")


@pytest.mark.parametrize("seed", range(5))
def test_random_windows_against_octree(tmp_path, seed):
    catalog = catalog_file(tmp_path, make_forests(forests=9, pixels=31, seed=seed))
    old, new = tmp_path / "old", tmp_path / "new"
    domain = dict(rp=3.7 + seed * 4.1, rt=6.3 + seed * 2.7, radius=9.1 + seed * 2.3)
    run(catalog, old, "lya-2pcf-3pcf-omp", **domain)
    run(catalog, new, "lya-los-tree-2pcf-3pcf-omp", threads=3, **domain)
    compare_products(new, old)
