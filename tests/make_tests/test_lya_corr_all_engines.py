"""DESI conversion, forest memory API, and multi-engine driver regressions."""
import json
import os
import shlex
import signal
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT), str(ROOT / "python")]
import lya_corr_all_engines as driver


@pytest.fixture
def desi_file(tmp_path):
    fits = pytest.importorskip("astropy.io.fits")
    ids = np.array([2**55+1, 2**55+2, 2**55+3], dtype=np.int64)
    meta = fits.BinTableHDU.from_columns([
        fits.Column(name="LOS_ID", format="K", array=ids),
        fits.Column(name="RA", format="D", unit="rad", array=[0., .002, .003]),
        fits.Column(name="DEC", format="D", unit="rad", array=[0., 0., .001]),
    ], name="METADATA")
    meta.header["BLINDING"] = "test-blind"
    wave = fits.ImageHDU(np.array([3700., 3800., 3900., 4000.]), name="LAMBDA")
    wave.header["BUNIT"] = "Angstrom"
    delta = np.array([[.2, .4, np.nan, .5], [-.1, .7, .5, .6], [.3, -.2, .8, .9]])
    weights = np.ones_like(delta)
    weights[0, 1] = 0
    weights[1, 2] = -1
    path = tmp_path / "delta.fits.gz"
    fits.HDUList([fits.PrimaryHDU(), wave, meta,
                  fits.ImageHDU(delta, name="DELTA_BLIND"),
                  fits.ImageHDU(weights, name="WEIGHT")]).writeto(path)
    return path, ids


def test_desi_layout_ids_weights_and_distances(desi_file):
    from astropy.cosmology import FlatLambdaCDM
    path, ids = desi_file
    cat = driver.read_desi([str(path)])
    assert cat.nbody == 9
    np.testing.assert_array_equal(np.unique(cat.forest_ids), ids)
    assert cat.forest_ids.dtype == np.int64
    assert cat.metadata["files"][0]["blinding"] == "test-blind"
    assert cat.delta[0] == .2
    chi = FlatLambdaCDM(H0=67.4, Om0=.315, Tcmb0=0).comoving_distance(
        3700/1215.67-1).value * .674
    np.testing.assert_allclose(np.linalg.norm(cat.positions[0]), chi, rtol=2e-14)
    np.testing.assert_array_equal(cat.weights, 1)


def test_desi_subsample_cache_and_duplicate(desi_file, tmp_path):
    path, ids = desi_file
    cat = driver.read_desi([str(path)], max_forests=2, pixel_stride=2)
    assert cat.nbody == 3
    cached = tmp_path / "cache.npz"
    driver.save_catalog(cached, cat)
    loaded = driver.read_npz(cached)
    np.testing.assert_array_equal(cat.forest_ids, loaded.forest_ids)
    np.testing.assert_array_equal(cat.positions, loaded.positions)
    with pytest.raises(ValueError, match="duplicate input"):
        driver.read_desi([str(path), str(path)])
    copied = tmp_path / "copied.fits.gz"
    copied.write_bytes(path.read_bytes())
    with pytest.raises(ValueError, match="duplicate LOS_ID"):
        driver.read_desi([str(path), str(copied)])


def test_desi_selection_and_column_failures(desi_file):
    path, _ = desi_file
    with pytest.raises(ValueError, match="no valid"):
        driver.read_desi([str(path)], z_min=4, z_max=5)
    with pytest.raises(ValueError, match="missing DELTA"):
        driver.read_desi([str(path)], delta_field="DELTA")
    with pytest.raises(ValueError, match="positive"):
        driver.read_desi([str(path)], pixel_stride=0)


def test_ascii_ids_are_not_floats(tmp_path):
    path = tmp_path / "pixels.txt"
    ids = [2**55+1, 2**55+2, 2**55+3]
    path.write_text("".join(f"{i+10} 0 0 0.25 1 {fid}\n" for i, fid in enumerate(ids)))
    cat = driver.read_ascii(path)
    np.testing.assert_array_equal(cat.forest_ids, ids)


@pytest.mark.parametrize("mutation", ["float_ids", "negative_weight", "zero_position", "nan_delta"])
def test_catalog_validation(mutation):
    cat = driver.synthetic_catalog(3, 2)
    if mutation == "float_ids":
        cat.forest_ids = cat.forest_ids.astype(float)
    if mutation == "negative_weight":
        cat.weights[0] = -1
    if mutation == "zero_position":
        cat.positions[0] = 0
    if mutation == "nan_delta":
        cat.delta[0] = np.nan
    with pytest.raises(ValueError):
        cat.normalized()


def test_engine_contracts():
    available = list(driver.LYA_ENGINES)
    assert len(available) == 14
    assert len(driver.resolve_engines(["all"], available, "both")) == 14
    assert len(driver.resolve_engines(["all-omp"], available, "2pcf")) == 3
    assert len(driver.resolve_engines(["all-mpi"], available, "3pcf")) == 2
    with pytest.raises(ValueError, match="same-quasar"):
        driver.resolve_engines(["octree-3pcf-3d-omp"], available)
    with pytest.raises(ValueError, match="--statistics"):
        driver.resolve_engines(["lya-2pcf-3pcf-omp"], available)
    with pytest.raises(ValueError, match="not compiled"):
        driver.resolve_engines(["lya-2pcf-omp"], [])


def test_memory_budget_guard(tmp_path):
    with pytest.raises(ValueError, match="histogram estimate"):
        driver.RunConfig(("lya-3pcf-omp",), tmp_path, r3_bins=100,
                         theta_bins=100, mu_bins=100).validate()


def test_relative_zero_and_sparse_keys():
    a = np.array([[0, .5, 0., 0., 2.], [1, 1.5, .2, 2., 10.]])
    b = np.array([[1, 1.5, .3, 3., 10.]])
    metrics, rows = driver.compare_products(a, b, "1d_2pcf")
    assert metrics["undefined_relative_bins"] == 1
    assert np.isnan(rows[0, -1])
    np.testing.assert_allclose(metrics["max_abs_relative"], .5)


def test_projection_sums_raw_weights(tmp_path):
    table = np.array([[0, 0, .5, .5, 1, 2, 2],
                      [0, 0, .5, .5, 3, 9, 3]])
    cfg = driver.RunConfig(("lya-1d-3pcf-omp",), tmp_path, r3_bins=1)
    projected = driver.projected_plot_data(table, "1d_3pcf", cfg)
    np.testing.assert_allclose(projected[0, 0], 11/5)
    assert np.isnan(projected[1, 1])


def _cyballs():
    cyballs = pytest.importorskip("cyballs")
    if not hasattr(cyballs.cballs, "set_forest_catalog"):
        pytest.fail("cyballs must be rebuilt with set_forest_catalog")
    return cyballs


def test_memory_api_rejects_bad_ids_and_recovers():
    cyballs = _cyballs()
    cat = driver.synthetic_catalog(3, 2)
    balls = cyballs.cballs()
    try:
        with pytest.raises(Exception, match="integer"):
            balls.set_forest_catalog(cat.positions, cat.delta, cat.weights,
                                     cat.forest_ids.astype(float))
        assert balls.set_forest_catalog(cat.positions, cat.delta, cat.weights, cat.forest_ids)
        assert balls.catalog_count == 1
    finally:
        balls.clear_catalogs()


@pytest.mark.parametrize("radial", [False, True])
def test_production_memory_oracles(tmp_path, radial):
    cyballs = _cyballs()
    import test_lya_forest_omp as three
    import test_lya_forest_1d_omp as one
    engines = tuple(n for n, s in driver.LYA_ENGINES.items()
                    if not s.mpi and s.radial == radial and cyballs.search_method_id(n) >= 0)
    if not engines:
        pytest.skip("OpenMP forest sibling not compiled")
    points = one.WIDE_ANGLE if radial else three.POINTS
    cat = driver.ForestCatalog(points[:, :3], points[:, 3], points[:, 4],
                              points[:, 5].astype(np.int64) + 2**55 + 1)
    cfg = driver.RunConfig(engines, tmp_path / "results", threads=2, plots=False,
                           rp_max=one.RP_MAX if radial else three.RP_MAX,
                           rp_bins=one.RP_BINS if radial else three.RP_BINS,
                           rt_max=three.RT_MAX, rt_bins=three.RT_BINS,
                           r3_max=one.R3_MAX if radial else three.R3_MAX,
                           r3_bins=one.R3_BINS_PER_SIDE if radial else three.R3_BINS,
                           theta_bins=three.THETA_BINS, mu_bins=three.MU_BINS)
    result = driver.run_engine_suite(cat, cfg)
    for name in result:
        root = cfg.output_dir / name
        orders = driver.LYA_ENGINES[name].orders
        if radial:
            if 2 in orders:
                one.assert_histogram_close(one.read_2pcf(root/"histXi2pcf_lya1d.txt"),
                                           one.oracle_2pcf(), name)
            if 3 in orders:
                one.assert_histogram_close(one.read_3pcf(root/"histZetaM_lya1d.txt"),
                                           one.oracle_3pcf(), name)
        else:
            if 2 in orders:
                three.assert_histogram_close(three.read_2pcf(root/"histXi2pcf_lya.txt"),
                                             three.oracle_2pcf(), name)
            if 3 in orders:
                three.assert_histogram_close(three.read_3pcf(root/"histZetaM_lya5d.txt"),
                                             three.oracle_3pcf(), name)
    summary = json.loads((cfg.output_dir / "summary.json").read_text())
    assert summary["catalog_registrations_per_rank"] == 1
    assert all(m["max_abs_correlation"] < 1e-12 for m in summary["comparisons"].values())
    with pytest.raises(RuntimeError, match="results already exist"):
        driver.run_engine_suite(cat, cfg)


def test_cli_help():
    p = subprocess.run([sys.executable, str(ROOT/"python/lya_corr_all_engines.py"), "--help"],
                       capture_output=True, text=True, timeout=30)
    assert p.returncode == 0, p.stderr
    assert "--mpi-ranks" in p.stdout and "--fits" in p.stdout


@pytest.mark.skipif(not os.environ.get("LYA_DRIVER_MPI_COMMAND"),
                    reason="set LYA_DRIVER_MPI_COMMAND to a two-rank launcher")
def test_mpi_multi_engine_driver(tmp_path):
    output = tmp_path / "mpi-results"
    command = shlex.split(os.environ["LYA_DRIVER_MPI_COMMAND"]) + [
        sys.executable, str(ROOT/"python/lya_corr_all_engines.py"),
        "--engine", "all", "--statistics", "both", "--threads", "2",
        "--synthetic-forests", "8", "--synthetic-pixels", "36", "--rp-max", "30",
        "--rt-max", "80", "--no-plots", "--output", str(output)]
    log_path = tmp_path / "mpi.log"
    with log_path.open("w") as log:
        proc = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT,
                                start_new_session=True)
        try:
            status = proc.wait(timeout=120)
        except subprocess.TimeoutExpired:
            os.killpg(proc.pid, signal.SIGKILL)
            proc.wait()
            pytest.fail("MPI driver timeout: " + log_path.read_text()[-4000:])
    assert status == 0, log_path.read_text()[-4000:]
    summary = json.loads((output/"summary.json").read_text())
    assert sum(name.endswith("-mpi") for name in summary["engines"]) == 7
    assert summary["catalog"]["pixels"] == 288
    assert summary["comparisons"]
    assert all(m["max_abs_correlation"] < 2e-12
               for m in summary["comparisons"].values())
