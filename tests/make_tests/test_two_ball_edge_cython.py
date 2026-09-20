#!/usr/bin/env python3
"""Cython lifetime, sparse-window, and driver edge-correction regressions."""

from pathlib import Path
import json
import sys
import tempfile

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tests" / "python"))
from cyballs import cballs, search_method_id
from kappa_corr_all_engines import KappaCatalog, RunConfig, engine_parameters, run_engine_suite
from test_two_ball_edge_corrections import ENGINES, MMAX, BINS, RMIN, RMAX, catalog, brute_force, edge_solution


@pytest.mark.parametrize("engine", ENGINES)
def test_edge_getters_and_recoverable_failures(engine):
    if search_method_id(engine) < 0:
        pytest.skip(f"{engine} is not built")
    data = catalog()
    probe = cballs()
    try:
        probe.set_catalog(data[0], kappa=data[1], weights=data[2])
    except Exception:
        data = catalog(dimension=2)
    finally:
        probe.clear_catalogs()
    signal, window = brute_force(data)
    expected = edge_solution(signal, window)
    with tempfile.TemporaryDirectory(prefix="two-ball-edge-cython-") as tmp:
        balls = cballs()
        balls.set_catalog(data[0], kappa=data[1], weights=data[2])
        params = dict(searchMethod=engine, rangeN=RMAX, rminHist=RMIN,
                      sizeHistN=BINS, mChebyshev=MMAX, sizeHistPhi=8,
                      lengthBox=2, numberThreads=2, useLogHist=False,
                      usePeriodic=False, verbose=0, verbose_log=0,
                      rootDir=tmp, nsmooth=2,
                      options="KKKCorrelation,weights-norm,only-3pcf,"
                              "edge-corrections,no-normalize-HistZeta,no-out-Hist,"
                              + ("no-one-ball,no-smooth-pivot"
                                 if engine in {"kdtree-omp", "kdtree-mpi"}
                                 else ("no-two-balls,no-smooth-pivot"
                                       if engine.startswith(("kdtree-2balls-",
                                                             "balltree-2balls-"))
                                       else "no-two-balls")))
        try:
            for _ in range(2):
                balls.set(params)
                balls.Run(level=["MainLoop"])
                assert balls.getEdgeCorrectionCPUTime() >= 0.0
                assert balls.getEdgeCorrectionWallTime() > 0.0
                actual = np.array([balls.getHistZetaM_EE_complex(m + 1)
                                   for m in range(MMAX + 1)])
                np.testing.assert_allclose(actual, expected, rtol=2e-10, atol=2e-10)
                assert not list(Path(tmp).glob("histZetaM*"))
                balls.struct_cleanup()
                with pytest.raises(Exception):
                    balls.getHistZetaM_EE_complex(1)
                bad = dict(params, options=params["options"] + ",only-2pcf")
                balls.set(bad)
                with pytest.raises(Exception, match="edge-corrections|mutually exclusive"):
                    balls.Run(level=["MainLoop"])
                balls.struct_cleanup()
            balls.set(dict(params, options=params["options"].replace("edge-corrections,", "")))
            balls.Run(level=["MainLoop"])
            assert balls.getEdgeCorrectionCPUTime() == 0.0
            assert balls.getEdgeCorrectionWallTime() == 0.0
            with pytest.raises(Exception, match="enable edge-corrections"):
                balls.getHistZetaM_EE_complex(1)
        finally:
            balls.struct_cleanup()
            balls.clear_catalogs()


def test_sparse_mask_windows_and_driver():
    engine = "octree-2balls-omp"
    if search_method_id(engine) < 0:
        pytest.skip(f"{engine} is not built")
    data = catalog()
    with tempfile.TemporaryDirectory(prefix="two-ball-edge-driver-") as tmp:
        config = RunConfig(engines=(engine,), output_dir=Path(tmp),
                           theta_min=RMIN, theta_max=RMAX, theta_scale="au",
                           bins=BINS, multipoles=MMAX, threads=2, use_log_bins=False,
                           result_type="edge_effects", verbose=0, verbose_log=0,
                           options=("weights-norm", "no-two-balls", "no-out-Hist"))
        params = engine_parameters(config, engine, True)
        assert {"read-mask", "edge-corrections", "no-normalize-HistZeta"}.issubset(
            params["options"].split(","))
        for valid in (1, 2, 3):
            mask = np.zeros(len(data[0]), dtype=np.uint8)
            mask[:valid] = 1
            result = run_engine_suite(KappaCatalog(positions=data[0], kappa=data[1],
                                                  weights=data[2], mask=mask), config)[engine]
            assert not result["warnings"], result["warnings"]
            for m in range(1, MMAX + 2):
                assert np.all(np.isnan(result[f"zeta_edge_complex_{m}"]))
            assert not np.any(result["scalar_window_valid"])
            assert np.all(np.isin(result["scalar_window_status"], (2, 3)))
        # A nonconstant field must retain both components in saved Python output.
        mask = np.ones(len(data[0]), dtype=np.uint8)
        mask[8::3] = 0
        result = run_engine_suite(KappaCatalog(positions=data[0], kappa=data[1],
                                              weights=data[2], mask=mask), config)[engine]
        signal, window = brute_force((*data[:3], mask))
        expected = edge_solution(signal, window)
        assert np.any(abs(expected.imag) > 1e-6)
        for m in range(MMAX + 1):
            np.testing.assert_allclose(result[f"zeta_edge_complex_{m+1}"], expected[m],
                                       rtol=2e-10, atol=2e-10)
        with np.load(Path(tmp) / engine / "python/histograms.npz") as saved:
            np.testing.assert_array_equal(saved["zeta_edge_complex_2"],
                                          result["zeta_edge_complex_2"])
            for name in ("status", "valid", "window_monopole", "pivot_ratio"):
                key = f"scalar_window_{name}"
                np.testing.assert_array_equal(saved[key], result[key])
        metadata = json.loads((Path(tmp) / engine / "python/result.json").read_text())
        assert metadata["run_settings"]["effective"]["sizeHistN"] == BINS
        assert metadata["run_settings"]["effective"]["catalog_sizes"] == [64]
        assert metadata["run_settings"]["requested"]["searchMethod"] == engine


def test_octree_two_ball_ggg_compatibility_with_mask_and_edge():
    engines = ("octree-ggg-omp", "octree-2balls-omp")
    if any(search_method_id(engine) < 0 for engine in engines):
        pytest.skip("octree GGG compatibility engines are not both registered")
    data = list(catalog())
    data[3][8::3] = 0
    options = (
        "KKKCorrelation,weights-norm,read-mask,edge-corrections,"
        "no-normalize-HistZeta,no-one-ball,no-smooth-pivot,no-out-Hist"
    )

    def run(engine, root):
        balls = cballs()
        balls.set_catalog(data[0], kappa=data[1], weights=data[2], mask=data[3])
        if engine == "octree-2balls-omp":
            engine_options = options + ",legacy-one-ball"
        else:
            engine_options = options
        balls.set(dict(searchMethod=engine, rangeN=RMAX, rminHist=RMIN,
                       sizeHistN=BINS, mChebyshev=MMAX, sizeHistPhi=8,
                       lengthBox=2, numberThreads=2, useLogHist=False,
                       usePeriodic=False, verbose=0, verbose_log=0,
                       rootDir=str(root), nsmooth=2, options=engine_options))
        try:
            balls.Run(level=["MainLoop"])
            result = {
                "nn": balls.getHistNN().copy(),
                "xi": balls.getHistXi2pcf().copy(),
            }
            for order in range(1, MMAX + 2):
                result[f"edge-{order}"] = (
                    balls.getHistZetaM_EE_complex(order).copy())
                for component in range(1, 5):
                    result[f"zeta-{order}-{component}"] = (
                        balls.getHistZetaMsincos(order, component).copy())
            return result
        finally:
            balls.struct_cleanup()
            balls.clear_catalogs()

    with tempfile.TemporaryDirectory(prefix="octree-ggg-compat-cython-") as tmp:
        reference = run(engines[0], Path(tmp) / "ggg")
        compatible = run(engines[1], Path(tmp) / "compat")
    assert reference.keys() == compatible.keys()
    for name in reference:
        np.testing.assert_array_equal(reference[name], compatible[name],
                                      err_msg=name)


def test_driver_rejects_unsupported_edge_requests():
    config = RunConfig(engines=("balltree-omp",), output_dir=Path("unused"),
                       result_type="edge_effects")
    with pytest.raises(ValueError, match="does not support"):
        engine_parameters(config, "balltree-omp", False)
    config.options = ("only-2pcf",)
    with pytest.raises(ValueError, match="require 3PCF"):
        engine_parameters(config, "octree-2balls-mpi", True)
