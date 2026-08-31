#!/usr/bin/env python3

import json
import importlib.util
import os
from pathlib import Path
import sys
import subprocess
import tempfile
from unittest.mock import patch

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))

from cyballs import cballs, search_method_id
from kappa_corr_all_engines import (
    KAPPA_ENGINES, KappaCatalog, RunConfig, engine_parameters, run_engine_suite,
    copy_engine_results, resolve_engines,
    INCOMPATIBLE_ENGINE_REASONS, METHOD_ALIASES,
)


EDGE_ENGINES = (
    "balltree-2balls-omp", "balltree-2balls-mpi",
    "octree-2balls-omp", "octree-2balls-mpi",
    "balltree-2balls-omp_3pcf", "balltree-2balls-mpi_3pcf",
    "octree-ggg-omp", "octree-ggg-mpi",
    "octree-balls4-omp", "octree-balls4-mpi",
)


def runnable_method(name):
    return (search_method_id(name) >= 0 and
            ("-mpi" not in name or importlib.util.find_spec("mpi4py") is not None))


def test_forest_engines_point_to_the_forest_driver():
    reasons = {name: reason for name, reason in INCOMPATIBLE_ENGINE_REASONS.items()
               if name.startswith("lya-")}
    assert len(reasons) == 14
    assert all("lya_corr_all_engines.py" in reason for reason in reasons.values())
    assert all("set_forest_catalog" in reason for reason in reasons.values())
    assert not set(reasons).intersection(KAPPA_ENGINES)


def test_physical_3d_aliases_do_not_become_kappa_engines():
    for alias, canonical in METHOD_ALIASES.items():
        assert alias in INCOMPATIBLE_ENGINE_REASONS
        assert canonical in INCOMPATIBLE_ENGINE_REASONS
        assert alias not in KAPPA_ENGINES


def test_renamed_balls4_driver_rejects_smooth_pivot():
    assert "octree-balls4-omp" in KAPPA_ENGINES
    assert "octree-kkk-balls4-omp" not in KAPPA_ENGINES
    config = RunConfig(engines=("octree-balls4-omp",), output_dir=Path("unused"))
    assert engine_parameters(config, "octree-balls4-omp", False)["useLogHist"]
    config.options = ("smooth-pivot",)
    with pytest.raises(ValueError, match="does not support options=smooth-pivot"):
        engine_parameters(config, "octree-balls4-omp", False)


def unit_sphere_catalog(count=32):
    indices = np.arange(count, dtype=np.float64)
    z = 1.0 - 2.0 * (indices + 0.5) / count
    phi = np.pi * (3.0 - np.sqrt(5.0)) * indices
    radius = np.sqrt(1.0 - z * z)
    positions = np.column_stack((radius * np.cos(phi), radius * np.sin(phi), z))
    kappa = positions[:, 0] - 0.4 * positions[:, 1] + 0.2 * positions[:, 2]
    kappa -= np.mean(kappa)
    return KappaCatalog(positions=positions, kappa=kappa)


def test_one_catalog_is_reused_across_enabled_engines():
    candidates = ("octree-ggg-omp", "kdtree-omp")
    engines = tuple(name for name in candidates if search_method_id(name) >= 0)
    if not engines:
        raise AssertionError("the build has no engine suitable for the kappa smoke test")

    with tempfile.TemporaryDirectory(prefix="ctreeballs-kappa-suite-") as tmp:
        output = Path(tmp)
        config = RunConfig(
            engines=engines,
            output_dir=output,
            theta_min=5.0,
            theta_max=100.0,
            bins=4,
            multipoles=2,
            threads=2,
            options=("only-2pcf", "no-out-Hist"),
            verbose=0,
            verbose_log=0,
        )
        results = run_engine_suite(unit_sphere_catalog(), config)

        if tuple(results) != engines:
            raise AssertionError(f"expected {engines}, got {tuple(results)}")
        for engine, result in results.items():
            if result["nbody"] != 32:
                raise AssertionError(f"{engine} published the wrong body count")
            if not {"r", "xi", "nn"}.issubset(result):
                raise AssertionError(f"{engine} did not publish its 2PCF arrays")

        summary = json.loads((output / "summary.json").read_text(encoding="utf-8"))
        if summary["catalog_reads"] != 1:
            raise AssertionError("the suite did not preserve the one-read contract")
        if summary["set_catalog_calls_per_process"] != 1:
            raise AssertionError("the suite registered the catalog more than once")
        if summary["failures"]:
            raise AssertionError(f"suite failures: {summary['failures']}")


def test_octree_two_ball_mask_parameters():
    for engine in ("octree-2balls-omp", "octree-2balls-mpi"):
        assert KAPPA_ENGINES[engine].supports_mask
        assert KAPPA_ENGINES[engine].mpi == engine.endswith("-mpi")
        config = RunConfig(engines=(engine,), output_dir=Path("unused"))
        assert "read-mask" not in engine_parameters(config, engine, False)["options"]
        for options in ((), ("read-mask",), ("read-mask,no-two-balls",)):
            config.options = options
            params = engine_parameters(config, engine, True)
            assert params["options"].split(",").count("read-mask") == 1
            assert params["iCatalogs"] == "1"
            assert params["searchMethod"] == engine


def test_masked_octree_suite_matches_filtered_catalog():
    engines = tuple(name for name in ("octree-2balls-omp", "octree-2balls-mpi")
                    if runnable_method(name))
    if not engines:
        print("SKIP: masked suite requires an octree two-ball engine")
        return

    catalog = unit_sphere_catalog(64)
    mask = (np.arange(catalog.nbody) % 3 != 0).astype(np.uint8)
    # Preserve the startup bounding-box reference frame in the filtered input.
    mask[np.argmin(catalog.positions, axis=0)] = 1
    mask[np.argmax(catalog.positions, axis=0)] = 1
    selected = mask.astype(bool)
    catalog.mask = mask
    catalog.weights = np.linspace(0.3, 1.7, catalog.nbody)
    filtered = KappaCatalog(positions=catalog.positions[selected],
                            kappa=catalog.kappa[selected],
                            weights=catalog.weights[selected])
    catalog.kappa[~selected] = 1.e6
    catalog.weights[~selected] = 1.e5
    originals = [array.copy() for array in
                 (catalog.positions, catalog.kappa, catalog.weights, mask)]
    registrations = []

    class RecordingBalls:
        def __init__(self):
            self.inner = cballs()

        def __getattr__(self, name):
            return getattr(self.inner, name)

        def set_catalog(self, *args, **kwargs):
            registrations.append((args, kwargs))
            return self.inner.set_catalog(*args, **kwargs)

    with tempfile.TemporaryDirectory(prefix="ctreeballs-kappa-mask-") as tmp:
        config = RunConfig(
            engines=engines, output_dir=Path(tmp) / "masked",
            theta_min=5.0, theta_max=70.0, bins=4, multipoles=2,
            threads=2, nsmooth=4, verbose=0, verbose_log=0,
            options=("no-two-balls", "no-out-Hist", "kappa-avg"),
        )
        with patch("cyballs.cballs", RecordingBalls):
            results = run_engine_suite(catalog, config)
        assert len(registrations) == 1, "engines registered the catalog repeatedly"
        args, kwargs = registrations[0]
        assert args[0] is catalog.positions
        assert kwargs["mask"] is catalog.mask
        assert kwargs["kappa"] is catalog.kappa
        assert kwargs["weights"] is catalog.weights

        config.output_dir = Path(tmp) / "filtered"
        reference = run_engine_suite(filtered, config)
        assert tuple(results) == engines
        for engine in engines:
            result = results[engine]
            assert not result["warnings"], result["warnings"]
            arrays = {name: value for name, value in result.items()
                      if isinstance(value, np.ndarray)}
            assert {"r", "xi", "nn", "zeta_cos_1"}.issubset(arrays)
            assert np.any(arrays["nn"] > 0)
            assert np.any(arrays["zeta_cos_1"] != 0)
            for name, value in arrays.items():
                assert np.all(np.isfinite(value)), (engine, name)
                np.testing.assert_allclose(value, reference[engine][name],
                                           rtol=3.e-12, atol=3.e-13,
                                           err_msg=f"{engine}: {name}")
                np.testing.assert_allclose(value, results[engines[0]][name],
                                           rtol=3.e-12, atol=3.e-13,
                                           err_msg=f"OMP/MPI: {name}")

    for original, current in zip(originals,
            (catalog.positions, catalog.kappa, catalog.weights, catalog.mask)):
        np.testing.assert_array_equal(original, current)


def test_unsupported_mask_engine_still_rejected():
    catalog = unit_sphere_catalog()
    catalog.mask = np.ones(catalog.nbody, dtype=np.uint8)
    with tempfile.TemporaryDirectory(prefix="ctreeballs-kappa-mask-error-") as tmp:
        config = RunConfig(engines=("balltree-2balls-omp",), output_dir=Path(tmp))
        with patch("cyballs.cballs") as constructor:
            try:
                run_engine_suite(catalog, config)
            except ValueError as exc:
                assert "read-mask" in str(exc)
                assert "octree-2balls-omp" in str(exc)
                assert "octree-2balls-mpi" in str(exc)
            else:
                raise AssertionError("an unsupported masked engine was accepted")
            constructor.assert_not_called()


@pytest.mark.parametrize("edge_request", (
    {"edge_corrections": True}, {"result_type": "edge_effects"},
    {"options": ("edge-corrections,no-out-Hist",)},
    {"options": "edge-corrections,no-out-Hist"},
))
def test_all_edge_request_forms_select_corrected_results(edge_request):
    config = RunConfig(engines=EDGE_ENGINES, output_dir=Path("unused"), **edge_request)
    assert config.wants_edge_corrections
    normalized = config.normalized()
    assert normalized.edge_corrections
    assert normalized.result_type == "edge_effects"
    assert normalized.normalized() == normalized
    for engine in EDGE_ENGINES:
        options = engine_parameters(config, engine, False)["options"].split(",")
        assert options.count("edge-corrections") == 1
        assert options.count("no-normalize-HistZeta") == 1
    config.options = ("edge-corrections,only-2pcf",)
    with pytest.raises(ValueError, match="require 3PCF"):
        config.normalized()


def test_edge_engine_groups_filter_capabilities_but_explicit_names_fail():
    available = list(KAPPA_ENGINES)
    assert tuple(resolve_engines(("all",), available, edge_corrections=True)) == EDGE_ENGINES
    for group, mpi in (("all-omp", False), ("all-mpi", True)):
        assert resolve_engines((group,), available, edge_corrections=True) == [
            name for name in EDGE_ENGINES if KAPPA_ENGINES[name].mpi == mpi
        ]
    assert resolve_engines(("all",), available, edge_corrections=True, masked=True) == [
        name for name in EDGE_ENGINES if KAPPA_ENGINES[name].supports_mask
    ]
    assert resolve_engines(("all",), ["octree-2balls-omp"], edge_corrections=True) == [
        "octree-2balls-omp"
    ]
    with pytest.raises(ValueError, match="does not support.*edge"):
        resolve_engines(("kdtree-omp",), available, edge_corrections=True)
    with pytest.raises(ValueError, match="read-mask"):
        resolve_engines(("balltree-2balls-omp",), available, edge_corrections=True, masked=True)
    with pytest.raises(ValueError, match="no requested.*edge corrections"):
        resolve_engines(("all",), ["kdtree-omp"], edge_corrections=True)


class EdgeResultStub:
    def getCPUTime(self):
        return 0.0

    def getNBody(self):
        return 32

    def getrBins(self):
        return np.arange(4)

    def getnMultipoles(self):
        return 2

    def getHistZetaM_EE(self, order):
        return np.full((4, 4), float(order))

    def getHistZetaM_EE_Im(self, order):
        return np.full((4, 4), -float(order))


def test_more_options_copies_complex_correction_not_raw_components():
    config = RunConfig(engines=(EDGE_ENGINES[0],), output_dir=Path("unused"),
                       bins=4, options=("edge-corrections,only-3pcf",))
    result = copy_engine_results(EdgeResultStub(), EDGE_ENGINES[0], config)
    assert result["edge_corrections"] is True
    assert result["result_type"] == "edge_effects"
    assert not result["warnings"]
    for order in range(1, 4):
        np.testing.assert_array_equal(result[f"zeta_edge_complex_{order}"], order * (1 - 1j))
        assert f"zeta_cos_{order}" not in result


@pytest.mark.parametrize("failure", ("missing", "nonfinite", "shape"))
def test_missing_or_invalid_corrected_results_fail(failure):
    balls = EdgeResultStub()
    if failure == "missing":
        def missing(order):
            raise AttributeError("extension has no imaginary edge getter")
        balls.getHistZetaM_EE_Im = missing
    elif failure == "nonfinite":
        balls.getHistZetaM_EE_Im = lambda order: np.full((4, 4), np.nan)
    else:
        balls.getHistZetaM_EE_Im = lambda order: np.zeros((3, 3))
    config = RunConfig(engines=(EDGE_ENGINES[0],), output_dir=Path("unused"),
                       bins=4, options=("only-3pcf",), edge_corrections=True)
    with pytest.raises(RuntimeError, match="requested edge-corrected 3PCF is unavailable"):
        copy_engine_results(balls, EDGE_ENGINES[0], config)


def test_one_catalog_edge_suite_publishes_all_available_engines(tmp_path):
    engines = tuple(name for name in EDGE_ENGINES if runnable_method(name))
    if not engines:
        pytest.skip("no edge-capable engine is built")
    registrations = []

    class RecordingBalls:
        def __init__(self):
            self.inner = cballs()

        def __getattr__(self, name):
            return getattr(self.inner, name)

        def set_catalog(self, *args, **kwargs):
            registrations.append((args, kwargs))
            return self.inner.set_catalog(*args, **kwargs)

    catalog = unit_sphere_catalog(64)
    catalog.weights = np.linspace(0.3, 1.7, catalog.nbody)
    config = RunConfig(
        engines=engines, output_dir=tmp_path, theta_min=5, theta_max=100,
        bins=4, multipoles=2, threads=2, nsmooth=2, verbose=0, verbose_log=0,
        options=("weights-norm,no-two-balls,no-one-ball,no-out-Hist",),
        edge_corrections=True,
    )
    with patch("cyballs.cballs", RecordingBalls):
        results = run_engine_suite(catalog, config)
    assert len(registrations) == 1
    assert registrations[0][0][0] is catalog.positions
    assert tuple(results) == engines
    summary = json.loads((tmp_path / "summary.json").read_text())
    assert summary["set_catalog_calls_per_process"] == 1
    assert not summary["failures"]
    for engine, result in results.items():
        assert not result["warnings"], (engine, result["warnings"])
        assert result["edge_corrections"] is True
        assert summary["engines"][engine]["result_type"] == "edge_effects"
        with np.load(tmp_path / engine / "python/histograms.npz") as saved:
            for order in range(1, 4):
                value = result[f"zeta_edge_complex_{order}"]
                assert value.shape == (4, 4)
                assert np.all(np.isfinite(value))
                np.testing.assert_array_equal(saved[f"zeta_edge_complex_{order}"], value)
        assert np.any(result["zeta_edge_complex_1"] != 0), engine


@pytest.mark.parametrize("edge_args", (
    ("--edge-corrections",), ("--type", "edge_effects"),
    ("--more-options", "edge-corrections"),
))
def test_cli_edge_groups_with_npz_mask(tmp_path, edge_args):
    available = [name for name in EDGE_ENGINES if search_method_id(name) >= 0
                 and not KAPPA_ENGINES[name].mpi and KAPPA_ENGINES[name].supports_mask]
    if not available:
        pytest.skip("no masked edge-capable OpenMP engine is built")
    catalog = unit_sphere_catalog(64)
    catalog.mask = (np.arange(catalog.nbody) % 3 != 0).astype(np.uint8)
    source = tmp_path / "catalog.npz"
    np.savez(source, positions=catalog.positions, kappa=catalog.kappa, mask=catalog.mask)
    proc = subprocess.run([
        sys.executable, str(ROOT / "python/kappa_corr_all_engines.py"),
        "--catalog-npz", str(source), "--engine", "all-omp", *edge_args,
        "--outdir", str(tmp_path / "results"), "--threads", "2",
        "--theta-min", "5", "--theta-max", "100", "--nbins", "4",
        "--multipoles", "2", "--nsmooth", "2", "--verbose", "0", "--verbose-log", "0",
        "--more-options", "no-two-balls,no-one-ball,no-out-Hist",
    ], text=True, capture_output=True, timeout=120, env=dict(os.environ))
    if (sys.platform == "darwin" and proc.returncode == -6
            and "OMP: Error #179: Function Can't open SHM failed" in proc.stderr):
        pytest.skip("host OpenMP runtime cannot initialize shared memory in the subprocess")
    assert proc.returncode == 0, proc.stdout + proc.stderr
    summary = json.loads((tmp_path / "results/summary.json").read_text())
    assert tuple(summary["engines"]) == tuple(available)
    assert not summary["failures"]
    assert summary["catalog_reads"] == 1
    assert all(item["edge_corrections"] for item in summary["engines"].values())


if __name__ == "__main__":
    test_one_catalog_is_reused_across_enabled_engines()
    print("PASS: kappa engines reuse one in-memory cyballs catalog")
    test_octree_two_ball_mask_parameters()
    test_masked_octree_suite_matches_filtered_catalog()
    test_unsupported_mask_engine_still_rejected()
    print("PASS: masked octree two-ball routing and available-engine 2PCF/3PCF")
