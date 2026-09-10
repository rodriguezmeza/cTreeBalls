#!/usr/bin/env python3
"""Contracts for the active convergence all-engines driver."""

import json
import math
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tests" / "python"))

from kappa_corr_all_engines import (  # noqa: E402
    INCOMPATIBLE_ENGINE_REASONS,
    KAPPA_ENGINES,
    KappaCatalog,
    RunConfig,
    discover_make_settings,
    engine_parameters,
    flatten_radial_matrix,
    make_plots,
    merge_statistics_options,
    print_engine_table,
    resolve_engines,
    run_engine_suite,
    smooth_pivot_mode,
    solve_scalar_mode_coupling,
)


ACTIVE_ENGINES = (
    "kdtree-2balls-omp", "kdtree-2balls-mpi",
    "balltree-2balls-omp", "balltree-2balls-mpi",
    "octree-2balls-omp", "octree-2balls-mpi",
)


def unit_sphere_catalog(count=32):
    index = np.arange(count, dtype=float)
    z = 1.0 - 2.0*(index + 0.5)/count
    phi = np.pi*(3.0 - np.sqrt(5.0))*index
    radius = np.sqrt(1.0-z*z)
    positions = np.column_stack((radius*np.cos(phi), radius*np.sin(phi), z))
    kappa = positions[:, 0] - 0.4*positions[:, 1] + 0.2*positions[:, 2]
    return KappaCatalog(positions, kappa-np.mean(kappa))


def test_registry_is_limited_to_active_scalar_addons():
    assert tuple(KAPPA_ENGINES) == ACTIVE_ENGINES
    assert resolve_engines(("all-omp",), ACTIVE_ENGINES) == list(ACTIVE_ENGINES[::2])
    assert resolve_engines(("all-mpi",), ACTIVE_ENGINES) == list(ACTIVE_ENGINES[1::2])
    assert resolve_engines(("all",), ACTIVE_ENGINES) == list(ACTIVE_ENGINES)


def test_other_field_families_are_rejected():
    for name in ("lya-2pcf-omp", "octree-shear-sphere-2balls-omp",
                 "octree-3pcf-3d-omp", "neighbor-boxes-omp"):
        assert name in INCOMPATIBLE_ENGINE_REASONS
        with pytest.raises(ValueError):
            resolve_engines((name,), (name,))


def test_statistics_selector_maps_to_native_options():
    assert merge_statistics_options("both", ()) == ()
    assert merge_statistics_options("2pcf", ()) == ("only-2pcf",)
    assert merge_statistics_options("3pcf", ("no-out-Hist",)) == (
        "no-out-Hist", "only-3pcf",
    )
    with pytest.raises(ValueError, match="conflicts"):
        merge_statistics_options("2pcf", ("only-3pcf",))


def test_active_engine_capabilities_and_parameters():
    for name, spec in KAPPA_ENGINES.items():
        assert spec.supports_mask and spec.supports_edge
        assert spec.supports_dual_node_bin_slop
        assert spec.mpi == name.endswith("-mpi")
        config = RunConfig(
            engines=(name,), output_dir=Path("unused"),
            dual_node_bin_slop=True,
        )
        options = engine_parameters(config, name, True)["options"].split(",")
        assert options.count("read-mask") == 1
        assert options.count("dual-node-bin-slop") == 1


def test_smooth_pivot_policy():
    for name in ACTIVE_ENGINES[:4]:
        config = RunConfig((name,), Path("unused"), smooth_pivot_compiled=True)
        assert smooth_pivot_mode(config, name) == "enabled-by-build-default"
        config.options = ("no-smooth-pivot",)
        assert smooth_pivot_mode(config, name) == "disabled-by-option"
    for name in ACTIVE_ENGINES[4:]:
        config = RunConfig((name,), Path("unused"), smooth_pivot_compiled=True)
        assert smooth_pivot_mode(config, name) == "unsupported"


def test_edge_correction_contract():
    config = RunConfig(
        ACTIVE_ENGINES, Path("unused"), edge_corrections=True,
        options=("only-3pcf",),
    ).normalized()
    for name in ACTIVE_ENGINES:
        options = engine_parameters(config, name, True)["options"].split(",")
        assert options.count("edge-corrections") == 1
        assert options.count("no-normalize-HistZeta") == 1
    with pytest.raises(ValueError, match="require 3PCF"):
        RunConfig(ACTIVE_ENGINES, Path("unused"), edge_corrections=True,
                  options=("only-2pcf",)).normalized()


def test_scalar_mode_coupling_identity_and_empty_bin_policy():
    signal = np.arange(12, dtype=float).reshape(2, 2, 3).astype(complex)
    window = np.zeros((2, 2, 5), dtype=complex)
    window[:, :, 2] = 4.0
    window[1, 1, 2] = 0.0
    corrected, diagnostics = solve_scalar_mode_coupling(signal, window, 1)
    np.testing.assert_array_equal(corrected[0], signal[0]/4.0)
    np.testing.assert_array_equal(corrected[1, 0], signal[1, 0]/4.0)
    np.testing.assert_array_equal(corrected[1, 1], 0.0)
    assert diagnostics == {"empty_bins": 1, "singular_bins": 0}


def test_flatten_and_plot_contract(tmp_path):
    matrix = np.arange(16, dtype=float).reshape(4, 4)
    np.testing.assert_array_equal(flatten_radial_matrix(matrix), matrix.ravel())
    with pytest.raises(ValueError, match="square"):
        flatten_radial_matrix(np.zeros((2, 3)))
    pytest.importorskip("matplotlib")
    results = {
        "a": {"r": np.arange(4), "xi": np.arange(4), "zeta_m_0": matrix},
        "b": {"r": np.arange(4), "xi": np.arange(4)+0.1,
              "zeta_m_0": matrix+0.1},
    }
    config = RunConfig(tuple(results), tmp_path, bins=4, multipoles=2).normalized()
    names = {Path(path).name for path in make_plots(results, config)}
    assert {"two_pcf.png", "three_pcf_m0_real.png",
            "three_pcf_flattened_real.png"}.issubset(names)


def test_make_info_status_one_is_parsed(monkeypatch):
    completed = subprocess.CompletedProcess(
        args=["cballs", "options=make-info"], returncode=1,
        stdout="  SMOOTHPIVOTON = 1\n  TPCFON = 1\n",
    )
    monkeypatch.setattr(subprocess, "run", lambda *args, **kwargs: completed)
    assert discover_make_settings(Path("/tmp/cballs")) == {
        "SMOOTHPIVOTON": "1", "TPCFON": "1",
    }


def test_engine_table_uses_dual_node_terminology(capsys):
    print_engine_table(ACTIVE_ENGINES, True)
    text = capsys.readouterr().out
    assert "frontier=dual-node" in text
    assert "smooth=default-on" in text
    assert "smooth=unsupported" in text


def test_one_catalog_is_reused_by_one_available_engine(tmp_path):
    try:
        from cyballs import search_method_id
    except ImportError:
        pytest.skip("cyballs is unavailable")
    engine = next((name for name in ACTIVE_ENGINES[::2]
                   if search_method_id(name) >= 0), None)
    if engine is None:
        pytest.skip("no active OpenMP convergence engine is compiled")
    config = RunConfig(
        (engine,), tmp_path, theta_min=5.0, theta_max=100.0,
        bins=4, multipoles=2, threads=1,
        options=("only-2pcf", "no-out-Hist", "no-smooth-pivot"),
        verbose=0, verbose_log=0, plots=False,
    )
    results = run_engine_suite(unit_sphere_catalog(), config)
    assert tuple(results) == (engine,)
    summary = json.loads((tmp_path / "summary.json").read_text())
    assert summary["catalog_reads"] == 1
    assert summary["set_catalog_calls_per_process"] == 1
    assert not summary["failures"]


if __name__ == "__main__":
    raise SystemExit(pytest.main([os.fspath(Path(__file__))]))
