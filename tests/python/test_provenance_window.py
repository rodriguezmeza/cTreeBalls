"""Permanent immutable-settings and scalar-window support regressions."""
import json
from pathlib import Path

import numpy as np
import pytest

from cyballs import cballs, search_method_id
from test_two_ball_edge_corrections import catalog, brute_force, edge_solution

ROOT = Path(__file__).resolve().parents[2]
RECORDED = json.loads((ROOT / "tests/fixtures/provenance_window/recorded_cases.json").read_text())
ENGINES = ("kdtree-2balls-omp", "balltree-2balls-omp", "octree-2balls-omp")
PATHS = [(name, False) for name in ENGINES] + [(name, True) for name in ENGINES]


def parameters(root, engine=ENGINES[0], legacy=False, edge=True, output=False):
    options = ["KKKCorrelation", "compute-HistN", "weights-norm", "read-mask",
               "no-smooth-pivot", "no-normalize-HistZeta"]
    options += ["legacy-one-ball", "no-one-ball"] if legacy else ["no-two-balls"]
    if edge:
        options += ["edge-corrections"]
    if not output:
        options += ["no-out-Hist"]
    return dict(searchMethod=engine, sizeHistN=4, mChebyshev=2, sizeHistPhi=8,
                rangeN=1.5, rminHist=.02, lengthBox=2, numberThreads=1,
                usePeriodic=False, useLogHist=False, nsmooth=2,
                verbose=0, verbose_log=0, rootDir=str(root), options=",".join(options))


def register(balls, data):
    balls.set_catalog(data[0], kappa=data[1], weights=data[2], mask=data[3])


@pytest.fixture
def balls():
    obj = cballs()
    yield obj
    obj.struct_cleanup()
    obj.clear_catalogs()


def test_immutable_settings_and_recorded_cache_invalidation(balls, tmp_path):
    case = RECORDED["provenance"]
    register(balls, catalog())
    balls.set(parameters(tmp_path))
    assert balls.run_settings is None
    pars = balls.pars
    balls.Run()
    recorded = balls.run_settings
    assert recorded["effective"]["sizeHistN"] == case["initial_bins"]
    assert recorded["effective"]["catalog_sizes"] == (64,)
    assert recorded["effective"]["numberThreads"] == 1
    for view in (pars, balls.pars, recorded["requested"], recorded["effective"]):
        with pytest.raises(TypeError):
            view["sizeHistN"] = case["next_bins"]
    with pytest.raises(TypeError):
        recorded["version"] = "changed"
    balls.set(sizeHistN=case["initial_bins"])
    balls.Run()
    assert balls.run_settings is recorded  # A true no-op can reuse the run.
    balls.set(sizeHistN=case["next_bins"])
    assert balls.run_settings is None
    with pytest.raises(Exception):
        balls.getHistNN()
    balls.Run()
    assert len(balls.getHistNN()) == case["next_bins"]
    assert balls.run_settings["effective"]["sizeHistN"] == case["next_bins"]
    assert pars["sizeHistN"] == recorded["requested"]["sizeHistN"] == case["initial_bins"]
    balls.struct_cleanup()
    assert balls.run_settings["effective"]["sizeHistN"] == case["next_bins"]
    balls.clean()
    assert balls.run_settings is None
    assert recorded["effective"]["sizeHistN"] == case["initial_bins"]


def test_registration_owns_arrays_and_parameter_values(balls, tmp_path):
    class MutableText:
        value = 4
        def __str__(self):
            return str(self.value)
    text = MutableText()
    data = tuple(array.copy() for array in catalog())
    reference = tuple(array.copy() for array in data)
    register(balls, data)
    balls.set(parameters(tmp_path))
    balls.set(sizeHistN=text)
    text.value = 6
    data[0][:] = 0  # Would otherwise become a degenerate catalog.
    data[1][:] = 0
    data[2][:] = 0
    data[3][:] = 0
    balls.Run()
    assert balls.run_settings["requested"]["sizeHistN"] == "4"
    expected = edge_solution(*brute_force(reference))
    actual = np.array([balls.getHistZetaM_EE_complex(m+1) for m in range(3)])
    np.testing.assert_allclose(actual, expected, rtol=2e-10, atol=2e-10)
    # Explicit registration replaces the owned snapshot and invalidates results.
    replacement = list(reference)
    replacement[1] = np.zeros(64)
    register(balls, replacement)
    assert balls.run_settings is None
    with pytest.raises(Exception):
        balls.getScalarWindowDiagnostics()
    balls.Run()
    valid = balls.getScalarWindowDiagnostics()["valid"]
    assert np.any(valid)
    np.testing.assert_array_equal(balls.getHistZetaM_EE_complex(1)[valid], 0)


def test_staged_run_failed_run_and_endrun_snapshots(balls, tmp_path):
    register(balls, catalog())
    balls.set(parameters(tmp_path))
    balls.Run(level=["SetNumberThreads"])
    assert balls.run_settings is None
    balls.Run(level=["MainLoop"])
    balls.struct_cleanup()
    balls.Run(level=["SetNumberThreads"])
    assert balls.run_settings is None
    balls.Run()
    snapshot = balls.run_settings
    balls.Run(level=["EndRun"])
    assert balls.run_settings is snapshot
    with pytest.raises(Exception):
        balls.getScalarWindowDiagnostics()
    balls.set(sizeHistN=0)
    with pytest.raises(Exception):
        balls.Run()
    assert balls.run_settings is None
    assert snapshot["effective"]["sizeHistN"] == 4
    balls.set(sizeHistN=4)
    balls.Run()
    assert balls.run_settings["effective"]["sizeHistN"] == 4


@pytest.mark.parametrize("engine,legacy", PATHS)
def test_support_and_valid_zero_across_shared_callers(balls, tmp_path, engine, legacy):
    if search_method_id(engine) < 0:
        pytest.skip(f"{engine} not built")
    data = list(catalog())
    data[1] = np.zeros(64)
    register(balls, data)
    balls.set(parameters(tmp_path, engine, legacy))
    if engine == "balltree-2balls-omp" and legacy:
        with pytest.raises(Exception, match="legacy-one-ball does not support edge-corrections"):
            balls.Run()
        assert balls.run_settings is None
        return
    balls.Run()
    diagnostic = balls.getScalarWindowDiagnostics()
    valid = diagnostic["valid"]
    assert np.any(valid)
    np.testing.assert_array_equal(valid, diagnostic["status"] == 1)
    assert np.all(diagnostic["window_monopole"][valid] > 0)
    assert np.all((diagnostic["pivot_ratio"][valid] > 0) & (diagnostic["pivot_ratio"][valid] <= 1))
    for m in range(1, 4):
        value = balls.getHistZetaM_EE_complex(m)
        np.testing.assert_array_equal(value[valid], 0)
        assert np.all(np.isnan(value[~valid]))
    diagnostic["status"][:] = 255
    assert not np.any(balls.getScalarWindowDiagnostics()["status"] == 255)
    # Only one visible point: no distinct triple can support any estimate.
    data[3][:] = 0
    data[3][0] = 1
    register(balls, data)
    balls.Run()
    diagnostic = balls.getScalarWindowDiagnostics()
    assert not np.any(diagnostic["valid"])
    np.testing.assert_array_equal(diagnostic["status"], 2)
    np.testing.assert_array_equal(diagnostic["window_monopole"], 0)
    assert np.all(np.isnan(diagnostic["pivot_ratio"]))
    assert np.all(np.isnan(balls.getHistZetaM_EE_complex(1)))


@pytest.mark.parametrize("legacy", (False, True))
def test_diagnostics_writer_failure_and_recovery(balls, tmp_path, legacy):
    register(balls, catalog())
    # The legacy octree exercises its independent Gauss-Jordan publication path.
    balls.set(parameters(tmp_path, "octree-2balls-omp", legacy, output=True))
    path = tmp_path / "histZetaM_window_diagnostics.txt"
    path.mkdir()
    with pytest.raises(Exception, match="window_diagnostics"):
        balls.Run()
    assert balls.run_settings is None
    with pytest.raises(Exception):
        balls.getScalarWindowDiagnostics()
    path.rmdir()
    balls.Run()
    diagnostic = balls.getScalarWindowDiagnostics()
    saved = np.loadtxt(path)
    np.testing.assert_array_equal(saved[:, :2], [(i+1, j+1) for i in range(4) for j in range(4)])
    for column, name in ((2, "status"), (3, "window_monopole"), (4, "pivot_ratio")):
        np.testing.assert_array_equal(saved[:, column].reshape(4, 4), diagnostic[name])
    balls.set(parameters(tmp_path, edge=False))
    balls.Run()
    with pytest.raises(Exception, match="unavailable"):
        balls.getScalarWindowDiagnostics()
    with pytest.raises(Exception, match="enable edge-corrections"):
        balls.getHistZetaM_EE_complex(1)


def test_two_live_settings_and_window_allocations_are_independent(balls, tmp_path):
    other = cballs()
    try:
        register(balls, catalog())
        balls.set(parameters(tmp_path / "first"))
        balls.Run()
        snapshot = balls.run_settings
        diagnostics = balls.getScalarWindowDiagnostics()
        register(other, catalog(count=32))
        other.set(dict(parameters(tmp_path / "second"), sizeHistN=6))
        other.Run()
        assert snapshot["effective"]["catalog_sizes"] == (64,)
        assert other.run_settings["effective"]["catalog_sizes"] == (32,)
        other.struct_cleanup()
        assert balls.run_settings is snapshot
        for name, value in balls.getScalarWindowDiagnostics().items():
            np.testing.assert_array_equal(value, diagnostics[name])
    finally:
        other.struct_cleanup()
        other.clear_catalogs()


@pytest.mark.parametrize("mmax", (2, 4, 7))
def test_legacy_octree_constant_field_has_analytic_nonzero_solution(balls, tmp_path, mmax):
    rng = np.random.default_rng(531)
    x, y = rng.uniform(-.45, .45, (2, 160))
    positions = np.column_stack((x, y, np.ones_like(x)))
    positions /= np.linalg.norm(positions, axis=1)[:, None]
    weights = rng.uniform(.7, 1.3, len(x))
    register(balls, (positions, np.full(len(x), .7), weights, np.ones(len(x), dtype=np.uint8)))
    balls.set(dict(parameters(tmp_path, "octree-2balls-omp", legacy=True),
                   rangeN=.75, rminHist=.03, sizeHistN=3, mChebyshev=mmax,
                   useLogHist=True, sizeHistPhi=16))
    balls.Run()
    assert np.all(balls.getScalarWindowDiagnostics()["valid"])
    np.testing.assert_allclose(balls.getHistZetaM_EE_complex(1), .7**3, rtol=0, atol=5e-13)
    for order in range(2, mmax+2):
        np.testing.assert_allclose(balls.getHistZetaM_EE_complex(order), 0, rtol=0, atol=5e-13)
