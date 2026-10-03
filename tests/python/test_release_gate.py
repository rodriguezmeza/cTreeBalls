"""Release identity, fail-closed runner, and persistent run metadata contracts."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import subprocess
import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'scripts'))
from active_release_gate import Gate, expected_registry, parse_registry
from cyballs import build_info, cballs, CosmoComputationError
from test_provenance_window import parameters, register
from test_two_ball_edge_corrections import catalog


def test_resolved_registry_is_explicit_and_complete():
    text = 'Searching methods registered in this executable (1):\n- octree-sincos-omp (id=24)\n'
    assert parse_registry(text) == expected_registry({})
    with pytest.raises(AssertionError): parse_registry(text.replace('(1):','(2):'))
    with pytest.raises(AssertionError): parse_registry('')
    assert parse_registry(text.replace('octree-sincos-omp','unknown')) != expected_registry({})


def test_current_native_registry_matches_the_resolved_gate_matrix(tmp_path):
    result = subprocess.run([str(ROOT/'cballs'), 'options=print-search-methods',
                             f'rootDir={tmp_path}', 'verbose=0', 'verbose_log=0'],
                            text=True, capture_output=True, timeout=30)
    assert result.returncode == 1, result.stdout + result.stderr
    assert parse_registry(result.stdout) == expected_registry(build_info()['resolved_settings'])


def test_runner_fails_on_command_failure_timeout_and_stale_output(tmp_path):
    args = argparse.Namespace(output=tmp_path/'gate', jobs=1, mpi_command='mpiexec')
    gate = Gate(args)
    with pytest.raises(RuntimeError):
        gate.command('failure', [sys.executable,'-c','raise SystemExit(3)'])
    with pytest.raises(RuntimeError):
        gate.command('timeout', [sys.executable,'-c','import time; time.sleep(10)'], timeout=.05)
    assert all(c['status']=='FAIL' for c in json.loads((args.output/'gate.json').read_text())['commands'])
    with pytest.raises(FileExistsError): Gate(args)


@pytest.mark.parametrize('logarithmic', [False, True])
def test_run_metadata_full_precision_immutable_and_saved(tmp_path, logarithmic):
    model = cballs()
    data = catalog()
    try:
        model.set(parameters(tmp_path/'résultat', output=True) | {'useLogHist': logarithmic})
        register(model, data)
        model.Run()
        result = model.getRunMetadata()
        native = json.loads((tmp_path/'résultat/run-metadata.json').read_text())
        assert result['build'] == native['build'] == build_info()
        assert len(result['build']['id']) == 64
        edges = np.geomspace(.02,1.5,5) if logarithmic else np.linspace(.02,1.5,5)
        np.testing.assert_allclose(result['bin_edges']['radial'], edges, rtol=1e-15, atol=1e-16)
        assert result['precision']['storage_bits'] == 64
        assert result['parallel']['threads_requested'] == result['parallel']['openmp_probe_threads'] == 1
        assert result['weights']['weights_norm_option'] and result['masks']['read_mask']
        assert not result['effective_smoothing']['enabled']
        assert result['inputs']['memory_catalogs'][0]['weights']['sha256'] == hashlib.sha256(data[2].tobytes()).hexdigest()
        with pytest.raises(TypeError): model.run_settings['provenance']['build']['id'] = 'corrupt'
        result['build']['id'] = 'detached'
        assert model.getRunMetadata()['build']['id'] == build_info()['id']
        model.struct_cleanup()
        assert model.getRunMetadata()['engine'] == 'kdtree-2balls-omp'
    finally:
        model.struct_cleanup()


def test_metadata_output_failure_propagates_and_recovers(tmp_path):
    model = cballs()
    root = tmp_path/'output'
    root.mkdir()
    (root/'run-metadata.json').mkdir()
    try:
        model.set(parameters(root, output=True))
        register(model, catalog())
        with pytest.raises(CosmoComputationError, match='run metadata'):
            model.Run()
        assert model.run_settings is None
        assert not model.getBodytableAllocated()
        (root/'run-metadata.json').rmdir()
        model.Run()
        assert json.loads((root/'run-metadata.json').read_text())['engine'] == 'kdtree-2balls-omp'
    finally:
        model.struct_cleanup()
