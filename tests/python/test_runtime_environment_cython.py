"""Successful informational startup must not try to fingerprint a catalog."""
import importlib.util
import os
from pathlib import Path
import subprocess
import sys

import pytest

ROOT = Path(__file__).resolve().parents[2]


@pytest.mark.parametrize('option', ['print-environment', 'print-env'])
def test_python_environment_report_stops_before_catalog_read(tmp_path, option):
    if importlib.util.find_spec('cyballs') is None:
        pytest.skip('requires a matching built/installed cyballs extension')
    environment = dict(os.environ, PYTHONPATH=str(ROOT), OMP_DYNAMIC='FALSE',
                       OMP_WAIT_POLICY='PASSIVE', OPENBLAS_NUM_THREADS='1')
    for name in ('CBALLS_MEMORY_BUDGET_MB', 'OMPI_COMM_WORLD_SIZE', 'PMI_SIZE', 'PMIX_RANK'):
        environment.pop(name, None)
    code = f'''from cyballs import cballs
model = cballs()
model.set({{'options': {option!r}, 'infile': 'missing-catalog.txt', 'verbose': 0}})
model.Run()
assert not model.getHistogramsAllocated()
model.clean_all()
print('PYTHON_REPORT_OK', flush=True)
'''
    result = subprocess.run([sys.executable, '-c', code], env=environment,
                            cwd=tmp_path, capture_output=True, text=True, timeout=45)
    assert result.returncode == 0, result.stdout + result.stderr
    assert 'cTreeBalls runtime environment' in result.stdout
    assert 'PYTHON_REPORT_OK' in result.stdout
    assert not list(tmp_path.iterdir())


def test_python_object_can_compute_after_an_environment_report(tmp_path):
    if importlib.util.find_spec('cyballs') is None:
        pytest.skip('requires a matching built/installed cyballs extension')
    environment = dict(os.environ, PYTHONPATH=str(ROOT), OMP_DYNAMIC='FALSE',
                       OMP_WAIT_POLICY='PASSIVE', OPENBLAS_NUM_THREADS='1')
    for name in ('CBALLS_MEMORY_BUDGET_MB', 'OMPI_COMM_WORLD_SIZE', 'PMI_SIZE', 'PMIX_RANK'):
        environment.pop(name, None)
    code = '''from cyballs import cballs
model = cballs()
model.set({'options': 'print-env', 'verbose': 0})
model.Run()
assert not model.getHistogramsAllocated()
model.set({'options': 'no-smooth-pivot,no-out-Hist', 'nbody': 8,
           'sizeHistN': 3, 'mChebyshev': 3, 'numberThreads': 1, 'rootDir': 'Output'})
model.Run()
assert model.getHistogramsAllocated()
assert len(model.getHistXi2pcf()) == 3
model.clean_all()
print('PYTHON_REUSE_OK', flush=True)
'''
    result = subprocess.run([sys.executable, '-c', code], env=environment,
                            cwd=tmp_path, capture_output=True, text=True, timeout=45)
    assert result.returncode == 0, result.stdout + result.stderr
    assert 'PYTHON_REUSE_OK' in result.stdout
