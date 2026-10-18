"""Exercise the environment diagnostic without initializing a science run."""
import os
from pathlib import Path
import re
import subprocess

import pytest

ROOT = Path(__file__).resolve().parents[2]
NATIVE = (
    'CBALLS_MEMORY_BUDGET_MB', 'CBALLS_SHEAR_PROFILE', 'CBALLS_SHEAR_PIVOT_TOL',
    'CBALLS_SHEAR_BIN_THETA', 'CBALLS_SCALAR_PIVOT_TOL', 'CBALLS_SCALAR_BIN_THETA',
    'CBALLS_SHEAR_TASKS', 'CBALLS_SHEAR_TASK_MEMORY_MB', 'CBALLS_TUTOR_PYTHON',
    'CBALLS_TUTOR_ROOT',
)
LAUNCH = ('OMPI_COMM_WORLD_SIZE', 'OMPI_COMM_WORLD_RANK', 'PMI_SIZE', 'PMI_RANK',
          'PMIX_RANK', 'MV2_COMM_WORLD_SIZE', 'SLURM_NTASKS', 'SLURM_PROCID')


@pytest.fixture
def binary():
    path = Path(os.environ.get('CBALLS', ROOT / 'cballs')).absolute()
    if not path.is_file():
        pytest.skip('build cballs first, or set CBALLS to the executable')
    return path


def run(binary, tmp_path, *args, overrides=None):
    environment = dict(os.environ)
    for name in (*NATIVE, *LAUNCH):
        environment.pop(name, None)
    environment.update(OMP_DYNAMIC='FALSE', OMP_WAIT_POLICY='PASSIVE')
    environment.update(overrides or {})
    return subprocess.run([str(binary), *args], cwd=tmp_path, env=environment,
                          capture_output=True, text=True, timeout=45)


@pytest.mark.parametrize('option', ['print-environment', 'print-env'])
def test_bare_report_is_read_only(binary, tmp_path, option):
    result = run(binary, tmp_path, 'options=' + option)
    assert result.returncode == 0, result.stderr
    assert result.stdout.startswith('cTreeBalls runtime environment\n')
    for name in NATIVE:
        assert f'  {name} = <unset>\n' in result.stdout
    assert 'Parsed native memory budget: 65536 MiB (68719476736 bytes) per rank.' in result.stdout
    assert 'unset CBALLS_MEMORY_BUDGET_MB' in result.stdout
    assert 'env -u CBALLS_MEMORY_BUDGET_MB' in result.stdout
    assert not list(tmp_path.iterdir())
    assert not result.stderr


@pytest.mark.parametrize('value', ['', '0', '-1', 'bad', '1.5', '99999999999999999999999'])
def test_report_survives_invalid_memory_override(binary, tmp_path, value):
    result = run(binary, tmp_path, 'options=print-env',
                 overrides={'CBALLS_MEMORY_BUDGET_MB': value, 'PMIX_RANK': 'bogus'})
    assert result.returncode == 0, result.stderr
    assert f'CBALLS_MEMORY_BUDGET_MB = "{value}"' in result.stdout
    assert 'INVALID CBALLS_MEMORY_BUDGET_MB' in result.stdout
    assert not list(tmp_path.iterdir())


def test_values_are_quoted_safe_and_allowlisted(binary, tmp_path):
    result = run(binary, tmp_path, 'options=print-env', overrides={
        'CBALLS_MEMORY_BUDGET_MB': '131072', 'CBALLS_SHEAR_PROFILE': '',
        'CBALLS_SHEAR_TASKS': '512', 'OMP_NUM_THREADS': '7',
        'CBALLS_SHEAR_PIVOT_TOL': 'bad\n\x1b[31m"\\',
        'UNRELATED_PRIVATE_VALUE': 'should_never_be_printed',
    })
    assert result.returncode == 0, result.stderr
    assert 'CBALLS_SHEAR_PROFILE = ""\n' in result.stdout
    assert 'ANY defined value enables profiling, including 0 and the empty string' in result.stdout
    assert 'CBALLS_SHEAR_TASKS = "512"\n' in result.stdout
    assert 'OMP_NUM_THREADS = "7"\n' in result.stdout
    assert '131072 MiB (137438953472 bytes)' in result.stdout
    assert r'CBALLS_SHEAR_PIVOT_TOL = "bad\n\x1b[31m\"\\"' in result.stdout
    assert '\x1b' not in result.stdout
    assert 'should_never_be_printed' not in result.stdout
    assert 'UNRELATED_PRIVATE_VALUE' not in result.stdout
    commands = [line for line in result.stdout.splitlines() if line.startswith('  unset ')]
    assert set(NATIVE) <= set(' '.join(commands).split())
    assert not set(LAUNCH) & set(' '.join(commands).split())


def test_large_environment_value_is_bounded(binary, tmp_path):
    result = run(binary, tmp_path, 'options=print-env',
                 overrides={'CBALLS_TUTOR_ROOT': 'x' * 4096})
    assert result.returncode == 0
    assert '[truncated after 1024 bytes]' in result.stdout
    assert 'x' * 1025 not in result.stdout


def test_no_stop_continues_normal_startup(binary, tmp_path):
    # A second informational command proves execution proceeds past the report
    # without needing an engine that may be absent from this build profile.
    result = run(binary, tmp_path, 'options=print-environment,no-stop,build-fingerprint',
                 'nbody=8', 'sizeHistN=3', 'mChebyshev=3', 'numberThreads=1',
                 'rootDir=Output', 'verbose=0', 'verbose_log=0')
    assert result.returncode == 0, result.stdout + result.stderr
    assert 'cTreeBalls runtime environment' in result.stdout
    assert '"source_sha256"' in result.stdout
    assert (tmp_path / 'Output').is_dir()


def test_extra_arguments_are_checked(binary, tmp_path):
    result = run(binary, tmp_path, 'options=print-env', 'notAParameter=1')
    assert result.returncode != 0
    assert 'notAParameter' in result.stdout + result.stderr


def test_normal_startup_report_stops_before_catalog_read(binary, tmp_path):
    result = run(binary, tmp_path, 'options=print-env', 'infile=missing-catalog.txt', 'verbose=0')
    assert result.returncode == 0, result.stdout + result.stderr
    assert 'cTreeBalls runtime environment' in result.stdout
    assert not (tmp_path / 'Output').exists()


def test_report_includes_compiled_reuse_defaults(binary, tmp_path):
    report = run(binary, tmp_path, 'options=print-env')
    make_info = run(binary, tmp_path, 'options=make-info')
    for name in ('CBALLS_SHEAR_PIVOT_TOL', 'CBALLS_SHEAR_BIN_THETA',
                 'CBALLS_SCALAR_PIVOT_TOL', 'CBALLS_SCALAR_BIN_THETA'):
        match = re.search(rf'^\s*{name}\s*=\s*(\S+)', make_info.stdout, re.M)
        assert match, make_info.stdout
        entry = re.search(rf'  {name} = [^\n]*\n    When unset: (\S+)', report.stdout)
        assert entry and float(entry[1]) == float(match[1])


def test_help_lists_both_spellings(binary, tmp_path):
    result = run(binary, tmp_path, 'options=print-options')
    assert '- print-environment [startup]' in result.stdout
    assert '- print-env [startup]' in result.stdout


def test_native_environment_inventory_stays_documented():
    # Native lookups may pass a name through a shared validator; collect both
    # direct getenv arguments and the CBALLS_* literals in those call sites.
    paths = ['general_libs/clib.c', 'source/treeload.c', 'include/ia_tutor_cli.h',
             'addons/balltree_2balls_omp/dual_node_pivot_reuse.h',
             'addons/octree_shear_sphere_2balls_omp/shear_pivot_reuse.h',
             'support/shear_kernel/shear_frontier_policy.h']
    read_controls = set()
    for path in paths:
        read_controls.update(re.findall(r'"(CBALLS_[A-Z0-9_]+)"', (ROOT / path).read_text()))
    report_source = (ROOT / 'source/runtime_environment.c').read_text()
    documented = set(re.findall(r'\{"(CBALLS_[A-Z0-9_]+)"', report_source))
    assert read_controls == documented == set(NATIVE)
