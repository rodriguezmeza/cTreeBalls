"""Reject incomplete, mislabeled or contaminated public source archives."""
import io
from pathlib import Path
import sys
import tarfile
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'scripts'))
from check_sdist import REQUIRED, inspect_archive
from release_profile import EXTERNAL_DEFAULTS, normalize_release_tree, check_external_defaults


def archive(tmp_path, *, omit=None, extra=None, distribution='cyballs', overrides=None):
    path = tmp_path/'cyballs-1.1.0.tar.gz'
    with tarfile.open(path, 'w:gz') as stream:
        for name in sorted((REQUIRED-{omit}) | ({extra} if extra else set())):
            data = (f'Name: {distribution}\nVersion: 1.1.0\n'.encode()
                    if name == 'PKG-INFO' else b'fixture\n')
            if name in EXTERNAL_DEFAULTS:
                key, value = EXTERNAL_DEFAULTS[name]
                data = f'{key} = {value}\n'.encode()
            if overrides and name in overrides:
                data = overrides[name].encode()
            member = tarfile.TarInfo('cyballs-1.1.0/'+name)
            member.size = len(data)
            stream.addfile(member, io.BytesIO(data))
    return path


def test_complete_public_source_policy(tmp_path):
    assert inspect_archive(archive(tmp_path))['status'] == 'PASS'


@pytest.mark.parametrize('missing', [
    'tests/python/shear_fits_catalog.py', 'tests/python/shear_products.py',
    'tests/python/run_desy3_shear_catalogs.py', 'tests/python/test_shear_fits_catalog.py',
    'include/cli_contracts.h', 'python/qualification_api.pxi',
    'scripts/qualify_results.py', 'scripts/benchmark_representative.py',
    'tests/python/test_clustered_contracts.py', 'tests/python/test_mpi_boundary_contracts.py',
    'ENGINEERING_COMPLETION.md',
    'addons/lya_forest_omp/search_lya_forest_los_tree_omp.c',
    'addons/lya_forest_omp/lya_pivot_frontier.h',
    'include/input_contracts.h',
    'tests/python/test_input_failure_contracts.py',
    'addons/balltree_2balls_omp/dual_node_pivot_reuse.h',
    'tests/python/test_scalar_pivot_reuse.py',
    'scripts/benchmark_scalar_pivot_reuse.py',
    'docs/SCALAR_HIERARCHICAL_REUSE.md',
    'addons/lya_forest_omp/lya_radial_moments.h',
    'tests/python/test_lya_hierarchy.py',
    'scripts/benchmark_lya_hierarchy.py',
    'docs/LYA_HIERARCHICAL_REUSE.md',
    'docs/LYA_LOS_HIERARCHICAL_REUSE.md',
    'docs/LYA_MULTIPOLE_HIERARCHICAL_REUSE.md',
    'tests/python/test_lya_multipole_hierarchy.py',
    'tests/python/test_lya_multipole_reconstruction.py',
    'docs/LYA_MULTIPOLE_RECONSTRUCTION.md',
    'addons/lya_forest_omp/lya_triplet_multipole.h',
    'addons/lya_forest_omp/lya_los_moments.h',
    'tests/python/test_lya_los_hierarchy.py',
    'python/ccyballs.pxd.in', 'requirements/build.txt'])
def test_missing_build_inputs_fail(tmp_path, missing):
    with pytest.raises(AssertionError, match='missing public source'):
        inspect_archive(archive(tmp_path, omit=missing))


@pytest.mark.parametrize('extra', ['addons/gsl/config.h',
    'addons/cfitsiolib/cfitsio_ver_4.6.3/fitscore.c', 'cyballs.so', '../escape'])
def test_contaminated_source_archives_fail(tmp_path, extra):
    with pytest.raises(AssertionError):
        inspect_archive(archive(tmp_path, extra=extra))


def test_distribution_name_matches_install_documentation(tmp_path):
    with pytest.raises(AssertionError, match='unexpected distribution'):
        inspect_archive(archive(tmp_path, distribution='cTreeBalls'))


@pytest.mark.parametrize('name', list(EXTERNAL_DEFAULTS))
@pytest.mark.parametrize('value', ['1', '-1', '$(LOCAL_DEFAULT)', '0 ', '0 # comment',
                                  '0\n{key} = 1', '0\n{key} += 0'])
def test_incompatible_or_ambiguous_packaged_defaults_fail(tmp_path, name, value):
    key, _ = EXTERNAL_DEFAULTS[name]
    replacement = f'{key} = {value.format(key=key)}\n'
    with pytest.raises(AssertionError, match='external dependency profile'):
        inspect_archive(archive(tmp_path, overrides={name: replacement}))


def test_release_normalization_does_not_modify_hardlinked_checkout(tmp_path):
    import os
    checkout, stage = tmp_path/'checkout', tmp_path/'stage'
    original = {}
    for name, (key, _) in EXTERNAL_DEFAULTS.items():
        src, dst = checkout/name, stage/name
        src.parent.mkdir(parents=True, exist_ok=True)
        dst.parent.mkdir(parents=True, exist_ok=True)
        original[name] = '# developer profile\n'+key+' = 1\nUNCHANGED = 17\n'
        src.write_text(original[name]); os.link(src, dst)
    normalize_release_tree(stage)
    once = {name: (stage/name).read_text() for name in EXTERNAL_DEFAULTS}
    assert check_external_defaults(once) == {'GSLINTERNAL': '0', 'CFITSIOLIBON': '0'}
    normalize_release_tree(stage)
    for name in EXTERNAL_DEFAULTS:
        assert (checkout/name).read_text() == original[name]
        assert (stage/name).read_text() == once[name]
        assert 'UNCHANGED = 17' in once[name]
        assert (checkout/name).stat().st_ino != (stage/name).stat().st_ino


def test_ambiguous_staging_rejected_before_either_file_changes(tmp_path):
    for name, (key, _) in EXTERNAL_DEFAULTS.items():
        p = tmp_path/name; p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(f'{key} = 1\n')
    broken = tmp_path/'addons/Makefile_addons_settings'
    broken.write_text(broken.read_text()+'CFITSIOLIBON = 0\n')
    with pytest.raises(RuntimeError, match='ambiguous'):
        normalize_release_tree(tmp_path)
    assert (tmp_path/'Makefile_settings').read_text() == 'GSLINTERNAL = 1\n'


def test_staged_defaults_select_external_branches_in_make(tmp_path):
    import subprocess
    for name, (key, _) in EXTERNAL_DEFAULTS.items():
        p = tmp_path/name; p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(f'{key} = 1\n')
    normalize_release_tree(tmp_path)
    (tmp_path/'Makefile').write_text(
        'include Makefile_settings\ninclude addons/Makefile_addons_settings\n'
        'ifeq ($(GSLINTERNAL),0)\nGSL_SELECTED = external\nendif\n'
        'ifeq ($(CFITSIOLIBON),0)\nFITS_SELECTED = external\nendif\n'
        'all:\n\t@echo "$(GSL_SELECTED) $(FITS_SELECTED)"\n')
    proc = subprocess.run(['make', '--no-print-directory', '-s'], cwd=tmp_path,
                          text=True, capture_output=True, timeout=30)
    assert proc.returncode == 0, proc.stdout+proc.stderr
    assert proc.stdout.strip() == 'external external'


@pytest.mark.parametrize('library,overrides', [
    ('GSL', ['GSL_CONFIG=false', 'GSL_LIB=', 'GSL_INCLUDE=', 'CFITSIO_LIB=/example/fits/lib', 'CFITSIO_INCLUDE=/example/fits/include']),
    ('CFITSIO', ['PKG_CONFIG=false', 'CFITSIO_LIB=', 'CFITSIO_INCLUDE=', 'GSL_LIB=/example/gsl/lib', 'GSL_INCLUDE=/example/gsl/include']),
])
def test_missing_native_dependency_has_actionable_error(library, overrides):
    import subprocess
    root = Path(__file__).resolve().parents[2]
    proc = subprocess.run(['make', '--no-print-directory', '-s', 'GSLINTERNAL=0',
        'CFITSIOLIBON=0', 'SLEEFON=0', *overrides, 'print-cyballs-build-env'],
        cwd=root, text=True, capture_output=True, timeout=30)
    assert proc.returncode != 0
    assert f'External {library} is required' in proc.stdout+proc.stderr


def test_resolved_dependency_flags_are_exported_for_cython():
    import subprocess
    root = Path(__file__).resolve().parents[2]
    proc = subprocess.run(['make', '--no-print-directory', '-s', 'GSLINTERNAL=0',
        'CFITSIOLIBON=0', 'SLEEFON=0', 'GSL_LIB=/example/gsl/lib',
        'GSL_INCLUDE=/example/gsl/include', 'CFITSIO_LIB=/example/fits/lib',
        'CFITSIO_INCLUDE=/example/fits/include', 'print-cyballs-build-env'],
        cwd=root, text=True, capture_output=True, timeout=30)
    assert proc.returncode == 0, proc.stdout+proc.stderr
    values = dict(line.split('=', 1) for line in proc.stdout.splitlines() if line.startswith('__CBALLS_'))
    assert values['__CBALLS_GSL_CFLAGS__'] == '-I/example/gsl/include'
    assert values['__CBALLS_GSL_LDFLAGS__'] == '-L/example/gsl/lib'
    assert values['__CBALLS_GSL_LIBS__'] == '-lgsl -lgslcblas'
    assert values['__CBALLS_CFITSIO_CFLAGS__'] == '-I/example/fits/include'
    assert values['__CBALLS_CFITSIO_LDFLAGS__'] == '-L/example/fits/lib'
    assert values['__CBALLS_CFITSIO_LIBS__'] == '-lcfitsio'
