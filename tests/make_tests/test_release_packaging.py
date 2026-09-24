"""Reject incomplete, mislabeled or contaminated public source archives."""
import io
from pathlib import Path
import sys
import tarfile
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'scripts'))
from check_sdist import REQUIRED, inspect_archive


def archive(tmp_path, *, omit=None, extra=None, distribution='cyballs'):
    path = tmp_path/'cyballs-1.1.0.tar.gz'
    with tarfile.open(path, 'w:gz') as stream:
        for name in sorted((REQUIRED-{omit}) | ({extra} if extra else set())):
            data = (f'Name: {distribution}\nVersion: 1.1.0\n'.encode()
                    if name == 'PKG-INFO' else b'fixture\n')
            member = tarfile.TarInfo('cyballs-1.1.0/'+name)
            member.size = len(data)
            stream.addfile(member, io.BytesIO(data))
    return path


def test_complete_public_source_policy(tmp_path):
    assert inspect_archive(archive(tmp_path))['status'] == 'PASS'


@pytest.mark.parametrize('missing', [
    'addons/lya_forest_omp/search_lya_forest_los_tree_omp.c',
    'python/ccyballs.pxd.in', 'requirements/build.txt'])
def test_missing_build_inputs_fail(tmp_path, missing):
    with pytest.raises(AssertionError, match='missing public source'):
        inspect_archive(archive(tmp_path, omit=missing))


@pytest.mark.parametrize('extra', ['addons/gsl/config.h',
    'addons/cfitsiolib/cfitsio_ver_4.6.3/fitscore.c', 'cyballs.so', '../escape',
    'addons/python_env/private.py', 'addons/octree_ggg_omp/search.c',
    'tests/python/kappa_corr_all_engines_TC.py'])
def test_contaminated_source_archives_fail(tmp_path, extra):
    with pytest.raises(AssertionError):
        inspect_archive(archive(tmp_path, extra=extra))


def test_distribution_name_matches_install_documentation(tmp_path):
    with pytest.raises(AssertionError, match='unexpected distribution'):
        inspect_archive(archive(tmp_path, distribution='cTreeBalls'))
