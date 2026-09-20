"""Publication boundaries for the maintained active-addon source snapshot."""
from pathlib import Path
import json
import re

ROOT = Path(__file__).resolve().parents[2]
ACTIVE = {
    'kdtree_2balls_omp', 'kdtree_2balls_mpi', 'balltree_2balls_omp', 'balltree_2balls_mpi',
    'octree_2balls_omp', 'octree_2balls_mpi', 'lya_forest_omp', 'lya_forest_mpi',
    'octree_shear_sphere_2balls_omp', 'kdtree_shear_sphere_2balls_omp',
    'balltree_shear_sphere_2balls_omp', 'kdtree_box_omp', 'neighbor_boxes_omp',
    'octree_3pcf_3d_omp', 'octree_3pcf_3d_mpi', 'gadget_io', 'class_lib', 'iolib',
    'cfitsio', 'pxd', 'cmdline_defs_settings',
}
SHARED = {'addons_include', 'balltree_shared', 'kdtree_shared', 'scalar_compat_shared',
          'shear_sphere_shared', 'shear_sphere_binary_2balls', 'native_octree_pair'}


def test_only_active_addons_and_shared_support_are_shipped():
    actual = {path.name for path in (ROOT/'addons').iterdir() if path.is_dir()}
    assert actual == ACTIVE | SHARED
    assert not list((ROOT/'python').glob('*.py'))
    for family in ('kappa', 'shear', 'lya'):
        assert (ROOT/f'tests/python/{family}_corr_all_engines.py').is_file()
        assert (ROOT/f'tests/python/README_{family}_corr_all_engines.md').is_file()
    assert not list((ROOT/'tests/python').glob('*_TC.py'))


def test_public_terminology_and_notebooks():
    forbidden = 'tree' + 'corr'
    for folder in ('addons', 'docs', 'include', 'source', 'main', 'tests', 'examples', 'python'):
        for path in (ROOT/folder).rglob('*'):
            if not path.is_file() or any(p in {'_build', '__pycache__'} for p in path.parts):
                continue
            if path.suffix not in {'.h', '.c', '.py', '.md', '.rst', '.txt', '.ipynb', '.m', '.1', '.html'}:
                continue
            assert forbidden.encode() not in path.read_bytes().lower(), path
            if path.suffix == '.ipynb':
                notebook = json.loads(path.read_text())
                assert notebook['nbformat'] == 4
                for cell in notebook['cells']:
                    if cell['cell_type'] == 'code':
                        source = ''.join(cell['source'])
                        if not re.search(r'^\s*[!%]', source, re.M):
                            compile(source, str(path), 'exec')
    readme = (ROOT/'README.md').read_text()
    before, acknowledgment = readme.split('## Acknowledgments', 1)
    assert forbidden not in before.lower()
    assert forbidden in acknowledgment.lower()
