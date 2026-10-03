"""Capabilities fail closed; shared changes select all dependent references."""
import copy
import json
from pathlib import Path
import subprocess
import sys
import pytest

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from generate_capabilities import validate, render
from capabilities_generated import CAPABILITIES, ENGINES, gate_plan
from affected_regressions import select


def test_generated_catalogue_is_current():
    subprocess.run([sys.executable,ROOT/'scripts/generate_capabilities.py','--check'],check=True)


def test_new_engine_requires_visible_test_declaration():
    data=copy.deepcopy(CAPABILITIES)
    e=copy.deepcopy(data['engines'][0]);e.update(name='new-engine',id=999,aliases=[])
    del e['gate']
    data['engines'].append(e)
    with pytest.raises(ValueError,match='missing gate'):validate(data)
    e['gate']=dict(oracle='scalar',tests=[],parallel='omp')
    with pytest.raises(ValueError,match='declared together'):validate(data)
    e['gate']=dict(oracle=None,tests=[],skip_reason='Awaiting an independent oracle')
    validate(data)
    # Legacy or undeclared engines cannot silently enter the release matrix.
    assert all(e['gate'].get('oracle') for e in ENGINES.values())
    active = next(iter(ENGINES.values()))
    with pytest.raises(ValueError,match='no matching retained'):
        gate_plan({active['name']:999999})


@pytest.mark.parametrize('path',['source/smooth_pivots.c','source/common_histogram.c',
    'include/scalar_moments.h','source/mpi_runtime.c','unknown/new_kernel.c'])
def test_shared_or_unowned_changes_select_full_public_regressions(path):
    plan=select([path])
    assert set(plan['engines'])=={n for n,e in ENGINES.items() if e['gate'].get('oracle')}
    assert plan['engines']
    assert 'tests/python/test_lya_forest_los_tree.py' in plan['tests']


def test_shared_shear_implementation_selects_all_tree_families():
    plan=select(['support/shear_binary/search_shear_sphere_binary_2balls_impl.h'])
    assert {f'{t}-shear-sphere-2balls-omp' for t in ('octree','kdtree','balltree')} <= set(plan['engines'])


def test_reference_numerics_are_not_generated():
    outputs=render(CAPABILITIES)
    assert not any(str(p).startswith('tests/') for p in outputs)
    assert Path('include/engine_registry_generated.h') in outputs


def test_removed_alternate_addons_are_not_selected():
    plan=select(['addons/balltree_2balls_omp/search_balltree_2balls_omp.c'])
    assert plan['alternate_regressions'] == []
    assert 'balltree-2balls-mpi_3pcf' not in plan['engines']


def test_relative_shared_kernel_includes_select_all_scalar_tree_families():
    plan=select(['addons/balltree_2balls_omp/search_balltree_2balls_omp.c'])
    assert {f'{tree}-2balls-{parallel}' for tree in ('kdtree','balltree','octree')
            for parallel in ('omp','mpi')} <= set(plan['engines'])
    assert 'addons/kdtree_2balls_omp/search_kdtree_2balls_omp.c' in plan['include_closure']
    assert not any('/../' in path for path in plan['include_closure'])

@pytest.mark.parametrize('name',['dual_node_radial_bins.h','dual_node_pair_acceptance.h','dual_node_task_schedule.h','dual_node_multipole.h','dual_node_search_policy.h'])
def test_extracted_scalar_contract_selects_all_six_engines(name):
    plan=select(['addons/balltree_2balls_omp/'+name])
    expected={f'{tree}-2balls-{parallel}' for tree in ('octree','kdtree','balltree') for parallel in ('omp','mpi')}
    assert expected <= set(plan['engines'])
    assert 'tests/python/test_two_ball_edge_corrections.py' in plan['tests']


@pytest.mark.parametrize('path,engines',[
    ('addons/kdtree_shear_sphere_2balls_omp/kdtree_shear_sphere_tree.c',
     {'kdtree-shear-sphere-2balls-mpi'}),
    ('support/pca_tree/fcfc_balltree.c',
     {'kdtree-shear-sphere-2balls-mpi','balltree-shear-sphere-2balls-mpi'}),
])
def test_mpi_shear_shared_tree_builders_select_mpi_regressions(path,engines):
    plan=select([path])
    if not engines <= set(plan['engines']):
        raise AssertionError('shared shear tree change omitted its MPI consumers')
    if 'tests/python/test_shear_sphere_mpi.py' not in plan['tests']:
        raise AssertionError('shared shear tree change omitted the MPI numerical regression')
