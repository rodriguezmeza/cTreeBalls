#!/usr/bin/env python3
"""Check the public cyballs source archive, without extracting untrusted paths."""
import argparse
from email.parser import BytesParser
import hashlib
import json
from pathlib import Path, PurePosixPath
import tarfile
from release_profile import EXTERNAL_DEFAULTS, check_external_defaults

REQUIRED = {
    'capabilities/build_profile.json', 'scripts/generate_profile_docs.py',
    'tests/python/benchmark_timing.py', 'tests/python/dual_node_compat.py',
    'support/pca_tree/Makefile_balltree_omp', 'support/kd_tree/Makefile_kdtree_omp',
    'docs/ACTIVE_PROFILE.md', 'docs/active_profile.rst',
    'scripts/release_profile.py', 'scripts/verify_sdist_install.py',
    'addons/Makefile_addons_settings', 'addons/Makefile_addons',
    'docs/SOURCE_DISTRIBUTION.md', 'docs/RESULT_AVAILABILITY_AND_IDENTITY.md',
    'source/scalar_window_io.c', 'include/scalar_window_solver.h',
    'tests/python/test_saved_scalar_edge.py', 'docs/saved_scalar_edge.md',
    'tests/python/shear_fits_catalog.py', 'tests/python/shear_products.py',
    'tests/python/run_desy3_shear_catalogs.py',
    'tests/python/test_shear_fits_catalog.py',
    'include/cli_contracts.h','python/qualification_api.pxi',
    'scripts/qualify_results.py','scripts/benchmark_representative.py',
    'tests/python/test_clustered_contracts.py','tests/python/test_mpi_boundary_contracts.py','ENGINEERING_COMPLETION.md',
    'include/input_contracts.h',
    'tests/python/test_input_failure_contracts.py',
    'RELEASE_WORK_SEQUENCE.md', 'python/resource_api.pxi', 'scripts/resource_plan.py',
    'scripts/workload_acceptance.py','scripts/benchmark_scaling.py',
    'tests/python/test_scientific_qualification.py',
    'addons/balltree_2balls_omp/dual_node_radial_bins.h',
    'addons/balltree_2balls_omp/dual_node_pair_acceptance.h',
    'addons/balltree_2balls_omp/dual_node_task_schedule.h',
    'addons/balltree_2balls_omp/dual_node_multipole.h',
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
    'addons/balltree_2balls_omp/dual_node_search_policy.h',
    'tests/python/test_release_verification_contracts.py',
    'tests/python/test_owned_results_resources.py',
    'tests/python/test_lya_pair_cells.py',

    'capabilities/engines.json','include/engine_registry_generated.h','include/mpi_runtime.h','include/tree_workspace.h',
    'source/engine_registry.c','source/runtime_context.c','source/memory_catalog.c',
    'source/common_histogram.c','source/smooth_pivots.c','source/mpi_runtime.c',
    'scripts/generate_capabilities.py','scripts/capabilities_generated.py',
    'scripts/affected_regressions.py','scripts/benchmark_contracts.py',
    'tests/test_runtime_context.c','tests/test_mpi_runtime.c',
    'tests/python/test_capability_contracts.py','ENGINE_CAPABILITIES.md','MAINTAINABILITY_CONTRACTS.md',

    'include/resource_contracts.h', 'scripts/accuracy_acceptance.py',
    'tests/test_resource_contracts.c', 'tests/python/test_resource_scientific_contracts.py',
    'RESOURCE_SCIENTIFIC_CONTRACTS.md', 'PKG-INFO', 'setup.py', 'pyproject.toml', 'MANIFEST.in', 'Makefile',
    'Makefile_machine', 'Makefile_settings', 'requirements/build.txt',
    'requirements/release.txt', 'include/tree_contracts.h',
    'python/cyballs.pyx', 'python/ccyballs.pxd.in',
    'scripts/build_fingerprint.py', 'scripts/active_release_gate.py',
    'scripts/release_gate_case.py', 'scripts/check_installed_package.py',
    'addons/Makefile_gsl', 'addons/cfitsiolib/Makefile_cfitsiolib',
    'addons/lya_forest_omp/search_lya_forest_los_tree_omp.c',
    'addons/lya_forest_omp/lya_forest_los_tree.h',
    'addons/lya_forest_omp/lya_pivot_frontier.h',
    'addons/lya_forest_omp/Makefile_lya_forest_omp',
    'addons/octree_3pcf_3d_omp/search_octree_3pcf_3d_omp.c',
    'tests/python/test_lya_forest_los_tree.py',
    'tests/python/test_lya_pivot_frontier.py',
    'tests/python/benchmark_lya_pivot_frontier.py',
    'tests/python/README_benchmark_lya_pivot_frontier.md',
    'tests/python/test_cython_in_memory_catalog.py',
    'tests/python/test_histogram_availability.py',
    'tests/python/test_effective_catalog_identity.py',
    'docs/RESULT_AVAILABILITY_AND_IDENTITY.md',
    'tests/python/test_release_packaging.py',
}


def inspect_archive(path):
    with tarfile.open(path, 'r:gz') as archive:
        members = archive.getmembers()
        roots = set()
        names = set()
        files = {}
        for member in members:
            parts = PurePosixPath(member.name)
            if parts.is_absolute() or '..' in parts.parts or not parts.parts:
                raise AssertionError(f'unsafe archive path: {member.name}')
            roots.add(parts.parts[0])
            name = '/'.join(parts.parts[1:])
            if member.issym() or member.islnk() or not (member.isfile() or member.isdir()):
                raise AssertionError(f'unsupported archive member: {member.name}')
            if member.isfile():
                if name in names:
                    raise AssertionError(f'duplicate archive member: {name}')
                names.add(name)
                files[name] = member
        if not (len(roots) == 1):
            raise AssertionError(f'expected one source root: {roots}')
        if not (REQUIRED <= names):
            raise AssertionError(f'missing public source files: {sorted(REQUIRED - names)}')
        forbidden = [name for name in names if name.startswith('addons/gsl/')
                     or name.startswith('addons/cfitsiolib/cfitsio_ver_')
                     or PurePosixPath(name).suffix in {'.o', '.a', '.so', '.dylib', '.pyc'}]
        if not (not forbidden):
            raise AssertionError(f'vendored libraries or build artifacts: {forbidden}')
        defaults = check_external_defaults({name: archive.extractfile(files[name]).read().decode('utf-8')
                                            for name in EXTERNAL_DEFAULTS})
        metadata = BytesParser().parsebytes(archive.extractfile(files['PKG-INFO']).read())
        if not (metadata['Name'] == 'cyballs'):
            raise AssertionError(f"unexpected distribution: {metadata['Name']}")
        if not (metadata['Version']):
            raise AssertionError('missing package version')
        root, = roots
        if not (root == f"cyballs-{metadata['Version']}"):
            raise AssertionError(f'inconsistent package root: {root}')
        if not (path.name == root + '.tar.gz'):
            raise AssertionError(f'inconsistent archive name: {path.name}')
    return dict(status='PASS', distribution=metadata['Name'], version=metadata['Version'],
                archive=str(path.resolve()), sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                file_count=len(names), packaged_defaults=defaults,
                native_dependencies='external GSL and CFITSIO')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('archive', type=Path)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = inspect_archive(args.archive)
    text = json.dumps(result, indent=2, sort_keys=True)+'\n'
    if args.output:
        args.output.write_text(text)
    print(text, end='')


if __name__ == '__main__':
    main()
