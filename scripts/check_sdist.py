#!/usr/bin/env python3
"""Check the public cyballs source archive, without extracting untrusted paths."""
import argparse
from email.parser import BytesParser
import hashlib
import json
from pathlib import Path, PurePosixPath
import tarfile

REQUIRED = {
    'capabilities/engines.json','include/engine_registry_generated.h','include/mpi_runtime.h','include/tree_workspace.h',
    'source/engine_registry.c','source/runtime_context.c','source/memory_catalog.c',
    'source/common_histogram.c','source/smooth_pivots.c','source/mpi_runtime.c',
    'scripts/generate_capabilities.py','scripts/capabilities_generated.py',
    'scripts/affected_regressions.py','scripts/benchmark_contracts.py',
    'tests/test_runtime_context.c','tests/test_mpi_runtime.c',
    'tests/make_tests/test_capability_contracts.py','ENGINE_CAPABILITIES.md','MAINTAINABILITY_CONTRACTS.md',

    'include/resource_contracts.h', 'scripts/accuracy_acceptance.py',
    'tests/test_resource_contracts.c', 'tests/make_tests/test_resource_scientific_contracts.py',
    'RESOURCE_SCIENTIFIC_CONTRACTS.md', 'PKG-INFO', 'setup.py', 'pyproject.toml', 'MANIFEST.in', 'Makefile',
    'Makefile_machine', 'Makefile_settings', 'requirements/build.txt',
    'requirements/release.txt', 'include/tree_contracts.h',
    'python/cyballs.pyx', 'python/ccyballs.pxd.in',
    'scripts/build_fingerprint.py', 'scripts/active_release_gate.py',
    'scripts/release_gate_case.py', 'scripts/check_installed_package.py',
    'addons/Makefile_gsl', 'addons/cfitsio/Makefile_cfitsio',
    'addons/lya_forest_omp/search_lya_forest_los_tree_omp.c',
    'addons/lya_forest_omp/lya_forest_los_tree.h',
    'addons/lya_forest_omp/Makefile_lya_forest_omp',
    'addons/octree_3pcf_3d_omp/search_octree_3pcf_3d_omp.c',
    'tests/make_tests/test_lya_forest_los_tree.py',
    'tests/make_tests/test_cython_in_memory_catalog.py',
    'tests/make_tests/test_release_packaging.py',
}

ALLOWED_ADDONS = {
    'kdtree_2balls_omp', 'kdtree_2balls_mpi', 'balltree_2balls_omp', 'balltree_2balls_mpi',
    'octree_2balls_omp', 'octree_2balls_mpi', 'lya_forest_omp', 'lya_forest_mpi',
    'octree_shear_sphere_2balls_omp', 'kdtree_shear_sphere_2balls_omp',
    'balltree_shear_sphere_2balls_omp', 'kdtree_box_omp', 'neighbor_boxes_omp',
    'octree_3pcf_3d_omp', 'octree_3pcf_3d_mpi', 'gadget_io', 'class_lib', 'iolib',
    'cfitsio', 'pxd', 'cmdline_defs_settings', 'addons_include', 'balltree_shared',
    'kdtree_shared', 'scalar_compat_shared', 'shear_sphere_shared',
    'shear_sphere_binary_2balls', 'native_octree_pair',
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
        assert len(roots) == 1, f'expected one source root: {roots}'
        assert REQUIRED <= names, f'missing public source files: {sorted(REQUIRED-names)}'
        forbidden = [name for name in names
                     if (name.startswith('addons/') and len(PurePosixPath(name).parts) > 2
                         and PurePosixPath(name).parts[1] not in ALLOWED_ADDONS)
                     or name.endswith('_TC.py') or name.startswith('run2/')
                     or PurePosixPath(name).suffix in {'.o', '.a', '.so', '.dylib', '.pyc'}]
        assert not forbidden, f'inactive/private addons or build artifacts: {forbidden}'
        metadata = BytesParser().parsebytes(archive.extractfile(files['PKG-INFO']).read())
        assert metadata['Name'] == 'cyballs', f'unexpected distribution: {metadata["Name"]}'
        assert metadata['Version'], 'missing package version'
        root, = roots
        assert root == f'cyballs-{metadata["Version"]}', f'inconsistent package root: {root}'
        assert path.name == root+'.tar.gz', f'inconsistent archive name: {path.name}'
    return dict(status='PASS', distribution=metadata['Name'], version=metadata['Version'],
                archive=str(path.resolve()), sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                file_count=len(names), native_dependencies='external GSL and CFITSIO')


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
