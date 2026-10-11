#!/usr/bin/env python3
"""MPI reuse parity and collective recovery. Launch with mpiexec -n 2 or -n 4.

Use an external timeout: a collective regression must fail rather than hang CI.
"""

import sys
from pathlib import Path

# Locate sibling analysis and benchmark helpers for direct script execution.
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python"))

import os
import ctypes
import re
import tempfile
import traceback
from unittest.mock import patch

sys.path[:0] = [str(Path(__file__).resolve().parents[2]), str(Path(__file__).resolve().parent)]
import numpy as np
from mpi4py import MPI
from cyballs import cballs, search_method_id
import test_shear_sphere_octree_omp as ref
from test_shear_pivot_reuse import clustered_catalog

COMM = MPI.COMM_WORLD
BASE = 'no-out-Hist,no-smooth-pivot,only-3pcf,shear-pivot-reuse'


def compare(actual, expected):
    for key in actual:
        if key != 'multipoles':
            np.testing.assert_allclose(actual[key], expected[key], rtol=3e-11, atol=3e-12, err_msg=key)
            continue
        if actual['window'].shape[0] == 1:
            np.testing.assert_allclose(actual[key], expected[key], rtol=3e-11, atol=3e-12)
            continue
        modes = np.arange(2*ref.NMAX+1)
        for i in range(ref.BINS):
            for j in range(ref.BINS):
                matrices = [v['window'][modes[:,None]-modes[None,:]+2*ref.NMAX,i,j]
                            for v in (actual, expected)]
                if max(np.linalg.cond(m) for m in matrices) < 1e8:
                    np.testing.assert_allclose(actual[key][:,:,i,j], expected[key][:,:,i,j], rtol=5e-9, atol=5e-12)
                else:
                    for result, matrix in zip((actual, expected), matrices):
                        x=result[key][:,:,i,j].T; b=result['upsilon'][:,:,i,j].T
                        if not np.any(x): continue  # Existing singular-window sentinel.
                        assert np.linalg.norm(matrix@x-b) <= 1e-11*(np.linalg.norm(matrix)*np.linalg.norm(x)+np.linalg.norm(b))+1e-24


def check_root(callback):
    error = None
    if COMM.rank == 0:
        try: callback()
        except Exception: error = traceback.format_exc()
    errors = COMM.allgather(error)
    if any(errors): raise AssertionError('\n'.join(e for e in errors if e))


def numerical(tree):
    engine = tree+'-shear-sphere-2balls-'
    first = clustered_catalog()
    second, third = ref.fixture(79), ref.fixture(63)
    masks = [np.ones(len(c[0]), dtype=np.uint8) for c in (first, second, third)]
    for mask in masks: mask[::5] = 0
    cases = [
        ('exact-limit', [ref.fixture(80)], {}, '1e-12', '0', 64, BASE, True),
        ('clustered', [first], {}, '.1', '0', 8, BASE, False),
        ('bin-slop', [first], {}, '3', '.2', 8, BASE, False),
        ('masked-cross', [first, second, third], dict(masks=masks), '.1', '.1', 8, BASE+',read-mask', False),
        ('repeated-role', [first, second, third], dict(masks=masks, catalog_order='1,2,1'), '.1', '0', 8, BASE+',read-mask', False),
        ('n0', [first], {}, '.1', '.2', 8, BASE+',no-count-multipoles', False),
        ('forced-window', [first], {}, '.1', '0', 8, BASE+',no-count-multipoles,edge-corrections', False),
        ('idle-ranks', [ref.fixture(3)], dict(nsmooth=1), '.1', '0', 1, BASE, False),
        ('zero-budget', [second], {}, '0', '0', 8, BASE, False),
        ('exact-fallback', [second], {}, '3', '.2', 8, BASE+',no-one-ball', False),
        ('combined', [first], {}, '.1', '0', 8, BASE.replace(',only-3pcf',''), False),
        ('pair-only', [first], {}, '.1', '0', 8, BASE.replace('only-3pcf','only-2pcf'), False),
        ('smooth-fallback', [ref.smooth_fixture()], dict(rsmooth=40), '.1', '.2', 8, BASE.replace('no-smooth-pivot','smooth-pivot'), False),
        ('expanded-frontier', [first], dict(nsmooth=1), '.1', '0', 320, BASE, False),
        ('memory-limited', [first], {}, '.1', '0', 4096, BASE, False),
    ]
    for name, catalogs, kwargs, phase, bin_slop, tasks, options, oracle in cases:
        env=dict(CBALLS_SHEAR_PIVOT_TOL=phase, CBALLS_SHEAR_BIN_THETA=bin_slop,
                 CBALLS_SHEAR_TASKS=str(tasks), CBALLS_SHEAR_TASK_MEMORY_MB='1' if name=='memory-limited' else '256')
        with patch.dict(os.environ, env):
            profiled = name in ('clustered', 'expanded-frontier', 'memory-limited')
            with tempfile.TemporaryFile(mode='w+b') as log:
                saved = os.dup(2)
                try:
                    if profiled: os.dup2(log.fileno(), 2)
                    with patch.dict(os.environ, {'CBALLS_SHEAR_PROFILE':'1'} if profiled else {}):
                        actual=ref.run_native_catalogs(catalogs, 2, engine=engine+'mpi', options=options,
                                                      collect_results=COMM.rank==0, **kwargs)
                    ctypes.CDLL(None).fflush(None)
                finally:
                    os.dup2(saved, 2); os.close(saved)
                log.seek(0); profile=log.read().decode()
            if profiled:
                rows=re.findall(r'frontier_tasks=(\d+) local_tasks=(\d+) histogram_bytes=(\d+)', profile)
                assert len(rows)==1, profile
                total, local, histbytes=map(int,rows[0])
                assert COMM.allreduce(local)==total
                if name=='expanded-frontier': assert total>256, profile
                if name=='memory-limited': assert histbytes<=1048576, profile
                if name=='clustered':
                    aggregates=sum(map(int,re.findall(r'aggregate_pivots=(\d+)',profile)))
                    assert COMM.allreduce(aggregates)>0, profile
            def verify():
                if name=='memory-limited':
                    # Rank budget changes the MPI partition; the exact-limit
                    # case separately establishes the physical estimator.
                    assert np.isfinite(actual['upsilon']).all()
                    return
                expected=ref.run_native_catalogs(catalogs, 2, engine=engine+'omp', options=options, **kwargs)
                compare(actual, expected)
                if oracle:
                    independent=ref.oracle(*catalogs[0])
                    compare(actual, {key: independent[key] for key in actual})
            check_root(verify)
        if COMM.rank==0: print('MPI_REUSE_NUMERICAL_OK', tree, name, flush=True)


def parallel_tree_build(tree="kdtree"):
    """Exercise rank-local build tasks above their population threshold."""
    catalog = ref.fixture(32771)
    catalog[2][::17] = 0.0
    mask = (np.arange(len(catalog[0])) % 11 != 0).astype(np.int32)
    assert mask.sum() >= 16384
    env = dict(CBALLS_SHEAR_PIVOT_TOL='3', CBALLS_SHEAR_BIN_THETA='.5',
               CBALLS_SHEAR_TASKS='128', CBALLS_SHEAR_TASK_MEMORY_MB='256')
    with patch.dict(os.environ, env), patch.object(ref, 'RMIN', .003), \
         patch.object(ref, 'RMAX', .04):
        for capacity in ((1, 8, 31) if tree in ('kdtree', 'balltree') else (1,)):
            kwargs = dict(options=BASE+',read-mask', masks=[mask],
                          nsmooth=capacity, use_log=True)
            actual = ref.run_native_catalogs([catalog], 2,
                engine=tree+'-shear-sphere-2balls-mpi',
                collect_results=COMM.rank == 0, **kwargs)
            def verify():
                expected = ref.run_native_catalogs([catalog], 1,
                    engine=tree+'-shear-sphere-2balls-omp', **kwargs)
                compare(actual, expected)
            check_root(verify)
    if COMM.rank == 0: print('MPI_PARALLEL_BUILD_OK', tree, flush=True)


def failures(tree):
    cases = [
        ('CBALLS_SHEAR_PIVOT_TOL', 'bad', 'MPI shear pivot-reuse controls'),
        ('CBALLS_SHEAR_BIN_THETA', 'nan', 'MPI shear pivot-reuse controls'),
        ('CBALLS_SHEAR_PIVOT_TOL', '.3', 'inconsistent run controls'),
        ('CBALLS_SHEAR_BIN_THETA', '.4', 'inconsistent run controls'),
        ('CBALLS_SHEAR_TASKS', 'bad', 'CBALLS_SHEAR_TASKS'),
        ('CBALLS_SHEAR_TASKS', '16', 'inconsistent run controls'),
        ('CBALLS_SHEAR_TASK_MEMORY_MB', '0', 'CBALLS_SHEAR_TASK_MEMORY_MB'),
        ('CBALLS_SHEAR_TASK_MEMORY_MB', '1', 'inconsistent run controls'),
    ]
    # Root and non-root application failures must stop before histogram MPI,
    # clean partial storage, preserve the same message, and permit another run.
    boundaries = ['MPI shear pivot-reuse controls', 'MPI shear frontier controls',
                  'MPI shear traversal workspace', 'MPI shear pivot-reuse traversal']
    for boundary in boundaries:
        for failing_rank in (0, COMM.size-1): cases.append(('inject', (boundary,failing_rank), boundary))
    with tempfile.TemporaryDirectory(prefix='mpi-reuse-fail-') as directory:
        for control, value, message in cases:
            model=cballs()
            model.set(dict(searchMethod=tree+'-shear-sphere-2balls-mpi',options=BASE,
                rootDir=directory,numberThreads=2,iCatalogs='1',usePeriodic=False,useLogHist=False,
                rangeN=ref.RMAX,rminHist=ref.RMIN,sizeHistN=ref.BINS,sizeHistPhi=ref.PHI_BINS,
                mChebyshev=ref.NMAX,theta=1,nsmooth=1,lengthBox=2.2,verbose=0,verbose_log=0))
            p,g,w=ref.fixture(48)
            def catalog(): model.set_catalog(p,gamma1=g.real,gamma2=g.imag,weights=w)
            env=dict(CBALLS_SHEAR_PIVOT_TOL='.1', CBALLS_SHEAR_BIN_THETA='0',
                     CBALLS_SHEAR_TASKS='8', CBALLS_SHEAR_TASK_MEMORY_MB='256')
            catalog()
            try:
                with patch.dict(os.environ, env):
                    with patch.dict(os.environ, {control:str(value)} if control!='inject' and COMM.rank==COMM.size-1 else {}):
                        if control=='inject':model._mpi_test_boundary(value[0] if COMM.rank==value[1] else '')
                        error=None
                        try: model.Run(level=['MainLoop'])
                        except Exception as exc: error=str(exc)
                    errors=COMM.allgather(error)
                    assert all(errors), errors
                    assert len(set(errors))==1, errors
                    # Parser failures carry their precise environment name.
                    assert message in errors[0] or (control in errors[0] and value in ('bad','nan')), errors
                    assert not model.getBodytableAllocated()
                    model._mpi_test_boundary()
                    catalog(); model.Run(level=['MainLoop'])
                    metadata = model.getRunMetadata()['shear_frontier']
                    assert metadata['target']==8 and metadata['global_tasks']>=8
                    assert metadata['local_tasks']>=0 and metadata['histogram_bytes']>0
                    assert COMM.allreduce(metadata['local_tasks'])==metadata['global_tasks']
                    with patch.dict(os.environ, CBALLS_SHEAR_TASKS='999'):
                        assert model.getRunMetadata()['shear_frontier']==metadata
                    if COMM.rank==0: assert np.isfinite(model.getShearUpsilonXMultipoles()).all()
            finally: model.struct_cleanup()
            if COMM.rank==0: print('MPI_REUSE_RECOVERY_OK', tree, control, value, flush=True)


def driver_publication(tree):
    import shear_corr_all_engines as driver
    p,g,w=ref.fixture(64)
    catalog=driver.ShearCatalog(positions=p,gamma1=g.real,gamma2=g.imag,
                               weights=w,geometry='sphere',metadata={})
    with tempfile.TemporaryDirectory(prefix='mpi-driver-rank-') as directory:
        config=driver.RunConfig(statistics='3pcf', min_sep=2, max_sep=100,
            sep_units='degree', bins=3,multipoles=2,phi_bins=16,threads=2,
            mpi_ranks=COMM.size,output_dir=Path(directory),options=tuple(BASE.split(',')))
        # An intentionally wrong hint must not decide publication or filenames.
        result=driver.run_ctreeballs(catalog,config,tree+'-shear-sphere-2balls-mpi',_mpi_rank=999)
        assert result['provenance']['parallel']['rank']==COMM.rank
        assert ('upsilon' in result)==(COMM.rank==0)
    COMM.Barrier()
    if COMM.rank==0:print('MPI_DRIVER_GLOBAL_RANK_OK',tree,flush=True)


def main():
    assert COMM.size >= 2, 'launch with at least two MPI ranks'
    for tree in ('octree','kdtree','balltree'):
        if search_method_id(tree+'-shear-sphere-2balls-mpi') < 0: continue
        numerical(tree)
        parallel_tree_build(tree)
        failures(tree)
        driver_publication(tree)
    if COMM.rank==0: print('MPI_REUSE_OK', flush=True)


if __name__=='__main__': main()
