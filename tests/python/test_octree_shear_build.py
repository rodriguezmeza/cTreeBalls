"""Native octree construction, parallel spherical moments and failure cleanup."""
import os
import tempfile
from unittest.mock import patch
import numpy as np
from cyballs import cballs, CosmoComputationError
import test_shear_sphere_octree_omp as ref

ENGINE = 'octree-shear-sphere-2balls-omp'
OPTIONS = 'no-out-Hist,no-smooth-pivot'


def compare_threads(catalogs, **kwargs):
    with patch.dict(os.environ, CBALLS_SHEAR_TASKS='128',
                    CBALLS_SHEAR_PIVOT_TOL='3', CBALLS_SHEAR_BIN_THETA='.5'):
        a = ref.run_native_catalogs(catalogs, 1, engine=ENGINE, **kwargs)
        b = ref.run_native_catalogs(catalogs, 4, engine=ENGINE, **kwargs)
    for key in a:
        np.testing.assert_array_equal(b[key], a[key], err_msg=key)
        assert np.isfinite(b[key]).all(), key


def test_parallel_octants_masks_and_exact():
    catalog = ref.fixture(32771)
    # Axis and split-plane points, including both poles, and zero weights.
    catalog[0][:6] = np.vstack((np.eye(3), -np.eye(3)))
    catalog[2][::17] = 0.0
    mask = (np.arange(len(catalog[0])) % 11 != 0).astype(np.int32)
    with patch.object(ref, 'RMIN', .003), patch.object(ref, 'RMAX', .04):
        for option in ['', ',no-one-ball', ',shear-pivot-reuse']:
            for masked in [False, True]:
                compare_threads([catalog], nsmooth=1, use_log=True,
                    masks=[mask] if masked else None,
                    options=OPTIONS+option+(',read-mask' if masked else ''))


def test_parallel_unbalanced_cross_catalogs():
    catalog = ref.fixture(32773)
    catalog[0][:] += [4., .5, .2]
    catalog[0][:] /= np.linalg.norm(catalog[0], axis=1)[:, None]
    other = ref.fixture(20003)
    catalogs = [catalog, other, ref.fixture(35)]
    masks = [(np.arange(len(c[0])) % 13 != 0).astype(np.int32) for c in catalogs]
    with patch.object(ref, 'RMIN', .003), patch.object(ref, 'RMAX', .015):
        compare_threads(catalogs, nsmooth=1, masks=masks, use_log=False,
            catalog_order='1,2,1', options=OPTIONS+',read-mask,shear-pivot-reuse')


def make_model(directory, catalog):
    model = cballs()
    model.set(dict(searchMethod=ENGINE, iCatalogs='1', usePeriodic=False,
        useLogHist=True, rminHist=.003, rangeN=.02, sizeHistN=3,
        sizeHistPhi=16, mChebyshev=2, theta=1, nsmooth=1, numberThreads=4,
        lengthBox=2.2, verbose=0, verbose_log=0, rootDir=directory,
        options=OPTIONS+',only-3pcf,shear-pivot-reuse'))
    p,g,w = catalog
    model.set_catalog(p, gamma1=g.real, gamma2=g.imag, weights=w)
    return model


def test_parallel_duplicate_failure_and_recovery():
    catalog = ref.fixture(20003)
    good = catalog[0].copy()
    catalog[0][1] = catalog[0][0]
    with tempfile.TemporaryDirectory(prefix='octree-build-duplicate-') as directory:
        model = make_model(directory, catalog)
        try:
            try: model.Run(level=['MainLoop'])
            except CosmoComputationError as exc:
                assert 'same position' in str(exc), str(exc)
            else: raise AssertionError('duplicate coordinates were accepted')
            assert not model.getTreeAllocated()
            assert not model.getBodytableAllocated()
            model.set_catalog(good, gamma1=catalog[1].real, gamma2=catalog[1].imag,
                              weights=catalog[2])
            model.Run(level=['MainLoop'])
            assert np.isfinite(model.getShearUpsilonXMultipoles()).all()
        finally: model.struct_cleanup()


def test_parallel_allocation_failure_and_recovery():
    catalog = ref.fixture(20003)
    for allocations in (0, 2, 32, 512):
        with tempfile.TemporaryDirectory(prefix='octree-build-allocation-') as directory:
            model = make_model(directory, catalog)
            try:
                model.Run(level=['StartRun_Common'])
                model._set_allocation_failure_after_for_tests(allocations)
                try:
                    try: model.Run(level=['MainLoop'])
                    except CosmoComputationError as exc:
                        assert 'allocation' in str(exc), str(exc)
                    else: raise AssertionError('injected allocation failure was ignored')
                finally: model._reset_allocation_failure_for_tests()
                assert not model.getTreeAllocated()
                assert not model.getBodytableAllocated()
                p,g,w=catalog
                model.set_catalog(p, gamma1=g.real, gamma2=g.imag, weights=w)
                model.Run(level=['MainLoop'])
                assert np.isfinite(model.getShearUpsilonXMultipoles()).all()
            finally: model.struct_cleanup()


if __name__ == '__main__':
    for test in [test_parallel_octants_masks_and_exact,
                 test_parallel_unbalanced_cross_catalogs,
                 test_parallel_duplicate_failure_and_recovery,
                 test_parallel_allocation_failure_and_recovery]:
        test()
        print('PASS', test.__name__, flush=True)
