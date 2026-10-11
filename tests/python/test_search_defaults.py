"""Fast shipped defaults, explicit overrides, and the retained exact limit."""
import os
import tempfile
from unittest.mock import patch

import numpy as np
import pytest
from cyballs import cballs
from test_shear_pivot_reuse import clustered_catalog

CONTROLS = ('CBALLS_SHEAR_PIVOT_TOL', 'CBALLS_SHEAR_BIN_THETA',
            'CBALLS_SCALAR_PIVOT_TOL', 'CBALLS_SCALAR_BIN_THETA')


@pytest.mark.parametrize('tree', ['octree', 'kdtree', 'balltree'])
@pytest.mark.parametrize('shear', [False, True])
def test_shipped_defaults_match_explicit_fast_controls(tree, shear):
    positions, gamma, weights = clustered_catalog()
    engine = tree + ('-shear-sphere-2balls-omp' if shear else '-2balls-omp')
    reuse = 'shear-pivot-reuse' if shear else 'scalar-pivot-reuse'
    def run(explicit=False, requested=False, unrelated_shear_controls=False):
        with patch.dict(os.environ), tempfile.TemporaryDirectory() as out:
            for name in CONTROLS:
                os.environ.pop(name, None)
            if explicit:
                os.environ.update(CBALLS_SHEAR_PIVOT_TOL='3', CBALLS_SHEAR_BIN_THETA='.5',
                                  CBALLS_SCALAR_PIVOT_TOL='.1', CBALLS_SCALAR_BIN_THETA='0')
            if unrelated_shear_controls:
                # Scalar engines must not even parse these shear-only values.
                os.environ.update(CBALLS_SHEAR_PIVOT_TOL='unused-by-kappa',
                                  CBALLS_SHEAR_BIN_THETA='unused-by-kappa')
            model = cballs()
            try:
                options = 'no-out-Hist,KKKCorrelation,weights-norm,no-normalize-HistZeta,edge-corrections'
                if explicit:
                    options += ',no-smooth-pivot'
                if requested:
                    options += ',' + reuse
                params = dict(searchMethod=engine, rootDir=out,
                    usePeriodic=False, useLogHist=False, rminHist=.03,
                    rangeN=1.75, sizeHistN=3, sizeHistPhi=16, mChebyshev=2,
                    numberThreads=2, verbose=0, verbose_log=0,
                    options=options)
                if explicit:
                    params['theta'] = 1
                    params['nsmooth'] = 8
                model.set(params)
                if shear:
                    model.set_catalog(positions, gamma1=gamma.real,
                                      gamma2=gamma.imag, weights=weights)
                else:
                    model.set_catalog(positions, kappa=gamma.real, weights=weights)
                model.Run(level=['MainLoop'])
                assert model.getTheta() == 1
                assert model.run_settings['effective']['nsmooth'] == 8
                if shear:
                    np.testing.assert_array_equal(model.getResults()['arrays']['xi_weight'],
                                                  model.getShearXiWeight())
                    info = model.getRunMetadata()['shear_hierarchical_reuse']
                    assert info['enabled'] == requested
                    if requested:
                        assert info['phase_budget_radians'] == 3
                        assert info['radial_bin_slop'] == .5
                    return [np.array(v, copy=True) for v in (
                        model.getShearXiPlus(), model.getShearXiMinus(),
                        model.getShearUpsilonXMultipoles(),
                        model.getShearWindowMultipoles(), model.getShearGammaXMultipoles())]
                return [np.array(v, copy=True) for v in (
                    model.getHistXi2pcf(),
                    [model.getHistZetaMsincos(i, j) for i in range(1,4) for j in range(1,5)],
                    [model.getHistZetaM_EE_complex(i) for i in range(1,4)])]
            finally:
                model.struct_cleanup()
    for requested in (False, True):
        explicit = run(explicit=True, requested=requested)
        for actual, target in zip(run(requested=requested), explicit):
            np.testing.assert_allclose(actual, target, rtol=2e-11, atol=1e-12,
                                       equal_nan=True)
        if not shear:
            for actual, target in zip(run(requested=requested,
                    unrelated_shear_controls=True), explicit):
                np.testing.assert_allclose(actual, target, rtol=2e-11,
                                           atol=1e-12, equal_nan=True)


@pytest.mark.parametrize('tree', ['octree', 'kdtree', 'balltree'])
def test_explicit_shear_overrides_enable_reuse(tree):
    positions, gamma, weights = clustered_catalog()
    with patch.dict(os.environ, CBALLS_SHEAR_PIVOT_TOL='.1',
                    CBALLS_SHEAR_BIN_THETA='0'), tempfile.TemporaryDirectory() as out:
        model = cballs()
        try:
            model.set(searchMethod=tree+'-shear-sphere-2balls-omp',
                rootDir=out, usePeriodic=False, useLogHist=False,
                rminHist=.03, rangeN=1.75, sizeHistN=3, sizeHistPhi=16,
                mChebyshev=2, nsmooth=1, numberThreads=2, theta=1,
                verbose=0, verbose_log=0, options='no-out-Hist,only-3pcf,shear-pivot-reuse')
            model.set_catalog(positions, gamma1=gamma.real, gamma2=gamma.imag,
                              weights=weights)
            model.Run(level=['MainLoop'])
            info = model.getRunMetadata()['shear_hierarchical_reuse']
            assert info['enabled']
            assert info['phase_budget_radians'] == .1
            assert info['radial_bin_slop'] == 0
        finally:
            model.struct_cleanup()


@pytest.mark.parametrize('log_bins', [False, True])
def test_core_zero_opening_cannot_accept_pruned_cells(log_bins):
    positions, gamma, weights = clustered_catalog()
    products = []
    for explicit_exact in (False, True):
        with tempfile.TemporaryDirectory() as out:
            model = cballs()
            try:
                options = 'KKKCorrelation,no-out-Hist,weights-norm,no-normalize-HistZeta'
                if explicit_exact:
                    options += ',no-one-ball,no-two-balls,no-smooth-pivot'
                model.set(searchMethod='octree-sincos-omp', rootDir=out, theta=0,
                    usePeriodic=False, useLogHist=log_bins, rminHist=.02,
                    rangeN=1.8, sizeHistN=6, mChebyshev=5, numberThreads=2,
                    verbose=0, verbose_log=0, options=options)
                model.set_catalog(positions, kappa=gamma.real, weights=weights)
                model.Run(level=['MainLoop'])
                assert model.getTheta() == 0
                products.append((model.getHistXi2pcf().copy(), np.array([
                    model.getHistZetaMsincos(i,j) for i in range(1,7)
                    for j in range(1,5)])))
            finally:
                model.struct_cleanup()
    for actual, exact in zip(*products):
        np.testing.assert_allclose(actual, exact, rtol=2e-11, atol=1e-12)
