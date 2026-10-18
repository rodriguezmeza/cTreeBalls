"""Unwindowed scalar moments against independent Decimal tuple enumeration."""
from pathlib import Path
import sys

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from cyballs import cballs
from test_exact_numerical_robustness import EDGES, POINTS, reference, compare

ENGINES = ['octree-2balls-omp', 'balltree-2balls-omp', 'kdtree-2balls-omp']


def run(engine, positions, field, weights, directory, order=7, threads=1):
    model = cballs()
    try:
        model.set(dict(searchMethod=engine, rootDir=str(directory),
            sizeHistN=4, mChebyshev=max(2,order), nsmooth=4, theta=0,
            rminHist=EDGES[0], rangeN=EDGES[-1], useLogHist=True,
            usePeriodic=False, numberThreads=threads, verbose=0, verbose_log=0,
            options='KKKCorrelation,weights-norm,no-normalize-HistZeta,'
                    'no-smooth-pivot,no-one-ball,no-two-balls,no-out-Hist'))
        model.set_catalog(positions, kappa=field, weights=weights)
        model.Run(level=['MainLoop'])
        components = np.array([[model.getHistZetaMsincos(m+1,c).copy()
                               for c in range(1,5)] for m in range(order+1)])
        return components
    finally:
        model.clean_all()


@pytest.mark.parametrize('engine', ENGINES)
@pytest.mark.parametrize('order', [0,7])
@pytest.mark.parametrize('threads', [1,3])
@pytest.mark.parametrize('heavy', [1.,1e16])
def test_signed_unwindowed_moments(tmp_path, engine, order, threads, heavy):
    rng = np.random.default_rng(18324)
    positions = rng.normal(size=(9,3))
    positions /= np.linalg.norm(positions,axis=1)[:,None]
    field = rng.uniform(-1,1,9)
    weights = rng.uniform(.5,2,9)
    weights[0] = heavy
    weights[2] = 0
    components = run(engine,positions,field,weights,tmp_path,order,threads)
    actual = components[:,0]+components[:,1]+1j*(components[:,2]-components[:,3])
    expected,_ = reference(positions,field,weights,order,False)
    compare(actual,expected)
    # A monopole has no sine component, including when its neighbors carry
    # zero weights and signed fields or dominate the other contributions.
    np.testing.assert_array_equal(components[0,1:],np.zeros_like(components[0,1:]))


@pytest.mark.parametrize('engine', ENGINES)
def test_zero_and_representable_subnormal_field(tmp_path, engine):
    zero = run(engine,POINTS,np.zeros(4),np.ones(4),tmp_path/'zero')
    np.testing.assert_array_equal(zero,np.zeros_like(zero))
    positions = np.random.default_rng(2451).normal(size=(4,3))
    positions /= np.linalg.norm(positions,axis=1)[:,None]
    field = np.full(4,1e-103)
    components = run(engine,positions,field,np.ones(4),tmp_path/'small',order=0)
    expected,_ = reference(positions,field,np.ones(4),0,False)
    assert np.any(components[0,0] > 0)
    np.testing.assert_allclose(components[:,0],expected.real,rtol=1e-12,atol=1e-320)


@pytest.mark.parametrize('engine', ENGINES)
@pytest.mark.parametrize('field', [1e-110,1e110])
def test_unrepresentable_unwindowed_products_fail(tmp_path, engine, field):
    with pytest.raises(Exception,match='numerical range'):
        run(engine,POINTS,np.full(4,field),np.ones(4),tmp_path)
