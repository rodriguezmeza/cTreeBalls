"""Computed products, not allocated buffers, define the scalar getter contract."""
import numpy as np
import pytest
from cyballs import cballs, CosmoSevereError, CosmoComputationError

GETTERS = ('getHistNN', 'getHistCF', 'getHistXi2pcf', 'getHistXi2pcf12', 'getHistXi2pcf13')
EXACT = 'no-smooth-pivot,no-one-ball,no-two-balls,no-out-Hist'


def make_model(tmp_path, method='octree-2balls-omp', options='only-2pcf'):
    rng = np.random.default_rng(2041)
    points = rng.uniform(-.8, .8, (36, 3))
    m = cballs()
    m.set(searchMethod=method, rootDir=str(tmp_path), rangeN=1.8, rminHist=.02,
          sizeHistN=4, mChebyshev=2, numberThreads=1, theta=0,
          lengthBox=2., useLogHist=False, verbose=0, verbose_log=0,
          options=EXACT+','+options)
    if method.startswith('lya-'):
        points[:, 0] += 10
        m.set_forest_catalog(points, 1+.1*points[:, 1], np.ones(36), np.arange(36, dtype=np.int64)//3)
    else:
        extra = dict(gamma1=np.full(36, .03), gamma2=np.full(36, .02)) if 'shear' in method else {}
        m.set_catalog(points, kappa=1+.1*points[:, 0], **extra)
    return m, points


def unavailable(model, names=GETTERS):
    for name in names:
        with pytest.raises(CosmoSevereError, match='not computed|live arrays|not allocated'):
            getattr(model, name)()


@pytest.mark.parametrize('method,options', [
    ('octree-2balls-omp', 'only-2pcf'),
    ('kdtree-2balls-omp', 'only-2pcf'),
    ('balltree-2balls-omp', 'KKKCorrelation'),
    ('octree-2balls-omp', 'only-3pcf'),
    ('octree-shear-sphere-2balls-omp', 'only-2pcf'),
    ('octree-3pcf-3d-omp', 'only-2pcf'),
    ('lya-2pcf-omp', 'only-2pcf'),
])
def test_uncomputed_products_reject_even_after_success(tmp_path, method, options):
    m, _ = make_model(tmp_path, method, options)
    try:
        unavailable(m)
        m.Run()
        unavailable(m, ('getHistCF', 'getHistXi2pcf12', 'getHistXi2pcf13'))
        if 'only-3pcf' in options or any(x in method for x in ('shear', '3pcf-3d', 'lya-')):
            unavailable(m)
        else:
            assert np.isfinite(m.getHistXi2pcf()).all()
        m.struct_cleanup()
        unavailable(m)
    finally:
        m.clean_all()


@pytest.mark.parametrize('method', ['octree-2balls-omp', 'kdtree-2balls-omp', 'balltree-2balls-omp', 'octree-sincos-omp'])
def test_count_cf_matches_independent_shell_normalization(tmp_path, method):
    mode = 'KKKCorrelation' if method == 'octree-sincos-omp' else 'only-2pcf'
    m, points = make_model(tmp_path, method, mode+',compute-HistN,and-CF')
    try:
        m.Run()
        edges = np.linspace(.02, 1.8, 5)
        d = np.linalg.norm(points[:, None]-points[None, :], axis=-1)
        distances = d[np.triu_indices(len(points), 1)]
        dd = np.histogram(distances[(distances > edges[0]) & (distances < edges[-1])], bins=edges)[0]
        expected = 2*dd*2.**3/(len(points)**2*(4*np.pi/3)*np.diff(edges**3))-1
        np.testing.assert_array_equal(m.getHistNN(), dd)
        np.testing.assert_allclose(m.getHistCF(), expected, rtol=2e-13, atol=2e-13)
        retained = m.getHistCF()
        m.set(options=EXACT+','+mode)
        unavailable(m)
        m.Run()
        unavailable(m, ('getHistCF',))
        np.testing.assert_allclose(retained, expected)
        m.set(rangeN=-1.)
        with pytest.raises((CosmoSevereError, CosmoComputationError)):
            m.Run()
        unavailable(m)
        m.set(rangeN=1.8, options=EXACT+','+mode+',compute-HistN,and-CF')
        m.Run()
        np.testing.assert_allclose(m.getHistCF(), expected, rtol=2e-13, atol=2e-13)
    finally:
        m.clean_all()


def test_startup_or_tree_only_is_not_a_completed_product(tmp_path):
    m, _ = make_model(tmp_path)
    try:
        m.Run(level=['SetNumberThreads'])
        unavailable(m)
        m.set(options=EXACT+',only-2pcf,make-tree')
        m.Run()
        unavailable(m)
    finally:
        m.clean_all()


def test_neighbor_box_empty_bins_are_initialized_products(tmp_path):
    m, _ = make_model(tmp_path, 'neighbor-boxes-omp')
    try:
        m.set(usePeriodic=True, lengthBox=10., rangeN=.4, rminHist=0.)
        m.set_catalog(np.array([[0.,0.,0.], [2.,0.,0.], [0.,2.,0.]]))
        m.Run()
        np.testing.assert_array_equal(m.getHistNN(), np.zeros(4))
        # Preserve this engine's documented file-output convention for empty bins.
        np.testing.assert_array_equal(m.getHistCF(), np.zeros(4))
        unavailable(m, ('getHistXi2pcf',))
    finally:
        m.clean_all()
