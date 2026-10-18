from libc.stdint cimport int64_t, uint64_t
import operator as _dm_operator

cdef extern from "lya_distortion.h":
    ctypedef struct cballs_lya_distortion_config:
        int rp_bins, rt_bins, model_rp_bins, model_rt_bins
        int projection_order, threads
        double rp_max, rt_max, model_rp_max, model_rt_max
        size_t max_bytes
    int cballs_lya_distortion_available()
    int cballs_lya_distortion_plan(size_t, const cballs_lya_distortion_config*, size_t*, char*, size_t)
    int cballs_lya_distortion_compute(const double*, const double*, const int64_t*,
        const double*, const double*, size_t, const cballs_lya_distortion_config*,
        double*, size_t, double*, size_t, uint64_t*, char*, size_t)


def _dm_array(value, name):
    array = np.asarray(value)
    if array.dtype.kind not in 'biuf':
        raise CosmoSevereError(name + ' must contain real numbers')
    array = np.ascontiguousarray(array, dtype=np.float64)
    if not np.isfinite(array).all():
        raise CosmoSevereError(name + ' must contain finite values')
    return array


def lya_distortion_matrix(positions, weights, forest_ids, log_lambda=None, *,
                          projection_weights=None, projection_order=1,
                          rp_max=200., rt_max=200., rp_bins=50, rt_bins=50,
                          model_rp_max=None, model_rt_max=None,
                          model_rp_bins=None, model_rt_bins=None,
                          number_threads=1, max_memory_mb=512.):
    """Native continuum-projection 2PCF response, xi_measured = D @ xi_model.

    positions are observer-centred (N,3); IDs are exact int64 forest labels.
    weights select/weight measured pixels. projection_weights (default weights)
    define the fit over the whole supplied forest, including pixels with zero
    data weight. Supply their true log wavelengths for mean+slope projection;
    a different logarithm base or affine shift gives the same projector.

    projection_order=1 removes mean+slope, 0 removes mean, -1 is an identity
    diagnostic. Constant-coordinate/singleton forests use mean-only projection.
    Data/model grids may differ; correlations outside MODEL support are zero
    in this forward model. Rows with no data weight are undefined (all NaN).
    Numerator and WDM are returned separately. No delta values enter D.

    Fixed primary-forest tasks commit in forest order, independently of thread
    assignment. Calls retain the GIL and keep no pointers after return.
    max_memory_mb bounds planned native workspace, input channels and outputs;
    CBALLS_MEMORY_BUDGET_MB also applies. MPI callers should compute once on
    the publishing rank; this routine parallelises with OpenMP, not MPI.
    """
    cdef cballs_lya_distortion_config c
    cdef char error[1024]
    cdef size_t planned = 0
    cdef uint64_t counts[8]
    cdef const double *lambda_ptr = NULL
    cdef const double *projection_ptr = NULL
    cdef np.ndarray pos, weight, ids, lam, pw, numerator, row_weights
    if not cballs_lya_distortion_available():
        raise CosmoSevereError('distortion computation requires a strict double 3D Ly-alpha build')
    pos = _dm_array(positions, 'positions')
    weight = _dm_array(weights, 'weights')
    labels = np.asarray(forest_ids)
    if pos.ndim != 2 or pos.shape[1] != 3 or pos.shape[0] < 1:
        raise CosmoSevereError('positions must have shape (N,3), N >= 1')
    n = pos.shape[0]
    if weight.ndim != 1 or len(weight) != n or np.any(weight < 0):
        raise CosmoSevereError('weights must be a non-negative (N,) array')
    if labels.ndim != 1 or len(labels) != n or labels.dtype.kind not in 'iu':
        raise CosmoSevereError('forest_ids must be an integer (N,) array')
    if labels.dtype.kind == 'u' and np.any(labels > np.iinfo(np.int64).max):
        raise CosmoSevereError('forest_ids exceed signed int64')
    ids = np.ascontiguousarray(labels, dtype=np.int64)
    if log_lambda is not None:
        lam = _dm_array(log_lambda, 'log_lambda')
        if lam.ndim != 1 or len(lam) != n:
            raise CosmoSevereError('log_lambda must have shape (N,)')
        lambda_ptr = <const double*>np.PyArray_DATA(lam)
    if projection_weights is not None:
        pw = _dm_array(projection_weights, 'projection_weights')
        if pw.ndim != 1 or len(pw) != n or np.any(pw < 0):
            raise CosmoSevereError('projection_weights must be a non-negative (N,) array')
        projection_ptr = <const double*>np.PyArray_DATA(pw)
    c.rp_bins = _dm_operator.index(rp_bins)
    c.rt_bins = _dm_operator.index(rt_bins)
    c.model_rp_bins = c.rp_bins if model_rp_bins is None else _dm_operator.index(model_rp_bins)
    c.model_rt_bins = c.rt_bins if model_rt_bins is None else _dm_operator.index(model_rt_bins)
    c.projection_order = _dm_operator.index(projection_order)
    c.threads = _dm_operator.index(number_threads)
    c.rp_max = rp_max; c.rt_max = rt_max
    c.model_rp_max = rp_max if model_rp_max is None else model_rp_max
    c.model_rt_max = rt_max if model_rt_max is None else model_rt_max
    budget = float(max_memory_mb)*2**20
    if not np.isfinite(budget) or budget <= 0 or budget > np.iinfo(np.intp).max:
        raise CosmoSevereError('max_memory_mb must be positive, finite and representable')
    c.max_bytes = int(budget)
    if c.projection_order == 1 and lambda_ptr == NULL:
        raise CosmoSevereError('mean+slope distortion requires true log_lambda metadata; it cannot be inferred from distance')
    if cballs_lya_distortion_plan(n, &c, &planned, error, sizeof(error)):
        raise CosmoSevereError(error)
    na, nb = int(c.rp_bins)*c.rt_bins, int(c.model_rp_bins)*c.model_rt_bins
    numerator = np.zeros((na, nb), dtype=np.float64)
    row_weights = np.zeros(na, dtype=np.float64)
    if cballs_lya_distortion_compute(<const double*>np.PyArray_DATA(pos),
        <const double*>np.PyArray_DATA(weight), <const int64_t*>np.PyArray_DATA(ids),
        lambda_ptr, projection_ptr, n, &c, <double*>np.PyArray_DATA(numerator),
        na*nb, <double*>np.PyArray_DATA(row_weights), na, counts, error, sizeof(error)):
        raise CosmoComputationError(error)
    valid = row_weights > 0
    matrix = np.full((na, nb), np.nan, dtype=np.float64)
    with np.errstate(over='ignore', invalid='ignore', divide='ignore'):
        np.divide(numerator, row_weights[:, None], out=matrix, where=valid[:, None])
    if np.any(valid & ~np.isfinite(matrix).all(axis=1)):
        raise CosmoComputationError('distortion row normalisation overflowed')
    names = ('pixels', 'forests', 'forest_pairs', 'pruned_forest_pairs',
             'forest_pairs_with_data', 'positive_weight_data_pairs',
             'degenerate_slope_forests', 'planned_bytes')
    return dict(matrix=matrix, numerator=numerator, weights=row_weights, valid_rows=valid,
        rp_max=c.rp_max, rt_max=c.rt_max, rp_bins=c.rp_bins, rt_bins=c.rt_bins,
        model_rp_max=c.model_rp_max, model_rt_max=c.model_rt_max,
        model_rp_bins=c.model_rp_bins, model_rt_bins=c.model_rt_bins,
        projection_order=c.projection_order, stats={name: int(counts[i]) for i,name in enumerate(names)},
        rp=np.repeat((np.arange(c.rp_bins)+.5)*c.rp_max/c.rp_bins,c.rt_bins),
        rt=np.tile((np.arange(c.rt_bins)+.5)*c.rt_max/c.rt_bins,c.rp_bins),
        model_rp=np.repeat((np.arange(c.model_rp_bins)+.5)*c.model_rp_max/c.model_rp_bins,c.model_rt_bins),
        model_rt=np.tile((np.arange(c.model_rt_bins)+.5)*c.model_rt_max/c.model_rt_bins,c.model_rp_bins))
