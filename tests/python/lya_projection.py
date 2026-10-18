"""Weighted per-forest continuum projection, matching lya2pcf fill_dw.

This pure Python source is also included in the cyballs extension. Engines
receive delta_projected, and multiply by the original weights themselves.
"""
import numpy as np


def _lya_projection_arrays(delta, weights, log_lambda):
    arrays = []
    for name, value in (('delta', delta), ('weights', weights), ('log_lambda', log_lambda)):
        a = np.asarray(value)
        if a.dtype.kind not in 'biuf' or a.ndim != 1:
            raise ValueError(name+' must be a one-dimensional real array')
        a = np.asarray(a, dtype=np.float64)
        if not np.isfinite(a).all():
            raise ValueError(name+' must contain finite values')
        arrays.append(a)
    if len({len(a) for a in arrays}) != 1:
        raise ValueError('delta, weights and log_lambda must have the same length')
    if np.any(arrays[1] < 0):
        raise ValueError('projection weights must be non-negative')
    return arrays


def _lya_project_one(delta, weights, lam, apply):
    positive = weights > 0
    if not positive.any():
        return delta.copy() if apply else None, None
    if not positive.all():
        active, stats = _lya_project_one(delta[positive], weights[positive], lam[positive], apply)
        projected = delta.copy() if apply else None
        if apply:
            projected[positive] = active
        return projected, stats
    # Scaling leaves the projector invariant and avoids overflow in weighted
    # sums. Translate wavelength before centering to preserve small spans.
    w = weights/weights.max()
    total = w.sum()
    with np.errstate(over='ignore', invalid='ignore'):
        x = lam-lam[np.flatnonzero(positive)[0]]
        x -= np.dot(w, x)/total
    if not np.isfinite(x).all():
        raise ValueError('log_lambda span is not representable')
    scale_x = np.max(np.abs(x[positive]), initial=0.)
    if scale_x:
        x = x/scale_x
    variance = np.dot(w, x*x)
    scale_delta = max(float(np.max(np.abs(delta), initial=0.)), np.finfo(float).tiny)
    value = delta/scale_delta
    mean = np.dot(w, value)/total
    moment = np.dot(w*x, value)
    result = value.copy() if apply else value
    if apply:
        result -= mean
        if variance > 0:
            result -= x*np.dot(w*x, result)/variance
    normalized_moment = np.dot(w*x, result)/np.sqrt(total*variance) if variance > 0 else 0.
    stats = dict(mean_only=not bool(variance > 0),
        mean_before=float(abs(mean)*scale_delta),
        normalized_first_moment_before=float(abs(moment)/np.sqrt(total*variance)*scale_delta) if variance > 0 else 0.,
        mean_after=float(abs(np.dot(w, result)/total)*scale_delta),
        normalized_first_moment_after=float(abs(normalized_moment)*scale_delta))
    with np.errstate(over='ignore', invalid='ignore'):
        projected = result*scale_delta if apply else None
    if apply and not np.isfinite(projected).all():
        raise ValueError('projected delta is not representable')
    return projected, stats


def project_forest_delta(delta, weights, log_lambda):
    """Return a new projected delta for one forest; leave weights unchanged.

    Remove the weighted mean and weighted first moment in true log wavelength.
    A constant-wavelength or one-positive-pixel forest uses mean removal only.
    Zero-weight pixels are unchanged because they cannot enter an estimator.
    """
    delta, weights, lam = _lya_projection_arrays(delta, weights, log_lambda)
    return _lya_project_one(delta, weights, lam, True)[0]


def _lya_project_groups(delta, weights, forest_ids, log_lambda, apply):
    delta, weights, lam = _lya_projection_arrays(delta, weights, log_lambda)
    ids = np.asarray(forest_ids)
    if ids.ndim != 1 or len(ids) != len(delta) or ids.dtype.kind not in 'iu':
        raise ValueError('forest_ids must be a same-length integer array')
    if ids.dtype.kind == 'u' and np.any(ids > np.iinfo(np.int64).max):
        raise ValueError('forest_ids exceed signed int64')
    starts = np.r_[0, np.flatnonzero(ids[1:] != ids[:-1])+1] if len(ids) else np.array([], dtype=int)
    # FITS inputs already group each forest. Avoid sorting the pixel catalog
    # unless a forest appears in several nonadjacent chunks.
    order = None
    if len(starts) != len(np.unique(ids[starts])):
        order = np.argsort(ids, kind='stable')
        sorted_ids = ids[order]
        starts = np.r_[0, np.flatnonzero(sorted_ids[1:] != sorted_ids[:-1])+1]
    stops = np.r_[starts[1:], len(ids)]
    out = delta.copy() if apply else None
    report = dict(forests=len(starts), pixels=len(ids), mean_only_forests=0,
        zero_weight_forests=0, max_abs_weighted_mean_before=0.,
        max_abs_normalized_first_moment_before=0., max_abs_weighted_mean_after=0.,
        max_abs_normalized_first_moment_after=0.,
        formula='delta_p = delta - <delta>_w - (log_lambda-<log_lambda>_w) * <(log_lambda-<log_lambda>_w)*delta>_w / <(log_lambda-<log_lambda>_w)^2>_w',
        coordinate='supplied true log wavelength; any logarithm base gives the same projector',
        weights='unchanged; engines multiply delta_p by w exactly once',
        applied=bool(apply))
    for start, stop in zip(starts, stops):
        index = slice(start, stop) if order is None else order[start:stop]
        projected, stats = _lya_project_one(delta[index], weights[index], lam[index], apply)
        if apply:
            out[index] = projected
        if stats is None:
            report['zero_weight_forests'] += 1
            continue
        report['mean_only_forests'] += int(stats['mean_only'])
        for name, field in (('mean_before', 'max_abs_weighted_mean_before'),
                            ('normalized_first_moment_before', 'max_abs_normalized_first_moment_before'),
                            ('mean_after', 'max_abs_weighted_mean_after'),
                            ('normalized_first_moment_after', 'max_abs_normalized_first_moment_after')):
            report[field] = max(report[field], stats[name])
    return out, report


def project_lya_delta(delta, weights, forest_ids, log_lambda, *, return_diagnostics=False):
    """Project each exact integer forest identity before 2PCF/3PCF searching.

    Inputs are never modified. The output is delta_p, not delta_p*w. Supply
    real log wavelength: Cartesian distance cannot substitute for it. IDs may
    be unsorted or exceed 2**53. Diagnostics report weighted constraint residuals.
    """
    projected, report = _lya_project_groups(delta, weights, forest_ids, log_lambda, True)
    return (projected, report) if return_diagnostics else projected


def forest_projection_diagnostics(delta, weights, forest_ids, log_lambda):
    """Measure weighted mean/slope residuals without changing the delta field."""
    return _lya_project_groups(delta, weights, forest_ids, log_lambda, False)[1]
