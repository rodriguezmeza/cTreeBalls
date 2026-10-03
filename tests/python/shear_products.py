"""Retain native shear arrays and attach an explicit analysis-validity policy."""
import numpy as np

SCHEMA='desy3-shear-v1'
PROJECTION='Porth-x; local-east-north; great-circle transport'
CONDITION_LIMIT=1.e10
RESIDUAL_LIMIT=1.e-8


def attach_products(result, config, catalog):
    result['schema']=SCHEMA if catalog.metadata.get('fits_format') == 'desy3' else 'spherical-shear-v1'
    result['catalog']={**catalog.metadata,'nbody':catalog.nbody,'geometry':catalog.geometry}
    result['shear_projection']=PROJECTION
    lo,hi=config.native_limits(catalog.geometry)
    edges=np.geomspace(lo,hi,config.bins+1) if config.use_log_bins else np.linspace(lo,hi,config.bins+1)
    result['chord_edges']=edges
    result['theta_edges_arcmin']=np.rad2deg(2*np.arcsin(edges/2))*60
    result['theta_arcmin']=np.rad2deg(2*np.arcsin(result['radius']/2))*60
    if config.statistics in ('2pcf','both'):
        for key in ('xi_plus','xi_minus','pair_weight'):
            if result[key].shape != (config.bins,) or not np.all(np.isfinite(result[key])):
                raise ValueError(f'malformed native {key}')
        result['pair_valid']=result['pair_weight']>0
    if config.statistics in ('3pcf','both'):
        n=config.multipoles;shape=(4,config.bins,config.bins,2*n+1)
        for key in ('upsilon','gamma'):
            if result[key].shape!=shape or not np.all(np.isfinite(result[key])):
                raise ValueError(f'malformed native {key}')
        if not np.array_equal(result['orders'],np.arange(-n,n+1)):
            raise ValueError('unexpected native shear multipole axis')
        if result['window'].shape!=(config.bins,config.bins,4*n+1) or not np.all(np.isfinite(result['window'])):
            raise ValueError('malformed native window')
        result.update(window_diagnostics(result['upsilon'],result['window'],result['gamma'],n))
        result['normalizations']={
            'upsilon':'raw distinct-triplet spin-2 sums',
            'normalized':'upsilon / window monopole; window-convolved',
            'gamma':'native finite-window mode-coupling solution',
            'gamma_valid_policy':f'positive window monopole; cond(C)<={CONDITION_LIMIT:g}; relative solve residual<={RESIDUAL_LIMIT:g}',
            'orders':'all signed orders -nmax..+nmax; no inferred symmetry or factor-of-two rescaling',
        }
    return result


def window_diagnostics(upsilon,window,gamma,n):
    n0=window[:,:,2*n]
    valid=np.isfinite(n0)&(n0.real>0)&(np.abs(n0.imag)<=1e-10*np.maximum(n0.real,1.))
    normalized=np.full_like(upsilon,np.nan+1j*np.nan)
    np.divide(upsilon,n0[None,:,:,None],out=normalized,where=valid[None,:,:,None])
    condition=np.full(n0.shape,np.inf);residual=np.full(n0.shape,np.inf)
    orders=np.arange(-n,n+1);indices=orders[:,None]-orders[None,:]+2*n
    for j,i in zip(*np.where(valid)):
        matrix=window[j,i,indices]/n0[j,i]
        try:condition[j,i]=np.linalg.cond(matrix)
        except np.linalg.LinAlgError:continue
        rhs=normalized[:,j,i,:].T;answer=gamma[:,j,i,:].T
        difference=np.max(np.abs(matrix@answer-rhs),axis=0)
        scale=np.maximum(np.max(np.abs(rhs),axis=0),np.finfo(float).tiny)
        residual[j,i]=np.max(difference/scale)
    gamma_valid=valid&np.isfinite(condition)&(condition<=CONDITION_LIMIT)&(residual<=RESIDUAL_LIMIT)
    return dict(normalized=normalized,window_monopole=n0,window_valid=valid,
                gamma_valid=gamma_valid,window_condition=condition,gamma_solve_relative_residual=residual)
