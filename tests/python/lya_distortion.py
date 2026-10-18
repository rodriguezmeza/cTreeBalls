"""Native 2PCF continuum response and portable matrix products.

Implements equations 2.5--2.7 of arXiv:2507.00129. This is a forward
operator, not an inverse correction, covariance, or 3PCF distortion model.
"""
from __future__ import annotations

import json
from pathlib import Path
import time

import numpy as np

DOMAIN_FIELDS = ('positions', 'forest_ids', 'projection_weights', 'data_indices')


def normalize_domain(domain, positions, ids, log_lambda):
    if domain is None:
        return None
    if set(domain) not in (set(DOMAIN_FIELDS),set(DOMAIN_FIELDS)|{'log_lambda'}):
        raise ValueError('distortion_domain requires ' + ', '.join(DOMAIN_FIELDS))
    result = {k: np.asarray(v) for k,v in domain.items()}
    n = len(result['positions'])
    if result['positions'].shape != (n, 3) or n < len(positions):
        raise ValueError('invalid distortion-domain positions')
    for key in ('forest_ids', 'projection_weights') + (('log_lambda',) if 'log_lambda' in result else ()):
        if result[key].shape != (n,):
            raise ValueError(f'invalid distortion-domain {key} shape')
    for key in ('forest_ids', 'data_indices'):
        a = result[key]
        if a.dtype.kind not in 'iu' or (a.dtype.kind == 'u' and np.any(a > np.iinfo(np.int64).max)):
            raise ValueError(f'distortion-domain {key} must be signed int64 compatible')
        result[key] = np.ascontiguousarray(a,dtype=np.int64)
    index = result['data_indices']
    if (index.shape != (len(positions),) or np.any(index < 0) or np.any(index >= n)
            or len(np.unique(index)) != len(index)):
        raise ValueError('invalid distortion-domain data_indices mapping')
    for key in ('positions','projection_weights') + (('log_lambda',) if 'log_lambda' in result else ()):
        result[key] = np.ascontiguousarray(result[key],dtype=np.float64)
        if not np.isfinite(result[key]).all():
            raise ValueError(f'non-finite distortion-domain {key}')
    if (np.any(result['projection_weights'] < 0)
            or np.any(np.linalg.norm(result['positions'],axis=1) <= 0)):
        raise ValueError('invalid distortion-domain weights or distances')
    if (not np.array_equal(result['positions'][index],positions)
            or not np.array_equal(result['forest_ids'][index],ids)
            or (log_lambda is None) != ('log_lambda' not in result)
            or (log_lambda is not None and not np.array_equal(result['log_lambda'][index],log_lambda))):
        raise ValueError('distortion-domain mapping does not match catalog pixels')
    if not np.array_equal(np.unique(result['forest_ids']),np.unique(ids)):
        raise ValueError('distortion domain must contain precisely the selected forests')
    return result


def concatenate_domains(pieces):
    """Combine per-forest full fitting domains with selected-pixel mappings."""
    if not pieces:
        return None
    offsets = np.cumsum([0]+[len(p['positions']) for p in pieces[:-1]])
    return {k: np.concatenate([p[k]+offset if k == 'data_indices' else p[k]
                              for p,offset in zip(pieces,offsets)]) for k in pieces[0]}


def compute_catalog_distortion(catalog, config):
    from cyballs import lya_distortion_matrix
    if catalog.metadata.get('project_delta') and config.distortion_projection_order != 1:
        raise ValueError('--project-delta removes mean+slope; its distortion order must be 1')
    domain = catalog.distortion_domain
    if domain is None:
        if (not catalog.metadata.get('project_delta')
                and (catalog.metadata.get('pixel_stride',1) != 1
                     or catalog.metadata.get('z_min',0) != 0 or catalog.metadata.get('z_max',10) != 10)):
            raise ValueError('cache lacks the full continuum-fitting domain; reload FITS with --compute-distortion --save-catalog')
        positions, ids, lam = catalog.positions, catalog.forest_ids, catalog.log_lambda
        weights, projection_weights = catalog.weights, catalog.weights
        scope = 'retained pixels with correlation weights'
    else:
        positions, ids = (domain[k] for k in ('positions','forest_ids'))
        lam = domain.get('log_lambda')
        projection_weights = domain['projection_weights']
        weights = np.zeros(len(positions))
        weights[domain['data_indices']] = catalog.weights
        scope = 'full valid input forest with original WEIGHT; data cuts/stride affect outer pairs only'
    started = time.perf_counter()
    result = lya_distortion_matrix(positions,weights,ids,lam,
        projection_weights=projection_weights,projection_order=config.distortion_projection_order,
        rp_max=config.rp_max,rt_max=config.rt_max,rp_bins=config.rp_bins,rt_bins=config.rt_bins,
        number_threads=config.threads,max_memory_mb=config.distortion_max_mib)
    result['metadata'] = dict(
        formula='arXiv:2507.00129 equations 2.5--2.7',
        convention='xi_measured[A] = sum_B D[A,B] xi_unprojected[B]',
        bin_order='rp major, rt minor; folded observer-bisector separations in Mpc/h',
        projection_domain=scope,
        projection_order=config.distortion_projection_order,
        model_support='zero correlation outside model grid; rows are not renormalized',
        scope='distinct-forest 3D 2PCF continuum mean/slope only; no survey mean-flux subtraction or redshift evolution of the model',
        wall_seconds=time.perf_counter()-started,
        stats=result['stats'])
    return result


def save_distortion(directory, result):
    """Save D, unnormalised sums and row weights; compatible with lya_analysis."""
    from astropy.io import fits
    directory = Path(directory)
    directory.mkdir(parents=True,exist_ok=True)
    if any(directory.iterdir()):
        raise FileExistsError(f'{directory}: distortion products already exist; choose a new output directory')
    arrays = {k: v for k,v in result.items() if k not in ('matrix','stats','metadata')}
    np.savez_compressed(directory/'distortion.npz',distortion=result['matrix'],
                        metadata=json.dumps(result.get('metadata',{})),**arrays)
    na,nb = result['matrix'].shape
    hdu = fits.BinTableHDU.from_columns([
        fits.Column(name='DM',format=f'{nb}D',array=result['matrix']),
        fits.Column(name='RAW_DM',format=f'{nb}D',array=result['numerator']),
        fits.Column(name='WDM',format='D',array=result['weights']),
        fits.Column(name='VALID',format='L',array=result['valid_rows']),
        fits.Column(name='RP',format='D',unit='Mpc/h',array=result['rp']),
        fits.Column(name='RT',format='D',unit='Mpc/h',array=result['rt']),
    ],name='DISTORTION')
    for key,value in dict(NP=result['rp_bins'],NT=result['rt_bins'],RPMAX=result['rp_max'],
        RTMAX=result['rt_max'],RPMIN=0.,MODNP=result['model_rp_bins'],MODNT=result['model_rt_bins'],
        MODRPMAX=result['model_rp_max'],MODRTMAX=result['model_rt_max'],PROJORD=result['projection_order']).items():
        hdu.header[key] = value
    hdu.header['DMCONV'] = 'xi_data = DM @ xi_model'
    hdu.header['BINORDER'] = 'rp-major rt-minor'
    hdu.header['HISTORY'] = 'arXiv:2507.00129 eqs 2.5-2.7; no row-sum renormalization'
    hdu.header['HISTORY'] = 'Undefined (zero WDM) rows are NaN; model is zero outside its grid'
    model = fits.BinTableHDU.from_columns([
        fits.Column(name='RP',format='D',unit='Mpc/h',array=result['model_rp']),
        fits.Column(name='RT',format='D',unit='Mpc/h',array=result['model_rt']),
    ],name='MODEL_BINS')
    fits.HDUList([fits.PrimaryHDU(),hdu,model]).writeto(directory/'distortion.fits')
    report = dict(result.get('metadata',{}), shape=[na,nb], valid_rows=int(result['valid_rows'].sum()),
                  npz=str(directory/'distortion.npz'),fits=str(directory/'distortion.fits'))
    (directory/'distortion.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    return report
