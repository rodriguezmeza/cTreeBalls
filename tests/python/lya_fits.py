"""eBOSS/PICCA delta products and explicit continuum preprocessing.

Projection follows equations (2.6)-(2.7) of arXiv:2507.00129. Raw spectra are
not delta products: continuum extraction must have been performed upstream.
"""
from __future__ import annotations

import math
import numbers
import numpy as np


def preprocess_forest(delta, weight, loglam, *, project_delta=False,
                      redshift_weight_exponent=0.0, weight_z_ref=2.25):
    delta, weight, loglam = (np.array(a, dtype=float, copy=True) for a in (delta, weight, loglam))
    if not math.isfinite(redshift_weight_exponent) or not math.isfinite(weight_z_ref) or weight_z_ref <= -1:
        raise ValueError("redshift weight exponent must be finite and weight_z_ref > -1")
    weight *= (10**loglam / 1215.67 / (1+weight_z_ref))**redshift_weight_exponent
    if not np.all(np.isfinite(weight)) or np.any(weight <= 0):
        raise ValueError("redshift weighting produced non-positive or non-finite weights")
    if project_delta:
        center = loglam - loglam[0]
        center -= np.average(center, weights=weight)
        delta -= np.average(delta, weights=weight)
        variance = np.dot(weight, center*center)
        if variance > 0:
            delta -= center * np.dot(weight*center, delta)/variance
    return delta, weight


def read_fits(paths, *, fits_layout="auto", eboss_angle_unit="rad", **kwargs):
    from astropy.io import fits
    from lya_corr_all_engines import expand_inputs, read_desi
    paths = expand_inputs(paths)
    if fits_layout == "auto":
        with fits.open(paths[0], memmap=False) as hdus:
            fits_layout = "desi" if {"LAMBDA", "METADATA", "WEIGHT"} <= {h.name for h in hdus} else "eboss"
    if fits_layout == "desi":
        return read_desi(paths, **kwargs)
    return read_eboss(paths, angle_unit=eboss_angle_unit, **kwargs)


def read_eboss(paths, *, omega_m=.315, h=.674, z_min=0., z_max=10.,
               max_forests=None, pixel_stride=1, delta_field="auto", angle_unit="rad",
               project_delta=False, redshift_weight_exponent=0., weight_z_ref=2.25):
    from astropy import units as u
    from astropy.cosmology import FlatLambdaCDM
    from astropy.io import fits
    from lya_corr_all_engines import ForestCatalog, expand_inputs, LYA_WAVELENGTH

    if not (math.isfinite(omega_m) and 0 < omega_m < 1 and math.isfinite(h) and h > 0):
        raise ValueError("require 0 < omega_m < 1 and finite h > 0")
    if not (math.isfinite(z_min) and math.isfinite(z_max) and 0 <= z_min < z_max):
        raise ValueError("require finite 0 <= z_min < z_max")
    if pixel_stride < 1 or (max_forests is not None and max_forests < 1):
        raise ValueError("pixel_stride and max_forests must be positive")
    if angle_unit not in ("rad", "deg"):
        raise ValueError("eBOSS header angle unit must be rad or deg")
    cosmology = FlatLambdaCDM(H0=100*h, Om0=omega_m, Tcmb0=0)
    pieces, provenance, seen = [], [], set()
    for path in expand_inputs(paths):
        if max_forests is not None and len(pieces) >= max_forests:
            break
        loaded, recognized, fields = 0, 0, set()
        with fits.open(path, memmap=False) as hdus:
            for hdu in hdus[1:]:
                if max_forests is not None and len(pieces) >= max_forests:
                    break
                if not isinstance(hdu, fits.BinTableHDU) or not {"LOGLAM", "WEIGHT"} <= set(hdu.columns.names):
                    continue
                recognized += 1
                field = delta_field
                if field == "auto":
                    field = "DELTA_BLIND" if "DELTA_BLIND" in hdu.columns.names else "DELTA"
                if field not in hdu.columns.names:
                    raise ValueError(f"{path}[{hdu.name}]: missing {field}")
                fields.add(field)
                identifier = next((hdu.header[k] for k in ("THING_ID", "LOS_ID", "TARGETID")
                                   if k in hdu.header), None)
                if not isinstance(identifier, numbers.Integral) or isinstance(identifier, bool):
                    raise ValueError(f"{path}[{hdu.name}]: expected integer THING_ID/LOS_ID/TARGETID header")
                identifier = int(identifier)
                if not np.iinfo(np.int64).min <= identifier <= np.iinfo(np.int64).max:
                    raise ValueError("forest ID exceeds signed int64")
                if identifier in seen:
                    raise ValueError(f"{path}: duplicate forest ID {identifier}")
                seen.add(identifier)
                ra, dec = [(float(hdu.header[k])*u.Unit(angle_unit)).to_value(u.rad) for k in ("RA", "DEC")]
                if not (math.isfinite(ra) and math.isfinite(dec) and abs(dec) <= np.pi/2):
                    raise ValueError(f"{path}: invalid sky coordinates for {identifier}")
                loglam, delta, weight = (np.asarray(hdu.data[k], dtype=float)
                                         for k in ("LOGLAM", field, "WEIGHT"))
                if any(a.ndim != 1 for a in (loglam, delta, weight)):
                    raise ValueError(f"{path}: expected scalar pixel rows, not vector-valued columns")
                with np.errstate(over="ignore", invalid="ignore"):
                    redshift = 10**loglam / LYA_WAVELENGTH - 1
                good = (np.isfinite(redshift) & (redshift > 0) & (redshift >= z_min) &
                        (redshift < z_max) & np.isfinite(delta) & np.isfinite(weight) & (weight > 0))
                indices = np.flatnonzero(good)[::pixel_stride]
                if not len(indices):
                    continue
                delta, weight = preprocess_forest(delta[indices], weight[indices], loglam[indices],
                    project_delta=project_delta, redshift_weight_exponent=redshift_weight_exponent,
                    weight_z_ref=weight_z_ref)
                chi = cosmology.comoving_distance(redshift[indices]).value*h
                los = np.array([np.cos(dec)*np.cos(ra), np.cos(dec)*np.sin(ra), np.sin(dec)])
                pieces.append((chi[:, None]*los, delta, weight, np.full(len(chi), identifier, dtype=np.int64)))
                loaded += 1
            if not recognized:
                raise ValueError(f"{path}: expected eBOSS/PICCA forest HDUs with LOGLAM, DELTA and WEIGHT")
            provenance.append(dict(path=str(path), delta_field=",".join(sorted(fields)),
                                   blinding="header/field as supplied; no unblinding", accepted_forests=loaded))
    if not pieces:
        raise ValueError("no valid forest pixels survived the selection")
    arrays = [np.concatenate([p[i] for p in pieces]) for i in range(4)]
    return ForestCatalog(*arrays, metadata=dict(
        input_format="eBOSS/PICCA forest HDUs", files=provenance, distance_unit="Mpc/h",
        cosmology=dict(model="FlatLambdaCDM", omega_m=omega_m, h=h, Tcmb0=0),
        wavelength_lya_angstrom=LYA_WAVELENGTH, z_min=z_min, z_max=z_max,
        max_forests=max_forests, pixel_stride=pixel_stride, angle_unit=angle_unit,
        project_delta=project_delta, projection_domain="retained pixels after cuts/stride",
        redshift_weight_exponent=redshift_weight_exponent, weight_z_ref=weight_z_ref,
        weights="input WEIGHT times ((1+z)/(1+z_ref))**exponent",
    )).normalized()
