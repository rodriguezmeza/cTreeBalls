"""Forest FITS readers (DESI, eBOSS/PICCA, Cartesian and RA-DEC-R).

Projection follows equations (2.6)-(2.7) of arXiv:2507.00129. Raw spectra are
not delta products: continuum extraction must have been performed upstream.
"""
from __future__ import annotations

import sys
from pathlib import Path

# Locate sibling analysis and benchmark helpers for direct script execution.
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python"))


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
        from lya_projection import project_forest_delta
        delta = project_forest_delta(delta, weight, loglam)
    return delta, weight


def read_fits(paths, *, fits_layout="auto", eboss_angle_unit="rad",
              radec_r_angle_unit="deg", radec_r_distance_unit="Mpc/h",
              radec_r_distance_source="input", **kwargs):
    from astropy.io import fits
    from lya_corr_all_engines import expand_inputs, read_desi
    paths = expand_inputs(paths)
    if fits_layout == "auto":
        with fits.open(paths[0], memmap=False) as hdus:
            if {"LAMBDA", "METADATA", "WEIGHT"} <= {h.name for h in hdus}:
                fits_layout = "desi"
            elif any(isinstance(h, fits.BinTableHDU) and
                     {"x", "y", "z", "los_id", "weight"} <= {k.lower() for k in h.columns.names}
                     for h in hdus[1:]):
                fits_layout = "cartesian"
            elif any(isinstance(h, fits.BinTableHDU) and
                     {"ra", "dec", "rcdis", "tid", "wdelta", "wobs", "zobs"} <= {k.lower() for k in h.columns.names}
                     for h in hdus[1:]):
                fits_layout = "radec-r"
            else:
                fits_layout = "eboss"
    if radec_r_distance_source != "input" and fits_layout != "radec-r":
        raise ValueError("--radec-r-distance-source requires RA-DEC-R FITS input")
    if fits_layout == "desi":
        return read_desi(paths, **kwargs)
    if fits_layout == "cartesian":
        return read_cartesian(paths, **kwargs)
    if fits_layout == "radec-r":
        return read_radec_r(paths, angle_unit=radec_r_angle_unit,
                           distance_unit=radec_r_distance_unit,
                           distance_source=radec_r_distance_source, **kwargs)
    if fits_layout != "eboss":
        raise ValueError(f"unsupported FITS layout: {fits_layout}")
    return read_eboss(paths, angle_unit=eboss_angle_unit, **kwargs)


def read_eboss(paths, *, omega_m=.315, h=.674, z_min=0., z_max=10.,
               max_forests=None, pixel_stride=1, delta_field="auto", angle_unit="rad",
               project_delta=False, redshift_weight_exponent=0., weight_z_ref=2.25,
               retain_projection_domain=False):
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
    pieces, provenance, seen, domains, wavelengths = [], [], set(), [], []
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
                full_good = (np.isfinite(redshift) & (redshift > 0) & np.isfinite(delta)
                             & np.isfinite(weight) & (weight > 0))
                good = full_good & (redshift >= z_min) & (redshift < z_max)
                indices = np.flatnonzero(good)[::pixel_stride]
                if not len(indices):
                    continue
                retained_delta, retained_weight = preprocess_forest(delta[indices], weight[indices], loglam[indices],
                    project_delta=project_delta, redshift_weight_exponent=redshift_weight_exponent,
                    weight_z_ref=weight_z_ref)
                full = np.flatnonzero(full_good) if retain_projection_domain and not project_delta else indices
                chi = cosmology.comoving_distance(redshift[full]).value*h
                los = np.array([np.cos(dec)*np.cos(ra), np.cos(dec)*np.sin(ra), np.sin(dec)])
                positions = chi[:,None]*los
                selected = np.searchsorted(full,indices)
                pieces.append((positions[selected], retained_delta, retained_weight,
                               np.full(len(indices),identifier,dtype=np.int64)))
                wavelengths.append(loglam[indices])
                if retain_projection_domain and not project_delta:
                    domains.append(dict(positions=positions,
                        forest_ids=np.full(len(full),identifier,dtype=np.int64),
                        log_lambda=loglam[full],projection_weights=weight[full],data_indices=selected))
                loaded += 1
            if not recognized:
                raise ValueError(f"{path}: expected eBOSS/PICCA forest HDUs with LOGLAM, DELTA and WEIGHT")
            provenance.append(dict(path=str(path), delta_field=",".join(sorted(fields)),
                                   blinding="header/field as supplied; no unblinding", accepted_forests=loaded))
    if not pieces:
        raise ValueError("no valid forest pixels survived the selection")
    arrays = [np.concatenate([p[i] for p in pieces]) for i in range(4)]
    from lya_distortion import concatenate_domains
    return ForestCatalog(*arrays, metadata=dict(
        input_format="eBOSS/PICCA forest HDUs", files=provenance, distance_unit="Mpc/h",
        cosmology=dict(model="FlatLambdaCDM", omega_m=omega_m, h=h, Tcmb0=0),
        wavelength_lya_angstrom=LYA_WAVELENGTH, z_min=z_min, z_max=z_max,
        max_forests=max_forests, pixel_stride=pixel_stride, angle_unit=angle_unit,
        project_delta=project_delta, projection_domain=("retained pixels after cuts/stride" if project_delta else "input projection unchanged; original valid forest assumed for response"),
        redshift_weight_exponent=redshift_weight_exponent, weight_z_ref=weight_z_ref,
        weights="input WEIGHT times ((1+z)/(1+z_ref))**exponent",
    ), log_lambda=np.concatenate(wavelengths),
        distortion_domain=concatenate_domains(domains)).normalized()


def read_cartesian(
    paths,
    *,
    max_forests=None,
    pixel_stride=1,
    delta_field="auto",
    omega_m=0.315,
    h=0.674,
    z_min=0.0,
    z_max=10.0,
    project_delta=False,
    redshift_weight_exponent=0.0,
    weight_z_ref=2.25,
    retain_projection_domain=False,
):
    """Read already-comoving x/y/z/delta/weight/los_id tables, without reprojection.

    As with ASCII/NPZ input, coordinates must already be in Mpc/h. Some private
    exported catalogs omit TUNIT; this assumption is recorded in provenance.
    """
    from astropy.io import fits
    from lya_corr_all_engines import ForestCatalog, expand_inputs

    if (
        omega_m != 0.315
        or h != 0.674
        or z_min != 0
        or z_max != 10
        or project_delta
        or redshift_weight_exponent != 0
        or weight_z_ref != 2.25
    ):
        raise ValueError(
            "Cartesian FITS coordinates are already comoving; cosmology, redshift cuts, "
            "projection and wavelength weighting require DESI/eBOSS wavelength data"
        )
    if pixel_stride < 1 or (max_forests is not None and max_forests < 1):
        raise ValueError("pixel_stride and max_forests must be positive")
    pieces, provenance, seen, domains = [], [], set(), []
    for path in expand_inputs(paths):
        if max_forests is not None and len(pieces) >= max_forests:
            break
        with fits.open(path, memmap=True) as hdus:
            recognized, loaded, fields = False, 0, set()
            for hdu in hdus[1:]:
                if not isinstance(hdu, fits.BinTableHDU):
                    continue
                columns = {k.lower(): k for k in hdu.columns.names}
                if not {"x", "y", "z", "weight", "los_id"} <= columns.keys():
                    continue
                recognized = True
                field = (
                    ("delta_blind" if "delta_blind" in columns else "delta")
                    if delta_field == "auto"
                    else delta_field.lower()
                )
                if field not in columns:
                    raise ValueError(f"{path}: missing {field}")
                fields.add(field)
                for axis in ("x", "y", "z"):
                    unit = hdu.columns[columns[axis]].unit
                    if unit and unit.lower().replace(" ", "") not in {
                        "mpc/h",
                        "mpch-1",
                        "mpch^-1",
                    }:
                        raise ValueError(
                            f"{path}: Cartesian coordinates require Mpc/h; found {unit}"
                        )
                ids = np.asarray(hdu.data[columns["los_id"]])
                if ids.ndim != 1 or ids.dtype.kind not in "iu":
                    raise ValueError(f"{path}: LOS_ID must be a scalar integer column")
                if ids.dtype.kind == "u" and np.any(ids > np.iinfo(np.int64).max):
                    raise ValueError("forest ID exceeds signed int64")
                # Preserve catalog forest order, then stride within each valid forest.
                unique, first = np.unique(ids, return_index=True)
                forest_order = unique[np.argsort(first)]
                remaining = None if max_forests is None else max_forests - len(pieces)
                # Only selected forests need Cartesian materialization.
                chosen = forest_order if remaining is None else forest_order[:remaining]
                rows = np.flatnonzero(np.isin(ids, chosen))
                positions = np.column_stack(
                    [hdu.data[columns[k]][rows] for k in ("x", "y", "z")]
                )
                delta, weight = [
                    np.asarray(hdu.data[columns[k]][rows], dtype=float)
                    for k in (field, "weight")
                ]
                selected_ids = ids[rows]
                valid = (
                    np.isfinite(positions).all(axis=1)
                    & (np.linalg.norm(positions, axis=1) > 0)
                    & np.isfinite(delta)
                    & np.isfinite(weight)
                    & (weight > 0)
                )
                order = np.argsort(selected_ids, kind="stable")
                ordered_ids = selected_ids[order]
                for identifier in chosen:
                    identifier = int(identifier)
                    if identifier in seen:
                        raise ValueError(f"{path}: duplicate forest ID {identifier}")
                    seen.add(identifier)
                    left, right = np.searchsorted(
                        ordered_ids, identifier, side="left"
                    ), np.searchsorted(ordered_ids, identifier, side="right")
                    group = order[left:right]
                    full = group[valid[group]]
                    take = full[::pixel_stride]
                    if not len(take):
                        continue
                    pieces.append(
                        (
                            positions[take],
                            delta[take],
                            weight[take],
                            np.full(len(take), identifier, dtype=np.int64),
                        )
                    )
                    if retain_projection_domain:
                        domains.append(dict(positions=positions[full],
                            forest_ids=np.full(len(full),identifier,dtype=np.int64),
                            projection_weights=weight[full],data_indices=np.arange(0,len(full),pixel_stride)))
                    loaded += 1
                if max_forests is not None and len(pieces) >= max_forests:
                    break
            if not recognized:
                raise ValueError(
                    f"{path}: expected Cartesian FITS columns x, y, z, delta, weight, los_id"
                )
            provenance.append(
                dict(
                    path=str(path),
                    delta_field=",".join(sorted(fields)),
                    blinding="as supplied; no unblinding",
                    accepted_forests=loaded,
                )
            )
    if not pieces:
        raise ValueError("no valid forest pixels survived the selection")
    arrays = [np.concatenate([p[i] for p in pieces]) for i in range(4)]
    from lya_distortion import concatenate_domains
    return ForestCatalog(
        *arrays,
        metadata=dict(
            input_format="Cartesian forest FITS",
            files=provenance,
            distance_unit="Mpc/h required; existing coordinates preserved (unlabelled exports assumed Mpc/h)",
            max_forests=max_forests,
            pixel_stride=pixel_stride,
            weights="input WEIGHT unchanged; zero/negative/nonfinite pixels excluded",
        ), distortion_domain=concatenate_domains(domains),
    ).normalized()


def _radec_r_ids(values, where):
    """Accept integer columns and provably safe integer-valued float exports."""
    values = np.asarray(values)
    if values.ndim != 1:
        raise ValueError(f"{where}: tid must be a scalar column")
    if values.dtype.kind in "iu":
        if values.dtype.kind == "u" and np.any(values > np.iinfo(np.int64).max):
            raise ValueError(f"{where}: tid exceeds signed int64")
    elif values.dtype.kind == "f":
        # The supplied catalog exports small integer target IDs as doubles.
        # Larger floating IDs could already have lost integer identity upstream.
        limit = min(2**53 - 1, 2**(np.finfo(values.dtype).nmant + 1) - 1)
        if (not np.isfinite(values).all() or np.any(values != np.trunc(values))
                or np.any(np.abs(values) > limit)):
            raise ValueError(f"{where}: floating tid must be exact integers within 2**53-1")
    else:
        raise ValueError(f"{where}: tid must contain integer forest identities")
    return values.astype(np.int64, copy=False)


def read_radec_r(paths, *, angle_unit="deg", distance_unit="Mpc/h", distance_source="input",
                 omega_m=.315, h=.674, z_min=0., z_max=10.,
                 max_forests=None, pixel_stride=1, delta_field="auto",
                 project_delta=False, redshift_weight_exponent=0., weight_z_ref=2.25,
                 retain_projection_domain=False):
    """Read vector-valued SPECTRA rows; join all chunks of each target ID.

    The default uses supplied rcdis as an absolute observer distance. Explicit
    distance_source="redshift" instead uses FlatLambdaCDM chi(zobs) in Mpc/h.
    TUNIT wins over the angle/distance fallbacks for unlabelled exports.
    Selection and stride apply to a complete joined sightline, not to chunks.
    """
    from astropy import units as u
    from astropy.io import fits
    from lya_corr_all_engines import ForestCatalog, expand_inputs, LYA_WAVELENGTH
    from lya_distortion import concatenate_domains

    if distance_source not in ("input", "redshift"):
        raise ValueError("RA-DEC-R distance_source must be input or redshift")
    if distance_source == "input" and omega_m != .315:
        raise ValueError("RA-DEC-R FITS supplies rcdis; --omega-m cannot recompute its distances")
    if not (math.isfinite(omega_m) and 0 < omega_m < 1):
        raise ValueError("require finite 0 < omega_m < 1")
    if not (math.isfinite(h) and h > 0):
        raise ValueError("require finite h > 0 for Mpc to Mpc/h conversion")
    if not (math.isfinite(z_min) and math.isfinite(z_max) and 0 <= z_min < z_max):
        raise ValueError("require finite 0 <= z_min < z_max")
    if pixel_stride < 1 or (max_forests is not None and max_forests < 1):
        raise ValueError("pixel_stride and max_forests must be positive")
    if angle_unit not in ("deg", "rad") or distance_unit not in ("Mpc", "Mpc/h"):
        raise ValueError("RA-DEC-R units must be deg/rad and Mpc/Mpc/h")
    paths = expand_inputs(paths)
    cosmology = None
    if distance_source == "redshift":
        from astropy.cosmology import FlatLambdaCDM
        cosmology = FlatLambdaCDM(H0=100*h, Om0=omega_m, Tcmb0=0)
    required = {"ra", "dec", "rcdis", "tid", "wdelta", "wobs", "zobs"}

    def tables(hdus, path):
        recognized = False
        for hdu in hdus[1:]:
            if not isinstance(hdu, fits.BinTableHDU):
                continue
            columns = {name.lower(): name for name in hdu.columns.names}
            if not required <= columns.keys():
                continue
            recognized = True
            field = ("delta_blind" if "delta_blind" in columns else "delta") if delta_field == "auto" else delta_field.lower()
            if field not in columns:
                raise ValueError(f"{path}[{hdu.name}]: missing {field}")
            yield hdu, columns, field
        if not recognized:
            raise ValueError(f"{path}: expected RA-DEC-R forest rows with {', '.join(sorted(required))}")

    def pixels(row, columns, field, where):
        arrays = [np.array(row[columns[k]], dtype=float, copy=True)
                  for k in ("rcdis", field, "wdelta", "wobs", "zobs")]
        if any(a.ndim != 1 or a.shape != arrays[0].shape for a in arrays):
            raise ValueError(f"{where}: rcdis/delta/wdelta/wobs/zobs must be equal-length pixel vectors")
        radius, delta, weight, wave, redshift = arrays
        good = (np.isfinite(delta)
                & np.isfinite(weight) & (weight > 0) & np.isfinite(wave)
                & (wave > 0) & np.isfinite(redshift) & (redshift > 0))
        if distance_source == "input":
            good &= np.isfinite(radius) & (radius > 0)
        return [a[good] for a in arrays]

    # First select complete forest identities. Even once the quota is full,
    # the collection pass below reads their chunks in ALL redshift extensions.
    # This also recovers earlier out-of-cut pixels for a full response domain.
    selected = None
    if max_forests is not None:
        selected = {}
        for path in paths:
            with fits.open(path, memmap=True) as hdus:
                for hdu, columns, field in tables(hdus, path):
                    for row in hdu.data:
                        identifier = int(_radec_r_ids([row[columns["tid"]]], str(path))[0])
                        if identifier in selected:
                            continue
                        *_, redshift = pixels(row, columns, field, f"{path}[{hdu.name}]")
                        if np.any((redshift >= z_min) & (redshift < z_max)):
                            selected[identifier] = None
                            if len(selected) == max_forests:
                                break
                    if len(selected) == max_forests:
                        break
            if len(selected) == max_forests:
                break
        if not selected:
            raise ValueError("no valid forest pixels survived the selection")

    forests = {} if selected is None else {identifier: None for identifier in selected}
    provenance = []
    for path in paths:
        accepted, fields, units = set(), set(), []
        with fits.open(path, memmap=True) as hdus:
            for hdu, columns, field in tables(hdus, path):
                where = f"{path}[{hdu.name}]"
                ids = _radec_r_ids(hdu.data[columns["tid"]], where)
                rows = range(len(ids)) if selected is None else np.flatnonzero(np.isin(ids, list(selected)))
                for key in ("ra", "dec"):
                    if np.asarray(hdu.data[columns[key]]).ndim != 1:
                        raise ValueError(f"{where}: {key} must be a scalar sky-coordinate column")
                ra_unit = u.Unit(hdu.columns[columns["ra"]].unit or angle_unit)
                dec_unit = u.Unit(hdu.columns[columns["dec"]].unit or angle_unit)
                ra_factor, dec_factor = ra_unit.to(u.rad), dec_unit.to(u.rad)
                radial_unit = hdu.columns[columns["rcdis"]].unit or distance_unit
                canonical = radial_unit.lower().replace(" ", "")
                if distance_source == "redshift":
                    radial_factor = None  # rcdis is not used to place pixels.
                elif canonical in {"mpc/h", "mpch-1", "mpch^-1", "h-1mpc", "h^-1mpc"}:
                    radial_factor = 1.
                elif canonical == "mpc":
                    radial_factor = h
                else:
                    raise ValueError(f"{where}: rcdis requires Mpc or Mpc/h; found {radial_unit}")
                wave_unit = u.Unit(hdu.columns[columns["wobs"]].unit or "Angstrom")
                wave_factor = wave_unit.to(u.AA)
                fields.add(field)
                units.append(dict(extension=hdu.name, ra=str(ra_unit), dec=str(dec_unit),
                                  rcdis=radial_unit, radius_scale_to_mpc_h=radial_factor,
                                  rcdis_unit_is_assumed=hdu.columns[columns["rcdis"]].unit is None,
                                  radius_source=distance_source,
                                  wobs=str(wave_unit)))
                for index in rows:
                    row = hdu.data[index]
                    identifier = int(ids[index])
                    radius, delta, weight, wave, redshift = pixels(row, columns, field, where)
                    if not len(radius):
                        continue
                    if distance_source == "redshift":
                        wave_z = wave*wave_factor/LYA_WAVELENGTH - 1.
                        if not np.allclose(redshift, wave_z, rtol=0, atol=1e-5):
                            raise ValueError(f"{where}: zobs disagrees with wobs/1215.67 - 1; "
                                             "cannot assign consistent absorption redshifts")
                        radius = cosmology.comoving_distance(redshift).value*h
                    else:
                        radius = radius*radial_factor
                    ra, dec = float(row[columns["ra"]])*ra_factor, float(row[columns["dec"]])*dec_factor
                    if not (math.isfinite(ra) and math.isfinite(dec) and abs(dec) <= np.pi/2):
                        raise ValueError(f"{where}: invalid sky coordinates for tid {identifier}")
                    los = np.array([np.cos(dec)*np.cos(ra), np.cos(dec)*np.sin(ra), np.sin(dec)])
                    if forests.get(identifier) is None:
                        forests[identifier] = dict(los=los, chunks=[])
                    elif not np.allclose(forests[identifier]["los"], los, rtol=0, atol=1e-12):
                        raise ValueError(f"{where}: inconsistent sky coordinates for repeated tid {identifier}")
                    forests[identifier]["chunks"].append((radius, delta, weight, wave*wave_factor, redshift))
                    accepted.add(identifier)
        provenance.append(dict(path=str(path), delta_field=",".join(sorted(fields)),
                               blinding="as supplied; no unblinding", accepted_forests=len(accepted), units=units))

    pieces, wavelengths, domains = [], [], []
    for identifier, forest in forests.items():
        if forest is None:
            continue
        radius, delta, weight, wave, redshift = [np.concatenate([chunk[i] for chunk in forest["chunks"]]) for i in range(5)]
        # Chunk boundaries do not create new forests or restart the stride.
        order = np.argsort(wave, kind="stable")
        radius, delta, weight, wave, redshift = [a[order] for a in (radius, delta, weight, wave, redshift)]
        if np.any(np.diff(wave) <= 0):
            raise ValueError(f"duplicate wavelength pixels for tid {identifier}; overlapping/repeated forest chunks")
        take = np.flatnonzero((redshift >= z_min) & (redshift < z_max))[::pixel_stride]
        if not len(take):
            forest.clear()
            continue
        positions = radius[:, None]*forest["los"]
        forest.clear()
        loglam = np.log10(wave)
        kept_delta, kept_weight = preprocess_forest(delta[take], weight[take], loglam[take],
            project_delta=project_delta, redshift_weight_exponent=redshift_weight_exponent,
            weight_z_ref=weight_z_ref)
        pieces.append((positions[take], kept_delta, kept_weight, np.full(len(take), identifier, dtype=np.int64)))
        wavelengths.append(loglam[take])
        if retain_projection_domain and not project_delta:
            domains.append(dict(positions=positions, forest_ids=np.full(len(radius),identifier,dtype=np.int64),
                                log_lambda=loglam, projection_weights=weight, data_indices=take))
    if not pieces:
        raise ValueError("no valid forest pixels survived the selection")
    arrays = [np.concatenate([p[i] for p in pieces]) for i in range(4)]
    metadata = dict(
        input_format="RA-DEC-R forest FITS", files=provenance, distance_unit="Mpc/h",
        distance_source=("input rcdis (not recomputed from redshift)" if distance_source == "input"
                         else "FlatLambdaCDM observer distance from zobs; input rcdis not used"),
        distance_source_mode=distance_source,
        unlabelled_angle_unit=angle_unit, unlabelled_distance_unit=distance_unit, h=h,
        forest_identity="tid; all redshift chunks joined; cid is not a forest identity",
        z_min=z_min, z_max=z_max, max_forests=max_forests, pixel_stride=pixel_stride,
        project_delta=project_delta, redshift_weight_exponent=redshift_weight_exponent,
        weight_z_ref=weight_z_ref, wavelength_lya_angstrom=LYA_WAVELENGTH,
        projection_domain=("retained pixels after cuts/stride" if project_delta else "input projection unchanged; original valid forest assumed for response"),
        weights="input wdelta times ((1+z)/(1+z_ref))**exponent; invalid pixels excluded jointly")
    if distance_source == "redshift":
        metadata["cosmology"] = dict(model="FlatLambdaCDM", omega_m=omega_m, h=h, Tcmb0=0)
    return ForestCatalog(*arrays, metadata=metadata, log_lambda=np.concatenate(wavelengths),
                         distortion_domain=concatenate_domains(domains)).normalized()
