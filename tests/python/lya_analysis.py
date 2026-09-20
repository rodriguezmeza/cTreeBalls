"""Estimator-aware comparisons and notebook-style Ly-alpha 2PCF analysis.

Wedges follow the geometric-bin overlap and inverse-diagonal-covariance
weighting in lya2pcf/plot_auxiliars.py, using NumPy with explicit empty-bin
handling. A distortion matrix is a forward operator (paper equation 2.5).
"""
from __future__ import annotations

from itertools import combinations
import math
from pathlib import Path
import numpy as np


def pyplot():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def write_comparisons(results, config):
    from lya_corr_all_engines import PRODUCTS, compare_products
    comparisons, paths = {}, []
    for key, (_, ndim, _) in PRODUCTS.items():
        members = [(name, row["products"][key]) for name, row in results.items() if key in row["products"]]
        # Reference first makes the direction of relative errors unambiguous.
        members.sort(key=lambda item: not item[0].startswith("lya2pcf-"))
        for (base, left), (name, right) in combinations(members, 2):
            metrics, rows = compare_products(left, right, key, relative_floor=config.relative_floor,
                                             rtol=config.rtol, atol=config.atol)
            tag = f"{key}__{name}__vs__{base}"
            comparisons[tag] = dict(reference=base, candidate=name, product=key, **metrics)
            np.savetxt(config.output_dir/(tag+".csv"), rows, delimiter=",", comments="",
                       header=",".join([f"bin{i}" for i in range(ndim)] +
                                       ["reference", "candidate", "difference", "relative"]))
        if not config.plots or len(members) < 2:
            continue
        plt = pyplot()
        fig, axes = plt.subplots(3, 1, figsize=(11, 9), sharex=True)
        base, left = members[0]
        all_keys = sorted({tuple(row[:ndim].astype(int)) for _, table in members for row in table})
        indices = {key: i for i, key in enumerate(all_keys)}
        for number, (name, table) in enumerate(members):
            values = np.full(len(all_keys), np.nan)
            for row in table:
                if row[-1] > 0:
                    values[indices[tuple(row[:ndim].astype(int))]] = row[-3]
            axes[0].plot(values, linewidth=.9, label=name, color=f"C{number % 10}")
        for number, (name, right) in enumerate(members[1:], start=1):
            _, rows = compare_products(left, right, key, relative_floor=config.relative_floor)
            x = [indices[tuple(row[:ndim].astype(int))] for row in rows]
            axes[1].plot(x, rows[:, -2], linewidth=.9, color=f"C{number % 10}")
            axes[2].plot(x, 100*rows[:, -1], linewidth=.9, color=f"C{number % 10}")
        axes[0].set(ylabel="xi" if "2pcf" in key else "zeta", title=f"{key}; reference: {base}")
        axes[1].set(ylabel="candidate - reference")
        axes[2].set(xlabel="lexicographically ordered histogram bin", ylabel="relative difference [%]")
        for ax in axes:
            ax.axhline(0, color="0.3", lw=.6)
            ax.grid(alpha=.2)
        axes[0].legend(fontsize=7, ncol=2)
        fig.text(.5, .01, f"Relative difference omitted where |reference| <= {config.relative_floor:g}; "
                 "these are numerical differences, not statistical uncertainties.", ha="center", fontsize=8)
        fig.tight_layout(rect=(0, .04, 1, 1))
        path = config.output_dir / f"{key}_differences.png"
        fig.savefig(path, dpi=140)
        plt.close(fig)
        paths.append(str(path))
    return comparisons, paths


def analysis_workspace_bytes(config):
    cells = config.rp_bins*config.rt_bins
    dense = 48*cells*cells if (config.covariance or config.distortion_matrix or config.reference_covariance) else 0
    return dense + 32*cells*config.wedge_bins*(len(config.wedge_mu_edges)-1)


def _validate_binning(metadata, config):
    expected = dict(rp_max=config.rp_max, rt_max=config.rt_max,
                    rp_bins=config.rp_bins, rt_bins=config.rt_bins)
    for key, value in expected.items():
        if key in metadata and not np.isclose(float(metadata[key]), value, rtol=1e-12, atol=0):
            raise ValueError(f"analysis input {key}={metadata[key]} does not match run value {value}")
    if "rp_min" in metadata and float(metadata["rp_min"]) != 0:
        raise ValueError("analysis input must have non-negative rp bins starting at zero")


def load_analysis_array(path, kind, config):
    """NPY, metadata-bearing NPZ, or lya2pcf/PICCA COR/DM FITS columns."""
    path = Path(path)
    size = config.rp_bins * config.rt_bins
    shape = (size, size) if kind in ("distortion", "covariance") else (size,)
    if math.prod(shape)*8 > config.max_hist_mib*2**20:
        raise ValueError(f"{kind} exceeds --max-hist-mib")
    metadata = {}
    fits_key = dict(distortion="DM", covariance="CO", model="DA")[kind]
    if path.suffix == ".npy":
        array = np.load(path, allow_pickle=False, mmap_mode="r")
    elif path.suffix == ".npz":
        with np.load(path, allow_pickle=False) as data:
            field = next((key for key in (kind, fits_key, "correlation" if kind == "model" else kind)
                          if key in data), None)
            if field is None:
                raise ValueError(f"{path}: missing {kind}/{fits_key} array")
            array = data[field]
            metadata = {key: data[key].item() for key in
                        ("rp_max", "rt_max", "rp_bins", "rt_bins", "rp_min") if key in data}
    else:
        from astropy.io import fits
        with fits.open(path, memmap=False) as hdus:
            hdu = next((h for h in hdus if isinstance(h, fits.BinTableHDU) and
                        fits_key in h.columns.names), None)
            if hdu is None:
                raise ValueError(f"{path}: no {fits_key} FITS column")
            array = np.array(hdu.data[fits_key])
            metadata = {name: hdu.header[key] for name, key in
                        (("rp_max", "RPMAX"), ("rt_max", "RTMAX"), ("rp_min", "RPMIN"),
                         ("rp_bins", "NP"), ("rt_bins", "NT")) if key in hdu.header}
    _validate_binning(metadata, config)
    if kind == "model" and array.shape == (config.rp_bins, config.rt_bins):
        array = array.ravel()
    if array.shape != shape:
        raise ValueError(f"{path}: {kind} shape {array.shape}; expected {shape} (rp-major, rt-minor)")
    if array.dtype.kind not in "biuf":
        raise ValueError(f"{kind} must contain real numeric values")
    array = np.asarray(array, dtype=float)
    if kind == "distortion":
        valid = np.isfinite(array).all(axis=1)
        if np.any(~valid & ~np.isnan(array).all(axis=1)):
            raise ValueError("distortion matrix contains partly non-finite rows")
    elif not np.isfinite(array).all():
        raise ValueError(f"{kind} contains non-finite values")
    if kind == "covariance":
        if not np.allclose(array, array.T, atol=1e-14, rtol=1e-8) or np.any(np.diag(array) < 0):
            raise ValueError("covariance must be symmetric with non-negative diagonal")
    return array


def apply_distortion(matrix, model):
    matrix, model = np.asarray(matrix), np.asarray(model)
    if matrix.ndim != 2 or model.ndim != 1 or matrix.shape[1] != len(model):
        raise ValueError("distortion/model shape mismatch")
    if not np.isfinite(model).all():
        raise ValueError("model must be finite, including unoccupied input bins")
    return matrix @ model


def wedge_geometry(config):
    """Area overlap in (rp,rt), sampled on a regular sub-bin grid."""
    mu_edges = np.array(config.wedge_mu_edges)
    rmax = config.wedge_r_max or min(config.rp_max, config.rt_max)
    radial_edges = np.linspace(0, rmax, config.wedge_bins+1)
    nr, nm = config.wedge_bins, len(mu_edges)-1
    geometry = np.zeros((nm*nr, config.rp_bins*config.rt_bins))
    offsets = (np.arange(config.wedge_subsamples)+.5)/config.wedge_subsamples
    for i in range(config.rp_bins):
        rp = (i+offsets)*config.rp_max/config.rp_bins
        for j in range(config.rt_bins):
            rt = (j+offsets)*config.rt_max/config.rt_bins
            radius = np.hypot(rp[:, None], rt[None, :])
            mu = rp[:, None]/radius
            rb = np.searchsorted(radial_edges, radius, side="right")-1
            mb = np.searchsorted(mu_edges, mu, side="right")-1
            valid = (rb >= 0) & (rb < nr) & (mb >= 0) & (mb < nm)
            geometry[:, i*config.rt_bins+j] = np.bincount(
                (mb[valid]*nr+rb[valid]).ravel(), minlength=nm*nr)/config.wedge_subsamples**2
    return geometry, .5*(radial_edges[:-1]+radial_edges[1:])


def project_wedges(xi, denominator, geometry, covariance=None):
    xi, denominator = np.asarray(xi).ravel(), np.asarray(denominator).ravel()
    active = (denominator > 0) & np.isfinite(xi)
    if covariance is None:
        weight = np.where(active, denominator, 0.)
    else:
        diagonal = np.diag(covariance)
        weight = np.divide(1., diagonal, out=np.zeros_like(diagonal), where=active & (diagonal > 0))
    projection = geometry*weight
    total = projection.sum(axis=1)
    projection = np.divide(projection, total[:, None], out=np.zeros_like(projection), where=total[:, None] > 0)
    value = projection @ np.nan_to_num(xi)
    value[total == 0] = np.nan
    covariance_out = None if covariance is None else projection @ covariance @ projection.T
    error = np.full(len(value), np.nan)
    if covariance_out is not None:
        variance = np.diag(covariance_out)
        if np.any(variance < -1e-10*np.max(np.abs(np.diag(covariance)), initial=0)):
            raise ValueError("covariance produces negative wedge variance; check the input matrix")
        error[total > 0] = np.sqrt(np.maximum(0, variance[total > 0]))
    return value, error, projection, covariance_out


def analyse_2pcf(results, config):
    members = [(name, row["products"]["3d_2pcf"]) for name, row in results.items()
               if "3d_2pcf" in row["products"]]
    if not members:
        return dict(status="no anisotropic 3D 2PCF selected"), []
    members.sort(key=lambda item: not item[0].startswith("lya2pcf-"))
    cells = config.rp_bins*config.rt_bins
    # Full matrices, validation copies, and projection workspaces, not just one input.
    if analysis_workspace_bytes(config) > config.max_hist_mib*2**20:
        raise ValueError("2PCF analysis workspace exceeds --max-hist-mib; reduce bins or increase budget")
    covariance, covariance_source = None, "unavailable: no covariance supplied/computed"
    if config.covariance:
        covariance = load_analysis_array(config.covariance, "covariance", config)
        covariance_source = str(config.covariance)
    elif "lya2pcf-cpu" in results and config.reference_covariance:
        covariance_source = results["lya2pcf-cpu"]["provenance"]["covariance"]
        with np.load(results["lya2pcf-cpu"]["reference_archive"], allow_pickle=False) as data:
            if "covariance" in data:
                covariance = data["covariance"]
                covariance_source = "shared lya2pcf weighted HEALPix subsampling covariance (unsmoothed)"
    model, distorted = None, None
    if config.model_correlation:
        model = load_analysis_array(config.model_correlation, "model", config)
        distorted = apply_distortion(load_analysis_array(config.distortion_matrix, "distortion", config), model)
    geometry, radii = wedge_geometry(config)
    analysis_dir = config.output_dir / "analysis_2pcf"
    analysis_dir.mkdir(exist_ok=True)
    report = dict(covariance=covariance_source,
                  covariance_scope="one shared covariance, not independently estimated for each native engine",
                  wedge_weighting="inverse covariance diagonal" if covariance is not None else "pair-weight sums",
                  wedge_mu_edges=list(config.wedge_mu_edges), wedge_subsamples=config.wedge_subsamples,
                  matrix_order="non-negative rp-major, rt-minor; square matching-bin matrices only",
                  matrix_metadata="NPY has no bin metadata; its order/limits must be supplied correctly by the user",
                  distortion="forward model D @ xi, never an inverse correction to measured xi",
                  model_residuals={}, archives={})
    curves, fields, paths = {}, [], []
    for name, table in members:
        num, den = np.zeros((config.rp_bins, config.rt_bins)), np.zeros((config.rp_bins, config.rt_bins))
        index = tuple(table[:, k].astype(int) for k in range(2))
        num[index], den[index] = table[:, -2], table[:, -1]
        xi = np.divide(num, den, out=np.full_like(num, np.nan), where=den > 0)
        value, error, projection, wedge_cov = project_wedges(xi, den, geometry, covariance)
        payload = dict(correlation=xi, numerator=num, denominator=den, wedge_r=radii,
                       wedge_mu_edges=config.wedge_mu_edges,
                       wedge_correlation=value.reshape(-1, len(radii)),
                       wedge_error=error.reshape(-1, len(radii)),
                       rp_max=config.rp_max, rt_max=config.rt_max,
                       rp_bins=config.rp_bins, rt_bins=config.rt_bins)
        if covariance is not None:
            payload.update(covariance=covariance, wedge_covariance=wedge_cov,
                           error=np.where(den > 0, np.sqrt(np.diag(covariance)).reshape(den.shape), np.nan))
        if distorted is not None:
            active = den.ravel() > 0
            if not np.isfinite(distorted[active]).all():
                raise ValueError("distortion matrix has undefined rows in measured occupied bins")
            model_wedge = projection @ np.nan_to_num(distorted)
            model_wedge[~np.isfinite(value)] = np.nan
            residual = xi.ravel()-distorted
            relative = np.divide(residual, distorted, out=np.full(cells, np.nan),
                                 where=active & (np.abs(distorted) > config.relative_floor))
            payload.update(model=model, distorted_model=distorted, model_residual=residual,
                           model_relative_residual=relative, model_wedge=model_wedge.reshape(-1, len(radii)))
            report["model_residuals"][name] = dict(
                max_abs=float(np.max(abs(residual[active]), initial=0)), occupied_bins=int(active.sum()))
        archive = analysis_dir / (name+".npz")
        np.savez_compressed(archive, **payload)
        report["archives"][name] = str(archive)
        curves[name] = payload
        fields.append(xi)
    if config.plots:
        paths.extend(plot_analysis(curves, fields, config, analysis_dir))
    return report, paths


def plot_analysis(curves, fields, config, output):
    plt = pyplot()
    names = list(curves)
    baseline = curves[names[0]]
    nr, nm = config.wedge_bins, len(config.wedge_mu_edges)-1
    fig, axes = plt.subplots(2, nm, figsize=(4.3*nm, 7), squeeze=False)
    radii = baseline["wedge_r"]
    for m in range(nm):
        base = baseline["wedge_correlation"][m]
        for number, (name, data) in enumerate(curves.items()):
            values = data["wedge_correlation"][m]
            color = f"C{number % 10}"
            axes[0, m].plot(radii, radii**2*values, label=name, lw=1, color=color)
            if name != names[0]:
                relative = np.divide(values-base, base, out=np.full(nr, np.nan),
                                     where=np.abs(base) > config.relative_floor)
                axes[1, m].plot(radii, relative*100, lw=1, color=color)
        if np.isfinite(baseline["wedge_error"][m]).any():
            axes[0, m].errorbar(radii, radii**2*base, yerr=radii**2*baseline["wedge_error"][m],
                                fmt="none", color="0.35", alpha=.5)
        if "model_wedge" in baseline:
            axes[0, m].plot(radii, radii**2*baseline["model_wedge"][m], "k--", label="D @ model")
        axes[0, m].set(title=f"{config.wedge_mu_edges[m]:g} <= mu < {config.wedge_mu_edges[m+1]:g}",
                       ylabel="r^2 xi [(Mpc/h)^2]")
        axes[1, m].set(xlabel="r [Mpc/h]", ylabel="difference to reference [%]")
        for ax in axes[:, m]:
            ax.axhline(0, color="0.3", lw=.6)
            ax.grid(alpha=.2)
            ax.set_xlim(0, config.wedge_r_max or min(config.rp_max, config.rt_max))
        if not any(np.isfinite(data["wedge_correlation"][m]).any() for data in curves.values()):
            axes[0, m].text(.5, .5, "No usable bins", transform=axes[0, m].transAxes, ha="center")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", fontsize=7, ncol=3)
    uncertainty = "error bars use shared covariance" if "covariance" in baseline else "no statistical error bars"
    fig.suptitle(f"3D 2PCF wedges; reference: {names[0]}; {uncertainty}")
    fig.tight_layout(rect=(0, .04+.025*math.ceil(len(labels)/3), 1, .96))
    path = output / "wedges.png"
    fig.savefig(path, dpi=140)
    plt.close(fig)
    paths = [str(path)]
    rp = (np.arange(config.rp_bins)+.5)*config.rp_max/config.rp_bins
    rt = (np.arange(config.rt_bins)+.5)*config.rt_max/config.rt_bins
    r2 = rp[:, None]**2 + rt[None, :]**2
    ncols = min(3, len(names))
    fig, axes = plt.subplots(math.ceil(len(names)/ncols), ncols,
                             figsize=(5*ncols, 4*math.ceil(len(names)/ncols)), squeeze=False)
    finite = np.concatenate([(v*r2)[np.isfinite(v)] for v in fields])
    limit = float(np.max(np.abs(finite), initial=1e-15))
    for ax, name, xi in zip(axes.flat, names, fields):
        plot = ax.imshow(r2*xi, origin="lower", aspect="auto", extent=(0, config.rt_max, 0, config.rp_max),
                         cmap="RdBu_r", vmin=-limit, vmax=limit)
        ax.set(xlabel="r transverse [Mpc/h]", ylabel="r parallel [Mpc/h]", title=name)
        ax.title.set_fontsize(9)
        fig.colorbar(plot, ax=ax, label="r^2 xi")
    for ax in list(axes.flat)[len(names):]:
        ax.set_visible(False)
    fig.tight_layout()
    path = output / "r2_correlation.png"
    fig.savefig(path, dpi=140)
    plt.close(fig)
    paths.append(str(path))
    if "distorted_model" in baseline:
        fig, axes = plt.subplots(2, 1, figsize=(11, 7), sharex=True)
        for name, data in curves.items():
            axes[0].plot(data["model_residual"], label=name, lw=.8)
            axes[1].plot(100*data["model_relative_residual"], lw=.8)
        axes[0].set(ylabel="measured xi - D @ model", title="Forward-model residuals (not inter-engine errors)")
        axes[1].set(xlabel="rp-major, rt-minor bin", ylabel="relative model residual [%]")
        axes[0].legend(fontsize=7, ncol=2)
        fig.tight_layout()
        path = output / "model_residuals.png"
        fig.savefig(path, dpi=140)
        plt.close(fig)
        paths.append(str(path))
    return paths
