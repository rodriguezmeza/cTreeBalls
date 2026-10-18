"""Keep anisotropic signed moments, finite-L bins and exact bins distinct."""
from pathlib import Path
import math
import numpy as np


def _grid_table(path, shape, columns):
    table = np.loadtxt(path, ndmin=2)
    if table.shape != (math.prod(shape), columns):
        raise ValueError(f"{path}: expected {math.prod(shape)} rows and {columns} columns")
    index = table[:, :len(shape)]
    if (not np.isfinite(index).all() or np.any(index != np.floor(index))
            or np.any(index < 0) or np.any(index >= np.asarray(shape))):
        raise ValueError(f"{path}: invalid histogram indices")
    flat = np.ravel_multi_index(index.astype(np.int64).T, shape)
    if len(np.unique(flat)) != len(flat):
        raise ValueError(f"{path}: duplicate histogram indices")
    if not np.array_equal(flat, np.arange(len(flat))):
        table = table[np.argsort(flat)]
    return table


def read_outputs(root, config):
    root = Path(root)
    shape = (config.r3_bins, config.r3_bins, config.theta_bins, config.theta_bins)
    raw_path = root/'histZetaM_lya_multipoles.txt'
    finite_path = root/'histZetaM_lya5d_multipole.txt'
    raw = _grid_table(raw_path, (*shape, config.lya3_lmax+1), 7)
    if not np.isfinite(raw).all():
        raise ValueError(f"{raw_path}: non-finite raw moments")
    # Higher Legendre denominator moments are signed. Only ell=0 is a
    # positive weight sum; a ratio of N_ell/D_ell is not a multipole estimator.
    if np.any(raw[raw[:, 4] == 0, -1] < 0):
        raise ValueError(f"{raw_path}: negative monopole denominator")
    finite = _grid_table(finite_path, (*shape, config.mu_bins), 13)
    if not np.isfinite(finite[:, :10]).all() or not np.isfinite(finite[:, -2:]).all():
        raise ValueError(f"{finite_path}: non-finite geometry or raw reconstructed sums")
    positive = finite[:, -1] > 0
    expected = finite[positive, -2]/finite[positive, -1]
    if (not np.isfinite(finite[positive, -3]).all()
            or not np.allclose(finite[positive, -3], expected, rtol=5e-12, atol=1e-14)
            or not np.isnan(finite[~positive, -3]).all()):
        raise ValueError(f"{finite_path}: invalid reconstruction ratio or NaN convention")
    index = finite[:, :5]
    centers = np.column_stack(((index[:, :2]+.5)*config.r3_max/config.r3_bins,
                              (index[:, 2:4]+.5)*np.pi/config.theta_bins,
                              -1+(index[:, 4]+.5)*2/config.mu_bins))
    if not np.allclose(finite[:, 5:10], centers, rtol=1e-12, atol=1e-13):
        raise ValueError(f"{finite_path}: bin coordinates disagree with the requested grid")
    return raw, finite


def reconstruction_diagnostics(finite, exact, config):
    shape = (config.r3_bins, config.r3_bins, config.theta_bins, config.theta_bins, config.mu_bins)
    num, den = np.zeros(len(finite)), np.zeros(len(finite))
    if len(exact):
        bins = np.ravel_multi_index(exact[:, :5].astype(np.int64).T, shape)
        num[bins], den[bins] = exact[:, -2], exact[:, -1]
    value = np.divide(num, den, out=np.zeros_like(num), where=den > 0)
    eligible = (den > 0) & (np.abs(value) > config.relative_floor)
    valid = np.isfinite(finite[:, -3]) & (finite[:, -1] > 0)
    selected = eligible & valid
    errors = np.abs((finite[selected, -3]-value[selected])/value[selected])

    def relative_l2(a, b):
        scale = float(np.max(np.abs(a), initial=0))
        norm = np.linalg.norm(a/scale) if scale else 0.
        return float(np.linalg.norm((b-a)/scale)/norm) if norm else None

    return dict(reference="exact hard-bin companion from this same run",
                relative_floor=config.relative_floor,
                eligible_bins=int(eligible.sum()),
                excluded_low_signal_bins=int(((den > 0) & ~eligible).sum()),
                invalid_eligible_bins=int((eligible & ~valid).sum()),
                negative_reconstructed_denominator_bins=int((finite[:, -1] < 0).sum()),
                numerator_relative_l2=relative_l2(num, finite[:, -2]),
                denominator_relative_l2=relative_l2(den, finite[:, -1]),
                max_relative_correlation_error=float(errors.max()) if len(errors) else None,
                p95_relative_correlation_error=float(np.percentile(errors, 95)) if len(errors) else None,
                empty_bin_denominator_leakage_l1=float(np.abs(finite[den == 0, -1]).sum()),
                scope="Approximate reconstruction diagnostics only; exact companion acceptance "
                      "does not certify finite-L bins, raw moments, or a different catalog")


def inspect_outputs(root, config, exact=None, plot_root=None):
    raw, finite = read_outputs(root, config)
    mode = config.lya3_mu_mode
    if (mode == "exact") != (exact is not None):
        raise ValueError("anisotropic exact-bin availability disagrees with requested mu mode")
    root = Path(root)
    report = dict(mode=mode, lmax=config.lya3_lmax, multipole_kernel=config.lya3_multipole_kernel,
        moments=dict(file=str(root/'histZetaM_lya_multipoles.txt'), rows=len(raw),
                     definition="Exact signed Legendre numerator and denominator moments; no ratio"),
        finite_reconstruction=dict(file=str(root/'histZetaM_lya5d_multipole.txt'), rows=len(finite),
                     approximate=True, negative_denominator_bins=int((finite[:, -1] < 0).sum()),
                     zero_denominator_bins=int((finite[:, -1] == 0).sum()),
                     undefined_ratio_bins=int(np.isnan(finite[:, -3]).sum()),
                     included_in_exact_comparisons=False),
        exact_bins=dict(available=exact is not None,
                        file=str(root/'histZetaM_lya5d.txt') if exact is not None else None),
        plots=[])
    if exact is not None:
        report['finite_reconstruction_vs_exact'] = reconstruction_diagnostics(finite, exact, config)
    if config.plots:
        destination = Path(plot_root) if plot_root is not None else root
        destination.mkdir(parents=True, exist_ok=True)
        report['plots'] = plot_outputs(destination, raw, finite, config)
        from lya_multipole_plots import write_anisotropic_plots
        paths, details = write_anisotropic_plots(destination, raw, config, exact)
        report['plots'].extend(paths)
        report['multipole_plotting'] = details
    return report


def plot_outputs(root, raw, finite, config):
    from lya_analysis import pyplot
    plt = pyplot()
    lmax = config.lya3_lmax
    orders = sorted({ell for ell in (0, 1, 2, lmax) if ell <= lmax})
    # Sum signed raw moments over the polar bins. No normalization by a
    # potentially negative/zero denominator moment is performed.
    fields = raw[:, -2:].reshape(config.r3_bins, config.r3_bins,
                                 config.theta_bins, config.theta_bins, lmax+1, 2).sum(axis=(2, 3))
    fig, axes = plt.subplots(len(orders), 2, figsize=(10, 3*len(orders)), squeeze=False)
    for row, ell in enumerate(orders):
        for component in range(2):
            value = fields[:, :, ell, component]
            limit = max(float(np.max(np.abs(value), initial=0)), 1e-15)
            ax = axes[row, component]
            artist = ax.imshow(value, origin='lower', cmap='RdBu_r', vmin=-limit, vmax=limit,
                               extent=(0, config.r3_max, 0, config.r3_max), aspect='auto')
            name = 'Numerator' if component == 0 else 'Denominator'
            ax.set(title=f'{name} moment, ell={ell}', xlabel='r2 [Mpc/h]', ylabel='r1 [Mpc/h]')
            fig.colorbar(artist, ax=ax)
    fig.suptitle('Anisotropic signed Legendre moments: sums over both polar bins; no ratios', fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, .96))
    path = root/'anisotropic_raw_moments.png'
    fig.savefig(path, dpi=140, bbox_inches='tight')
    plt.close(fig)
    paths = [str(path)]

    fig, axes = plt.subplots(3, 1, figsize=(11, 8), sharex=True)
    sample = np.unique(np.linspace(0, len(finite)-1, min(6000, len(finite))).astype(int))
    for ax, column, label in zip(axes, (-3, -2, -1), ('Reconstructed zeta', 'Signed numerator', 'Signed denominator')):
        ax.plot(sample, finite[sample, column], lw=.7)
        ax.set_ylabel(label)
        ax.axhline(0, color='.4', lw=.6)
        ax.grid(alpha=.2)
    axes[-1].set_xlabel('Lexicographic native bin index (plot samples at most 6000 bins)')
    fig.suptitle(f'APPROXIMATE finite-L reconstruction, L={lmax}; undefined ratios remain NaN', fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, .95))
    path = root/'anisotropic_finite_reconstruction.png'
    fig.savefig(path, dpi=140, bbox_inches='tight')
    plt.close(fig)
    paths.append(str(path))
    return paths
