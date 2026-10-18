"""Radial views of native moments and binned angular correlation projections."""
from pathlib import Path
import math

import numpy as np
from numpy.polynomial.legendre import Legendre


def selected_orders(config, available):
    requested = getattr(config, 'plot_multipole_orders', None)
    if requested is None:
        return list(available)
    if (not requested or any(not isinstance(x, int) or not 0 <= x <= 32 for x in requested)
            or len(set(requested)) != len(requested)):
        raise ValueError('plot_multipole_orders must contain distinct integers in [0,32]')
    return [ell for ell in sorted(requested) if ell in available]


def anisotropic_radial_moments(raw, config):
    """Sum polar cells before normalizing by the positive monopole weight."""
    fields = raw[:, -2:].reshape(config.r3_bins, config.r3_bins,
        config.theta_bins, config.theta_bins, config.lya3_lmax+1, 2).sum(axis=(2, 3))
    num, den = fields[..., 0], fields[..., 1]
    d0 = den[..., :1]
    factor = 2*np.arange(config.lya3_lmax+1)+1
    q = np.divide(num*factor, d0, out=np.full_like(num, np.nan), where=d0 > 0)
    window = np.divide(den*factor, d0, out=np.full_like(den, np.nan), where=d0 > 0)
    return num, den, q, window


def binned_correlation_multipoles(exact, config):
    """Integrate the piecewise-constant N(mu)/D(mu), requiring full support.

    These are projections of the binned correlation, not unbinned moments or
    a deconvolution of the survey window. Empty angular bins stay undefined.
    """
    shape = (config.r3_bins, config.r3_bins, config.mu_bins)
    num, den = np.zeros(shape), np.zeros(shape)
    if len(exact):
        indices = tuple(exact[:, i].astype(int) for i in (0, 1, 4))
        np.add.at(num, indices, exact[:, -2])
        np.add.at(den, indices, exact[:, -1])
    supported = den > 0
    correlation = np.divide(num, den, out=np.full(shape, np.nan), where=supported)
    edges = np.linspace(-1., 1., config.mu_bins+1)
    coverage = np.sum(supported*np.diff(edges), axis=-1)/2
    # Higher orders from a coarse angular grid would imply unsupported
    # angular resolution. Even these retained projections have binning error.
    orders = np.arange(min(config.lya3_lmax, config.mu_bins-1)+1)
    weights = np.array([(2*ell+1)/2*np.diff(Legendre.basis(int(ell)).integ()(edges))
                        for ell in orders])
    coefficients = np.einsum('ijm,lm->ijl', correlation, weights)
    return coefficients, coverage


def radial_views(root, stem, values, centers, orders, title, label, bounds=None):
    """Paginate all requested orders into radial maps and fixed-r1 slices."""
    from lya_analysis import pyplot
    plt = pyplot()
    root = Path(root)
    root.mkdir(parents=True, exist_ok=True)
    paths = []
    if not orders:
        return paths
    step = centers[1]-centers[0] if len(centers) > 1 else 2*centers[0]
    lo, hi = bounds if bounds is not None else (centers[0]-step/2, centers[-1]+step/2)
    for start in range(0, len(orders), 9):
        page = orders[start:start+9]
        columns = min(3, len(page))
        rows = math.ceil(len(page)/columns)
        suffix = f'_page{start//9+1}' if len(orders) > 9 else ''
        for kind in ('maps', 'slices'):
            fig, axes = plt.subplots(rows, columns, figsize=(5*columns, 3.9*rows), squeeze=False)
            for ax, ell in zip(axes.flat, page):
                field = values[..., ell]
                if kind == 'maps':
                    finite = field[np.isfinite(field)]
                    limit = max(float(np.max(np.abs(finite), initial=0)), 1e-15)
                    cmap = plt.get_cmap('RdBu_r').with_extremes(bad='0.85')
                    artist = ax.imshow(np.ma.masked_invalid(field), origin='lower', cmap=cmap,
                        vmin=-limit, vmax=limit, extent=(lo, hi, lo, hi), aspect='equal')
                    fig.colorbar(artist, ax=ax, label=label)
                    ax.set(xlabel='r2 [Mpc/h]', ylabel='r1 [Mpc/h]')
                else:
                    occupied = np.flatnonzero(np.isfinite(field).any(axis=1))
                    picks = np.unique(np.linspace(0, len(occupied)-1, min(3, len(occupied))).astype(int))
                    for i in occupied[picks]:
                        ax.plot(centers, field[i], marker='.', label=f'r1={centers[i]:.3g}')
                    diagonal = np.diag(field)
                    if np.isfinite(diagonal).any():
                        ax.plot(centers, diagonal, 'k--', lw=1, label='r1=r2')
                    if ax.lines:
                        ax.legend(fontsize=8)
                    else:
                        ax.text(.5, .5, 'No supported bins', ha='center', transform=ax.transAxes)
                    ax.axhline(0, color='.5', lw=.6)
                    ax.grid(alpha=.25)
                    ax.set(xlabel='r2 [Mpc/h]', ylabel=label)
                ax.set_title(f'ell={ell}')
            for ax in list(axes.flat)[len(page):]:
                ax.set_visible(False)
            fig.suptitle(title, fontsize=11)
            fig.tight_layout(rect=(0, 0, 1, .94))
            path = root/f'{stem}_{kind}{suffix}.png'
            fig.savefig(path, dpi=140, bbox_inches='tight')
            plt.close(fig)
            paths.append(str(path))
    return paths


def write_anisotropic_plots(root, raw, config, exact=None):
    root = Path(root)
    root.mkdir(parents=True, exist_ok=True)
    num, den, q, window = anisotropic_radial_moments(raw, config)
    centers = (np.arange(config.r3_bins)+.5)*config.r3_max/config.r3_bins
    orders = selected_orders(config, range(config.lya3_lmax+1))
    if not orders:
        raise ValueError('none of --plot-multipole-orders are available in the anisotropic output')
    indices = np.indices(num.shape).reshape(3, -1).T
    d0 = np.broadcast_to(den[..., :1], den.shape)
    csv = root/'anisotropic_multipoles.csv'
    np.savetxt(csv, np.column_stack((indices, centers[indices[:, 0]], centers[indices[:, 1]],
        num.ravel(), den.ravel(), d0.ravel(), q.ravel(), window.ravel())), delimiter=',',
        header='b1,b2,ell,r1_mpc_h,r2_mpc_h,N_ell,D_ell,D_0,Q_ell,window_ell', comments='')
    paths = radial_views(root, 'anisotropic_multipoles', q, centers, orders,
        'Window-weighted Legendre moments: Q_l=(2l+1) N_l / D_0\n'
        'Polar bins summed first; survey-window coupling remains; gray bins undefined', 'Q_l')
    paths += radial_views(root, 'anisotropic_window_multipoles', window, centers, orders,
        'Angular sampling moments: W_l=(2l+1) D_l / D_0\n'
        'Polar bins summed first; gray bins undefined', 'W_l')
    report = dict(csv=str(csv), plotted_orders=orders,
        definition='Q_l=(2l+1) sum_polar N_l / sum_polar D_0; window_l=(2l+1) sum_polar D_l / sum_polar D_0',
        window_corrected=False, undefined_radial_cells=int(np.isnan(q[..., 0]).sum()),
        scope='Exact sampled weighted moments; Q_l is not a survey-window-corrected correlation multipole. No N_l/D_l ratio is used.')
    if exact is not None:
        binned, coverage = binned_correlation_multipoles(exact, config)
        binned_orders = selected_orders(config, range(binned.shape[-1]))
        indices = np.indices(binned.shape).reshape(3, -1).T
        binned_csv = root/'anisotropic_binned_multipoles.csv'
        np.savetxt(binned_csv, np.column_stack((indices, centers[indices[:, 0]], centers[indices[:, 1]],
            binned.ravel(), coverage[indices[:, 0], indices[:, 1]])), delimiter=',',
            header='b1,b2,ell,r1_mpc_h,r2_mpc_h,zeta_ell_binned,mu_coverage', comments='')
        paths += radial_views(root, 'anisotropic_binned_multipoles', binned, centers, binned_orders,
            'Correlation multipoles projected from exact angular bins\n'
            'Finite mu-bin resolution; full mu coverage required; gray bins undefined', 'zeta_l (binned)')
        report['binned_correlation'] = dict(csv=str(binned_csv), plotted_orders=binned_orders,
            definition='(2l+1)/2 sum_mu [sum_polar N_mu / sum_polar D_mu] integral_bin P_l(mu) dmu',
            complete_radial_cells=int(np.all(np.isfinite(binned), axis=-1).sum()),
            orders_available=list(range(binned.shape[-1])),
            scope='Piecewise-constant projection of exact-bin correlation; finite mu resolution, not unbinned coefficients or finite-L reconstruction. Missing angular support gives NaN.')
    return paths, report


def write_scalar_plots(root, table, config):
    """Use the scalar engine's native normalization without an extra factor."""
    nb = config.r3_bins
    values = np.full((nb, nb, config.multipole_lmax+1), np.nan)
    table = table[table[:, -1] > 0]
    if len(table):
        indices = (table[:, 1].astype(int)-1, table[:, 2].astype(int)-1, table[:, 0].astype(int))
        values[indices] = table[:, -3]
    orders = selected_orders(config, range(config.multipole_lmax+1))
    if not orders:
        raise ValueError('none of --plot-multipole-orders are available in the scalar output')
    centers = config.multipole_rmin+(np.arange(nb)+.5)*(config.r3_max-config.multipole_rmin)/nb
    # Native centers are authoritative; the driver requests linear bins.
    for row in table:
        centers[int(row[1])-1] = row[3]
        centers[int(row[2])-1] = row[4]
    return radial_views(root, 'scalar_multipoles', values, centers, orders,
        'Scalar 3PCF multipoles: native numerator / denominator\n'
        'Only independent native shell pairs are shown; gray bins undefined', 'zeta_l (native)',
        bounds=(config.multipole_rmin, config.r3_max))
