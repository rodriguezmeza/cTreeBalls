#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Acknowledgements: Axel Romero Tisnado
#                   Python script by Tisnado
#                   adapted to use cballys module
#
#
import os, argparse, numpy as np
import matplotlib as mpl; mpl.use("Agg")
import matplotlib.pyplot as plt

# Aesthetic (Times/STIX)
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman","Nimbus Roman","TeX Gyre Termes","DejaVu Serif"],
    "mathtext.fontset": "stix",
    "axes.titlesize": 18, "axes.labelsize": 16,
    "xtick.labelsize": 12, "ytick.labelsize": 12,
    "legend.fontsize": 12,
})

def load_xy(path):
    arr = np.loadtxt(path)
    if arr.ndim == 1:
        arr = arr.reshape(1, -1)
    if arr.shape[1] < 2:
        raise ValueError(f"{path}: necesito ≥2 columnas (theta, y)")
    th = arr[:, 0].astype(float)
    y  = arr[:, 1].astype(float)
    return th, y

def crop_theta(th, y, tmin=None, tmax=None):
    m = np.isfinite(th) & np.isfinite(y)
    if tmin is not None: m &= (th >= tmin)
    if tmax is not None: m &= (th <= tmax)
    return th[m], y[m]

def interp_to_ref(th_src, y_src, th_ref):
    
    y_i = np.interp(th_ref, th_src, y_src, left=np.nan, right=np.nan)
    return y_i

def error_stats(y, yref, eps=1e-30):
    
    valid = np.isfinite(y) & np.isfinite(yref)
    if not np.any(valid):
        return np.nan, np.nan, np.nan
    dy = y[valid] - yref[valid]
    denom = np.maximum(np.abs(yref[valid]), eps)

    mae = 100.0 * np.mean(np.abs(dy) / denom)
    rmse = 100.0 * np.sqrt(np.mean((dy / denom) ** 2))
    l1 = 100.0 * (np.sum(np.abs(dy)) / max(np.sum(np.abs(yref[valid])), eps))
    return mae, rmse, l1

def residual_percent(y, yref, eps=1e-30):
    valid = np.isfinite(y) & np.isfinite(yref)
    out = np.full_like(yref, np.nan, dtype=float)
    out[valid] = 100.0 * (y[valid] - yref[valid]) / np.maximum(np.abs(yref[valid]), eps)
    return out

def main():
    ap = argparse.ArgumentParser(description="Compare 2 curves (convergence Xi2pcf) with residuals and metrics.")
    ap.add_argument("--file-a", required=True, help="TXT A (p.ej. DES TreeCorr)")
    ap.add_argument("--file-b", required=True, help="TXT B (p.ej. DES Corrfunc)")
    ap.add_argument("--label-a", default="cBalls_a")
    ap.add_argument("--label-b", default="cBalls_b")
    ap.add_argument("--ref", choices=["a","b"], default="a",
        help="Reference curve for errors/residuals")
    ap.add_argument("--theta-min", type=float, default=None)
    ap.add_argument("--theta-max", type=float, default=None)
    ap.add_argument("--scale", choices=["loglog","linear","semilogx","semilogy"], default="loglog")
    ap.add_argument("--outdir", default=None)
    ap.add_argument("--out", default=None,
        help="Base name for the output (no extension)")
    # In case some TXT already have Y=θ*w:
    ap.add_argument("--y-is-thetax-a", action="store_true")
    ap.add_argument("--y-is-thetax-b", action="store_true")
    # If you want to multiply by θ in the final plot
    ap.add_argument("--plot-mul-theta", action="store_true",
                    help="Multiply Y by θ (only for main figure; do not affect metrics if already are in θ*w)")
    ap.add_argument("--xscale", choices=["radian","arcmin","degree"],
                    default="radians",
                    help="x units: arcmin (default), radian or degree.")
    args = ap.parse_args()

    #B Load:
    thA, yA = load_xy(args.file-a if False else args.file_a)
    thB, yB = load_xy(args.file_b)
    #E

    #B Cut in θ:
    thA, yA = crop_theta(thA, yA, args.theta_min, args.theta_max)
    thB, yB = crop_theta(thB, yB, args.theta_min, args.theta_max)

    if args.y_is_thetax_a:
        yA = np.where(thA != 0, yA / thA, np.nan)
    if args.y_is_thetax_b:
        yB = np.where(thB != 0, yB / thB, np.nan)
    #E

    #B Choose reference:
    if args.ref == "a":
        th_ref, y_ref, lab_ref = thA, yA, args.label_a
        thX, yX, labX = thB, yB, args.label_b
    else:
        th_ref, y_ref, lab_ref = thB, yB, args.label_b
        thX, yX, labX = thA, yA, args.label_a
    #E

    # Interpole B to the grid of reference for metrics and residuals
    yX_i = interp_to_ref(thX, yX, th_ref)

    # Metrics:
    mae_X, rmse_X, l1_X = error_stats(yX_i, y_ref)

    # Residuals (%):
    resX = residual_percent(yX_i, y_ref)

    #B Figure: upper panel superimposed; lower panel, residuals (%)
    fig = plt.figure(figsize=(7.6, 7.2), dpi=140)
    gs = fig.add_gridspec(nrows=2, ncols=1, height_ratios=[2.6, 1.4], hspace=0.06)
    ax = fig.add_subplot(gs[0, 0])
    axr = fig.add_subplot(gs[1, 0], sharex=ax)

    ax.xaxis.label.set_visible(False)

    plotA = thA*yA if args.plot_mul_theta else yA
    plotB = thB*yB if args.plot_mul_theta else yB
    plot_ref = th_ref*y_ref if args.plot_mul_theta else y_ref

    #B Curves:
    ax.plot(thA, plotA, lw=1.9, label=args.label_a)
    ax.plot(thB, plotB, lw=1.9, label=args.label_b)
    #E

    #B Scales:
    if args.scale == "loglog":
        ax.set_xscale("log"); ax.set_yscale("log")
    elif args.scale == "semilogx":
        ax.set_xscale("log")
    elif args.scale == "semilogy":
        ax.set_yscale("log")
    #E

    ax.tick_params(axis='x', which='both', labelbottom=False, bottom=False)

    ax.set_ylabel(r"$\theta\, \zeta_{\kappa\kappa}(\theta)$" if args.plot_mul_theta else r"$w_{\kappa\kappa}(\theta)$")
    ax.grid(alpha=0.25, ls=":")
    ax.legend(frameon=True, ncols=1, loc="best")
    ax.set_title(f"Comparison between different outputs", pad=6)

    txtX = f"{labX}: MAE={mae_X:.3f}%  RMSE={rmse_X:.3f}%  L1={l1_X:.3f}%"

    # Residuals (%), on reference's grid:
    axr.axhline(0.0, color="k", lw=1.0, alpha=0.5)
    axr.plot(th_ref, resX, lw=1.6, label=f"{labX} - {lab_ref}")

    #  axr.set_xlabel(r"$\theta$ [deg]")
    if args.xscale == "degree":
        axr.set_xlabel(r"$\theta\ [{\rm deg}]$")
    else:
        if args.xscale == "arcmin":
            axr.set_xlabel(r"$\theta\ [{\rm arcmin}]$")
        else:
            axr.set_xlabel(r"$\theta$ [radians]")

    axr.set_ylabel(r"$\Delta(\%)$")
    axr.grid(alpha=0.25, ls=":")
    if args.scale in ("loglog","semilogx"):
        axr.set_xscale("log")

    #B Output:
    outdir = args.outdir or os.path.dirname(os.path.abspath(args.file_a))
    os.makedirs(outdir, exist_ok=True)
    base = args.out or f"compare_vs_{ {'a':args.label_a,'b':args.label_b}[args.ref] }"
    fig.tight_layout()
    pdf = os.path.join(outdir, f"{base}.pdf")
    fig.savefig(pdf)
    plt.close(fig)
    #E

    #B TXT with metrics:
    with open(os.path.join(outdir, f"{base}_metrics.txt"), "w") as f:
        f.write(f"Reference: {lab_ref}\n")
        f.write(f"{txtX}\n")
    #E

    print("[OK] Figura:", pdf)
    print("[OK] Métricas:")
    print("   ", txtX)

if __name__ == "__main__":
    main()
