#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Acknowledgements: Axel Romero Tisnado
#                   Python script by Tisnado
#                   adapted to use cballys module
#
import os, argparse, numpy as np
import matplotlib as mpl; mpl.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

#  (Times/STIX)
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman","Nimbus Roman","TeX Gyre Termes","DejaVu Serif"],
    "mathtext.fontset": "stix",
    "axes.titlesize": 18, "axes.labelsize": 16,
    "xtick.labelsize": 12, "ytick.labelsize": 12,
    "legend.fontsize": 12,
})

def load_z_matrix(path,type=None,nmonopoles=8):
    if type=="sincos":
        yfileCos = path+"/"+"histZetaM_cos_"
        yfileSin = path+"/"+"histZetaM_sin_"
    else:
        if type=="edge_effects":
            yfileEE = path+"/"+"histZetaM_EE_"
    expData = []
    if type=="sincos":
        for j in range(1,nmonopoles+1):
            expData.append(np.loadtxt(yfileCos + str(j) +".txt")
                    + np.loadtxt(yfileSin + str(j) +".txt"))
    else:
        if type=="edge_effects":
            for j in range(1,nmonopoles+1):
                expData.append(np.loadtxt(yfileEE + str(j) +".txt"))

    return expData

def crop_bins(y, tmin=None, tmax=None):
    m = np.isfinite(y)
    if tmin is not None: m &= (th >= tmin)
    if tmax is not None: m &= (th <= tmax)
    return y[m]

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
    ap = argparse.ArgumentParser(description="Compare 2 set of curves (convergence monopoles of Xi3pcf) with residuals and metrics.")
    ap.add_argument("--file-a", required=True, help="TXT A (p.ej. DES cBalls)")
    ap.add_argument("--file-b", default=None, help="TXT B (p.ej. DES cBalls)")
    ap.add_argument("--label-a", default="cBalls_a")
    ap.add_argument("--label-b", default="cBalls_b")
    ap.add_argument("--ref", choices=["a","b"], default="a",
        help="Reference curves for errors/residuals")
    ap.add_argument("--type", choices=["sincos","edge_effects"], default="sincos")
    ap.add_argument("--monopoles", type=int, default=8)
    ap.add_argument("--bin-min", type=float, default=None)
    ap.add_argument("--bin-max", type=float, default=None)
    ap.add_argument("--scale",
        choices=["loglog","linear","semilogx","semilogy"], default="linear")
    ap.add_argument("--outdir", default="./")
    ap.add_argument("--out", default="Xi3pcf", help="Output base name (no extension)")
    # In case some TXT already have Y=θ*w:
    ap.add_argument("--y-is-thetax-a", action="store_true")
    ap.add_argument("--y-is-thetax-b", action="store_true")
    # If you want to multiply by θ in the final plot
    ap.add_argument("--plot-mul-theta", action="store_true",
                    help="Multiplica la Y mostrada por θ (solo para la figura; no afecta métricas si ya están en θ*w)")
    ap.add_argument("--xscale", choices=["radian","arcmin","degree","bins"],
                    default="bins",
                    help="x units: number of bin (default), arcmin, radian or degree.")
    args = ap.parse_args()

    if args.file_b == None:
        file_b = args.file-a if False else args.file_a
    else:
        file_b = args.file_b

    #B Output:
    outdir = args.outdir or os.path.dirname(os.path.abspath(args.file_a))
    os.makedirs(outdir, exist_ok=True)
    base = args.out or f"compare_vs_{ {'a':args.label_a,'b':args.label_b}[args.ref] }"
    ofilename = os.path.join(outdir, f"{base}_flatten.pdf")
    ofilename_err = os.path.join(outdir, f"{base}_flatten_err.pdf")
    #E

    xlabel = r"triangle index"
    ylabel0 = r"$ \zeta_0$"
    ylabel1 = r"$ \zeta_1$"
    ylabel2 = r"$ \zeta_2$"
    ylabel3 = r"$ \zeta_3$"
    ylabel4 = r"$ \zeta_4$"
    dylabel0 = r"$\Delta \zeta_0 (\%)$"
    dylabel1 = r"$\Delta \zeta_1 (\%)$"
    dylabel2 = r"$\Delta \zeta_2 (\%)$"
    dylabel3 = r"$\Delta \zeta_3 (\%)$"
    dylabel4 = r"$\Delta \zeta_4 (\%)$"
    label1 = args.label_a
    label2 = args.label_b

    lw = 1.75
    color0 = "blue"
    color1 = "red"
    linestyle1 = "-"
    linestyle2 = "dashed"

    #B Load:
    yA = load_z_matrix(args.file-a if False else args.file_a, args.type, args.monopoles)
    yB = load_z_matrix(file_b, args.type, args.monopoles)
    #E

    print('triangle indexes yA:',yA[0].shape[0]*yA[0].shape[0])
    print('triangle indexes yB:',yB[0].shape[0]*yB[0].shape[0])
    xmin = args.bin_min
    xmax = args.bin_max

    #B Scales
    xscale="linear"
    yscale="linear"
    if args.scale == "loglog":
        xscale="log"; yscale="log"
    elif args.scale == "semilogx":
        xscale="log"
    elif args.scale == "semilogy":
        yscale="log"
    #E

    #B start making plots:
    fig = plt.figure(tight_layout=True,figsize=(15, 15))
    gs = gridspec.GridSpec(5, 3)
    formatter=mpl.ticker.ScalarFormatter(useMathText=True)
    formatter.set_scientific(True)

    # m=0
    ax = fig.add_subplot(gs[0, :])
    ax.plot(range(len(yA[0].flatten())),yA[0].flatten(),
        color=color0,linestyle=linestyle1,linewidth=lw,alpha=1,label=label1)
    ax.plot(range(len(yB[0].flatten())),yB[0].flatten(),
        color=color1,linestyle=linestyle2,linewidth=lw,alpha=1,label=label2)
    ax.set_ylabel(ylabel0,fontsize=35)
    ax.set_xlim(xmin,xmax)
    ax.tick_params(labelsize = 25)
    ax.legend(fontsize=25)
    #ax.set_ylim(ylim1,ylim2)
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)
    ax.axes.xaxis.set_ticklabels([])

    # m=1
    ax = fig.add_subplot(gs[1, :])
    ax.plot(range(len(yA[1].flatten())),yA[1].flatten(),
        color=color0,linestyle=linestyle1,linewidth=lw,alpha=1)
    ax.plot(range(len(yB[1].flatten())),yB[1].flatten(),
        color=color1,linestyle=linestyle2,linewidth=lw,alpha=1)
    ax.set_ylabel(ylabel1,fontsize=35)
    ax.set_xlim(xmin,xmax)
    ax.tick_params(labelsize = 25)
    #ax.legend(fontsize=25)
    #ax.set_ylim(ylim1,ylim2)
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)
    ax.axes.xaxis.set_ticklabels([])
    # m=2
    ax = fig.add_subplot(gs[2, :])
    ax.plot(range(len(yA[2].flatten())),yA[2].flatten(),
        color=color0,linestyle=linestyle1,linewidth=lw,alpha=1)
    ax.plot(range(len(yB[2].flatten())),yB[2].flatten(),
        color=color1,linestyle=linestyle2,linewidth=lw,alpha=1)
    ax.set_ylabel(ylabel2,fontsize=35)
    ax.set_xlim(xmin,xmax)
    ax.tick_params(labelsize = 25)
    #ax.legend(fontsize=25)
    #ax.set_ylim(ylim1,ylim2)
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)
    ax.axes.xaxis.set_ticklabels([])
    # m=3
    ax = fig.add_subplot(gs[3, :])
    ax.plot(range(len(yA[3].flatten())),yA[3].flatten(),
        color=color0,linestyle=linestyle1,linewidth=lw,alpha=1)
    ax.plot(range(len(yB[3].flatten())),yB[3].flatten(),
        color=color1,linestyle=linestyle2,linewidth=lw,alpha=1)
    ax.set_ylabel(ylabel3,fontsize=35)
    ax.set_xlim(xmin,xmax)
    ax.tick_params(labelsize = 25)
    #ax.legend(fontsize=25)
    #ax.set_ylim(ylim1,ylim2)
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)
    ax.axes.xaxis.set_ticklabels([])
    # m=4
    ax = fig.add_subplot(gs[4, :])
    ax.plot(range(len(yA[4].flatten())),yA[4].flatten(),
        color=color0,linestyle=linestyle1,linewidth=lw,alpha=1)
    ax.plot(range(len(yB[4].flatten())),yB[4].flatten(),
        color=color1,linestyle=linestyle2,linewidth=lw,alpha=1)
    ax.set_xlabel(xlabel,fontsize=30)
    ax.set_ylabel(ylabel4,fontsize=35)
    ax.set_xlim(xmin,xmax)
    ax.tick_params(labelsize = 25)
    #ax.legend(fontsize=25)
    #ax.set_ylim(ylim1,ylim2)
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)

    fig.align_labels()
    plt.savefig(ofilename,dpi=300)

    # plotting 3pcf err rel
    fig = plt.figure(tight_layout=True,figsize=(15, 15))
    gs = gridspec.GridSpec(5, 3)
    formatter=mpl.ticker.ScalarFormatter(useMathText=True)
    formatter.set_scientific(True)

    # m=0
    ax = fig.add_subplot(gs[0, :])
    ax.plot(range(len(yA[0].flatten())),(yA[0].flatten()-yB[0].flatten())*100/yB[0].flatten(),
            color=color0,linestyle=linestyle1,linewidth=lw,alpha=1,label=label1)
    ax.set_ylabel(dylabel0,fontsize=35)
    ax.set_xlim(xmin,xmax)
    ax.tick_params(labelsize = 25)
    #ax.legend(fontsize=25)
    ax.axes.xaxis.set_ticklabels([])
    # m=1
    ax = fig.add_subplot(gs[1, :])
    ax.plot(range(len(yA[1].flatten())),(yA[1].flatten()-yB[1].flatten())*100/yB[1].flatten(),
            color=color0,linestyle=linestyle1,linewidth=lw,alpha=1)
    ax.set_ylabel(dylabel1,fontsize=35)
    ax.set_xlim(xmin,xmax)
    ax.tick_params(labelsize = 25)
    ax.axes.xaxis.set_ticklabels([])
    # m=2
    ax = fig.add_subplot(gs[2, :])
    ax.plot(range(len(yA[2].flatten())),(yA[2].flatten()-yB[2].flatten())*100/yB[2].flatten(),
            color=color0,linestyle=linestyle1,linewidth=lw,alpha=1)
    ax.set_ylabel(dylabel2,fontsize=35)
    ax.set_xlim(xmin,xmax)
    ax.tick_params(labelsize = 25)
    ax.axes.xaxis.set_ticklabels([])
    # m=3
    ax = fig.add_subplot(gs[3, :])
    ax.plot(range(len(yA[3].flatten())),(yA[3].flatten()-yB[3].flatten())*100/yB[3].flatten(),
            color=color0,linestyle=linestyle1,linewidth=lw,alpha=1)
    ax.set_ylabel(dylabel3,fontsize=35)
    ax.set_xlim(xmin,xmax)
    ax.tick_params(labelsize = 25)
    ax.axes.xaxis.set_ticklabels([])
    # m=4
    ax = fig.add_subplot(gs[4, :])
    ax.plot(range(len(yA[4].flatten())),(yA[4].flatten()-yB[4].flatten())*100/yB[4].flatten(),
            color=color0,linestyle=linestyle1,linewidth=lw,alpha=1)
    ax.set_xlabel(xlabel,fontsize=30)
    ax.set_ylabel(dylabel4,fontsize=35)
    ax.set_xlim(xmin,xmax)
    ax.tick_params(labelsize = 25)

    fig.align_labels()
    plt.savefig(ofilename_err,dpi=300)

    # Elige referencia
    if args.ref == "a":
        y_ref, lab_ref = yA, args.label_a
        yX, labX = yB, args.label_b
    else:
        y_ref, lab_ref = yB, args.label_b
        yX, labX = yA, args.label_a

    y_ref_flatten = y_ref[0].flatten()
    yXflatten = yX[0].flatten()

    # Métricas
    mae_X, rmse_X, l1_X = error_stats(yXflatten, y_ref_flatten)

    # Residuales (%)
    resX = residual_percent(yXflatten, y_ref_flatten)

    txtX = f"{labX}: MAE={mae_X:.3f}%  RMSE={rmse_X:.3f}%  L1={l1_X:.3f}%"

    plt.close(fig)
    #E end making plots

    # TXT with metricss
    with open(os.path.join(outdir, f"{base}_metrics.txt"), "w") as f:
        f.write(f"Reference: {lab_ref}\n")
        f.write(f"{txtX}\n")

    print("[OK] Figure:", ofilename)
    print("[OK] Errors figure:", ofilename_err)
    print("[OK] Métricas:")
    print("   ", txtX)

if __name__ == "__main__":
    main()
