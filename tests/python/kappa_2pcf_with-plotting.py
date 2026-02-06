#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Acknowledgements: Axel Romero Tisnado
#                   Python script by Tisnado
#                   adapted to use cballys module
#
#
import os, time, argparse, numpy as np
from astropy.io import fits
import healpy as hp
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
#B for cBalls
import sys
# Determine the absolute path of the target (cballys) directory
#   these two lines won´t be necessary if cballys is in searching path
#target_directory = os.path.abspath('/opt/homebrew/anaconda3/lib/python3.13/site-packages/')
# Append the directory to sys.path
#sys.path.append(target_directory)
from cballys import cballs
#E

def set_fonts():
    try:
        fm.findfont('Times New Roman', fallback_to_default=False)
        serif = ['Times New Roman','Nimbus Roman No9 L',
                'TeX Gyre Termes','Liberation Serif','DejaVu Serif']
    except Exception:
        serif = ['Nimbus Roman No9 L','TeX Gyre Termes',
                'Liberation Serif','DejaVu Serif']
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": serif,
        "mathtext.fontset": "stix",
        "axes.titlesize": 16, "axes.labelsize": 16,
        "xtick.labelsize": 12, "ytick.labelsize": 12,
        "legend.fontsize": 12
    })

def theta_centers(res):
    names = set(res.dtype.names or [])
    if {'theta_min','theta_max'} <= names:
        return 0.5*(res['theta_min']+res['theta_max'])
    if {'thetamin','thetamax'} <= names:
        return 0.5*(res['thetamin']+res['thetamax'])
    if {'rmin','rmax'} <= names:
        return 0.5*(res['rmin']+res['rmax'])
    if 'meanth' in names:
        return res['meanth']
    return np.arange(res.size, dtype=float)

def weight_avg(res):
    names = set(res.dtype.names or [])
    for key in ('weightavg','wavg','w'):
        if key in names:
            return res[key]
    return np.zeros(res.size, dtype=float)

# -------- cBalls --------
def wtheta_cballs(tmin_deg, tmax_deg, nbins, nthreads=16,
                    outdir=None, fits=None, mask=None, xscale="arcmin"):
    rmin=np.radians(tmin_deg)
    rmax=np.radians(tmax_deg)

    Balls = cballs()
    if mask == None:
        Balls.set({'infile':fits})
        Balls.set({'infileformat':'fits-healpix'})
        Balls.set({'iCatalogs':'1'})
    else:
        Balls.set({'infile':fits+','+mask})
        Balls.set({'infileformat':'fits-healpix,fits-healpix'})
        Balls.set({'iCatalogs':'1,2'})
    Balls.set({'rangeN':rmax})
    Balls.set({'rminHist':rmin})
    Balls.set({'sizeHistN':nbins})
    Balls.set({'numberThreads':nthreads})
    if outdir == None:
        Balls.set({'rootDir':'Output'})
    else:
        Balls.set({'rootDir':outdir})
    if mask == None:
        Balls.set({'options':'compute-HistN'})
    else:
        Balls.set({'options':'compute-HistN,read-mask'})

    cputime=Balls.Run()

    rr = Balls.getrBins()
    xi = Balls.getHistXi2pcf()
    if xscale == "degree":
        th = np.degrees(rr)             # rr en rad → deg
    else:
        if xscale == "arcmin":
            th = np.degrees(rr)*60.0    # rr en rad → arcmin
        else:
            th = rr                     # rr en rad

    print('Searching cputime=',cputime,' sec.')

    return th, xi

def load_kappa(path, field=0):
    k = hp.read_map(path, field=field, dtype=float)
    nside = hp.get_nside(k)
    return k, nside

def prepare_samples(kappa, nside, nside_down=None):
    """Always use ALL pixels in the map catalog."""
    if nside_down and (nside_down < nside):
        k_down = hp.ud_grade(kappa, nside_out=nside_down,
                order_in='RING', order_out='RING', power=0.0)
        mask = np.isfinite(k_down)
        dk = k_down[mask].astype(np.float64)
        dk -= np.nanmean(dk)
        pix = np.where(mask)[0]
        nside_eff = nside_down
    else:
        mask = np.isfinite(kappa)
        dk = kappa[mask].astype(np.float64)
        dk -= np.nanmean(dk)
        pix = np.where(mask)[0]
        nside_eff = nside

    return dk, nside_eff

def main():
    ap = argparse.ArgumentParser(description="κ-κ angular correlation with cBalls.")
    ap.add_argument("--fits", required=True, help="Path to the HEALPix's map FITS.")
    ap.add_argument("--mask", default=None,
            help="Path to HEALPix's mask map FITS. FITS and Mask files need to have same Nside")
    ap.add_argument("--field", type=int, default=0, help="FITS' field (default 0).")
    ap.add_argument("--outdir", default=None,
                    help="Output directory (default, near fits file).")
    ap.add_argument("--threads", type=int, default=max(1, (os.cpu_count() or 2)-1))
    ap.add_argument("--nside-down", type=int, default=1024,
                    help="NSIDE to up/down-grade (0=no degrade). Does not work with --mask")
    ap.add_argument("--theta-min", type=float, default=0.122505, help="r minimum [deg].")
    ap.add_argument("--theta-max", type=float, default=3.627997,  help="r maximum [deg].")
    ap.add_argument("--nbins", type=int, default=20)
    # Figure's options
    ap.add_argument("--y-mul-theta", action="store_true",
                    help="Multiplica Y por θ (en grados) al graficar.")
    ap.add_argument("--scale", choices=["loglog","linear"], default="loglog",
                    help="Escala de la figura: loglog (default) o linear.")
    ap.add_argument("--xscale", choices=["radian","arcmin","degree"], default="arcmin",
                    help="x units: arcmin (default), radian or degree.")
    args = ap.parse_args()

    set_fonts()

    engine = "cBalls"

    outdir = args.outdir or os.path.dirname(os.path.abspath(args.fits))
    os.makedirs(outdir, exist_ok=True)

    #B Setting kappa catalog:
    if args.nside_down == 0:
        print('nside_down: nothing to do...')
        k2save = args.fits
        print('fits file: ', k2save)
        nside_eff = 0
        base = os.path.splitext(os.path.basename(args.fits))[0]
    else:
        print('nside_down: doing down sizing...')
        kappa, nside_in = load_kappa(args.fits, field=args.field)
        # map to target NSIDE
        k2, nside_eff = prepare_samples(kappa, nside_in,
                    nside_down=(args.nside_down if args.nside_down>0 else None))
        if args.outdir == None:
            outdir = "Output"
            os.makedirs(outdir, exist_ok=True)
        else:
            outdir = args.outdir
            os.makedirs(outdir, exist_ok=True)
        base = os.path.splitext(os.path.basename(args.fits))[0]
        k2save=os.path.join(outdir, f"{base}_degraded.fits")
        hp.fitsfunc.write_map(k2save, k2, overwrite=True)
        #B here if needed downsize mask file
        #E
    #E

    th, w = wtheta_cballs(args.theta_min, args.theta_max,
            args.nbins, nthreads=args.threads, outdir=outdir,
            fits=k2save, mask=args.mask, xscale=args.xscale)
    err = None

    # 2pcf ASCII saving (wtheta “pure”)
    txt = os.path.join(outdir, f"wtheta_{base}_{engine}.txt")
    if err is None:
        np.savetxt(txt, np.column_stack([th, w]), header="theta_deg  wtheta(delta_k)")
    else:
        np.savetxt(txt, np.column_stack([th, w, err]), header="theta_deg  wtheta  sigma")

    # Figure
    yplot = th * w if args.y_mul_theta else w
    ylabel = r"$\theta\, w_{\kappa\kappa}(\theta)$" if args.y_mul_theta else r"$w_{\kappa\kappa}(\theta)$"

    fig, ax = plt.subplots(figsize=(6.8, 5.2), dpi=140)
    if args.scale == "loglog":
        mask = (th > 0) & (yplot > 0)
        if err is not None:
            ylo = th*w - err*th if args.y_mul_theta else (w - err)
            yhi = th*w + err*th if args.y_mul_theta else (w + err)
            m2 = (th > 0) & (ylo > 0) & (yhi > 0)
            ax.fill_between(th[m2], ylo[m2], yhi[m2], alpha=0.25, linewidth=0)
        ax.plot(th[mask], yplot[mask], lw=1.8)
        ax.set_xscale("log"); ax.set_yscale("log")
    else:
        if err is not None:
            ylo = th*w - err*th if args.y_mul_theta else (w - err)
            yhi = th*w + err*th if args.y_mul_theta else (w + err)
            ax.fill_between(th, ylo, yhi, alpha=0.25, linewidth=0)
        ax.plot(th, yplot, lw=1.8)

    if args.xscale == "degree":
        ax.set_xlabel(r"$\theta\ [{\rm deg}]$")
    else:
        ax.set_xlabel(r"$\theta\ [{\rm arcmin}]$")
    ax.set_ylabel(ylabel)
    ax.set_title(f"{base}  |  NSIDE={args.nside_down}  |  {engine}", pad=6)  # título corto
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, f"wtheta_{base}_{engine}.png"))
    plt.close(fig)

if __name__ == "__main__":
    main()
