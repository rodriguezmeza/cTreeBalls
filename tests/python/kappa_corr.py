#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Acknowledgements: Axel Romero Tisnado
#                   Python script by Tisnado
#                   adapted to use cballys module
#
#
import os, argparse, numpy as np
from pathlib import Path
import healpy as hp
#B for cBalls
import sys
# Determine the absolute path of the target (cballys) directory
#   these two lines won´t be necessary if cballys is in searching path
#target_directory = os.path.abspath('/opt/homebrew/anaconda3/lib/python3.13/site-packages/')
# Append the directory to sys.path
#sys.path.append(target_directory)
from cballys import cballs
#E

# -------- cBalls --------
def cballs_corr(tmin_deg, tmax_deg, nbins, nthreads=16, outdir='Output',
                fitsfile=None, mask=None, fileformat=None, type='sincos'):
    rmin=np.radians(tmin_deg)
    rmax=np.radians(tmax_deg)
    print('rmin, rmax in degrees', np.degrees(rmin), np.degrees(rmax))
    print('rmin, rmax in arcmin', np.degrees(rmin)*60, np.degrees(rmax)*60)

    Balls = cballs()
    if mask == None:
        if fileformat==None:
            Balls.set({'infile':fitsfile})
            Balls.set({'infileformat':'fits-healpix'})
            Balls.set({'iCatalogs':'1'})
        else:
            Balls.set({'infile':fitsfile})
            Balls.set({'infileformat':fileformat})
            Balls.set({'iCatalogs':'1'})
    else:
        Balls.set({'infile':fitsfile+','+mask})
        Balls.set({'infileformat':'fits-healpix,fits-healpix'})
        Balls.set({'iCatalogs':'1,2'})
    Balls.set({'rangeN':rmax})
    Balls.set({'rminHist':rmin})
    Balls.set({'sizeHistN':nbins})
    Balls.set({'numberThreads':nthreads})
    Balls.set({'rootDir':outdir})
    if mask == None:
        Balls.set({'options':'compute-HistN'})
    else:
        Balls.set({'options':'compute-HistN,read-mask'})
    Balls.set({'verbose':0})
    Balls.set({'verbose_log':0})

    print('Running...')
    cputime=Balls.Run()

    nmonopoles = Balls.getnMonopoles()
    
    # getHistZetaM_sincos(m, type): m multipole,
    #   type: 1 - cos; 2 - sin; 3 - sincos; 4 - cossin
    # edge_effects not programmed yet
    monopolesData = []
    if type=="sincos":
        for j in range(1,nmonopoles+1):
            zcos = Balls.getHistZetaMsincos(j, 1)
            zsin = Balls.getHistZetaMsincos(j, 2)
            monopolesData.append(zcos+zsin)
    else:
        if type=="edge_effects":
            for j in range(1,nmonopoles+1):
                zcos = Balls.getHistZetaMsincos(j, 1)
                zsin = Balls.getHistZetaMsincos(j, 2)
                monopolesData.append(zcos+zsin)
    rr = Balls.getrBins()
    xi = Balls.getHistXi2pcf()
    nn = Balls.getHistNN()
    print('Searching cputime=',cputime,' sec.')
    return rr, xi, nn, monopolesData

def load_kappa(path, field=0):
    k = hp.read_map(path, field=0, dtype=np.float64, verbose=False)
    nside = hp.get_nside(k)
    return k, nside

def prepare_samples(kappa, nside, nside_down=None):
    """Always use all valid píxels."""
    if nside_down and (nside_down < nside):
        print('nside_down < nside: ', nside_down, nside)
        k_down = hp.ud_grade(kappa, nside_out=nside_down, order_in='RING', order_out='RING', power=0.0)
        mask = np.isfinite(k_down)
        dk = k_down[mask].astype(np.float64)
        dk -= np.nanmean(dk)
        pix = np.where(mask)[0]
        nside_eff = nside_down
    else:
        print('nside_down and nside: ', nside_down, nside)
        mask = np.isfinite(kappa)
        dk = kappa[mask].astype(np.float64)
        dk -= np.nanmean(dk)
        pix = np.where(mask)[0]
        nside_eff = nside

    return dk, nside_eff

def main():
    ap = argparse.ArgumentParser(description="κ-κ angular correlation with cBalls.")
    ap.add_argument("--fits", type=Path, required=True,
                    help="Path to the Healpix (fits) map.")
    ap.add_argument("--file-format", default=None, help="Format of the input map.")
    ap.add_argument("--mask", default=None,
                    help="Path to healpix mask map FITS.")
    ap.add_argument("--field", type=int, default=0,
                    help="Fits' field map (default 0).")
    ap.add_argument("--outdir", default="Output",
                    help="Output directory (default: 'Output').")
    ap.add_argument("--type", default="sincos", choices=["sincos","edge_effects"],
                    help="Type of convergence matrices (default: 'sincos').")
    ap.add_argument("--threads", type=int, default=max(1, (os.cpu_count() or 2)-1))
    ap.add_argument("--nside-down", type=int, default=0,
                    help="NSIDE to up/down-grade (default 0 = no up/down-grade).")
    #B radians:
    # rmin = 0.00213811
    # rmax = 0.0633205
    ap.add_argument("--theta-min", type=float, default=0.12250467913471644,
                    help="r minimum [deg].")
    ap.add_argument("--theta-max", type=float, default=3.6279974066581295,
                    help="r maximum [deg].")
    #E
    ap.add_argument("--nbins", type=int, default=20)
    args = ap.parse_args()

    #B Setting kappa catalog:
    if args.file_format==None:
        if args.nside_down == 0:
            print('nside_down: nothing to do...')
            k2save = args.fits
            print('fits file: ', k2save)
            base = os.path.splitext(os.path.basename(args.fits))[0]
        else:
            kappa, nside_in = load_kappa(args.fits, field=args.field)
            # map to target NSIDE
            k2, nside_eff = prepare_samples(kappa, nside_in,
                    nside_down=(args.nside_down if args.nside_down>0 else None))
            outdir = args.outdir
            os.makedirs(outdir, exist_ok=True)
            base = os.path.splitext(os.path.basename(args.fits))[0]
            k2save=os.path.join(outdir, f"{base}_degraded.fits")
            hp.fitsfunc.write_map(k2save, k2, overwrite=True)
    else:
        k2save = args.fits
        print('map file: ', k2save)
        base = os.path.splitext(os.path.basename(args.fits))[0]
    #E

    # Computing convergence 2pcf/3pcf
    rr, xi, nn, zeta = cballs_corr(args.theta_min, args.theta_max,
            args.nbins, nthreads=args.threads, outdir=args.outdir,
            fitsfile=k2save, mask=args.mask,
            fileformat=args.file_format, type=args.type)
    print('rr:',rr)
    print('xi:',xi)
    print('nn',nn)
    print()
    print()
    print('zeta monopole:')
    print(zeta[0])

if __name__ == "__main__":
    main()
