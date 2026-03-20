#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Acknowledgements: Axel Romero Tisnado
#
import os, argparse, gc, numpy as np
from pathlib import Path
import healpy as hp
#B for cBalls
import sys
# Determine the absolute path of the target (cyballs) directory
#   these two lines won´t be necessary if cballys is in searching path
#target_directory = os.path.abspath('/opt/homebrew/anaconda3/lib/python3.13/site-packages/')
# Append the directory to sys.path
#sys.path.append(target_directory)
from cyballs import cballs
#E

# -------- cBalls --------
def cballs_corr(tmin_deg, tmax_deg, nbins, nthreads=16, outdir='Output',
                fitsfile1=None, fitsfile2=None, fitsfile3=None,
                mask=None, fileformat=None, type='sincos'):
    rmin=np.radians(tmin_deg)
    rmax=np.radians(tmax_deg)
    print('rmin, rmax in degrees', np.degrees(rmin), np.degrees(rmax))
    print('rmin, rmax in arcmin', np.degrees(rmin)*60, np.degrees(rmax)*60)

    Balls = cballs()
    Balls.set({'searchMethod':'octree-g1g2g3-omp'})
    if mask == None:
        if fileformat==None:
            Balls.set({'infile':fitsfile1+','+fitsfile2+','+fitsfile3})
            Balls.set({'infileformat':'fits-healpix,fits-healpix,fits-healpix'})
            Balls.set({'iCatalogs':'1,2,3'})
        else:
            Balls.set({'infile':fitsfile1+fitsfile2+fitsfile3})
            Balls.set({'infileformat':fileformat+fileformat+fileformat})
            Balls.set({'iCatalogs':'1,2,3'})
    else:
        Balls.set({'infile':fitsfile1+fitsfile2+fitsfile3+','+mask})
        Balls.set({'infileformat':'fits-healpix,fits-healpix,fits-healpix,fits-healpix'})
        Balls.set({'iCatalogs':'1,2,3,4'})
    Balls.set({'rangeN':rmax})
    Balls.set({'rminHist':rmin})
    Balls.set({'sizeHistN':nbins})
    Balls.set({'numberThreads':nthreads})
    Balls.set({'rootDir':outdir})
    if mask == None:
        Balls.set({'options':'compute-HistN,out-m-HistZeta'})
    else:
        Balls.set({'options':'compute-HistN,out-m-HistZeta,read-mask'})
    Balls.set({'verbose':2})
    Balls.set({'verbose_log':2})
    Balls.set({'theta':1.0})

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
    xi12 = Balls.getHistXi2pcf12()
    xi13 = Balls.getHistXi2pcf13()
    nn = Balls.getHistNN()

    print('Searching cputime=',cputime,' sec.')
    Balls.clean_all()
    gc.collect()

    return rr, xi12, xi13, nn, monopolesData

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
    ap.add_argument("--fits1", type=Path, required=True, help="Path to the Healpix (fits) map.")
    ap.add_argument("--fits2", type=Path, required=True, help="Path to the Healpix (fits) map.")
    ap.add_argument("--fits3", type=Path, required=True, help="Path to the Healpix (fits) map.")
    ap.add_argument("--file-format", default=None, help="Format of the input map.")
    ap.add_argument("--mask", default=None,
                    help="Path to healpix mask map FITS.")
    ap.add_argument("--field", type=int, default=0, help="Fits' field map (default 0).")
    ap.add_argument("--outdir", default="Output", help="Output directory (default: 'Output').")
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

    #B Setting kappa catalog 1:
    if args.file_format==None:
        if args.nside_down == 0:
            print('nside_down: nothing to do...')
            k2save1 = str(args.fits1)
            print('fits file1: ', k2save1)
            base1 = os.path.splitext(os.path.basename(args.fits1))[0]
        else:
            kappa1, nside_in1 = load_kappa(args.fits1, field=args.field)
            # map to target NSIDE
            k21, nside_eff1 = prepare_samples(kappa1, nside_in1,
                    nside_down=(args.nside_down if args.nside_down>0 else None))
            outdir = args.outdir
            os.makedirs(outdir, exist_ok=True)
            base1 = os.path.splitext(os.path.basename(args.fits1))[0]
            k2save1=os.path.join(outdir, f"{base1}_degraded.fits")
            hp.fitsfunc.write_map(k2save1, k21, overwrite=True)
    else:
        k2save1 = args.fits1
        print('map file1: ', k2save1)
        base1 = os.path.splitext(os.path.basename(args.fits1))[0]
    #E

    #B Setting kappa catalog 2:
    if args.file_format==None:
        if args.nside_down == 0:
            print('nside_down: nothing to do...')
            k2save2 = str(args.fits2)
            print('fits file2: ', k2save2)
            base2 = os.path.splitext(os.path.basename(args.fits2))[0]
        else:
            kappa2, nside_in2 = load_kappa(args.fits2, field=args.field)
            # map to target NSIDE
            k22, nside_eff2 = prepare_samples(kappa2, nside_in2,
                    nside_down=(args.nside_down if args.nside_down>0 else None))
            outdir = args.outdir
            os.makedirs(outdir, exist_ok=True)
            base2 = os.path.splitext(os.path.basename(args.fits2))[0]
            k2save2=os.path.join(outdir, f"{base2}_degraded.fits")
            hp.fitsfunc.write_map(k2save2, k22, overwrite=True)
    else:
        k2save2 = args.fits2
        print('map file2: ', k2save2)
        base2 = os.path.splitext(os.path.basename(args.fits2))[0]
    #E

    #B Setting kappa catalog 3:
    if args.file_format==None:
        if args.nside_down == 0:
            print('nside_down: nothing to do...')
            k2save3 = str(args.fits3)
            print('fits file3: ', k2save3)
            base3 = os.path.splitext(os.path.basename(args.fits3))[0]
        else:
            kappa3, nside_in3 = load_kappa(args.fits3, field=args.field)
            # map to target NSIDE
            k23, nside_eff3 = prepare_samples(kappa3, nside_in3,
                    nside_down=(args.nside_down if args.nside_down>0 else None))
            outdir = args.outdir
            os.makedirs(outdir, exist_ok=True)
            base3 = os.path.splitext(os.path.basename(args.fits3))[0]
            k2save3=os.path.join(outdir, f"{base3}_degraded.fits")
            hp.fitsfunc.write_map(k2save3, k23, overwrite=True)
    else:
        k2save3 = args.fits3
        print('map file3: ', k2save3)
        base3 = os.path.splitext(os.path.basename(args.fits3))[0]
    #E

    # Computing convergence 2pcf/3pcf
    rr, xi12, xi13, nn, zeta = cballs_corr(args.theta_min, args.theta_max,
            args.nbins, nthreads=args.threads, outdir=args.outdir,
            fitsfile1=k2save1, fitsfile2=k2save2, fitsfile3=k2save3,
            mask=args.mask, fileformat=args.file_format, type=args.type)
    print('rr:',rr)
    print('xi12:',xi12)
    print('xi13:',xi13)
    print('nn',nn)
    print()
    print()
    print('zeta monopole:')
    print(zeta[0])

if __name__ == "__main__":
    main()
