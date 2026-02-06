#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os, argparse, numpy as np
from pathlib import Path
import healpy as hp
#B for cBalls
#import sys


def load_kappa(path, field=0):
    k = hp.read_map(path, field=0, dtype=np.float64, verbose=False)
    nside = hp.get_nside(k)
    return k, nside

def prepare_samples(kappa, nside, nside_down=None):
    """Always use all valid píxels."""
#    if nside_down and (nside_down < nside):
    if nside_down:
        print('nside_down, nside: ', nside_down, nside)
        k_down = hp.ud_grade(kappa, nside_out=nside_down,
                order_in='RING', order_out='RING', power=0.0, pess=True)
#                order_in='RING', order_out='RING', power=0.0)
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
    ap = argparse.ArgumentParser(description="Down/up grade HEALPix file.")
    ap.add_argument("--fits", type=Path, required=True,
                    help="Path to the Healpix (fits) map.")
    ap.add_argument("--field", type=int, default=0,
                    help="Fits' field map (default 0).")
    ap.add_argument("--outdir", default='./',
                    help="Output directory (default: './').")
    ap.add_argument("--nside-down", type=int, default=0,
                    help="NSIDE to up/down-grade (default 0 = no up/down-grade).")
    args = ap.parse_args()

    if args.nside_down == 0:
        print('nside_down: nothing to do...')
    else:
        kappa, nside_in = load_kappa(args.fits, field=args.field)
        # map to target NSIDE
        k2, nside_eff = prepare_samples(kappa, nside_in,
                    nside_down=(args.nside_down if args.nside_down>0 else None))
        outdir = args.outdir
        os.makedirs(outdir, exist_ok=True)
        base = os.path.splitext(os.path.basename(args.fits))[0]
        k2save=os.path.join(outdir, f"{base}_{args.nside_down}.fits")
        hp.fitsfunc.write_map(k2save, k2, overwrite=True)

if __name__ == "__main__":
    main()
