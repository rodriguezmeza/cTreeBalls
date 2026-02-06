#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Acknowledgements: Axel Romero Tisnado
#                   Python script by Tisnado
#                   adapted to use cballys module
#
import os, re, time, glob, argparse, numpy as np
import matplotlib as mpl; mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
#from Corrfunc.theory.xi import xi as corrfunc_xi
#B for cBalls
import sys
# Determine the absolute path of the target (cballys) directory
#   these two lines won´t be necessary if cballys is in searching path
#target_directory = os.path.abspath('/opt/homebrew/anaconda3/lib/python3.13/site-packages/')
# Append the directory to sys.path
#sys.path.append(target_directory)
from cballys import cballs
#E

# ---------- Estilo Times + tamaños ----------
BASE_FONTSIZE = 14
def set_fonts():
    try:
        fm.findfont('Times New Roman', fallback_to_default=False)
        serif = ['Times New Roman','Times','Nimbus Roman No9 L','DejaVu Serif']
    except Exception:
        serif = ['Nimbus Roman No9 L','TeX Gyre Termes','Liberation Serif','DejaVu Serif']
    mpl.rcParams.update({
        "font.family": "serif",
        "font.serif": serif,
        "mathtext.fontset": "stix",
        "axes.titlesize": BASE_FONTSIZE + 2,
        "axes.labelsize":  BASE_FONTSIZE + 0,
        "xtick.labelsize": BASE_FONTSIZE - 2,
        "ytick.labelsize": BASE_FONTSIZE - 2,
        "legend.fontsize": BASE_FONTSIZE - 2,
        "savefig.dpi": 300,
        "savefig.facecolor": "white",
        "axes.facecolor": "white",
        "figure.figsize": (7.2, 5.8),
    })

# ---------- helpers de I/O ----------
def find_ascii(inpath: str, tag: str):
    if os.path.isfile(inpath): return inpath
    if not os.path.isdir(inpath): raise FileNotFoundError(f"No existe: {inpath}")
    cands = sorted(glob.glob(os.path.join(inpath, "*.ascii")))
    if not cands: raise FileNotFoundError(f"No encontré .ascii dentro de: {inpath}")
    tag_hits   = [p for p in cands if tag.replace(".","p") in os.path.basename(p)]
    halos_hits = [p for p in cands if "halos" in os.path.basename(p).lower()]
    if tag_hits:   return tag_hits[0]
    if halos_hits: return halos_hits[0]
    return cands[0]

def parse_header(path):
    tokens=None; z_val=None; box=None
    with open(path,"r",encoding="utf-8",errors="ignore") as f:
        for line in f:
            if not line.startswith("#"): break
            L=line.strip()
            if L.startswith("#id "): tokens=L[1:].split()
            if "Redshift" in L:
                m=re.search(r"Redshift\s*=\s*([0-9.eE+\-]+)",L)
                if m:
                    try: z_val=float(m.group(1))
                    except: pass
            if re.match(r'#\s*a\s*=',L):
                m=re.search(r"a\s*=\s*([0-9.eE+\-]+)",L)
                if m:
                    try:
                        a=float(m.group(1))
                        if a>0: z_val=1.0/a - 1.0
                    except: pass
            if "Box size" in L or "box size" in L:
                m=re.search(r"Box size:\s*([0-9.eE+\-]+)",L)
                if m:
                    try: box=float(m.group(1))
                    except: pass
    return tokens, z_val, box

# ---------- selección de columnas (robusta) ----------
def _canon(t:str)->str:
    return t.lower().replace(" ", "").replace("[","(").replace("]",")").replace("{","(").replace("}",")").replace("\\","")
def _find_pos_col(tokens, axis):  # axis in {'x','y','z'}
    if not tokens: return None
    cand = {
        'x': ('x','xpos','posx','x_mpc','x(kpc/h)','x(mpc/h)','x[h^-1mpc]'),
        'y': ('y','ypos','posy','y_mpc','y(kpc/h)','y(mpc/h)','y[h^-1mpc]'),
        'z': ('z','zpos','posz','z_mpc','z(kpc/h)','z(mpc/h)','z[h^-1mpc]')
    }[axis]
    toks=[_canon(t) for t in tokens]
    for i,t in enumerate(toks):
        if t in cand or any(t.startswith(cc.split('(')[0]+"(") for cc in cand):
            return i+1
    return None

def _find_contains(tokens, keys):
    if not tokens: return None
    toks=[_canon(t) for t in tokens]
    for i,t in enumerate(toks):
        for k in keys:
            if _canon(k) in t: return i+1
    return None

# ---------- reading of position and masses ----------
def read_xyz_and_mass(ascii_fn, only_centrals=True, x_col=9, y_col=10, z_col=11,
                      mass_col=3, pid_col=None, type_col=None):
    tokens, z_val, Lh = parse_header(ascii_fn)
    arr = np.genfromtxt(ascii_fn, comments="#")
    if arr.size==0: raise RuntimeError("Catálogo vacío")
    if arr.ndim==1: arr=arr.reshape(1,-1)
    ncol=arr.shape[1]

    if tokens:
        _x = _find_pos_col(tokens,'x'); _y = _find_pos_col(tokens,'y'); _z = _find_pos_col(tokens,'z')
        _m = _find_contains(tokens, ["mvir","m_vir","m200","m200c","mass"])
        _pid  = _find_contains(tokens, ["pid","upid"])
        _type = _find_contains(tokens, ["type"])
        x_col = _x or x_col; y_col = _y or y_col; z_col = _z or z_col
        mass_col = _m or mass_col
        pid_col = _pid if _pid is not None else pid_col
        type_col = _type if _type is not None else type_col

    def ok(c): return c is not None and 1<=c<=ncol
    X = arr[:, (x_col-1) if ok(x_col) else 8].astype(np.float32)
    Y = arr[:, (y_col-1) if ok(y_col) else 9].astype(np.float32)
    Z = arr[:, (z_col-1) if ok(z_col) else 10].astype(np.float32)
    M = arr[:, (mass_col-1) if ok(mass_col) else 2].astype(np.float64)

    if only_centrals:
        if ok(pid_col):
            mask_cen = (arr[:, pid_col-1].astype(int)==-1)
        elif ok(type_col):
            mask_cen = (arr[:, type_col-1].astype(int)==0)
        else:
            mask_cen = np.ones(len(arr), dtype=bool)
    else:
        mask_cen = np.ones(len(arr), dtype=bool)

    good = np.isfinite(X)&np.isfinite(Y)&np.isfinite(Z)&np.isfinite(M)&(M>0)&mask_cen
    X,Y,Z,M = X[good],Y[good],Z[good],M[good]
    xyz = np.column_stack([X,Y,Z]).astype(np.float32, copy=False)
    return dict(xyz=xyz, mass=M, z=z_val, box=Lh)

# ---------- BAO helpers ----------
def gaussian_smooth(y, sigma_bins=2):
    sigma=max(0.5,float(sigma_bins)); half=int(3*sigma)
    x=np.arange(-half,half+1,1.0); w=np.exp(-0.5*(x/sigma)**2); w/=w.sum()
    return np.convolve(y,w,mode="same")

def detect_bao_peak(r, r2xi, rmin=80.0, rmax=140.0, smooth_bins=3):
    y=gaussian_smooth(r2xi, smooth_bins)
    sel=(r>=rmin)&(r<=rmax)
    if not np.any(sel): return None,None
    i=np.argmax(y[sel]); j=np.where(sel)[0][0]+i
    return r[j], y[j]

#B ---------- cBalls ----------
def run_bao_cballs(inpath, outdir, tag, box, rmin, rmax, nbins, nthreads,
            all_halos=False, masalogmin=None, masalogmax=None):
    set_fonts()
    os.makedirs(outdir, exist_ok=True)

    engine = "cBalls"

    ascii_fn = find_ascii(inpath, tag)
    print(f"[info] catalog: {ascii_fn}")

    d = read_xyz_and_mass(ascii_fn, only_centrals=(not all_halos))
    L = float(d["box"]) if (d["box"] is not None and np.isfinite(d["box"])) else float(box)

    # counting and cutting masses (log10)
    N0 = d["xyz"].shape[0]
    if (masalogmin is not None) or (masalogmax is not None):
        logM = np.log10(d["mass"])
        m = np.ones_like(logM, dtype=bool)
        if masalogmin is not None: m &= (logM >= float(masalogmin))
        if masalogmax is not None: m &= (logM <= float(masalogmax))
        d["xyz"]  = d["xyz"][m]
        d["mass"] = d["mass"][m]
    Nsel = d["xyz"].shape[0]
    print(f"[info] N_halos: total={N0}  after cutting={Nsel}  (centrals={'no' if all_halos else 'yes'})")
    if Nsel < 10: raise RuntimeError("To few halos to estimate xi(r).")

    np.savetxt(os.path.join(outdir, f"xyz_{tag}.txt"),
                np.column_stack([d["xyz"][:,0], d["xyz"][:,1], d["xyz"][:,2]]),
                header=f"x[h^-1 Mpc] y[h^-1 Mpc] z[h^-1 Mpc]", fmt="%.6f")

    infile = os.path.join(outdir, f"xyz_{tag}.txt")
    Balls = cballs()
    Balls.set({'searchMethod':'neighbor-boxes-omp'})
    Balls.set({'infile':infile})
    Balls.set({'infileformat':'multi-columns-ascii'})
    Balls.set({'columns':'1,2,3'})
    Balls.set({'iCatalogs':'1'})
    Balls.set({'rootDir':outdir})
    Balls.set({'rangeN':rmax})
    Balls.set({'rminHist':rmin})
    Balls.set({'sizeHistN':nbins})
    Balls.set({'stepState':100000})
    Balls.set({'verbose':2})
    Balls.set({'verbose_log':2})
    Balls.set({'numberThreads':nthreads})
    Balls.set({'useLogHist':False})
    Balls.set({'usePeriodic':True})
    Balls.set({'lengthBox':box})
    Balls.set({'options':'compute-HistN,and-CF,only-pos,cute-box-fmt'})

    print('Running...')
    cputime=Balls.Run()
    print('Searching cputime=',cputime,' sec.')

    rmid = Balls.getrBins()
    xir = Balls.getHistCF()
    dd = Balls.getHistNN()
    r2xi = (rmid**2) * xir

    # saving
    np.savetxt(os.path.join(outdir, f"xi_{tag}_{engine}.txt"),
               np.column_stack([rmid, xir]),
               header=f"r_mid[h^-1 Mpc]  xi(r)   N_total={N0}   N_sel={Nsel}", fmt="%.6f")
    np.savetxt(os.path.join(outdir, f"xi_dd_{tag}_{engine}.txt"),
               np.column_stack([rmid, xir, dd.astype(np.int64)]),
               header=f"r_mid  xi  DD_pairs   N_total={N0}   N_sel={Nsel}",
               fmt=["%.6f","%.6e","%d"])

    # plots
    ttl = rf"$z={d['z']:.3f}$" if d['z'] is not None else f"{tag}"
    # xi(r)
    fig1,ax1=plt.subplots(figsize=(7.2,5.8), dpi=160)
    ax1.plot(rmid, xir, lw=1.8, color="tab:blue")
    ax1.set_xlabel(r"$r\,[h^{-1}\mathrm{Mpc}]$"); ax1.set_ylabel(r"$\xi(r)$"); ax1.set_title(ttl, pad=6)
    fig1.tight_layout(); fig1.savefig(os.path.join(outdir, f"xi_{tag}.png")); plt.close(fig1)
    # r^2 xi(r)
    fig2,ax2=plt.subplots(figsize=(7.2,5.8), dpi=160)
    ax2.plot(rmid, r2xi, lw=1.8, color="tab:blue")
    rbao,_ = detect_bao_peak(rmid, r2xi, rmin=80.0, rmax=140.0, smooth_bins=3)
    if rbao is not None:
        ax2.axvline(rbao, color="crimson", ls="--", lw=1.2,
                    label=rf"$r_{{\rm BAO}}\approx {rbao:.1f}\ \mathrm{{Mpc}}\ h^{{-1}}$")
        ax2.legend(loc="lower left", frameon=True)
    ax2.set_xlabel(r"$r\,[h^{-1}\mathrm{Mpc}]$"); ax2.set_ylabel(r"$r^2\,\xi(r)$"); ax2.set_title(ttl, pad=6)
    fig2.tight_layout(); fig2.savefig(os.path.join(outdir, f"xi_r2_{tag}_{engine}.png")); plt.close(fig2)
#E cBalls

def main():
    ap=argparse.ArgumentParser(description="BAO from halos with Corrfunc (periodic box).")
    ap.add_argument("--in", dest="inpath", required=True)
    ap.add_argument("--out",dest="outdir", required=True)
    ap.add_argument("--tag",dest="tag", required=True)
    ap.add_argument("--box",dest="box", type=float, required=True)
    ap.add_argument("--threads", type=int, default=max(1,(os.cpu_count() or 2)-1))
    ap.add_argument("--rmin", type=float, default=5.0)
    ap.add_argument("--rmax", type=float, default=200.0)
    ap.add_argument("--nbins", type=int, default=80)
    # selection
    ap.add_argument("--all-halos", action="store_true",
                    help="Usar centrales + satélites (default: solo centrales)")
    # Mass cut (log10 M/Msun h^-1)
    ap.add_argument("--masamin", type=float, default=None, help="log10(Mvir_min)")
    ap.add_argument("--masamax", type=float, default=None, help="log10(Mvir_max)")
    # columns fallbacks (1-based) if NO header
    ap.add_argument("--x-col", type=int, default=9)
    ap.add_argument("--y-col", type=int, default=10)
    ap.add_argument("--z-col", type=int, default=11)
    ap.add_argument("--mass-col", type=int, default=3)
    ap.add_argument("--pid-col",  type=int, default=None)
    ap.add_argument("--type-col", type=int, default=None)
    args=ap.parse_args()

    set_fonts()

    run_bao_cballs(inpath=args.inpath, outdir=args.outdir, tag=args.tag, box=args.box,
                rmin=args.rmin, rmax=args.rmax, nbins=args.nbins, nthreads=args.threads,
                all_halos=args.all_halos, masalogmin=args.masamin,
                masalogmax=args.masamax)

if __name__ == "__main__":
    main()
