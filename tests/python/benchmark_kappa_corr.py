#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Acknowledgements: Axel Romero Tisnado
#                   Python script by Tisnado
#                   adapted to use cballys module
#
"""
Benchmark of convergence, kappa(theta) using several NSIDEs,
with repetitions to have estable CPU times.

- CSV summary per NSIDE with mean/std/median/min/max.
- CPU time versus Npix log-log plots.

Output:
  - TXT: wtheta_<tag>_<engine>_ns<NSIDE>_N<N>.txt (or with _repXX if --save-each)
  - CSV per repetition: cputime_<tag>.csv  (long format)
  - CSV summary:        cputime_<tag>_meanstd.csv
  - PNG summary:        cputime_<tag>_<engine>_meanstd.png
"""
import os, time, argparse, gc, numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
#B settings for cBalls
import sys
# if cballys not in path
#   Determine the absolute path to its target directory
#target_directory = os.path.abspath('/opt/homebrew/anaconda3/lib/python3.13/site-packages/')
# Append the directory to sys.path
#sys.path.append(target_directory)
from cballys import cballs
#E

# ---------- estilo ----------
def set_fonts():
    try:
        fm.findfont('Times New Roman', fallback_to_default=False)
        serif = ['Times New Roman','Nimbus Roman No9 L','TeX Gyre Termes','Liberation Serif','DejaVu Serif']
    except Exception:
        serif = ['Nimbus Roman No9 L','TeX Gyre Termes','Liberation Serif','DejaVu Serif']
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": serif,
        "mathtext.fontset": "stix",
        "axes.titlesize": 16, "axes.labelsize": 16,
        "xtick.labelsize": 12, "ytick.labelsize": 12,
        "legend.fontsize": 12
    })

def guess_tag(path):
    b = os.path.basename(path)
    return os.path.splitext(b)[0]

def load_kappa(path):
    if path.lower().endswith(".npy"):
        k = np.load(path)
    elif path.lower().endswith((".fits", ".fit")):
        k = hp.read_map(path, field=0, dtype=np.float64, verbose=False)
    else:
        raise SystemExit(f"[error] not valid format: {path}")
    return k

def degrade_map(kappa, ns_target):
    ns_in = hp.get_nside(kappa)
    if ns_target >= ns_in:
        return kappa.copy(), ns_in, (ns_target != ns_in)
    k2 = hp.ud_grade(kappa, nside_out=ns_target, power=0.0,
                     order_in='RING', order_out='RING')
    return k2, ns_target, True

def catalog_from_map(kappa):
    """Use ALL valid pixels of the map"""
    ns = hp.get_nside(kappa)
    valid = np.isfinite(kappa)
    pix = np.where(valid)[0]
    th, ph = hp.pix2ang(ns, pix)
    ra  = np.degrees(ph)
    dec = 90.0 - np.degrees(th)
    w   = kappa[pix].astype(np.float64)
    w  -= np.nanmean(w)
    return ra, dec, w

# -------- cBalls --------
def wtheta_cballs(tmin_deg, tmax_deg, nbins, nthreads=16, outdir=None, fits=None):

    rmin=np.radians(tmin_deg)
    rmax=np.radians(tmax_deg)

    Balls = cballs()
    Balls.set({'infile':fits})
    Balls.set({'infileformat':'fits-healpix'})
    Balls.set({'rangeN':rmax})
    Balls.set({'rminHist':rmin})
    Balls.set({'sizeHistN':nbins})
    Balls.set({'numberThreads':nthreads})
    if outdir == None:
        Balls.set({'rootDir':'Output'})
    else:
        Balls.set({'rootDir':outdir})
    Balls.set({'verbose':0})
    Balls.set({'verbose_log':0})
    Balls.set({'theta':1.05})

    dt=Balls.Run()

    rr = Balls.getrBins()
    xi = Balls.getHistXi2pcf()
    th = np.degrees(rr)  # rr en rad → deg
    print('cleaning all...')
    Balls.clean_all()
    print('done.')
    gc.collect()

    return th, xi, dt
# ------------------------

def parse_nsides(ns_str):
    ns = []
    for chunk in ns_str.split(','):
        chunk = chunk.strip()
        if not chunk: continue
        ns.append(int(chunk))
    return sorted(set(ns))

def ensure_csv(path):
    if not os.path.isfile(path):
        with open(path, "w") as f:
            f.write("tag,engine,ns_target,repeat,N_points,theta_min,theta_max,nbins,threads,seconds\n")

def write_summary_csv(path, rows):
    with open(path, "w") as f:
        f.write("engine,ns_target,repeats,N_points_mean,time_mean_s,time_std_s,time_median_s,time_min_s,time_max_s\n")
        for r in rows:
            f.write(f"{r['engine']},{r['ns']},{r['repeats']},{int(round(r['N_mean']))},"
                    f"{r['t_mean']:.6f},{r['t_std']:.6f},{r['t_median']:.6f},{r['t_min']:.6f},{r['t_max']:.6f}\n")

def _nside_to_npix(ns): return 12*ns*ns
def _npix_to_nside(npix): return (np.sqrt(npix/12.0)).astype(int)

def make_speed_plot_per_engine(path_png, engine, stats_by_ns, tag, theta_min, theta_max, nbins, threads):
    """
    stats_by_ns: dict[ns] -> {'t_mean':, 't_std':}
    Figure: X = Npix (log), upper axis with NSIDE; Y = time [s] (log).
    """
    set_fonts()
    fig, ax = plt.subplots(figsize=(7.0, 5.0), dpi=140)

    ns = np.array(sorted(stats_by_ns.keys()))
    Npix = _nside_to_npix(ns).astype(float)
    mu = np.array([stats_by_ns[n]['t_mean'] for n in ns], dtype=float)
    sd = np.array([stats_by_ns[n]['t_std']  for n in ns], dtype=float)

    ax.plot(Npix, mu, marker='o', lw=1.8, label=engine.capitalize())
    ax.fill_between(Npix, np.maximum(mu-sd, 1e-12), mu+sd, alpha=0.25)

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel("Number of pixels")
    ax.set_ylabel("Time [s]")
    ax.set_title(f"{engine.capitalize()}", pad=6)
    ax.grid(True, which='both', ls=':', alpha=0.5)
    ax.legend(loc='upper left', frameon=False)

    # eje superior con NSIDE
    ax2 = ax.twiny()
    ax2.set_xscale('log')
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks(Npix)
    ax2.set_xticklabels([f"{n}" for n in ns])
    ax2.set_xlabel("NSIDE")

    note = rf"$\theta\in[{theta_min},\,{theta_max}]$ deg, nbins={nbins}, threads={threads}"
    ax.text(0.98, 0.02, note, transform=ax.transAxes, ha='right', va='bottom', fontsize=11)  # ← antes estaba (0.02,0.02) y ha='left'

    fig.tight_layout()
    fig.savefig(path_png)
    plt.close(fig)

def main():
    ap = argparse.ArgumentParser(description="Benchmark cBalls on kappa HEALPix (several NSIDEs) with repetitions.")
    ap.add_argument("--kappa", required=True,
            help="κ map .fits o .npy (ideally NSIDE=8192).")
    ap.add_argument("--outdir", required=True,
            help="Directory of TXT/CSV/PNG (i.e. ./simulations/cputime)")
    ap.add_argument("--nsides",  default="128,256,512,1024")
    ap.add_argument("--theta-min", type=float, default=0.122505)
    ap.add_argument("--theta-max", type=float, default=3.627997)
    ap.add_argument("--nbins",     type=int,   default=20)
    ap.add_argument("--threads",   type=int,   default=max(1, (os.cpu_count() or 2)//2))
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--tag",  type=str, default=None)

    # Repetitions and modes
    ap.add_argument("--repeats", type=int, default=10,
            help="N repetitions per NSIDE×engine.")
    ap.add_argument("--csv-split", action="store_true",
            help="Save one CSV per repetition (instead of only one).")
    ap.add_argument("--resample-each", action="store_true",
            help="Regenerate randoms in each repetition (pixel catalog keep fix).")
    ap.add_argument("--save-each", action="store_true",
            help="Save TXT of wtheta in each repetition (default only first one).")

    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    tag = args.tag or guess_tag(args.kappa)

    engine = "cBalls"

    # CSV(s)
    base_csv = os.path.join(args.outdir, f"cputime_{tag}.csv")
    if not args.csv_split:
        ensure_csv(base_csv)

    # load base κ
    kappa = load_kappa(args.kappa)

    ns_list = parse_nsides(args.nsides)
    bins = np.linspace(args.theta_min, args.theta_max, args.nbins+1)

    # --- storage for summary/plot
    recs = []   # each element: dict(engine, ns, N, seconds)

    for ns_tgt in ns_list:
        # Target NSIDE map
        k2, ns_eff, _ = degrade_map(kappa, ns_tgt)
        k2save=args.kappa+"_degrade"
        hp.fitsfunc.write_map(k2save, k2, overwrite=True)

        # base catalog (ALL valid pixels)
        ra, dec, w = catalog_from_map(k2)
        N = len(ra)
        if N < 2:
            print(f"[warn] NSIDE={ns_eff}: muy pocos puntos (N={N}), lo salto.")
            # liberar temporales
            del k2
            gc.collect()
            continue

        # randoms base para Corrfunc (mismo tamaño)
        rng_base = np.random.default_rng(args.seed)
        u = rng_base.uniform(-1.0, 1.0, N)
        ra_r_base  = rng_base.uniform(0.0, 360.0, N)
        dec_r_base = np.degrees(np.arcsin(u))

        for rep in range(1, args.repeats+1):
            # re-generar randoms en cada repetición si se pide
            if args.resample_each:
                rng_rep = np.random.default_rng(args.seed + rep)
                u = rng_rep.uniform(-1.0, 1.0, N)
                ra_r = rng_rep.uniform(0.0, 360.0, N)
                dec_r = np.degrees(np.arcsin(u))
            else:
                ra_r, dec_r = ra_r_base, dec_r_base
            th, wth, dt = wtheta_cballs(args.theta_min, args.theta_max,
                            args.nbins, nthreads=args.threads,
                            outdir=args.outdir, fits=k2save)

            # TXT (guarda solo en la primera repetición por defecto)
            rep_suf = f"_rep{rep:02d}" if args.save_each else ""
            if (rep == 1) or args.save_each:
                txt_name = f"wtheta_{tag}_{engine}_ns{ns_eff}_N{N}{rep_suf}.txt"
                np.savetxt(os.path.join(args.outdir, txt_name),
                            np.column_stack([th, wth]),
                            header="theta_deg  w_kappa(theta)")

            # CSV large (single or split)
            csv_path = (os.path.join(args.outdir, f"cputime_{tag}_rep{rep:02d}.csv")
                        if args.csv_split else base_csv)
            ensure_csv(csv_path)
            with open(csv_path, "a") as f:
                f.write(f"{tag},{engine},{ns_eff},{rep},{N},{args.theta_min},"
                        f"{args.theta_max},{args.nbins},{args.threads},{dt:.6f}\n")

            recs.append({"engine":engine, "ns":ns_eff, "N":N, "seconds":dt})
            print(f"[OK] {engine}  NSIDE={ns_eff}  rep={rep}  N={N:,}  time={dt:.3f}s")

            # ---- liberar memoria de la repetición ----
            del th, wth
            gc.collect()

        # liberar memoria de este NSIDE
        del k2, ra, dec, w, ra_r_base, dec_r_base
        gc.collect()

    # ---------- SUMMARY ----------
    if len(recs) == 0:
        print("[warn] there is no registers to summary."); return

    # agrupar por engine, ns
    groups = {}
    for r in recs:
        eng = r["engine"]; ns = r["ns"]
        groups.setdefault(eng, {}).setdefault(ns, {"times":[], "Ns":[]})
        groups[eng][ns]["times"].append(r["seconds"])
        groups[eng][ns]["Ns"].append(r["N"])

    # CSV resumen + figuras por engine
    sum_rows = []
    for eng, dd in groups.items():
        stats_by_ns = {}
        for ns, d in dd.items():
            t = np.array(d["times"], dtype=float)
            Ns = np.array(d["Ns"], dtype=float)
            t_mean = float(np.mean(t))
            t_std  = float(np.std(t, ddof=1)) if len(t)>1 else 0.0
            stats_by_ns[ns] = {"t_mean": t_mean, "t_std": t_std}
            sum_rows.append({
                "engine": eng, "ns": ns, "repeats": len(t),
                "N_mean": float(np.mean(Ns)),
                "t_mean": t_mean, "t_std": t_std,
                "t_median": float(np.median(t)),
                "t_min": float(np.min(t)), "t_max": float(np.max(t)),
            })
        png = os.path.join(args.outdir, f"cputime_{tag}_{eng}_meanstd.png")
        make_speed_plot_per_engine(png, eng, stats_by_ns,
                                   tag, args.theta_min, args.theta_max,
                                   args.nbins, args.threads)
        print(f"[OK] plot saved: {png}")

    sum_csv = os.path.join(args.outdir, f"cputime_{tag}_meanstd.csv")
    write_summary_csv(sum_csv, sorted(sum_rows, key=lambda r:(r["engine"], r["ns"])))
    print(f"[OK] summary saved: {sum_csv}")
    print()

if __name__ == "__main__":
    main()
