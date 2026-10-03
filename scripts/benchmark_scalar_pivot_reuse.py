#!/usr/bin/env python3
"""Benchmark scalar hierarchical reuse; writes native products and run provenance.

Use --reuse for the opt-in traversal, --exact for the reference, or neither
for the original traversal. Run with --mpi under mpiexec for MPI engines.
Timings exclude catalog generation/loading and product extraction.
"""
import os, sys, time, json, tempfile, argparse
from pathlib import Path
import numpy as np

p = argparse.ArgumentParser()
p.add_argument("--module", default=str(Path(__file__).resolve().parents[1]))
p.add_argument("--engine", choices=("octree", "kdtree", "balltree"), default="balltree")
p.add_argument("--n", type=int, default=32768)
p.add_argument("--reuse", action="store_true")
p.add_argument("--exact", action="store_true")
p.add_argument("--leaf", type=int, default=16)
p.add_argument("--repeat", type=int, default=1)
p.add_argument("--threads", type=int, default=16)
p.add_argument("--tree-theta", type=float, default=1.0)
p.add_argument("--max-n", type=int, default=3)
p.add_argument("--bins", type=int, default=20)
p.add_argument("--mpi", action="store_true")
p.add_argument("--geometry", choices=("fullsky", "clustered"), default="fullsky")
p.add_argument("--catalog-npz", type=Path)
p.add_argument("--sampling-seed", type=int, default=8675309)
p.add_argument("--min-sep", type=float, default=0.02)
p.add_argument("--max-sep", type=float, default=0.5)
p.add_argument(
    "--phase-tol",
    type=float,
    default=float(os.environ.get("CBALLS_SCALAR_PIVOT_TOL", ".1")),
)
p.add_argument(
    "--bin-theta",
    type=float,
    default=float(os.environ.get("CBALLS_SCALAR_BIN_THETA", "0")),
)
p.add_argument("--statistic", choices=("2pcf", "3pcf", "both"), default="3pcf")
p.add_argument("--no-edge", action="store_true")
p.add_argument("--profile", action="store_true")
p.add_argument("--output", required=True)
a = p.parse_args()
if min(a.n, a.repeat, a.threads, a.leaf, a.bins) <= 0 or a.max_n < 2:
    p.error("n, repeat, threads, leaf and bins must be positive; max-n must be >=2")
if a.reuse and a.exact:
    p.error("choose --reuse or --exact")
os.environ["CBALLS_SCALAR_PIVOT_TOL"] = str(a.phase_tol)
os.environ["CBALLS_SCALAR_BIN_THETA"] = str(a.bin_theta)
if a.mpi:
    from mpi4py import MPI

    comm = MPI.COMM_WORLD
else:
    comm = None
sys.path.insert(0, a.module)
import cyballs

if a.catalog_npz:
    z = np.load(a.catalog_npz)
    if a.n > len(z["positions"]):
        raise ValueError("--n exceeds catalog population")
    selected = np.random.default_rng(a.sampling_seed).permutation(len(z["positions"]))[
        : a.n
    ]
    pos = z["positions"][selected].copy()
    field = z["kappa"][selected].copy()
    weights = z["weights"][selected].copy() if "weights" in z else np.ones(len(pos))
    minimum = a.min_sep
    maximum = a.max_sep
else:
    rng = np.random.default_rng(214729)
    pos = rng.normal(size=(a.n, 3))
    pos /= np.linalg.norm(pos, axis=1)[:, None]
    if a.geometry == "clustered":
        centers = rng.normal(size=(128, 3))
        centers /= np.linalg.norm(centers, axis=1)[:, None]
        pos = centers[np.arange(a.n) % 128] + 0.001 * rng.normal(size=(a.n, 3))
        pos /= np.linalg.norm(pos, axis=1)[:, None]
    field = 0.05 * rng.normal(size=a.n)
    weights = rng.uniform(0.3, 1.7, a.n)
    minimum = a.min_sep
    maximum = a.max_sep
Path(a.output).resolve().parent.mkdir(parents=True, exist_ok=True)
records = []
for rep in range(a.repeat):
    with tempfile.TemporaryDirectory(prefix="scalar-hierarchy-benchmark-") as out:
        out = comm.bcast(out if comm.rank == 0 else None, root=0) if comm else out
        m = cyballs.cballs()
        opts = "KKKCorrelation,no-smooth-pivot,no-normalize-HistZeta,weights-norm,no-balltree-tree-cache,no-native-tree-cache,no-out-Hist"
        if a.statistic != "both":
            opts += ",only-" + a.statistic
        if a.statistic != "2pcf" and not a.no_edge:
            opts += ",edge-corrections"
        if a.profile:
            opts += ",dual-node-profile"
        if a.reuse:
            opts += ",scalar-pivot-reuse"
        if a.exact:
            opts += ",no-one-ball,no-two-balls"
        m.set(
            dict(
                searchMethod=a.engine + "-2balls-" + ("mpi" if comm else "omp"),
                usePeriodic=False,
                useLogHist=True,
                rminHist=minimum,
                rangeN=maximum,
                sizeHistN=a.bins,
                sizeHistPhi=max(32, 2 * a.max_n + 2),
                mChebyshev=a.max_n,
                theta=a.tree_theta,
                nsmooth=a.leaf,
                numberThreads=a.threads,
                verbose=1 if a.profile else 0,
                verbose_log=0,
                rootDir=out,
                options=opts,
            )
        )
        m.set_catalog(pos, kappa=field, weights=weights)
        if comm:
            comm.Barrier()
        t = time.perf_counter()
        cpu = time.process_time()
        m.Run(level=["MainLoop"])
        wall = time.perf_counter() - t
        cpu = time.process_time() - cpu
        if comm:
            wall = comm.allreduce(wall, op=MPI.MAX)
            cpu = comm.allreduce(cpu, op=MPI.SUM)
        if comm is None or comm.rank == 0:
            result = {}
            if a.statistic != "2pcf":
                parts = np.array(
                    [
                        [m.getHistZetaMsincos(order, c).copy() for c in range(1, 5)]
                        for order in range(1, a.max_n + 2)
                    ]
                )
                result.update(
                    components=parts,
                    signal=parts[:, 0] + parts[:, 1] + 1j * (parts[:, 2] - parts[:, 3]),
                )
                if not a.no_edge:
                    result.update(
                        corrected=np.array(
                            [
                                m.getHistZetaM_EE_complex(order).copy()
                                for order in range(1, a.max_n + 2)
                            ]
                        ),
                        w0=m.getScalarWindowDiagnostics()["window_monopole"],
                    )
            if a.statistic != "3pcf":
                result.update(pairs=m.getHistNN().copy(), xi=m.getHistXi2pcf().copy())
            np.savez(a.output + ".npz", **result)
            records.append(dict(wall=wall, cpu=cpu, metadata=m.getRunMetadata()))
            print(a.engine, a.n, a.reuse, wall, flush=True)
        m.struct_cleanup()
        if comm:
            comm.Barrier()
if comm is None or comm.rank == 0:
    Path(a.output + ".json").write_text(
        json.dumps(dict(args=vars(a), samples=records), indent=2, default=str)
    )
