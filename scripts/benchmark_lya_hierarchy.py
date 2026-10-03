#!/usr/bin/env python3
"""Before/after Ly-alpha hierarchy benchmark, retaining raw products and provenance.

Fresh Python processes; optional MPI; identical catalogs/bins/zero-slop defaults.
Native MainLoop and complete Python Run timings are reported separately.
Both exclude Python input loading, catalog registration and result extraction.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shlex
import statistics
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
METHODS = [
    f"lya-{prefix}{stat}-{parallel}"
    for prefix in ("", "los-tree-")
    for parallel in ("omp", "mpi")
    for stat in ("2pcf", "3pcf", "2pcf-3pcf")
] + ["lya-anisotropic-multipole-3pcf-omp"]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def worker(case):
    import numpy as np
    import resource

    mpi = case["parameters"]["searchMethod"].endswith("-mpi")
    if mpi:
        from mpi4py import MPI

        comm = MPI.COMM_WORLD
    else:
        comm = None
    module = Path(case["module"]).resolve()
    sys.path.insert(0, str(module))
    import cyballs

    if Path(cyballs.__file__).resolve().parent != module:
        raise RuntimeError("wrong cyballs extension loaded")
    with np.load(case["catalog"], allow_pickle=False) as z:
        arrays = {k: z[k] for k in ("positions", "delta", "weights", "forest_ids")}
    out = Path(case["output"])
    for iteration in range(case["warmups"] + 1):
        m = cyballs.cballs()
        try:
            m.set(case["parameters"] | {"rootDir": str(out / f"run-{iteration}")})
            m.set_forest_catalog(**arrays)
            if comm:
                comm.Barrier()
            start = time.perf_counter()
            cpu = time.process_time()
            m.Run(level=["MainLoop"])
            run_wall = time.perf_counter() - start
            run_cpu = time.process_time() - cpu
            native = m.getTimings()
            wall = native["wall_seconds"]
            cpu = native["process_cpu_seconds"]
            rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * (
                1 if sys.platform == "darwin" else 1024
            )
            if comm:
                wall = comm.allreduce(wall, op=MPI.MAX)
                cpu = comm.allreduce(cpu, op=MPI.SUM)
                run_wall = comm.allreduce(run_wall, op=MPI.MAX)
                run_cpu = comm.allreduce(run_cpu, op=MPI.SUM)
                max_rss = comm.allreduce(rss, op=MPI.MAX)
                sum_rss = comm.allreduce(rss, op=MPI.SUM)
            else:
                max_rss = sum_rss = rss
            if iteration == case["warmups"] and (comm is None or comm.rank == 0):
                np.savez_compressed(
                    out / "products.npz", **m.getForestResults()["arrays"]
                )
                result = dict(
                    wall_seconds=wall,
                    process_cpu_seconds=cpu,
                    run_wall_seconds=run_wall,
                    run_process_cpu_seconds=run_cpu,
                    maximum_rank_peak_rss_bytes=max_rss,
                    sum_rank_peak_rss_bytes=sum_rss,
                    metadata=m.getRunMetadata(),
                    parameters=case["parameters"],
                    extension=str(cyballs.__file__),
                    extension_sha256=sha(cyballs.__file__),
                )
                (out / "sample.json").write_text(json.dumps(result, indent=2) + "\n")
        finally:
            m.struct_cleanup()
        if comm:
            comm.Barrier()


def compare(reference, candidate):
    import numpy as np

    if set(reference) != set(candidate):
        raise AssertionError("product names differ")
    result = {}

    def measure(name, a, b):
        if a.shape != b.shape:
            raise AssertionError((name, "shape mismatch"))
        finite = np.isfinite(a) & np.isfinite(b)
        diff = b[finite] - a[finite]
        result[name] = dict(
            passed=bool(np.allclose(a, b, rtol=3e-11, atol=1e-10, equal_nan=False)),
            same_finite_support=bool(np.array_equal(np.isfinite(a), np.isfinite(b))),
            relative_l2=float(
                np.linalg.norm(diff) / max(np.linalg.norm(a[finite]), 1e-300)
            ),
            max_absolute=float(np.max(np.abs(diff), initial=0)),
        )

    for name in reference:
        measure(name, reference[name], candidate[name])
    for stat in ("pair", "triple"):
        n, d = stat + "_numerator", stat + "_denominator"
        if n not in reference:
            continue
        result[stat + "_occupancy"] = dict(
            passed=bool(np.array_equal(reference[d] > 0, candidate[d] > 0))
        )
        measure(
            stat + "_ratio",
            np.divide(
                reference[n],
                reference[d],
                out=np.zeros_like(reference[n]),
                where=reference[d] > 0,
            ),
            np.divide(
                candidate[n],
                candidate[d],
                out=np.zeros_like(candidate[n]),
                where=candidate[d] > 0,
            ),
        )
    return result


def main():
    if len(sys.argv) == 3 and sys.argv[1] == "--worker":
        worker(json.loads(Path(sys.argv[2]).read_text()))
        return
    import numpy as np

    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--reference-module", type=Path, default=ROOT)
    p.add_argument("--candidate-module", type=Path, default=ROOT)
    p.add_argument("--reference-parameters", default='{"lya2Kernel":0,"lya3Kernel":0}')
    p.add_argument("--candidate-parameters", help="JSON overrides; default kernels 1/5 for hard bins, 0/0 for anisotropic multipoles")
    p.add_argument("--methods", nargs="+", choices=METHODS, default=METHODS[:3])
    p.add_argument(
        "--catalog", type=Path, help="NPZ: positions, delta, weights, forest_ids"
    )
    p.add_argument(
        "--geometry", choices=("forests", "clustered", "bent"), default="forests"
    )
    p.add_argument("--forests", type=int, default=32)
    p.add_argument("--pixels", type=int, default=64)
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--ranks", type=int, default=2)
    p.add_argument("--mpi-command", default="mpiexec")
    p.add_argument("--radius", type=float, default=160)
    p.add_argument("--pair-bins", type=int, default=30)
    p.add_argument("--radial-bins", type=int, default=4)
    p.add_argument("--polar-bins", type=int, default=4)
    p.add_argument("--mu-bins", type=int, default=4)
    p.add_argument("--warmups", type=int, default=1)
    p.add_argument("--repeats", type=int, default=3)
    p.add_argument("--timeout", type=float, default=600)
    p.add_argument("--outdir", type=Path, required=True)
    a = p.parse_args()
    if (
        min(
            a.forests,
            a.pixels,
            a.threads,
            a.ranks,
            a.pair_bins,
            a.radial_bins,
            a.polar_bins,
            a.mu_bins,
            a.repeats,
        )
        < 1
        or a.warmups < 0
    ):
        p.error("invalid counts")
    if not np.isfinite(a.radius) or a.radius <= 0:
        p.error("radius must be positive and finite")
    out = a.outdir.resolve()
    out.mkdir(parents=True, exist_ok=False)
    if a.catalog:
        with np.load(a.catalog, allow_pickle=False) as z:
            data = {k: z[k] for k in ("positions", "delta", "weights", "forest_ids")}
    else:
        rng = np.random.default_rng(187321)
        directions = np.column_stack(
            (np.ones(a.forests), rng.uniform(-0.13, 0.13, (a.forests, 2)))
        )
        directions /= np.linalg.norm(directions, axis=1)[:, None]
        chi = rng.uniform(80, 260, (a.forests, a.pixels, 1))
        xyz = directions[:, None, :] * chi
        if a.geometry == "clustered":
            xyz = directions[:, None, :] * rng.uniform(
                100, 140, (a.forests, 1, 1)
            ) + 0.0001 * rng.normal(size=xyz.shape)
        if a.geometry == "bent":
            xyz += 2 * rng.normal(size=xyz.shape)
        data = dict(
            positions=xyz.reshape(-1, 3),
            delta=rng.normal(size=a.forests * a.pixels),
            weights=rng.uniform(0.1, 2, a.forests * a.pixels),
            forest_ids=np.repeat(np.arange(a.forests, dtype=np.int64), a.pixels),
        )
    cat = out / "catalog.npz"
    np.savez_compressed(cat, **data)
    report = dict(
        arguments=vars(a),
        pixels=len(data["delta"]),
        forest_count=int(len(np.unique(data["forest_ids"]))),
        catalog_sha256=sha(cat),
        timing_scope="native MainLoop; maximum rank wall and sum rank process CPU",
        memory_scope="per-process high-water RSS including Python and imported libraries; sum of rank maxima is not simultaneous peak",
        samples=[],
        comparisons=[],
        summary=[],
    )

    def save():
        (out / "summary.json").write_text(
            json.dumps(report, indent=2, default=str, allow_nan=False) + "\n"
        )

    env = dict(os.environ, OMP_DYNAMIC="FALSE", OMP_WAIT_POLICY="PASSIVE")
    for method in a.methods:
        result = {"reference": [], "candidate": []}
        params = dict(
            searchMethod=method,
            numberThreads=a.threads,
            verbose=1,
            verbose_log=0,
            usePeriodic=False,
            useLogHist=False,
            rangeN=a.radius,
            rminHist=0.1,
            sizeHistN=4,
            lya2RpMax=a.radius,
            lya2RtMax=a.radius,
            lya2RpBins=a.pair_bins,
            lya2RtBins=a.pair_bins,
            lya3RMax=a.radius,
            lya3RBins=a.radial_bins,
            lya3ThetaBins=a.polar_bins,
            lya3MuBins=a.mu_bins,
            options="no-out-Hist,no-smooth-pivot",
        )
        for rep in range(a.repeats):
            for mode in (
                ("reference", "candidate")
                if rep % 2 == 0
                else ("candidate", "reference")
            ):
                folder = out / f"{method}-{mode}-{rep}"
                folder.mkdir()
                overrides = getattr(a, mode + "_parameters")
                if overrides is None:
                    overrides = '{"lya2Kernel":0,"lya3Kernel":0}' if "anisotropic-multipole" in method else '{"lya2Kernel":1,"lya3Kernel":5}'
                settings = params | json.loads(overrides)
                if "-2pcf-" in method and "-3pcf-" not in method:
                    settings["lya3Kernel"] = 0
                if "-3pcf-" in method and "-2pcf-" not in method:
                    settings["lya2Kernel"] = 0
                case = dict(
                    module=str(getattr(a, mode + "_module").resolve()),
                    catalog=str(cat),
                    output=str(folder),
                    parameters=settings,
                    warmups=a.warmups,
                )
                spec = folder / "case.json"
                spec.write_text(json.dumps(case, indent=2))
                command = [
                    sys.executable,
                    str(Path(__file__).resolve()),
                    "--worker",
                    str(spec),
                ]
                if method.endswith("-mpi"):
                    command = (
                        shlex.split(a.mpi_command) + ["-n", str(a.ranks)] + command
                    )
                with (folder / "process.log").open("w") as log:
                    subprocess.run(
                        command,
                        env=env,
                        stdout=log,
                        stderr=subprocess.STDOUT,
                        check=True,
                        timeout=a.timeout,
                    )
                sample = json.loads((folder / "sample.json").read_text())
                sample.update(method=method, mode=mode, repeat=rep, command=command)
                report["samples"].append(sample)
                result[mode].append((folder, sample))
                save()
                print(method, mode, rep, round(sample["wall_seconds"], 6), flush=True)
        reference = np.load(result["reference"][0][0] / "products.npz")
        for mode in ("reference", "candidate"):
            first = np.load(result[mode][0][0] / "products.npz")
            for folder, sample in result[mode]:
                current = np.load(folder / "products.npz")
                for key in first:
                    np.testing.assert_array_equal(first[key], current[key])
                metrics = compare(reference, current)
                report["comparisons"].append(
                    dict(
                        method=method,
                        mode=mode,
                        repeat=sample["repeat"],
                        metrics=metrics,
                    )
                )
        times = {
            mode: statistics.median(s["wall_seconds"] for _, s in result[mode])
            for mode in result
        }
        passed = all(
            m["passed"]
            for x in report["comparisons"]
            if x["method"] == method
            for m in x["metrics"].values()
        )
        report["summary"].append(
            dict(
                method=method,
                median_wall=times,
                speedup=times["reference"] / times["candidate"],
                passed=passed,
            )
        )
        save()
    if not all(x["passed"] for x in report["summary"]):
        raise SystemExit(
            "FAILED numerical comparison; timings retained but not accepted"
        )


if __name__ == "__main__":
    main()
