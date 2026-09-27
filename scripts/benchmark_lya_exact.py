#!/usr/bin/env python3
"""Compare two exact Ly-alpha executables with retained accuracy, CPU, wall and RSS.

Requires Python 3.10+ and NumPy. No cyballs installation is required. The seeded
catalog recipe is the one in tests/python/lya_corr_all_engines.py. Binaries must
be built with identical compiler settings. Run on an otherwise idle host.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import os
from pathlib import Path
import platform
import re
import statistics
import subprocess
import threading
import time

import numpy as np


METHODS = ("lya-2pcf-omp", "lya-3pcf-omp", "lya-2pcf-3pcf-omp")
PRODUCTS = (("histXi2pcf_lya.txt", 2, 5), ("histZetaM_lya5d.txt", 5, 11))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def catalog(path, shape, seed):
    forests, pixels = map(int, shape.lower().split("x"))
    if forests < 3 or pixels < 1:
        raise ValueError("catalog dimensions must be >=3 forests x >=1 pixel")
    rng = np.random.default_rng(seed)
    ra, dec = rng.uniform(-.008, .008, (2, forests))
    los = np.column_stack((np.cos(dec)*np.cos(ra), np.cos(dec)*np.sin(ra), np.sin(dec)))
    chi = 4000 + np.arange(pixels)*2.1 + rng.uniform(-.4, .4, (forests, 1))
    pos = (chi[..., None]*los[:, None, :]).reshape(-1, 3)
    rows = np.column_stack((pos, rng.normal(0, .2, len(pos)), rng.uniform(.5, 2, len(pos)),
                            np.repeat(np.arange(forests), pixels)))
    np.savetxt(path, rows, fmt=("%.17g",)*5 + ("%.0f",))
    return dict(forests=forests, pixels_per_forest=pixels, nbody=len(pos), seed=seed,
                sha256=digest(path))


def command(binary, input_path, directory, method, threads, write):
    return [str(binary), f"search={method}", f"infile={input_path}",
            "infileformat=lya-ascii", "iCatalogs=1", f"rootDir={directory}",
            f"numberThreads={threads}", "usePeriodic=false", "useLogHist=false",
            "rangeN=80", "rminHist=0.1", "sizeHistN=8", "theta=0",
            "lya2RpMax=40", "lya2RtMax=40", "lya2RpBins=20", "lya2RtBins=20",
            "lya3RMax=25", "lya3RBins=8", "lya3ThetaBins=8", "lya3MuBins=10",
            "verbose=2", "verbose_log=0", "options=" + ("" if write else "no-out-Hist")]


def run(binary, input_path, directory, method, threads, write):
    directory.mkdir()
    native_command = command(binary, input_path, directory, method, threads, write)
    mac = platform.system() == "Darwin"
    start = time.perf_counter()
    with (directory / "process.log").open("w") as stream:
        proc = subprocess.Popen(native_command, stdout=stream, stderr=subprocess.STDOUT,
                                env=dict(os.environ, OMP_DYNAMIC="FALSE", OMP_WAIT_POLICY="PASSIVE"))
        timeout = threading.Timer(7200, proc.kill)
        timeout.start()
        try:
            _, status, usage = os.wait4(proc.pid, 0)
            proc.returncode = os.waitstatus_to_exitcode(status)
        finally:
            timeout.cancel()
    wall = time.perf_counter() - start
    log = (directory / "process.log").read_text()
    if proc.returncode:
        raise RuntimeError(f"exit {proc.returncode}: {directory / 'process.log'}")
    match = re.search(re.escape(method) + r": accepted=(\d+) pairs=(\d+) ordered_triplets=(\d+) CPU=([\d.eE+-]+)", log)
    if not match:
        raise RuntimeError(f"missing search CPU/counters: {directory / 'process.log'}")
    return dict(command=native_command, directory=str(directory), wall_seconds=wall,
                search_cpu_seconds=float(match[4]), process_cpu_seconds=usage.ru_utime + usage.ru_stime,
                accepted=int(match[1]), pairs=int(match[2]), ordered_triplets=int(match[3]),
                peak_rss_bytes=int(usage.ru_maxrss) * (1 if mac else 1024))


def compare(left, right, method):
    results = {}
    for name, keys, col in PRODUCTS:
        if (name.startswith("histXi") and method == "lya-3pcf-omp") or (
            name.startswith("histZeta") and method == "lya-2pcf-omp"):
            continue
        def load(path):
            lines = [line for line in path.read_text().splitlines()
                     if line.strip() and not line.startswith("#")]
            return (np.atleast_2d(np.loadtxt(io.StringIO("\n".join(lines)))) if lines
                    else np.empty((0, col+2)))
        a, b = (load(p / name) for p in (left, right))
        np.testing.assert_array_equal(a[:, :keys], b[:, :keys])
        # Numerator, denominator AND normalized correlation: cancellation can
        # make per-bin relative correlation errors misleading near zero.
        metrics = {}
        for label, index in (("numerator", col), ("denominator", col+1), ("correlation", col-1)):
            x, y = a[:, index], b[:, index]
            np.testing.assert_allclose(y, x, rtol=3e-12, atol=3e-12)
            diff = y-x
            norm = float(np.linalg.norm(x))
            metrics[label] = dict(max_absolute=float(np.max(np.abs(diff), initial=0)),
                                 relative_l2=float(np.linalg.norm(diff)) / max(norm, 1e-300))
        counts = lambda p: re.findall(r"(?:pairs|triplets): (\d+)", (p / name).read_text())
        assert counts(left) == counts(right), f"count mismatch: {name}"
        results[name] = dict(rows=len(a), occupied_rows=int(np.count_nonzero(a[:, col+1])), metrics=metrics,
                             byte_identical=(left/name).read_bytes() == (right/name).read_bytes())
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", required=True, type=Path)
    parser.add_argument("--candidate", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--pair-catalog", default="128x256", help="forests x pixels")
    parser.add_argument("--triple-catalog", default="32x64", help="forests x pixels")
    parser.add_argument("--catalog", type=Path, help="use the same existing six-column lya-ascii file for all methods")
    parser.add_argument("--pair-input", type=Path, help="existing six-column catalog for pair-only runs")
    parser.add_argument("--triple-input", type=Path, help="existing six-column catalog for 3PCF and combined runs")
    parser.add_argument("--seed", type=int, default=9252026)
    parser.add_argument("--threads", type=int, nargs="+", default=[1, 4])
    parser.add_argument("--repeats", type=int, default=3)
    args = parser.parse_args()
    if args.repeats < 1 or min(args.threads) < 1:
        parser.error("repeats and threads must be positive")
    if args.catalog and (args.pair_input or args.triple_input):
        parser.error("use --catalog or the per-statistic input options, not both")
    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=False)
    binaries = {k: getattr(args, k).resolve() for k in ("baseline", "candidate")}
    report = dict(status="RUNNING", platform=platform.platform(), cpu_count=os.cpu_count(),
                  numpy=np.__version__, binaries={k: dict(path=str(v), sha256=digest(v)) for k, v in binaries.items()},
                  threads=args.threads, repeats=args.repeats, catalogs={}, accuracy=[], samples=[], summary=[],
                  methodology="Fresh processes, paired alternating execution order, OMP_WAIT_POLICY=PASSIVE. "
                  "Accuracy/warmup runs write histograms; timed runs suppress them. Search CPU is summed process "
                  "CPU from clock(); wall includes input, tree construction and process startup. Peak RSS is whole process.")
    paths = {}
    for kind, shape in (("pair", args.pair_catalog), ("triple", args.triple_catalog)):
        paths[kind] = root / (kind + ".txt")
        source = args.catalog or getattr(args, kind + "_input")
        if source:
            paths[kind].write_bytes(source.read_bytes())
            nbody = sum(bool(line.strip()) and not line.lstrip().startswith("#")
                        for line in paths[kind].read_text().splitlines())
            report["catalogs"][kind] = dict(source=str(source.resolve()), nbody=nbody, sha256=digest(paths[kind]))
        else:
            report["catalogs"][kind] = catalog(paths[kind], shape, args.seed)
    def save():
        (root / "benchmark.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    save()
    try:
        for method in METHODS:
            input_path = paths["pair" if method == METHODS[0] else "triple"]
            first_accuracy = {}
            for threads in args.threads:
                prefix = f"{method}-t{threads}"
                dirs = {}
                validations = {}
                for label, binary in binaries.items():
                    dirs[label] = root / f"{prefix}-{label}-accuracy"
                    validations[label] = run(binary, input_path, dirs[label], method, threads, True)
                    if label in first_accuracy:
                        for name, _, _ in PRODUCTS:
                            if (dirs[label]/name).exists():
                                assert (dirs[label]/name).read_bytes() == (first_accuracy[label]/name).read_bytes(), (
                                    f"{label} output depends on thread count: {method}")
                    first_accuracy[label] = dirs[label]
                # Accepted visits measure implementation work and may fall
                # after earlier ownership pruning; scientific counts may not.
                for counter in ("pairs", "ordered_triplets"):
                    assert validations["baseline"][counter] == validations["candidate"][counter], counter
                report["accuracy"].append(dict(method=method, threads=threads,
                                                products=compare(dirs["baseline"], dirs["candidate"], method)))
                save()
                for repeat in range(args.repeats):
                    for label in (("baseline", "candidate") if repeat % 2 == 0 else ("candidate", "baseline")):
                        row = run(binaries[label], input_path, root / f"{prefix}-{label}-{repeat}", method, threads, False)
                        for counter in ("accepted", "pairs", "ordered_triplets"):
                            assert row[counter] == validations[label][counter], counter
                        row.update(method=method, threads=threads, label=label, repeat=repeat,
                                   host_load=os.getloadavg() if hasattr(os, "getloadavg") else None)
                        report["samples"].append(row)
                        save()
                summary = dict(method=method, threads=threads, timings={})
                for label in binaries:
                    rows = [r for r in report["samples"] if (r["method"], r["threads"], r["label"]) == (method, threads, label)]
                    summary["timings"][label] = {
                        metric: dict(median=statistics.median(r[metric] for r in rows),
                                     minimum=min(r[metric] for r in rows), maximum=max(r[metric] for r in rows))
                        for metric in ("search_cpu_seconds", "wall_seconds", "peak_rss_bytes")}
                for metric in ("search_cpu_seconds", "wall_seconds"):
                    a, b = (summary["timings"][label][metric]["median"] for label in binaries)
                    summary[metric + "_speedup"] = a / b if b > 0 else None
                report["summary"].append(summary)
                save()
                print(f"PASS {prefix}: CPU speedup={summary['search_cpu_seconds_speedup']:.3f}; "
                      f"wall speedup={summary['wall_seconds_speedup']:.3f}", flush=True)
        report["status"] = "PASS"
    except Exception as error:
        report.update(status="FAIL", error=str(error))
        raise
    finally:
        save()


if __name__ == "__main__":
    main()
