#!/usr/bin/env python3
"""GGG window-policy regressions and a bounded masked-sky benchmark."""

import argparse
import json
import re
import tempfile
import time
from pathlib import Path

import numpy as np

from cyballs import cballs


def make_catalog(count):
    rng = np.random.default_rng(87314)
    positions = rng.normal(size=(count, 3))
    positions /= np.linalg.norm(positions, axis=1)[:, None]
    positions = np.ascontiguousarray(positions[np.argsort(positions[:, 2])])
    kappa = 0.4 + 0.3*positions[:, 0] - 0.2*positions[:, 1]
    weights = rng.uniform(0.7, 1.3, count)
    mask = ((positions[:, 2] > -0.25) & (positions[:, 0] > 0)).astype(np.uint8)
    return positions, kappa, weights, mask


def run_case(catalog, threads=1, mmax=7, full_window=False,
             raw=False, exact=False, write_windows=False, profile=False):
    with tempfile.TemporaryDirectory(prefix="ctreeballs-ggg-fast-") as tmp:
        root = Path(tmp)
        options = ["read-mask", "compute-HistN", "and-CF", "KKKCorrelation"]
        if full_window:
            options.append("ggg-full-window")
        if raw:
            options.append("no-normalize-HistZeta")
        if exact:
            options.append("no-one-ball")
        if profile:
            options.append("ggg-profile")
        if not write_windows:
            options.append("no-out-Hist")
        balls = cballs()
        balls.set({
            "searchMethod": "octree-ggg-omp",
            "rangeN": 0.0633205,
            "rminHist": 0.00213811,
            "sizeHistN": 20,
            "mChebyshev": mmax,
            "lengthBox": 2,
            "theta": 1,
            "numberThreads": threads,
            "verbose": 0,
            "verbose_log": int(profile),
            "rootDir": str(root),
            "options": ",".join(options),
        })
        positions, kappa, weights, mask = catalog
        balls.set_catalog(positions, kappa=kappa, weights=weights, mask=mask)
        try:
            started = time.perf_counter()
            balls.Run(level=["MainLoop"])
            elapsed = time.perf_counter() - started
            result = {
                "NN": balls.getHistNN().copy(),
                "CF": balls.getHistCF().copy(),
                "KK": balls.getHistXi2pcf().copy(),
            }
            for order in range(1, mmax + 2):
                for component in range(1, 5):
                    result[f"KKK_{order}_{component}"] = (
                        balls.getHistZetaMsincos(order, component).copy()
                    )
        finally:
            balls.struct_cleanup()
        windows = {
            path.name: np.loadtxt(path)
            for path in root.glob("histZetaM_*_N_*.txt")
        }
        stats = None
        log = root / "tmp" / "cballs.log"
        if profile and log.exists():
            match = re.search(
                r"GGG profile: window_orders=(\d+) work=([\d.]+) "
                r"wait=([\d.]+) merge=([\d.]+) wall=([\d.]+)",
                log.read_text(),
            )
            if match:
                stats = dict(zip(
                    ("window_orders", "work", "wait", "merge", "wall"),
                    map(float, match.groups()),
                ))
        return result, windows, elapsed, stats


def assert_identical(reference, candidate):
    assert reference.keys() == candidate.keys()
    for name in reference:
        np.testing.assert_array_equal(reference[name], candidate[name], err_msg=name)


def test_window_policy_preserves_histograms():
    catalog = make_catalog(4099)
    for mmax, raw, exact in ((2, False, False), (7, True, True)):
        short, short_n, _, _ = run_case(
            catalog, mmax=mmax, raw=raw, exact=exact, write_windows=True
        )
        full, full_n, _, _ = run_case(
            catalog, mmax=mmax, raw=raw, exact=exact,
            write_windows=True, full_window=True,
        )
        assert_identical(short, full)
        assert len(short_n) == 4*(mmax + 1)
        assert len(full_n) == 4*(2*mmax + 1)
        for name in short_n:
            np.testing.assert_array_equal(short_n[name], full_n[name])
        assert np.any(short["NN"] > 0)
        assert np.any(short["KKK_1_1"] != 0)


def test_masked_chunk_determinism():
    # More than four 4096-pivot chunks, including fully masked chunks.
    catalog = make_catalog(20483)
    reference, _, _, _ = run_case(catalog, threads=1, mmax=3)
    for _ in range(2):
        candidate, _, _, stats = run_case(
            catalog, threads=4, mmax=3, profile=True
        )
        assert_identical(reference, candidate)
        assert stats is not None
        assert stats["window_orders"] == 4
        assert all(np.isfinite(value) and value >= 0 for value in stats.values())


def benchmark(args):
    catalog = make_catalog(args.nbody)
    measurements = []
    for _ in range(args.repeat):
        result, _, elapsed, stats = run_case(
            catalog, threads=args.threads, full_window=args.full_window,
            exact=args.no_one_ball, profile=True,
        )
        measurements.append({"run_wall_seconds": elapsed, "profile": stats})
    np.savez(args.benchmark.with_suffix(".npz"), **result)
    args.benchmark.write_text(json.dumps({
        "nbody": args.nbody,
        "unmasked": int(np.count_nonzero(catalog[3])),
        "threads": args.threads,
        "full_window": args.full_window,
        "exact": args.no_one_ball,
        "measurements": measurements,
    }, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--benchmark", type=Path, help="write benchmark JSON and NPZ")
    parser.add_argument("--nbody", type=int, default=131072)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--repeat", type=int, default=3)
    parser.add_argument("--full-window", action="store_true")
    parser.add_argument("--no-one-ball", action="store_true")
    args = parser.parse_args()
    if min(args.nbody, args.threads, args.repeat) <= 0:
        parser.error("nbody, threads, and repeat must be positive")
    if args.benchmark:
        benchmark(args)
    else:
        test_window_policy_preserves_histograms()
        test_masked_chunk_determinism()
        print("PASS: GGG fast/full window equivalence and masked OpenMP determinism")
