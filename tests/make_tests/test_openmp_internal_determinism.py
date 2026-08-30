#!/usr/bin/env python3

import os
import tempfile
from pathlib import Path

import numpy as np

from cyballs import cballs


TEST_DIR = Path(__file__).resolve().parent.parent
CATALOG = TEST_DIR / "catalogs/Abraham/kappa_nres12_zs9NS256r000.bin"
THREADS = int(os.environ.get("CBALLS_DETERMINISM_THREADS", "4"))
REPEATS = int(os.environ.get("CBALLS_DETERMINISM_REPEATS", "3"))


def run_once(root_dir, threads):
    balls = cballs()
    balls.set(
        {
            "searchMethod": "octree-ggg-omp",
            "infile": str(CATALOG),
            "infileformat": "binary",
            "rangeN": 0.0633205,
            "rminHist": 0.00313811,
            "sizeHistN": 20,
            "mChebyshev": 7,
            "numberThreads": threads,
            "verbose": 0,
            "verbose_log": 0,
            "rootDir": root_dir,
            "options": "compute-HistN,and-CF,out-m-HistZeta,KKKCorrelation",
        }
    )
    try:
        balls.Run(level=["MainLoop"])
        result = {
            "histNN": balls.getHistNN().copy(),
            "histCF": balls.getHistCF().copy(),
            "histXi2pcf": balls.getHistXi2pcf().copy(),
        }
        for m in range(1, 9):
            for component in range(1, 5):
                result[f"histZetaM[{m},{component}]"] = (
                    balls.getHistZetaMsincos(m, component).copy()
                )
        return result
    finally:
        balls.struct_cleanup()


def assert_identical(reference, candidate, run):
    if reference.keys() != candidate.keys():
        raise AssertionError(f"run {run} returned a different histogram set")
    for name in reference:
        if not np.array_equal(reference[name], candidate[name]):
            changed = np.count_nonzero(reference[name] != candidate[name])
            raise AssertionError(
                f"run {run} changed {changed} values in {name}"
            )


def test_openmp_runs_are_bitwise_identical_across_thread_counts():
    if THREADS <= 1:
        raise AssertionError("CBALLS_DETERMINISM_THREADS must exceed one")
    if REPEATS < 2:
        raise AssertionError("CBALLS_DETERMINISM_REPEATS must be at least two")
    if not CATALOG.is_file():
        raise AssertionError(f"determinism catalog is missing: {CATALOG}")

    with tempfile.TemporaryDirectory(prefix="ctreeballs-omp-determinism-") as root:
        reference = run_once(root, 1)
        for run in range(1, REPEATS + 1):
            assert_identical(reference, run_once(root, THREADS), run)


if __name__ == "__main__":
    test_openmp_runs_are_bitwise_identical_across_thread_counts()
    print(
        f"PASS: 1-thread reference and {REPEATS} {THREADS}-thread runs "
        "are internally bitwise identical"
    )
