#!/usr/bin/env python3

import tempfile

import numpy as np

from cyballs import cballs


def run_case(positions, kappa, weights, threads, only_2pcf):
    with tempfile.TemporaryDirectory(prefix="ctreeballs-ggg-2pcf-") as root:
        options = [
            "compute-HistN",
            "and-CF",
            "KKKCorrelation",
            "no-one-ball",
            "no-out-Hist",
        ]
        if only_2pcf:
            options.append("only-2pcf")

        balls = cballs()
        balls.set(
            {
                "searchMethod": "octree-ggg-omp",
                "rangeN": 0.7,
                "rminHist": 0.025,
                "sizeHistN": 8,
                "mChebyshev": 4,
                "lengthBox": 1.2,
                "numberThreads": threads,
                "verbose": 0,
                "verbose_log": 0,
                "rootDir": root,
                "options": ",".join(options),
            }
        )
        balls.set_catalog(positions, kappa=kappa, weights=weights)
        try:
            balls.Run(level=["MainLoop"])
            return {
                "histNN": balls.getHistNN().copy(),
                "histCF": balls.getHistCF().copy(),
                "histXi2pcf": balls.getHistXi2pcf().copy(),
            }
        finally:
            balls.struct_cleanup()


def assert_identical(reference, candidate):
    for name in reference:
        if not np.array_equal(reference[name], candidate[name]):
            changed = np.count_nonzero(reference[name] != candidate[name])
            raise AssertionError(f"{name} changed in {changed} bins")


def test_only_2pcf_matches_full_estimator_and_is_deterministic():
    rng = np.random.default_rng(7341)
    positions = rng.uniform(-0.45, 0.45, size=(192, 3))
    kappa = 0.4 + positions[:, 0] - 0.3 * positions[:, 1]
    weights = rng.uniform(0.7, 1.3, len(positions))

    full = run_case(positions, kappa, weights, threads=1, only_2pcf=False)
    only_one = run_case(positions, kappa, weights, threads=1, only_2pcf=True)
    only_many = run_case(positions, kappa, weights, threads=4, only_2pcf=True)

    assert_identical(full, only_one)
    assert_identical(only_one, only_many)


if __name__ == "__main__":
    test_only_2pcf_matches_full_estimator_and_is_deterministic()
    print("PASS: octree-GGG only-2pcf equivalence and OpenMP determinism")
