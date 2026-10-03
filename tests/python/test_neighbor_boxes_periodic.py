#!/usr/bin/env python3

import tempfile

import numpy as np

from cyballs import cballs


LBOX = 10.0
RANGE = 2.4
NBINS = 6


def expected_pair_counts(positions):
    positions = np.mod(positions, LBOX)
    delta = np.abs(positions[:, None, :] - positions[None, :, :])
    delta = np.minimum(delta, LBOX - delta)
    distance = np.sqrt(np.sum(delta * delta, axis=2))
    valid = (~np.eye(len(positions), dtype=bool)) & (distance < RANGE)
    radial_bin = (distance[valid] * NBINS / RANGE).astype(int)
    return np.bincount(radial_bin, minlength=NBINS)


def run_case(positions, threads, *, cute_box_format=False):
    options = ["compute-HistN", "no-out-Hist"]
    if cute_box_format:
        options.append("cute-box-fmt")

    with tempfile.TemporaryDirectory(prefix="ctreeballs-neighbor-boxes-") as root:
        balls = cballs()
        balls.set(
            {
                "searchMethod": "neighbor-boxes-omp",
                "iCatalogs": "1",
                "usePeriodic": "true",
                "useLogHist": "false",
                "lengthBox": LBOX,
                "rangeN": RANGE,
                "rminHist": 0.0,
                "sizeHistN": NBINS,
                "mChebyshev": 2,
                "numberThreads": threads,
                "verbose": 0,
                "verbose_log": 0,
                "rootDir": root,
                "options": ",".join(options),
            }
        )
        balls.set_catalog(
            np.ascontiguousarray(positions),
            kappa=np.zeros(len(positions), dtype=np.float64),
        )
        try:
            balls.Run(level=["MainLoop"])
            return balls.getHistNN().copy()
        finally:
            balls.struct_cleanup()


def test_periodic_coordinate_conventions_and_openmp_determinism():
    rng = np.random.default_rng(83519)
    positions = rng.uniform(0.0, LBOX, size=(96, 3))
    positions[:4] = np.array(
        [
            [0.05, 0.10, 0.15],
            [9.95, 0.10, 0.15],
            [0.05, 9.90, 0.15],
            [0.05, 0.10, 9.85],
        ]
    )

    expected = expected_pair_counts(positions)
    reference = run_case(positions, 1)
    centered = run_case(positions - 0.5 * LBOX, 1)
    translated = run_case(positions + np.array([2.0, -3.0, 5.0]) * LBOX, 1)
    threaded = run_case(positions, 4)
    formatted = run_case(positions, 1, cute_box_format=True)

    for label, result in (
        ("reference", reference),
        ("centered", centered),
        ("translated", translated),
        ("four-thread", threaded),
        ("cute-box-fmt", formatted),
    ):
        if not np.array_equal(result, expected):
            raise AssertionError(
                f"{label} counts differ: expected {expected}, received {result}"
            )


if __name__ == "__main__":
    test_periodic_coordinate_conventions_and_openmp_determinism()
    print("PASS: neighbor-boxes periodic wrapping and OpenMP determinism")
