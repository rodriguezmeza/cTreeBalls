#!/usr/bin/env python3
"""Run cyballs directly from NumPy arrays, without a catalog file."""

import tempfile

import numpy as np

from cyballs import cballs


def make_catalog(nbody=128, seed=1729):
    rng = np.random.default_rng(seed)
    positions = rng.uniform(-0.45, 0.45, size=(nbody, 3))
    kappa = 1.0 + positions[:, 0] - 0.5 * positions[:, 1]
    weights = rng.uniform(0.8, 1.2, size=nbody)
    mask = np.ones(nbody, dtype=bool)
    return positions, kappa, weights, mask


def main():
    positions, kappa, weights, mask = make_catalog()
    positions_before = positions.copy()

    with tempfile.TemporaryDirectory(prefix="cyballs-memory-example-") as root:
        balls = cballs()
        balls.set(
            {
                "searchMethod": "octree-ggg-omp",
                "rangeN": 0.8,
                "rminHist": 0.02,
                "sizeHistN": 12,
                "mChebyshev": 3,
                "lengthBox": 1.0,
                "numberThreads": 2,
                "verbose": 0,
                "verbose_log": 0,
                "rootDir": root,
                "options": "no-out-Hist",
            }
        )

        # C receives an owned copy. Tree centering/reordering cannot alter NumPy.
        balls.set_catalog(
            positions,
            kappa=kappa,
            weights=weights,
            mask=mask,
        )

        try:
            balls.Run(level=["MainLoop"])
            radius = balls.getrBins().copy()
            xi2pcf = balls.getHistXi2pcf().copy()
            zeta_m1 = balls.getHistZetaMsincos(1, 1).copy()
        finally:
            balls.struct_cleanup()

    np.testing.assert_array_equal(positions, positions_before)
    print("r bins:", radius)
    print("2PCF:", xi2pcf)
    print("3PCF m=1, component=1 shape:", zeta_m1.shape)


if __name__ == "__main__":
    main()
