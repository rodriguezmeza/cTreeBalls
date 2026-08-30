#!/usr/bin/env python3

import tempfile

import numpy as np

from cyballs import CosmoSevereError, cballs


def catalog_positions(x, y, kappa, weights):
    probe = cballs()
    positions_2d = np.column_stack((x, y))
    try:
        probe.set_catalog(positions_2d, kappa=kappa, weights=weights)
    except CosmoSevereError:
        return np.column_stack((x, y, np.zeros_like(x)))
    probe.clear_catalogs()
    return positions_2d


def run_ggg(positions, kappa, weights, mmax, theta=None):
    with tempfile.TemporaryDirectory(prefix="ctreeballs-ggg-edge-") as root:
        parameters = {
            "searchMethod": "octree-ggg-omp",
            "rangeN": 0.75,
            "rminHist": 0.03,
            "sizeHistN": 3,
            "mChebyshev": mmax,
            "lengthBox": 1.2,
            "numberThreads": 1,
            "verbose": 0,
            "verbose_log": 0,
            "rootDir": root,
            "options": (
                "no-normalize-HistZeta,edge-corrections,"
                "no-one-ball,no-out-Hist"
            ),
        }
        if theta is not None:
            parameters["theta"] = theta

        balls = cballs()
        balls.set(parameters)
        balls.set_catalog(positions, kappa=kappa, weights=weights)
        try:
            balls.Run(level=["MainLoop"])
            raw = []
            edge = []
            for order in range(mmax + 1):
                components = [
                    balls.getHistZetaMsincos(order + 1, component)
                    for component in range(1, 5)
                ]
                raw.append(
                    components[0] + components[1]
                    + 1j*(components[2] - components[3])
                )
                edge.append(balls.getHistZetaM_EE_complex(order + 1))
            return np.asarray(raw), np.asarray(edge)
        finally:
            balls.struct_cleanup()


def test_constant_field_edge_solution_and_theta_invariance():
    rng = np.random.default_rng(531)
    x = rng.uniform(-0.45, 0.45, 160)
    y = rng.uniform(-0.45, 0.45, 160)
    weights = rng.uniform(0.7, 1.3, len(x))
    constant = 0.7
    kappa = np.full(len(x), constant)
    positions = catalog_positions(x, y, kappa, weights)

    for mmax in (2, 3, 4, 7):
        _, edge = run_ggg(positions, kappa, weights, mmax)
        np.testing.assert_allclose(
            edge[0], constant**3, rtol=0.0, atol=5.0e-13
        )
        np.testing.assert_allclose(
            edge[1:], 0.0, rtol=0.0, atol=5.0e-13
        )

    varying = 0.4 + 0.7*x - 0.3*y + 0.2*x*y
    default_raw, default_edge = run_ggg(
        positions, varying, weights, mmax=3
    )
    exact_raw, exact_edge = run_ggg(
        positions, varying, weights, mmax=3, theta=1.0
    )
    np.testing.assert_allclose(
        default_raw, exact_raw, rtol=2.0e-13, atol=2.0e-13
    )
    np.testing.assert_allclose(
        default_edge, exact_edge, rtol=2.0e-12, atol=2.0e-13
    )


if __name__ == "__main__":
    test_constant_field_edge_solution_and_theta_invariance()
    print("PASS: octree-GGG edge correction and theta invariance")
