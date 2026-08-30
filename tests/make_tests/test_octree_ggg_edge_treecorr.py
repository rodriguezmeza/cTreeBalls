#!/usr/bin/env python3

import tempfile
from pathlib import Path

import numpy as np

from cyballs import CosmoSevereError, cballs

try:
    import treecorr
except ImportError:
    treecorr = None


def require_2d_build(x, y, kappa, weights):
    probe = cballs()
    try:
        probe.set_catalog(
            np.column_stack((x, y)), kappa=kappa, weights=weights
        )
    except CosmoSevereError:
        try:
            import pytest
        except ImportError as exc:
            raise RuntimeError(
                "TreeCorr comparison requires a DEFDIMENSION=2 cyballs build"
            ) from exc
        pytest.skip("TreeCorr flat-sky comparison requires DEFDIMENSION=2")
    probe.clear_catalogs()


def read_window_mode(root, order):
    index = order + 1

    def component(name):
        return np.loadtxt(root / f"histZetaM_{name}_N_{index}.txt")

    return (
        component("cos") + component("sin")
        + 1j*(component("sincos") - component("cossin"))
    )


def solve_treecorr_edge(raw, window, mmax):
    center = window.shape[2] // 2
    orders = np.arange(-mmax, mmax + 1)
    result = np.zeros((mmax + 1,) + raw.shape[:2], dtype=np.complex128)

    for radial_1 in range(raw.shape[0]):
        for radial_2 in range(raw.shape[1]):
            nzero = window[radial_1, radial_2, center]
            if nzero == 0.0:
                continue
            matrix = np.empty((len(orders), len(orders)), dtype=np.complex128)
            rhs = np.empty(len(orders), dtype=np.complex128)
            for row, ell in enumerate(orders):
                rhs[row] = raw[radial_1, radial_2, center + ell]/nzero
                for column, order in enumerate(orders):
                    matrix[row, column] = window[
                        radial_1, radial_2, center + ell - order
                    ]/nzero
            solution = np.linalg.solve(matrix, rhs)
            result[:, radial_1, radial_2] = solution[mmax:]
    return result


def assert_scaled_close(actual, expected, relative, absolute):
    error = np.max(np.abs(actual - expected))
    scale = np.max(np.abs(expected))
    if error > absolute + relative*scale:
        raise AssertionError(
            f"maximum error {error:.6g} exceeds "
            f"{absolute + relative*scale:.6g} at scale {scale:.6g}"
        )


def test_octree_ggg_raw_window_and_edge_match_treecorr():
    if treecorr is None:
        try:
            import pytest
        except ImportError as exc:
            raise RuntimeError("TreeCorr is required for this comparison") from exc
        pytest.skip("TreeCorr is not installed")

    rng = np.random.default_rng(9137)
    count = 72
    x = rng.uniform(-0.43, 0.43, count)
    y = rng.uniform(-0.30, 0.44, count)
    kappa = 0.4 + 0.7*x - 0.3*y + 0.2*x*y
    weights = rng.uniform(0.7, 1.3, count)
    require_2d_build(x, y, kappa, weights)

    minimum = 0.05
    maximum = 0.7
    bins = 3
    mmax = 3
    with tempfile.TemporaryDirectory(prefix="ctreeballs-treecorr-ggg-") as tmp:
        root = Path(tmp)
        balls = cballs()
        balls.set(
            {
                "searchMethod": "octree-ggg-omp",
                "rangeN": maximum,
                "rminHist": minimum,
                "sizeHistN": bins,
                "mChebyshev": mmax,
                "lengthBox": 1.2,
                "numberThreads": 1,
                "verbose": 0,
                "verbose_log": 0,
                "rootDir": str(root),
                "options": (
                    "no-normalize-HistZeta,edge-corrections,no-one-ball"
                ),
            }
        )
        balls.set_catalog(
            np.column_stack((x, y)), kappa=kappa, weights=weights
        )
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
            raw = np.asarray(raw)
            edge = np.asarray(edge)
        finally:
            balls.struct_cleanup()

        window = np.asarray(
            [read_window_mode(root, order) for order in range(2*mmax + 1)]
        )

    reference = treecorr.KKKCorrelation(
        min_sep=minimum,
        max_sep=maximum,
        nbins=bins,
        bin_type="LogMultipole",
        max_n=2*mmax,
        brute=True,
        verbose=0,
    )
    reference.process(treecorr.Catalog(x=x, y=y, k=kappa, w=weights))
    center = reference.zeta.shape[2] // 2

    for order in range(mmax + 1):
        assert_scaled_close(
            raw[order], reference.zeta[:, :, center + order],
            relative=5.0e-9, absolute=2.0e-8
        )
    for order in range(2*mmax + 1):
        assert_scaled_close(
            window[order], reference.weight[:, :, center + order],
            relative=5.0e-8, absolute=2.0e-4
        )

    reference_edge = solve_treecorr_edge(
        reference.zeta, reference.weight, mmax
    )
    assert_scaled_close(
        edge, reference_edge, relative=1.0e-7, absolute=2.0e-9
    )


if __name__ == "__main__":
    if treecorr is None:
        raise SystemExit("SKIP: TreeCorr is not installed")
    test_octree_ggg_raw_window_and_edge_match_treecorr()
    print("PASS: octree-GGG raw, window, and edge modes match TreeCorr")
