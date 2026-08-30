#!/usr/bin/env python3
"""Cython in-memory regression for octree-3pcf-3d-omp."""

from pathlib import Path
import tempfile

import numpy as np

from cyballs import cballs
import test_octree_3pcf_3d_survey as survey_oracle


RMIN = 0.05
RMAX = 0.80
NBINS = 4


def read_rows(path):
    return np.asarray([
        [float(token) for token in line.split()]
        for line in path.read_text(encoding="ascii").splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ])


def direct_xi(positions, kappa, weights):
    numerator = np.zeros(NBINS)
    denominator = np.zeros(NBINS)
    scale = NBINS / (RMAX - RMIN)
    for pivot in range(len(positions)):
        for neighbor in range(len(positions)):
            if pivot == neighbor:
                continue
            radius = np.linalg.norm(positions[neighbor] - positions[pivot])
            if radius <= RMIN or radius >= RMAX:
                continue
            radial_bin = int((radius - RMIN) * scale)
            pair_weight = weights[pivot] * weights[neighbor]
            numerator[radial_bin] += (
                pair_weight * kappa[pivot] * kappa[neighbor]
            )
            denominator[radial_bin] += pair_weight
    return numerator, denominator


def run_once(positions, kappa, weights, threads, root):
    balls = cballs()
    balls.set({
        "searchMethod": "octree-3pcf-3d-omp",
        "rangeN": RMAX,
        "rminHist": RMIN,
        "sizeHistN": NBINS,
        "sizeHistPhi": 8,
        "mChebyshev": 3,
        "lengthBox": 1.0,
        "useLogHist": False,
        "usePeriodic": False,
        "numberThreads": threads,
        "scanLevel": 2,
        "verbose": 0,
        "verbose_log": 0,
        "rootDir": str(root),
        "options": "compute-2pcf-3d,compute-3pcf-3d",
    })
    balls.set_catalog(positions, kappa=kappa, weights=weights)
    try:
        balls.Run(level=["MainLoop"])
        if balls.getNBody() != len(positions):
            raise AssertionError("Cython did not publish the memory-catalog size")
        xi = read_rows(root / "histXi2pcf_3d.txt")
        zeta = read_rows(root / "histZetaM_3d.txt")
    finally:
        balls.struct_cleanup()
    return xi, zeta


def run_survey_once(threads, root):
    data = np.asarray(survey_oracle.DATA, dtype=np.float64)
    randoms = np.asarray(survey_oracle.RANDOMS, dtype=np.float64)
    balls = cballs()
    balls.set({
        "searchMethod": "octree-3pcf-3d-omp",
        "rangeN": survey_oracle.RMAX,
        "rminHist": survey_oracle.RMIN,
        "sizeHistN": survey_oracle.NBINS,
        "sizeHistPhi": 8,
        "mChebyshev": survey_oracle.LMAX,
        "lengthBox": 2.0,
        "useLogHist": False,
        "usePeriodic": False,
        "numberThreads": threads,
        "scanLevel": 2,
        "verbose": 0,
        "verbose_log": 0,
        "rootDir": str(root),
        "iCatalogs": "1,2",
        "options": (
            "survey-estimator-3d,compute-2pcf-3d,compute-3pcf-3d"
        ),
    })
    balls.set_catalog(data[:, :3], weights=data[:, 3], catalog=0)
    balls.set_catalog(randoms[:, :3], weights=randoms[:, 3], catalog=1)
    try:
        balls.Run(level=["MainLoop"])
        survey_oracle.check_oracle(root)
    finally:
        balls.struct_cleanup()


def main():
    rng = np.random.default_rng(80317)
    positions = rng.uniform(-0.4, 0.4, size=(32, 3))
    kappa = rng.normal(size=32)
    weights = rng.uniform(0.5, 1.5, size=32)
    expected_num, expected_den = direct_xi(positions, kappa, weights)

    with tempfile.TemporaryDirectory(prefix="ctreeballs-3d-cython-") as tmp:
        tmp = Path(tmp)
        one_root = tmp / "one-thread"
        many_root = tmp / "three-threads"
        one_root.mkdir()
        many_root.mkdir()
        one_xi, one_zeta = run_once(
            positions, kappa, weights, 1, one_root
        )
        many_xi, many_zeta = run_once(
            positions, kappa, weights, 3, many_root
        )

        np.testing.assert_allclose(one_xi[:, 3], expected_num,
                                   rtol=3.0e-11, atol=3.0e-11)
        np.testing.assert_allclose(one_xi[:, 4], expected_den,
                                   rtol=3.0e-11, atol=3.0e-11)
        expected_xi = np.divide(
            expected_num, expected_den,
            out=np.zeros_like(expected_num), where=expected_den != 0,
        )
        np.testing.assert_allclose(one_xi[:, 2], expected_xi,
                                   rtol=3.0e-11, atol=3.0e-11)
        np.testing.assert_allclose(many_xi, one_xi,
                                   rtol=3.0e-11, atol=3.0e-11)
        np.testing.assert_allclose(many_zeta, one_zeta,
                                   rtol=3.0e-11, atol=3.0e-11)

        survey_one_root = tmp / "survey-one-thread"
        survey_many_root = tmp / "survey-three-threads"
        survey_one_root.mkdir()
        survey_many_root.mkdir()
        run_survey_once(1, survey_one_root)
        run_survey_once(3, survey_many_root)
        survey_oracle.compare_outputs(survey_one_root, survey_many_root)

    print(
        "PASS: cyballs in-memory legacy/survey 3D estimators and "
        "OpenMP comparison"
    )


if __name__ == "__main__":
    main()
