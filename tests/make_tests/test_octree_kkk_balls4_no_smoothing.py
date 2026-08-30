#!/usr/bin/env python3

import os
import subprocess
import tempfile
from pathlib import Path

import numpy as np

from cyballs import cballs


def unit_sphere_catalog(count=48):
    index = np.arange(count, dtype=np.float64)
    z = 1.0 - 2.0*(index + 0.5)/count
    phi = index*np.pi*(3.0 - np.sqrt(5.0))
    radius = np.sqrt(1.0 - z*z)
    positions = np.column_stack((radius*np.cos(phi), radius*np.sin(phi), z))
    kappa = 1.0 + 0.2*positions[:, 0] - 0.1*positions[:, 1]
    weights = 0.8 + 0.15*(1.0 + positions[:, 2])
    return positions, kappa, weights


def write_catalog(path, positions, kappa):
    with path.open("w", encoding="ascii") as stream:
        stream.write("# nbody NDIM Lx Ly Lz\n")
        stream.write(f"# {len(positions)} 3 2.0 2.0 2.0\n")
        for point, field in zip(positions, kappa):
            stream.write(
                "{:.17g} {:.17g} {:.17g} {:.17g}\n".format(
                    point[0], point[1], point[2], field
                )
            )


EXACT_OPTIONS = (
    "no-out-Hist,no-normalize-HistZeta,KKKCorrelation,"
    "no-one-ball,no-two-balls,compute-HistN,and-CF"
)
ACCELERATED_OPTIONS = (
    "no-out-Hist,no-normalize-HistZeta,KKKCorrelation,compute-HistN,and-CF"
)


def parameters(root_dir, threads, search_method="octree-kkk-balls4-omp",
               options=EXACT_OPTIONS):
    return {
        "searchMethod": search_method,
        "rangeN": 1.5,
        "rminHist": 0.05,
        "sizeHistN": 8,
        "mChebyshev": 3,
        "lengthBox": 2.0,
        "numberThreads": threads,
        "useLogHist": True,
        "usePeriodic": False,
        "verbose": 0,
        "verbose_log": 0,
        "rootDir": root_dir,
        "options": options,
    }


def run_cython(positions, kappa, weights, root_dir, threads,
               search_method="octree-kkk-balls4-omp", options=EXACT_OPTIONS,
               mask=None):
    balls = cballs()
    balls.set(parameters(root_dir, threads, search_method, options))
    balls.set_catalog(positions, kappa=kappa, weights=weights, mask=mask)
    try:
        balls.Run(level=["MainLoop"])
        result = {
            "histNN": balls.getHistNN().copy(),
            "histCF": balls.getHistCF().copy(),
            "histXi2pcf": balls.getHistXi2pcf().copy(),
        }
        if search_method == "octree-kkk-balls4-omp":
            for multipole in range(1, 5):
                for component in range(1, 5):
                    name = f"zeta[{multipole},{component}]"
                    result[name] = balls.getHistZetaMsincos(
                        multipole, component
                    ).copy()
        return result
    finally:
        balls.struct_cleanup()


def run_executable(executable, catalog, root_dir):
    command = [
        executable,
        "search=octree-kkk-balls4-omp",
        f"in={catalog}",
        "infmt=columns-ascii",
        f"rootDir={root_dir}",
        "rangeN=1.5",
        "rminHist=0.05",
        "sizeHistN=8",
        "mChebyshev=3",
        "lengthBox=2.0",
        "nthreads=2",
        "useLogHist=true",
        "usePeriodic=false",
        "verb=0",
        "verblog=0",
        "options=no-normalize-HistZeta,KKKCorrelation,no-one-ball,"
        "no-two-balls,compute-HistN,and-CF",
    ]
    completed = subprocess.run(
        command, check=False, capture_output=True, text=True
    )
    if completed.returncode != 0:
        raise AssertionError(
            "standalone octree-kkk-balls4-omp failed:\n"
            + completed.stdout
            + completed.stderr
        )
    for filename in ("histNN.txt", "histCF.txt", "histXi2pcf.txt"):
        output = root_dir / filename
        if not output.is_file():
            raise AssertionError(f"standalone output is missing {output}")
        values = np.loadtxt(output)
        if values.shape != (8, 2) or not np.all(np.isfinite(values)):
            raise AssertionError(f"invalid standalone 2PCF output in {output}")


def assert_results(reference, candidate):
    nonzero = False
    for name in reference:
        if not np.all(np.isfinite(reference[name])):
            raise AssertionError(f"non-finite values in {name}")
        nonzero = nonzero or np.any(reference[name] != 0.0)
        np.testing.assert_allclose(
            candidate[name], reference[name], rtol=2.0e-13, atol=2.0e-13,
            err_msg=name,
        )
    if not nonzero:
        raise AssertionError("the BALLS4 search returned only zero 3PCF matrices")


def assert_2pcf_matches_ggg(reference, candidate):
    for name in ("histNN", "histCF", "histXi2pcf"):
        if not np.all(np.isfinite(candidate[name])):
            raise AssertionError(f"non-finite BALLS4 values in {name}")
        np.testing.assert_allclose(
            candidate[name], reference[name], rtol=2.0e-13, atol=2.0e-13,
            err_msg=f"exact BALLS4 versus octree-GGG {name}",
        )


def test_no_smoothing_profile_runs_c_and_cython():
    executable = os.environ.get("CBALLS")
    if not executable:
        raise AssertionError("CBALLS must name the standalone executable")

    positions, kappa, weights = unit_sphere_catalog()
    with tempfile.TemporaryDirectory(prefix="ctreeballs-balls4-normal-tree-") as tmp:
        root = Path(tmp)
        catalog = root / "unit_sphere.txt"
        write_catalog(catalog, positions, kappa)
        run_executable(executable, catalog, root / "c-output")
        one_thread = run_cython(
            positions, kappa, weights, str(root / "cy-one"), 1
        )
        three_threads = run_cython(
            positions, kappa, weights, str(root / "cy-three"), 3
        )
        assert_results(one_thread, three_threads)
        ggg = run_cython(
            positions, kappa, weights, str(root / "cy-ggg"), 1,
            "octree-ggg-omp"
        )
        assert_2pcf_matches_ggg(ggg, one_thread)
        accelerated_one = run_cython(
            positions, kappa, weights, str(root / "cy-accelerated-one"), 1,
            options=ACCELERATED_OPTIONS,
        )
        accelerated_three = run_cython(
            positions, kappa, weights, str(root / "cy-accelerated-three"), 3,
            options=ACCELERATED_OPTIONS,
        )
        assert_results(accelerated_one, accelerated_three)


def test_masked_b4_partition_matches_octree_ggg():
    positions, kappa, weights = unit_sphere_catalog(256)
    mask = positions[:, 2] > -0.15
    options = EXACT_OPTIONS + ",read-mask"

    with tempfile.TemporaryDirectory(prefix="ctreeballs-balls4-mask-") as tmp:
        root = Path(tmp)
        ggg = run_cython(
            positions, kappa, weights, str(root / "ggg"), 2,
            search_method="octree-ggg-omp", options=options, mask=mask,
        )
        balls4 = run_cython(
            positions, kappa, weights, str(root / "balls4"), 2,
            options=options, mask=mask,
        )

    assert_2pcf_matches_ggg(ggg, balls4)
    if not np.any(balls4["histNN"] != 0.0):
        raise AssertionError("masked BALLS4 regression produced no pairs")


if __name__ == "__main__":
    test_no_smoothing_profile_runs_c_and_cython()
    test_masked_b4_partition_matches_octree_ggg()
    print(
        "PASS: octree-kkk-balls4-omp 2PCF matches exact octree-GGG "
        "in C, cyballs, and masked B4 partitions"
    )
