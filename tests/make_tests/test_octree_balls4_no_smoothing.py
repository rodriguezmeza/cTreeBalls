#!/usr/bin/env python3

import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from cyballs import cballs, search_method_id


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


def parameters(root_dir, threads, search_method="octree-balls4-omp",
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
               search_method="octree-balls4-omp", options=EXACT_OPTIONS,
               mask=None, rsmooth=None):
    balls = cballs()
    params = parameters(root_dir, threads, search_method, options)
    if rsmooth is not None:
        params["rsmooth"] = str(rsmooth)
    balls.set(params)
    balls.set_catalog(positions, kappa=kappa, weights=weights, mask=mask)
    try:
        balls.Run(level=["MainLoop"])
        result = {
            "histNN": balls.getHistNN().copy(),
            "histCF": balls.getHistCF().copy(),
            "histXi2pcf": balls.getHistXi2pcf().copy(),
        }
        if search_method == "octree-balls4-omp":
            for multipole in range(1, 5):
                for component in range(1, 5):
                    name = f"zeta[{multipole},{component}]"
                    result[name] = balls.getHistZetaMsincos(
                        multipole, component
                    ).copy()
        return result
    finally:
        balls.struct_cleanup()
        balls.clear_catalogs()


def run_executable(executable, catalog, root_dir):
    command = [
        executable,
        "search=octree-balls4-omp",
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
        command, check=False, capture_output=True, text=True, timeout=120
    )
    if completed.returncode != 0:
        raise AssertionError(
            "standalone octree-balls4-omp failed:\n"
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
    for rejected_name, options in (
        ("octree-balls4-omp", "KKKCorrelation,smooth-pivot"),
        ("octree-kkk-balls4-omp", "KKKCorrelation"),
    ):
        rejected = [arg for arg in command
                    if not arg.startswith(("search=", "options="))]
        rejected += [f"search={rejected_name}", f"options={options}"]
        completed = subprocess.run(rejected, capture_output=True, text=True, timeout=120)
        assert completed.returncode != 0, completed.stdout
        if rejected_name == "octree-balls4-omp":
            assert "does not support options=smooth-pivot" in completed.stdout + completed.stderr


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
    if not any(np.any(value != 0.0) for name, value in reference.items()
               if name.startswith("zeta[")):
        raise AssertionError("the BALLS4 search returned only zero 3PCF matrices")


def assert_2pcf_matches_ggg(reference, candidate, check_cf=True):
    names = ("histNN", "histCF", "histXi2pcf") if check_cf else ("histNN", "histXi2pcf")
    for name in names:
        if not np.all(np.isfinite(candidate[name])):
            raise AssertionError(f"non-finite BALLS4 values in {name}")
        np.testing.assert_allclose(
            candidate[name], reference[name], rtol=2.0e-13, atol=2.0e-13,
            err_msg=f"exact BALLS4 versus octree-GGG {name}",
        )


def assert_2pcf_matches_pairs(positions, kappa, weights, candidate, mask=None):
    i, j = np.triu_indices(len(positions), 1)
    if mask is not None:
        keep = mask[i] & mask[j]
        i, j = i[keep], j[keep]
    distance = np.linalg.norm(positions[i] - positions[j], axis=1)
    edges = np.geomspace(0.05, 1.5, 9)
    counts = np.histogram(distance, bins=edges)[0]
    field = weights * kappa
    numerator = np.histogram(distance, bins=edges, weights=field[i]*field[j])[0]
    xi = np.divide(numerator, counts, out=np.zeros(8), where=counts > 0)
    np.testing.assert_array_equal(candidate["histNN"], counts)
    np.testing.assert_allclose(candidate["histXi2pcf"], xi, rtol=2e-13, atol=2e-13)


def test_no_smoothing_profile_runs_c_and_cython():
    executable = os.environ.get("CBALLS", str(ROOT / "cballs"))
    assert search_method_id("octree-balls4-omp") == 69
    assert search_method_id("octree-kkk-balls4-omp") == -1

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
        assert_2pcf_matches_pairs(positions, kappa, weights, one_thread)
        # GGG's SMOOTHPIVOT build uses a different weighted estimator.
        # Keep that reference exact and unweighted; check weights directly above.
        unit_weights = np.ones_like(weights)
        unweighted = run_cython(
            positions, kappa, unit_weights, str(root / "cy-unweighted"), 1
        )
        ggg = run_cython(
            positions, kappa, unit_weights, str(root / "cy-ggg"), 1,
            "octree-ggg-omp", rsmooth=0.0,
        )
        assert_2pcf_matches_ggg(ggg, unweighted)
        accelerated_one = run_cython(
            positions, kappa, weights, str(root / "cy-accelerated-one"), 1,
            options=ACCELERATED_OPTIONS,
        )
        accelerated_three = run_cython(
            positions, kappa, weights, str(root / "cy-accelerated-three"), 3,
            options=ACCELERATED_OPTIONS,
        )
        assert_results(accelerated_one, accelerated_three)
        large_radius = run_cython(
            positions, kappa, weights, str(root / "cy-large-rsmooth"), 3,
            rsmooth=10.0,
        )
        assert_results(one_thread, large_radius)
        snapshot = os.environ.get("BALLS4_PROFILE_OUTPUT")
        if snapshot:
            np.savez(snapshot, **one_thread,
                     **{f"accelerated_{key}": value for key, value in accelerated_one.items()})


def test_smooth_pivot_rejection_is_recoverable():
    positions, kappa, weights = unit_sphere_catalog()
    with tempfile.TemporaryDirectory(prefix="ctreeballs-balls4-rejection-") as tmp:
        balls = cballs()
        balls.set_catalog(positions, kappa=kappa, weights=weights)
        try:
            for _ in range(3):
                balls.set(parameters(tmp, 2, options=EXACT_OPTIONS + ",smooth-pivot"))
                try:
                    balls.Run(level=["MainLoop"])
                except Exception as exc:
                    assert "does not support options=smooth-pivot" in str(exc), str(exc)
                else:
                    raise AssertionError("octree-balls4-omp accepted smooth-pivot")
                finally:
                    balls.struct_cleanup()
            balls.set(parameters(tmp, 2))
            balls.Run(level=["MainLoop"])
            assert np.any(balls.getHistXi2pcf() != 0)
        finally:
            balls.struct_cleanup()
            balls.clear_catalogs()


def test_masked_b4_partition_matches_octree_ggg():
    positions, kappa, weights = unit_sphere_catalog(256)
    mask = positions[:, 2] > -0.15
    options = EXACT_OPTIONS + ",read-mask"

    with tempfile.TemporaryDirectory(prefix="ctreeballs-balls4-mask-") as tmp:
        root = Path(tmp)
        unit_weights = np.ones_like(weights)
        ggg = run_cython(
            positions, kappa, unit_weights, str(root / "ggg"), 2,
            search_method="octree-ggg-omp", options=options, mask=mask,
            rsmooth=0.0,
        )
        unweighted = run_cython(
            positions, kappa, unit_weights, str(root / "unweighted"), 2,
            options=options, mask=mask,
        )
        balls4 = run_cython(
            positions, kappa, weights, str(root / "balls4"), 2,
            options=options, mask=mask,
        )

    assert_2pcf_matches_ggg(ggg, unweighted, check_cf=False)
    assert_2pcf_matches_pairs(positions, kappa, weights, balls4, mask=mask)
    # BALLS4 density normalization uses the selected population, not the full map.
    delta = np.log10(1.5/0.05)/8
    expected_cf = np.zeros(8)
    selected = unweighted["histNN"] > 0
    expected_cf[selected] = (unweighted["histNN"][selected] * 8.0 /
        (2*np.pi*delta**3*np.count_nonzero(mask)**2*(np.arange(8)[selected]+0.5)**2)-1)
    np.testing.assert_allclose(unweighted["histCF"], expected_cf, rtol=2e-13, atol=2e-13)
    if not np.any(balls4["histNN"] != 0.0):
        raise AssertionError("masked BALLS4 regression produced no pairs")


if __name__ == "__main__":
    test_no_smoothing_profile_runs_c_and_cython()
    test_masked_b4_partition_matches_octree_ggg()
    test_smooth_pivot_rejection_is_recoverable()
    print(
        "PASS: octree-balls4-omp C/Cython smoke, exact and masked 2PCF oracles, "
        "2PCF/3PCF thread comparisons, and smooth-pivot rejection"
    )
