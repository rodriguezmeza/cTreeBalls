#!/usr/bin/env python3
"""Oracle and mode regression for octree-3pcf-3d-omp."""

import argparse
import math
from pathlib import Path
import subprocess
import tempfile


CATALOG = (
    (0.00, 0.00, 0.00, 0.50, 1.00),
    (0.41, 0.13, 0.22, -0.70, 2.00),
    (-0.31, 0.53, 0.11, 1.20, 0.80),
    (0.23, -0.43, 0.61, 0.90, 1.50),
    (-0.52, -0.21, -0.34, -1.10, 0.60),
    (0.72, 0.39, -0.12, 0.30, 1.10),
)
RMIN = 0.05
RMAX = 1.60
NBINS = 4
LMAX = 3


def legendre(ell, value):
    if ell == 0:
        return 1.0
    if ell == 1:
        return value
    previous = 1.0
    current = value
    for order in range(2, ell + 1):
        following = ((2 * order - 1) * value * current
                     - (order - 1) * previous) / order
        previous, current = current, following
    return current


def separation(left, right):
    delta = tuple(right[index] - left[index] for index in range(3))
    radius = math.sqrt(sum(value * value for value in delta))
    return radius, tuple(value / radius for value in delta)


def radial_bin(radius):
    if radius <= RMIN or radius >= RMAX:
        return None
    index = int((radius - RMIN) * NBINS / (RMAX - RMIN))
    return index if 0 <= index < NBINS else None


def direct_oracle(los_ids=None):
    xi_num = [0.0] * NBINS
    xi_den = [0.0] * NBINS
    zeta_num = [[[0.0] * NBINS for _ in range(NBINS)]
                for _ in range(LMAX + 1)]
    zeta_den = [[[0.0] * NBINS for _ in range(NBINS)]
                for _ in range(LMAX + 1)]

    for pivot_index, pivot in enumerate(CATALOG):
        pivot_field = pivot[3]
        pivot_weight = pivot[4]
        neighbors = [[] for _ in range(NBINS)]
        for neighbor_index, neighbor in enumerate(CATALOG):
            if neighbor_index == pivot_index:
                continue
            if los_ids is not None and los_ids[neighbor_index] == los_ids[pivot_index]:
                continue
            radius, direction = separation(pivot, neighbor)
            bin_index = radial_bin(radius)
            if bin_index is None:
                continue
            neighbor_field = neighbor[3]
            neighbor_weight = neighbor[4]
            xi_num[bin_index] += (
                pivot_weight * pivot_field * neighbor_weight * neighbor_field
            )
            xi_den[bin_index] += pivot_weight * neighbor_weight
            neighbors[bin_index].append(
                (neighbor_index, neighbor_weight, neighbor_field, direction)
            )

        for ell in range(LMAX + 1):
            for first_bin in range(NBINS):
                for second_bin in range(NBINS):
                    for first in neighbors[first_bin]:
                        for second in neighbors[second_bin]:
                            if first[0] == second[0]:
                                continue
                            cosine = sum(a * b for a, b in zip(first[3], second[3]))
                            cosine = min(1.0, max(-1.0, cosine))
                            zeta_num[ell][first_bin][second_bin] += (
                                pivot_weight * pivot_field
                                * first[1] * first[2]
                                * second[1] * second[2]
                                * (2 * ell + 1) * legendre(ell, cosine)
                            )
                            zeta_den[ell][first_bin][second_bin] += (
                                pivot_weight * first[1] * second[1]
                            )
    return xi_num, xi_den, zeta_num, zeta_den


def read_rows(path):
    return [
        [float(token) for token in line.split()]
        for line in path.read_text(encoding="ascii").splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]


def run_case(executable, catalog, root, options, threads, expect_success=True,
             infile_format="multi-columns-ascii", columns="1,2,3,4,5",
             theta=1.25):
    root.mkdir()
    command = [
        str(executable),
        "search=octree-3pcf-3d-omp",
        f"in={catalog}",
        f"infmt={infile_format}",
        f"columns={columns}",
        "iCatalogs=1",
        f"rootDir={root}",
        f"numberThreads={threads}",
        "usePeriodic=false",
        "useLogHist=false",
        f"rangeN={RMAX}",
        f"rminHist={RMIN}",
        f"sizeHistN={NBINS}",
        f"mChebyshev={LMAX}",
        f"theta={theta}",
        "stepState=1000000",
        "verbose=0",
        "verbose_log=0",
        f"options={options}",
    ]
    result = subprocess.run(command, text=True, capture_output=True, check=False)
    (root / "run.log").write_text(result.stdout + result.stderr, encoding="utf-8")
    if expect_success and result.returncode != 0:
        raise AssertionError(f"cballs failed:\n{result.stdout}\n{result.stderr}")
    if not expect_success and result.returncode == 0:
        raise AssertionError("cballs accepted incompatible 3D mode options")
    return result


def assert_close(actual, expected, label, tolerance=3.0e-6):
    if not math.isclose(actual, expected, rel_tol=tolerance, abs_tol=tolerance):
        raise AssertionError(f"{label}: {actual:.17g} != {expected:.17g}")


def check_oracle(root):
    xi_num, xi_den, zeta_num, zeta_den = direct_oracle()
    xi_rows = read_rows(root / "histXi2pcf_3d.txt")
    zeta_rows = read_rows(root / "histZetaM_3d.txt")
    if len(xi_rows) != NBINS:
        raise AssertionError("unexpected 2PCF row count")
    if len(zeta_rows) != (LMAX + 1) * NBINS * NBINS:
        raise AssertionError("unexpected 3PCF row count")

    for row in xi_rows:
        bin_index = int(row[0]) - 1
        assert_close(row[3], xi_num[bin_index], f"xi numerator bin {bin_index}")
        assert_close(row[4], xi_den[bin_index], f"xi denominator bin {bin_index}")
        expected = xi_num[bin_index] / xi_den[bin_index] if xi_den[bin_index] else 0.0
        assert_close(row[2], expected, f"xi bin {bin_index}")

    for row in zeta_rows:
        ell = int(row[0])
        first_bin = int(row[1]) - 1
        second_bin = int(row[2]) - 1
        expected_num = zeta_num[ell][first_bin][second_bin]
        expected_den = zeta_den[ell][first_bin][second_bin]
        assert_close(row[6], expected_num,
                     f"zeta numerator ({ell},{first_bin},{second_bin})")
        assert_close(row[7], expected_den,
                     f"zeta denominator ({ell},{first_bin},{second_bin})")
        expected = expected_num / expected_den if expected_den else 0.0
        assert_close(row[5], expected,
                     f"zeta ({ell},{first_bin},{second_bin})")


def compare_outputs(left, right):
    for filename in ("histXi2pcf_3d.txt", "histZetaM_3d.txt"):
        left_rows = read_rows(left / filename)
        right_rows = read_rows(right / filename)
        if len(left_rows) != len(right_rows):
            raise AssertionError(f"{filename}: thread-count row mismatch")
        for row_index, (left_row, right_row) in enumerate(zip(left_rows, right_rows)):
            for column, (left_value, right_value) in enumerate(zip(left_row, right_row)):
                assert_close(left_value, right_value,
                             f"{filename}:{row_index}:{column}", tolerance=2.0e-11)


def check_fits_los(executable, tmp):
    try:
        import numpy as np
        from astropy.io import fits
    except ImportError:
        print("SKIP: FITS LOS regression requires NumPy and Astropy")
        return

    los_ids = (10, 10, 20, 20, 30, 30)
    table = np.array(
        [tuple(row) + (los_id,) for row, los_id in zip(CATALOG, los_ids)],
        dtype=[("X", "f8"), ("Y", "f8"), ("Z", "f8"),
               ("DELTA", "f8"), ("WEIGHT", "f8"), ("LOS_ID", "i8")],
    )
    catalog = tmp / "weighted_xyz_los.fits"
    fits.BinTableHDU(table).writeto(catalog)
    root = tmp / "fits-los"
    run_case(
        executable, catalog, root,
        "with-weight,exclude-same-los,only-2pcf-3d", 2,
        infile_format="fits", columns="1,2,3,4,5,6",
    )
    expected_num, expected_den, _, _ = direct_oracle(los_ids)
    for row in read_rows(root / "histXi2pcf_3d.txt"):
        bin_index = int(row[0]) - 1
        assert_close(row[3], expected_num[bin_index],
                     f"FITS LOS xi numerator bin {bin_index}")
        assert_close(row[4], expected_den[bin_index],
                     f"FITS LOS xi denominator bin {bin_index}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--cballs", type=Path,
        default=Path(__file__).resolve().parents[2] / "cballs",
    )
    args = parser.parse_args()
    executable = args.cballs.resolve()
    if not executable.is_file():
        raise SystemExit(f"cballs executable not found: {executable}")

    with tempfile.TemporaryDirectory(prefix="ctreeballs-3d-") as tmp_name:
        tmp = Path(tmp_name)
        catalog = tmp / "weighted_xyz.txt"
        catalog.write_text(
            "".join(" ".join(str(value) for value in row) + "\n" for row in CATALOG),
            encoding="ascii",
        )
        combined = "pos-and-convergence-weight,compute-2pcf-3d,compute-3pcf-3d"
        one_thread = tmp / "combined-1"
        three_threads = tmp / "combined-3"
        run_case(executable, catalog, one_thread, combined, 1)
        run_case(executable, catalog, three_threads, combined, 3)
        check_oracle(one_thread)
        compare_outputs(one_thread, three_threads)

        divided_global_cutoff = tmp / "rcut-over-theta"
        run_case(
            executable, catalog, divided_global_cutoff,
            f"{combined},Rcut/theta", 2, theta=1.25,
        )
        check_oracle(divided_global_cutoff)
        compare_outputs(one_thread, divided_global_cutoff)

        only_xi = tmp / "only-xi"
        run_case(executable, catalog, only_xi,
                 "pos-and-convergence-weight,only-2pcf-3d", 2)
        if not (only_xi / "histXi2pcf_3d.txt").is_file():
            raise AssertionError("only-2pcf-3d did not write its output")
        if (only_xi / "histZetaM_3d.txt").exists():
            raise AssertionError("only-2pcf-3d unexpectedly wrote 3PCF output")

        default_mode = tmp / "default"
        run_case(executable, catalog, default_mode,
                 "pos-and-convergence-weight", 2)
        if not (default_mode / "histZetaM_3d.txt").is_file():
            raise AssertionError("default mode did not write 3PCF output")
        if (default_mode / "histXi2pcf_3d.txt").exists():
            raise AssertionError("default mode unexpectedly wrote 2PCF output")

        run_case(executable, catalog, tmp / "conflict",
                 "pos-and-convergence-weight,only-2pcf-3d,only-3pcf-3d",
                 1, expect_success=False)
        check_fits_los(executable, tmp)

    print("PASS: 3D scalar 2PCF/3PCF oracle, Rcut/theta, modes, and OpenMP comparison")


if __name__ == "__main__":
    main()
