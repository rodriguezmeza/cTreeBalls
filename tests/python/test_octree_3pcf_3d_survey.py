#!/usr/bin/env python3
"""Direct survey-estimator oracle for octree-3pcf-3d-omp."""

import argparse
import math
from pathlib import Path
import subprocess
import tempfile


RMIN = 0.08
RMAX = 1.45
NBINS = 3
LMAX = 3
FOURPI = 4.0 * math.pi

DATA = (
    (0.10, 0.16, 0.22, 1.00),
    (0.34, 0.24, 0.18, 0.80),
    (0.63, 0.13, 0.31, 1.20),
    (0.20, 0.57, 0.28, 0.90),
    (0.51, 0.48, 0.12, 1.10),
    (0.82, 0.39, 0.41, 0.70),
    (0.17, 0.29, 0.78, 1.30),
    (0.47, 0.66, 0.69, 0.85),
    (0.78, 0.72, 0.56, 1.05),
    (1.03, 0.52, 0.27, 0.95),
)

RANDOMS = (
    (0.04, 0.09, 0.12, 0.70),
    (0.29, 0.11, 0.36, 1.10),
    (0.57, 0.07, 0.15, 0.90),
    (0.91, 0.19, 0.33, 1.20),
    (1.16, 0.08, 0.52, 0.80),
    (0.13, 0.43, 0.19, 1.05),
    (0.42, 0.39, 0.47, 0.95),
    (0.73, 0.51, 0.23, 1.15),
    (1.08, 0.46, 0.61, 0.75),
    (0.07, 0.76, 0.38, 1.25),
    (0.36, 0.81, 0.72, 0.85),
    (0.69, 0.87, 0.44, 1.00),
    (1.01, 0.79, 0.17, 0.90),
    (0.18, 0.31, 0.94, 1.10),
    (0.55, 0.58, 1.03, 0.80),
    (0.89, 0.68, 0.88, 1.20),
    (1.19, 0.37, 0.91, 0.95),
    (1.24, 0.91, 0.63, 1.05),
)


def legendre(ell, value):
    if ell == 0:
        return 1.0
    if ell == 1:
        return value
    previous = 1.0
    current = value
    for order in range(2, ell + 1):
        following = (
            (2 * order - 1) * value * current
            - (order - 1) * previous
        ) / order
        previous, current = current, following
    return current


def radial_bin(radius):
    if radius <= RMIN or radius >= RMAX:
        return None
    index = int((radius - RMIN) * NBINS / (RMAX - RMIN))
    return index if 0 <= index < NBINS else None


def neighbors(catalog, pivot_index):
    pivot = catalog[pivot_index]
    shells = [[] for _ in range(NBINS)]
    for index, point in enumerate(catalog):
        if index == pivot_index:
            continue
        vector = tuple(point[axis] - pivot[axis] for axis in range(3))
        radius = math.sqrt(sum(value * value for value in vector))
        bin_index = radial_bin(radius)
        if bin_index is None:
            continue
        shells[bin_index].append(
            (index, point[3], tuple(value / radius for value in vector))
        )
    return shells


def direct_counts(catalog):
    xi = [0.0] * NBINS
    multipoles = [
        [[0.0] * NBINS for _ in range(NBINS)]
        for _ in range(LMAX + 1)
    ]
    for pivot_index, pivot in enumerate(catalog):
        shells = neighbors(catalog, pivot_index)
        for bin_index, shell in enumerate(shells):
            for _, weight, _ in shell:
                xi[bin_index] += pivot[3] * weight
        for ell in range(LMAX + 1):
            for first_bin in range(NBINS):
                for second_bin in range(NBINS):
                    for first in shells[first_bin]:
                        for second in shells[second_bin]:
                            if first[0] == second[0]:
                                continue
                            cosine = sum(
                                left * right
                                for left, right in zip(first[2], second[2])
                            )
                            cosine = min(1.0, max(-1.0, cosine))
                            multipoles[ell][first_bin][second_bin] += (
                                pivot[3] * first[1] * second[1]
                                * (2 * ell + 1) * legendre(ell, cosine)
                            )
    return xi, multipoles


def wigner_zero_squared(left, middle, right):
    total = left + middle + right
    if (
        total % 2
        or left + middle < right
        or left + right < middle
        or middle + right < left
    ):
        return 0.0
    half = total // 2
    log_value = (
        2.0 * math.lgamma(half + 1.0)
        - 2.0 * math.lgamma(half - left + 1.0)
        - 2.0 * math.lgamma(half - middle + 1.0)
        - 2.0 * math.lgamma(half - right + 1.0)
        + math.lgamma(total - 2 * left + 1.0)
        + math.lgamma(total - 2 * middle + 1.0)
        + math.lgamma(total - 2 * right + 1.0)
        - math.lgamma(total + 2.0)
    )
    return math.exp(log_value)


def solve(matrix, vector):
    size = len(vector)
    work = [row[:] for row in matrix]
    result = vector[:]
    for column in range(size):
        pivot = max(range(column, size), key=lambda row: abs(work[row][column]))
        if abs(work[pivot][column]) < 1.0e-14:
            raise AssertionError("survey oracle generated a singular window")
        work[column], work[pivot] = work[pivot], work[column]
        result[column], result[pivot] = result[pivot], result[column]
        for row in range(column + 1, size):
            factor = work[row][column] / work[column][column]
            for inner in range(column + 1, size):
                work[row][inner] -= factor * work[column][inner]
            result[row] -= factor * result[column]
    for row in range(size - 1, -1, -1):
        result[row] = (
            result[row]
            - sum(work[row][column] * result[column]
                  for column in range(row + 1, size))
        ) / work[row][row]
    return result


def corrected_oracle(numerator, randoms, first_bin, second_bin):
    n_encore = []
    r_encore = []
    for ell in range(LMAX + 1):
        sign = -1.0 if ell % 2 else 1.0
        conversion = FOURPI * sign * math.sqrt(2 * ell + 1)
        n_encore.append(numerator[ell][first_bin][second_bin] / conversion)
        r_encore.append(randoms[ell][first_bin][second_bin] / conversion)
    window = [value / r_encore[0] for value in r_encore]
    matrix = []
    for signal in range(LMAX + 1):
        row = []
        for output in range(LMAX + 1):
            coupling = 0.0
            for window_ell in range(LMAX + 1):
                coupling += (
                    window[window_ell]
                    * wigner_zero_squared(signal, window_ell, output)
                    * math.sqrt(
                        (2 * signal + 1)
                        * (2 * window_ell + 1)
                        * (2 * output + 1)
                    )
                    / FOURPI
                )
            row.append(coupling)
        matrix.append(row)
    encore = solve(
        matrix,
        [value / r_encore[0] for value in n_encore],
    )
    legendre_values = [
        value * (-1.0 if ell % 2 else 1.0)
        * math.sqrt(2 * ell + 1) / FOURPI
        for ell, value in enumerate(encore)
    ]
    return encore, legendre_values, n_encore, r_encore


def read_rows(path):
    return [
        [float(token) for token in line.split()]
        for line in path.read_text(encoding="ascii").splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]


def write_catalog(path, rows, weight_scale=1.0):
    path.write_text(
        "".join(
            f"{x:.17g} {y:.17g} {z:.17g} 1 "
            f"{weight * weight_scale:.17g}\n"
            for x, y, z, weight in rows
        ),
        encoding="ascii",
    )


def run_case(executable, data_file, random_file, root, threads):
    root.mkdir()
    command = [
        str(executable),
        "search=octree-3pcf-3d-omp",
        f"in={data_file},{random_file}",
        "infmt=multi-columns-ascii,multi-columns-ascii",
        "columns=1,2,3,4,5",
        "iCatalogs=1,2",
        f"rootDir={root}",
        f"numberThreads={threads}",
        "usePeriodic=false",
        "useLogHist=false",
        f"rangeN={RMAX}",
        f"rminHist={RMIN}",
        f"sizeHistN={NBINS}",
        f"mChebyshev={LMAX}",
        "theta=1.25",
        "stepState=1000000",
        "verbose=0",
        "verbose_log=0",
        "options=pos-and-convergence-weight,survey-estimator-3d,"
        "compute-2pcf-3d,compute-3pcf-3d",
    ]
    result = subprocess.run(command, text=True, capture_output=True, check=False)
    (root / "run.log").write_text(
        result.stdout + result.stderr, encoding="utf-8"
    )
    if result.returncode != 0:
        raise AssertionError(f"cballs survey run failed:\n{result.stdout}\n{result.stderr}")


def assert_close(actual, expected, label, tolerance=4.0e-10):
    if not math.isclose(actual, expected, rel_tol=tolerance, abs_tol=tolerance):
        raise AssertionError(f"{label}: {actual:.17g} != {expected:.17g}")


def compare_outputs(left, right):
    for filename in ("histXi2pcf_3d_survey.txt", "histZetaM_3d_survey.txt"):
        left_rows = read_rows(left / filename)
        right_rows = read_rows(right / filename)
        if len(left_rows) != len(right_rows):
            raise AssertionError(f"{filename}: thread-count row mismatch")
        for row_index, (left_row, right_row) in enumerate(zip(left_rows, right_rows)):
            for column, (left_value, right_value) in enumerate(
                zip(left_row, right_row)
            ):
                assert_close(
                    left_value, right_value,
                    f"{filename}:{row_index}:{column}",
                    tolerance=2.0e-10,
                )


def compare_scale_invariance(reference, scaled):
    reference_xi = read_rows(reference / "histXi2pcf_3d_survey.txt")
    scaled_xi = read_rows(scaled / "histXi2pcf_3d_survey.txt")
    for index, (reference_row, scaled_row) in enumerate(
        zip(reference_xi, scaled_xi)
    ):
        assert_close(reference_row[2], scaled_row[2], f"scaled xi {index}")
        if int(scaled_row[5]) != 1:
            raise AssertionError("scaled survey xi was marked invalid")

    reference_zeta = read_rows(reference / "histZetaM_3d_survey.txt")
    scaled_zeta = read_rows(scaled / "histZetaM_3d_survey.txt")
    for index, (reference_row, scaled_row) in enumerate(
        zip(reference_zeta, scaled_zeta)
    ):
        assert_close(
            reference_row[5], scaled_row[5], f"scaled Legendre zeta {index}"
        )
        assert_close(
            reference_row[6], scaled_row[6], f"scaled ENCORE zeta {index}"
        )
        if int(scaled_row[13]) != 1:
            raise AssertionError("scaled survey 3PCF was marked invalid")


def check_oracle(root, weight_scale=1.0):
    data = tuple(
        (x, y, z, weight * weight_scale) for x, y, z, weight in DATA
    )
    random_catalog = tuple(
        (x, y, z, weight * weight_scale)
        for x, y, z, weight in RANDOMS
    )
    alpha = sum(point[3] for point in data) / sum(
        point[3] for point in random_catalog
    )
    combined = data + tuple(
        (x, y, z, -alpha * weight)
        for x, y, z, weight in random_catalog
    )
    normalized_randoms = tuple(
        (x, y, z, alpha * weight)
        for x, y, z, weight in random_catalog
    )
    numerator_xi, numerator_zeta = direct_counts(combined)
    random_xi, random_zeta = direct_counts(normalized_randoms)

    xi_rows = read_rows(root / "histXi2pcf_3d_survey.txt")
    zeta_rows = read_rows(root / "histZetaM_3d_survey.txt")
    if len(xi_rows) != NBINS:
        raise AssertionError("unexpected survey 2PCF row count")
    if len(zeta_rows) != LMAX * NBINS * (NBINS - 1) // 2:
        raise AssertionError("unexpected survey 3PCF row count")

    for row in xi_rows:
        radial_bin_index = int(row[0]) - 1
        expected = (
            numerator_xi[radial_bin_index] / random_xi[radial_bin_index]
            if random_xi[radial_bin_index] else 0.0
        )
        assert_close(row[2], expected, f"survey xi bin {radial_bin_index}")
        assert_close(row[3], numerator_xi[radial_bin_index], "survey xi N")
        assert_close(row[4], random_xi[radial_bin_index], "survey xi R")
        if int(row[5]) != 1:
            raise AssertionError("survey xi unexpectedly marked invalid")

    expected_by_pair = {}
    for first_bin in range(NBINS):
        for second_bin in range(first_bin + 1, NBINS):
            expected_by_pair[first_bin, second_bin] = corrected_oracle(
                numerator_zeta, random_zeta, first_bin, second_bin
            )
    for row in zeta_rows:
        ell = int(row[0])
        first_bin = int(row[1]) - 1
        second_bin = int(row[2]) - 1
        encore, legendre_values, n_encore, r_encore = expected_by_pair[
            first_bin, second_bin
        ]
        assert_close(row[5], legendre_values[ell], "survey zeta Legendre")
        assert_close(row[6], encore[ell], "survey zeta ENCORE")
        assert_close(
            row[7], numerator_zeta[ell][first_bin][second_bin],
            "survey raw N",
        )
        assert_close(
            row[8], random_zeta[ell][first_bin][second_bin],
            "survey raw R",
        )
        assert_close(row[9], n_encore[ell], "survey ENCORE N")
        assert_close(row[10], r_encore[ell], "survey ENCORE R")
        assert_close(row[11], r_encore[0], "survey R0")
        if int(row[13]) != 1:
            raise AssertionError("survey 3PCF unexpectedly marked invalid")


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

    with tempfile.TemporaryDirectory(prefix="ctreeballs-3d-survey-") as tmp_name:
        tmp = Path(tmp_name)
        data_file = tmp / "data.txt"
        random_file = tmp / "random.txt"
        write_catalog(data_file, DATA)
        write_catalog(random_file, RANDOMS)
        one_thread = tmp / "one-thread"
        three_threads = tmp / "three-threads"
        run_case(executable, data_file, random_file, one_thread, 1)
        run_case(executable, data_file, random_file, three_threads, 3)
        check_oracle(one_thread)
        compare_outputs(one_thread, three_threads)

        scaled_data_file = tmp / "data-scaled.txt"
        scaled_random_file = tmp / "random-scaled.txt"
        weight_scale = 1.0e-8
        write_catalog(scaled_data_file, DATA, weight_scale)
        write_catalog(scaled_random_file, RANDOMS, weight_scale)
        scaled = tmp / "scaled"
        run_case(executable, scaled_data_file, scaled_random_file, scaled, 2)
        check_oracle(scaled, weight_scale)
        compare_scale_invariance(one_thread, scaled)

        invalid_random_file = tmp / "random-negative-weight.txt"
        invalid_randoms = list(RANDOMS)
        invalid_randoms[0] = (*invalid_randoms[0][:3], -invalid_randoms[0][3])
        write_catalog(invalid_random_file, invalid_randoms)
        invalid = tmp / "invalid-negative-weight"
        try:
            run_case(executable, data_file, invalid_random_file, invalid, 1)
        except AssertionError:
            log = (invalid / "run.log").read_text(encoding="utf-8")
            if "negative or non-finite" not in log:
                raise AssertionError(
                    "negative-weight startup failed for the wrong reason"
                )
        else:
            raise AssertionError("negative random weight was accepted")

    print("PASS: ENCORE-style 3D survey 2PCF/3PCF oracle and OpenMP comparison")


if __name__ == "__main__":
    main()
