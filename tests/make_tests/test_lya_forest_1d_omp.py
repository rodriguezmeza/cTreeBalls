#!/usr/bin/env python3
"""Oracle, transverse-invariance, and OpenMP tests for radial Ly-alpha."""

from __future__ import annotations

import itertools
import os
from pathlib import Path
import re
import subprocess
import tempfile

import numpy as np


RADII = np.array([10.0, 12.25, 9.0, 11.5, 8.25, 13.0, 7.5, 10.75])
DELTA = np.array([0.50, -0.20, 0.30, 0.70, -0.40, 0.10, 0.60, -0.30])
WEIGHT = np.array([1.20, 0.80, 1.10, 0.90, 1.30, 0.70, 1.40, 1.00])
FOREST = np.array([0, 0, 1, 1, 2, 2, 3, 3])

RP_MAX = 4.75
RP_BINS = 5
R3_MAX = 5.75
R3_BINS_PER_SIDE = 4


def catalog(directions: np.ndarray) -> np.ndarray:
    positions = RADII[:, None] * directions
    return np.column_stack((positions, DELTA, WEIGHT, FOREST))


COLLINEAR = catalog(np.tile(np.array([1.0, 0.0, 0.0]), (RADII.size, 1)))
WIDE_ANGLE = catalog(
    np.array(
        [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [-1.0, 0.0, 0.0],
            [0.0, -1.0, 0.0],
            [0.0, 0.0, -1.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
        ]
    )
)


def positive_bin(value: float, maximum: float, bins: int) -> int | None:
    if not 0.0 <= value < maximum:
        return None
    result = int(value / maximum * bins)
    return result if result < bins else None


def signed_bin(value: float, maximum: float, bins_per_side: int) -> int | None:
    if not -maximum < value < maximum:
        return None
    result = int((value + maximum) / maximum * bins_per_side)
    return result if 0 <= result < 2 * bins_per_side else None


def oracle_2pcf() -> dict[int, tuple[float, float]]:
    result: dict[int, list[float]] = {}
    for i, j in itertools.combinations(range(RADII.size), 2):
        if FOREST[i] == FOREST[j]:
            continue
        bin_index = positive_bin(abs(RADII[j] - RADII[i]), RP_MAX, RP_BINS)
        if bin_index is None:
            continue
        numerator = DELTA[i] * WEIGHT[i] * DELTA[j] * WEIGHT[j]
        denominator = WEIGHT[i] * WEIGHT[j]
        cell = result.setdefault(bin_index, [0.0, 0.0])
        cell[0] += numerator
        cell[1] += denominator
    return {key: tuple(value) for key, value in result.items()}


def oracle_3pcf() -> dict[tuple[int, int], tuple[float, float]]:
    result: dict[tuple[int, int], list[float]] = {}
    for pivot in range(RADII.size):
        neighbors = [
            index
            for index in range(RADII.size)
            if index != pivot
            and FOREST[index] != FOREST[pivot]
            and abs(RADII[index] - RADII[pivot]) < R3_MAX
        ]
        for q, r in itertools.combinations(neighbors, 2):
            if FOREST[q] == FOREST[r]:
                continue
            bq = signed_bin(RADII[q] - RADII[pivot], R3_MAX, R3_BINS_PER_SIDE)
            br = signed_bin(RADII[r] - RADII[pivot], R3_MAX, R3_BINS_PER_SIDE)
            assert bq is not None and br is not None
            numerator = np.prod(
                [DELTA[index] * WEIGHT[index] for index in (pivot, q, r)]
            )
            denominator = np.prod([WEIGHT[index] for index in (pivot, q, r)])
            for key in ((bq, br), (br, bq)):
                cell = result.setdefault(key, [0.0, 0.0])
                cell[0] += numerator
                cell[1] += denominator
    return {key: tuple(value) for key, value in result.items()}


def read_2pcf(path: Path) -> dict[int, tuple[float, float]]:
    result = {}
    for line in path.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split()
        if float(fields[4]) != 0.0:
            result[int(fields[0])] = (float(fields[3]), float(fields[4]))
    return result


def read_3pcf(path: Path) -> dict[tuple[int, int], tuple[float, float]]:
    result = {}
    for line in path.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split()
        result[(int(fields[0]), int(fields[1]))] = (
            float(fields[5]),
            float(fields[6]),
        )
    return result


def assert_histogram_close(actual: dict, expected: dict, label: str) -> None:
    if actual.keys() != expected.keys():
        raise AssertionError(f"{label}: occupied bins differ: {actual.keys() ^ expected.keys()}")
    for key in actual:
        if not np.allclose(actual[key], expected[key], rtol=2e-13, atol=2e-13):
            raise AssertionError(f"{label} bin {key}: {actual[key]} != {expected[key]}")


def run(
    binary: Path,
    input_path: Path,
    output: Path,
    threads: int,
    method: str = "lya-1d-2pcf-3pcf-omp",
) -> None:
    output.mkdir()
    command = [
        os.fspath(binary),
        f"search={method}",
        f"infile={input_path}",
        "infileformat=lya-ascii",
        "iCatalogs=1",
        f"rootDir={output}",
        f"numberThreads={threads}",
        "usePeriodic=false",
        "useLogHist=false",
        "rangeN=30",
        "rminHist=0.1",
        "sizeHistN=4",
        f"lya2RpMax={RP_MAX}",
        f"lya2RpBins={RP_BINS}",
        f"lya3RMax={R3_MAX}",
        f"lya3RBins={R3_BINS_PER_SIDE}",
        "verbose=0",
        "verbose_log=0",
        "options=",
    ]
    completed = subprocess.run(command, text=True, capture_output=True, check=False)
    if completed.returncode:
        raise RuntimeError(
            f"cballs failed with {threads} thread(s)\n"
            f"{completed.stdout}\n{completed.stderr}"
        )


def main() -> None:
    repository = Path(__file__).resolve().parents[2]
    binary = Path(os.environ.get("CBALLS", repository / "cballs")).resolve()
    if not binary.is_file():
        raise SystemExit(f"cballs executable not found: {binary}")
    build_env = subprocess.check_output(
        ["make", "--no-print-directory", "-s", "print-cyballs-build-env"],
        cwd=repository,
        text=True,
    )
    if "-DLYAFORESTOMP" not in build_env.split("__CBALLS_CPPFLAGS__=", 1)[-1]:
        raise SystemExit("the active build does not enable LYAFORESTOMP")
    block_match = re.search(
        r"-DLYA1D_OMP_PIVOT_BLOCK_SIZE=(\d+)", build_env
    )
    if block_match is None or int(block_match.group(1)) < 1:
        raise SystemExit("the active build has no valid radial pivot block size")
    block_size = int(block_match.group(1))

    with tempfile.TemporaryDirectory(prefix="ctreeballs-lya1d-") as temporary:
        root = Path(temporary)
        collinear_input = root / "collinear.txt"
        wide_angle_input = root / "wide_angle.txt"
        block_input = root / "multiple_blocks.txt"
        one_forest_input = root / "one_forest.txt"
        np.savetxt(collinear_input, COLLINEAR, fmt=("%.17g",) * 5 + ("%.0f",))
        np.savetxt(wide_angle_input, WIDE_ANGLE, fmt=("%.17g",) * 5 + ("%.0f",))
        block_count = 3 * block_size + 17
        indices = np.arange(block_count)
        block_radii = 100.0 + np.mod(indices * 0.3125, 60.0)
        block_catalog = np.column_stack(
            (
                block_radii,
                np.zeros(block_count),
                np.zeros(block_count),
                0.2 * np.sin(0.1 * indices),
                1.0 + 0.05 * np.mod(indices, 7),
                indices // 2,
            )
        )
        np.savetxt(block_input, block_catalog, fmt=("%.17g",) * 5 + ("%.0f",))
        one_forest_count = 129
        one_forest_indices = np.arange(one_forest_count)
        one_forest_catalog = np.column_stack(
            (
                200.0 + 0.03125 * one_forest_indices,
                np.zeros(one_forest_count),
                np.zeros(one_forest_count),
                np.cos(0.2 * one_forest_indices),
                0.5 + 0.01 * one_forest_indices,
                np.zeros(one_forest_count),
            )
        )
        np.savetxt(
            one_forest_input,
            one_forest_catalog,
            fmt=("%.17g",) * 5 + ("%.0f",),
        )

        one = root / "one"
        many = root / "many"
        transverse = root / "transverse"
        two_only = root / "two_only"
        three_only = root / "three_only"
        block_one = root / "block_one"
        block_many = root / "block_many"
        tree_one = root / "tree_one"
        tree_many = root / "tree_many"
        tree_block_one = root / "tree_block_one"
        tree_block_many = root / "tree_block_many"
        tree_one_forest = root / "tree_one_forest"
        run(binary, collinear_input, one, 1)
        run(binary, collinear_input, many, min(4, os.cpu_count() or 1))
        run(binary, wide_angle_input, transverse, 1)
        run(binary, collinear_input, two_only, 1, "lya-1d-2pcf-omp")
        run(binary, collinear_input, three_only, 1, "lya-1d-3pcf-omp")
        run(binary, block_input, block_one, 1, "lya-1d-2pcf-omp")
        run(
            binary,
            block_input,
            block_many,
            min(4, os.cpu_count() or 1),
            "lya-1d-2pcf-omp",
        )
        run(binary, collinear_input, tree_one, 1, "lya-1d-tree-2pcf-omp")
        run(
            binary,
            collinear_input,
            tree_many,
            min(4, os.cpu_count() or 1),
            "lya-1d-tree-2pcf-omp",
        )
        run(binary, block_input, tree_block_one, 1, "lya-1d-tree-2pcf-omp")
        run(
            binary,
            block_input,
            tree_block_many,
            min(4, os.cpu_count() or 1),
            "lya-1d-tree-2pcf-omp",
        )
        run(
            binary,
            one_forest_input,
            tree_one_forest,
            min(4, os.cpu_count() or 1),
            "lya-1d-tree-2pcf-omp",
        )

        output2 = "histXi2pcf_lya1d.txt"
        output3 = "histZetaM_lya1d.txt"
        assert_histogram_close(read_2pcf(one / output2), oracle_2pcf(), "1D 2PCF")
        assert_histogram_close(read_3pcf(one / output3), oracle_3pcf(), "1D 3PCF")
        assert_histogram_close(
            read_2pcf(tree_one / output2), oracle_2pcf(), "1D tree 2PCF"
        )
        for name in (output2, output3):
            baseline = (one / name).read_bytes()
            if baseline != (many / name).read_bytes():
                raise AssertionError(f"OpenMP output is not deterministic: {name}")
            if baseline != (transverse / name).read_bytes():
                raise AssertionError(f"transverse coordinates changed radial output: {name}")
        if (two_only / output2).read_bytes() != (one / output2).read_bytes():
            raise AssertionError("2PCF-only output differs from the combined method")
        if (three_only / output3).read_bytes() != (one / output3).read_bytes():
            raise AssertionError("3PCF-only output differs from the combined method")
        if (two_only / output3).exists() or (three_only / output2).exists():
            raise AssertionError("a radial-only method emitted a disabled product")
        if (block_one / output2).read_bytes() != (block_many / output2).read_bytes():
            raise AssertionError("multi-block radial 2PCF is not deterministic")
        if (tree_one / output2).read_bytes() != (tree_many / output2).read_bytes():
            raise AssertionError("radial tree 2PCF is not deterministic")
        if (tree_block_one / output2).read_bytes() != (
            tree_block_many / output2
        ).read_bytes():
            raise AssertionError("multi-task radial tree 2PCF is not deterministic")
        assert_histogram_close(
            read_2pcf(tree_block_one / output2),
            read_2pcf(block_one / output2),
            "tree versus range-scan 2PCF",
        )
        if read_2pcf(tree_one_forest / output2):
            raise AssertionError("same-forest subtraction left nonzero 2PCF bins")
        if (tree_one / output3).exists():
            raise AssertionError("radial tree 2PCF emitted a disabled 3PCF product")

    print("PASS: radial Ly-alpha scan/tree oracles, forest exclusion, and deterministic OpenMP reductions")


if __name__ == "__main__":
    main()
