#!/usr/bin/env python3
"""Independent oracle and OpenMP determinism test for the Ly-alpha addon."""

from __future__ import annotations

import itertools
import math
import os
from pathlib import Path
import subprocess
import tempfile

import numpy as np


POINTS = np.array(
    [
        (10.0, 0.0, 0.0, 0.50, 1.20, 0),
        (12.0, 0.0, 0.0, -0.20, 0.80, 0),
        (9.0, 3.0, 0.0, 0.30, 1.10, 1),
        (11.0, 4.0, 0.0, 0.70, 0.90, 1),
        (8.0, 0.0, 4.0, -0.40, 1.30, 2),
        (13.0, 0.0, 3.0, 0.10, 0.70, 2),
        (7.0, -2.0, 3.0, 0.60, 1.40, 3),
        (10.0, -3.0, 2.0, -0.30, 1.00, 3),
    ],
    dtype=np.float64,
)

RP_MAX = 30.0
RT_MAX = 30.0
RP_BINS = 5
RT_BINS = 6
R3_MAX = 30.0
R3_BINS = 4
THETA_BINS = 5
MU_BINS = 6


def positive_bin(value: float, maximum: float, bins: int) -> int | None:
    if not 0.0 <= value < maximum:
        return None
    result = int(value / maximum * bins)
    return result if result < bins else None


def theta_bin(value: float) -> int:
    return min(int(np.clip(value, 0.0, math.pi) / math.pi * THETA_BINS), THETA_BINS - 1)


def mu_bin(value: float) -> int:
    return min(int((np.clip(value, -1.0, 1.0) + 1.0) * 0.5 * MU_BINS), MU_BINS - 1)


def oracle_2pcf() -> dict[tuple[int, int], tuple[float, float]]:
    result: dict[tuple[int, int], list[float]] = {}
    for left, right in itertools.combinations(POINTS, 2):
        if int(left[5]) == int(right[5]):
            continue
        p = left[:3]
        q = right[:3]
        dp = np.linalg.norm(p)
        dq = np.linalg.norm(q)
        cosine = np.clip(np.dot(p / dp, q / dq), -1.0, 1.0)
        rp = abs(dp - dq) * math.sqrt(max(0.0, 0.5 * (1.0 + cosine)))
        rt = (dp + dq) * math.sqrt(max(0.0, 0.5 * (1.0 - cosine)))
        bp = positive_bin(rp, RP_MAX, RP_BINS)
        bt = positive_bin(rt, RT_MAX, RT_BINS)
        if bp is None or bt is None:
            continue
        numerator = left[3] * left[4] * right[3] * right[4]
        denominator = left[4] * right[4]
        cell = result.setdefault((bp, bt), [0.0, 0.0])
        cell[0] += numerator
        cell[1] += denominator
    return {key: tuple(value) for key, value in result.items()}


def oracle_3pcf() -> dict[tuple[int, int, int, int, int], tuple[float, float]]:
    result: dict[tuple[int, int, int, int, int], list[float]] = {}
    for pivot_index, pivot in enumerate(POINTS):
        los = pivot[:3] / np.linalg.norm(pivot[:3])
        neighbors = [
            candidate
            for index, candidate in enumerate(POINTS)
            if index != pivot_index
            and int(candidate[5]) != int(pivot[5])
            and np.linalg.norm(candidate[:3] - pivot[:3]) < R3_MAX
        ]
        for q, r in itertools.combinations(neighbors, 2):
            if int(q[5]) == int(r[5]):
                continue
            dq = q[:3] - pivot[:3]
            dr = r[:3] - pivot[:3]
            rq = np.linalg.norm(dq)
            rr = np.linalg.norm(dr)
            bq = positive_bin(rq, R3_MAX, R3_BINS)
            br = positive_bin(rr, R3_MAX, R3_BINS)
            assert bq is not None and br is not None
            tq = theta_bin(math.acos(np.clip(np.dot(dq, los) / rq, -1.0, 1.0)))
            tr = theta_bin(math.acos(np.clip(np.dot(dr, los) / rr, -1.0, 1.0)))
            bm = mu_bin(np.dot(dq, dr) / (rq * rr))
            numerator = np.prod((pivot[3] * pivot[4], q[3] * q[4], r[3] * r[4]))
            denominator = pivot[4] * q[4] * r[4]
            for key in ((bq, br, tq, tr, bm), (br, bq, tr, tq, bm)):
                cell = result.setdefault(key, [0.0, 0.0])
                cell[0] += numerator
                cell[1] += denominator
    return {key: tuple(value) for key, value in result.items()}


def read_2pcf(path: Path) -> dict[tuple[int, int], tuple[float, float]]:
    result = {}
    for line in path.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split()
        if float(fields[6]) != 0.0:
            result[(int(fields[0]), int(fields[1]))] = (float(fields[5]), float(fields[6]))
    return result


def read_3pcf(path: Path) -> dict[tuple[int, int, int, int, int], tuple[float, float]]:
    result = {}
    for line in path.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split()
        result[tuple(map(int, fields[:5]))] = (float(fields[11]), float(fields[12]))
    return result


def assert_histogram_close(actual: dict, expected: dict, label: str) -> None:
    if actual.keys() != expected.keys():
        raise AssertionError(f"{label}: occupied bins differ: {actual.keys() ^ expected.keys()}")
    for key in actual:
        if not np.allclose(actual[key], expected[key], rtol=2e-13, atol=2e-13):
            raise AssertionError(f"{label} bin {key}: {actual[key]} != {expected[key]}")


def run(binary: Path, catalog: Path, output: Path, threads: int) -> None:
    output.mkdir()
    command = [
        os.fspath(binary),
        "search=lya-2pcf-3pcf-omp",
        f"infile={catalog}",
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
        f"lya2RtMax={RT_MAX}",
        f"lya2RpBins={RP_BINS}",
        f"lya2RtBins={RT_BINS}",
        f"lya3RMax={R3_MAX}",
        f"lya3RBins={R3_BINS}",
        f"lya3ThetaBins={THETA_BINS}",
        f"lya3MuBins={MU_BINS}",
        "verbose=0",
        "verbose_log=0",
        "options=",
    ]
    completed = subprocess.run(command, text=True, capture_output=True, check=False)
    if completed.returncode:
        raise RuntimeError(
            f"cballs failed with {threads} thread(s)\n{completed.stdout}\n{completed.stderr}"
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

    with tempfile.TemporaryDirectory(prefix="ctreeballs-lya-") as temporary:
        root = Path(temporary)
        catalog = root / "forest.txt"
        np.savetxt(catalog, POINTS, fmt=("%.17g",) * 5 + ("%.0f",))
        one = root / "one"
        many = root / "many"
        run(binary, catalog, one, 1)
        run(binary, catalog, many, min(4, os.cpu_count() or 1))

        assert_histogram_close(
            read_2pcf(one / "histXi2pcf_lya.txt"), oracle_2pcf(), "2PCF"
        )
        assert_histogram_close(
            read_3pcf(one / "histZetaM_lya5d.txt"), oracle_3pcf(), "3PCF"
        )
        for name in ("histXi2pcf_lya.txt", "histZetaM_lya5d.txt"):
            if (one / name).read_bytes() != (many / name).read_bytes():
                raise AssertionError(f"OpenMP output is not deterministic: {name}")

    print("PASS: Ly-alpha 2PCF/3PCF oracle and OpenMP determinism")


if __name__ == "__main__":
    main()
