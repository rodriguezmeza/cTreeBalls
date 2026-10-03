#!/usr/bin/env python3
"""Compare the production Ly-alpha 2PCF with the original lya2pcf CPU code."""

from __future__ import annotations

import importlib
import io
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np


RP_MAX = 64.0
RT_MAX = 64.0
RP_BINS = 8
RT_BINS = 8

# Each item is ra, dec, comoving distances, weights, and unweighted deltas.
FORESTS = (
    (0.03, -0.02, (80.0, 94.0, 111.0), (1.20, 0.80, 1.10), (0.50, -0.20, 0.30)),
    (0.11, 0.04, (82.0, 101.0, 119.0), (0.90, 1.30, 0.70), (-0.40, 0.60, 0.10)),
    (0.22, -0.01, (77.0, 96.0, 124.0), (1.40, 1.00, 0.60), (0.20, -0.50, 0.70)),
    (0.34, 0.07, (85.0, 105.0, 128.0), (0.75, 1.25, 0.95), (-0.30, 0.40, -0.10)),
    # This distant sightline exercises angular/r_transverse rejection.
    (1.40, -0.03, (88.0, 108.0), (1.05, 0.85), (0.35, -0.25)),
)


def reference_source() -> Path:
    value = os.environ.get("LYA2PCF_SOURCE", "")
    if not value:
        raise SystemExit(
            "set LYA2PCF_SOURCE to the original lya2pcf checkout "
            "(the directory containing correlation_procedures_cpu.py)"
        )
    source = Path(value).expanduser().resolve()
    required = (
        source / "parameters.py",
        source / "forest_class.py",
        source / "correlation_procedures_cpu.py",
    )
    missing = [path.name for path in required if not path.is_file()]
    if missing:
        raise SystemExit(f"invalid LYA2PCF_SOURCE={source}: missing {', '.join(missing)}")
    return source


def load_reference(source: Path):
    sys.dont_write_bytecode = True
    sys.path.insert(0, os.fspath(source))
    try:
        forest_module = importlib.import_module("forest_class")
        correlation_module = importlib.import_module("correlation_procedures_cpu")
    except ModuleNotFoundError as exc:
        raise SystemExit(
            f"the lya2pcf reference dependencies are incomplete: {exc.name!r} is missing"
        ) from exc
    return forest_module, correlation_module


def make_reference_data(forest_module) -> dict[int, list[object]]:
    data: dict[int, list[object]] = {}
    for index, (ra, dec, distances, weights, deltas) in enumerate(FORESTS):
        distance = np.asarray(distances, dtype=np.float64)
        weight = np.asarray(weights, dtype=np.float64)
        delta = np.asarray(deltas, dtype=np.float64)
        forest = forest_module.quasar(
            fib=index, pl=index, name=1000 + index, ra=ra, dec=dec,
            lenght=distance.size,
        )
        forest.dc = distance
        forest.we = weight
        forest.dw = weight * delta
        forest.rx = distance * forest.x
        forest.ry = distance * forest.y
        forest.rz = distance * forest.z
        data.setdefault(int(forest.pix), []).append(forest)
    return data


def run_reference(correlations, data: dict[int, list[object]]) -> tuple[np.ndarray, np.ndarray]:
    shape = (RP_BINS, RT_BINS)
    minimum_distance = min(float(forest.dc[0]) for forests in data.values() for forest in forests)
    angular_limit = 2.0 * math.asin(0.5 * RT_MAX / minimum_distance)

    # pair_correlation is JIT-compiled on its first call, after these production
    # module globals have been set to the same histogram contract as cTreeBalls.
    correlations.rpmax = RP_MAX
    correlations.rtmax = RT_MAX
    correlations.numpix_rp = RP_BINS
    correlations.numpix_rt = RT_BINS
    correlations.init(data, io.StringIO(), shape, angular_limit)

    denominator = np.zeros(shape, dtype=np.float64)
    numerator = np.zeros(shape, dtype=np.float64)
    for pixel in sorted(data):
        weight_histogram, weighted_delta_histogram = correlations.two_point_per_pixel(pixel)
        denominator += weight_histogram
        numerator += weighted_delta_histogram
    return numerator, denominator


def convert_catalog(source: Path, data: dict[int, list[object]], root: Path) -> Path:
    saved = root / "data.npy"
    catalog = root / "forest.txt"
    np.save(saved, data, allow_pickle=True)
    converter = (
        Path(__file__).resolve().parents[2]
        / "addons/lya_forest_omp/lya2pcf_to_cballs.py"
    )
    completed = subprocess.run(
        [
            sys.executable,
            os.fspath(converter),
            "--lya2pcf-source",
            os.fspath(source),
            "-o",
            os.fspath(catalog),
            os.fspath(saved),
        ],
        text=True,
        capture_output=True,
        check=False,
        env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"},
    )
    if completed.returncode:
        raise RuntimeError(f"lya2pcf converter failed\n{completed.stdout}\n{completed.stderr}")
    return catalog


def run_cballs(binary: Path, catalog: Path, output: Path) -> None:
    output.mkdir()
    command = [
        os.fspath(binary),
        "search=lya-2pcf-omp",
        f"infile={catalog}",
        "infileformat=lya-ascii",
        "iCatalogs=1",
        f"rootDir={output}",
        "numberThreads=2",
        "usePeriodic=false",
        "useLogHist=false",
        "rangeN=64",
        "rminHist=0.1",
        "sizeHistN=4",
        f"lya2RpMax={RP_MAX}",
        f"lya2RtMax={RT_MAX}",
        f"lya2RpBins={RP_BINS}",
        f"lya2RtBins={RT_BINS}",
        "verbose=0",
        "verbose_log=0",
        "options=",
    ]
    completed = subprocess.run(command, text=True, capture_output=True, check=False)
    if completed.returncode:
        raise RuntimeError(f"cballs failed\n{completed.stdout}\n{completed.stderr}")


def read_cballs(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    rows = np.loadtxt(path, comments="#", ndmin=2)
    if rows.shape != (RP_BINS * RT_BINS, 7):
        raise AssertionError(f"unexpected cTreeBalls 2PCF shape: {rows.shape}")
    numerator = np.zeros((RP_BINS, RT_BINS), dtype=np.float64)
    denominator = np.zeros_like(numerator)
    correlation = np.zeros_like(numerator)
    seen: set[tuple[int, int]] = set()
    for row in rows:
        key = int(row[0]), int(row[1])
        if key in seen or not (0 <= key[0] < RP_BINS and 0 <= key[1] < RT_BINS):
            raise AssertionError(f"invalid or duplicate cTreeBalls bin: {key}")
        seen.add(key)
        correlation[key] = row[4]
        numerator[key] = row[5]
        denominator[key] = row[6]
    return numerator, denominator, correlation


def assert_same(actual: np.ndarray, expected: np.ndarray, label: str) -> None:
    np.testing.assert_allclose(actual, expected, rtol=3e-13, atol=3e-13, err_msg=label)


def main() -> None:
    source = reference_source()
    repository = Path(__file__).resolve().parents[2]
    binary = Path(os.environ.get("CBALLS", repository / "cballs")).resolve()
    if not binary.is_file():
        raise SystemExit(f"cballs executable not found: {binary}")

    forest_module, correlations = load_reference(source)
    data = make_reference_data(forest_module)
    reference_numerator, reference_denominator = run_reference(correlations, data)
    if np.count_nonzero(reference_denominator) < 4:
        raise AssertionError("reference fixture does not exercise enough occupied bins")
    reference_correlation = np.divide(
        reference_numerator,
        reference_denominator,
        out=np.zeros_like(reference_numerator),
        where=reference_denominator != 0.0,
    )

    with tempfile.TemporaryDirectory(prefix="ctreeballs-lya2pcf-reference-") as temporary:
        root = Path(temporary)
        catalog = convert_catalog(source, data, root)
        output = root / "cballs"
        run_cballs(binary, catalog, output)
        actual_numerator, actual_denominator, actual_correlation = read_cballs(
            output / "histXi2pcf_lya.txt"
        )

    assert_same(actual_numerator, reference_numerator, "weighted-delta numerator")
    assert_same(actual_denominator, reference_denominator, "weight denominator")
    assert_same(actual_correlation, reference_correlation, "normalized correlation")
    print(f"PASS: cTreeBalls 2PCF matches lya2pcf CPU reference at {source}")


if __name__ == "__main__":
    main()
