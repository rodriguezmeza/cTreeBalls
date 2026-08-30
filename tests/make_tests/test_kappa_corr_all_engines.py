#!/usr/bin/env python3

import json
from pathlib import Path
import sys
import tempfile

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))

from cyballs import search_method_id
from kappa_corr_all_engines import KappaCatalog, RunConfig, run_engine_suite


def unit_sphere_catalog(count=32):
    indices = np.arange(count, dtype=np.float64)
    z = 1.0 - 2.0 * (indices + 0.5) / count
    phi = np.pi * (3.0 - np.sqrt(5.0)) * indices
    radius = np.sqrt(1.0 - z * z)
    positions = np.column_stack((radius * np.cos(phi), radius * np.sin(phi), z))
    kappa = positions[:, 0] - 0.4 * positions[:, 1] + 0.2 * positions[:, 2]
    kappa -= np.mean(kappa)
    return KappaCatalog(positions=positions, kappa=kappa)


def test_one_catalog_is_reused_across_enabled_engines():
    candidates = ("octree-ggg-omp", "kdtree-omp")
    engines = tuple(name for name in candidates if search_method_id(name) >= 0)
    if not engines:
        raise AssertionError("the build has no engine suitable for the kappa smoke test")

    with tempfile.TemporaryDirectory(prefix="ctreeballs-kappa-suite-") as tmp:
        output = Path(tmp)
        config = RunConfig(
            engines=engines,
            output_dir=output,
            theta_min=5.0,
            theta_max=100.0,
            bins=4,
            multipoles=2,
            threads=2,
            options=("only-2pcf", "no-out-Hist"),
            verbose=0,
            verbose_log=0,
        )
        results = run_engine_suite(unit_sphere_catalog(), config)

        if tuple(results) != engines:
            raise AssertionError(f"expected {engines}, got {tuple(results)}")
        for engine, result in results.items():
            if result["nbody"] != 32:
                raise AssertionError(f"{engine} published the wrong body count")
            if not {"r", "xi", "nn"}.issubset(result):
                raise AssertionError(f"{engine} did not publish its 2PCF arrays")

        summary = json.loads((output / "summary.json").read_text(encoding="utf-8"))
        if summary["catalog_reads"] != 1:
            raise AssertionError("the suite did not preserve the one-read contract")
        if summary["set_catalog_calls_per_process"] != 1:
            raise AssertionError("the suite registered the catalog more than once")
        if summary["failures"]:
            raise AssertionError(f"suite failures: {summary['failures']}")


if __name__ == "__main__":
    test_one_catalog_is_reused_across_enabled_engines()
    print("PASS: kappa engines reuse one in-memory cyballs catalog")
