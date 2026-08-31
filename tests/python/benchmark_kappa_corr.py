#!/usr/bin/env python3
"""Compatibility entry point for the maintained multi-backend benchmark.

Original convergence benchmark by Axel Romero Tisnado.
The maintained runner preserves the common --kappa/--nsides/--threads options.
Output now uses timings.csv, summary.csv, comparisons.csv and values.csv.
"""

import os
from pathlib import Path
import runpy
import sys


def main():
    root = Path(__file__).resolve().parents[2]
    script = root / "addons/python_env/cputime_comparison/benchmark_kappa_corr.py"
    if not script.is_file():
        raise SystemExit(
            "The maintained benchmark is not bundled in this source installation. "
            "Run benchmark_kappa_corr.py from your cputime_comparison environment."
        )
    os.environ.setdefault("CTREEBALLS_ROOT", str(root))
    sys.path.insert(0, str(script.parent))
    if "--kappa" in sys.argv and "--scenarios" not in sys.argv:
        sys.argv.extend(["--scenarios", "sphere-convergence"])
    runpy.run_path(str(script), run_name="__main__")


if __name__ == "__main__":
    main()
