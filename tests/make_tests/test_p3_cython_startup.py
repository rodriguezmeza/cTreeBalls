#!/usr/bin/env python3

import gc
import os
import resource
import subprocess
import sys
import tempfile

from cyballs import CosmoComputationError, cballs


def current_rss_bytes():
    try:
        with open("/proc/self/statm", encoding="ascii") as statm:
            resident_pages = int(statm.read().split()[1])
        return resident_pages * os.sysconf("SC_PAGE_SIZE")
    except (FileNotFoundError, IndexError, OSError, ValueError):
        try:
            rss_kib = subprocess.check_output(
                ["ps", "-o", "rss=", "-p", str(os.getpid())],
                encoding="ascii",
            )
            return int(rss_kib.strip()) * 1024
        except (OSError, subprocess.SubprocessError, ValueError):
            maximum_rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
            if sys.platform == "darwin":
                return maximum_rss
            return maximum_rss * 1024


def base_parameters(root_dir):
    return {
        "searchMethod": "octree-ggg-omp",
        "testmodel": "simple-cubic",
        "nbody": 8,
        "rootDir": root_dir,
        "numberThreads": 2,
        "verbose": 0,
        "verbose_log": 0,
        "options": "no-out-Hist",
    }


def invalid_parameters(index, root_dir):
    parameters = base_parameters(root_dir)
    scenario = index % 3

    if scenario == 0:
        parameters["searchMethod"] = "not-a-compiled-search-method"
    elif scenario == 1:
        parameters.update(
            {
                "useLogHist": False,
                "rangeN": 0.1,
                "rminHist": 0.1,
            }
        )
    else:
        parameters.update(
            {
                "infile": os.path.join(root_dir, "missing-catalog.dat"),
                "infileformat": "ascii",
            }
        )

    return parameters


def assert_no_c_ownership(balls, iteration):
    flags = {
        "command": balls.getCMDAllocated(),
        "global": balls.getGDAllocated(),
        "global-2": balls.getAllocated2(),
        "histograms": balls.getHistogramsAllocated(),
        "tree": balls.getTreeAllocated(),
        "body tables": balls.getBodytableAllocated(),
    }
    retained = [name for name, value in flags.items() if value]
    if retained:
        raise AssertionError(
            f"failure {iteration} retained C ownership: {', '.join(retained)}"
        )


def test_repeated_failures_have_bounded_rss():
    repeats = int(os.environ.get("CBALLS_FAILURE_REPEATS", "120"))
    warmup = int(os.environ.get("CBALLS_FAILURE_WARMUP", "24"))
    limit_mb = float(os.environ.get("CBALLS_FAILURE_RSS_LIMIT_MB", "32"))
    if repeats <= warmup:
        raise ValueError("CBALLS_FAILURE_REPEATS must exceed the warmup count")

    with tempfile.TemporaryDirectory(prefix="ctreeballs-p3-failures-") as root_dir:
        baseline_rss = None
        peak_rss = 0

        for index in range(repeats):
            balls = cballs()
            balls.set(invalid_parameters(index, root_dir))
            try:
                balls.Run(level=["StartRun_Common"])
            except CosmoComputationError:
                pass
            else:
                raise AssertionError(f"failure scenario {index % 3} succeeded")

            assert_no_c_ownership(balls, index)
            balls.struct_cleanup()
            del balls

            if index == warmup - 1 or index >= warmup and index % 8 == 0:
                gc.collect()
                rss = current_rss_bytes()
                if baseline_rss is None:
                    baseline_rss = rss
                peak_rss = max(peak_rss, rss)

        gc.collect()
        peak_rss = max(peak_rss, current_rss_bytes())
        growth = peak_rss - baseline_rss
        limit = int(limit_mb * 1024 * 1024)
        if growth > limit:
            raise AssertionError(
                "repeated failed startups grew RSS by "
                f"{growth / (1024 * 1024):.1f} MiB "
                f"after warmup (limit {limit_mb:.1f} MiB)"
            )

        print(
            "PASS: repeated Cython startup failures retained no C ownership; "
            f"post-warmup RSS growth={growth / (1024 * 1024):.1f} MiB"
        )


def test_repeated_partial_tree_failures_have_bounded_rss():
    repeats = int(os.environ.get("CBALLS_TREE_FAILURE_REPEATS", "40"))
    warmup = int(os.environ.get("CBALLS_TREE_FAILURE_WARMUP", "8"))
    limit_mb = float(os.environ.get("CBALLS_TREE_FAILURE_RSS_LIMIT_MB", "16"))
    if repeats <= warmup:
        raise ValueError("CBALLS_TREE_FAILURE_REPEATS must exceed the warmup count")

    with tempfile.TemporaryDirectory(prefix="ctreeballs-p3-tree-failures-") as root_dir:
        baseline_rss = None
        peak_rss = 0

        for index in range(repeats):
            balls = cballs()
            parameters = base_parameters(root_dir)
            parameters.update(
                {
                    "testmodel": "simple-cubic",
                    "nbody": 32768,
                    "numberThreads": 1,
                    "theta": 0.5,
                }
            )
            balls.set(parameters)

            try:
                balls.Run(level=["MainLoop"])
            except CosmoComputationError as error:
                if "cellRadius index out of range" not in str(error):
                    raise AssertionError(
                        f"partial-tree failure {index} failed for the wrong reason: {error}"
                    ) from error
            else:
                raise AssertionError(
                    f"partial-tree failure {index} unexpectedly succeeded"
                )

            assert_no_c_ownership(balls, index)
            balls.struct_cleanup()
            del balls

            if index == warmup - 1 or index >= warmup and index % 4 == 0:
                gc.collect()
                rss = current_rss_bytes()
                if baseline_rss is None:
                    baseline_rss = rss
                peak_rss = max(peak_rss, rss)

        gc.collect()
        peak_rss = max(peak_rss, current_rss_bytes())
        growth = peak_rss - baseline_rss
        limit = int(limit_mb * 1024 * 1024)
        if growth > limit:
            raise AssertionError(
                "repeated partial-tree failures grew RSS by "
                f"{growth / (1024 * 1024):.1f} MiB "
                f"after warmup (limit {limit_mb:.1f} MiB)"
            )

        print(
            "PASS: repeated partial-tree failures retained no C ownership; "
            f"post-warmup RSS growth={growth / (1024 * 1024):.1f} MiB"
        )


if __name__ == "__main__":
    test_repeated_failures_have_bounded_rss()
    test_repeated_partial_tree_failures_have_bounded_rss()
