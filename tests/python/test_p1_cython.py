#!/usr/bin/env python3

import tempfile

from cyballs import CosmoComputationError, cballs


def small_cubic_parameters(root_dir):
    return {
        "testmodel": "simple-cubic",
        "nbody": 8,
        "rootDir": root_dir,
        "numberThreads": 2,
        "verbose": 0,
        "verbose_log": 0,
        "options": "no-out-Hist",
    }


def test_failed_startup_releases_body_tables():
    with tempfile.TemporaryDirectory(prefix="ctreeballs-failed-startup-") as root_dir:
        for _ in range(8):
            balls = cballs()
            parameters = small_cubic_parameters(root_dir)
            parameters["searchMethod"] = "not-a-compiled-search-method"
            balls.set(parameters)

            try:
                balls.Run(level=["StartRun_Common"])
            except CosmoComputationError:
                pass
            else:
                raise AssertionError("invalid search method unexpectedly succeeded")

            if balls.getBodytableAllocated():
                raise AssertionError("failed startup retained body-table ownership")


def test_compiled_default_search_runs():
    with tempfile.TemporaryDirectory(prefix="ctreeballs-default-search-") as root_dir:
        balls = cballs()
        balls.set(small_cubic_parameters(root_dir))
        try:
            balls.Run(level=["MainLoop"])
            if not balls.getBodytableAllocated():
                raise AssertionError("successful run did not publish body-table ownership")
        finally:
            balls.struct_cleanup()

        if balls.getBodytableAllocated():
            raise AssertionError("explicit cleanup retained body-table ownership")


if __name__ == "__main__":
    test_failed_startup_releases_body_tables()
    test_compiled_default_search_runs()
    print("PASS: Cython startup ownership and compiled default search")
