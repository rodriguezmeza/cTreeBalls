#!/usr/bin/env python3

import tempfile

import numpy as np

from cyballs import cballs


def parameters(root_dir, nbody):
    return {
        "searchMethod": "octree-ggg-omp",
        "testmodel": "simple-cubic",
        "nbody": nbody,
        "rootDir": root_dir,
        "numberThreads": 2,
        "verbose": 0,
        "verbose_log": 0,
        "options": "no-out-Hist",
    }


def test_two_live_instances_keep_independent_results_and_ownership():
    with tempfile.TemporaryDirectory(prefix="ctreeballs-p0-instances-") as root_dir:
        first = cballs()
        second = cballs()
        first.set(parameters(root_dir, 8))
        second.set(parameters(root_dir, 27))

        first.Run(level=["MainLoop"])
        first_histogram = first.getHistNN().copy()
        second.Run(level=["MainLoop"])
        second_histogram = second.getHistNN().copy()

        first_bodytable = first._runtime_bodytable_address()
        second_bodytable = second._runtime_bodytable_address()
        if not first_bodytable or not second_bodytable:
            raise AssertionError("a live object has no runtime body table")
        if first_bodytable == second_bodytable:
            raise AssertionError("live objects share one runtime body table")

        np.testing.assert_array_equal(first.getHistNN(), first_histogram)
        first.struct_cleanup()

        if first._runtime_bodytable_address():
            raise AssertionError("cleaning the first object retained its body table")
        if second._runtime_bodytable_address() != second_bodytable:
            raise AssertionError("cleaning the first object changed the second body table")

        if not second.getBodytableAllocated() or not second.getTreeAllocated():
            raise AssertionError("cleaning the first object changed second ownership")
        np.testing.assert_array_equal(second.getHistNN(), second_histogram)
        second.struct_cleanup()
        if second._runtime_bodytable_address():
            raise AssertionError("cleaning the second object retained its body table")


if __name__ == "__main__":
    test_two_live_instances_keep_independent_results_and_ownership()
    print("PASS: multiple live cballs objects retain independent C state")
