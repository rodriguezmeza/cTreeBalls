#!/usr/bin/env python3

import ctypes
import gc
import os
import resource
import subprocess
import sys
import tempfile

from cyballs import CosmoComputationError, cballs


def current_memory_bytes():
    # Darwin can retain hundreds of MiB of freed malloc pages in RSS. Measure
    # live allocations there so allocator caching is not mistaken for a leak.
    # These constructors allocate through malloc/calloc, covered by all zones.
    if sys.platform == "darwin":
        class MallocStatistics(ctypes.Structure):
            _fields_ = [("blocks_in_use", ctypes.c_uint),
                        ("size_in_use", ctypes.c_size_t),
                        ("max_size_in_use", ctypes.c_size_t),
                        ("size_allocated", ctypes.c_size_t)]
        libc = ctypes.CDLL("/usr/lib/libSystem.B.dylib")
        libc.malloc_zone_statistics.argtypes = [
            ctypes.c_void_p, ctypes.POINTER(MallocStatistics)]
        libc.malloc_zone_statistics.restype = None
        statistics = MallocStatistics()
        libc.malloc_zone_statistics(None, ctypes.byref(statistics))
        return statistics.size_in_use

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
        "searchMethod": "octree-2balls-omp",
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


def test_repeated_failures_have_bounded_memory():
    repeats = int(os.environ.get("CBALLS_FAILURE_REPEATS", "120"))
    warmup = int(os.environ.get("CBALLS_FAILURE_WARMUP", "24"))
    limit_mb = float(os.environ.get("CBALLS_FAILURE_RSS_LIMIT_MB", "32"))
    if repeats <= warmup:
        raise ValueError("CBALLS_FAILURE_REPEATS must exceed the warmup count")

    with tempfile.TemporaryDirectory(prefix="ctreeballs-p3-failures-") as root_dir:
        baseline_memory = None
        peak_memory = 0

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
                rss = current_memory_bytes()
                if baseline_memory is None:
                    baseline_memory = rss
                peak_memory = max(peak_memory, rss)

        gc.collect()
        peak_memory = max(peak_memory, current_memory_bytes())
        growth = peak_memory - baseline_memory
        limit = int(limit_mb * 1024 * 1024)
        if growth > limit:
            raise AssertionError(
                "repeated failed startups grew retained memory by "
                f"{growth / (1024 * 1024):.1f} MiB "
                f"after warmup (limit {limit_mb:.1f} MiB)"
            )

        print(
            "PASS: repeated Cython startup failures retained no C ownership; "
            f"post-warmup memory growth={growth / (1024 * 1024):.1f} MiB"
        )


def test_repeated_large_tree_runs_have_bounded_memory():
    repeats = int(os.environ.get("CBALLS_TREE_FAILURE_REPEATS", "40"))
    warmup = int(os.environ.get("CBALLS_TREE_FAILURE_WARMUP", "8"))
    limit_mb = float(os.environ.get("CBALLS_TREE_FAILURE_RSS_LIMIT_MB", "16"))
    if repeats <= warmup:
        raise ValueError("CBALLS_TREE_FAILURE_REPEATS must exceed the warmup count")

    with tempfile.TemporaryDirectory(prefix="ctreeballs-p3-tree-failures-") as root_dir:
        baseline_memory = None
        peak_memory = 0

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

            # setradius now clamps large radii into the overflow bin. This
            # formerly failing fixture must complete and release its full tree.
            balls.Run(level=["MainLoop"])
            if balls.getNBody() <= 0:
                raise AssertionError("large-tree run lost its input catalog")
            balls.struct_cleanup()
            assert_no_c_ownership(balls, index)
            del balls

            if index == warmup - 1 or index >= warmup and index % 4 == 0:
                gc.collect()
                rss = current_memory_bytes()
                if baseline_memory is None:
                    baseline_memory = rss
                peak_memory = max(peak_memory, rss)

        gc.collect()
        peak_memory = max(peak_memory, current_memory_bytes())
        growth = peak_memory - baseline_memory
        limit = int(limit_mb * 1024 * 1024)
        if growth > limit:
            raise AssertionError(
                "repeated large-tree runs grew retained memory by "
                f"{growth / (1024 * 1024):.1f} MiB "
                f"after warmup (limit {limit_mb:.1f} MiB)"
            )

        print(
            "PASS: repeated large-tree runs released all C ownership; "
            f"post-warmup memory growth={growth / (1024 * 1024):.1f} MiB"
        )


if __name__ == "__main__":
    test_repeated_failures_have_bounded_memory()
    test_repeated_large_tree_runs_have_bounded_memory()
