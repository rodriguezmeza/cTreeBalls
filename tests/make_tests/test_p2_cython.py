#!/usr/bin/env python3

import math
import os
import tempfile

from cyballs import CosmoComputationError, CosmoSevereError, cballs


def assert_no_c_ownership(balls, context):
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
            f"{context} retained C ownership: {', '.join(retained)}"
        )


def test_allocation_failures_raise_and_recover():
    repeats = int(os.environ.get("CBALLS_ALLOCATION_FAILURE_REPEATS", "24"))
    sweep_limit = int(
        os.environ.get("CBALLS_ALLOCATION_FAILURE_SWEEP_LIMIT", "256")
    )

    with tempfile.TemporaryDirectory(prefix="ctreeballs-p2-allocation-") as root_dir:
        balls = cballs()
        parameters = {
            "testmodel": "simple-cubic",
            "nbody": 8,
            "rootDir": root_dir,
            "options": "no-out-Hist",
            "numberThreads": 2,
            "verbose": 0,
            "verbose_log": 0,
        }

        for iteration in range(repeats):
            balls.set_default(**parameters)
            balls._set_allocation_failure_after_for_tests(0)
            try:
                balls.Run(level=["StartRun_Common"])
            except CosmoComputationError as exc:
                if "memory allocation failed" not in str(exc):
                    raise AssertionError(
                        f"allocation failure {iteration} returned the wrong diagnostic: {exc}"
                    ) from exc
            else:
                raise AssertionError(
                    f"injected allocation failure {iteration} unexpectedly succeeded"
                )
            finally:
                balls._reset_allocation_failure_for_tests()

            assert_no_c_ownership(balls, f"allocation failure {iteration}")

        for failure_point in range(sweep_limit):
            balls.set_default(**parameters)
            balls._set_allocation_failure_after_for_tests(failure_point)
            try:
                balls.Run(level=["StartRun_Common"])
            except CosmoComputationError as exc:
                if "memory allocation failed" not in str(exc):
                    raise AssertionError(
                        f"allocation point {failure_point} returned the wrong "
                        f"diagnostic: {exc}"
                    ) from exc
            else:
                balls.struct_cleanup()
                assert_no_c_ownership(balls, "allocation sweep success")
                break
            finally:
                balls._reset_allocation_failure_for_tests()

            assert_no_c_ownership(
                balls, f"allocation sweep point {failure_point}"
            )
        else:
            raise AssertionError(
                f"allocation sweep did not reach a successful startup in "
                f"{sweep_limit} attempts"
            )

        if failure_point == 0:
            raise AssertionError("startup performed no fault-injectable allocations")

        balls.set_default(**parameters)
        balls.Run(level=["StartRun_Common"])
        balls.struct_cleanup()
        assert_no_c_ownership(balls, "post-failure recovery run")


def test_openmp_allocation_failures_raise_and_recover():
    sweep_limit = int(
        os.environ.get("CBALLS_OPENMP_ALLOCATION_FAILURE_SWEEP_LIMIT", "384")
    )

    with tempfile.TemporaryDirectory(prefix="ctreeballs-p2-openmp-allocation-") as root_dir:
        parameters = {
            "searchMethod": "octree-ggg-omp",
            "testmodel": "simple-cubic",
            "nbody": 8,
            "rootDir": root_dir,
            "options": "no-out-Hist",
            "numberThreads": 2,
            "verbose": 0,
            "verbose_log": 0,
        }

        for method in ("octree-ggg-omp", "kdtree-omp", "kdtree-box-omp"):
            parameters["searchMethod"] = method
            for failure_point in range(sweep_limit):
                balls = cballs()
                balls.set_default(**parameters)
                balls._set_allocation_failure_after_for_tests(failure_point)
                try:
                    balls.Run(level=["MainLoop"])
                except CosmoComputationError as exc:
                    if "allocation failed" not in str(exc):
                        raise AssertionError(
                            f"{method} allocation point {failure_point} "
                            f"returned the wrong diagnostic: {exc}"
                        ) from exc
                else:
                    balls.struct_cleanup()
                    assert_no_c_ownership(
                        balls, f"{method} allocation sweep success"
                    )
                    break
                finally:
                    balls._reset_allocation_failure_for_tests()

                assert_no_c_ownership(
                    balls, f"{method} allocation sweep point {failure_point}"
                )
            else:
                raise AssertionError(
                    f"{method} allocation sweep did not reach a successful "
                    f"run in {sweep_limit} attempts"
                )


def test_degenerate_histogram_domain_is_recoverable():
    with tempfile.TemporaryDirectory(prefix="ctreeballs-p2-cython-") as root_dir:
        balls = cballs()
        balls.set(
            {
                "searchMethod": "octree-ggg-omp",
                "testmodel": "simple-cubic",
                "nbody": 8,
                "rootDir": root_dir,
                "rangeN": 0.1,
                "rminHist": 0.1,
                "useLogHist": False,
                "numberThreads": 2,
                "verbose": 0,
                "verbose_log": 0,
            }
        )

        try:
            balls.Run(level=["StartRun_Common"])
        except CosmoComputationError as exc:
            if "rminHist must be finite" not in str(exc):
                raise AssertionError(f"unexpected Cython diagnostic: {exc}") from exc
        else:
            raise AssertionError("degenerate Cython histogram domain succeeded")

        if balls.getBodytableAllocated():
            raise AssertionError("failed Cython startup retained body-table ownership")


def test_nonfinite_theta_is_recoverable():
    with tempfile.TemporaryDirectory(prefix="ctreeballs-p2-theta-") as root_dir:
        balls = cballs()
        balls.set(
            {
                "searchMethod": "kdtree-omp",
                "testmodel": "simple-cubic",
                "nbody": 8,
                "rootDir": root_dir,
                "theta": math.nan,
                "numberThreads": 2,
                "verbose": 0,
                "verbose_log": 0,
            }
        )

        try:
            balls.Run(level=["StartRun_Common"])
        except (CosmoComputationError, CosmoSevereError) as exc:
            message = str(exc)
            if ("theta must be finite and non-negative" not in message
                    and "invalid numeric value for theta" not in message):
                raise AssertionError(f"unexpected Cython diagnostic: {exc}") from exc
        else:
            raise AssertionError("non-finite Cython theta succeeded")

        if balls.getBodytableAllocated():
            raise AssertionError("failed Cython startup retained body-table ownership")


def test_small_cubic_lattice_is_recoverable():
    with tempfile.TemporaryDirectory(prefix="ctreeballs-p2-cubic-") as root_dir:
        balls = cballs()
        balls.set(
            {
                "searchMethod": "octree-ggg-omp",
                "testmodel": "simple-cubic",
                "nbody": 3,
                "rootDir": root_dir,
                "numberThreads": 2,
                "verbose": 0,
                "verbose_log": 0,
            }
        )

        try:
            balls.Run(level=["StartRun_Common"])
        except CosmoComputationError as exc:
            if "at least two lattice points per dimension" not in str(exc):
                raise AssertionError(f"unexpected Cython diagnostic: {exc}") from exc
        else:
            raise AssertionError("one-point Cython cubic lattice succeeded")

        if balls.getBodytableAllocated():
            raise AssertionError("failed Cython startup retained body-table ownership")


def test_exact_theta_and_short_progress_interval_run():
    for method in ("octree-ggg-omp", "kdtree-omp", "kdtree-box-omp"):
        with tempfile.TemporaryDirectory(prefix="ctreeballs-p2-exact-") as root_dir:
            balls = cballs()
            balls.set(
                {
                    "searchMethod": method,
                    "testmodel": "simple-cubic",
                    "nbody": 8,
                    "rootDir": root_dir,
                    "theta": 0.0,
                    "stepState": 1,
                    "options": "no-out-Hist",
                    "numberThreads": 2,
                    "verbose": 0,
                    "verbose_log": 0,
                }
            )

            try:
                balls.Run(level=["MainLoop"])
            finally:
                balls.struct_cleanup()


if __name__ == "__main__":
    test_allocation_failures_raise_and_recover()
    test_openmp_allocation_failures_raise_and_recover()
    test_degenerate_histogram_domain_is_recoverable()
    test_nonfinite_theta_is_recoverable()
    test_small_cubic_lattice_is_recoverable()
    test_exact_theta_and_short_progress_interval_run()
    print(
        "PASS: Cython recovers allocation failures, guards zero-divisor "
        "domains, and preserves exact mode"
    )
