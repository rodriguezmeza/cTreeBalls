#!/usr/bin/env python3

import os
import re
import tempfile
from pathlib import Path

from cyballs import CosmoComputationError, cballs


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


def balls0357_parameters(root_dir, **overrides):
    parameters = {
        "searchMethod": "balls-omp-0357",
        "testmodel": "simple-cubic",
        "nbody": 8,
        "rootDir": root_dir,
        "options": "no-out-Hist",
        "useLogHist": True,
        "numberThreads": 2,
        "verbose": 0,
        "verbose_log": 0,
    }
    parameters.update(overrides)
    return parameters


def test_addon_contains_no_process_exit():
    addon_dir = Path(__file__).resolve().parents[2] / "addons" / "balls_omp_0357"
    exit_call = re.compile(r"\bexit\s*\(")

    for path in sorted((*addon_dir.glob("*.c"), *addon_dir.glob("*.h"))):
        if exit_call.search(path.read_text(encoding="utf-8")):
            raise AssertionError(f"BALLS0357 still contains a process exit: {path}")


def test_validation_failure_is_recoverable():
    with tempfile.TemporaryDirectory(prefix="ctreeballs-balls0357-status-") as root_dir:
        balls = cballs()
        balls.set_default(
            **balls0357_parameters(root_dir, useLogHist=False)
        )

        try:
            balls.Run(level=["MainLoop"])
        except CosmoComputationError as exc:
            diagnostic = str(exc)
            if "BALLS0357" not in diagnostic or "useLogHist" not in diagnostic:
                raise AssertionError(
                    f"BALLS0357 returned the wrong diagnostic: {diagnostic}"
                ) from exc
        else:
            raise AssertionError("invalid BALLS0357 histogram mode succeeded")

        assert_no_c_ownership(balls, "BALLS0357 validation failure")


def test_worker_allocation_failure_is_recoverable():
    sweep_limit = int(
        os.environ.get("CBALLS_BALLS0357_ALLOCATION_SWEEP_LIMIT", "512")
    )

    with tempfile.TemporaryDirectory(prefix="ctreeballs-balls0357-worker-") as root_dir:
        balls = cballs()
        parameters = balls0357_parameters(root_dir)
        saw_worker_failure = False

        for failure_point in range(sweep_limit):
            balls.set_default(**parameters)
            balls._set_allocation_failure_after_for_tests(failure_point)
            try:
                balls.Run(level=["MainLoop"])
            except CosmoComputationError as exc:
                diagnostic = str(exc)
                if "allocation failed" not in diagnostic:
                    raise AssertionError(
                        f"BALLS0357 allocation point {failure_point} returned "
                        f"the wrong diagnostic: {diagnostic}"
                    ) from exc
                if "OpenMP histogram allocation failed" in diagnostic:
                    saw_worker_failure = True
            else:
                balls.struct_cleanup()
                assert_no_c_ownership(balls, "BALLS0357 allocation sweep success")
                break
            finally:
                balls._reset_allocation_failure_for_tests()

            assert_no_c_ownership(
                balls, f"BALLS0357 allocation sweep point {failure_point}"
            )
        else:
            raise AssertionError(
                "BALLS0357 allocation sweep did not reach a successful run in "
                f"{sweep_limit} attempts"
            )

        if not saw_worker_failure:
            raise AssertionError(
                "BALLS0357 allocation sweep never exercised an OpenMP worker failure"
            )

        balls.set_default(**parameters)
        balls.Run(level=["MainLoop"])
        balls.struct_cleanup()
        assert_no_c_ownership(balls, "BALLS0357 post-failure recovery")


if __name__ == "__main__":
    test_addon_contains_no_process_exit()
    test_validation_failure_is_recoverable()
    test_worker_allocation_failure_is_recoverable()
    print("PASS: BALLS0357 failures propagate and recover without process exit")
