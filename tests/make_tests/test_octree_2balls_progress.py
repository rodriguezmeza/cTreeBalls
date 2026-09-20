#!/usr/bin/env python3
"""Completed-pivot progress, output routing, and histogram invariance."""

import argparse
from pathlib import Path
import re
import subprocess
import tempfile

import numpy as np

from test_octree_2balls_mask import catalog, parameters, write_catalog


ENGINES = ("octree-2balls-omp", "kdtree-2balls-omp", "balltree-2balls-omp")


def progress_pattern(engine):
    return re.compile(re.escape(engine) + r": 3PCF progress: completed pivots (\d+)"
                      r" / (\d+) \(([\d.]+)%\); elapsed ([\d.]+) s")


PROGRESS = progress_pattern("octree-2balls-omp")
FINAL = "3PCF pivots complete; reducing histograms and finalizing outputs"


def check_progress(text, total, interval, engine="octree-2balls-omp"):
    matches = progress_pattern(engine).findall(text)
    assert matches, text
    done = [int(row[0]) for row in matches]
    assert done[0] == 0 and done[-1] == total, done
    assert all(int(row[1]) == total for row in matches), matches
    assert all(left < right for left, right in zip(done, done[1:])), done
    assert all(right - left >= interval or right == total
               for left, right in zip(done, done[1:])), done
    elapsed = [float(row[3]) for row in matches]
    assert elapsed == sorted(elapsed), elapsed
    for completed, _, percent, _ in matches:
        assert abs(float(percent) - 100 * int(completed) / total) <= 0.051
    assert matches[-1][2] == "100.0"
    assert text.count(FINAL) == 1, text
    assert text.index(FINAL) > text.rindex("3PCF progress:"), text
    return matches


def suite(binary, reference_binary=None, engine="octree-2balls-omp"):
    progress = progress_pattern(engine)

    def check(text, total, interval):
        return check_progress(text, total, interval, engine)

    with tempfile.TemporaryDirectory(prefix="two-ball-pivot-progress-") as tmp:
        root = Path(tmp)
        data = catalog(257)
        write_catalog(root / "catalog.txt", data)
        serial = 0

        def run(options, threads=4, verbose=1, verbose_log=1, step=37,
                executable=binary, overrides=None, input_name="catalog.txt",
                neighbor_name=None):
            nonlocal serial
            serial += 1
            out = root / str(serial)
            params = parameters(out, threads, options)
            params.update(verbose=verbose, verbose_log=verbose_log, stepState=step)
            params.update(overrides or {})
            command = [str(executable), f"search={engine}"]
            if neighbor_name:
                command += [f"in={root/input_name},{root/neighbor_name}",
                            "infmt=columns-ascii-all,columns-ascii-all", "iCatalogs=1,2"]
            else:
                command += [f"in={root/input_name}", "infmt=columns-ascii-all"]
            command.extend(
                f"{key}={str(value).lower() if isinstance(value, bool) else value}"
                for key, value in params.items()
            )
            result = subprocess.run(command, capture_output=True, text=True, timeout=120)
            assert result.returncode == 0, result.stdout + result.stderr
            logs = list(out.rglob("cballs.log"))
            assert len(logs) == (1 if verbose_log else 0), logs
            histograms = {p.name: p.read_bytes() for p in out.glob("hist*.txt")}
            assert histograms, result.stdout
            return result.stdout, logs[0].read_text() if logs else "", histograms

        for approximate in (False, True):
            for masked in (False, True):
                options = ["KKKCorrelation", "only-3pcf", "no-smooth-pivot",
                           "no-normalize-HistZeta"]
                if not approximate:
                    options.append("no-two-balls")
                if masked:
                    options.append("read-mask")
                total = int(data[3].sum()) if masked else len(data[0])
                quiet_stdout, quiet_log, expected = run(options, verbose=0, verbose_log=0)
                assert not progress.search(quiet_stdout + quiet_log)
                if reference_binary:
                    _, _, reference = run(options, verbose=0, verbose_log=0,
                                          executable=reference_binary)
                    assert expected == reference, "progress changed the original histograms"
                for threads in (1, 4):
                    stdout, log, histograms = run(options, threads=threads)
                    assert check(stdout, total, 37) == check(log, total, 37)
                    assert histograms == expected, "progress changed the histograms"
                stdout, log, histograms = run(options, verbose=0)
                assert not progress.search(stdout)
                check(log, total, 37)
                assert histograms == expected
                stdout, log, histograms = run(options, verbose_log=0, step=1000000)
                assert len(check(stdout, total, 1000000)) == 2
                assert not progress.search(log)
                assert histograms == expected

        options = ["KKKCorrelation", "only-3pcf", "no-two-balls", "no-smooth-pivot"]
        stdout, log, _ = run(options, step=1)
        assert len(check(stdout, len(data[0]), 1)) == len(data[0]) + 1
        assert progress.findall(stdout) == progress.findall(log)
        for options in (["only-2pcf", "no-two-balls", "no-smooth-pivot"],
                        ["only-3pcf", "no-two-balls", "no-smooth-pivot", "dual-node-direct-triples"]):
            stdout, log, _ = run(options)
            assert not progress.search(stdout + log)
            assert FINAL not in stdout + log

        if engine != "octree-2balls-omp":
            # Exercise partial-frontier completion with several pivots per leaf,
            # including bins completed collectively before children are visited.
            for capacity in (1, 16, 1024):
                options = ["KKKCorrelation", "only-3pcf", "no-smooth-pivot",
                           "no-normalize-HistZeta", "dual-node-bin-slop"]
                stdout, log, histograms = run(options, overrides=dict(nsmooth=capacity), step=1)
                assert check(stdout, len(data[0]), 1) == check(log, len(data[0]), 1)
                _, _, quiet = run(options, overrides=dict(nsmooth=capacity), verbose=0, verbose_log=0)
                assert histograms == quiet
                if reference_binary:
                    _, _, reference = run(options, overrides=dict(nsmooth=capacity),
                                          executable=reference_binary, verbose=0, verbose_log=0)
                    assert histograms == reference

            # An out-of-range neighbor group lets KD finish the entire pivot
            # leaf at once; a larger search range exercises body completions.
            pivot = list(catalog(257))
            neighbor = list(catalog(19))
            pivot[0] = pivot[0]*.001 + [1., 0., 0.]
            neighbor[0] = neighbor[0]*.001 + [-1., 0., 0.]
            write_catalog(root/"pivots.txt", pivot)
            write_catalog(root/"neighbors.txt", neighbor)
            options = ["KKKCorrelation", "only-3pcf", "no-smooth-pivot",
                       "no-normalize-HistZeta"]
            for limit in (.1, 2.5):
                settings = dict(nsmooth=1024, rangeN=limit)
                stdout, log, histograms = run(options, overrides=settings, step=1,
                    input_name="pivots.txt", neighbor_name="neighbors.txt")
                matches = check(stdout, 257, 1)
                assert matches == check(log, 257, 1)
                if engine == "kdtree-2balls-omp" and limit == .1:
                    assert len(matches) == 2, "fully resolved KD group was not counted collectively"
                if reference_binary:
                    _, _, reference = run(options, overrides=settings, executable=reference_binary,
                        input_name="pivots.txt", neighbor_name="neighbors.txt", verbose=0, verbose_log=0)
                    assert histograms == reference

            for extra_options in ([], ["only-3pcf", "no-balltree-persistent-frontier"]):
                options = ["KKKCorrelation", "no-smooth-pivot", *extra_options]
                stdout, log, histograms = run(options)
                assert check(stdout, 257, 37) == check(log, 257, 37)
                if reference_binary:
                    _, _, reference = run(options, verbose=0, verbose_log=0,
                                          executable=reference_binary)
                    assert histograms == reference

            info = subprocess.run([str(binary), "options=make-info"], text=True,
                                  capture_output=True, timeout=30)
            if re.search(r"SMOOTHPIVOTON\s*=\s*1", info.stdout):
                grouped = tuple(np.repeat(a, 3, axis=0) for a in catalog(17))
                write_catalog(root/"grouped.txt", grouped)
                options = ["KKKCorrelation", "only-3pcf", "smooth-pivot",
                           "no-normalize-HistZeta"]
                stdout, log, histograms = run(options, input_name="grouped.txt",
                                              overrides=dict(rsmooth=1e-6), step=1)
                assert check(stdout, 17, 1) == check(log, 17, 1)
                if reference_binary:
                    _, _, reference = run(options, input_name="grouped.txt",
                        overrides=dict(rsmooth=1e-6), verbose=0, verbose_log=0, executable=reference_binary)
                    assert histograms == reference
    print(f"PASS: {engine} completed pivots, masks, threads, logging, and unchanged histograms")


def cython_suite(engines):
    from cyballs import cballs

    positions, kappa, weights, _ = catalog(257)
    with tempfile.TemporaryDirectory(prefix="two-ball-pivot-progress-cython-") as tmp:
        for engine in engines:
            histograms = []
            for iteration in range(2):
                out = Path(tmp) / engine / str(iteration)
                model = cballs()
                params = parameters(out, 4, ["KKKCorrelation", "only-3pcf",
                                            "no-smooth-pivot", "no-out-Hist"])
                params.update(searchMethod=engine, verbose=0, verbose_log=1, stepState=37)
                model.set(params)
                model.set_catalog(positions, kappa=kappa, weights=weights)
                try:
                    model.Run(level=["MainLoop"])
                    histograms.append(model.getHistZetaMsincos(1, 1).copy())
                finally:
                    model.struct_cleanup()
                logs = list(out.rglob("cballs.log"))
                assert len(logs) == 1, logs
                check_progress(logs[0].read_text(), len(positions), 37, engine)
            np.testing.assert_array_equal(*histograms)
            print(f"PASS: {engine} Python progress resets and repeated multipoles agree")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cballs", type=Path, default=Path("./cballs"))
    parser.add_argument("--reference-cballs", type=Path)
    parser.add_argument("--engine", choices=ENGINES, action="append",
                        help="repeat to test multiple methods; default: octree-2balls-omp")
    parser.add_argument("--cython", action="store_true",
                        help="also check repeated calls through the local Python extension")
    args = parser.parse_args()
    engines = args.engine or ["octree-2balls-omp"]
    for engine in engines:
        suite(args.cballs.resolve(),
              args.reference_cballs.resolve() if args.reference_cballs else None, engine)
    if args.cython:
        cython_suite(engines)
