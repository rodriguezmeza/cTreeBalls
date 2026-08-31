#!/usr/bin/env python3
"""Independent Ly-alpha oracles and OpenMP/MPI comparisons for all seven modes."""
import argparse
import json
import os
from pathlib import Path
import re
import shlex
import subprocess
import sys
import tempfile

import numpy as np
import test_lya_forest_omp as three
import test_lya_forest_1d_omp as radial

METHODS = ("lya-2pcf", "lya-3pcf", "lya-2pcf-3pcf", "lya-1d-2pcf",
           "lya-1d-3pcf", "lya-1d-2pcf-3pcf", "lya-1d-tree-2pcf")


def params(kind, catalog, output, threads=1):
    p = dict(infile=str(catalog), infileformat="lya-ascii", iCatalogs="1",
             rootDir=str(output), numberThreads=threads, usePeriodic=False,
             useLogHist=False, rangeN=30.0, rminHist=0.1, sizeHistN=4,
             verbose=0, verbose_log=0, options="")
    if kind < 3:
        p.update(lya2RpMax=three.RP_MAX, lya2RtMax=three.RT_MAX,
                 lya2RpBins=three.RP_BINS, lya2RtBins=three.RT_BINS,
                 lya3RMax=three.R3_MAX, lya3RBins=three.R3_BINS,
                 lya3ThetaBins=three.THETA_BINS, lya3MuBins=three.MU_BINS)
    else:
        p.update(lya2RpMax=radial.RP_MAX, lya2RpBins=radial.RP_BINS,
                 lya3RMax=radial.R3_MAX, lya3RBins=radial.R3_BINS_PER_SIDE)
    return p


def products(kind):
    return tuple(name for condition, name in (
        (kind in (0, 2), "histXi2pcf_lya.txt"),
        (kind in (1, 2), "histZetaM_lya5d.txt"),
        (kind in (3, 5, 6), "histXi2pcf_lya1d.txt"),
        (kind in (4, 5), "histZetaM_lya1d.txt")) if condition)


def check_oracle(kind, output):
    if kind in (0, 2):
        three.assert_histogram_close(three.read_2pcf(output/products(kind)[0]),
                                     three.oracle_2pcf(), "3D 2PCF")
    if kind in (1, 2):
        three.assert_histogram_close(three.read_3pcf(output/"histZetaM_lya5d.txt"),
                                     three.oracle_3pcf(), "3D 3PCF")
    if kind in (3, 5, 6):
        radial.assert_histogram_close(radial.read_2pcf(output/"histXi2pcf_lya1d.txt"),
                                      radial.oracle_2pcf(), "radial 2PCF")
    if kind in (4, 5):
        radial.assert_histogram_close(radial.read_3pcf(output/"histZetaM_lya1d.txt"),
                                      radial.oracle_3pcf(), "radial 3PCF")
    for name in products(kind):
        a = np.atleast_2d(np.loadtxt(output/name))
        np.testing.assert_allclose(a[:, -3], a[:, -2]/np.where(a[:, -1] != 0,
                                                            a[:, -1], 1),
                                   rtol=2e-13, atol=2e-13)
        assert np.all(np.isfinite(a))


def compare(kind, left, right, exact=False):
    for name in products(kind):
        if exact:
            assert (left/name).read_bytes() == (right/name).read_bytes(), name
        else:
            np.testing.assert_allclose(np.loadtxt(left/name), np.loadtxt(right/name),
                                       rtol=3e-12, atol=3e-12, err_msg=name)
            # Counts in file headers must agree exactly.
            count = lambda p: re.findall(r"(?:pairs|triplets): (\d+)", p.read_text())
            assert count(left/name) == count(right/name)


def write_catalog(path, data):
    np.savetxt(path, data, fmt=("%.17g",)*5 + ("%.0f",))


def c_tests(binary, mpi, mpi_only):
    with tempfile.TemporaryDirectory(prefix="lya-mpi-") as tmp:
        root = Path(tmp)
        c3, c1 = root/"three.txt", root/"radial.txt"
        write_catalog(c3, three.POINTS)
        write_catalog(c1, radial.WIDE_ANGLE)
        sequence = 0

        def run(kind, suffix="mpi", threads=1, catalog=None, extra=None,
                fail=False, overrides=None):
            nonlocal sequence
            sequence += 1
            out = root/str(sequence)
            p = params(kind, catalog or (c3 if kind < 3 else c1), out, threads)
            p.update(extra or {})
            p["searchMethod"] = METHODS[kind]+"-"+suffix
            command = [str(binary)] + [f"{k}={str(v).lower() if isinstance(v, bool) else v}"
                                      for k, v in p.items()]
            if overrides:
                command = [sys.executable, str(Path(__file__).resolve()),
                           "--rank-launch", json.dumps([command, overrides])]
            command = (mpi if suffix == "mpi" else []) + command
            proc = subprocess.run(command, capture_output=True, text=True, timeout=120)
            assert (proc.returncode != 0) if fail else (proc.returncode == 0), proc.stdout+proc.stderr
            if fail:
                return proc.stdout+proc.stderr
            expected = set() if "no-out-Hist" in p["options"] else set(products(kind))
            assert {x.name for x in out.glob("*lya*.txt")} == expected
            return out

        for kind, name in enumerate(METHODS):
            one = run(kind)
            many = run(kind, threads=3)
            check_oracle(kind, one)
            compare(kind, one, many, exact=True)
            if not mpi_only:
                omp = run(kind, suffix="omp", threads=3)
                compare(kind, omp, one)
            for extra in (dict(lya2RpMax=0) if kind not in (1, 4) else dict(lya3RMax=0),
                          dict(usePeriodic=True)):
                run(kind, extra=extra, fail=True)
            if kind < 3:
                for theta in (0.0, 3.0):
                    check_oracle(kind, run(kind, extra=dict(theta=theta)))
            print(f"PASS: {name}-mpi oracle and thread/rank agreement", flush=True)

        # More than two radial blocks, plus enough interval-tree tasks for both ranks.
        rng = np.random.default_rng(8903)
        radii = rng.uniform(20, 70, 4096)
        large = np.column_stack((radii, np.zeros((len(radii), 2)),
                                 rng.normal(size=len(radii)), rng.uniform(.3, 2, len(radii)),
                                 np.arange(len(radii)) % 53))
        large_file = root/"large.txt"
        write_catalog(large_file, large)
        scan = run(3, catalog=large_file, threads=3)
        tree = run(6, catalog=large_file, threads=3)
        np.testing.assert_allclose(np.loadtxt(scan/"histXi2pcf_lya1d.txt"),
                                   np.loadtxt(tree/"histXi2pcf_lya1d.txt"), rtol=3e-12, atol=3e-10)
        if not mpi_only:
            for kind, actual in ((3, scan), (6, tree)):
                compare(kind, actual, run(kind, suffix="omp", catalog=large_file, threads=1))
        triple_file = root/"triples.txt"
        write_catalog(triple_file, large[:640])
        for kind in (4, 5):
            one = run(kind, catalog=triple_file, threads=1)
            many = run(kind, catalog=triple_file, threads=3)
            compare(kind, one, many, exact=True)
            if not mpi_only:
                compare(kind, one, run(kind, suffix="omp", catalog=triple_file, threads=3))
        one_forest = large[:64].copy()
        one_forest[:, 5] = 0
        write_catalog(root/"same.txt", one_forest)
        for kind in (2, 5, 6):
            empty = run(kind, catalog=root/"same.txt", extra={"options": "lya-output-empty-bins"})
            for name in products(kind):
                np.testing.assert_array_equal(np.atleast_2d(np.loadtxt(empty/name))[:, -3:], 0)

        for kind in (2, 5, 6):
            nonroot = root/f"nonroot-{kind}"
            run(kind, overrides={"1": {"rootDir": str(nonroot)}})
            assert not nonroot.exists(), "non-root rank wrote output"
            run(kind, extra={"options": "no-out-Hist"})
            run(kind, fail=True, overrides={"1": {"lya2RpMax": "0"}})
            # Only rank 1 fails input. Other ranks must return, not wait in a later reduction.
            message = run(kind, fail=True, overrides={"1": {"infile": str(root/"missing.txt")}})
            assert "Ly-alpha catalog input" in message or "lya-ascii" in message
            bad_out = root/f"bad-{kind}"
            bad_out.mkdir()
            (bad_out/products(kind)[0]).mkdir()
            run(kind, extra={"rootDir": str(bad_out)}, fail=True)
        print("PASS: radial scan/tree, same-quasar exclusion, rank-local input and root output failures", flush=True)


def cython_worker(root):
    from mpi4py import MPI
    import cyballs
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    root = Path(root)
    b = cyballs.cballs()
    for kind, name in enumerate(METHODS):
        assert cyballs.search_method_id(name+"-mpi") == 185+kind
        out = root/f"py-{kind}"
        p = params(kind, root/("three.txt" if kind < 3 else "radial.txt"), out, threads=2)
        p["searchMethod"] = name+"-mpi"
        try:
            b.set(p)
            b.Run(level=["MainLoop"])
            if rank == 0:
                check_oracle(kind, out)
        finally:
            b.struct_cleanup()
        comm.Barrier()
    for _ in range(3):
        p = params(2, root/("missing.txt" if rank == 1 else "three.txt"), root/"py-bad")
        p["searchMethod"] = "lya-2pcf-3pcf-mpi"
        failed = False
        try:
            b.set(p)
            b.Run(level=["MainLoop"])
        except Exception:
            failed = True
        finally:
            b.struct_cleanup()
        assert all(comm.allgather(failed))
    p = params(2, root/"three.txt", root/"py-recovered", threads=2)
    p["searchMethod"] = "lya-2pcf-3pcf-mpi"
    b.set(p)
    b.Run(level=["MainLoop"])
    b.struct_cleanup()
    if rank == 0:
        check_oracle(2, root/"py-recovered")
        print("PASS: all seven Cython MPI methods and repeated failure recovery", flush=True)
    comm.Barrier()
    MPI.Finalize()
    try:
        b.set(p)
        b.Run(level=["MainLoop"])
    except Exception as exc:
        assert "MPI_Finalize" in str(exc), str(exc)
    else:
        raise AssertionError("cyballs accepted a run after MPI_Finalize")
    finally:
        b.struct_cleanup()
    if rank == 0:
        print("PASS: Cython rejects reuse after MPI_Finalize", flush=True)


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "--rank-launch":
        command, overrides = json.loads(sys.argv[2])
        rank = next((os.environ[k] for k in ("OMPI_COMM_WORLD_RANK", "PMI_RANK", "PMIX_RANK")
                     if k in os.environ), None)
        if rank is None:
            raise RuntimeError("MPI launcher did not provide a rank")
        changes = overrides.get(rank, {})
        command = [arg for arg in command if arg.split("=", 1)[0] not in changes]
        command.extend(f"{k}={v}" for k, v in changes.items())
        os.execv(command[0], command)
    if len(sys.argv) > 1 and sys.argv[1] == "--cython-worker":
        cython_worker(sys.argv[2])
        return
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cballs", type=Path, default=Path("./cballs"))
    parser.add_argument("--mpi-command", default="mpiexec -n 2")
    parser.add_argument("--mpi-only", action="store_true")
    parser.add_argument("--cython", action="store_true")
    args = parser.parse_args()
    mpi = shlex.split(args.mpi_command)
    c_tests(args.cballs.resolve(), mpi, args.mpi_only)
    if args.cython:
        with tempfile.TemporaryDirectory(prefix="lya-mpi-python-") as tmp:
            root = Path(tmp)
            write_catalog(root/"three.txt", three.POINTS)
            write_catalog(root/"radial.txt", radial.WIDE_ANGLE)
            subprocess.run([*mpi, sys.executable, str(Path(__file__).resolve()),
                            "--cython-worker", tmp], check=True, timeout=120)


if __name__ == "__main__":
    main()
