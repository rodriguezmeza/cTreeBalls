#!/usr/bin/env python3
"""MPI regressions for scalar and ENCORE-style survey 3D multipoles."""
import argparse
import json
import os
from pathlib import Path
import re
import shlex
import signal
import subprocess
import sys
import tempfile

import numpy as np
import test_octree_3pcf_3d_omp as scalar
import test_octree_3pcf_3d_survey as survey

METHOD = "octree-3pcf-3d-mpi"
BOTH = "compute-2pcf-3d,compute-3pcf-3d"

def launch(command, timeout=120):
    with subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                          text=True, start_new_session=True) as process:
        try:
            stdout, stderr = process.communicate(timeout=timeout)
        except subprocess.TimeoutExpired:
            os.killpg(process.pid, signal.SIGKILL)
            stdout, stderr = process.communicate()
            raise AssertionError(f"MPI command timed out: {command}\n{stdout}\n{stderr}")
        return subprocess.CompletedProcess(command, process.returncode, stdout, stderr)


def parameters(root, survey_mode=False, threads=1):
    reference = survey if survey_mode else scalar
    return dict(searchMethod=METHOD, rootDir=str(root), numberThreads=threads,
                usePeriodic=False, useLogHist=False, rangeN=reference.RMAX,
                rminHist=reference.RMIN, sizeHistN=reference.NBINS,
                mChebyshev=reference.LMAX, sizeHistPhi=8, lengthBox=3.0,
                stepState=1000000, verbose=0, verbose_log=0,
                iCatalogs="1,2" if survey_mode else "1",
                options=BOTH + (",survey-estimator-3d" if survey_mode else ""))


def products(root):
    return {p.name: p for p in root.glob("*_3d*.txt")}


def compare(left, right, exact=False):
    a, b = products(left), products(right)
    assert a.keys() == b.keys(), (a.keys(), b.keys())
    for name in a:
        if exact:
            assert a[name].read_bytes() == b[name].read_bytes(), name
        else:
            np.testing.assert_allclose(np.loadtxt(a[name]), np.loadtxt(b[name]),
                                       rtol=3e-10, atol=3e-10, err_msg=name)
            counts = lambda p: re.findall(
                r"(?:npivots|accepted_neighbor_visits) (\d+)", p.read_text())
            assert counts(a[name]) == counts(b[name])


def check_scalar(root, los_ids=None):
    if los_ids is None:
        scalar.check_oracle(root)
    xn, xd, zn, zd = scalar.direct_oracle(los_ids)
    xi = np.loadtxt(root/"histXi2pcf_3d.txt")
    zeta = np.loadtxt(root/"histZetaM_3d.txt")
    for actual, expected in ((xi[:, 3], xn), (xi[:, 4], xd),
                             (zeta[:, 6], np.asarray(zn).ravel()),
                             (zeta[:, 7], np.asarray(zd).ravel())):
        np.testing.assert_allclose(actual, expected, rtol=3e-11, atol=3e-11)


def fixtures(root):
    np.savetxt(root/"scalar.txt", scalar.CATALOG, fmt="%.17g")
    survey.write_catalog(root/"data.txt", survey.DATA)
    survey.write_catalog(root/"random.txt", survey.RANDOMS)


def c_tests(binary, mpi, mpi_only, fits=False):
    with tempfile.TemporaryDirectory(prefix="cb3d-mpi-") as tmp:
        root = Path(tmp)
        fixtures(root)
        sequence = 0

        def run(survey_mode=False, threads=1, extra=None, fail=False, overrides=None):
            nonlocal sequence
            sequence += 1
            out = root/str(sequence)
            p = parameters(out, survey_mode, threads)
            p.update(infile=(str(root/"data.txt")+","+str(root/"random.txt")
                             if survey_mode else str(root/"scalar.txt")),
                     infileformat=("multi-columns-ascii,multi-columns-ascii"
                                   if survey_mode else "multi-columns-ascii"),
                     columns="1,2,3,4,5")
            p["options"] += ",pos-and-convergence-weight"
            p.update(extra or {})
            out = Path(p["rootDir"])
            command = [str(binary)] + [f"{k}={str(v).lower() if isinstance(v, bool) else v}"
                                      for k, v in p.items()]
            if overrides:
                command = [sys.executable, str(Path(__file__).resolve()), "--rank-launch",
                           json.dumps([command, overrides])]
            if p["searchMethod"].endswith("-mpi"):
                command = mpi + command
            proc = launch(command)
            assert (proc.returncode != 0) if fail else proc.returncode == 0, proc.stdout+proc.stderr
            return out

        for mode, oracle in ((False, check_scalar), (True, survey.check_oracle)):
            one = run(mode, 1)
            oracle(one)
            compare(one, run(mode, 3), exact=True)
            compare(one, run(mode, 2, {"searchMethod": "octree-ggg-3d-mpi"}), exact=True)
            if not mpi_only:
                compare(one, run(mode, 3, {"searchMethod": "octree-3pcf-3d-omp"}))
            for option, prefix in (("only-2pcf-3d", "histXi2pcf"),
                                   ("only-3pcf-3d", "histZetaM")):
                options = option+",pos-and-convergence-weight"
                if mode:
                    options += ",survey-estimator-3d"
                out = run(mode, extra={"options": options})
                assert len(products(out)) == 1
                assert next(iter(products(out))).startswith(prefix)
            run(mode, extra={"options": "only-2pcf-3d,only-3pcf-3d"}, fail=True)
            suffix = ",survey-estimator-3d" if mode else ""
            silent = run(mode, extra={"options": BOTH+suffix+",no-out-Hist,pos-and-convergence-weight"})
            assert not products(silent)
            empty = run(mode, extra={"rminHist": 10.0, "rangeN": 20.0})
            for name, path in products(empty).items():
                values = np.loadtxt(path)
                if mode and name.startswith("histZetaM"):
                    assert np.isfinite(np.delete(values, 12, axis=1)).all()
                    assert np.isposinf(values[:, 12]).all()
                else:
                    assert np.isfinite(values).all()
                if name.startswith("histXi2pcf"):
                    np.testing.assert_array_equal(values[:, 2:5], 0)
                    if mode:
                        np.testing.assert_array_equal(values[:, 5], 0)
                elif mode:
                    np.testing.assert_array_equal(values[:, 5:12], 0)
                    np.testing.assert_array_equal(values[:, 13], 0)
                else:
                    np.testing.assert_array_equal(values[:, 5:], 0)
            nonroot = root/f"nonroot-{mode}"
            run(mode, overrides={"1": {"rootDir": str(nonroot)}})
            assert not nonroot.exists(), "non-root rank wrote output"
            run(mode, fail=True, overrides={"1": {"mChebyshev": "33"}})
            run(mode, fail=True, overrides={"1": {"infile": str(root/"missing.txt")}})
            bad = root/f"bad-output-{mode}"
            bad.mkdir()
            (bad/("histZetaM_3d_survey.txt" if mode else "histZetaM_3d.txt")).mkdir()
            run(mode, extra={"rootDir": str(bad)}, fail=True)
            print(f"PASS: {'survey edge-corrected' if mode else 'scalar'} MPI oracles, modes, errors and threads", flush=True)

        for theta in (0.0, 3.0):
            check_scalar(run(extra={"theta": theta, "options": BOTH+",Rcut/theta,pos-and-convergence-weight"}))

        if fits:
            from astropy.io import fits as fits_io
            los_ids = (10, 10, 20, 20, 30, 30)
            table = np.array([tuple(row)+(los,) for row, los in zip(scalar.CATALOG, los_ids)],
                             dtype=[("X", "f8"), ("Y", "f8"), ("Z", "f8"),
                                    ("DELTA", "f8"), ("WEIGHT", "f8"), ("LOS_ID", "i8")])
            fits_io.BinTableHDU(table).writeto(root/"los.fits")
            out = run(extra={"infile": str(root/"los.fits"), "infileformat": "fits",
                             "columns": "1,2,3,4,5,6",
                             "options": BOTH+",with-weight,exclude-same-los"})
            check_scalar(out, los_ids)
            print("PASS: MPI FITS 2PCF/3PCF and LOS exclusion oracle", flush=True)

        # Enough pivot blocks to put real work on every rank in the two-rank test.
        rng = np.random.default_rng(90211)
        large = np.column_stack((rng.uniform(-.7, .7, (257, 3)),
                                 rng.normal(size=257), rng.uniform(.5, 1.5, 257)))
        np.savetxt(root/"large.txt", large, fmt="%.17g")
        randoms = np.column_stack((rng.uniform(-.9, .9, (311, 3)),
                                   np.ones(311), rng.uniform(.5, 1.5, 311)))
        np.savetxt(root/"large-random.txt", randoms, fmt="%.17g")
        for mode in (False, True):
            extra = {"infile": str(root/"large.txt")}
            if mode:
                extra["infile"] = str(root/"large.txt")+","+str(root/"large-random.txt")
            one = run(mode, 1, extra)
            compare(one, run(mode, 3, extra), exact=True)
            if not mpi_only:
                compare(one, run(mode, 3, dict(extra, searchMethod="octree-3pcf-3d-omp")))
            if mode:
                zeta = np.loadtxt(one/"histZetaM_3d_survey.txt")
                assert np.isfinite(zeta).all()
                assert (zeta[:, 13] == 1).all()
        if not mpi_only:
            for extra in ({"useLogHist": True}, {"usePeriodic": True, "lengthBox": 3.0}):
                extra["infile"] = str(root/"large.txt")
                compare(run(extra=extra),
                        run(extra=dict(extra, searchMethod="octree-3pcf-3d-omp"), threads=3))
        print("PASS: multi-block MPI/OpenMP agreement, logarithmic/periodic bins, Rcut/theta", flush=True)


def cython_worker(root):
    from mpi4py import MPI
    import cyballs
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    root = Path(root)
    assert cyballs.search_method_id(METHOD) == 192
    assert cyballs.search_method_id("octree-ggg-3d-mpi") == 192
    for mode, oracle in ((False, check_scalar), (True, survey.check_oracle)):
        b = cyballs.cballs()
        p = parameters(root/f"python-{mode}", mode, threads=2)
        b.set(p)
        if mode:
            for index, rows in enumerate((survey.DATA, survey.RANDOMS)):
                rows = np.asarray(rows)
                b.set_catalog(rows[:, :3], weights=rows[:, 3], catalog=index)
        else:
            rows = np.asarray(scalar.CATALOG)
            b.set_catalog(rows[:, :3], kappa=rows[:, 3], weights=rows[:, 4])
        try:
            b.Run(level=["MainLoop"])
            if rank == 0:
                oracle(Path(p["rootDir"]))
        finally:
            b.struct_cleanup()
        comm.Barrier()

    # Masking in memory must equal a filtered auto catalog.
    rows = np.asarray(scalar.CATALOG)
    mask = np.array([1, 0, 1, 1, 0, 1], dtype=np.uint8)
    for masked in (True, False):
        b = cyballs.cballs()
        p = parameters(root/f"mask-{masked}", threads=2)
        if masked:
            p["options"] += ",read-mask"
            data = rows
        else:
            data = rows[mask.astype(bool)]
        b.set(p)
        b.set_catalog(data[:, :3], kappa=data[:, 3], weights=data[:, 4],
                      mask=mask if masked else None)
        try:
            b.Run(level=["MainLoop"])
        finally:
            b.struct_cleanup()
    if rank == 0:
        compare(root/"mask-True", root/"mask-False")

    b = cyballs.cballs()
    for attempt in range(4):
        b.clear_catalogs()
        p = parameters(root/f"recover-{attempt}", True, threads=2)
        b.set(p)
        for index, rows in enumerate((survey.DATA, survey.RANDOMS)):
            data = np.array(rows)
            if attempt < 3 and rank == 1 and index == 1:
                data[0, 3] = -1
            b.set_catalog(data[:, :3], weights=data[:, 3], catalog=index)
        failed = False
        try:
            b.Run(level=["MainLoop"])
        except Exception:
            failed = True
        finally:
            b.struct_cleanup()
        assert all(value == (attempt < 3) for value in comm.allgather(failed))
    if rank == 0:
        survey.check_oracle(Path(p["rootDir"]))
        print("PASS: MPI Cython in-memory scalar/survey, masks, repeated failure recovery", flush=True)
    comm.Barrier()
    MPI.Finalize()
    try:
        b.Run(level=["MainLoop"])
    except Exception as exc:
        assert "MPI_Finalize" in str(exc), str(exc)
    else:
        raise AssertionError("cyballs accepted reuse after MPI_Finalize")
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
        command = [x for x in command if x.split("=", 1)[0] not in changes]
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
    parser.add_argument("--fits", action="store_true",
                        help="also check FITS LOS input; requires CFITSIO and Astropy")
    args = parser.parse_args()
    mpi = shlex.split(args.mpi_command)
    c_tests(args.cballs.resolve(), mpi, args.mpi_only, args.fits)
    if args.cython:
        with tempfile.TemporaryDirectory(prefix="cb3d-mpi-python-") as tmp:
            result = launch([*mpi, sys.executable, str(Path(__file__).resolve()),
                             "--cython-worker", tmp], timeout=180)
            assert result.returncode == 0, result.stdout+result.stderr
            print(result.stdout, end="")


if __name__ == "__main__":
    main()
