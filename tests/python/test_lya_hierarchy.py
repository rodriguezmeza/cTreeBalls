#!/usr/bin/env python3
"""Exact hierarchy contracts for all six 3D Ly-alpha engines, active under -O.

Native independent triangle oracles, adversarial catalogs, MPI partitioning,
range sums, certified mixed-forest moments, and explicit approximation checks.
"""
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
import test_lya_forest_omp as oracle

ROOT = Path(__file__).resolve().parents[2]
PRODUCTS = {2: "histXi2pcf_lya.txt", 3: "histZetaM_lya5d.txt"}


def require(value, message):
    if not value:
        raise AssertionError(message)


def forests(nf=6, npix=12, seed=989, clustered=False):
    rng = np.random.default_rng(seed)
    if clustered:
        centers = rng.uniform([90, -10, -20], [150, 30, 20], (nf, 3))
        xyz = centers[:, None, :] + 0.0001 * rng.normal(size=(nf, npix, 3))
    else:
        directions = np.column_stack((np.ones(nf), rng.uniform(-0.13, 0.13, (nf, 2))))
        directions /= np.linalg.norm(directions, axis=1)[:, None]
        xyz = directions[:, None, :] * rng.uniform(80, 260, (nf, npix, 1))
    result = np.column_stack(
        (
            xyz.reshape(-1, 3),
            rng.normal(size=nf * npix),
            rng.uniform(0.1, 2, nf * npix),
            np.repeat(np.arange(nf) * 37 - 2**40, npix),
        )
    )
    rng.shuffle(result)
    return result


def compare(left, right, orders, exact=False):
    for order in orders:
        a, b = left[order], right[order]
        np.testing.assert_array_equal(
            a[:, : 4 if order == 2 else 10], b[:, : 4 if order == 2 else 10]
        )
        np.testing.assert_array_equal(a[:, -1] > 0, b[:, -1] > 0)
        if exact:
            np.testing.assert_array_equal(a, b)
        else:
            np.testing.assert_allclose(a[:, -2:], b[:, -2:], rtol=3e-12, atol=3e-12)
            np.testing.assert_allclose(a[:, -3], b[:, -3], rtol=3e-11, atol=3e-12)
        require(
            left["counts"][order] == right["counts"][order], "geometric count mismatch"
        )


def cython_mpi_contract(los_tree=False):
    from mpi4py import MPI

    sys.path.insert(0, str(ROOT))
    import cyballs

    comm = MPI.COMM_WORLD
    with tempfile.TemporaryDirectory(prefix="lya-hierarchy-mpi-") as tmp:
        tmp = comm.bcast(tmp if comm.rank == 0 else None, root=0)
        params = dict(
            searchMethod="lya-los-tree-2pcf-3pcf-mpi" if los_tree else "lya-2pcf-3pcf-mpi",
            rootDir=tmp,
            numberThreads=2,
            usePeriodic=False,
            useLogHist=False,
            rangeN=30.0,
            rminHist=0.1,
            sizeHistN=4,
            lya2RpMax=30.0,
            lya2RtMax=30.0,
            lya2RpBins=5,
            lya2RtBins=6,
            lya3RMax=30.0,
            lya3RBins=4,
            lya3ThetaBins=5,
            lya3MuBins=6,
            lya2Kernel=1,
            lya3Kernel=5,
            options="no-smooth-pivot,no-out-Hist",
            verbose=0,
            verbose_log=0,
        )
        for case in ("rank-controls", "rank-budget", "recover"):
            m = cyballs.cballs()
            settings = dict(params)
            if case == "rank-controls" and comm.rank == 1:
                settings["lya3Kernel"] = 4
            old = os.environ.get("CBALLS_MEMORY_BUDGET_MB")
            if case == "rank-budget":
                os.environ["CBALLS_MEMORY_BUDGET_MB"] = (
                    ".01" if comm.rank == 1 else "2048"
                )
            try:
                m.set(settings)
                p = oracle.POINTS
                m.set_forest_catalog(
                    p[:, :3], p[:, 3], p[:, 4], p[:, 5].astype(np.int64)
                )
                error = None
                try:
                    m.Run(level=["MainLoop"])
                except Exception as exc:
                    error = str(exc)
                errors = comm.allgather(error)
                if case != "recover":
                    require(all(errors), (case, errors))
                else:
                    require(not any(errors), errors)
                    if comm.rank == 0:
                        require(
                            m.getRunMetadata()["lya_hierarchy"]["enabled"],
                            "missing hierarchy metadata",
                        )
            finally:
                if old is None:
                    os.environ.pop("CBALLS_MEMORY_BUDGET_MB", None)
                else:
                    os.environ["CBALLS_MEMORY_BUDGET_MB"] = old
                m.struct_cleanup()
            comm.Barrier()
        if comm.rank == 0:
            print(
                "PASS rank-local control/budget failures and Cython MPI recovery",
                flush=True,
            )


def main():
    if "--cython-mpi-worker" in sys.argv:
        cython_mpi_contract("--los-tree" in sys.argv)
        return
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--los-tree", action="store_true")
    parser.add_argument("--cballs", type=Path, default=ROOT / "cballs")
    parser.add_argument(
        "--mpi-command", help="launcher including -n, e.g. mpiexec -n 2"
    )
    args = parser.parse_args()
    serial = 0
    with tempfile.TemporaryDirectory(prefix="lya-hierarchy-test-") as tmp:
        root = Path(tmp)

        def run(
            points, stat="2pcf-3pcf", kernel=5, pair=1, threads=1, mpi=False, extra=None
        ):
            nonlocal serial
            serial += 1
            out = root / str(serial)
            out.mkdir()
            path = out / "forest.txt"
            np.savetxt(path, points, fmt=("%.17g",) * 5 + ("%.0f",))
            settings = dict(
                search=("lya-los-tree-" if args.los_tree else "lya-") + stat + "-" + ("mpi" if mpi else "omp"),
                infile=path,
                infileformat="lya-ascii",
                iCatalogs=1,
                rootDir=out,
                numberThreads=threads,
                usePeriodic="false",
                useLogHist="false",
                rangeN=160,
                rminHist=0.1,
                sizeHistN=4,
                lya2RpMax=30,
                lya2RtMax=30,
                lya2RpBins=5,
                lya2RtBins=6,
                lya3RMax=160,
                lya3RBins=4,
                lya3ThetaBins=5,
                lya3MuBins=6,
                lya3Kernel=0 if stat == "2pcf" else kernel,
                lya2Kernel=0 if stat == "3pcf" else pair,
                verbose=2,
                verbose_log=0,
                options="no-smooth-pivot,lya-output-empty-bins",
            )
            settings.update(extra or {})
            command = ([*shlex.split(args.mpi_command)] if mpi else []) + [
                str(args.cballs.resolve())
            ]
            command += [f"{k}={v}" for k, v in settings.items()]
            res = subprocess.run(command, capture_output=True, text=True, timeout=120)
            require(res.returncode == 0, res.stdout + res.stderr)
            result = {
                "counts": {},
                "log": res.stdout,
                "metadata": json.loads((out / "run-metadata.json").read_text()),
            }
            for order in (2, 3):
                f = out / PRODUCTS[order]
                if f.exists():
                    result[order] = np.loadtxt(f, ndmin=2)
                    result["counts"][order] = int(
                        re.search(r"(?:pairs|triplets): (\d+)", f.read_text())[1]
                    )
            return result

        # Independent NumPy oracle, including permutations and all forest vetoes.
        for stat in ("2pcf", "3pcf", "2pcf-3pcf"):
            orders = (2, 3) if stat == "2pcf-3pcf" else (int(stat[0]),)
            a = run(oracle.POINTS, stat, extra={"lya3RMax": 30})
            for order in orders:
                reader = oracle.oracle_2pcf if order == 2 else oracle.oracle_3pcf
                expected = reader()
                table = a[order]
                actual = {
                    tuple(row[: 2 if order == 2 else 5].astype(int)): tuple(row[-2:])
                    for row in table
                    if row[-1] > 0
                }
                oracle.assert_histogram_close(
                    actual, expected, "hierarchy independent oracle"
                )
            if args.mpi_command:
                for kernel in (3, 4, 5):
                    b = run(
                        oracle.POINTS,
                        stat,
                        kernel=kernel,
                        mpi=True,
                        extra={"lya3RMax": 30},
                    )
                    compare(a, b, orders)
        print("PASS six methods / independent pair and triangle oracles", flush=True)

        for case in (
            "ordinary",
            "clustered",
            "bent",
            "closed",
            "boundary",
            "zero",
            "dominant",
            "one-forest",
            "two-forests",
            "translated",
            "small-weight",
            "split-forest-collisions",
        ):
            points = forests(clustered=case == "clustered")
            if case == "bent":
                points[::3, :3] += np.array([2, 8, -9])
            if case == "closed":
                points[:, :3] = [100, 20, 0] + 20 * np.column_stack(
                    (
                        np.cos(np.arange(len(points))),
                        np.sin(np.arange(len(points))),
                        np.zeros(len(points)),
                    )
                )
            if case == "boundary":
                points[:, :3] = [100, 0, 0] + np.arange(len(points))[
                    :, None
                ] * np.array([10, 0, 0])
            if case == "zero":
                points[::3, 4] = 0
            if case == "dominant":
                points[points[:, 5] == points[0, 5], 4] = 1e12
            if case == "one-forest":
                points[:, 5] = 17
            if case == "two-forests":
                points[:, 5] = np.arange(len(points)) % 2
            if case == "translated":
                points[:, :3] += np.array([1e8, -2e8, 3e8])
            if case == "split-forest-collisions":
                points = forests(clustered=True, nf=6, npix=16)
                ids = np.unique(points[:, 5])
                colliding = []
                for fid in range(100000):
                    key = (fid * 0x9E3779B97F4A7C15) & ((1 << 64) - 1)
                    if (key ^ (key >> 33)) & 255 == 0:
                        colliding.append(fid)
                    if len(colliding) == len(ids):
                        break
                for j, fid in enumerate(ids):
                    indices = np.flatnonzero(points[:, 5] == fid)
                    points[indices[:8], :3] += np.array([50.0, 20.0, 10.0])
                    points[indices, 5] = colliding[j]
            if case == "small-weight":
                points[:, 4] = 1e-110
            ref = run(points, kernel=1, pair=0)
            one = run(points)
            many = run(points, threads=3)
            compare(ref, one, (2, 3))
            compare(one, many, (2, 3), exact=True)
            if args.mpi_command:
                mpi = run(points, threads=2, mpi=True)
                compare(one, mpi, (2, 3))
            if case == "ordinary":
                require(
                    one["metadata"]["lya_hierarchy"]["pixel_fallback"],
                    "sparse pivot fallback not exercised",
                )
            if case == "split-forest-collisions":
                require(
                    not one["metadata"]["lya_hierarchy"]["pixel_fallback"],
                    "collision fixture bypassed the hierarchy",
                )
            if case == "clustered":
                info = one["metadata"]["lya_hierarchy"]
                require(
                    info["enabled"]
                    and info["certificates"] > 0
                    and info["represented"] > info["certificates"],
                    info,
                )
                pair = run(points, "2pcf")
                triple = run(points, "3pcf")
                compare(one, pair, (2,))
                compare(one, triple, (3,))
            require(
                not one["metadata"]["lya_geometry"]["three_point_approximate"],
                "zero slop changed estimator",
            )
        print(
            "PASS adversarial geometry, zero/dominant weights, forest exclusions and threads",
            flush=True,
        )

        if args.los_tree:
            # Many small forests force the sparse fallback while retaining
            # enough roots to exercise the mixed-forest per-pivot hierarchy.
            for case in ("ordinary", "colliding", "dominant", "translated", "rotated"):
                points=forests(nf=48,npix=3)
                if case=="rotated": points[:,:3]=points[:,[2,0,1]]
                if case=="dominant": points[::9,4]=1e12
                if case=="translated": points[:,:3]+=np.array([1e8,-2e8,3e8])
                if case=="colliding":
                    ids=np.unique(points[:,5]); colliding=[]
                    for fid in range(100000):
                        key=(fid*0x9E3779B97F4A7C15)&((1<<64)-1)
                        if (key^(key>>33))&255==0: colliding.append(fid)
                        if len(colliding)==len(ids): break
                    for old,new in zip(ids,colliding): points[points[:,5]==old,5]=new
                ref=run(points,kernel=1,pair=0,extra={"search":"lya-2pcf-3pcf-omp"})
                fast=run(points)
                compare(ref,fast,(2,3));compare(fast,run(points,threads=4),(2,3),exact=True)
                info=fast["metadata"]["lya_hierarchy"]
                require(info["pixel_fallback"] and info["nodes"]>0,(case,info))
                if case=="ordinary": require(info["certificates"]>0,info)
                if args.mpi_command: compare(fast,run(points,mpi=True,threads=2),(2,3))
            # Many leg bins per forest should retain the established loop.
            # A one-pixel pivot cap makes the fallback decision explicit.
            points=forests(nf=12,npix=40)
            extra={"lya3RMax":300,"lya3RBins":16,"lya3PivotCellMax":1}
            ref=run(points,"3pcf",kernel=1,extra=extra|{"search":"lya-3pcf-omp"})
            fast=run(points,"3pcf",extra=extra)
            compare(ref,fast,(3,))
            require(fast["metadata"]["lya_hierarchy"]["pixel_fallback"],"long-frontier fallback bypassed")
            require(fast["metadata"]["lya_hierarchy"]["nodes"]==0,"repeated forest groups must retain segment loop")
            if args.mpi_command: compare(fast,run(points,"3pcf",mpi=True,extra=extra),(3,))
            # A compact block includes pivots from many forests. Its cached
            # frontier must retain the first pivot's forest for later pivots.
            points=forests(nf=9,npix=16)
            points[:,:3]=100+.01*(points[:,:3]-100)
            ref=run(points,kernel=1,pair=0,extra={"search":"lya-2pcf-3pcf-omp"})
            fast=run(points,kernel=0,pair=0)
            compare(ref,fast,(2,3));compare(fast,run(points,kernel=0,pair=0,threads=4),(2,3),exact=True)
            require(fast["metadata"]["lya_los_tree"]["reuses"]>0,"LOS discovery reuse not exercised")
            if args.mpi_command: compare(fast,run(points,kernel=0,pair=0,mpi=True,threads=2),(2,3))
            # Near the expanded discovery boundary, actual pivot cuts must
            # decide inclusion. A block's first forest is a source for others.
            center=np.array([100.,0.,0.]); rows=[]
            for j in range(64): rows.append([*(center+np.array([0.,j*.005,0.])),.1+j*.01,1.,j%3])
            for j,d in enumerate([30.,30.-1e-5,30.+1e-5,29.9,30.1]):
                rows.append([*(center+np.array([d,0.,0.])),.2,0. if j==0 else 1.,j+10])
            points=np.array(rows)
            for stat in ("2pcf","3pcf","2pcf-3pcf"):
                extra={"lya3RMax":30,"lya2RpMax":30,"lya2RtMax":2}
                ref=run(points,stat,kernel=1,pair=0,extra=extra|{"search":f"lya-{stat}-omp"})
                fast=run(points,stat,kernel=0,pair=0,extra=extra)
                orders=(2,3) if stat=="2pcf-3pcf" else (int(stat[0]),)
                compare(ref,fast,orders)
                require(fast["metadata"]["lya_los_tree"]["reuses"]>0,"boundary cache bypassed")
                if args.mpi_command: compare(fast,run(points,stat,kernel=0,pair=0,mpi=True,extra=extra),orders)
            print("PASS LOS mixed-forest hierarchy, hash collisions, discovery reuse and ownership",flush=True)

        # Straight forests with dense radial patches activate certified range sums.
        points = forests(nf=8, npix=80)
        for fid in np.unique(points[:, 5]):
            idx = np.flatnonzero(points[:, 5] == fid)[:16]
            direction = points[idx[0], :3] / np.linalg.norm(points[idx[0], :3])
            points[idx, :3] = (100 + 0.0001 * np.arange(len(idx)))[:, None] * direction
        extra = {"lya2RpMax": 100, "lya2RtMax": 100, "lya2RpBins": 30, "lya2RtBins": 30}
        ref = run(points, "2pcf", pair=0, extra=extra)
        fast = run(points, "2pcf", extra=extra)
        compare(ref, fast, (2,))
        require(
            fast["metadata"]["lya_hierarchy"]["pair_ranges"] > 0,
            "range sums not exercised",
        )
        if args.mpi_command:
            compare(fast, run(points, "2pcf", mpi=True, extra=extra), (2,))
        print("PASS certified 2PCF range moments", flush=True)

        # Positive slop remains explicit and may move bins, never physical cuts.
        points = forests(nf=5, npix=15)
        ref = run(points, "3pcf", kernel=1, extra={"lya3RMax": 73})
        extra = dict(lya3RMax=73, lya3MuSlop=0.2, lya3RadialSlop=0.1, lya3PolarSlop=0.1)
        fast = run(points, "3pcf", extra=extra)
        require(ref["counts"] == fast["counts"], "slop changed domain")
        np.testing.assert_allclose(
            ref[3][:, -2:].sum(axis=0),
            fast[3][:, -2:].sum(axis=0),
            rtol=1e-11,
            atol=1e-10,
        )
        require(
            fast["metadata"]["lya_geometry"]["three_point_approximate"],
            "slop provenance",
        )
        if args.mpi_command:
            compare(fast, run(points, "3pcf", mpi=True, extra=extra), (3,))
        print("PASS explicit slop / strict physical cutoff / MPI agreement", flush=True)


if __name__ == "__main__":
    main()
