#!/usr/bin/env python3
"""Independent raw scalar oracle, rotations, inactive smoothing, and MPI checks."""
import os
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
MPI_MODE = os.environ.get("CBALLS_TEST_MPI") == "1"
if MPI_MODE:
    from mpi4py import MPI
    COMM = MPI.COMM_WORLD
else:
    COMM = None

from cyballs import cballs, search_method_id

METHODS = (
    "octree-ggg-omp", "kdtree-omp", "balltree-omp", "octree-balls4-omp",
    "balltree-2balls-omp", "octree-2balls-omp", "balltree-2balls-omp_3pcf",
)
EDGES = np.geomspace(2*np.sin(np.pi/360), np.sqrt(3), 5)


def catalog(patch=False):
    rng = np.random.default_rng(9428)
    positions = rng.normal(size=(36, 3))
    if patch:
        positions[:, 2] = np.abs(positions[:, 2]) + 0.4
    positions /= np.linalg.norm(positions, axis=1)[:, None]
    return positions, rng.normal(size=36), rng.uniform(0.5, 1.8, 36)


def oracle(positions, field, weights):
    """Enumerate ordered, distinct legs. No tree, phase basis, or moment products."""
    pairs = np.zeros(4)
    result = np.zeros((4, 4, 4), dtype=complex)
    wk = weights*field
    for i, pivot in enumerate(positions):
        legs = positions-pivot
        distances = np.linalg.norm(legs, axis=1)
        bins = np.searchsorted(EDGES, distances, side="right")-1
        valid = np.flatnonzero((distances > EDGES[0]) & (distances < EDGES[-1]))
        np.add.at(pairs, bins[valid], 0.5)
        radius = np.linalg.norm(pivot)
        if radius == 0:
            continue
        normal = pivot/radius
        tangent = legs - np.outer(legs @ normal, normal)
        lengths = np.linalg.norm(tangent, axis=1)
        for j in valid:
            for k in valid:
                if j == k or min(lengths[j], lengths[k]) < 1e-14:
                    continue
                u, v = tangent[j]/lengths[j], tangent[k]/lengths[k]
                phase = np.arctan2(normal @ np.cross(u, v), u @ v)
                result[:, bins[j], bins[k]] += wk[i]*wk[j]*wk[k]*np.exp(
                    1j*np.arange(4)*phase)
    return pairs, result


def run(method, positions, field, weights, threads=1, rsmooth=None,
        exact=True, theta=1.0):
    root = tempfile.mkdtemp(prefix="scalar-contract-") if COMM is None or COMM.rank == 0 else None
    if COMM is not None:
        root = COMM.bcast(root, root=0)
    options = ["compute-HistN", "no-out-Hist", "KKKCorrelation",
               "weights-norm", "no-normalize-HistZeta"]
    if exact:
        options += ["no-one-ball", "no-two-balls"]
    params = dict(searchMethod=method, rangeN=EDGES[-1], rminHist=EDGES[0],
                  sizeHistN=4, mChebyshev=3, lengthBox=2.0, numberThreads=threads,
                  useLogHist=True, usePeriodic=False, theta=theta,
                  verbose=0, verbose_log=0, rootDir=root, options=",".join(options))
    if rsmooth is not None:
        params["rsmooth"] = str(rsmooth)
    model = cballs()
    try:
        model.set(params)
        model.set_catalog(positions, kappa=field, weights=weights)
        model.Run(level=["MainLoop"])
        if COMM is not None and COMM.rank != 0:
            return None
        modes = []
        for m in range(1, 5):
            cc, ss, sc, cs = [model.getHistZetaMsincos(m, c).copy() for c in range(1, 5)]
            modes.append(cc+ss+1j*(sc-cs))
        pairs = None if method.endswith("_3pcf") else model.getHistNN().copy()
        return pairs, np.asarray(modes)
    finally:
        model.struct_cleanup()
        if COMM is not None:
            COMM.Barrier()
        if COMM is None or COMM.rank == 0:
            import shutil
            shutil.rmtree(root)


def check_close(actual, expected, label):
    error = None
    if COMM is None or COMM.rank == 0:
        try:
            if actual[0] is not None:
                np.testing.assert_array_equal(actual[0], expected[0], err_msg=label+" pairs")
            np.testing.assert_allclose(actual[1], expected[1], rtol=2e-11, atol=2e-10, err_msg=label)
        except AssertionError as exc:
            error = str(exc)
    if COMM is not None:
        error = COMM.bcast(error, root=0)
    if error:
        raise AssertionError(error)


def test_scalar_contract():
    methods = [m.replace("-omp", "-mpi") for m in METHODS if m != "kdtree-omp"] if MPI_MODE else METHODS
    rotation, _ = np.linalg.qr(np.random.default_rng(23).normal(size=(3, 3)))
    rotation[:, 0] *= np.linalg.det(rotation)
    checked = 0
    for patch in (False, True):
        positions, field, weights = catalog(patch)
        expected = oracle(positions, field, weights)
        for method in methods:
            if search_method_id(method) < 0:
                continue
            check_close(run(method, positions, field, weights), expected, method)
            check_close(run(method, positions @ rotation, field, weights, threads=3),
                        expected, method+" rotated/threads=3")
            check_close(run(method, positions, field, weights, rsmooth=0), expected,
                        method+" explicit zero smoothing")
            check_close(run(method, positions, field, weights, rsmooth=20), expected,
                        method+" radius without opt-in")
            # Aggregation must converge to the same estimator at small theta.
            check_close(run(method, positions, field, weights, exact=False, theta=1e-6),
                        expected, method+" converged aggregation")
            checked += 1
    assert checked, "No scalar engines were tested"


def test_undefined_bearings():
    positions = np.array([[0, 0, 0], [0, 0, 0.5], [0, 0, 1], [0, 0, -1],
                          [0.6, 0, 0.8], [0, 0.8, 0.6], [-0.8, 0.6, 0]])
    field = np.linspace(-0.5, 1.0, len(positions))
    weights = np.linspace(0.6, 1.4, len(positions))
    expected = oracle(positions, field, weights)
    methods = ([m.replace("-omp", "-mpi") for m in METHODS if m != "kdtree-omp"]
               if MPI_MODE else METHODS)
    for method in methods:
        if search_method_id(method) >= 0:
            check_close(run(method, positions, field, weights), expected,
                        method+" undefined bearings retain ordinary pairs")


def test_treecorr_reference():
    if MPI_MODE:
        return
    try:
        import treecorr
    except ImportError:
        import pytest
        pytest.skip("TreeCorr is not installed")
    positions, field, weights = catalog(True)
    expected = oracle(positions, field, weights)[1]
    corr = treecorr.KKKCorrelation(min_sep=EDGES[0], max_sep=EDGES[-1], nbins=4,
        max_n=3, bin_type="LogMultipole", bin_slop=0, angle_slop=0, brute=True)
    corr.process(treecorr.Catalog(x=positions[:, 0], y=positions[:, 1], z=positions[:, 2],
                                 k=field, w=weights), num_threads=1)
    np.testing.assert_allclose(np.moveaxis(corr.zeta[:, :, 3:], -1, 0), expected,
                               rtol=2e-6, atol=2e-6)


def test_c_executable():
    if MPI_MODE:
        return
    positions, field, _ = catalog(True)
    expected = oracle(positions, field, np.ones(len(field)))
    with tempfile.TemporaryDirectory(prefix="scalar-contract-c-") as directory:
        root = Path(directory)
        for fmt in ("columns-ascii", "multi-columns-ascii"):
            path = root/f"{fmt}.txt"
            header = (f"nbody NDIM Lx Ly Lz\n{len(field)} 3 2 2 2"
                      if fmt == "columns-ascii" else "")
            np.savetxt(path, np.column_stack((positions, field)), header=header)
            for method in METHODS:
                if search_method_id(method) < 0:
                    continue
                output = root/f"{fmt}-{method}"
                output.mkdir()
                subprocess.run([str(ROOT/"cballs"), f"search={method}", f"in={path}",
                    f"infmt={fmt}", "columns=1,2,3,4", f"rootDir={output}",
                    f"rangeN={EDGES[-1]}", f"rminHist={EDGES[0]}", "sizeHistN=4",
                    "mChebyshev=3", "lengthBox=2", "nthreads=2", "useLogHist=true",
                    "usePeriodic=false", "verb=0", "verblog=0",
                    "options=pos-and-convergence,compute-HistN,no-normalize-HistZeta,"
                    "no-one-ball,no-two-balls,KKKCorrelation,weights-norm,out-m-HistZeta"],
                    check=True, capture_output=True, text=True, timeout=60)
                if not method.endswith("_3pcf"):
                    np.testing.assert_array_equal(
                        np.loadtxt(output/"histNN.txt")[:, 1], expected[0])
                modes = []
                for m in range(1, 5):
                    cc, ss, sc, cs = [np.loadtxt(output/f"histZetaM_{name}_{m}.txt")
                                      for name in ("cos", "sin", "sincos", "cossin")]
                    modes.append(cc+ss+1j*(sc-cs))
                np.testing.assert_allclose(modes, expected[1], rtol=2e-8, atol=2e-7,
                                           err_msg=f"{method}: {fmt}")


if __name__ == "__main__":
    test_scalar_contract()
    test_undefined_bearings()
    if not MPI_MODE:
        test_c_executable()
    if COMM is None or COMM.rank == 0:
        print("PASS: scalar oracle, rotations, inactive smoothing, and parallel consistency")
