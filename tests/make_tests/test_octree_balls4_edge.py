#!/usr/bin/env python3
"""Independent masked/distinct-triplet oracle and C/OpenMP/MPI comparisons."""
import argparse
from pathlib import Path
import re
import shlex
import subprocess
import tempfile

import numpy as np

from test_two_ball_edge_corrections import (
    BINS, MMAX, RMIN, RMAX, catalog, polar, edge_solution, write_catalog, load_results,
)


def oracle(data, neighbor_data=None):
    p, k, w, mask = data
    same = neighbor_data is None
    q, qk, qw, qm = data if same else neighbor_data
    signal = np.zeros((MMAX+1, BINS, BINS), complex)
    window = np.zeros((2*MMAX+1, BINS, BINS), complex)
    for pivot in np.flatnonzero(mask):
        indices = np.flatnonzero(qm)
        if same:
            indices = indices[indices != pivot]
        indices = indices[np.linalg.norm(p[pivot] - q[indices], axis=1) > 0]
        distance, phase = polar(p[pivot], q[indices])
        keep = (distance > RMIN) & (distance < RMAX)
        indices, distance, phase = indices[keep], distance[keep], phase[keep]
        bins = np.floor(np.log(distance/RMIN)*BINS/np.log(RMAX/RMIN)).astype(int)
        rows, cols = np.meshgrid(bins, bins, indexing="ij")
        product = w[pivot]*np.outer(qw[indices], qw[indices])
        np.fill_diagonal(product, 0)
        field = k[pivot]*np.outer(qk[indices], qk[indices])
        for m in range(2*MMAX+1):
            moment = product*np.outer(phase**m, phase.conj()**m)
            np.add.at(window[m], (rows, cols), moment)
            if m <= MMAX:
                np.add.at(signal[m], (rows, cols), moment*field)
    return signal, window, edge_solution(signal, window)


def assert_histograms(actual, expected):
    for a, b in zip(actual, expected):
        assert np.all(np.isfinite(a))
        np.testing.assert_allclose(a, b, rtol=2e-9, atol=2e-8)


def cli_tests(executable, mpi_command=()):
    with tempfile.TemporaryDirectory(prefix="balls4-edge-regression-") as tmp:
        root = Path(tmp)
        data = catalog(count=160)
        data[3][9::3] = 0
        data[1][data[3] == 0] = 1e6
        data[2][data[3] == 0] = 1e5
        write_catalog(root/"input.txt", data)
        expected = oracle(data)
        serial = 0

        def run(engine, threads=1, *, edge=True, exact=True, output=None,
                failure=False, input_name="input.txt", cross=False,
                require_cells=False, theta=1.0):
            nonlocal serial
            serial += 1
            output = root/str(serial) if output is None else output
            output.mkdir(exist_ok=True)
            options = ["KKKCorrelation", "weights-norm", "read-mask",
                       "no-normalize-HistZeta", "compute-HistN", "and-CF"]
            if edge:
                options += ["edge-corrections"]
            if exact:
                options += ["no-one-ball", "no-two-balls"]
            command = [str(executable), f"search={engine}",
                       f"in={root/input_name}", "infmt=columns-ascii-all",
                       "useLogHist=true", "usePeriodic=false", f"rootDir={output}",
                       f"rangeN={RMAX}", f"rminHist={RMIN}", f"sizeHistN={BINS}",
                       f"mChebyshev={MMAX}", "sizeHistPhi=8", "nsmooth=2",
                       f"theta={theta}",
                       f"numberThreads={threads}", f"verbose={2 if require_cells else 0}",
                       "verbose_log=0",
                       f"options={','.join(options)}"]
            if cross:
                command = [s for s in command if not s.startswith(("in=", "infmt=", "options="))]
                command += [f"in={root/input_name},{root/'neighbor.txt'}",
                            "infmt=columns-ascii-all,columns-ascii-all", "iCatalogs=1,2",
                            "options=" + ",".join(s for s in options if s != "read-mask")]
            if engine.endswith("-mpi"):
                command = [*mpi_command, *command]
            proc = subprocess.run(command, text=True, capture_output=True, timeout=120)
            if failure:
                assert proc.returncode != 0, proc.stdout
                return proc.stdout+proc.stderr
            assert proc.returncode == 0, proc.stdout+proc.stderr
            if require_cells:
                accepted = re.findall(r"nbccalc\s*=\s*(\d+)", proc.stdout)
                assert accepted and max(map(int, accepted)) > 0, proc.stdout
            pair = tuple(np.loadtxt(output/f"{name}.txt") for name in
                         ("histNN", "histCF", "histXi2pcf"))
            if edge:
                return load_results(output), pair
            def read(component, m):
                return np.loadtxt(output/f"histZetaM_{component}_{m+1}.txt")
            matrices = tuple(read(c, m) for c in ("cos", "sin", "sincos", "cossin")
                             for m in range(MMAX+1))
            return matrices, pair

        base, pairs = run("octree-balls4-omp")
        assert_histograms(base, expected)
        for theta in (0.0, 3.0):
            bounded, _ = run("octree-balls4-omp", theta=theta)
            assert_histograms(bounded, expected)
        three, three_pairs = run("octree-balls4-omp", 3)
        for a, b in zip((*base, *pairs), (*three, *three_pairs)):
            np.testing.assert_array_equal(a, b)
        approximate, _ = run("octree-balls4-omp", 3, exact=False)
        np.testing.assert_allclose(approximate[2], edge_solution(*approximate[:2]),
                                   rtol=2e-9, atol=2e-8)
        # Normal (non-edge) masked traversal is preserved in the MPI addon.
        normal, normal_pairs = run("octree-balls4-omp", 3, edge=False)
        if mpi_command:
            for threads in (1, 3):
                distributed, distributed_pairs = run("octree-balls4-mpi", threads)
                for a, b in zip((*base, *pairs), (*distributed, *distributed_pairs)):
                    np.testing.assert_array_equal(a, b)
            distributed, distributed_pairs = run("octree-balls4-mpi", 3, edge=False)
            for a, b in zip((*normal, *normal_pairs), (*distributed, *distributed_pairs)):
                np.testing.assert_allclose(a, b, rtol=2e-8, atol=2e-8)
            other, _ = run("octree-balls4-mpi", 3, exact=False)
            for a, b in zip(approximate, other):
                np.testing.assert_array_equal(a, b)
        for engine in ("octree-balls4-omp", "octree-balls4-mpi") if mpi_command else ("octree-balls4-omp",):
            bad = root/("bad-"+engine)
            bad.mkdir()
            (bad/"histZetaM_window_Re_1.txt").mkdir()
            message = run(engine, output=bad, failure=True)
            assert "histZetaM_window_Re_1.txt" in message, message
            empty = tuple(x.copy() for x in data)
            empty[3][:] = 0
            empty[3][:2] = 1
            write_catalog(root/"sparse.txt", empty)
            sparse, _ = run(engine, input_name="sparse.txt")
            np.testing.assert_array_equal(sparse[2], 0)
            constant = tuple(x.copy() for x in data)
            constant[1][:] = 2.0
            write_catalog(root/"constant.txt", constant)
            const, _ = run(engine, input_name="constant.txt", exact=False)
            expected_const = oracle(constant)
            assert_histograms(const, expected_const)
        cross_data = tuple(x.copy() for x in data)
        cross_data[3][:] = 1
        cross_data[1][:] = 0.4 + cross_data[0][:, 0]
        cross_data[2][:] = 1.0
        neighbor_data = catalog(count=96)
        neighbor_data[0][:] = neighbor_data[0][:, [1, 2, 0]]
        write_catalog(root/"cross.txt", cross_data)
        write_catalog(root/"neighbor.txt", neighbor_data)
        cross_expected = oracle(cross_data, neighbor_data)
        clustered = catalog(count=136)
        centers = np.array(((0.3, 0.25, 0.4), (-0.3, -0.2, 0.4),
                            (0.3, -0.3, -0.2), (-0.3, 0.3, -0.4)))
        clustered[0][8:] = np.repeat(centers, 32, axis=0) + clustered[0][8:]*1.e-6
        write_catalog(root/"clusters.txt", clustered)
        cluster_expected = oracle(clustered)
        for engine in ("octree-balls4-omp", "octree-balls4-mpi") if mpi_command else ("octree-balls4-omp",):
            cross_result, _ = run(engine, input_name="cross.txt", cross=True)
            assert_histograms(cross_result, cross_expected)
            cluster_exact, _ = run(engine, input_name="clusters.txt")
            assert_histograms(cluster_exact, cluster_expected)
            clustered_result, _ = run(engine, input_name="clusters.txt", exact=False,
                                      require_cells=True, theta=0.001)
            for a, b in zip(clustered_result, cluster_expected):
                np.testing.assert_allclose(a, b, rtol=1.e-4, atol=2.e-5)
        print("PASS: BALLS4 mask, weighted complex edge oracle, C/OpenMP"
              + ("/MPI" if mpi_command else "")
              + ", sparse windows, output failures, cross catalogs, and cell aggregation")


def cython_tests():
    from cyballs import cballs, search_method_id
    from kappa_corr_all_engines import KappaCatalog, RunConfig, run_engine_suite
    data = catalog(count=160)
    data[3][9::3] = 0
    expected = oracle(data)[2]
    with tempfile.TemporaryDirectory(prefix="balls4-edge-cython-") as tmp:
        for engine in ("octree-balls4-omp", "octree-balls4-mpi"):
            assert search_method_id(engine) >= 0
            config = RunConfig(engines=(engine,), output_dir=Path(tmp)/engine,
                               theta_min=RMIN, theta_max=RMAX, theta_scale="au",
                               bins=BINS, multipoles=MMAX, threads=3,
                               edge_corrections=True, verbose=0, verbose_log=0,
                               options=("weights-norm", "no-one-ball", "no-two-balls", "no-out-Hist"))
            result = run_engine_suite(KappaCatalog(positions=data[0], kappa=data[1],
                                                   weights=data[2], mask=data[3]), config)[engine]
            assert not result["warnings"], result["warnings"]
            for m in range(MMAX+1):
                np.testing.assert_allclose(result[f"zeta_edge_complex_{m+1}"],
                                           expected[m], rtol=2e-9, atol=2e-8)
        balls = cballs()
        balls.set_catalog(data[0], kappa=data[1], weights=data[2], mask=data[3])
        params = dict(searchMethod="octree-balls4-omp", rangeN=RMAX, rminHist=RMIN,
                      sizeHistN=BINS, mChebyshev=MMAX, useLogHist=True,
                      rootDir=tmp, verbose=0, verbose_log=0, numberThreads=2,
                      options="KKKCorrelation,read-mask,no-normalize-HistZeta,no-out-Hist")
        try:
            for _ in range(3):
                balls.set(dict(params, options=params["options"]+",smooth-pivot,edge-corrections"))
                try:
                    balls.Run(level=["MainLoop"])
                except Exception as exc:
                    assert "smooth-pivot" in str(exc)
                else:
                    raise AssertionError("smooth-pivot accepted")
                balls.struct_cleanup()
            balls.set(params)
            balls.Run(level=["MainLoop"])
            try:
                balls.getHistZetaM_EE_complex(1)
            except Exception as exc:
                assert "enable edge-corrections" in str(exc)
            else:
                raise AssertionError("uncorrected histogram exposed as corrected")
        finally:
            balls.struct_cleanup()
            balls.clear_catalogs()
    print("PASS: BALLS4 Cython/driver edge getters, in-memory mask, and failure recovery")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cballs", type=Path)
    parser.add_argument("--mpi-command", default="")
    parser.add_argument("--cython", action="store_true")
    args = parser.parse_args()
    if args.cballs:
        cli_tests(args.cballs.resolve(), shlex.split(args.mpi_command))
    if args.cython:
        cython_tests()


if __name__ == "__main__":
    main()
