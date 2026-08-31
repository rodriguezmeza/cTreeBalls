#!/usr/bin/env python3
"""Independent complex-window oracle for the six TreeCorr-style addons."""

import argparse
import itertools
from pathlib import Path
import shlex
import subprocess
import tempfile

import numpy as np


ENGINES = (
    "balltree-2balls-omp", "balltree-2balls-mpi",
    "balltree-2balls-omp_3pcf", "balltree-2balls-mpi_3pcf",
    "octree-2balls-omp", "octree-2balls-mpi",
)
MMAX, BINS, RMIN, RMAX = 2, 4, 0.02, 1.5


def catalog(dimension=3, count=64):
    rng = np.random.default_rng(81142)
    positions = rng.uniform(-0.5, 0.5, (count, dimension))
    anchors = list(itertools.product((-0.6, 0.6), repeat=dimension))
    positions[:len(anchors)] = anchors
    kappa = 0.5 + positions[:, 0] - 0.3 * positions[:, 1]
    weights = rng.uniform(0.4, 1.7, count)
    mask = np.ones(count, dtype=np.uint8)
    return positions, kappa, weights, mask


def write_catalog(path, data):
    positions, kappa, weights, mask = data
    dimension = positions.shape[1]
    np.savetxt(path, np.column_stack((positions, kappa, weights, mask)),
               fmt="%.17g", comments="",
               header=f"# nbody NDIM box\n# {len(positions)} {dimension} "
                      + " ".join(["2"] * dimension))


def polar(pivot, neighbors):
    dr = pivot - neighbors
    distance = np.linalg.norm(dr, axis=1)
    if len(pivot) == 2:
        return distance, (-dr[:, 0] - 1j * dr[:, 1]) / distance
    normal = pivot / np.linalg.norm(pivot)
    reference = np.array([0.0, 0.0, 1.0])
    if abs(normal[2]) > 0.99:
        reference = np.array([1.0, 0.0, 0.0])
    north = reference - normal * (reference @ normal)
    north /= np.linalg.norm(north)
    side = np.cross(north, normal)
    phase = (-dr @ north) + 1j * (-dr @ side)
    length = np.abs(phase)
    return distance, np.divide(phase, length, out=np.full_like(phase, np.nan),
                               where=length > 1e-14 * distance)


def brute_force(data, neighbor_data=None):
    positions, kappa, weights, mask = data
    selected = np.flatnonzero(mask)
    same_catalog = neighbor_data is None
    qpos, qkappa, qweights, qmask = data if same_catalog else neighbor_data
    qselected = np.flatnonzero(qmask)
    signal = np.zeros((MMAX + 1, BINS, BINS), dtype=complex)
    window = np.zeros((2 * MMAX + 1, BINS, BINS), dtype=complex)
    for pivot in selected:
        neighbors = qselected[qselected != pivot] if same_catalog else qselected
        neighbors = neighbors[np.linalg.norm(positions[pivot] - qpos[neighbors], axis=1) > 0]
        distance, phase = polar(positions[pivot], qpos[neighbors])
        valid = (distance > RMIN) & (distance < RMAX) & np.isfinite(phase)
        neighbors, distance, phase = neighbors[valid], distance[valid], phase[valid]
        bins = ((distance - RMIN) * BINS / (RMAX - RMIN)).astype(int)
        rows, cols = np.meshgrid(bins, bins, indexing="ij")
        weight = weights[pivot] * np.outer(qweights[neighbors], qweights[neighbors])
        np.fill_diagonal(weight, 0)
        field = kappa[pivot] * np.outer(qkappa[neighbors], qkappa[neighbors])
        for order in range(2 * MMAX + 1):
            moment = weight * np.outer(phase**order, (phase.conj())**order)
            np.add.at(window[order], (rows, cols), moment)
            if order <= MMAX:
                np.add.at(signal[order], (rows, cols), moment * field)
    return signal, window


def edge_solution(signal, window):
    result = np.zeros_like(signal)
    modes = np.arange(-MMAX, MMAX + 1)
    for i in range(BINS):
        for j in range(BINS):
            if window[0, i, j].real <= 0:
                continue
            difference = modes[:, None] - modes[None, :]
            matrix = window[np.abs(difference), i, j]
            matrix = np.where(difference < 0, matrix.conj(), matrix)
            rhs = signal[np.abs(modes), i, j]
            rhs = np.where(modes < 0, rhs.conj(), rhs)
            if np.linalg.cond(matrix) < 1.e12:
                result[:, i, j] = np.linalg.solve(matrix, rhs)[MMAX:]
    return result


def load_results(root):
    def read(suffix, order):
        return np.loadtxt(root / f"histZetaM_{suffix}_{order + 1}.txt")
    signal = np.array([read("cos", m) + read("sin", m)
                       + 1j * (read("sincos", m) - read("cossin", m))
                       for m in range(MMAX + 1)])
    window = np.array([read("window_Re", m) + 1j * read("window_Im", m)
                       for m in range(2 * MMAX + 1)])
    edge = np.array([read("EE", m) + 1j * read("EE_Im", m)
                     for m in range(MMAX + 1)])
    return signal, window, edge


def assert_results(actual, expected, raw_tolerance=2.e-8):
    for index, (a, b) in enumerate(zip(actual, expected)):
        assert np.all(np.isfinite(a))
        np.testing.assert_allclose(a, b, rtol=raw_tolerance if index == 0 else 2e-10,
                                   atol=raw_tolerance if index == 0 else 2e-10)


def cli_tests(executable, engines, dimension, mpi_command):
    with tempfile.TemporaryDirectory(prefix="ctreeballs-two-ball-edge-") as tmp:
        root = Path(tmp)
        data = catalog(dimension)
        signal, window = brute_force(data)
        reference = signal, window, edge_solution(signal, window)
        write_catalog(root / "input.txt", data)
        serial = 0

        def run(engine, datafile="input.txt", threads=1, extra=(), exact=True,
                expect_failure=False, output=None, weighted=True):
            nonlocal serial
            serial += 1
            output = root / str(serial) if output is None else output
            output.mkdir(exist_ok=True)
            options = ["KKKCorrelation", "only-3pcf",
                       "edge-corrections", "no-normalize-HistZeta", *extra]
            if weighted:
                options.append("weights-norm")
            if exact:
                options.append("no-two-balls")
            files = datafile if isinstance(datafile, tuple) else (datafile,)
            command = [str(executable), f"search={engine}",
                       "in=" + ",".join(str(root / name) for name in files),
                       "infmt=" + ",".join(["columns-ascii-all"] * len(files)),
                       "iCatalogs=1,2" if len(files) == 2 else "iCatalogs=1",
                       f"rootDir={output}",
                       f"numberThreads={threads}", "useLogHist=false",
                       "usePeriodic=false", f"rangeN={RMAX}", f"rminHist={RMIN}",
                       f"sizeHistN={BINS}", f"mChebyshev={MMAX}",
                       "sizeHistPhi=8", "nsmooth=2", "theta=1", "verbose=2",
                       "verbose_log=0", f"options={','.join(options)}"]
            if "-mpi" in engine and mpi_command:
                command = [*mpi_command, *command]
            proc = subprocess.run(command, text=True, capture_output=True, timeout=120)
            if expect_failure:
                assert proc.returncode != 0, proc.stdout
                return proc.stdout + proc.stderr
            assert proc.returncode == 0, proc.stdout + proc.stderr
            assert "edge correction uses window modes 0..4" in proc.stdout, proc.stdout
            return load_results(output)

        for engine in engines:
            one = run(engine)
            assert_results(one, reference)
            unweighted = (data[0], data[1], np.ones(len(data[0])), data[3])
            su, wu = brute_force(unweighted)
            assert_results(run(engine, weighted=False), (su, wu, edge_solution(su, wu)))
            neighbor_data = tuple(array.copy() for array in data)
            neighbor_data[0][:, :2] = data[0][:, [1, 0]] * np.array([-1, 1])
            write_catalog(root / "neighbor.txt", neighbor_data)
            sc, wc = brute_force(data, neighbor_data)
            assert_results(run(engine, ("input.txt", "neighbor.txt")),
                           (sc, wc, edge_solution(sc, wc)))
            many = run(engine, threads=4)
            for a, b in zip(one, many):
                np.testing.assert_array_equal(a, b)
            aggregate = run(engine, exact=False)
            repeated = run(engine, exact=False, threads=4)
            for a, b in zip(aggregate, repeated):
                np.testing.assert_array_equal(a, b)
            # Test the solver with the actual approximated signal/window too.
            np.testing.assert_allclose(aggregate[2], edge_solution(*aggregate[:2]),
                                       rtol=2e-7, atol=2e-8)
            if "-mpi" not in engine and engine != "balltree-2balls-omp":
                assert_results(run(engine, extra=("treecorr-direct-triples",)), reference)
            if engine.startswith("octree"):
                masked = tuple(array.copy() for array in data)
                masked[3][2**dimension::2] = 0
                masked[1][masked[3] == 0] = 1e6
                masked[2][masked[3] == 0] = 1e5
                write_catalog(root / "mask.txt", masked)
                sm, wm = brute_force(masked)
                assert_results(run(engine, "mask.txt", extra=("read-mask",)),
                               (sm, wm, edge_solution(sm, wm)))
            constant = tuple(array.copy() for array in data)
            constant[1][:] = 2.0
            write_catalog(root / "constant.txt", constant)
            _, constant_window, corrected = run(engine, "constant.txt", exact=False)
            modes = np.arange(-MMAX, MMAX + 1)
            difference = modes[:, None] - modes[None, :]
            for i in range(BINS):
                for j in range(BINS):
                    matrix = constant_window[np.abs(difference), i, j]
                    matrix = np.where(difference < 0, matrix.conj(), matrix)
                    # Narrow outer-bin angles can make an invertible window
                    # ill-conditioned; allow its predicted roundoff amplification.
                    tolerance = 1e-10 + 128 * np.finfo(float).eps * 8 * np.linalg.cond(matrix)
                    target = np.zeros(MMAX + 1)
                    target[0] = 8
                    np.testing.assert_allclose(corrected[:, i, j], target,
                                               rtol=0, atol=tolerance)
            bad_output = root / f"bad-{serial}"
            bad_output.mkdir()
            (bad_output / "histZetaM_window_Re_1.txt").mkdir()
            message = run(engine, output=bad_output, expect_failure=True)
            assert "histZetaM_window_Re_1.txt" in message, message
            print(f"PASS: {engine}: raw/window/edge oracle, determinism, failures", flush=True)


def treecorr_tests(executable, engines, mpi_command):
    import treecorr

    data = catalog(dimension=2)
    positions, kappa, weights, _ = data
    # TreeCorr stores leaf w and w*k as float. Dyadic fields and weights isolate
    # the estimator comparison from that different input-storage precision.
    kappa[:] = np.rint(8 * kappa) / 8
    weights[:] = np.rint(8 * weights) / 8
    reference = treecorr.KKKCorrelation(
        min_sep=RMIN, max_sep=RMAX, nbins=BINS, max_n=2 * MMAX,
        bin_type="LogMultipole", brute=True, bin_slop=0, angle_slop=0, verbose=0)
    reference.process(treecorr.Catalog(x=positions[:, 0], y=positions[:, 1],
                                      k=kappa, w=weights), num_threads=1)
    center = 2 * MMAX
    signal = np.moveaxis(reference.zeta[:, :, center:center + MMAX + 1], -1, 0)
    window = np.moveaxis(reference.weight[:, :, center:], -1, 0)
    expected = signal, window, edge_solution(signal, window)
    with tempfile.TemporaryDirectory(prefix="ctreeballs-two-ball-treecorr-") as tmp:
        root = Path(tmp)
        write_catalog(root / "input.txt", data)
        for engine in engines:
            output = root / engine
            output.mkdir()
            command = [str(executable), f"search={engine}",
                       f"in={root / 'input.txt'}", "infmt=columns-ascii-all",
                       "iCatalogs=1", f"rootDir={output}", "numberThreads=2",
                       "useLogHist=true", "usePeriodic=false", f"rangeN={RMAX}",
                       f"rminHist={RMIN}", f"sizeHistN={BINS}", f"mChebyshev={MMAX}",
                       "sizeHistPhi=8", "nsmooth=2", "verbose=0", "verbose_log=0",
                       "options=KKKCorrelation,weights-norm,only-3pcf,"
                       "edge-corrections,no-normalize-HistZeta,no-two-balls"]
            if "-mpi" in engine and mpi_command:
                command = [*mpi_command, *command]
            proc = subprocess.run(command, text=True, capture_output=True, timeout=120)
            assert proc.returncode == 0, proc.stdout + proc.stderr
            actual = load_results(output)
            assert_results(actual, expected)
            print(f"PASS: {engine}: TreeCorr {treecorr.__version__} raw/window/edge", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cballs", type=Path, required=True)
    parser.add_argument("--engine", action="append", choices=ENGINES)
    parser.add_argument("--dimension", type=int, choices=(2, 3), default=3)
    parser.add_argument("--mpi-command", default="")
    parser.add_argument("--treecorr", action="store_true", help="compare a 2D build to TreeCorr")
    args = parser.parse_args()
    executable = args.cballs.resolve()
    available = subprocess.run([str(executable), "options=print-search-methods"],
                               text=True, capture_output=True, check=False).stdout
    engines = args.engine or [engine for engine in ENGINES if f"- {engine} (id=" in available]
    if not engines:
        raise SystemExit("no two-ball engines selected; pass --engine explicitly")
    if args.treecorr:
        if args.dimension != 2:
            parser.error("--treecorr requires --dimension 2 and a 2D executable")
        treecorr_tests(executable, engines, shlex.split(args.mpi_command))
    else:
        cli_tests(executable, engines, args.dimension, shlex.split(args.mpi_command))


if __name__ == "__main__":
    main()
