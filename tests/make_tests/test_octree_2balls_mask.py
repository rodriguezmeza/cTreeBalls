#!/usr/bin/env python3
"""Masked native-octree 2PCF/3PCF: CLI, Cython, and MPI regressions."""

import argparse
import itertools
import os
from pathlib import Path
import re
import shlex
import subprocess
import tempfile

import numpy as np


def catalog(count=64):
    rng = np.random.default_rng(28461)
    positions = rng.uniform(-0.45, 0.45, (count, 3))
    # Keep bounding-box anchors valid, so filtered and full inputs use the
    # same reference frame after the normal startup recentering.
    positions[:8] = list(itertools.product((-0.7, 0.7), repeat=3))
    kappa = rng.uniform(-1.0, 2.0, count)
    weights = rng.uniform(0.3, 1.7, count)
    mask = np.ones(count, dtype=np.uint8)
    mask[8:] = (np.arange(count - 8) % 3 == 0)
    return positions, kappa, weights, mask


def write_catalog(path, data):
    positions, kappa, weights, mask = data
    values = np.column_stack((positions, kappa, weights, mask))
    np.savetxt(path, values, fmt="%.17g", comments="",
               header=f"# nbody NDIM Lx Ly Lz\n# {len(positions)} 3 2 2 2")


def parameters(root, threads=1, options=()):
    return {
        "rangeN": 2.5, "rminHist": 0.001, "sizeHistN": 4,
        "mChebyshev": 2, "lengthBox": 2, "theta": 1,
        "nsmooth": 4, "usePeriodic": False, "useLogHist": False,
        "numberThreads": threads, "verbose": 1, "verbose_log": 0,
        "stepState": 1000000, "rootDir": str(root),
        "options": ",".join(("compute-HistN", "out-m-HistZeta", *options)),
    }


def assert_same(reference, candidate, exact=False):
    assert reference.keys() == candidate.keys(), "histogram sets differ"
    for name in reference:
        assert np.all(np.isfinite(candidate[name])), name
        if exact:
            np.testing.assert_array_equal(reference[name], candidate[name],
                                          err_msg=name)
        else:
            np.testing.assert_allclose(reference[name], candidate[name],
                                       rtol=3e-12, atol=3e-13, err_msg=name)


def pair_oracle(data, weighted):
    positions, kappa, weights, mask = data
    selected = np.flatnonzero(mask)
    counts = np.zeros(4)
    numerator = np.zeros(4)
    denominator = np.zeros(4)
    for i, j in itertools.combinations(selected, 2):
        distance = np.linalg.norm(positions[i] - positions[j])
        if not 0.001 < distance < 2.5:
            continue
        radial = int((distance - 0.001) / ((2.5 - 0.001) / 4))
        weight = weights[i] * weights[j] if weighted else 1.0
        counts[radial] += 1
        numerator[radial] += weight * kappa[i] * kappa[j]
        denominator[radial] += weight
    return counts, np.divide(numerator, denominator,
                             out=np.zeros(4), where=denominator != 0)


def monopole_oracle(data, raw=False):
    positions, kappa, weights, mask = data
    selected = np.flatnonzero(mask)
    numerator = np.zeros((4, 4))
    denominator = np.zeros((4, 4))
    for pivot in selected:
        neighbors = []
        for neighbor in selected:
            distance = np.linalg.norm(positions[pivot] - positions[neighbor])
            if neighbor != pivot and 0.001 < distance < 2.5:
                radial = int((distance - 0.001) / ((2.5 - 0.001) / 4))
                neighbors.append((neighbor, radial))
        for (q, nq), (r, nr) in itertools.permutations(neighbors, 2):
            weight = weights[pivot] * weights[q] * weights[r]
            numerator[nq, nr] += weight * kappa[pivot] * kappa[q] * kappa[r]
            denominator[nq, nr] += weight
    return numerator if raw else np.divide(
        numerator, denominator, out=np.zeros((4, 4)), where=denominator != 0)


def cli_suite(executable, mpi_command=None):
    with tempfile.TemporaryDirectory(prefix="ctreeballs-2balls-mask-") as tmp:
        root = Path(tmp)
        data = catalog()
        selected = data[3].astype(bool)
        filtered = tuple(array[selected] for array in data)
        write_catalog(root / "masked.txt", data)
        write_catalog(root / "filtered.txt", filtered)
        poison = tuple(array.copy() for array in data)
        poison[1][~selected] = 1.e6
        poison[2][~selected] = 1.e5
        write_catalog(root / "poison.txt", poison)
        empty = tuple(array.copy() for array in data)
        empty[3][:] = 0
        write_catalog(root / "empty.txt", empty)
        serial = 0

        def run(input_name, options, threads=1, mpi=False, fails=False,
                companion=False, expected_mask=None, require_aggregation=False):
            nonlocal serial
            serial += 1
            out = root / str(serial)
            out.mkdir()
            method = "octree-2balls-mpi" if mpi else "octree-2balls-omp"
            params = parameters(out, threads, options)
            if companion:
                params["iCatalogs"] = "1,1"
            name_map = {"searchMethod": "search"}
            args = [str(executable), f"search={method}"]
            path = str(root / f"{input_name}.txt")
            args += [f"in={path},{path}" if companion else f"in={path}",
                     "infmt=columns-ascii-all,columns-ascii-all" if companion
                     else "infmt=columns-ascii-all"]
            args += [f"{name_map.get(key, key)}={str(value).lower() if isinstance(value, bool) else value}"
                     for key, value in params.items()]
            if mpi:
                args = [*mpi_command, *args]
            result = subprocess.run(args, text=True, capture_output=True,
                                    timeout=120)
            (out / "run.log").write_text(result.stdout + result.stderr)
            if fails:
                assert result.returncode != 0, "empty mask was accepted"
                assert "no unmasked bodies" in result.stdout + result.stderr, (
                    result.stdout + result.stderr)
                return {}
            assert result.returncode == 0, result.stdout + result.stderr
            if "read-mask" in options:
                validity = data[3] if expected_mask is None else expected_mask
                expected = f"mask keeps {int(validity.sum())} of {len(validity)} bodies"
                assert expected in result.stdout, result.stdout
                assert f"binary capacity={2 * int(validity.sum())} nodes" in result.stdout
            if require_aggregation:
                counts = re.findall(r"(?:nbccalc|accepted cell orientations) = (\d+)",
                                    result.stdout)
                assert any(int(count) > 0 for count in counts), result.stdout
            outputs = {p.name: np.loadtxt(p) for p in out.glob("hist*.txt")}
            assert outputs, "no histograms were written"
            return outputs

        for weighted in (False, True):
            for raw in (False, True):
                options = ["no-two-balls"]
                if weighted:
                    options.append("weights-norm")
                if raw:
                    options.append("no-normalize-HistZeta")
                reference = run("filtered", options)
                masked = run("masked", [*options, "read-mask"])
                assert_same(reference, masked)
                assert_same(masked, run("poison", [*options, "read-mask"]), True)
                assert_same(masked, run("masked", [*options, "read-mask"], 4), True)
                direct = run("filtered", [*options, "treecorr-direct-triples"])
                assert_same(direct, masked)
                counts, correlation = pair_oracle(data, weighted)
                np.testing.assert_allclose(masked["histNN.txt"][:, -1], counts)
                # The CLI writer publishes Xi with %16.8e precision.
                np.testing.assert_allclose(masked["histXi2pcf.txt"][:, -1],
                                           correlation, rtol=5e-9, atol=5e-12)
                if mpi_command:
                    distributed = run("masked", [*options, "read-mask"], 2, True)
                    assert_same(masked, distributed)
                    assert_same(distributed, run("masked", [*options, "read-mask"],
                                                 2, True), True)

        large = catalog(1024)
        large_poison = tuple(array.copy() for array in large)
        large_poison[1][large[3] == 0] = 1.e6
        large_poison[2][large[3] == 0] = 1.e5
        write_catalog(root / "large.txt", large)
        write_catalog(root / "large-poison.txt", large_poison)
        for mode in ("only-2pcf", "only-3pcf"):
            options = [mode, "weights-norm", "read-mask"]
            reference = run("large", options, expected_mask=large[3],
                            require_aggregation=True)
            assert_same(reference, run("large-poison", options, 4,
                                        expected_mask=large[3]), True)
            if mpi_command:
                assert_same(reference, run("large", options, 2, True,
                                            expected_mask=large[3]))

        options = ["only-2pcf", "no-two-balls", "weights-norm", "and-CF"]
        reference = run("filtered", options)
        masked = run("masked", [*options, "read-mask"])
        assert_same(reference, masked)
        assert_same(masked, run("masked", [*options, "read-mask"], companion=True))
        run("empty", ["read-mask"], fails=True)
        if mpi_command:
            assert_same(masked, run("masked", [*options, "read-mask"], 2, True))
            run("empty", ["read-mask"], 2, True, fails=True)
    print("PASS: masked octree two-ball CLI pairs, multipoles, normalization, and determinism")


def cython_suite():
    from cyballs import cballs

    data = catalog()
    selected = data[3].astype(bool)
    with tempfile.TemporaryDirectory(prefix="ctreeballs-2balls-mask-cython-") as tmp:
        balls = cballs()

        def run(values, masked, threads=1, raw=False):
            options = ["no-two-balls", "weights-norm", "no-out-Hist"]
            if raw:
                options.append("no-normalize-HistZeta")
            if masked:
                options.append("read-mask")
            balls.set(dict(parameters(tmp, threads, options),
                           searchMethod="octree-2balls-omp", verbose=0))
            balls.set_catalog(values[0], kappa=values[1], weights=values[2],
                              mask=values[3])
            try:
                balls.Run(level=["MainLoop"])
                result = {"NN": balls.getHistNN().copy(),
                          "KK": balls.getHistXi2pcf().copy()}
                for order in range(1, 4):
                    for component in range(1, 5):
                        result[f"KKK_{order}_{component}"] = (
                            balls.getHistZetaMsincos(order, component).copy())
                return result
            finally:
                balls.struct_cleanup()

        reference = run(tuple(array[selected] for array in data), False)
        masked = run(data, True)
        assert_same(reference, masked)
        np.testing.assert_allclose(masked["KKK_1_1"], monopole_oracle(data),
                                   rtol=2e-12, atol=2e-13)
        np.testing.assert_allclose(run(data, True, raw=True)["KKK_1_1"],
                                   monopole_oracle(data, raw=True),
                                   rtol=2e-12, atol=2e-13)
        assert_same(masked, run(data, True, 4), True)
        all_valid = tuple(array.copy() for array in data)
        all_valid[3][:] = 1
        assert_same(run(all_valid, True), run(data, False), True)
        for live in (1, 2):
            sparse = tuple(array.copy() for array in data)
            sparse[3][:] = 0
            sparse[3][:live] = 1
            for raw in (False, True):
                result = run(sparse, True, raw=raw)
                for name, values in result.items():
                    assert np.all(np.isfinite(values)), name
                    if name.startswith("KKK") or (live == 1 and name == "NN"):
                        assert not np.any(values), name
        empty = tuple(array.copy() for array in data)
        empty[3][:] = 0
        for _ in range(3):
            try:
                run(empty, True)
            except Exception as error:
                assert "no unmasked bodies" in str(error), str(error)
            else:
                raise AssertionError("empty mask was accepted")
            assert_same(masked, run(data, True), True)
    print("PASS: masked octree two-ball Cython input, sparse masks, and failure recovery")


def fits_suite(executable, mpi_command=None):
    import healpy as hp

    nside = 4
    count = hp.nside2npix(nside)
    positions = np.array(hp.pix2vec(nside, np.arange(count))).T
    kappa = (0.3 + positions[:, 0] - 0.2 * positions[:, 1]).astype(np.float32)
    mask = ((positions[:, 2] > -0.3) & (positions[:, 0] > -0.4)).astype(np.uint8)
    with tempfile.TemporaryDirectory(prefix="ctreeballs-2balls-mask-fits-") as tmp:
        root = Path(tmp)
        hp.write_map(root / "data.fits", kappa, dtype=np.float32, nest=False)
        hp.write_map(root / "mask.fits", mask.astype(np.float32),
                     dtype=np.float32, nest=False)
        write_catalog(root / "inline.txt", (positions, kappa, np.ones(count), mask))
        results = []
        cases = [(False, False), (True, False)]
        if mpi_command:
            cases.append((True, True))
        for index, (fits, mpi) in enumerate(cases):
            out = root / str(index)
            out.mkdir()
            method = "octree-2balls-mpi" if mpi else "octree-2balls-omp"
            params = parameters(out, 2, ["read-mask", "no-two-balls", "weights-norm"])
            # HEALPix antipodal directions have an ambiguous angular phase;
            # use the local angular range intended for these sky estimators.
            params["rangeN"] = 0.713
            if fits:
                params.update(infile=f"{root}/data.fits,{root}/mask.fits",
                              infileformat="fits-healpix,fits-healpix", iCatalogs="1,1")
            else:
                params.update(infile=str(root / "inline.txt"),
                              infileformat="columns-ascii-all")
            aliases = {"infile": "in", "infileformat": "infmt"}
            command = [str(executable), f"search={method}", *[
                f"{aliases.get(key, key)}={str(value).lower() if isinstance(value, bool) else value}"
                for key, value in params.items()]]
            result = subprocess.run([*(mpi_command if mpi else []), *command],
                                    text=True, capture_output=True, timeout=120)
            assert result.returncode == 0, result.stdout + result.stderr
            results.append({p.name: np.loadtxt(p) for p in out.glob("hist*.txt")})
        for candidate in results[1:]:
            assert_same(results[0], candidate)
    print("PASS: HEALPix companion masks match embedded masks in octree two-ball searches")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cballs", type=Path)
    parser.add_argument("--mpi", action="store_true")
    parser.add_argument("--cython", action="store_true")
    parser.add_argument("--fits", action="store_true",
                        help="also exercise CFITSIO HEALPix input (requires healpy)")
    args = parser.parse_args()
    if args.fits and not args.cballs:
        parser.error("--fits requires --cballs")
    if args.cballs:
        mpi = None
        if args.mpi:
            launcher = os.environ.get("MPIEXEC", "mpiexec")
            extra = shlex.split(os.environ.get("MPIEXEC_ARGS", ""))
            if not extra:
                version = subprocess.run([launcher, "--version"], text=True,
                                         capture_output=True, timeout=10)
                if "open mpi" in (version.stdout + version.stderr).lower():
                    extra = ["--oversubscribe", "--bind-to", "none"]
            mpi = [launcher, *extra, "-n", os.environ.get("MPI_RANKS", "2")]
        cli_suite(args.cballs.resolve(), mpi)
        if args.fits:
            fits_suite(args.cballs.resolve(), mpi)
    if args.cython:
        cython_suite()
    if not args.cballs and not args.cython:
        parser.error("select --cballs and/or --cython")
