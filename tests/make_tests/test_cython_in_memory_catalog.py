#!/usr/bin/env python3

import tempfile
from pathlib import Path

import numpy as np

from cyballs import CosmoSevereError, cballs


def parameters(root_dir):
    return {
        "searchMethod": "octree-ggg-omp",
        "rangeN": 0.8,
        "rminHist": 0.05,
        "sizeHistN": 8,
        "mChebyshev": 3,
        "lengthBox": 1.0,
        "numberThreads": 2,
        "verbose": 0,
        "verbose_log": 0,
        "rootDir": root_dir,
        "options": "no-out-Hist",
    }


def write_ascii_catalog(path, positions, kappa):
    with path.open("w", encoding="ascii") as stream:
        stream.write("# nbody NDIM Lx Ly Lz\n")
        stream.write(f"# {len(positions)} 3 1.0 1.0 1.0\n")
        for point, field in zip(positions, kappa):
            stream.write(
                "{:.17g} {:.17g} {:.17g} {:.17g}\n".format(
                    point[0], point[1], point[2], field
                )
            )


def collect(balls):
    balls.Run(level=["MainLoop"])
    result = {
        "rBins": balls.getrBins().copy(),
        "histNN": balls.getHistNN().copy(),
        "histCF": balls.getHistCF().copy(),
        "histXi2pcf": balls.getHistXi2pcf().copy(),
    }
    for m in range(1, 5):
        for component in range(1, 5):
            result[f"zeta[{m},{component}]"] = (
                balls.getHistZetaMsincos(m, component).copy()
            )
    return result


def assert_same(reference, candidate):
    if reference.keys() != candidate.keys():
        raise AssertionError("file and memory runs returned different outputs")
    for name in reference:
        np.testing.assert_allclose(
            candidate[name], reference[name], rtol=1.0e-14, atol=1.0e-14,
            equal_nan=True, err_msg=name
        )


def test_in_memory_catalog_matches_file_reader_and_reloads_cleanly():
    rng = np.random.default_rng(1729)
    positions = rng.uniform(-0.4, 0.4, size=(32, 3))
    kappa = 1.0 + positions[:, 0] - 0.5 * positions[:, 1]
    original_positions = positions.copy()

    with tempfile.TemporaryDirectory(prefix="ctreeballs-memory-catalog-") as tmp:
        root = Path(tmp)
        catalog_path = root / "catalog.txt"
        write_ascii_catalog(catalog_path, positions, kappa)

        file_balls = cballs()
        file_parameters = parameters(str(root))
        file_parameters.update(
            {"infile": str(catalog_path), "infileformat": "columns-ascii"}
        )
        file_balls.set(file_parameters)
        reference = collect(file_balls)
        file_balls.struct_cleanup()

        memory_balls = cballs()
        memory_balls.set(parameters(str(root)))
        memory_balls.set_catalog(positions, kappa=kappa)
        if memory_balls.catalog_count != 1:
            raise AssertionError("catalog was not registered")
        first = collect(memory_balls)
        if memory_balls.getNBody() != len(positions):
            raise AssertionError("C did not publish the in-memory body count")
        if not memory_balls._runtime_bodytable_address():
            raise AssertionError("C did not publish ownership of the body table")
        if not np.array_equal(positions, original_positions):
            raise AssertionError("the C run modified the NumPy positions")
        assert_same(reference, first)

        memory_balls.struct_cleanup()
        if memory_balls._runtime_bodytable_address():
            raise AssertionError("cleanup retained the in-memory body table")
        second = collect(memory_balls)
        assert_same(first, second)
        memory_balls.struct_cleanup()


def test_in_memory_catalog_validation_is_recoverable():
    balls = cballs()
    try:
        balls.set_catalog(np.zeros((3, 2)))
    except CosmoSevereError:
        pass
    else:
        raise AssertionError("wrong-dimensional positions were accepted")

    bad = np.zeros((3, 3))
    bad[1, 0] = np.nan
    try:
        balls.set_catalog(bad)
    except CosmoSevereError:
        pass
    else:
        raise AssertionError("non-finite positions were accepted")


if __name__ == "__main__":
    test_in_memory_catalog_matches_file_reader_and_reloads_cleanly()
    test_in_memory_catalog_validation_is_recoverable()
    print("PASS: cyballs in-memory catalogs match file input and cleanly reload")
