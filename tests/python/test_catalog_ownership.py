"""No-copy registration ownership, recovery, and numerical equivalence."""
import gc
import inspect
import json
import weakref

import numpy as np
import pytest

import cyballs
from cyballs import cballs, CosmoSevereError


def arrays(count=48):
    rng = np.random.default_rng(791)
    positions = rng.normal(size=(count, 3))
    positions /= np.linalg.norm(positions, axis=1)[:, None]
    return dict(positions=positions, kappa=1 + .2 * positions[:, 0],
                weights=rng.uniform(.8, 1.2, count),
                mask=np.array(np.arange(count) % 7 != 0, dtype=np.uint8),
                gamma1=rng.normal(0, .02, count), gamma2=rng.normal(0, .02, count))


def configured(tmp_path, engine="octree-2balls-omp"):
    model = cballs()
    model.set(searchMethod=engine, rootDir=str(tmp_path), rminHist=.02, rangeN=1.8,
              sizeHistN=4, mChebyshev=2, numberThreads=1, theta=0., lengthBox=2.,
              useLogHist=False, verbose=0, verbose_log=0,
              options="no-out-Hist,no-smooth-pivot,no-one-ball,no-two-balls,read-mask,weights-norm")
    return model


def test_transfer_retains_original_arrays_and_release_drops_references():
    values = arrays()
    refs = {key: weakref.ref(value) for key, value in values.items()}
    model = cballs()
    model.set_catalog(**values, take_ownership=True)
    assert all(not value.flags.writeable for value in values.values())
    with pytest.raises(ValueError, match="read-only"):
        values["positions"][0, 0] = 7.
    del values
    gc.collect()
    assert all(ref() is not None for ref in refs.values())
    model.struct_cleanup()
    assert all(ref() is not None for ref in refs.values())
    model.clear_catalogs()
    gc.collect()
    assert all(ref() is None for ref in refs.values())


@pytest.mark.parametrize("release", ["replace", "clean_all", "destroy"])
def test_other_release_paths(release):
    values = arrays()
    refs = [weakref.ref(value) for value in values.values()]
    model = cballs()
    model.set_catalog(**values, take_ownership=True)
    del values
    if release == "replace":
        model.set_catalog(**arrays(), take_ownership=True)
    elif release == "clean_all":
        model.clean_all()
    else:
        del model
    gc.collect()
    assert all(ref() is None for ref in refs)


@pytest.mark.parametrize("bad", ["view", "fortran", "float32", "byteorder", "list", "nonfinite", "mask-dtype", "mask-value", "shape", "gamma-shape", "gamma-missing"])
def test_rejected_transfer_is_atomic_and_keeps_previous_registration(tmp_path, bad):
    model = configured(tmp_path)
    old = arrays()
    old_ref = weakref.ref(old["positions"])
    model.set_catalog(**old, take_ownership=True)
    model.Run(level=["MainLoop"])
    native = model._runtime_bodytable_address()
    before = model.getHistXi2pcf().copy()
    values = arrays()
    if bad == "view": values["positions"] = values["positions"].view()
    elif bad == "fortran": values["positions"] = np.asfortranarray(values["positions"])
    elif bad == "float32": values["positions"] = values["positions"].astype(np.float32)
    elif bad == "byteorder": values["positions"] = values["positions"].astype(">f8")
    elif bad == "list": values["positions"] = values["positions"].tolist()
    elif bad == "nonfinite": values["kappa"][0] = np.nan
    elif bad == "mask-dtype": values["mask"] = values["mask"].astype(bool)
    elif bad == "mask-value": values["mask"][0] = 2
    elif bad == "shape": values["positions"] = np.zeros((48, 2))
    elif bad == "gamma-shape": values["gamma2"] = np.zeros(47)
    elif bad == "gamma-missing": values.pop("gamma2")
    with pytest.raises(CosmoSevereError):
        model.set_catalog(**values, take_ownership=True)
    assert all(value.flags.writeable for value in values.values() if isinstance(value, np.ndarray))
    assert model.catalog_count == 1
    assert model._runtime_bodytable_address() == native
    np.testing.assert_array_equal(model.getHistXi2pcf(), before)
    del old
    assert old_ref() is not None
    model.clean_all()


def test_default_still_copies_and_isolates_caller_mutation(tmp_path):
    values = arrays()
    frozen = {key: value.copy() for key, value in values.items()}
    model = configured(tmp_path / "copy")
    model.set_catalog(**values)
    assert all(value.flags.writeable for value in values.values())
    for value in values.values(): value[:] = 0
    model.Run(level=["MainLoop"])
    result = model.getResults()
    reference = configured(tmp_path / "reference")
    reference.set_catalog(**frozen, take_ownership=True)
    reference.Run(level=["MainLoop"])
    for key, value in result["arrays"].items():
        np.testing.assert_array_equal(value, reference.getResults()["arrays"][key])
    model.clean_all(); reference.clean_all()


@pytest.mark.parametrize("engine", [f"{tree}{field}-2balls-omp"
    for field in ("", "-shear-sphere") for tree in ("octree", "balltree", "kdtree")])
def test_copy_and_transfer_match_across_engines_and_reload(tmp_path, engine):
    if cyballs.search_method_id(engine) < 0:
        pytest.skip(f"{engine} is not built")
    reference = None
    for transfer in (False, True):
        values = arrays()
        model = configured(tmp_path / str(transfer), engine)
        model.set_catalog(**values, take_ownership=transfer)
        del values
        gc.collect()
        for repeat in range(2):
            model.Run(level=["SetNumberThreads"])
            model.Run(level=["MainLoop"])
            result = model.getResults()
            if reference is None:
                reference = result
            assert result["arrays"].keys() == reference["arrays"].keys()
            for key in result["arrays"]:
                np.testing.assert_array_equal(result["arrays"][key], reference["arrays"][key])
            assert result["metadata"]["inputs"]["effective_catalogs"] == reference["metadata"]["inputs"]["effective_catalogs"]
            model.struct_cleanup()
        model.clean_all()


def test_transfer_supports_multiple_catalogs_and_readonly_owned_inputs():
    model = cballs()
    first, second = arrays(), arrays()
    first["positions"].setflags(write=False)
    model.set_catalog(**first, take_ownership=True)
    model.set_catalog(**second, catalog=1, take_ownership=True)
    assert model.catalog_count == 2
    assert all(not value.flags.writeable for value in second.values())
    model.clear_catalogs()
    assert not first["positions"].flags.writeable


def test_driver_option_and_shared_transfer_catalog(tmp_path):
    import kappa_corr_all_engines as driver
    assert not driver.parse_arguments([]).transfer_catalog_ownership
    assert driver.parse_arguments(["--transfer-catalog-ownership"]).transfer_catalog_ownership
    assert inspect.signature(cballs.set_catalog).parameters["take_ownership"].default is False
    engines = ("octree-2balls-omp", "balltree-2balls-omp", "kdtree-2balls-omp")
    catalog = driver.synthetic_healpix_catalog(2)
    run = driver.RunConfig(engines, tmp_path, bins=4, multipoles=2, threads=1,
        theta_min=1., theta_max=100., tree_theta=0., verbose=0, verbose_log=0,
        plots=False, transfer_catalog_ownership=True,
        options=("no-out-Hist", "no-smooth-pivot", "no-one-ball", "no-two-balls"))
    results = driver.run_engine_suite(catalog, run)
    assert set(results) == set(engines)
    assert not catalog.positions.flags.writeable
    assert not catalog.kappa.flags.writeable
    summary = json.loads((tmp_path / "summary.json").read_text())
    assert summary["catalog_registration"] == "transfer"
    assert summary["set_catalog_calls_per_process"] == 1
    assert not summary["failures"]
