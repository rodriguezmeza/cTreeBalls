"""Memory controls preserve the selected catalog and native science outputs."""
import argparse
import gc
import weakref

import numpy as np
import pytest

from catalog_memory import (add_memory_arguments, gather_field, healpix_positions, register_catalog, register_forest_catalog)


@pytest.mark.parametrize('nested', [False, True])
def test_chunked_coordinates_and_endian_gather(nested):
    hp = pytest.importorskip('healpy')
    pixels = np.array([31, 4, 765, 0, 15, 201, 55], dtype=np.int64)
    calls = []
    class BoundedHealpy:
        @staticmethod
        def pix2vec(nside, selected, nest):
            calls.append(len(selected))
            return hp.pix2vec(nside, selected, nest=nest)
    result = healpix_positions(BoundedHealpy, 8, pixels, nested=nested, chunk_pixels=3)
    np.testing.assert_array_equal(result, np.column_stack(hp.pix2vec(8, pixels, nest=nested)))
    assert calls == [3, 3, 1] and result.flags.owndata
    field = np.arange(768, dtype='>f8')
    np.testing.assert_array_equal(gather_field(field, pixels, 3), field[pixels])


def test_transfer_prepares_views_but_preserves_owned_buffers():
    class Model:
        def set_catalog(self, positions, *, kappa, take_ownership=False):
            self.positions, self.kappa, self.transfer = positions, kappa, take_ownership
    model = Model()
    pos = np.ones((10, 3))
    complex_field = np.arange(10) + 1j
    register_catalog(model, True, pos, kappa=complex_field.real)
    assert model.transfer and model.positions is pos
    assert model.kappa.flags.owndata and not np.shares_memory(model.kappa, complex_field)
    np.testing.assert_array_equal(model.kappa, complex_field.real)
    register_catalog(model, False, pos, kappa=complex_field.real)
    assert not model.transfer and np.shares_memory(model.kappa, complex_field)


def test_parser_defaults_and_invalid_chunk():
    parser = argparse.ArgumentParser()
    add_memory_arguments(parser, healpix=True)
    args = parser.parse_args([])
    assert not args.transfer_catalog_ownership and args.fits_chunk_pixels == 1048576
    with pytest.raises(SystemExit):
        parser.parse_args(['--fits-chunk-pixels', '0'])


def test_forest_transfer_lifetime_and_atomic_validation():
    cyballs = pytest.importorskip('cyballs')
    model = cyballs.cballs()
    pos = np.ones((6, 3))
    delta = np.linspace(-.2, .2, 6).copy()
    weight = np.ones(6)
    ids = np.repeat(np.array([2**53+1, 2**53+3, 2**53+5], dtype=np.int64), 2)
    model.set_forest_catalog(pos, delta, weight, ids, take_ownership=True)
    assert not pos.flags.writeable and not delta.flags.writeable and not weight.flags.writeable
    assert ids.flags.writeable  # IDs have private compact storage.
    ref = weakref.ref(pos)
    del pos
    gc.collect()
    assert ref() is not None
    model.struct_cleanup()
    assert ref() is not None
    replacement = np.ones((6, 3))
    bad_weight = np.ones(6)
    bad_weight[-1] = -1
    with pytest.raises(Exception, match='weights'):
        model.set_forest_catalog(replacement, delta.copy(), bad_weight, ids, take_ownership=True)
    assert replacement.flags.writeable and bad_weight.flags.writeable and model.catalog_count == 1
    with pytest.raises(Exception, match='project delta before'):
        model.set_forest_catalog(replacement, delta.copy(), weight.copy(), ids,
                                 take_ownership=True, project_delta=True)
    model.clear_catalogs()
    gc.collect()
    assert ref() is None


@pytest.mark.parametrize('engine', ['lya-2pcf-omp', 'lya-1d-tree-2pcf-omp', 'lya-los-tree-3pcf-omp'])
def test_forest_copy_transfer_products(engine, tmp_path):
    cyballs = pytest.importorskip('cyballs')
    from lya_corr_all_engines import synthetic_catalog, RunConfig, engine_parameters, read_products
    from kappa_corr_all_engines import discover_cython_methods
    if engine not in discover_cython_methods([engine]):
        pytest.skip('engine absent from build')
    catalog = synthetic_catalog(4, 5, 17)
    products = []
    for transfer in (False, True):
        model = cyballs.cballs()
        output = tmp_path / str(transfer)
        output.mkdir()
        config = RunConfig((engine,), output, threads=1, plots=False)
        try:
            register_forest_catalog(model, transfer, catalog.positions.copy(),
                                    catalog.delta.copy(), catalog.weights.copy(), catalog.forest_ids)
            model.set(engine_parameters(config, engine, output))
            model.Run()
            products.append(read_products(output, engine))
        finally:
            model.clean_all()
    assert products[0].keys() == products[1].keys()
    for key in products[0]:
        np.testing.assert_array_equal(products[0][key], products[1][key])
