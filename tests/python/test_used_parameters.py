"""Used-value records survive driver cleanup and describe actual native controls."""
from pathlib import Path
from types import SimpleNamespace
import json
import sys

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'tests/python'))
sys.path.insert(0, str(ROOT/'tests/python'))
import kappa_corr_all_engines as kappa
import shear_corr_all_engines as shear
from used_parameters import save_native_parameters, _write


def read_values(path):
    native, extra = {}, {}
    for line in path.read_text().splitlines():
        if '=' not in line:
            continue
        key, value = line.split('=', 1)
        if key.lstrip().startswith('%'):
            extra[key.lstrip()[1:].strip()] = json.loads(value)
        else:
            native[key.strip()] = value.strip()
    return native, extra


@pytest.mark.parametrize('tree', ['octree', 'kdtree', 'balltree'])
@pytest.mark.parametrize('field', ['shear', 'kappa'])
def test_native_record_survives_cleanup_with_defaults_and_overrides(tmp_path, monkeypatch, tree, field):
    for prefix in ('SHEAR', 'SCALAR'):
        for suffix in ('PIVOT_TOL', 'BIN_THETA'):
            monkeypatch.delenv(f'CBALLS_{prefix}_{suffix}', raising=False)
    rng = np.random.default_rng(761)
    positions = rng.normal(size=(96, 3))
    positions /= np.linalg.norm(positions, axis=1)[:, None]
    g1, g2 = .01*positions[:, 0], .02*positions[:, 1]
    engine = tree + ('-shear-sphere-2balls-omp' if field == 'shear' else '-2balls-omp')
    from cyballs import search_method_id
    if search_method_id(engine) < 0:
        pytest.skip(engine + ' is not compiled')
    out = tmp_path/field
    for override in (False, True):
        prefix = 'SHEAR' if field == 'shear' else 'SCALAR'
        if override:
            monkeypatch.setenv(f'CBALLS_{prefix}_PIVOT_TOL', '.2')
            monkeypatch.setenv(f'CBALLS_{prefix}_BIN_THETA', '.1')
        phase = .2 if override else (3.0 if field == 'shear' else .1)
        slop = .1 if override else (.5 if field == 'shear' else 0.0)
        if field == 'shear':
            catalog = shear.ShearCatalog(positions, g1, g2, geometry='sphere').normalized()
            config = shear.RunConfig(output_dir=out, bins=4, multipoles=2,
                min_sep=.1, max_sep=1.8, sep_units='radian', threads=1,
                statistics='3pcf', phi_bins=16, plots=False,
                options=('shear-pivot-reuse', 'no-smooth-pivot')).normalized()
            if override:
                config.nsmooth, config.tree_theta = 11, .7
            result = shear.run_ctreeballs(catalog, config, engine)
            limits = config.native_limits('sphere')
        else:
            catalog = kappa.KappaCatalog(positions, g1)
            config = kappa.RunConfig((engine,), out, bins=4, multipoles=2,
                theta_min=np.degrees(.1), theta_max=np.degrees(1.8), threads=1,
                plots=False, options=('only-3pcf', 'scalar-pivot-reuse',
                                      'no-smooth-pivot', 'no-out-Hist')).normalized()
            if override:
                config.nsmooth, config.tree_theta = 11, .7
            result = kappa.run_engine_suite(catalog, config)[engine]
            params = kappa.engine_parameters(config, engine, False)
            limits = (params['rminHist'], params['rangeN'])
        path = out/engine/'parameters-usedvalues'
        assert Path(result['parameters_file']) == path
        native, extra = read_values(path)
        assert extra['run_stage'] == 'MainLoop complete'
        assert native['searchMethod'] == engine
        assert int(native['nsmooth']) == (11 if override else 8)
        assert float(native['theta']) == (.7 if override else 1.0)
        assert (float(native['rminHist']), float(native['rangeN'])) == limits
        assert extra[f'used.CBALLS_{prefix}_PIVOT_TOL'] == phase
        assert extra[f'used.CBALLS_{prefix}_BIN_THETA'] == slop
        assert extra['native.catalog_sizes'] == [len(positions)]
        assert extra['parallel.rank'] == 0
        assert 'verbose' in native and 'preScript' in native  # Full native file.
        if field == 'shear':
            assert not Path(native['rootDir']).exists()  # Native scratch was removed.
        assert not list(path.parent.glob('.parameters-usedvalues-*'))


def test_nonroot_rank_does_not_read_or_publish_parameters(tmp_path):
    model = SimpleNamespace(getRunMetadata=lambda: {'parallel': {'rank': 3}})
    assert save_native_parameters(model, SimpleNamespace(output_dir=tmp_path)) is None
    assert not list(tmp_path.iterdir())


def test_failed_search_retains_startup_record(tmp_path, monkeypatch):
    import cyballs
    native_class = cyballs.cballs
    engine = 'balltree-shear-sphere-2balls-omp'
    if cyballs.search_method_id(engine) < 0:
        pytest.skip(engine + ' is not compiled')
    path = tmp_path/engine/'parameters-usedvalues'
    class FailingSearch:
        def __init__(self):
            self.native = native_class()
        def __getattr__(self, name):
            return getattr(self.native, name)
        def Run(self, level):
            if level == ['MainLoop']:
                assert path.is_file()  # Visible before counting starts.
                raise RuntimeError('intentional search failure')
            return self.native.Run(level=level)
    monkeypatch.setattr(cyballs, 'cballs', FailingSearch)
    config = shear.RunConfig(output_dir=tmp_path, bins=4, multipoles=2,
                             threads=1, plots=False).normalized()
    catalog = shear.synthetic_spherical_shear_catalog(32).normalized()
    with pytest.raises(RuntimeError, match='intentional search failure'):
        shear.run_ctreeballs(catalog, config, engine)
    native, extra = read_values(path)
    assert extra['run_stage'] == 'initialized'
    assert int(native['nsmooth']) == 8


def test_atomic_output_failure_keeps_previous_record(tmp_path, monkeypatch):
    import used_parameters
    old = _write(tmp_path, ['previous record\n'])
    def fail(*args):
        raise OSError('write denied')
    monkeypatch.setattr(used_parameters.os, 'replace', fail)
    with pytest.raises(OSError, match='write denied'):
        _write(tmp_path, ['new record\n'])
    assert Path(old).read_text() == 'previous record\n'
    assert list(tmp_path.iterdir()) == [Path(old)]
