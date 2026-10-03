"""Keep optional external comparisons working across the naming boundary."""
from pathlib import Path
import sys
from types import SimpleNamespace
import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'python'))
import dual_node_compat as compat


@pytest.mark.parametrize('kind', ['NN', 'KK', 'GG', 'NNN', 'KKK', 'GGG'])
def test_translates_only_backend_tolerances(monkeypatch, kind):
    seen = {}
    def construct(**kwargs):
        seen.update(kwargs)
        return seen
    monkeypatch.setattr(compat, '_backend', lambda: SimpleNamespace(**{kind+'Correlation': construct}))
    result = compat.correlation(kind, bin_theta=.12, angle_theta=.03, nbins=4)
    assert result['nbins'] == 4
    assert sorted(v for k,v in result.items() if k != 'nbins') == [.03, .12]
    assert 'bin_theta' not in result and 'angle_theta' not in result


@pytest.mark.parametrize('name', ['bin_theta', 'angle_theta'])
@pytest.mark.parametrize('value', [-1, float('nan'), float('inf')])
def test_invalid_tolerances_fail_before_import(monkeypatch, name, value):
    monkeypatch.setattr(compat, '_backend', lambda: pytest.fail('must validate before import'))
    with pytest.raises(ValueError, match=name):
        compat.correlation('KK', **{name: value})


def test_unknown_estimator_fails():
    with pytest.raises(ValueError, match='Unsupported'):
        compat.correlation('forest-anisotropic')


def test_installed_backend_pair_and_triplet_computation():
    try:
        compat.version()
    except compat.BackendUnavailable as exc:
        pytest.skip(str(exc))
    # Weighted raw scalar pairs have an independent, direct numerator oracle.
    x = np.array([0., 1., 2.5, 4., 6.])
    k = np.array([1., -2., .5, 3., -1.])
    w = np.array([1., 2., 3., 4., 5.])
    cat = compat.catalog(x=x, y=x*.1, k=k, g1=k*.01, g2=k*.02, w=w)
    for kind in ('KK','GG'):
        result = compat.correlation(kind, min_sep=.1, max_sep=10., nbins=1)
        result.process(cat, num_threads=1)
        assert result.npairs.sum() == 10
        assert np.isfinite(result.weight).all()
        if kind=='KK':
            num=sum(w[i]*w[j]*k[i]*k[j] for i in range(5) for j in range(i+1,5))
            den=sum(w[i]*w[j] for i in range(5) for j in range(i+1,5))
            np.testing.assert_allclose(result.xi[0], num/den, rtol=1e-13)
    for kind in ('KKK','GGG'):
        result = compat.correlation(kind, bin_type='LogMultipole', min_sep=.1, max_sep=10., nbins=2, max_n=2)
        result.process(cat, num_threads=1)
        assert result.ntri.sum() > 0
        assert np.isfinite(result.weight).all()
