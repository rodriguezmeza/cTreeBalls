"""Profile equality and execution coverage must fail closed."""
from pathlib import Path
import sys
import pytest
sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'scripts'))
from check_installed_package import check_profile
from capabilities_generated import ENGINES, expected_registry
from active_release_gate import regression_coverage, ROOT


def test_installed_profile_is_exact_not_a_count():
    intended={'resolved_settings':{'CCFLAG':'-DADDONS -DLYAFORESTOMP'}}
    expected=expected_registry(intended['resolved_settings'])
    assert check_profile(intended,intended,lambda n:expected.get(n,-1))==expected
    with pytest.raises(AssertionError,match='profile mismatch'):
        check_profile({'resolved_settings':{'CCFLAG':''}},intended,lambda n:expected.get(n,-1))
    extra=next(n for n in ENGINES if n not in expected)
    with pytest.raises(AssertionError,match='registry mismatch'):
        check_profile(intended,intended,lambda n:55 if n==extra else expected.get(n,-1))


def test_declared_module_requires_successful_execution():
    registry={'lya-2pcf-omp':169}
    plan=regression_coverage(registry,[])
    assert 'tests/make_tests/test_lya_pair_cells.py' in plan['missing']
    rows=[dict(status='PASS',argv=['python','-m','pytest',str(ROOT/p)]) for p in plan['required']]
    assert not regression_coverage(registry,rows)['missing']
    rows[-1]['status']='FAIL'
    assert regression_coverage(registry,rows)['missing']
