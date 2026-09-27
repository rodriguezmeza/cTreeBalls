import sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
from workload_acceptance import compare


def test_nan_candidate_is_rejected_and_json_safe():
    import json
    result=compare(np.array([1.,2.,0.]),np.array([np.nan,2.,0.]))
    assert not result['passed'] and not result['same_finite_mask']
    json.dumps(result,allow_nan=False)


def test_weak_signed_observable_has_fixed_absolute_and_relative_limits():
    target=np.array([-1e-5,1e-5,0.])
    assert compare(target,target.copy())['passed']
    result=compare(target,target+1e-5)
    assert not result['passed'] and result['weak_bins']==1


def test_scaling_selection_is_nested_and_preserves_within_forest_order():
    from benchmark_scaling import select
    ids=np.array([9,2,9,2,7,9,7,2,7,9])
    data=dict(forest_ids=ids,positions=np.arange(30).reshape(10,3),delta=np.arange(10),weights=np.ones(10))
    small,rows=select(data,2,2);large,larger_rows=select(data,3,2)
    np.testing.assert_array_equal(rows,[1,7,4,8])
    np.testing.assert_array_equal(larger_rows[:len(rows)],rows)
    np.testing.assert_array_equal(small['positions'],data['positions'][rows])
