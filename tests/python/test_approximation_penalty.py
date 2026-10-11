"""Measured approximation penalties must expose weak and undefined signals."""
from copy import deepcopy
import json
from pathlib import Path

import numpy as np
import pytest
from cyballs import qualification_report, publish_result_packet, save_result_packet
from test_clustered_catalogs import fixture, run


@pytest.fixture
def reference(tmp_path):
    return run(fixture(41), tmp_path/'reference', True, False, 3, 'exact')


def products(reference, actual, exact, name='signal'):
    target=deepcopy(reference); candidate=deepcopy(reference)
    target['arrays']={name:np.asarray(exact)}
    candidate['arrays']={name:np.asarray(actual)}
    return candidate,target


def test_relative_error_has_no_absolute_denominator_floor(reference):
    candidate,target=products(reference,[2e-30],[1e-30])
    row=qualification_report(candidate,target,rtol=.05,atol=1e-10)['observables']['signal']
    assert row['relative_l2']==1
    assert row['relative_l2_percent']==100
    assert not row['relative_target_met']
    # The legacy absolute-plus-relative acceptance is a separate statistic.
    assert row['status']=='QUALIFIED'


def test_zero_reference_and_invalid_bins_are_not_zero_percent(reference):
    for candidate_values,target_values in [([1.],[0.]),([np.nan],[1.])]:
        candidate,target=products(reference,candidate_values,target_values)
        report=qualification_report(candidate,target,rtol=.05)
        row=report['observables']['signal']
        assert row['relative_l2_percent'] is None
        assert not row['relative_target_met']
        json.dumps(report,allow_nan=False)


def test_weak_shear_component_is_reported_separately(reference):
    target=np.ones((4,3,2,2),dtype=complex);target[3]*=1e-12
    actual=target.copy();actual[3]*=1.2
    candidate,reference=products(reference,actual,target,'upsilon_raw')
    report=qualification_report(candidate,reference,rtol=.05,atol=0)
    assert report['observables']['upsilon_raw']['relative_target_met']
    row=report['observables']['upsilon_raw_component_3']
    assert row['relative_l2_percent']==pytest.approx(20)
    assert not row['relative_target_met']
    assert report['status']=='REJECTED'


def test_reuse_is_an_approximation_control_not_a_different_estimator(reference):
    candidate=deepcopy(reference)
    candidate['metadata']['options']+=',shear-pivot-reuse,scalar-pivot-reuse'
    assert qualification_report(candidate,reference)['status']=='QUALIFIED'


class PacketModel:
    def __init__(self,packet):self.packet=deepcopy(packet)
    def qualifyAgainst(self,reference,**kwargs):
        self.packet['metadata']['qualification']=qualification_report(self.packet,reference,**kwargs)
    def getResults(self):return deepcopy(self.packet)


def test_report_only_retains_results_above_target(reference,tmp_path,capsys):
    candidate,target=products(reference,[2.],[1.])
    exact=tmp_path/'exact';save_result_packet(target,exact)
    path=publish_result_packet(PacketModel(candidate),tmp_path/'report',exact,
                               rtol=.05,atol=0,fail_on_rejection=False)
    assert Path(path).is_file()
    report=json.loads((tmp_path/'report/approximation-penalty.json').read_text())
    assert report['status']=='MEASURED'
    assert report['observables']['signal']['relative_l2_percent']==100
    assert 'EXCEEDS TARGET' in capsys.readouterr().out
    with pytest.raises(ValueError,match='qualification failed'):
        publish_result_packet(PacketModel(candidate),tmp_path/'strict',exact,rtol=.05,atol=0)


def test_no_reference_is_unmeasured(reference,tmp_path,capsys):
    reference['metadata']['qualification']={'status':'UNMEASURED'}
    publish_result_packet(PacketModel(reference),tmp_path/'none',rtol=.05)
    report=json.loads((tmp_path/'none/approximation-penalty.json').read_text())
    assert report['status']=='UNMEASURED'
    assert report['observables']=={}
    assert 'UNMEASURED' in capsys.readouterr().out
