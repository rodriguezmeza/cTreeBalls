"""Missing, stale, failed or incomplete component receipts cannot approve release."""
import json
from pathlib import Path
import subprocess
import sys
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
from release_readiness import REQUIRED, collect
from acceptance_metadata import acceptance_policy


def receipts(tmp_path):
    paths={}
    for name,checks in REQUIRED.items():
        data=dict(schema_version=1,component=name,source_id='current',status='PASS',acceptance_policy=acceptance_policy(),
                  checks={c:'PASS' for c in checks},commands=[dict(returncode=0)])
        path=tmp_path/(name+'.json');path.write_text(json.dumps(data));paths[name]=path
    return paths


def test_requires_all_components(tmp_path):
    paths=receipts(tmp_path)
    assert collect(paths,'current')['status']=='PASS'
    for name in REQUIRED:
        partial=dict(paths);partial.pop(name)
        assert collect(partial,'current')['status']=='INCOMPLETE'


@pytest.mark.parametrize('state', ['FAIL','RUNNING','UNAVAILABLE','SKIPPED'])
@pytest.mark.parametrize('component', REQUIRED)
def test_unsuccessful_component_never_approves(tmp_path,state,component):
    paths=receipts(tmp_path);data=json.loads(paths[component].read_text())
    data['status']=state;paths[component].write_text(json.dumps(data))
    assert collect(paths,'current')['status']==('FAIL' if state=='FAIL' else 'INCOMPLETE')


@pytest.mark.parametrize('component', REQUIRED)
def test_missing_required_check_and_failed_command_rejected(tmp_path,component):
    paths=receipts(tmp_path);data=json.loads(paths[component].read_text())
    for required in REQUIRED[component]:
        changed=dict(data,checks=dict(data['checks']));changed['checks'].pop(required)
        paths[component].write_text(json.dumps(changed))
        assert collect(paths,'current')['status']=='INCOMPLETE'
    data['commands'][0]['returncode']=3;paths[component].write_text(json.dumps(data))
    assert collect(paths,'current')['status']=='FAIL'


@pytest.mark.parametrize('damage', ['old-source','wrong-component','invalid-json'])
def test_mixed_or_corrupted_evidence_fails(tmp_path,damage):
    paths=receipts(tmp_path);p=paths['package'];data=json.loads(p.read_text())
    if damage=='old-source':data['source_id']='old'
    if damage=='wrong-component':data['component']='active'
    p.write_text('{' if damage=='invalid-json' else json.dumps(data))
    assert collect(paths,'current')['status']=='FAIL'


def test_collect_exit_status_survives_optimized_python(tmp_path):
    root=Path(__file__).resolve().parents[2]
    result=subprocess.run([sys.executable,'-OO',str(root/'scripts/release_readiness.py'),
                           'collect','--output',str(tmp_path/'result.json')],capture_output=True,text=True)
    assert result.returncode!=0
    assert json.loads((tmp_path/'result.json').read_text())['status']=='INCOMPLETE'


def test_software_pass_keeps_science_rejection_visible(tmp_path):
    paths=receipts(tmp_path);p=paths['active'];data=json.loads(p.read_text())
    data['science_acceptance']='REJECTED';p.write_text(json.dumps(data))
    result=collect(paths,'current')
    assert result['status']=='PASS' and result['science_acceptance']=='REJECTED'
