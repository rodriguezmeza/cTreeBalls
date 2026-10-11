"""Reject stale evidence; validate benchmark scheduling and interprocess locking."""
import copy,json,os,subprocess,sys
from pathlib import Path
import pytest
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from acceptance_metadata import REQUIRED,acceptance_policy,source_identity
from benchmark_dedicated import configurations,schedule,summary,node_lock,validate_prerequisites
from platform_support import matrix,collect
from release_readiness import collect as release_collect

def complete_platforms(tmp_path):
    paths=[]
    for row in matrix():
        p=tmp_path/(row['id']+'.json');p.write_text(json.dumps(dict(schema_version=1,row=row['id'],status='PASS',source_id='current',acceptance_policy=acceptance_policy())));paths.append(p)
    return paths

def test_platform_requires_each_unique_current_row(tmp_path):
    paths=complete_platforms(tmp_path)
    assert collect(paths,'current')['status']=='PASS'
    assert collect(paths[:-1],'current')['status']=='INCOMPLETE'
    assert collect(paths+[paths[0]],'current')['status']=='FAIL'
    assert collect(paths,'other')['status']=='FAIL'

@pytest.mark.parametrize('key,value',[('status','SKIPPED'),('status','FAIL'),('acceptance_policy',{}),('schema_version',9)])
def test_platform_damage_cannot_approve(tmp_path,key,value):
    paths=complete_platforms(tmp_path);data=json.loads(paths[0].read_text());data[key]=value;paths[0].write_text(json.dumps(data))
    assert collect(paths,'current')['status']=='FAIL'

def test_policy_is_copied_and_sensitive_to_case_controls(monkeypatch):
    import accuracy_cases
    first=acceptance_policy();damaged=copy.deepcopy(first);damaged['bounded_accuracy']['cases'][0]['theta']=99
    assert acceptance_policy()==first
    monkeypatch.setitem(accuracy_cases.SCALAR[0],'theta',.987)
    assert acceptance_policy()['id']!=first['id']

def test_release_rejects_stale_policy(tmp_path):
    paths={}
    for name,checks in REQUIRED.items():
        p=tmp_path/(name+'.json');paths[name]=p
        p.write_text(json.dumps(dict(schema_version=1,component=name,source_id='current',status='PASS',acceptance_policy=acceptance_policy(),checks={c:'PASS' for c in checks},commands=[dict(returncode=0)])))
    assert release_collect(paths,'current')['status']=='PASS'
    data=json.loads(paths['package'].read_text());data['acceptance_policy']['id']='stale';paths['package'].write_text(json.dumps(data))
    assert release_collect(paths,'current')['status']=='FAIL'

def test_identity_covers_ci_examples_and_manual(tmp_path):
    for name in ('Makefile','Makefile_tests','Makefile_settings','Makefile_machine','setup.py','pyproject.toml','MANIFEST.in'):(tmp_path/name).write_text('')
    for name in ('.github/test.yml','examples/test.ipynb','docs/man/cballs.m','capabilities/support_matrix.json'):
        path=tmp_path/name;path.parent.mkdir(parents=True,exist_ok=True);before=source_identity(tmp_path)
        path.write_text('first');after=source_identity(tmp_path);assert before!=after
        path.write_text('second');assert source_identity(tmp_path)!=after

def test_schedule_has_every_repeat_in_reproducible_random_order():
    configs=configurations(['kdtree-2balls-omp','octree-2balls-omp']);rows=schedule(configs,7,91)
    assert len(rows)==7*len(configs) and len(set(rows))==len(rows)
    assert rows==schedule(configs,7,91) and rows!=schedule(configs,7,92)
    for c in configs:assert sum(tuple(r[:3])==c for r in rows)==7
    with pytest.raises(ValueError):configurations(['bogus'])
    with pytest.raises(ValueError):configurations(['kdtree-2balls-omp']*2)

def test_summary_retains_outlier_and_actual_count():
    s=summary([1,2,3,4,100]);assert s['n']==5 and s['median']==3 and s['max']==100 and s['mad']==1
    with pytest.raises(ValueError):summary([1,float('nan')])

def test_lock_excludes_independent_process(tmp_path):
    path=tmp_path/'benchmark.lock';code='import sys;from pathlib import Path;from benchmark_dedicated import node_lock\nwith node_lock(Path(sys.argv[1])):pass'
    env=dict(os.environ,PYTHONPATH=str(ROOT/'scripts'))
    with node_lock(path):
        child=subprocess.run([sys.executable,'-c',code,str(path)],env=env,capture_output=True,text=True)
        assert child.returncode!=0 and 'node lock' in child.stderr
    child=subprocess.run([sys.executable,'-c',code,str(path)],env=env,capture_output=True,text=True)
    assert child.returncode==0,child.stderr

@pytest.mark.parametrize('change',['source','build','policy','missing-component','failed-accuracy'])
def test_timing_requires_matching_prior_correctness(change):
    release=dict(status='PASS',source_id='source',acceptance_policy=acceptance_policy(),components={k:dict(status='PASS') for k in REQUIRED})
    accuracy=dict(status='PASS',acceptance_policy=acceptance_policy(),build=dict(id='build'),cases=[dict(acceptance_status='PASS')],exact_cross_engine=[dict(passed=True)])
    validate_prerequisites(release,accuracy,'build','source')
    if change=='source':release['source_id']='old'
    if change=='build':accuracy['build']['id']='old'
    if change=='policy':accuracy['acceptance_policy']['id']='old'
    if change=='missing-component':release['components'].pop('sanitizers')
    if change=='failed-accuracy':accuracy['cases'][0]['acceptance_status']='FAIL'
    with pytest.raises(ValueError):validate_prerequisites(release,accuracy,'build','source')

def test_hosted_release_requires_complete_platform_receipt(tmp_path):
    paths={}
    for name,checks in REQUIRED.items():
        p=tmp_path/(name+'.json');paths[name]=p
        p.write_text(json.dumps(dict(schema_version=1,component=name,source_id='current',status='PASS',acceptance_policy=acceptance_policy(),checks={c:'PASS' for c in checks},commands=[dict(returncode=0)])))
    missing=tmp_path/'support.json'
    assert release_collect(paths,'current')['platform_support']['status']=='NOT_RUN'
    assert release_collect(paths,'current',missing)['status']=='INCOMPLETE'
    support=collect(complete_platforms(tmp_path),'current');missing.write_text(json.dumps(support))
    assert release_collect(paths,'current',missing)['status']=='PASS'
    support['rows'].pop(next(iter(support['rows'])));missing.write_text(json.dumps(support))
    assert release_collect(paths,'current',missing)['status']=='INCOMPLETE'

def test_platform_rejects_modified_policy_with_unchanged_id(tmp_path):
    paths=complete_platforms(tmp_path);data=json.loads(paths[0].read_text())
    data['acceptance_policy']['survey_defaults']['rtol']=10
    paths[0].write_text(json.dumps(data))
    assert collect(paths,'current')['status']=='FAIL'


def test_platform_receipt_rejects_failed_installed_case(tmp_path,monkeypatch):
    import platform_support as support
    monkeypatch.setattr(support.platform,'system',lambda:'Darwin')
    monkeypatch.setattr(support.platform,'machine',lambda:'arm64')
    monkeypatch.setattr(support.platform,'python_version',lambda:'3.13.0')
    path=tmp_path/'installed.json';path.write_text(json.dumps(dict(status='PASS',cases=[dict(status='FAIL')]*11)))
    with pytest.raises(ValueError,match='incomplete'):support.receipt('darwin-3.13-current',path)


def test_setup_queries_use_backend_interpreter_despite_conflicting_environment(tmp_path):
    fake=tmp_path/'wrong-python';fake.write_text('#!/bin/sh\necho WRONG_BUILD_INTERPRETER >&2\nexit 99\n');fake.chmod(0o755)
    env=dict(os.environ,PYTHON=str(fake),MAKEFLAGS='-e')
    result=subprocess.run([sys.executable,str(ROOT/'setup.py'),'--name'],cwd=ROOT,env=env,
                          capture_output=True,text=True,timeout=120)
    assert result.returncode==0,result.stdout+result.stderr
    assert 'WRONG_BUILD_INTERPRETER' not in result.stdout+result.stderr
    assert result.stdout.strip().endswith('cyballs')
