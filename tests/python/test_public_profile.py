"""The selected branch must match its buildable, documented source layout."""
import json
from pathlib import Path
import subprocess
import sys

import pytest

ROOT=Path(__file__).resolve().parents[2]
PROFILE=json.loads((ROOT/'capabilities/build_profile.json').read_text())


def test_addon_directories_match_the_selected_profile():
    actual={p.name for p in (ROOT/'addons').iterdir() if p.is_dir()}
    assert actual==set(PROFILE['enabled_addon_directories'])
    assert 'python_env' not in actual
    assert {p.name for p in (ROOT/'support').iterdir() if p.is_dir()}==set(PROFILE['shared_support_directories'])


def test_python_tests_and_binding_sources_are_separated():
    assert not list((ROOT/'python').glob('*.py'))
    assert all(p.is_relative_to(ROOT/'tests/python') for p in (ROOT/'tests').rglob('*.py'))
    for name in ('cyballs.pyx','ccyballs.pxd.in','resource_api.pxi','qualification_api.pxi'):
        assert (ROOT/'python'/name).is_file()


def test_generated_profile_reference_is_current():
    subprocess.run([sys.executable,ROOT/'scripts/generate_profile_docs.py','--check'],check=True)


def test_every_capability_source_pattern_has_a_shipped_match():
    data=json.loads((ROOT/'capabilities/engines.json').read_text())
    for engine in data['engines']:
        for pattern in engine['sources']:
            assert any(ROOT.glob(pattern)),(engine['name'],pattern)
        for test in engine['gate']['tests']:
            assert (ROOT/test).is_file(),test


def test_inactive_addon_override_fails_explicitly():
    p=subprocess.run(['make','--no-print-directory','-s','BALLSON=1','print-cyballs-build-env'],
                     cwd=ROOT,text=True,capture_output=True,timeout=30)
    assert p.returncode!=0
    assert 'does not ship these inactive addons' in p.stdout+p.stderr


def test_live_help_and_registry_match_shipped_profile(tmp_path):
    from cyballs import build_info, search_method_id
    data=json.loads((ROOT/'capabilities/engines.json').read_text())
    sys.path.insert(0,str(ROOT/'scripts'))
    from capabilities_generated import expected_registry
    expected=expected_registry(build_info()['resolved_settings'])
    assert len(expected)==41
    assert all(search_method_id(e['name'])==e['id'] for e in data['engines'])
    outputs={}
    for option in ('make-info','print-options','print-search-methods'):
        p=subprocess.run([ROOT/'cballs',f'options={option}',f'rootDir={tmp_path/option}'],
                         text=True,capture_output=True,timeout=30)
        assert p.returncode in (0,1),p.stdout+p.stderr
        outputs[option]=p.stdout+p.stderr
    import re
    methods=set(re.findall(r'^- (\S+) \(id=',outputs['print-search-methods'],re.M))
    assert methods==set(expected)
    for name in PROFILE['inactive_switches']:
        assert not re.search(r'^\s*'+name+r'\s*=',outputs['make-info'],re.M),name
    for option in ('NNLandySzalay1','NNStandard','smooth-min-cell','pivot-loop'):
        assert not re.search(r'^- '+option+r' \[',outputs['print-options'],re.M)
    gsl=build_info()['resolved_settings']['GSLINTERNAL']
    assert ('using internal GSL' in outputs['make-info'])==(gsl=='1')
