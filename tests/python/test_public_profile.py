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
    programs = list((ROOT/'tests').rglob('*.py'))
    assert programs
    assert all(path.parent == ROOT/'tests/python' for path in programs)
    assert list((ROOT/'tests/python').glob('test_*.py'))
    assert list((ROOT/'tests/python').glob('benchmark_*.py'))
    assert (ROOT/'tests/python/test_benchmark_timing.py').is_file()
    assert (ROOT/'tests/python/kappa_corr_all_engines.py').is_file()
    assert (ROOT/'tests/python/shear_fits_catalog.py').is_file()
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
    assert len(expected)==42
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


def test_public_documentation_uses_project_terminology():
    import re
    forbidden = re.compile('treecorr', re.I)
    for directory in ('docs', 'tutor', 'tests/python', 'support', 'addons'):
        for path in (ROOT/directory).rglob('*'):
            if not path.is_file() or any(x in ('_build', 'build', '__pycache__') for x in path.parts):
                continue
            if path.suffix in ('.md', '.rst') or path.name.lower().startswith('readme'):
                text = path.read_text(errors='replace')
                assert not forbidden.search(text), path.relative_to(ROOT)
                for name in re.findall(r'(?<![\w/])python/([\w.-]+\.py)\b', text):
                    assert not (ROOT/'tests/python'/name).is_file(), (path.relative_to(ROOT), name)
    body, acknowledgments = (ROOT/'README.md').read_text().split('## Acknowledgments', 1)
    assert not forbidden.search(body)
    assert forbidden.search(acknowledgments)


def test_man_inventory_and_projection_include_match_the_profile():
    import re
    engines = json.loads((ROOT/'capabilities/engines.json').read_text())['engines']
    manual = (ROOT/'docs/man/cballs.m').read_text()
    assert manual == (ROOT/'docs/man/cballs.1').read_text()
    inventory = manual.split('.SH ACTIVE SEARCH METHODS\n', 1)[1].split('.SH COMMON PARAMETERS', 1)[0]
    assert set(re.findall(r'^\.B (\S+)$', inventory, re.M)) == {e['name'] for e in engines}
    assert 'include "../tests/python/lya_projection.py"' in (ROOT/'python/cyballs.pyx').read_text()
    assert not (ROOT/'python/lya_projection.py').exists()
