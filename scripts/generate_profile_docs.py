#!/usr/bin/env python3
"""Generate the selected-profile directory and runtime-method reference."""
import argparse
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]


def render():
    profile=json.loads((ROOT/'capabilities/build_profile.json').read_text())
    engines=json.loads((ROOT/'capabilities/engines.json').read_text())['engines']
    intro=('This testing branch contains the enabled Makefile addon set. '
           'The core octree method is always available. Shared implementation '
           'dependencies live in support/ and do not register standalone methods. '
           'The Cython binding sources stay in python/; Python analysis, benchmark '
           'and regression scripts live in tests/python/. Inactive addons and '
           'addons/python_env are excluded from the branch.')
    md=['# Active source profile\n\n',intro+'\n\n',
        'The development checkout includes the enabled bundled GSL and CFITSIO '
        'sources. A public source distribution excludes those libraries and sets '
        '`GSLINTERNAL=0, CFITSIOLIBON=0` in its staged Makefiles.\n\n',
        '## Shipped addon directories\n\n']
    md+=['- `addons/'+name+'`\n' for name in profile['enabled_addon_directories']]
    md+=['\n## Shared implementation dependencies\n\n']
    md+=['- `support/'+name+'`\n' for name in profile['shared_support_directories']]
    md+=['\n## Enabled runtime methods\n\n', '| Name | ID | Geometry and estimator |\n| --- | ---: | --- |\n']
    md += [f'| `{e["name"]}` | {e["id"]} | {e["geometry"]}; {e["correlations"]} |\n' for e in engines]
    md+=['\nUse `options=make-info`, `options=print-options` and '
         '`options=print-search-methods` to inspect a particular compiled executable. '
         'The latter is authoritative after any supported build overrides.\n']
    rst=['Active Source Profile\n=====================\n\n',intro+'\n\n',
         'Runtime discovery\n-----------------\n\n.. code-block:: sh\n\n'
         '   ./cballs options=make-info\n   ./cballs options=print-options\n'
         '   ./cballs options=print-search-methods\n\n',
         'The selected profile has '+str(len(engines))+' registered methods. '
         'The executable reports only methods actually compiled; unsupported '
         'standalone-addon overrides fail explicitly. The development profile uses '
         'bundled GSL/CFITSIO, while the public source archive uses external libraries.\n\n',
         'Enabled methods\n---------------\n\n']
    for e in engines:
        rst += ['``'+e['name']+'`` (ID '+str(e['id'])+')\n',
                '    '+e['geometry']+'. '+e['correlations']+'.\n\n']
    rst+=['Source organization\n-------------------\n\n',
          '* ``addons/`` contains enabled engines and their I/O/binding dependencies.\n'
          '* ``support/`` contains shared builders and compatibility kernels.\n'
          '* ``tests/python/`` contains all Python testing and benchmark scripts.\n'
          '* ``tests/make_tests/`` contains shell/native-test launchers.\n'
          '* ``python/`` contains Cython source and declaration files only.\n']
    return {ROOT/'docs/ACTIVE_PROFILE.md':''.join(md),ROOT/'docs/active_profile.rst':''.join(rst)}


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--check',action='store_true')
    args=parser.parse_args();stale=[]
    for p,text in render().items():
        if not p.exists() or p.read_text()!=text:
            if args.check:stale.append(str(p.relative_to(ROOT)))
            else:p.write_text(text)
    if stale:raise SystemExit('stale profile documentation: '+', '.join(stale))
    print('PASS: active-profile documentation')

if __name__=='__main__':main()
