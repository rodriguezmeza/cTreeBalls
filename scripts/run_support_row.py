#!/usr/bin/env python3
"""Build one declared dependency row and verify its installed numerical API."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
from platform_support import matrix, receipt, ROOT

def run(row_id, output):
    row=next(r for r in matrix() if r['id']==row_id)
    output=output.resolve();output.mkdir(parents=True,exist_ok=False)
    report=dict(status='RUNNING',row=row_id,commands=[])
    env=dict(os.environ,PYTHONNOUSERSITE='1',PYTHONOPTIMIZE='0')
    for key in ('PYTHONPATH','PYTHONHOME'):env.pop(key,None)
    def call(label,argv,cwd=ROOT):
        log=output/(label+'.log')
        with log.open('w') as f:
            result=subprocess.run(list(map(str,argv)),cwd=cwd,env=env,stdout=f,stderr=subprocess.STDOUT,timeout=1800)
        report['commands'].append(dict(label=label,argv=list(map(str,argv)),returncode=result.returncode))
        if result.returncode:raise RuntimeError(label+' failed; see '+str(log))
        return log.read_text()
    try:
        call('build-environment',[sys.executable,'-m','pip','install',*row['build_requirements']])
        call('native-build',['make','-j2','cballs','cyballs-static-lib',f'PYTHON={sys.executable}'])
        raw=call('fingerprint',['make','--no-print-directory','print-build-fingerprint',f'PYTHON={sys.executable}'])
        build=next(json.loads(v) for v in raw.splitlines() if v.startswith('{'))
        expected=output/'expected-build.json';expected.write_text(json.dumps(build)+'\n')
        call('wheel',[sys.executable,'-m','pip','wheel','--no-build-isolation','--no-deps','--wheel-dir',output/'wheel',ROOT])
        call('runtime-environment',[sys.executable,'-m','venv',output/'runtime'])
        py=output/'runtime/bin/python'
        call('runtime-dependencies',[py,'-m','pip','install',*row['runtime_requirements'],'packaging'])
        wheels=list((output/'wheel').glob('*.whl'))
        if len(wheels)!=1:raise RuntimeError('expected one wheel')
        call('install',[py,'-m','pip','install','--no-deps',wheels[0]])
        outside=output/'outside';outside.mkdir()
        evidence=output/'installed.json'
        call('installed-oracles',[py,ROOT/'scripts/check_installed_package.py','--source-root',ROOT,
            '--expected-build',expected,'--output',evidence],cwd=outside)
        call('receipt',[py,ROOT/'scripts/platform_support.py','receipt','--row',row_id,
            '--evidence',evidence,'--output',output/'platform-result.json'],cwd=outside)
        report['status']='PASS'
    except Exception as e:
        report.update(status='FAIL',error=str(e))
    finally:(output/'run.json').write_text(json.dumps(report,indent=2)+'\n')
    print(report['status'],output)
    return 0 if report['status']=='PASS' else 1
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--row',required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    raise SystemExit(run(a.row,a.output))
