#!/usr/bin/env python3
"""Declared platform matrix, machine-readable receipts and complete coverage checks."""
import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform
import subprocess
import sys
from acceptance_metadata import source_identity, acceptance_policy
ROOT=Path(__file__).resolve().parents[1]

def matrix():
    data=json.loads((ROOT/'capabilities/support_matrix.json').read_text())
    if data['schema_version']!=1 or len({r['id'] for r in data['rows']})!=len(data['rows']):
        raise ValueError('invalid support matrix')
    return data['rows']

def receipt(row_id, evidence):
    row=next(r for r in matrix() if r['id']==row_id)
    actual=dict(system=platform.system(),machine=platform.machine(),python=platform.python_version())
    expected_machine={'arm64':{'arm64','aarch64'},'x86_64':{'x86_64','AMD64'}}[row['machine']]
    if actual['system']!=row['system'] or actual['machine'] not in expected_machine or '.'.join(actual['python'].split('.')[:2])!=row['python']:
        raise ValueError('host does not match declared row: '+str(actual))
    data=json.loads(evidence.read_text())
    if (data.get('status')!='PASS' or len(data.get('cases',[]))<11
            or any(c.get('status')!='PASS' for c in data['cases'])):
        raise ValueError('installed scientific verification is incomplete')
    from packaging.requirements import Requirement
    versions={}
    for spec in row['runtime_requirements']:
        req=Requirement(spec);version=importlib.metadata.version(req.name);versions[req.name]=version
        if version not in req.specifier:raise ValueError('runtime version outside row: '+spec)
    built=data['build']['toolchain']['build_packages']
    for spec in row['build_requirements']:
        req=Requirement(spec)
        if req.name not in built or built[req.name] not in req.specifier:
            raise ValueError('build version outside row: '+spec)
    return dict(schema_version=1,status='PASS',row=row_id,source_id=source_identity(),
        acceptance_policy=acceptance_policy(),host=actual,versions=versions,
        installed_build=data['build'],evidence=str(evidence.resolve()),
        evidence_sha256=hashlib.sha256(evidence.read_bytes()).hexdigest())

def collect(paths, identity):
    rows={r['id']:dict(status='INCOMPLETE',reason='missing row evidence') for r in matrix()}
    for path in paths:
        try:
            data=json.loads(path.read_text());key=data['row']
            if key not in rows:raise ValueError('unknown support row '+key)
            if rows[key].get('path'):raise ValueError('duplicate support row '+key)
            ok=(data.get('schema_version')==1 and data.get('status')=='PASS'
                and data.get('source_id')==identity
                and data.get('acceptance_policy')==acceptance_policy())
            rows[key]=dict(status='PASS' if ok else 'FAIL',path=str(path))
        except (OSError,ValueError,KeyError,TypeError,AttributeError) as e:
            return dict(status='FAIL',error=str(e),rows=rows,source_id=identity)
    status='FAIL' if any(v['status']=='FAIL' for v in rows.values()) else ('PASS' if all(v['status']=='PASS' for v in rows.values()) else 'INCOMPLETE')
    return dict(schema_version=1,status=status,source_id=identity,rows=rows,
                acceptance_policy=acceptance_policy())

def main():
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='action',required=True)
    sub.add_parser('matrix')
    r=sub.add_parser('receipt');r.add_argument('--row',required=True);r.add_argument('--evidence',type=Path,required=True);r.add_argument('--output',type=Path,required=True)
    c=sub.add_parser('collect');c.add_argument('--directory',type=Path,required=True);c.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    if a.action=='matrix':print(json.dumps({'include':matrix()},separators=(',',':')));return 0
    if a.action=='receipt':report=receipt(a.row,a.evidence)
    else:report=collect(sorted(a.directory.rglob('platform-result.json')),source_identity())
    a.output.write_text(json.dumps(report,indent=2)+'\n');print(report['status'],a.output)
    return 0 if report['status']=='PASS' else 1
if __name__=='__main__':raise SystemExit(main())
