#!/usr/bin/env python3
"""Select regression owners from changed paths and transitive C includes.

Unknown paths select the complete public matrix. Shared top-level contracts
also select it. Selection is conservative; the release gate always stays full.
"""
import argparse
from collections import defaultdict
from fnmatch import fnmatch
import hashlib
import json
from pathlib import Path
import re
import subprocess

from capabilities_generated import CAPABILITIES, ENGINES, expected_registry, gate_plan
from build_fingerprint import source_files
ROOT = Path(__file__).resolve().parents[1]


def include_dependents(root=ROOT):
    root=Path(root).resolve()
    paths=[p for p in source_files(root) if p.suffix in ('.c','.h','.cpp')]
    by_name=defaultdict(list)
    for p in paths: by_name[p.name].append(str(p.relative_to(root)))
    reverse=defaultdict(set)
    for p in paths:
        owner=str(p.relative_to(root))
        for included in re.findall(r'^\s*#\s*include\s+"([^"]+)"',p.read_text(errors='replace'),re.M):
            local=(p.parent/included).resolve()
            targets=[str(local.relative_to(root))] if local.is_file() and local.is_relative_to(root) else by_name[Path(included).name]
            for target in targets: reverse[target].add(owner)
    return reverse


def select(changed, registry=None, root=ROOT):
    registry=registry if registry is not None else {n:e['id'] for n,e in ENGINES.items() if e['gate'].get('oracle')}
    gate_plan(registry)  # A newly enabled engine must have an executable oracle.
    reverse=include_dependents(root)
    pending=list(dict.fromkeys(changed)); closure=set(pending)
    while pending:
        for parent in reverse.get(pending.pop(), ()):
            if parent not in closure: closure.add(parent); pending.append(parent)
    selected=set(); reasons={}
    for path in sorted(closure):
        shared=any(fnmatch(path,pattern) for pattern in CAPABILITIES['shared_sources'])
        owners={name for name in registry if any(fnmatch(path,p) for p in ENGINES[name]['sources'])}
        if shared or not owners:
            owners=set(registry)
            reasons[path]='shared contract' if shared else 'unknown ownership: full matrix'
        else: reasons[path]='declared owners'
        selected.update(owners)
    alternate=[]
    for name,e in ENGINES.items():
        if e['gate'].get('oracle') or not e['gate'].get('test_profile'):continue
        if any(any(fnmatch(path,p) for p in CAPABILITIES['shared_sources']+e['sources']) or reasons[path].startswith('unknown') for path in closure):
            alternate.append(dict(engine=name,tests=e['gate']['tests'],test_profile=e['gate']['test_profile']))
    cases=gate_plan({name:registry[name] for name in sorted(selected)})
    return dict(schema_version=1,changed=sorted(set(changed)),include_closure=sorted(closure),
                reasons=reasons,engines=sorted(selected),cases=cases,alternate_regressions=alternate,
                tests=sorted({test for c in cases for test in c['tests']}),
                policy='conservative include closure; unknown/shared paths select all active engines')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--changed',nargs='*',default=[])
    parser.add_argument('--base',help='Git base ref; includes committed, staged and unstaged tracked changes against it')
    parser.add_argument('--baseline',type=Path,help='retained build-fingerprint.json (also works without Git)')
    parser.add_argument('--output',required=True,type=Path)
    args=parser.parse_args();changed=list(args.changed)
    if args.base:
        result=subprocess.run(['git','diff','--name-only','-z',args.base,'--'],cwd=ROOT,check=True,capture_output=True)
        changed.extend(p for p in result.stdout.decode().split('\0') if p)
    if args.baseline:
        old=json.loads(args.baseline.read_text())['source_files']
        new={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in source_files(ROOT)}
        changed.extend(p for p in set(old)|set(new) if old.get(p)!=new.get(p))
    if not (args.base or args.baseline or args.changed):changed=['Makefile']
    try:
        from cyballs import build_info
        registry=expected_registry(build_info()['resolved_settings'])
    except ImportError: registry=None
    report=select(changed,registry)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,indent=2)+'\n')
    print(f'Selected {len(report["engines"])} engines and {len(report["tests"])} reference modules: {args.output}')


if __name__=='__main__': main()
