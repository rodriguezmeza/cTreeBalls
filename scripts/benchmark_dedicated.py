#!/usr/bin/env python3
"""Repeated, randomized, cold-process measurements after fixed correctness gates.

--smoke validates this workflow locally and never certifies node performance.
Dedicated runs require an explicit exclusive-allocation attestation, at least
five measured repetitions and a host-wide cooperative lock.
"""
import argparse
from contextlib import contextmanager
import fcntl
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import random
import statistics
import subprocess
import sys
import tempfile
import time
from acceptance_metadata import acceptance_policy, source_identity, ACCURACY_ATOL, EXACT_RTOL
from accuracy_cases import SCALAR, SHEAR
ROOT=Path(__file__).resolve().parents[1]

def validate_prerequisites(release, accuracy, build_id, identity):
    from acceptance_metadata import REQUIRED
    if release.get('status')!='PASS' or release.get('source_id')!=identity:
        raise ValueError('requires completed release evidence for the current source')
    if any(release.get('components',{}).get(k,{}).get('status')!='PASS' for k in REQUIRED):
        raise ValueError('release components are incomplete')
    for evidence in (release,accuracy):
        if evidence.get('acceptance_policy')!=acceptance_policy():
            raise ValueError('correctness evidence uses a different acceptance policy')
    if accuracy.get('status')!='PASS' or accuracy.get('build',{}).get('id')!=build_id:
        raise ValueError('bounded-accuracy evidence must pass on the measured build')
    if not accuracy.get('cases') or not accuracy.get('exact_cross_engine'):
        raise ValueError('bounded-accuracy evidence is incomplete')
    if any(c.get('passed') is not True for c in accuracy['exact_cross_engine']):
        raise ValueError('exact-control evidence contains a failure')
    if any(c.get('acceptance_status')!='PASS' for c in accuracy['cases'] if not c.get('diagnostic')):
        raise ValueError('bounded accuracy contains a rejected configuration')

def configurations(engines):
    cases={c['engine']:c for c in SCALAR+SHEAR if not c.get('compatibility') and not c.get('diagnostic')}
    if len(engines)!=len(set(engines)) or not engines or any(e not in cases for e in engines):
        raise ValueError('choose unique engines from the bounded accuracy policy')
    result=[]
    for engine in engines:
        for fixture in ('signed_sky','clustered_smooth'):
            modes=['exact','approx']
            if fixture=='clustered_smooth' and cases[engine]['smoothing']:modes+=['smooth']
            result += [(engine,fixture,mode) for mode in modes]
    return result

def schedule(configs,repeats,seed):
    rows=[(*c,repeat) for repeat in range(repeats) for c in configs]
    random.Random(seed).shuffle(rows)
    return rows

def summary(values):
    if not values or any(not math.isfinite(v) or v<0 for v in values):raise ValueError('invalid timing')
    ordered=sorted(values);median=statistics.median(values)
    # Linear quantiles, defined also for the two-sample local smoke.
    def quantile(q):
        x=(len(ordered)-1)*q;i=int(x);j=min(i+1,len(ordered)-1)
        return ordered[i]+(x-i)*(ordered[j]-ordered[i])
    return dict(n=len(values),median=median,min=min(values),max=max(values),
                q25=quantile(.25),q75=quantile(.75),mad=statistics.median(abs(v-median) for v in values))

@contextmanager
def node_lock(path):
    # Keep this inode in place after close: unlinking a held lock could let a
    # second process lock a different inode under the same pathname.
    with path.open('a+') as stream:
        try:fcntl.flock(stream,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError as e:raise RuntimeError('another cTreeBalls benchmark holds the node lock') from e
        try:yield
        finally:fcntl.flock(stream,fcntl.LOCK_UN)

def run(args):
    import numpy as np
    sys.path.insert(0,str(ROOT));import cyballs
    if args.threads<1 or args.warmups<1 or args.repeats<(2 if args.smoke else 5):
        raise ValueError('need positive threads/warmups and >=2 smoke or >=5 dedicated repeats')
    if not args.smoke and (not args.exclusive_node or not args.allocation_id):
        raise ValueError('dedicated runs require --exclusive-node and --allocation-id')
    identity=source_identity();build=cyballs.build_info();policy=acceptance_policy()
    release=json.loads(args.release.read_text());accuracy=json.loads(args.accuracy.read_text())
    validate_prerequisites(release,accuracy,build['id'],identity)
    configs=configurations(args.engines)
    for engine,fixture,mode in configs:
        if mode!='exact' and not any(c['engine']==engine and c['fixture']==fixture and c['mode']==mode and c.get('acceptance_status')=='PASS' for c in accuracy['cases']):
            raise ValueError('measured configuration lacks prior correctness evidence')
    output=args.output.resolve();output.mkdir(parents=True,exist_ok=False)
    report=dict(schema_version=1,status='RUNNING',performance_status='NOT_MEASURED',
        source_id=identity,build=build,acceptance_policy=policy,
        methodology='cold process per sample; sequential workers; discarded warmups; shuffled measured configurations; fixed correctness limits',
        host=dict(system=platform.platform(),machine=platform.machine(),node=platform.node(),logical_cpus=os.cpu_count()),
        allocation=dict(exclusive_attested=args.exclusive_node,identifier=args.allocation_id,
            slurm_job=os.environ.get('SLURM_JOB_ID'),slurm_nodes=os.environ.get('SLURM_JOB_NUM_NODES')),
        threads=args.threads,repeats=args.repeats,warmups=args.warmups,seed=args.seed,
        inputs={str(p.resolve()):hashlib.sha256(p.read_bytes()).hexdigest() for p in (args.release,args.accuracy)},
        warmup_samples=[],samples=[],summary=[])
    def save():
        tmp=output/'benchmark.json.tmp';tmp.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n');tmp.replace(output/'benchmark.json')
    def sample(config,index,warmup):
        engine,fixture,mode=config;label=('warmup' if warmup else 'sample')+f'-{index:04d}'
        dest=output/label;log=output/(label+'.log')
        command=[sys.executable,ROOT/'scripts/benchmark_accuracy_and_timing.py','--worker','--module-root',ROOT,
                 '--engine',engine,'--fixture',fixture,'--mode',mode,'--threads',str(args.threads),'--output',dest]
        with log.open('w') as stream:
            done=subprocess.run(list(map(str,command)),cwd=ROOT,stdout=stream,stderr=subprocess.STDOUT,timeout=args.timeout)
        if done.returncode:raise RuntimeError('worker failed; inspect '+str(log))
        row=json.loads((dest/'sample.json').read_text());row.update(directory=label,sequence=index)
        if row['build']['id']!=build['id']:raise ValueError('worker build changed')
        (report['warmup_samples'] if warmup else report['samples']).append(row);save()
    save()
    try:
        # The lock prevents cooperating runs by this user; scheduler exclusivity
        # is an operator responsibility, not inferred from low load averages.
        with node_lock(Path(tempfile.gettempdir())/f'ctreeballs-benchmark-{os.getuid()}.lock'):
            for i,c in enumerate(configs*args.warmups):sample(c,i,True)
            for i,(*c,repeat) in enumerate(schedule(configs,args.repeats,args.seed)):
                sample(c,i,False);report['samples'][-1]['repeat']=repeat
            reference={}
            for row in report['samples']:
                if row['mode']=='exact':reference.setdefault((row['engine'],row['fixture']),row)
            for row in report['samples']:
                ref=reference[(row['engine'],row['fixture'])]
                if ref['fixture_sha256']!=row['fixture_sha256']:raise ValueError('fixture changed across repetitions')
                metrics={}
                with np.load(output/row['directory']/'products.npz') as values,np.load(output/ref['directory']/'products.npz') as targets:
                    if set(values.files)!=set(targets.files):raise ValueError('observable set changed')
                    for name in targets.files:
                        target=targets[name];value=values[name]
                        error=float(np.linalg.norm(value-target));norm=float(np.linalg.norm(target))
                        limit=ACCURACY_ATOL+(EXACT_RTOL if row['mode']=='exact' else row['rtol'])*norm
                        metrics[name]=dict(error_l2=error if math.isfinite(error) else None,limit_l2=limit,
                                           passed=math.isfinite(error) and error<=limit)
                row['accuracy']=metrics;row['status']='PASS' if all(v['passed'] for v in metrics.values()) else 'FAIL'
            for c in configs:
                rows=[v for v in report['samples'] if (v['engine'],v['fixture'],v['mode'])==c]
                vals=[next(p['wall_seconds'] for p in v['phases'] if p['name']=='main_loop') for v in rows]
                report['summary'].append(dict(engine=c[0],fixture=c[1],mode=c[2],main_loop_wall_seconds=summary(vals),
                    maximum_peak_rss_bytes=max(v['peak_rss_bytes'] for v in rows)))
            if source_identity()!=identity:raise ValueError('source changed during measurement')
            report['status']='PASS' if all(v['status']=='PASS' for v in report['samples']) else 'FAIL'
            report['performance_status']='LOCAL_SMOKE_ONLY' if args.smoke else ('MEASURED' if report['status']=='PASS' else 'REJECTED')
    except Exception as error:report.update(status='FAIL',error=str(error))
    finally:save()
    print(report['status'],report['performance_status'],output/'benchmark.json')
    return 0 if report['status']=='PASS' else 1

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('release','accuracy','output'):p.add_argument('--'+key,required=True,type=Path)
    p.add_argument('--engines',nargs='+',default=['kdtree-2balls-omp','balltree-2balls-omp','octree-2balls-omp'])
    p.add_argument('--repeats',type=int,default=7);p.add_argument('--warmups',type=int,default=1)
    p.add_argument('--threads',type=int,default=2);p.add_argument('--seed',type=int,default=20261011)
    p.add_argument('--timeout',type=int,default=600);p.add_argument('--smoke',action='store_true')
    p.add_argument('--exclusive-node',action='store_true');p.add_argument('--allocation-id')
    return run(p.parse_args())
if __name__=='__main__':raise SystemExit(main())
