#!/usr/bin/env python3
"""Cold-process phase timings, peak RSS and retained exact/approximate accuracy.

Each measured sample runs alone in a fresh process. Accuracy limits come from
accuracy_acceptance.py, never from observed errors. No speedup is accepted
without the products, memory record, configuration and source identity.
"""
import argparse
from contextlib import contextmanager
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import resource
import statistics
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parents[1]


def peak_bytes():
    value=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if sys.platform=='darwin':return int(value)
    if sys.platform.startswith('linux'):return int(value*1024)
    raise RuntimeError('peak RSS unit must be defined for this platform')


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def worker(args):
    # Use the same workload implementation for before/after builds, changing
    # only the selected extension. Its actual identity is recorded below.
    sys.path[:0]=[str(args.module_root),str(ROOT/'scripts')]
    import numpy as np
    import accuracy_acceptance as accuracy
    import cyballs
    case=next(c for c in accuracy.SCALAR+accuracy.SHEAR if c['engine']==args.engine and not c.get('diagnostic'))
    if cyballs.search_method_id(args.engine)<0:raise RuntimeError('requested engine is unavailable')
    args.output.mkdir(parents=True,exist_ok=False)
    phases=[]
    @contextmanager
    def phase(name):
        rss0=peak_bytes();wall=time.perf_counter();cpu=time.process_time()
        try:yield
        finally:phases.append(dict(name=name,wall_seconds=time.perf_counter()-wall,
                                  process_cpu_seconds=time.process_time()-cpu,
                                  process_high_water_before_bytes=rss0,
                                  process_high_water_after_bytes=peak_bytes()))
    data=accuracy.fixture(args.fixture)
    np.savez(args.output/'fixture.npz',positions=data[0],kappa=data[1],weights=data[2],gamma=data[3])
    imported_peak=peak_bytes()
    values,metadata,timings=accuracy.run(case,data,args.mode,args.output/'run',phase=phase,threads=args.threads)
    workload_peak=peak_bytes()
    np.savez_compressed(args.output/'products.npz',**values)
    record=dict(engine=args.engine,fixture=args.fixture,mode=args.mode,threads=args.threads,
                phases=phases,peak_rss_bytes=workload_peak,import_and_fixture_peak_rss_bytes=imported_peak,
                peak_rss_semantics='process high-water mark including interpreter, imports and workload; phase snapshots are cumulative, not independent phase peaks',
                metadata=metadata,timings=timings,extension=str(Path(cyballs.__file__).resolve()),
                extension_sha256=sha(Path(cyballs.__file__)),build=cyballs.build_info(),
                fixture_sha256=sha(args.output/'fixture.npz'),products_sha256=sha(args.output/'products.npz'),
                rtol=case['rtol'],atol=1e-10)
    (args.output/'sample.json').write_text(json.dumps(record,indent=2,allow_nan=False)+'\n')


def run(args):
    import numpy as np
    # The manifest drives public engine eligibility, and the accepted accuracy
    # workload limits the benchmark settings that may claim a passing result.
    sys.path.insert(0,str(ROOT/'scripts'))
    from capabilities_generated import ENGINES
    args.output.mkdir(parents=True,exist_ok=False)
    report=dict(schema_version=1,status='RUNNING',host=dict(system=platform.platform(),machine=platform.machine(),
                python=sys.version,logical_cpus=os.cpu_count()),repeats=args.repeats,
                methodology='fresh process per sample; one process at a time; cold tree cache; 512-point retained exact subsamples; no speedup threshold',
                samples=[],summary=[])
    def save(): (args.output/'benchmark.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    save()
    try:
        for engine in args.engines:
            assert ENGINES[engine]['gate'].get('oracle'),engine
            for fixture in ('signed_sky','clustered_smooth'):
                exact=None;fixture_digest=None
                modes=['exact','approx']
                if fixture=='clustered_smooth' and engine!='octree-2balls-omp':modes+=['smooth']
                for mode in modes:
                    for repeat in range(args.repeats):
                        label=f'{engine}-{fixture}-{mode}-{repeat}'
                        directory=args.output/label;log=args.output/(label+'.log')
                        command=[sys.executable,Path(__file__).resolve(),'--worker','--engine',engine,'--fixture',fixture,
                                 '--mode',mode,'--threads',str(args.threads),'--module-root',args.module_root,'--output',directory]
                        with log.open('w') as stream:
                            completed=subprocess.run([str(v) for v in command],stdout=stream,stderr=subprocess.STDOUT,
                                                     cwd=args.module_root,timeout=240)
                        if completed.returncode:raise RuntimeError(label+' failed; inspect '+str(log))
                        row=json.loads((directory/'sample.json').read_text());arrays=np.load(directory/'products.npz')
                        if exact is None:
                            exact={k:arrays[k].copy() for k in arrays.files};fixture_digest=row['fixture_sha256']
                        assert row['fixture_sha256']==fixture_digest,'changed fixture across samples'
                        row.update(directory=label,repeat=repeat,accuracy={})
                        row['native_phases_seconds']=[{key:float(value) for key,value in re.findall(r'(\w+)\s*=\s*([0-9.eE+-]+)',line)}
                                                     for line in log.read_text().splitlines() if 'phase-timers:' in line]
                        row['native_phase_semantics']='*_wall are elapsed phase seconds; *_thread sum instrumented per-thread work and must not be added to elapsed wall time'
                        for key,target in exact.items():
                            error=float(np.linalg.norm(arrays[key]-target));norm=float(np.linalg.norm(target))
                            limit=(1e-11 if mode=='exact' else row['rtol'])*norm+row['atol']
                            row['accuracy'][key]=dict(absolute_l2=error,reference_l2=norm,relative_l2=error/max(norm,row['atol']),
                                                       limit_l2=limit,passed=bool(np.isfinite(error) and error<=limit))
                        row['status']='PASS' if all(v['passed'] for v in row['accuracy'].values()) else 'FAIL'
                        report['samples'].append(row);save()
                    rows=[r for r in report['samples'] if (r['engine'],r['fixture'],r['mode'])==(engine,fixture,mode)]
                    timing=[next(p['wall_seconds'] for p in r['phases'] if p['name']=='main_loop') for r in rows]
                    report['summary'].append(dict(engine=engine,fixture=fixture,mode=mode,samples=len(rows),
                        main_loop_wall_median_seconds=statistics.median(timing),main_loop_wall_min_seconds=min(timing),
                        main_loop_wall_max_seconds=max(timing),peak_rss_max_bytes=max(r['peak_rss_bytes'] for r in rows),
                        maximum_relative_l2={k:max(r['accuracy'][k]['relative_l2'] for r in rows) for k in exact},
                        status='PASS' if all(r['status']=='PASS' for r in rows) else 'FAIL'))
                    save()
        report['status']='PASS' if report['samples'] and all(r['status']=='PASS' for r in report['samples']) else 'FAIL'
    except Exception as error:
        report.update(status='FAIL',error=str(error));raise
    finally:save()
    print(report['status']+': '+str(args.output/'benchmark.json'))
    if report['status']!='PASS':raise SystemExit(1)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',required=True,type=Path)
    parser.add_argument('--module-root',type=Path,default=ROOT)
    parser.add_argument('--engines',nargs='+',default=['kdtree-2balls-omp','balltree-2balls-omp','octree-2balls-omp'])
    parser.add_argument('--repeats',type=int,default=3)
    parser.add_argument('--threads',type=int,default=2)
    parser.add_argument('--worker',action='store_true',help=argparse.SUPPRESS)
    parser.add_argument('--engine',help=argparse.SUPPRESS)
    parser.add_argument('--fixture',choices=('signed_sky','clustered_smooth'),help=argparse.SUPPRESS)
    parser.add_argument('--mode',choices=('exact','approx','smooth'),help=argparse.SUPPRESS)
    args=parser.parse_args();args.output=args.output.resolve();args.module_root=args.module_root.resolve()
    if args.repeats<1 or args.threads<1:parser.error('repeats and threads must be positive')
    (worker if args.worker else run)(args)


if __name__=='__main__': main()
