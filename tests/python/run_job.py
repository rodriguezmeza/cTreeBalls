#!/usr/bin/env python3
"""Plan or run one registered engine in a Perlmutter CPU allocation."""
from __future__ import annotations
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[2]


def engines(root):
    return {e['name']: e for e in json.loads((root/'capabilities/engines.json').read_text())['engines']}


def integer_env(env, name, default):
    text = env.get(name, str(default))
    if not re.fullmatch(r'[0-9]+', text):
        raise ValueError(f'{name} must be an integer, got {text!r}')
    return int(text)


def layout(engine, env, nodes=None, ranks_per_node=None, threads=None):
    mpi = engine.endswith('-mpi')
    nodes = nodes if nodes is not None else integer_env(env, 'SLURM_JOB_NUM_NODES', 2 if mpi else 1)
    rpn = ranks_per_node if ranks_per_node is not None else integer_env(env, 'SLURM_NTASKS_PER_NODE', 2 if mpi else 1)
    cpus = integer_env(env, 'SLURM_CPUS_PER_TASK', 256//max(rpn,1))
    threads = threads if threads is not None else cpus//2
    if min(nodes,rpn,threads) < 1 or rpn*threads > 128:
        raise ValueError('require positive counts and at most 128 physical cores per CPU node')
    if not mpi and (nodes != 1 or rpn != 1):
        raise ValueError('OpenMP engines require one node and one rank')
    if 'SLURM_JOB_ID' in env:
        if nodes > integer_env(env,'SLURM_JOB_NUM_NODES',nodes) or 2*threads > cpus:
            raise ValueError('requested placement exceeds this Slurm allocation')
        if nodes*rpn > integer_env(env,'SLURM_NTASKS',nodes*rpn):
            raise ValueError('requested MPI ranks exceed this Slurm allocation')
        if rpn > integer_env(env,'SLURM_NTASKS_PER_NODE',rpn):
            raise ValueError('requested ranks per node exceed this Slurm allocation')
    return dict(nodes=nodes,ranks_per_node=rpn,ranks=nodes*rpn,threads_per_rank=threads,
                logical_cpus_per_rank=2*threads,physical_cores=nodes*rpn*threads)


def plan(args, env):
    root = args.root.resolve()
    table = engines(root)
    if args.engine not in table:
        raise ValueError('unknown registered engine: '+args.engine)
    shape = layout(args.engine,env,args.nodes,args.ranks_per_node,args.threads)
    native = args.mode == 'native'
    shear = '-shear-' in args.engine
    forest = args.engine.startswith('lya-')
    if args.mode == 'shear' and not shear:
        raise ValueError('shear mode requires a shear engine')
    if args.mode == 'kappa' and (shear or forest or '-3pcf-3d-' in args.engine):
        raise ValueError('kappa mode requires an angular scalar engine')
    if args.mode == 'lya' and not (forest or '-3pcf-3d-' in args.engine):
        raise ValueError('lya mode requires a forest or physical-3D engine')
    if native and args.params is None:
        raise ValueError('native mode needs --params (or CBALLS_PARAMS) with file-backed input')
    extra = list(args.arguments)
    if extra[:1] == ['--']: extra.pop(0)
    if native and any(x.split('=',1)[0] in {'search','searchMethod','numberThreads','rootDir','paramfile'} for x in extra):
        raise ValueError('engine, threads, parameter file and output are managed by this launcher')
    managed = {'--engine','--engines','--threads','--mpi-ranks','--mpiexec',
               '--mpi-extra-arg','--mpi-worker','--outdir','--output','--cballs',
               '--list-engines','--help','-h'}
    if not native and any(x.split('=',1)[0] in managed for x in extra):
        raise ValueError('driver engine/launcher/output arguments are managed by this launcher')
    if not native and not any(x.split('=',1)[0] in {'--fits','--catalog','--catalog-npz','--ascii','--synthetic-nbody','--synthetic-nside','--synthetic'} for x in extra):
        raise ValueError('supply an input catalogue after -- (synthetic input must be explicit)')
    executable = Path(env.get('CBALLS_EXECUTABLE', str(root/'cballs'))).resolve()
    python = env.get('CBALLS_PYTHON', sys.executable)
    base = Path(env.get('CBALLS_OUTPUT_ROOT', str(Path(env.get('SCRATCH', str(root/'results')))/'ctreeballs-jobs')))
    label = env.get('SLURM_JOB_ID','preview')+'-'+args.engine
    out = (args.outdir or base/(label if args.dry_run else label+'-'+str(time.time_ns()))).resolve()
    working = Path(env.get('CBALLS_WORKDIR',str(root))).resolve()
    runenv = dict(env,MPICH_GPU_SUPPORT_ENABLED='0',OMP_NUM_THREADS=str(shape['threads_per_rank']),OMP_DYNAMIC='FALSE',
                  OMP_PLACES='cores',OMP_PROC_BIND='spread',OPENBLAS_NUM_THREADS='1',
                  MKL_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1',PYTHONUNBUFFERED='1',
                  PYTHONDONTWRITEBYTECODE='1',PYTHONNOUSERSITE='1',MPLBACKEND='Agg',
                  PYTHONPATH=str(root)+os.pathsep+str(root/'tests/python')+os.pathsep+str(root/'tests/python'),
                  CTREEBALLS_SOURCE_ROOT=str(root),CTREEBALLS_CBALLS=str(executable),
                  TMPDIR=str(out/'shared-tmp'),
                  MPLCONFIGDIR=env.get('MPLCONFIGDIR',str(out/'python-cache/matplotlib')),
                  XDG_CACHE_HOME=env.get('XDG_CACHE_HOME',str(out/'python-cache')))
    # No automatic task-count tuning: preserve the native adaptive policy.
    # Explicit user overrides are exported identically by srun to all ranks.
    defaults = dict(CBALLS_SHEAR_PIVOT_TOL='3',CBALLS_SHEAR_BIN_THETA='0.5',
                    CBALLS_SCALAR_PIVOT_TOL='0.1',CBALLS_SCALAR_BIN_THETA='0',
                    CBALLS_SHEAR_TASK_MEMORY_MB='256')
    for k,v in defaults.items():runenv.setdefault(k,v)
    launch = [env.get('CBALLS_SRUN','srun'),'--export=ALL','--kill-on-bad-exit=1',
              '--cpu-bind=cores','--nodes='+str(shape['nodes']),
              '--ntasks-per-node='+str(shape['ranks_per_node']),
              '--cpus-per-task='+str(shape['logical_cpus_per_rank'])]
    effective_parameters = None
    if native:
        params=args.params.expanduser().resolve()
        if not params.is_file():raise ValueError('parameter file not found: '+str(params))
        text = params.read_text()
        # Parameter files take precedence over CLI defaults in the native reader.
        # Write a separate effective file rather than silently losing overrides.
        overrides = {}
        for item in extra:
            key, sep, value = item.partition('=')
            if not sep or not re.fullmatch(r'[A-Za-z][A-Za-z0-9_]*', key) or '\n' in value or '\r' in value:
                raise ValueError('native overrides must be single-line key=value arguments')
            overrides[key] = value
        overrides.update(searchMethod=args.engine, numberThreads=str(shape['threads_per_rank']),
                         rootDir=str(out/'native'))
        lines = [line for line in text.splitlines()
                 if not any(re.match(r'^\s*'+re.escape(key)+r'\s*=',line) for key in overrides)]
        effective_parameters = '\n'.join(lines+[f'{k} = {v}' for k,v in overrides.items()])+'\n'
        match = re.search(r'^\s*infile\s*=\s*([^\n#%]*)',effective_parameters,re.M)
        if not match or not match[1].strip() or re.search(r'(?:^|[,\s])m[0-9]+(?:[,\s]|$)',match[1]):
            raise ValueError('native mode requires file-backed infile; memory catalogues such as m0 cannot be replayed')
        command=launch+['--ntasks='+str(shape['ranks']),str(executable),
                        'paramfile='+str(out/'effective.parameters')]
    else:
        driver=root/'tests/python'/f'{args.mode}_corr_all_engines.py'
        common=[python,str(driver),'--engine',args.engine,'--threads',str(shape['threads_per_rank']),
                '--mpi-ranks',str(shape['ranks'])]
        if args.mode=='shear':
            # The shear parent writes a shared packet and invokes srun itself.
            common += ['--outdir',str(out/'products'),'--cballs',str(executable)]
            command=common+['--mpiexec',launch[0]]+['--mpi-extra-arg='+x for x in launch[1:]]+extra
            if not args.engine.endswith('-mpi'):
                command=launch+['--ntasks=1']+common+extra
        else:
            # Scalar/forest drivers are MPI-aware: enter them once per rank.
            common += ['--output' if args.mode=='lya' else '--outdir',str(out/'products')]
            if args.mode=='kappa':common += ['--cballs',str(executable)]
            command=launch+['--ntasks='+str(shape['ranks'])]+common+extra
    settings={k:v for k,v in runenv.items() if k.startswith(('OMP_','CBALLS_SHEAR_','CBALLS_SCALAR_','CBALLS_LYA_')) or k in ['MPICH_GPU_SUPPORT_ENABLED','TMPDIR','MPLCONFIGDIR','XDG_CACHE_HOME','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','LOADEDMODULES','CONDA_PREFIX','CTREEBALLS_SOURCE_ROOT']}
    return dict(engine=args.engine,mode=args.mode,layout=shape,command=command,working_directory=str(working),
                output=str(out),environment=settings,executable=str(executable),
                effective_parameters=effective_parameters,
                scientific_status='UNMEASURED unless this input is compared with an exact reference'),runenv


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('engine');p.add_argument('--root',type=Path,default=Path(os.environ.get('CTREEBALLS_ROOT',ROOT)))
    p.add_argument('--mode',choices=['native','shear','kappa','lya'],default='native')
    p.add_argument('--params',type=Path,default=os.environ.get('CBALLS_PARAMS'))
    p.add_argument('--outdir',type=Path);p.add_argument('--nodes',type=int)
    p.add_argument('--ranks-per-node',type=int);p.add_argument('--threads',type=int)
    p.add_argument('--dry-run',action='store_true')
    # Split explicitly: REMAINDER after a positional would swallow our options.
    words=list(sys.argv[1:] if argv is None else argv)
    tail=words[words.index('--')+1:] if '--' in words else []
    words=words[:words.index('--')] if '--' in words else words
    args=p.parse_args(words);args.arguments=tail
    try:record,env=plan(args,os.environ)
    except (ValueError,OSError) as exc:p.error(str(exc))
    if args.dry_run:
        print(json.dumps(record,indent=2));return 0
    if not os.environ.get('SLURM_JOB_ID'):p.error('run in a Slurm allocation, or use --dry-run')
    out=Path(record['output']);out.mkdir(parents=True,exist_ok=False)
    (out/'shared-tmp').mkdir()
    if args.params:
        shutil.copy2(args.params.expanduser().resolve(),out/'input.parameters')
        (out/'effective.parameters').write_text(record['effective_parameters'])
    (out/'job-settings.json').write_text(json.dumps(record,indent=2)+'\n')
    start=time.monotonic();code=1;error=None
    try:
        # cballs intentionally returns 1 for informational options.
        probe=subprocess.run([record['executable'],'options=print-search-methods'],
            cwd=record['working_directory'],env=env,stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,text=True,timeout=60)
        (out/'compiled-engines.txt').write_text(probe.stdout)
        if probe.returncode not in (0,1) or not re.search(r'^- '+re.escape(args.engine)+r' \(id=',probe.stdout,re.M):
            raise RuntimeError('requested engine is unavailable; see '+str(out/'compiled-engines.txt'))
        record['executable_sha256']=hashlib.sha256(Path(record['executable']).read_bytes()).hexdigest()
        record['command_shell_display']=shlex.join(record['command'])
        (out/'job-settings.json').write_text(json.dumps(record,indent=2)+'\n')
        print('Output:',out,flush=True);print(shlex.join(record['command']),flush=True)
        with (out/'run.log').open('w') as log, subprocess.Popen(record['command'],cwd=record['working_directory'],env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,bufsize=1) as proc:
            for line in proc.stdout:
                sys.stdout.write(line);sys.stdout.flush();log.write(line);log.flush()
            code=proc.wait()
        if args.mode=='native':
            used=list((out/'native').rglob('*-usedvalues'))
            if len(used)==1:shutil.copy2(used[0],out/'parameters-usedvalues')
    except (OSError,RuntimeError,subprocess.TimeoutExpired) as exc:
        error=str(exc);print(error,file=sys.stderr)
    finally:
        (out/'job-status.json').write_text(json.dumps(dict(exit_code=code,error=error,
            launcher_wall_seconds=time.monotonic()-start),indent=2)+'\n')
    return code

if __name__=='__main__':raise SystemExit(main())
