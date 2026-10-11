#!/usr/bin/env python3
"""Print a CPU strong-scaling sweep; submit only with explicit --submit."""
import argparse
import json
from pathlib import Path
import shlex
import subprocess
import sys

HERE=Path(__file__).resolve().parents[1]/'slurm_scripts'

def commands(engine,account,nodes,rpn,repeats,tail,qos,walltime):
    if not engine.endswith('-mpi') or not (HERE/(engine+'.slurm')).is_file():
        raise ValueError('select a registered MPI engine')
    if repeats<1 or any(n<1 for n in nodes):raise ValueError('positive nodes/repeats required')
    for n in nodes:
        for r in rpn:
            if r<1 or 128%r:raise ValueError('ranks per node must divide 128')
            for _ in range(repeats):
                yield ['sbatch','--parsable','--account='+account,'--qos='+qos,
                       '--time='+walltime,'--nodes='+str(n),'--ntasks-per-node='+str(r),
                       '--cpus-per-task='+str(256//r),str(HERE/(engine+'.slurm')),*tail]

def main(argv=None):
    words=list(sys.argv[1:] if argv is None else argv)
    cut=words.index('--') if '--' in words else len(words)
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('engine');p.add_argument('--account',required=True)
    p.add_argument('--nodes',nargs='+',type=int,default=[1,2,4])
    p.add_argument('--ranks-per-node',nargs='+',type=int,default=[2,4])
    p.add_argument('--repeats',type=int,default=3);p.add_argument('--qos',default='regular')
    p.add_argument('--time',default='01:00:00');p.add_argument('--submit',action='store_true')
    a=p.parse_args(words[:cut]);tail=words[cut+1:]
    if not tail:p.error('pass job arguments after --, including an input catalogue or parameter file')
    if any(x.split('=',1)[0] in {'--nodes','--ranks-per-node','--threads','--outdir'} for x in tail[:tail.index('--') if '--' in tail else len(tail)]):
        p.error('let the sweep choose placement and unique output directories')
    try:jobs=list(commands(a.engine,a.account,a.nodes,a.ranks_per_node,a.repeats,tail,a.qos,a.time))
    except ValueError as exc:p.error(str(exc))
    for command in jobs:
        print(shlex.join(command),flush=True)
        if a.submit:
            result=subprocess.run(command,check=True,capture_output=True,text=True)
            print(json.dumps({'job_id':result.stdout.strip(),'command':command}),flush=True)

if __name__=='__main__':main()
