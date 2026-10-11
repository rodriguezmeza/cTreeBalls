#!/usr/bin/env python3
"""Validate the Python MPI environment inside the two-node validation job."""
import json
from pathlib import Path
import socket
from mpi4py import MPI
import cyballs

ROOT=Path(__file__).resolve().parents[2]
comm=MPI.COMM_WORLD
expected=json.loads((ROOT/'capabilities/engines.json').read_text())['engines']
row=dict(rank=comm.rank,host=socket.gethostname(),mpi=MPI.Get_library_version(),
         build_id=cyballs.build_info()['id'],extension=cyballs.__file__,
         missing=[e['name'] for e in expected if cyballs.search_method_id(e['name'])<0])
rows=comm.allgather(row)
errors=[]
if comm.size!=2:errors.append('validation requires exactly two ranks')
if len({r['host'] for r in rows})!=comm.size:errors.append('validation requires different compute nodes')
if len({r['build_id'] for r in rows})!=1:errors.append('ranks imported different cTreeBalls builds')
if len({r['mpi'] for r in rows})!=1:errors.append('ranks imported different MPI libraries')
if any(r['missing'] for r in rows):errors.append('the validation profile needs every registered engine')
if comm.rank==0:print(json.dumps(dict(ranks=rows,errors=errors),indent=2),flush=True)
if errors:raise SystemExit(1)
