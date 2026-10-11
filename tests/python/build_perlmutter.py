#!/usr/bin/env python3
"""Build all registered engines and install cyballs in the active environment."""
import argparse
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[2]

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--jobs',type=int,default=8)
    p.add_argument('--dry-run',action='store_true')
    a=p.parse_args()
    if a.jobs<1:p.error('--jobs must be positive')
    settings=json.loads((ROOT/'capabilities/build_profile.json').read_text())['enabled_settings'].copy()
    settings.update(CC='cc',MPICC='cc',PYTHON=sys.executable,SINGLEPON='0',
                    DEFDIMENSION='3',OPTFLAG='-O3 -fno-fast-math',OMPFLAG='-fopenmp',
                    GSLINTERNAL='0',CFITSIOLIBON='0',SLEEFON='0',THETA='1',
                    CBALLS_SHEAR_PIVOT_TOL='3',CBALLS_SHEAR_BIN_THETA='0.5',
                    CBALLS_SCALAR_PIVOT_TOL='0.1',CBALLS_SCALAR_BIN_THETA='0')
    command=['make','-j'+str(a.jobs),*[k+'='+v for k,v in sorted(settings.items())],'cballs','cyballs']
    print(shlex.join(command),flush=True)
    if a.dry_run:return
    if os.environ.get('PE_ENV')!='GNU' or not os.environ.get('CRAY_MPICH_DIR'):
        p.error('load PrgEnv-gnu and cray-mpich; use build_perlmutter.sh')
    # Fail before an expensive build if Python/scientific or external C dependencies are missing.
    for module in ('numpy','Cython','setuptools','scipy','matplotlib','astropy'):
        __import__(module)
    # Check mpi4py availability without initializing MPI on the login node.
    import importlib.util
    if importlib.util.find_spec('mpi4py') is None:p.error('activate an environment with NERSC-compatible mpi4py')
    subprocess.run(['gsl-config','--version'],check=True)
    subprocess.run(['pkg-config','--modversion','cfitsio'],check=True)
    subprocess.run(command,cwd=ROOT,check=True)

if __name__=='__main__':main()
