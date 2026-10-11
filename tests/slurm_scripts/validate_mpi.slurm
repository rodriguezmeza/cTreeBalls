#!/bin/bash
#SBATCH --job-name=cballs-mpi-validation
#SBATCH --constraint=cpu
#SBATCH --qos=regular
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=00:30:00
#SBATCH --output=cballs-%x-%j.out
#SBATCH --export=ALL

# Slurm copies this script to a spool directory. Resolve the checkout explicitly.
set -euo pipefail
export CTREEBALLS_ROOT="${CTREEBALLS_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
if [[ -n "${CBALLS_ENV_FILE:-}" ]]; then source "$CBALLS_ENV_FILE"; fi
if [[ ! -f "$CTREEBALLS_ROOT/tests/python/run_job.py" ]]; then
    echo 'Set CTREEBALLS_ROOT to the shared cTreeBalls checkout.' >&2
    exit 2
fi
cd "$CTREEBALLS_ROOT"
export MPICH_GPU_SUPPORT_ENABLED=0
export OMP_NUM_THREADS=2 OMP_PLACES=cores OMP_PROC_BIND=spread OMP_DYNAMIC=FALSE
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg
export PYTHONPATH="$CTREEBALLS_ROOT:$CTREEBALLS_ROOT/tests/python:$CTREEBALLS_ROOT/tests/python:$CTREEBALLS_ROOT/tests/python:$CTREEBALLS_ROOT/scripts"
export TMPDIR="${CBALLS_OUTPUT_ROOT:-${SCRATCH:?Set SCRATCH or CBALLS_OUTPUT_ROOT}/ctreeballs-jobs}/validation-${SLURM_JOB_ID}"
mkdir -p "$TMPDIR"
py="${CBALLS_PYTHON:-python3}"
launch=(srun --export=ALL --kill-on-bad-exit=1 --cpu-bind=cores --nodes=2 --ntasks=2 --ntasks-per-node=1 --cpus-per-task=4)
"${launch[@]}" "$py" tests/python/check_mpi_environment.py
"${launch[@]}" "$py" tests/python/test_shear_pivot_reuse_mpi.py
"${launch[@]}" "$py" tests/python/test_scalar_pivot_reuse.py --mpi
"$py" tests/python/test_lya_forest_mpi.py --mpi-only --cython --mpi-command "${launch[*]}"
"$py" tests/python/test_lya_hierarchy.py --mpi-command "${launch[*]}"
"$py" tests/python/test_lya_hierarchy.py --los-tree --mpi-command "${launch[*]}"
"$py" tests/python/test_octree_3pcf_3d_mpi.py --mpi-only --cython --mpi-command "${launch[*]}"
echo 'PASS: cross-node MPI numerical and recovery checks'
