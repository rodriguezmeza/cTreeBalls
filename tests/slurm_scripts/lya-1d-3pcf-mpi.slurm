#!/bin/bash
#SBATCH --job-name=lya-1d-3pcf-mpi
#SBATCH --constraint=cpu
#SBATCH --qos=regular
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=128
#SBATCH --time=01:00:00
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
exec "${CBALLS_PYTHON:-python3}" "$CTREEBALLS_ROOT/tests/python/run_job.py" lya-1d-3pcf-mpi "$@"
