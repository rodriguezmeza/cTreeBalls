#!/bin/bash
#SBATCH --job-name=ctreeballs-acceptance-benchmark
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --exclusive
#SBATCH --time=02:00:00
# Supply account/partition/cpus-per-task to sbatch for the target site.
set -euo pipefail
if [ "$#" -ne 3 ]; then
    echo 'Usage: sbatch --cpus-per-task=N benchmark_dedicated.slurm RELEASE_JSON ACCURACY_JSON NEW_OUTPUT_DIRECTORY' >&2
    exit 2
fi
: "${SLURM_JOB_ID:?Submit this script with sbatch}"
: "${SLURM_CPUS_PER_TASK:?Set --cpus-per-task}"
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export OMP_PLACES=cores
export OMP_PROC_BIND=close
# Run in the accepted checkout, with its matching Python extension already built.
"${CBALLS_PYTHON:-python3}" scripts/benchmark_dedicated.py \
    --release "$1" --accuracy "$2" --output "$3" \
    --threads "$SLURM_CPUS_PER_TASK" --repeats 7 --warmups 1 \
    --exclusive-node --allocation-id "slurm:$SLURM_JOB_ID"
