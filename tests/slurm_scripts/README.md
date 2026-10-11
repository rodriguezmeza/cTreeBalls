# Perlmutter CPU jobs for cTreeBalls

Reviewed against the current source on 2026-10-11. This directory contains one
job for each of the 41 registered engines: **18 MPI+OpenMP and 23 OpenMP**.
The OpenMP jobs provide a one-node comparison. There are no MPI implementations
of the remaining OpenMP-only methods; a job script cannot create one.
dual node is an optional external reference, not a registered cTreeBalls engine.
`engine_jobs.json` records names, aliases, geometries and default placement.

These scripts have local planning, failure-path and MPI numerical checks.
**Execution and scaling on Perlmutter remain unmeasured.** Nothing here submits
jobs automatically. Account names, catalogues and the Python environment belong
to the user; they are not hardcoded.

## What changed since the earlier MPI report

The report that MPI shear rejects `shear-pivot-reuse` is outdated. All three
MPI shear specializations now define `SHEAR_SPHERE_PIVOT_REUSE`, use the shared
hierarchical search, check approximation controls collectively, and reduce raw
histograms before root publication. The octree, k-d tree and ball tree also use
their parallel rank-local construction paths. This does not distribute a tree
across nodes: **each rank retains a catalogue and its own tree**.

The shared 3PCF scheduling policy now aims for at least 64 tasks and normally
8 times the number of ranks times the threads per rank. It applies a per-rank
task-histogram memory limit, then assigns tasks to ranks by estimated cost.
The resulting task count also depends on the tree: inspect `shear_frontier`
metadata in the Python output instead of assuming the target was reached.
The remaining constants of 256 in binary pair/neighbor frontier construction
are separate from this adaptive global 3PCF task target.

The scheduler is not unlimited. Histogram storage scales with radial bins and
multipoles, and one histogram is still required even if it exceeds the requested
budget. `CBALLS_SHEAR_TASK_MEMORY_MB=256` bounds task histograms per rank; it does
not bound total process memory. Catalogue/tree replicas, thread scratch and
other histograms are additional. Shear reuse scratch can use up to 64 MiB per
thread. Task creation and cost assignment are replicated; binary frontier
selection scans the growing frontier, and rank assignment scans the rank load
array. These preparation costs deserve profiling at large rank counts.
For example, a double-precision 3PCF with 20 radial bins, maximum order 7 and
all count-window modes needs 569,600 bytes per task histogram. A 256 MiB budget
allows 471 task histograms per rank, below the 512-task target for 64 threads.
Increasing `CBALLS_SHEAR_TASK_MEMORY_MB` can permit more tasks, if total rank
memory fits; it is not automatically faster because splitting also changes reuse
and scheduling overhead. With more radial bins this limit becomes tighter.

Assignment is static between ranks and dynamic among local OpenMP threads;
there is no cross-rank work stealing.

Scalar two-ball MPI engines also support optional `scalar-pivot-reuse`, with
separate scalar controls. Its benefit depends strongly on catalogue geometry;
previous DES-like scalar tests did not show a general speed improvement.
Ly-alpha MPI engines distribute body/pivot blocks, forest rows or radial tree
work, depending on the estimator. Their moments/hierarchies can accelerate
zero-slop calculations without intentional bin approximation. The physical-3D
multipole engine has its own pivot-block scheduler. These estimators and their
normalizations are not interchangeable.

Source entry points: `support/shear_kernel/shear_frontier_policy.h`,
`support/shear_kernel/search_octree_shear_omp.c`, `support/shear_binary/`,
`source/mpi_runtime.c`, `docs/SCALAR_HIERARCHICAL_REUSE.md`,
`docs/LYA_HIERARCHICAL_REUSE.md`, `docs/LYA_LOS_HIERARCHICAL_REUSE.md`,
`docs/LYA_MULTIPOLE_HIERARCHICAL_REUSE.md`, and the MPI addon READMEs.

## Approximation controls and numerical error

| Control | Job default / meaning | Accuracy implication |
|---|---|---|
| `theta` / `--tree-theta` | Compiled/driver default 1; parameter files retain their own value | Cell acceptance; not a relative error percentage |
| `CBALLS_SHEAR_PIVOT_TOL` | 3; phase budget in radians for shear hierarchical reuse | Larger values may accept more reuse; zero disables reuse, not every other approximation |
| `CBALLS_SHEAR_BIN_THETA` | 0.5; radial bin-position allowance | Independent of theta; does not overwrite it |
| `nsmooth` / `--nsmooth` | Compiled/driver default 8 | Leaf capacity where supported; with smoothing disabled this is not a request to average the catalogue |
| `CBALLS_SCALAR_PIVOT_TOL` | 0.1 radians | Only scalar hierarchical reuse; shear environment values do not control kappa |
| `CBALLS_SCALAR_BIN_THETA` | 0 | Conservative radial containment in scalar reuse |
| Ly-alpha slops | Existing engine defaults, zero for certified exact hierarchy paths | Nonzero angular/radial/polar slops permit bin approximations; they are not percentage error bounds |
| `CBALLS_SHEAR_TASKS` | Unset: adaptive | Positive value requests a target, still subject to memory/tree limits |
| `CBALLS_SHEAR_TASK_MEMORY_MB` | 256 per rank | Task-histogram budget, not total memory |

The launcher exports the same environment to all ranks, preserves explicit user
overrides, and saves these values. It does not add reuse options automatically.
Request `--more-options shear-pivot-reuse` or, for the scalar engines,
`--more-options scalar-pivot-reuse`. Shear reuse needs 3PCF, positive theta,
no pivot smoothing and compatible traversal options; a 2PCF-only job cannot
benefit from it. MPI shear rejects `legacy-one-ball`.

For an angular shear reference, use `--tree-theta 0 --exact-tree
--no-smooth-pivot --more-options no-one-ball,no-two-balls`, with the same sample,
weights, bins and multipoles. For scalar two-ball references use theta zero,
no smoothing and `no-two-balls`, without scalar reuse. Ly-alpha reference
controls depend on the selected estimator: use its zero-slop/reference-kernel
checks, not angular-shear flags. See the family documentation above.

For every observable X, report
`relative_L2 = ||X_approx - X_reference||_2 / ||X_reference||_2`.
Use matching finite, supported bins, report coverage, and handle a zero reference
norm explicitly. Compare raw shear numerators, count/window modes and corrected
multipoles separately; an ill-conditioned window can amplify a small raw error.
A single aggregate can hide a poor individual multipole. Use identical masks,
separation units/cutoffs, catalogue identities and order selection.

There is **no universal 5% guarantee** for theta=1, phase=3, bin allowance=0.5
and leaf capacity=8. No new DES error measurement is supplied by these jobs.
The existing driver `--accuracy-reference /shared/exact-job/products` can
compare a matching saved reference; add `--accuracy-report-only` if a failed
threshold should be reported without making the whole job fail. An exact small
sample checks that sample, not automatically the full survey.

MPI parity at a fixed policy and scientific approximation error are different
checks. Changing ranks/threads can change the adaptive frontier, reuse grouping
and floating-point reduction order, even with the same approximation values.
Use `CBALLS_SHEAR_TASKS` to request a consistent frontier when studying this,
verify actual metadata, and still compare numerical products.

## Prepare a shared checkout and environment

Use a shared path accessible to every compute node, for example `$SCRATCH/cTreeBalls`.
Set the checkout and an existing CPU allocation account:

```bash
export CTREEBALLS_ROOT="$SCRATCH/cTreeBalls"
export CPU_ACCOUNT=YOUR_CPU_ACCOUNT
cd "$CTREEBALLS_ROOT"
module load PrgEnv-gnu cray-mpich
# Activate your prepared Python environment here.
export CBALLS_PYTHON="$(command -v python3)"
```

Python needs NumPy, SciPy, Cython, setuptools, matplotlib, astropy and mpi4py;
HEALPix kappa inputs additionally need healpy. MPI tests need the test
requirements. Use NERSC-compatible mpi4py linked to Cray MPICH. Do not mix it
with an Open MPI Python installation. Check MPI imports in a compute allocation,
not as an MPI run on a login node. Follow
[NERSC parallel Python guidance](https://docs.nersc.gov/development/languages/python/parallel-python/).

The supplied build uses external GSL (`gsl-config`) and CFITSIO (`pkg-config
cfitsio`), matching the public source package. Make these development libraries
available in the active environment. Then:

```bash
bash tests/slurm_scripts/build_perlmutter.sh --jobs 8
```

This builds all current registered engines and **installs the matching cyballs
extension into the active Python environment**. It sets `CC=cc MPICC=cc`, GNU
OpenMP, double precision, `-O3 -fno-fast-math`, and disables optional SLEEF for
this initial baseline. It preserves the requested approximation defaults.
Use `python3 tests/python/build_perlmutter.py --dry-run` to inspect the
Make command first. Do not build concurrently with production jobs using the
same checkout. Native and Python builds must match.

The build helpers now understand Open MPI, MPICH and the Cray `cc` flag-query
interfaces. This fixes the former assumption that every MPI wrapper implements
`-show`. See [NERSC compiler guidance](https://docs.nersc.gov/development/compilers/)
and [HPE cc query options](https://cpe.ext.hpe.com/docs/latest/craype/cc.html).

Jobs inherit the submission environment. If activation must be repeated in the
batch shell, set `CBALLS_ENV_FILE` to an absolute path to your shell environment
setup file. All job wrappers source it before launching Python. The file should
load the same modules and environment used to build the code.

## Run the DES shear calculation

```bash
unset CBALLS_SHEAR_PROFILE
export CBALLS_SHEAR_PIVOT_TOL=3
export CBALLS_SHEAR_BIN_THETA=0.5
sbatch --account="$CPU_ACCOUNT" \
  tests/slurm_scripts/balltree-shear-sphere-2balls-mpi.slurm \
  --mode shear -- \
  --fits "$SCRATCH/catalogues/DESY3_Takahashi_r002_bin1_region1.fits" \
  --fits-format desy3 --des-shear-convention local-east-north \
  --geometry sphere --binning sofia-fig1 --statistics 3pcf \
  --min-sep 8 --max-sep 200 --sep-units arcmin \
  --multipoles 7 --phi-bins 32 --tree-theta 1 --nsmooth 8 \
  --no-smooth-pivot --more-options shear-pivot-reuse \
  --max-points 0 --verbose 1
```

Substitute the octree or kdtree job filename to compare the other MPI shear
engines. Substitute the corresponding `-omp.slurm` filename for its one-node
OpenMP baseline. Do not add an outer `srun` to the shear Python launcher: the
single controller writes a shared catalogue packet and launches MPI workers
with `srun` itself. `TMPDIR` is set inside the shared output directory so all
nodes can read that packet. The kappa and Ly-alpha drivers instead run once per
rank under `srun`. Native mode also uses one executable per rank.

For longer work override `sbatch --time=...` and, when relevant, the Python
driver's `--timeout` (seconds). `--verbose 1` enables the engine's live progress
counters. `CBALLS_SHEAR_PROFILE=1` adds diagnostic instrumentation; unset it for
production timings. Setting it to `0` does not disable a presence-based switch.

## Kappa, Ly-alpha and native parameter files

```bash
sbatch --account="$CPU_ACCOUNT" tests/slurm_scripts/kdtree-2balls-mpi.slurm \
  --mode kappa -- --catalog-npz "$SCRATCH/catalogues/kappa.npz" \
  --statistics 3pcf --tree-theta 1 --nsmooth 8 \
  --no-smooth-pivot --verbose 1

sbatch --account="$CPU_ACCOUNT" tests/slurm_scripts/lya-los-tree-2pcf-3pcf-mpi.slurm \
  --mode lya -- --catalog "$SCRATCH/catalogues/forest_pixels.npz" \
  --statistics both --verbose 1

sbatch --account="$CPU_ACCOUNT" tests/slurm_scripts/octree-3pcf-3d-mpi.slurm \
  --params "$SCRATCH/configs/physical-3d.params"
```

Input schemas and scientific bin controls remain those of the corresponding
`tests/python/*_corr_all_engines.py` driver. Set them for your science;
angular kappa, physical-3D multipoles and the several Ly-alpha estimators have
different geometries. Use the driver `--help` before submission.

Every engine accepts native mode, the default. Provide a file-backed native
`infile`; Python's `infile=m0` cannot be replayed by the executable. Extra native
arguments go after `--` as `key=value`. The launcher copies the original file,
writes `effective.parameters`, replaces engine/threads/output, and applies
explicit extra overrides there. Use canonical native parameter names, such as
`searchMethod` rather than `search`, in parameter files. Input paths are resolved
from `CBALLS_WORKDIR` (default: the checkout); absolute paths are preferable.
This avoids native parameter-file values silently overriding launch settings.

Preview any job without submitting:

```bash
python3 tests/python/run_job.py balltree-shear-sphere-2balls-mpi \
  --dry-run --mode shear -- --fits /shared/catalog.fits --statistics 3pcf
```

## Placement, scaling and saved outputs

Perlmutter CPU nodes have 128 physical cores and 256 logical CPUs. The default
MPI job uses **2 nodes, 2 ranks per node, 64 OpenMP threads per rank** and requests
128 logical CPUs per rank. OpenMP-only jobs use one node and 128 threads.
Both use `OMP_PLACES=cores`, `OMP_PROC_BIND=spread` and `srun --cpu-bind=cores`.
This follows the layout in [NERSC hybrid job examples](https://docs.nersc.gov/jobs/examples/)
and [CPU affinity guidance](https://docs.nersc.gov/jobs/affinity/).
The scripts request CPU nodes because these engines contain no GPU offload.
They set `MPICH_GPU_SUPPORT_ENABLED=0` for host-only MPI communication.
If your mpi4py was built with GPU-aware Cray MPICH, keep its required CUDA
runtime modules available even for CPU work, as described in the
[NERSC Python instructions](https://docs.nersc.gov/development/languages/python/using-python-perlmutter/).

Override placement with sbatch, not driver arguments, for example:
`--nodes=2 --ntasks-per-node=4 --cpus-per-task=64` gives 32 threads per rank.
The launcher rejects physical-core oversubscription and multi-rank OpenMP jobs.
More ranks replicate more memory and construction work; four or eight ranks per
node may improve locality or balance but can lose to two ranks on a large tree.

Print a sweep (no jobs are submitted by default):

```bash
python3 tests/python/submit_scaling.py \
  balltree-shear-sphere-2balls-mpi --account "$CPU_ACCOUNT" \
  --nodes 1 2 4 --ranks-per-node 2 4 --repeats 3 -- \
  --mode shear -- --fits /shared/catalog.fits --fits-format desy3 \
  --des-shear-convention local-east-north --geometry sphere \
  --statistics 3pcf --binning sofia-fig1 --multipoles 7 --phi-bins 32 \
  --min-sep 8 --max-sep 200 --sep-units arcmin --tree-theta 1 \
  --nsmooth 8 --no-smooth-pivot --more-options shear-pivot-reuse --verbose 1
```

Inspect the printed commands, then add `--submit` before the separator to submit
them. Compare medians of at least three fresh runs, using the same catalogue,
sample seed, bins, multipoles and approximation settings. Record total wall time,
search time, tree construction, MPI rank count, per-rank memory and L2 error.
Do not use resumed outputs as timing measurements. Strong-scaling speedup is
`T(one-node baseline)/T(N nodes)` and node efficiency is `speedup/N` when the
per-node resource count is fixed. Sum of process CPU time is not wall time.

Output defaults to `$SCRATCH/ctreeballs-jobs/JOBID-ENGINE-UNIQUE_ID`; set
`CBALLS_OUTPUT_ROOT` to another shared filesystem location if needed. Each job
saves `job-settings.json`, the exact argument vector, environment controls,
compiled-engine listing, executable SHA256, `run.log` and `job-status.json`.
Python driver products retain their normal summaries, plots and parameter files.
Shear/kappa products include their usual `parameters-usedvalues`. Native runs
retain the original/effective inputs and copy their used-values file to the job
root when one is present. Temporary MPI packets use `shared-tmp/` and successful
Python temporary-directory cleanup removes them. Failed jobs return failure to
Slurm; check `sacct` as well, since a time-limit kill may prevent final status
publication. For allocation memory/cost use, for example,
`sacct -j JOBID --format=JobID,State,Elapsed,AllocCPUS,TotalCPU,MaxRSS`.

Local script checks can be rerun with `make test-perlmutter-jobs`.

## Validate on Perlmutter before production scaling

```bash
sbatch --account="$CPU_ACCOUNT" tests/slurm_scripts/validate_mpi.slurm
```

This runs two MPI ranks on two different nodes, prints MPI-library/build/host
information, then runs shear reuse parity, rank-mismatch/failure recovery,
parallel tree construction, scalar reuse, Ly-alpha and physical-3D multipole numerical tests.
It does not calibrate a full DES catalogue or prove scaling. All nodes must
see the checkout, input files, output directory and Python environment.

## Engine jobs

Each name below is an individual `.slurm` file. The JSON inventory includes
aliases and estimator descriptions.

| MPI+OpenMP | OpenMP |
|---|---|
| `kdtree-2balls-mpi` | `octree-sincos-omp` |
| `balltree-2balls-mpi` | `kdtree-2balls-omp` |
| `octree-2balls-mpi` | `balltree-2balls-omp` |
| `octree-shear-sphere-2balls-mpi` | `octree-2balls-omp` |
| `kdtree-shear-sphere-2balls-mpi` | `octree-shear-sphere-2balls-omp` |
| `balltree-shear-sphere-2balls-mpi` | `kdtree-shear-sphere-2balls-omp` |
| `octree-3pcf-3d-mpi` | `balltree-shear-sphere-2balls-omp` |
| `lya-2pcf-mpi` | `kdtree-box-omp` |
| `lya-3pcf-mpi` | `neighbor-boxes-omp` |
| `lya-2pcf-3pcf-mpi` | `octree-3pcf-3d-omp` |
| `lya-1d-2pcf-mpi` | `lya-2pcf-omp` |
| `lya-1d-3pcf-mpi` | `lya-3pcf-omp` |
| `lya-1d-2pcf-3pcf-mpi` | `lya-2pcf-3pcf-omp` |
| `lya-1d-tree-2pcf-mpi` | `lya-los-tree-2pcf-omp` |
| `lya-1d-tree-3pcf-mpi` | `lya-los-tree-3pcf-omp` |
| `lya-los-tree-2pcf-mpi` | `lya-los-tree-2pcf-3pcf-omp` |
| `lya-los-tree-3pcf-mpi` | `lya-1d-2pcf-omp` |
| `lya-los-tree-2pcf-3pcf-mpi` | `lya-1d-3pcf-omp` |
|  | `lya-1d-2pcf-3pcf-omp` |
|  | `lya-1d-tree-2pcf-omp` |
|  | `lya-1d-tree-3pcf-omp` |
|  | `lya-1d-tree-same-los-2pcf-omp` |
|  | `lya-anisotropic-multipole-3pcf-omp` |
