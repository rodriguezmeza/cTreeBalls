# balltree-shear-sphere-2balls-mpi

MPI+OpenMP counterpart of `balltree-shear-sphere-2balls-omp`, enabled in the
active Makefile profile. Both 2PCF and 3PCF use the optimized shared shear kernel.

```sh
make BALLTREESHEARSPHERE2BALLSMPION=1 cballs
# Set searchMethod = balltree-shear-sphere-2balls-mpi and numberThreads = 2 in parameters.ini.
mpiexec -n 4 ./cballs paramfile=parameters.ini
```

Use `options=only-2pcf,no-smooth-pivot` for pairs,
`options=only-3pcf,no-smooth-pivot` for triplets, or
`options=no-smooth-pivot` for both. Preserve other required options in your
parameter file. For exact unsmoothed results use `theta=0`.

With `BALLS4SCANLEVON=1`, work-estimated pivot tasks are assigned largest first
to the least-loaded rank; OpenMP dynamically schedules that rank's tasks.
Catalogs and trees are replicated on every rank; task histograms are local.
This distributes computation, not catalog storage. Use the same thread count
and settings on every rank, and load the same catalog on every rank.

The pair-only route uses the symmetric dual-tree frontier. Exact combined
runs (no-one-ball, no-two-balls, or theta=0) fuse pairs with the 3PCF ring walk.
Approximate unsmoothed combined runs retain independent pair acceptance;
smoothed combined runs share transport in the pivot walk. The default 3PCF
uses exact body pivots with optional accepted neighbor cells. Small-log-bin
lookup, spin-2 transport and signed-mode products share the OpenMP optimizations.
Binary trees also use direct great-circle pair projection.

Raw histograms are packed into a real reduction and counters into an integer
reduction. Only rank zero normalizes, solves mode coupling, reconstructs angular
Gamma, writes products, and publishes `getResults()`. Invoke Python `Run()`
collectively on all ranks; read global results on rank zero. Tree preparation,
smoothing, result and traversal allocations have collective error boundaries;
returned application errors permit cleanup and same-object recovery. Process
death and a failed MPI communicator are outside this recovery guarantees.

`theta`, `rsmooth`, `nsmooth`, `SMOOTHPIVOTON`, `BALLS4SCANLEVON` and `THETA`
retain the shared OpenMP kernel's meanings. Approximation requires qualification
against an exact reference for each observable and data geometry. The
`shear-pivot-reuse` option is supported when `BALLS4SCANLEVON=1`, with
consistent controls across ranks; the MPI engine rejects `legacy-one-ball`. MPI reduction ordering can change
roundoff across rank/thread counts; bitwise identity is not promised.

Set `CBALLS_SHEAR_PROFILE=1` for rank/thread timings. Regression coverage:

```sh
python tests/python/test_shear_sphere_mpi.py \
    --output /tmp/shear-mpi-check --mpi-command 'mpiexec'
```

This compares 1/2/4 ranks and 1/2 threads against OpenMP and independent shear
reference calculations, covers both orders and combined runs, masks, cross catalogs, repeated
roles, poles, sparse work, smoothing and approximation. It
checks raw histograms at strict floating-point tolerance. Corrected bins with
window condition number above 1e8 are marked unstable in the retained evidence
and checked by backward residual; they are not forward-accuracy qualifications.
The release checks also faults every reached named application boundary on either rank under owned
and borrowed MPI, followed by recovery.


Tree construction shares the OpenMP ball-tree builder on each rank. Catalogue
selection, pointer preparation, source frames, large sibling branches and
upper-node statistics/moments run in parallel. PCA projections are cached in
the member pass and reused by the unchanged median partition. Set
`numberThreads` for workers per rank; no catalogue or tree is distributed
between ranks during construction. Additional MPI ranks therefore replicate
this work and its memory. The projection cache is temporary and capped at
256 MiB per rank; `no-balltree-shear-projection-cache` forces the direct
calculation and `no-balltree-parallel-build` forces serial construction.
The parallel-build MPI regression includes odd, masked catalogues above the
construction threshold and leaf capacities 1, 8 and 31.

For Perlmutter CPU jobs, shared inputs, effective parameter files and scaling
checks, see [the batch-job guide](../../tests/slurm_scripts/README.md).
