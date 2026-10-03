# kdtree-shear-sphere-2balls-mpi

MPI+OpenMP counterpart of `kdtree-shear-sphere-2balls-omp`, enabled in the
active Makefile profile. Both 2PCF and 3PCF use the optimized shared shear kernel.

```sh
make KDTREESHEARSPHERE2BALLSMPION=1 cballs
mpiexec -n 4 ./cballs parameters.ini \
    searchMethod=kdtree-shear-sphere-2balls-mpi numberThreads=2
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
death and a failed MPI communicator are outside this recovery contract.

`theta`, `rsmooth`, `nsmooth`, `SMOOTHPIVOTON`, `BALLS4SCANLEVON` and `THETA`
retain the shared OpenMP kernel's meanings. Approximation requires qualification
against an exact reference for each observable and data geometry. The
experimental `shear-pivot-reuse` option and `legacy-one-ball` compatibility
option are rejected in these MPI engines. MPI reduction ordering can change
roundoff across rank/thread counts; bitwise identity is not promised.

Set `CBALLS_SHEAR_PROFILE=1` for rank/thread timings. Regression coverage:

```sh
python tests/python/test_shear_sphere_mpi.py \
    --output /tmp/shear-mpi-check --mpi-command 'mpiexec'
```

This compares 1/2/4 ranks and 1/2 threads against OpenMP and independent shear
oracles, covers both orders and combined runs, masks, cross catalogs, repeated
roles, poles, sparse work, smoothing and approximation. It
checks raw histograms at strict floating-point tolerance. Corrected bins with
window condition number above 1e8 are marked unstable in the retained evidence
and checked by backward residual; they are not forward-accuracy qualifications.
The release gate also faults every reached named application boundary on either rank under owned
and borrowed MPI, followed by recovery.
