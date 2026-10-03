# balltree-2balls-mpi

This addon is the deterministic MPI+OpenMP counterpart of
`balltree-2balls-omp`. It uses the same FCFC PCA ball tree, dual-node-style
dual-node 2PCF traversal, and the production body-pivot LogMultipole 3PCF
traversal. Neighbor moments and repeated-neighbor subtraction share the
OpenMP implementation through `BALLTREE_2BALLS_PRIMARY_FEATURES=1`.

The shared [scalar numerical contract](../../docs/3pcf.rst) defines
observer-centered tangent angles, chord bins, weighted raw distinct triplets,
and the undefined-bearing policy. Use `no-normalize-HistZeta,weights-norm`
for raw comparisons. Repeated neighbors are already excluded natively.

Every rank builds the same fixed tree frontier. Frontier slot `i` belongs to
rank `i % nranks`, OpenMP processes the owned slots, and task-indexed
histograms are reduced to rank 0. Rank 0 publishes tasks in the same order as
the OpenMP method and is the only rank that writes output files.

`options=legacy-one-ball` selects the privately linked distributed FCFC
ball-tree compatibility kernel. This replaces the default-profile need for the
removed standalone one-ball search while retaining its runtime controls.

Build and run with:

```text
make BALLTREE2BALLSMPION=1 cballs
mpiexec -n 4 ./cballs parameters.ini \
    search=balltree-2balls-mpi numberThreads=4
```

`TWOPCFON` and `TPCFON` select the compiled correlation orders.
`only-2pcf` and `only-3pcf` select one order at runtime when both are built.
The remaining options match `balltree-2balls-omp`, including `theta`,
`nsmooth`, `no-two-balls`, `weights-norm`, `no-normalize-HistZeta`,
`compute-HistN`, and `and-CF`.

The default 3PCF forms per-pivot neighbor multipoles and contracts radial-bin
moments; it does not enumerate every distinct triple. No separate 3PCF addon
is required for this production path. `dual-node-direct-triples` selects the
small-catalog validation traversal, whose worst-case work is cubic. MPI
partitions work without changing either estimator's asymptotic cost.

Approximation and smoothing must be qualified for the requested observable
and catalog geometry. Raw-multipole agreement alone does not qualify a
window-corrected result; retain the workload qualification and validity masks.

Run the focused regression with:

```text
make test-balltree-2balls-mpi
```

## Optional hierarchical scalar 3PCF

`options=scalar-pivot-reuse,no-smooth-pivot` enables shared neighbor moments
and completion of resolved radial-bin pairs at parent pivots. It supports
both OpenMP and MPI, preserves strict radial cutoffs, and leaves the independent
2PCF path in place. Runtime controls are `CBALLS_SCALAR_PIVOT_TOL` (default 0.1)
and `CBALLS_SCALAR_BIN_THETA` (default 0). Numerical qualification and speed
depend on the observable, controls, and geometry. See the
[hierarchy contract, MPI rules, tests, and benchmark commands](../../docs/SCALAR_HIERARCHICAL_REUSE.md).
