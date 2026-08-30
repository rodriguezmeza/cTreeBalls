# balltree-2balls-mpi

This addon is the deterministic MPI+OpenMP counterpart of
`balltree-2balls-omp`. It uses the same FCFC PCA ball tree, TreeCorr-style
dual-node 2PCF traversal, and genuine `process3`/`process21`/`process111`
triple-node 3PCF traversal.

Every rank builds the same fixed tree frontier. Frontier slot `i` belongs to
rank `i % nranks`, OpenMP processes the owned slots, and task-indexed
histograms are reduced to rank 0. Rank 0 publishes tasks in the same order as
the OpenMP method and is the only rank that writes output files.

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

The direct triple-node estimator has cubic worst-case work. MPI distributes
that work but does not change its asymptotic cost; use the LogMultipole addons
for production 3PCF calculations on large catalogs.

Run the focused regression with:

```text
make test-balltree-2balls-mpi
```
