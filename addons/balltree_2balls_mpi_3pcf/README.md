# balltree-2balls-mpi_3pcf

This addon is the deterministic MPI+OpenMP counterpart of
`balltree-2balls-omp_3pcf`. It uses the same FCFC PCA ball tree and TreeCorr
LogMultipole estimator. A fixed pivot-node frontier is built independently on
every rank, frontier slot `i` is owned by rank `i % nranks`, and OpenMP scans
the owned slots. Task-indexed histograms are reduced to rank 0 and published
in the same fixed order as the OpenMP implementation.

The shared kernel transports inherited moments between pivot tangent bases
and bounds projected angular error. Preserve observer-centered positions
and chord bins in 3D. Raw `no-normalize-HistZeta,weights-norm` products
are weighted distinct-triplet sums, not normalized correlations.
See [the scalar contract](../../docs/3pcf.rst).

Build and run with:

```text
make BALLTREE2BALLSMPI3PCFON=1 TPCFON=1 cballs
mpiexec -n 4 ./cballs parameters.ini \
    search=balltree-2balls-mpi_3pcf numberThreads=4
```

`numberThreads` is the OpenMP thread count per MPI rank. The production
options and normalization policy are identical to
`balltree-2balls-omp_3pcf`: `theta`, `nsmooth`,
`treecorr-singleton-leaves`, `no-two-balls`, `weights-norm`, and
`no-normalize-HistZeta` are supported. The cubic
`treecorr-direct-triples` validation traversal is intentionally not
distributed; run that oracle with `balltree-2balls-omp_3pcf`.

Only rank 0 writes histogram files. Run the focused regression with:

```text
make test-balltree-2balls-mpi-3pcf
```
