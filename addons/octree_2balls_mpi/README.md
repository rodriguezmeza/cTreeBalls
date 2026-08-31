# octree-2balls-mpi

`octree-2balls-mpi` is the MPI+OpenMP counterpart of
`octree-2balls-omp`. It uses the same native-octree binary view, dual-node
2PCF traversal, and LogMultipole 3PCF estimator.

It shares the repaired observer-frame, tangent-angle, chord-bin, and distinct
raw-moment contract described in [the scalar guide](../../docs/3pcf.rst).
`read-mask` and complex `edge-corrections,no-normalize-HistZeta` are
supported; `weights-norm` applies catalog weights. Raw and corrected
products are not interchangeable.

The deterministic frontier is reproduced on every rank. Frontier slot `i`
is owned by rank `i % nranks`; OpenMP distributes the rank-local slots. Each
slot has an independent histogram, and MPI reduces those slot arrays before
rank 0 publishes them in the same fixed order as the OpenMP implementation.
This makes rank-count changes reproducible without serializing the search.

Build and run:

```sh
make OCTREE2BALLSMPION=1 cballs
mpiexec -n 4 ./cballs parameters.ini \
    search=octree-2balls-mpi numberThreads=4
```

`only-2pcf`, `only-3pcf`, `no-two-balls`, `weights-norm`, and
`no-normalize-HistZeta` have the same meaning as in `octree-2balls-omp`.
The cubic `treecorr-direct-triples` validation path is intentionally not
distributed; run it with `octree-2balls-omp` on a small catalog.

Run the regression test with:

```sh
make test-octree-2balls-mpi
```
