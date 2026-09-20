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
is owned by rank `i % nranks`; OpenMP distributes the rank-local slots. The
3PCF keeps independent task histograms and the same adaptive inherited
neighbor frontier as the OpenMP implementation. The independent 2PCF instead
uses one reusable histogram per OpenMP thread, combines them with a binary
tree reduction, and performs a single aggregate MPI reduction before rank 0
publishes the result. This removes dense per-task pair storage and serial
publication while preserving the numerical contract. Each catalog also uses
the OpenMP implementation's density/bin-width/`theta` leaf-capacity policy;
exact mode keeps the configured `nsmooth` value.

Each rank also uses the deterministic compact-view topology and parallel
post-order geometry build described by `octree-2balls-omp`. Linux vector-log
selection is shared through `NATIVE_PAIR_VECTOR_LOG=auto|sleef|compiler`.

Build and run:

```sh
make OCTREE2BALLSMPION=1 cballs
mpiexec -n 4 ./cballs parameters.ini \
    search=octree-2balls-mpi numberThreads=4
```

`only-2pcf`, `only-3pcf`, `no-two-balls`, `weights-norm`, and
`no-normalize-HistZeta` have the same meaning as in `octree-2balls-omp`.
The cubic `dual-node-direct-triples` validation path is intentionally not
distributed; run it with `octree-2balls-omp` on a small catalog.

## One-ball compatibility

`options=legacy-one-ball` switches this search name to the privately linked
distributed one-ball implementation. It builds the full threaded native octree and
adaptive frontier and uses the dynamic MPI scheduler, histogram
reductions, root-only output, and finalizer. The native two-ball MPI runtime is
not initialized in this mode.

`behavior-ball`/`no-one-ball`, `compute-HistN`, masks, normalization, complex
edge correction, `ggg-full-window`, and `ggg-profile` retain their GGG
meanings. With `SMOOTHPIVOTON=1`, smoothing is enabled by default and
`no-smooth-pivot` disables it. `only-2pcf` is supported; `only-3pcf` is
rejected because GGG does not provide a true skip-2PCF execution path. Do not
combine compatibility mode with `no-two-balls`, `dual-node-bin-slop`, or
`dual-node-direct-triples`.

The compatibility objects are private dependencies of this addon, so a narrow
build with `OCTREEGGGOMPON=0 OCTREEGGGMPION=0` still supports the alias without
registering either public GGG search name.

Run the regression test with:

```sh
make test-octree-2balls-mpi
```
