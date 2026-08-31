# Octree BALLS4 MPI

Build with `OCTREEBALLS4MPION=1`; select `search=octree-balls4-mpi`
(method ID 184). The OpenMP addon need not also be enabled.
Requires an MPI C compiler, OpenMP, 3D, nonperiodic geometry, and log bins.
Use `TWOPCFON=1` for 2PCF and `TPCFON=1` for 3PCF.

Each rank loads the same catalog/tree. With the all-engines Python driver,
rank 0 reads the catalog once and broadcasts it using mpi4py before cyballs
copies it into C memory. MPI must provide at least MPI_THREAD_FUNNELED.
MPI initialized by Python is not finalized by this addon.

Normal runs assign disjoint B4 frontier nodes to ranks. Raw histograms are
reduced before normalization. Raw and edge runs use fixed blocks of body pivots from
that partition; block reduction order is independent of thread/rank count.
Only rank 0 owns published results and writes output.

`read-mask` and `edge-corrections,no-normalize-HistZeta` are supported.
`weights-norm` selects weighted complex signal/window moments.
Window modes extend through twice `mChebyshev`; empty/singular bins are zero.
Edge runs and uncorrected `no-normalize-HistZeta` runs exclude repeated
neighbors and use body pivots. Only legacy normalized non-edge runs retain
the original cell-pivot estimator. Either `no-one-ball` or `no-two-balls`
makes the raw/edge path exact. `smooth-pivot` is rejected independently of
`SMOOTHPIVOTON`. Observer-frame, chord-bin, and complex-mode conventions
are described in [the scalar contract](../../docs/3pcf.rst).

```sh
mpiexec -n 2 ./cballs parameters.ini search=octree-balls4-mpi numberThreads=4
make test-octree-balls4-mpi
```

See `docs/addons.rst` and `python/README_kappa_corr_all_engines.md`.
