# FCFC-style octree-GGG MPI addon

This addon keeps the production `octree-ggg-omp` estimator and adds hybrid
MPI/OpenMP execution as `search=octree-ggg-mpi`. MPI ranks claim pivot batches
dynamically from a rank-0 RMA counter, following FCFC's scheduler. Each rank
uses its local octree and OpenMP workers; raw histograms and counters are
reduced to rank 0 before normalization and output.

The OMP/MPI siblings share the repaired observer-tangent angular contract,
chord bins, runtime-only smoothing, and computed-histogram publication.
Use `no-normalize-HistZeta,weights-norm` for raw distinct-triplet
comparisons. `read-mask` is supported; complex edge correction requires
`edge-corrections,no-normalize-HistZeta`, `NMultipolesON=1`, and
`NONORMHISTON=1`. Only rank 0 owns published output.
See [the scalar guide](../../docs/3pcf.rst) and
[the OpenMP README](../octree_ggg_omp/Readme.txt).

Build and run with:

```bash
make OCTREEGGGMPION=1 cballs
mpiexec -n 4 ./cballs parameters.ini \
  search=octree-ggg-mpi numberThreads=2
```

The selected MPI implementation must provide `MPI_THREAD_FUNNELED`. Only rank
0 creates output files or runs preprocessing and post-processing commands.

Run the OMP/MPI regression with:

```bash
make OCTREEGGGMPION=1 test-octree-ggg-mpi
```

Rebuild Cython with the same profile before Python tests. Import mpi4py.MPI
first; every rank must enter the same run and cleanup sequence.
