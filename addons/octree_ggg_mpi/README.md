# FCFC-style octree-GGG MPI addon

This addon keeps the production `octree-ggg-omp` estimator and adds hybrid
MPI/OpenMP execution as `search=octree-ggg-mpi`. MPI ranks claim pivot batches
dynamically from a rank-0 RMA counter, following FCFC's scheduler. Each rank
uses its local octree and OpenMP workers; raw histograms and counters are
reduced to rank 0 before normalization and output.

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
