# FCFC-style Ball Tree OpenMP/MPI Addon

This addon adds `search=balltree-omp` for cTreeBalls 2PCF and 3PCF runs.
It uses the balanced PCA split and enclosing-sphere strategy described by FCFC,
while retaining cTreeBalls' pivot-local histogram and angular-multipole kernels.

The leaf capacity is controlled by `nsmooth` (FCFC recommends `8`). Exact
point traversal is the default. `behavior-ball` enables cTreeBalls' approximate
cell aggregation; `no-one-ball` explicitly restores exact traversal.

The adapted FCFC construction code is distributed under the MIT license. See
the notice in `fcfc_balltree.c`.

## Hybrid MPI and OpenMP

Build the optional MPI profile with an MPI C compiler available as `mpicc`:

```sh
make BALLTREEMPION=1 cballs
mpiexec -n 4 ./cballs parameters.ini search=balltree-mpi numberThreads=2
```

`balltree-mpi` builds each catalog tree on rank 0 and broadcasts it. Rank 0
also creates a same-level ball-tree frontier. MPI ranks claim frontier blocks
dynamically through an MPI-3 RMA counter, while OpenMP threads process the
pivots in each block. Histograms and counters are reduced to rank 0, which is
the only rank that writes output files or runs pre/post-processing scripts.

For 3PCF calculations, all neighbors of a pivot remain on the same rank. This
is required because the angular multipole product is pivot-local; MPI
parallelism is therefore over pivot blocks rather than over partial neighbor
lists. Runs assume a homogeneous MPI job because the internal tree-node
representation is broadcast between ranks.

Run the regression with:

```sh
make test-balltree-mpi
```
