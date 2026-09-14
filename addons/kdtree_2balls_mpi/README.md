# KD-tree two-ball MPI addon

`search=kdtree-2balls-mpi` is the MPI+OpenMP counterpart of
`kdtree-2balls-omp`. Every rank builds the same read-only median KD tree and
owns a deterministic cyclic subset of the pair or pivot frontier. Histograms,
normalization planes, edge-window modes, and counters are reduced to rank zero;
only rank zero publishes files or Cython-visible results.

The MPI engine supports the same mask, smooth-pivot, two-ball, exact,
`only-2pcf`, `only-3pcf`, weighting, and edge-correction options as the OpenMP
engine.
`dual-node-direct-triples` distributes the direct validation frontier and must
be combined with `no-smooth-pivot`, as in the OpenMP engine. It requires
`MPI_THREAD_FUNNELED`; MPI calls remain outside OpenMP worker regions.

The MPI partner compiles the same adaptive 2PCF leaves, batched vector-log
kernel, persistent unresolved-neighbor 3PCF frontier, and deterministic
large-subtree OpenMP builder as the OMP method. `dual-node-profile` reports the
root rank's build/frontier/traversal/reduction wall times and rank-reduced
transport, scratch-clear, and multipole-product thread-seconds.
