# KD-tree two-ball MPI addon

`search=kdtree-2balls-mpi` is the MPI+OpenMP counterpart of
`kdtree-2balls-omp`. Every rank builds the same read-only median KD tree and
owns a deterministic cyclic subset of the pair or pivot frontier. Histograms,
normalization planes, edge-window modes, and counters are reduced to rank zero;
only rank zero publishes files or Cython-visible results.

The MPI engine supports the same mask, smooth-pivot, two-ball, exact,
`only-2pcf`, `only-3pcf`, weighting, and edge-correction options as the OpenMP
engine. `legacy-one-ball` dispatches to the actual legacy KD search while
retaining this addon's communicator and deterministic reductions.
`dual-node-direct-triples` distributes the direct validation frontier and must
be combined with `no-smooth-pivot`, as in the OpenMP engine. It requires
`MPI_THREAD_FUNNELED`; MPI calls remain outside OpenMP worker regions.

The MPI partner compiles the same adaptive 2PCF leaves, batched vector-log
kernel, persistent unresolved-neighbor 3PCF frontier, and deterministic
large-subtree OpenMP builder as the OMP method. `dual-node-profile` reports the
root rank's build/frontier/traversal/reduction wall times and rank-reduced
transport, scratch-clear, and multipole-product thread-seconds.

## Optional hierarchical scalar 3PCF

`options=scalar-pivot-reuse,no-smooth-pivot` enables shared neighbor moments
and completion of resolved radial-bin pairs at parent pivots. It supports
both OpenMP and MPI, preserves strict radial cutoffs, and leaves the independent
2PCF path in place. Runtime controls are `CBALLS_SCALAR_PIVOT_TOL` (default 0.1)
and `CBALLS_SCALAR_BIN_THETA` (default 0). Numerical qualification and speed
depend on the observable, controls, and geometry. See the
[hierarchy contract, MPI rules, tests, and benchmark commands](../../docs/SCALAR_HIERARCHICAL_REUSE.md).
