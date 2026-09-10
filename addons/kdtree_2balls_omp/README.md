# KD-tree two-ball OpenMP addon

`search=kdtree-2balls-omp` builds a median-split KD tree and applies the
dual-node cell-scanning strategy. Its 2PCF is a genuine dual-node traversal:
two nodes are accepted only when their summed radii fit the requested radial
bin tolerance; otherwise the larger or comparable nodes are split with
dual-node's `0.585` rule. The 3PCF uses exact pivots and a dual-node-style
LogMultipole scan of the neighbor KD tree, including second-moment subtraction
for repeated neighbors.

`nsmooth` is the maximum leaf occupancy and `theta` controls radial and angular
bin slop. `options=no-two-balls` is the exact validation path, while
`dual-node-bin-slop` selects dual-node-compatible near-bin acceptance.
`only-2pcf` and `only-3pcf` avoid all work and output for the other order.

`options=legacy-one-ball` is a compatibility mode that dispatches to the
privately linked one-ball KD implementation before the two-ball tree is built. In that
mode, one-ball node acceptance is the default and `no-one-ball` selects exact
body traversal. Smooth pivots, masks, `only-2pcf`, `only-3pcf`, normalization,
and edge correction therefore retain their one-ball compatibility meanings.
Two-ball-only controls (`no-two-balls`, `dual-node-bin-slop`, and
`dual-node-direct-triples`) are rejected. The same compatibility option is
available through `search=kdtree-2balls-mpi`, where the legacy frontier and
histogram reductions use that method's active MPI communicator.

Masks are applied while building both trees. With `SMOOTHPIVOTON=1`, smoothing
is enabled by default and `no-smooth-pivot` disables it. The smoothing prepass
is deterministic. It constructs a tree of active smoothed pivots and a separate
tree of raw neighbors; original body identity is retained to remove self-pairs
exactly. `dual-node-direct-triples` is available with OpenMP or MPI, but only
with `no-smooth-pivot`; MPI distributes its deterministic task frontier and
reduces task-local histograms before root-only normalization and publication.

`edge-corrections,no-normalize-HistZeta` accumulates the mask/window modes
through order `2*mChebyshev` and applies the shared complex edge-correction
solve. The implementation inherits `weights-norm`, `compute-HistN`,
`out-m-HistZeta`, in-memory catalogs, and recoverable error handling from the
dual-node-style estimator contract.

The traversal is adapted from dual-node by Mike Jarvis under its BSD license;
the full notice is in `addons/balltree_2balls_omp/dual-node_LICENSE`.
