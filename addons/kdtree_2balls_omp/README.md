# KD-tree two-ball OpenMP addon

`search=kdtree-2balls-omp` builds a median-split KD tree and applies the
dual-node cell-scanning strategy. Its 2PCF is a genuine dual-node traversal:
two nodes are accepted only when their summed radii fit the requested radial
bin tolerance; otherwise the larger or comparable nodes are split with
dual-node's `0.585` rule. The approximate 3PCF recursively splits aggregate
pivot nodes and carries their unresolved neighbor frontiers into each child;
leaf pivots complete the remaining LogMultipole work. Second-moment subtraction
removes repeated neighbors.

`nsmooth` is the maximum 3PCF leaf occupancy and `theta` controls radial and
angular error acceptance. `options=no-two-balls` is the exact validation path, while
`dual-node-bin-slop` selects dual-node-compatible near-bin acceptance.
`only-2pcf` and `only-3pcf` avoid all work and output for the other order.
For an isolated 2PCF run, the KD implementation selects 4-, 8-, or 16-point
leaves from the full-sky catalog density, relative radial-bin width, and the
requested `theta` error budget. The thresholds were calibrated against exact
pair histograms; dense catalogs and narrow bins select smaller leaves. Add
`dual-node-bucket-leaves` to force `nsmooth` leaves, or
`dual-node-singleton-leaves` for singleton leaves.

The KD specialization uses natural-log scaling for logarithmic radial-bin
indices. This is algebraically equivalent to the base-10 expression and keeps
the same bin assignments. Terminal leaf pairs are gathered into fixed-size
batches so distance filtering and logarithmic bin lookup are no longer mixed
with histogram branches. On macOS, batches of at least eight logarithms use
Accelerate/vForce; other platforms and short tails use the scalar libm path.
The exact batched kernel is covered by an independent-oracle and OpenMP
determinism test.

The LogMultipole path caches each pivot's spherical tangent basis and norm,
reuses them for angular extent tests, and reads exact neighbor bodies from the
contiguous packed-point array. In approximate 3PCF mode, splitting a pivot
carries only its unresolved neighbor-node frontier into each child.
Completed outer radial rings are transported to the child's tangent basis and
reused; the scan no longer restarts every body pivot at the neighbor root.

KD construction remains serial below 262144 selected bodies. At and above that
cutoff, subtrees of at least 65536 bodies are OpenMP tasks whose disjoint point
ranges and preorder node-index ranges are assigned before execution. This
keeps topology and output deterministic while avoiding task overhead on small
catalogs.

Add `dual-node-profile` to emit native phase diagnostics for tree build,
top-level frontier construction, pair/pivot traversal, tangent-basis
transport, scratch clearing, multipole products, and reduction/publication.
Wall phases are elapsed seconds; the three hot-kernel fields are summed
thread-seconds (and rank-reduced by the MPI partner). Profiling is opt-in.

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

The traversal follows the dual-node method by Mike Jarvis under its BSD license;
the full notice is in `addons/balltree_2balls_omp/DUAL_NODE_LICENSE`.
