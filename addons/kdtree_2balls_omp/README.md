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

## Completed-pivot progress

The production 3PCF reports completed active pivots in exact and approximate
mode. Add `stepState=10000 verbose=1 verbose_log=1` to the cballs command or
parameter file. For example:

```text
kdtree-2balls-omp: 3PCF progress: completed pivots 10000 / 100000 (10.0%); elapsed 2.35 s
```

Counts start at zero and end at the pivot-tree population, after masking and
any smoothing. A smoothed representative counts as one active pivot, not the
number of original bodies it represents. Counts describe completed work, not
input catalog indices or partial-ring operations.

The approximate KD traversal may finish all remaining radial bins for an
entire pivot node at once. Its active pivots are then counted together; a node
with unresolved bins is not counted until its descendants finish. Consequently
progress can jump by a whole group. Ordinary body completions are published in
local batches of at most 64 (smaller for smaller `stepState`), with any remainder
published when a task finishes. Output is serialized and flushed so counts
remain monotonic under dynamic OpenMP scheduling.

`verbose=0 verbose_log=1` sends progress only to `cballs.log`;
`verbose=1 verbose_log=0` sends it only to the terminal; setting both to zero
disables counting/publication overhead. `stepState=1` reports each completed
body or group and can generate substantial output. The final 100% line means
pivot traversal is complete; histogram reduction, normalization, edge correction
and output may still take time.

This counter is for the OpenMP production 3PCF only. A combined 2PCF/3PCF run
starts it after the independent pair traversal. `only-2pcf`, the
`dual-node-direct-triples` validation path, and MPI are unchanged. Compatibility
mode continues to use the legacy kernel's own reporting.

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
the full notice is in `addons/balltree_2balls_omp/DUAL_NODE_LICENSE`.
