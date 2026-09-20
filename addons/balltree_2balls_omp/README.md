# FCFC PCA-tree dual-node and LogMultipole OpenMP addon

`search=balltree-2balls-omp` uses the shared FCFC PCA ball tree for two
independent correlation engines:

* `TWOPCFON = 1` enables the scalar 2PCF dual-node traversal.
* `TPCFON = 1` enables the scalar angular-multipole 3PCF LogMultipole traversal.

Either feature can be enabled alone, or both can be enabled in the same build.
When both are built, `options=only-2pcf` or `options=only-3pcf` runs just one
engine. Supplying both selectors is an error, as is selecting an engine that
was compiled out.
The production 3PCF follows the same body-pivot LogMultipole contract as
`kdtree-2balls-omp`: it accumulates radial field moments over accepted neighbor
nodes, subtracts repeated-neighbor second moments, and forms every radial-bin
pair without enumerating physical triples. `options=dual-node-direct-triples`
retains the older `process3`/`process21`/`process111` traversal as a small-catalog
validation oracle.

Node pairs retain exact radial-bin containment. Oriented node triples use
dual-node-style `bin_slop` and the LogMultipole default angular tolerance,
scaled by `theta`; this is the approximation that permits useful cell
aggregation. The exact reference remains available with `no-two-balls`.
`options=no-two-balls` disables aggregation and forces exact body neighbors,
providing the reference path used by the regression tests.
`nsmooth` sets leaf capacity, `theta` scales the split tolerances, and
`weights-norm` selects weighted scalar normalization.

For an ordinary auto-correlation 2PCF, the production engine chooses a
deterministic leaf capacity from catalog size, radial-bin width, and `theta`.
Catalogs below 262144 active points keep at least eight bodies per leaf;
denser catalogs may use four-body leaves when their smaller cells are likely
to repay the extra construction work.  Terminal body pairs are evaluated in
batches of 256.  Logarithmic bins use the platform vector-log backend on
Apple builds. Linux builds use explicit SLEEF AVX2, SSE2, or AArch64 Advanced
SIMD entry points when `pkg-config` can resolve the library and fall back to
scalar `log` otherwise. Set `SLEEFON=1` to require SLEEF or `SLEEFON=0` to
disable it; `make test-sleef-vector-log` checks both SIMD and scalar-tail
results against libm and rejects a build that selected only one-wide lanes.

The OpenMP auto-2PCF and LogMultipole auto-3PCF paths keep a two-entry,
process-local cache of compact
trees.  Cache keys include catalog contents and all tree-shaping field,
weight, mask, smoothing, and leaf settings, so modifying a catalog causes a
rebuild. Cached trees own packed positions, scalar values, weights, and moments;
they retain no body or Python-model pointers after installation. This is useful
for repeated 2PCFs, repeated 3PCFs, or thread sweeps over the same catalog.
Cross-correlations continue to build independent trees. Add
`no-balltree-tree-cache` when measuring cold construction or when the
application prefers immediate release. With `dual-node-profile`, the search
log reports `compact-tree cache = hit` or `miss`.

Large ball trees use deterministic preassigned node ranges and OpenMP subtree
tasks. Range statistics combine covariance, aggregate moments, and field
moments in deterministic chunks; the enclosing-sphere and aggregate-radius
pass is also fused and parallel above a large cutoff. Sufficiently large
disjoint children build concurrently without allocator contention. Add
`no-balltree-parallel-build` for a serial-construction diagnostic run.

The production 3PCF passes each pivot child a bounded sparse list of unresolved
neighbor nodes. It prunes geometrically irrelevant entries and refines useful
large entries without copying any multipole scratch level between pivot nodes.
Leaf pivots accumulate into their own cleared scratch arrays, preserving the
original arithmetic order. Add `no-balltree-persistent-frontier` to compare
against a neighbor-root restart for every body pivot.

## Completed-pivot progress

The production 3PCF reports completed active pivots in exact and approximate
mode, with either the persistent neighbor frontier or body-root restart.
Add `stepState=10000 verbose=1 verbose_log=1` to the cballs command or parameter
file. Example output:

```text
balltree-2balls-omp: 3PCF progress: completed pivots 10000 / 100000 (10.0%); elapsed 2.35 s
```

The counter starts at zero and ends at the pivot-tree population, excluding
masked bodies and pivots absorbed by smoothing. Each active smoothed
representative counts once. Counts measure finished pivots, not input indices
or node visits. Workers publish local batches of at most 64, smaller when
`stepState` is smaller, and flush any remainder at task completion. A displayed
count can therefore pass the requested interval. Updates and flushed output
are serialized so dynamic OpenMP scheduling cannot reorder progress lines.

Use `verbose=0 verbose_log=1` for log-only reporting and
`verbose=1 verbose_log=0` for terminal-only reporting. Setting both to zero
disables counting/publication overhead. `stepState=1` prints each pivot and can
be expensive. The 100% line precedes histogram reduction, normalization, edge
correction and output; it does not mean the entire application has exited.

Only the OpenMP production 3PCF has this new counter. Combined 2PCF/3PCF runs
start it after the independent pair traversal. `only-2pcf`,
`dual-node-direct-triples`, MPI, and legacy-kernel reporting are unchanged.

`read-mask` removes masked bodies before either role-specific tree is built.
With `SMOOTHPIVOTON=1`, deterministic smooth pivots are enabled by default and
the pivot tree stores the same grouped field and normalization sums as the KD
two-ball engine; `no-smooth-pivot` restores ordinary body pivots. Complex
`edge-corrections,no-normalize-HistZeta` uses window modes through twice the
requested signal order.

`options=legacy-one-ball` dispatches to the actual legacy ball-tree
implementation. Its controls remain unchanged: `behavior-ball` enables
one-ball node aggregation, while `no-one-ball` forces exact traversal. Use it
with `search=balltree-2balls-omp` for OpenMP or
`search=balltree-2balls-mpi` for MPI. Do not combine compatibility mode with
`no-two-balls`, `dual-node-bin-slop`, or `dual-node-direct-triples`. The old
standalone search names are disabled in the default build profile. The legacy
ball-tree kernel has no angular-window solver: combining `legacy-one-ball`
with `edge-corrections` is rejected. Remove `legacy-one-ball` to solve the window.

Three-dimensional angular phases and acceptance use projected tangent
bearings in the original observer frame, with chord-distance bins.
See [the scalar contract](../../docs/3pcf.rst) for the positive complex-mode
convention and undefined-bearing policy. Raw `no-normalize-HistZeta` runs
already exclude repeated neighbors; do not remove them again in Python.

The 3PCF multipoles are normalized at runtime by default. Add
`no-normalize-HistZeta` to retain raw distinct-triplet sums. With
`weights-norm`, the denominator is the distinct-triplet weight sum; otherwise
it is the distinct-triplet count. This runtime choice applies even when the
build has `NONORMHISTON=1`.

With `BALLS4SCANLEVON=1`, pair and direct-triple scans use a balanced spatial
frontier with at least 64 tasks and grow it with the available OpenMP/MPI
worker pool. Task-local histograms are still reduced in frontier order, so the
optimization preserves deterministic results. The production LogMultipole
variant already uses a stronger work-estimated frontier of up to 256 tasks;
the flag enables the same documented frontier contract without replacing that
scheduler.

The direct validation engine explicitly visits distinct triples and has cubic
worst-case work. Do not use `dual-node-direct-triples` for a full-sky catalog;
the default LogMultipole path is the production algorithm.

The node recursion follows dual-node by Mike Jarvis, distributed under its
BSD-style license. The PCA ball-tree construction is adapted from FCFC by
Cheng Zhao under the MIT license; see the notices in
`addons/balltree_shared/fcfc_balltree.c`. dual-node's redistribution terms are in
`dual-node_LICENSE`.
