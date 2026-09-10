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
standalone search names are disabled in the default build profile.

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
`addons/balltree_omp/fcfc_balltree.c`. dual-node's redistribution terms are in
`dual-node_LICENSE`.
