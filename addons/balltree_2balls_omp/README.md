# TreeCorr-style dual/triple-node OpenMP addon

`search=balltree-2balls-omp` uses the shared FCFC PCA ball tree for two
independent correlation engines:

* `TWOPCFON = 1` enables the scalar 2PCF dual-node traversal.
* `TPCFON = 1` enables the scalar angular-multipole 3PCF triple-node traversal.

Either feature can be enabled alone, or both can be enabled in the same build.
When both are built, `options=only-2pcf` or `options=only-3pcf` runs just one
engine. Supplying both selectors is an error, as is selecting an engine that
was compiled out.
The 3PCF recursion follows TreeCorr's `process3`, `process21`, and `process111`
decomposition: triples contained in one node, split across two nodes, and split
across three nodes are visited separately. Auto-correlations visit every
unordered physical triple once and accumulate all six pivot/leg orientations.
For a two-catalog run, catalog 1 supplies the pivot and catalog 2 supplies the
two legs.

Node pairs retain exact radial-bin containment. Oriented node triples use
TreeCorr-style `bin_slop` and the LogMultipole default angular tolerance,
scaled by `theta`; this is the approximation that permits useful cell
aggregation. The exact reference remains available with `no-two-balls`.
`options=no-two-balls` disables both kinds of aggregation and forces exact body
pairs and triplets, providing the reference path used by the regression tests.
`nsmooth` sets leaf capacity, `theta` scales the split tolerances, and
`weights-norm` selects weighted scalar normalization.

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

The triple-node engine explicitly visits distinct triples. Its worst-case work
is cubic, so it is intended for validation and medium-sized catalogs, not a
full-sky `NSIDE=1024` map. TreeCorr's production LogMultipole implementation
avoids this cost with pair-multipole accumulation, which is a different
algorithm. For a large catalog use `options=only-2pcf`, downsample before a
direct 3PCF run, or use one of cTreeBalls' fast pivot-neighbor multipole
searches such as `octree-ggg-omp`.

The node recursion follows TreeCorr by Mike Jarvis, distributed under its
BSD-style license. The PCA ball-tree construction is adapted from FCFC by
Cheng Zhao under the MIT license; see the notices in
`addons/balltree_omp/fcfc_balltree.c`. TreeCorr's redistribution terms are in
`TreeCorr_LICENSE`.
