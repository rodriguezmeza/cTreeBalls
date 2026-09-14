# octree-shear-sphere-2balls-omp

This independent addon computes full-sky spin-2 `xi+`, `xi-`, and Gamma-x
multipole estimators. Neighbor shears are
parallel transported along great circles into the pivot tangent frame.

Its 2PCF uses a dual-node traversal over the native
cTreeBalls octree. For nodes of radii `s1` and `s2` at chord distance `d`, the
production path accepts their aggregate shear when their combined angular
extent is within the radial-bin and spin-phase tolerance. Otherwise it splits
the larger node; comparable nodes may both split using dual-node's empirical
`0.585` rule. `options=no-two-balls` or `options=no-one-ball` selects the exact
body-pair fallback.

The 3PCF retains the spherical LogMultipole tree scan. This preserves the
existing cTreeBalls radial-multipole result contract rather than changing to
the paper's `(r,u,v)` triangle histogram. In a combined 2PCF+3PCF run, the
dual-node pair pass is performed separately and the pivot scan computes only
the 3PCF. With `BALLS4SCANLEVON=1`, an adaptive work-estimated pivot-cell
frontier replaces fixed pivot blocks. Each task retains exact body pivots,
conservatively prefilters disjoint neighbor roots, and publishes a task-local
histogram in spatial order under a 256 MiB memory cap.

The hot path reuses radial-bin edges, bin-major ring storage, nonnegative
weight multipoles, SIMD-friendly complex accumulation, and normalized body
positions. Native aggregate-cell centers are still normalized when used;
unlike binary-tree centers, they are not assumed to lie on the unit sphere.
Set `CBALLS_SHEAR_PROFILE=1` to report per-thread clearing, tree walk,
transport, radial lookup, ring accumulation, and reduction timings.

Build with `OCTREESHEARSPHERE2BALLSOMPON=1`. The addon supports `only-2pcf`,
`only-3pcf`, masks, shear edge correction, `SMOOTHPIVOTON`, and
`BALLS4SCANLEVON`. Smoothing is enabled by default when compiled; use
`options=no-smooth-pivot` to enable the dual-node 2PCF path. A smoothed run
uses representative body pivots because dual-node traversal has no unique
pivot ownership. Its literal radius must satisfy `2*rsmooth <= rminHist`; an
automatically selected radius is capped at this bound so smooth groups cannot
contain measured pairs.

The independent symmetric dual-tree 2PCF uses chunked frontier tasks. It
performs one transport solve per unordered pair, reuses the conjugate reverse
rotation, and merges tasks deterministically.

The split heuristic is adapted from dual-node, copyright Mike Jarvis, under
dual-node's BSD-style license. See `dual-node_LICENSE`.
