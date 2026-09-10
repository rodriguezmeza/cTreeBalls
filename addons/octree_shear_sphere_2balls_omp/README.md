# octree-shear-sphere-2balls-omp

This independent addon computes full-sky spin-2 `xi+`, `xi-`, and Gamma-x
multipole estimators. Neighbor shears are
parallel transported along great circles into the pivot tangent frame.

Its 2PCF uses a dual-node/Jarvis-style dual-node traversal over the native
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
the 3PCF.

`options=legacy-one-ball` dispatches this search name to the privately linked
one-ball compatibility kernel. It uses a combined body-pivot 2PCF/3PCF traversal
and preserves its order
selectors, accepted-cell policy, spherical transport, smoothing, masks,
mode-coupling correction, normalization, output, and deterministic reduction.
The native two-ball traversal remains the default when the option is absent.
The addon privately links the ordinary spherical kernel when
`OCTREESHEARSPHEREOMPON=0`, so this compatibility mode also works in a narrow
two-ball-only build.

Build with `OCTREESHEARSPHERE2BALLSOMPON=1`. The addon supports `only-2pcf`,
`only-3pcf`, masks, shear edge correction, `SMOOTHPIVOTON`, and
`BALLS4SCANLEVON`. Smoothing is enabled by default when compiled; use
`options=no-smooth-pivot` to enable the dual-node 2PCF path. A smoothed run
uses representative body pivots because dual-node traversal has no unique
pivot ownership. Its literal radius must satisfy `2*rsmooth <= rminHist`; an
automatically selected radius is capped at this bound so smooth groups cannot
contain measured pairs.

The split heuristic is adapted from dual-node, copyright Mike Jarvis, under
dual-node's BSD-style license. See `dual-node_LICENSE`.
