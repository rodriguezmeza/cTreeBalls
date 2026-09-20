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
the 3PCF. With `BALLS4SCANLEVON=1`, an adaptive work-estimated pivot-cell
frontier replaces fixed pivot blocks. Each task retains exact body pivots,
conservatively prefilters disjoint neighbor roots, and publishes a task-local
histogram in spatial order under a 256 MiB memory cap.

## Opt-In Pivot-Cell Reuse

`options=shear-pivot-reuse,no-smooth-pivot,only-3pcf` enables an experimental
3PCF traversal in this OpenMP addon with `BALLS4SCANLEVON=1`. It is off by
default. The ball-tree OpenMP addon also implements this option; the KD-tree
and MPI addons do not. Calibrate each tree implementation separately.

Each pivot cell resolves neighbors that are safe for every member pivot and
passes only the unresolved nodes to its children. Accepted partial multipoles
stay at their original tangent frames, without copying a full scratch level
to each child. A completed branch combines those rings with direct
ancestor-to-final-frame transport. When all neighbors are resolved, the branch
uses the cell's weighted pivot shear and weight instead of visiting its bodies.
Distinct-neighbor second moments and window multipoles follow the same reuse.
The separate 2PCF path is unchanged, including in a combined run.

Set `CBALLS_SHEAR_PIVOT_TOL` to a finite value in `[0,3]` (default `0.1`). This
is a **phase budget in radians**, not a percentage tolerance on the result.
Zero disables the new traversal. Acceptance requires full radial-bin
containment for the combined pivot/neighbor extent. The phase budget includes
the highest stored ring order, the spherical bearing derivative, a holonomy
allowance for changing the transport path, and the stored cell transport
errors. Caps touching zero or antipodal separation are subdivided. Accepted
rings are transported directly from their original frames, avoiding untracked
errors from repeatedly transporting them along the pivot hierarchy.

For example, add the following to an otherwise unchanged native 3PCF command:

```bash
CBALLS_SHEAR_PIVOT_TOL=0.1 ./cballs your-run.params \
  searchMethod=octree-shear-sphere-2balls-omp \
  options=GGGCorrelation,no-smooth-pivot,only-3pcf,shear-pivot-reuse
```

Calibrate on the same catalog, weights, mask, bins, multipole order, and edge
correction as the intended production run. Compare against
`options=GGGCorrelation,no-smooth-pivot,no-one-ball,only-3pcf`, and inspect both
raw and corrected multipoles, including weak coefficients. A small phase
budget does not guarantee a fixed relative error after cancellations or an
ill-conditioned window correction. Nor does it guarantee a speedup: sparse
catalogs and strict budgets can make inherited-list work more expensive than
the existing body-pivot scan. Retain the default traversal unless calibration
demonstrates the required accuracy and a measured performance benefit.

`no-one-ball`, `no-two-balls`, smoothing, nonpositive `theta`, partial pivot
intervals, inactive unmasked pivots, and pair-only runs disable reuse. As before, `no-two-balls` alone
does not make the existing 3PCF neighbor scan exact: use `no-one-ball`.
`legacy-one-ball` and unsupported shear engines reject the new option.
Thread-local reuse storage is bounded to 64 MiB per worker; allocation failure
aborts the correlation without publishing a partial result.

`CBALLS_SHEAR_PROFILE=1` additionally reports `SHEAR_REUSE` counters for
aggregate pivots, represented bodies, accepted interactions, ancestor-ring
merges, unresolved-list sizes, depth, and scratch bytes. Reuse `walk` timing
includes its clearing and reduction work, so these reported subphases overlap.
Run `make test-shear-sphere-2balls` for the standard and reuse regression suites.

The hot path reuses radial-bin edges, bin-major ring storage, nonnegative
weight multipoles, SIMD-friendly complex accumulation, and normalized body
positions. Native aggregate-cell centers are still normalized when used;
unlike binary-tree centers, they are not assumed to lie on the unit sphere.
Set `CBALLS_SHEAR_PROFILE=1` to report per-thread clearing, tree walk,
transport, radial lookup, ring accumulation, and reduction timings.

The pair traversal reuses normalized centers and radii for acceptance,
splitting, and accumulation. Tangent bases and angular radii are evaluated
only when needed. The pivot scan rejects cells before computing bearings
and caches its angular acceptance constants. Logarithmic bin assignment
retains the original expression and boundary behavior. These are implementation
optimizations; they do not loosen acceptance thresholds or change histogram
normalization.

For exact unsmoothed 3PCF validation use `options=no-smooth-pivot,no-one-ball,only-3pcf`.
The `no-two-balls` switch alone disables pair-node acceptance in the 2PCF
path; it is not a substitute for `no-one-ball` in the 3PCF pivot scan.

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

The independent symmetric dual-tree 2PCF uses chunked frontier tasks. It
performs one transport solve per unordered pair, reuses the conjugate reverse
rotation, and merges tasks deterministically.

The split heuristic is adapted from dual-node, copyright Mike Jarvis, under
dual-node's BSD-style license. See `dual-node_LICENSE`.
