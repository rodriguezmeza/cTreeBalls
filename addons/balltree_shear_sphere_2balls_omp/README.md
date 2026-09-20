# balltree-shear-sphere-2balls-omp

This addon computes full-sky weak-lensing shear `xi+`, `xi-`, and all four
natural Gamma-x 3PCF multipole components over an FCFC-style PCA ball tree.
Input positions are observer-centered 3D vectors and shear components are
defined in each sample's local east/north frame.

Nodes are split at the median of their dominant principal axis and carry a
conservative spherical chord-radius bound. Their weighted first and second
spin-2 moments are parallel transported into the tangent frame at the node
center. The estimator shares the symmetric dual-node 2PCF and accepted
neighbor-node LogMultipole 3PCF traversal with the KD-tree variant, including
the combined-radius opening criterion and 0.585 split heuristic.

`only-2pcf` and `only-3pcf` skip unused work. Either `no-two-balls` or
`no-one-ball` disables cell aggregation in this binary-tree implementation,
for both 2PCF and 3PCF. `dual-node-bin-slop` enables the looser dual-node
radial criterion, and `nsmooth` sets leaf capacity. Masks, shear
mode-coupling edge correction, deterministic `BALLS4SCANLEVON`, and the
`SMOOTHPIVOTON`/`no-smooth-pivot` contract match the spherical octree and KD
engines. Explicit smoothing radii must satisfy `2*rsmooth <= rminHist`.

## Construction and Pair Projection

The builder skips PCA on terminal leaves, where no split direction is needed.
It allocates exactly the required node count and uses deterministic preorder
node ranges. For at least 16,384 selected points, large sibling subtrees can
build concurrently; tasks stop below an 8,192-point subtree. Point order within
each node, moment summation order, and the resulting tree are independent of
the worker schedule. Upper-node center and PCA statistics share one member
pass. At 16,384 members and above, statistics and direct spin-2/radius passes
use up to 32 fixed chunks (target 4,096 members), reduced in index order.
The same chunks are used with one thread, so changing the worker count does
not change the tree or its moments. Smaller nodes create no range tasks.
`no-balltree-parallel-build` selects the serial construction diagnostic while
retaining the leaf-PCA and allocation improvements.

Normalized source directions, tangent bases, weights, and weighted shears are
prepared once per build and reused at all levels. This temporary cache is
freed before searching and capped at 256 MiB (80 bytes per input body in the
double-precision build). If the cap or an allocation fails, construction uses
the identical direct source-frame calculation. `no-balltree-shear-member-cache`
forces that diagnostic path. It changes neither moment definitions nor
transport accuracy. Bounds are rounded outward and include displacement
between normalized and stored centers.

Unsmoothed repeated catalog roles share a tree within one correlation call.
There is no cross-call shear-tree cache. Smoothed pivot and neighbor roles
retain separate trees. First and second moments are still summed directly from
node members, so parallel construction introduces no new transport approximation.
The shared builder also serves the MPI partner; this addon's new pair projection
is enabled only in the OpenMP wrapper.

The symmetric 2PCF projects both shears onto their connecting great circle.
Squared bearings are normalized directly to avoid cancellation at small
separations. One update accounts for both orientations of an auto pair;
cross-catalog pairs retain complex xi-plus. This removes explicit pair transport
and a duplicate reverse update without changing opening criteria, radial bins,
or estimator definitions. Floating-point rounding can differ from the previous
transport-then-project implementation.

## Optional Pivot-Cell Reuse

The default 3PCF still uses individual body pivots. With `BALLS4SCANLEVON=1`,
`options=shear-pivot-reuse,no-smooth-pivot` enables a separate, calibrated
approximation: pivot cells inherit sparse unresolved neighbor lists and
accepted multipole rings. Rings remain in their original tangent frame and
are rotated directly into the final pivot frame once; no full scratch array
is copied to each child. Fully resolved cells use their summed pivot shear
and weight. Repeated-neighbor second-moment subtraction is retained.

Every accepted interaction must fit entirely inside one radial bin. A
conservative angular-cap test budgets bearing variation, spherical transport
holonomy, and stored moment error, rejecting coincident and antipodal caps.
`CBALLS_SHEAR_PIVOT_TOL` sets the phase budget in radians (default `0.1`, range
`[0,3]`). While reuse is enabled, this budget replaces `theta` as the 3PCF
angular-acceptance control; positive `theta` is required to enable it. The
ordinary body-pivot and independent 2PCF paths retain their usual `theta` tests.
This is **not a percentage-error bound**: cancellation and
mode-coupling inversion can amplify small phase errors. Calibrate each catalog,
mask, separation range, binning, and multipole order against its own exact run.
Check every raw, window, and corrected complex coefficient, not just RMS error.

Zero tolerance, `theta=0`, `no-one-ball`, `no-two-balls`, smoothing, partial
pivot ranges, or inactive valid pivots fall back to the body-pivot path. Zero
tolerance alone does not disable ordinary neighbor-cell approximation.
The option has no effect on the independent 2PCF. Unsupported engines,
including the MPI partner, and legacy mode reject it. Per-worker reuse scratch
is bounded by 64 MiB; allocation or transport failure aborts before publication.
No run-to-run body/Python pointers are cached.

For example, with a pre-existing parameter file (calibrate the budget first):

```bash
CBALLS_SHEAR_PIVOT_TOL=0.03 ./cballs your-run.params \
  searchMethod=balltree-shear-sphere-2balls-omp numberThreads=8 \
  options=GGGCorrelation,no-smooth-pivot,only-3pcf,shear-pivot-reuse
```

Profile with `CBALLS_SHEAR_PROFILE=1`: `binary_tree_build`
reports build wall time and `unique_trees`; per-worker counters report radial,
transport, ring, clearing, traversal, and product work. Sampled kernel estimates
overlap with traversal time and must not be summed as independent CPU phases.
`SHEAR_REUSE` additionally reports aggregate pivots, represented bodies,
ancestor-ring merges, unresolved-list peaks, and allocated scratch bytes.

Construction and pair optimizations require no new option. For an exact,
unsmoothed 3PCF reference with an existing parameter file:

```bash
./cballs your-run.params searchMethod=balltree-shear-sphere-2balls-omp \
  numberThreads=8 \
  options=GGGCorrelation,no-smooth-pivot,no-one-ball,only-3pcf
```

Validate approximate settings against this reference on the relevant catalog,
binning, and multipole order. Neither `theta` nor the unchanged defaults imply
a universal percentage-error guarantee.

The PCA split and initial enclosing-sphere strategy are adapted from FCFC by
Cheng Zhao under the MIT license. The two-node traversal strategy is adapted
from dual-node by Mike Jarvis under its BSD-style license; see
`dual-node_LICENSE`.
