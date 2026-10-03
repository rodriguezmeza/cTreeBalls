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
for both 2PCF and 3PCF. `dual-node-bin-theta` enables the looser dual-node
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
from dual-node by Mike Jarvis under its BSD-style license.

## Fast shared kernels and combined orders

This OpenMP wrapper enables the same qualified fast arithmetic used by the
spherical octree and KD-tree: direct normalization of squared spin-2
orientations, shared positive/negative ring products, and guarded lookup in
small logarithmic histograms. The lookup shortcut requires double precision,
at most 32 bins, `deltaR >= 1e-8`, `rminHist >= 1e-100`, and
`rangeN <= 1e100`. Edge neighborhoods within `1e-10*distance`, linear bins,
and fallback domains retain the original radial expression.

Exact combined runs with `no-one-ball`, `no-two-balls`, or `theta=0` now
collect pairs and 3PCF rings in the same body-pivot visits. Radial bin, bearing,
and transported neighbor shear are computed once; pivot weight and weighted
shear are cached per pivot. Node second moments still rotate by the squared
spin-2 rotation. Smoothed combined visits also share this arithmetic while
retaining the existing smoothing estimator and safety bound.

Approximate unsmoothed combined runs retain independent pair and ring
acceptance. Pivot reuse remains opt-in with its own phase budget; exact and
smoothing controls still select its body-pivot fallback. The existing PCA
builder, member cache, great-circle pair kernel, frontier scheduler, and
reduction order are unchanged. These private fast switches do not enable a
new path in the spherical ball-tree MPI wrapper. No opening or smoothing
control was relaxed, and floating-point roundoff can differ from older builds.

Run `make test-shear-sphere-balltree-2balls` for independent spin oracles,
log/linear bins and edge ULPs, masked one/two/three catalog roles, exact/reuse
fallbacks, thread determinism, serial/parallel construction, member-cache
fallback, near-polar pairs, and reuse calibration. Scientific/profile checks
in the ball-tree and shared reuse suites use explicit failures and remain
active when these tests are run directly with Python `-O`. The active release
gate now includes both ball-tree construction and reuse suites when this
engine is enabled.

`scripts/benchmark_shear_sphere.py` already supports this engine. Example:

```bash
python scripts/benchmark_shear_sphere.py \
  --engine balltree-shear-sphere-2balls-omp --output /tmp/ball-exact \
  --n 8192 --geometry clustered --order both --threads 4 --exact
python scripts/benchmark_shear_sphere.py \
  --engine balltree-shear-sphere-2balls-omp --output /tmp/ball-approx \
  --n 8192 --geometry clustered --order both --threads 4 --theta .1 \
  --reference /tmp/ball-exact/products-0.npz
```

For reuse, pass `--pivot-reuse` and set `CBALLS_SHEAR_PIVOT_TOL` explicitly.
Use `--nsmooth` to compare leaf capacities and `--smooth` for an explicit
arcminute radius. Match fixture, bins, modes, masks, and weights to the exact
reference; an NPZ alone does not certify provenance. The example 2% L2 plus
`1e-12` absolute policy is checked per raw/window/corrected observable and
finite mask. Record per-bin errors too. Passing an A/B implementation check
is separate from qualifying an approximation for a science application.

Runtime `theta`, `rsmooth`, and `nsmooth` keep their existing meanings.
`BALLS4SCANLEVON=1` retains adaptive task scheduling; compile-time `THETA`
remains a legacy scan-level parameter, not a percentage-error bound. The
active profile's default leaf capacity (8) and all accuracy defaults are
unchanged by this optimization.


### Hierarchical 3PCF reuse

All three OpenMP spherical two-ball shear engines support `shear-pivot-reuse`.
Partial multipole rings retain their original spherical acceptance frames.
An unresolved neighbor marks only the radial bins that its enclosing distance
interval can intersect. Once both legs of a radial pair are complete, that pair
is accumulated at the current pivot cell. Descendants inherit a completion mask
and never accumulate that pair again. Mixed resolved/unresolved pairs wait until
both legs are complete; diagonal self-neighbor subtraction is retained.

The option is off by default and requires `no-smooth-pivot`, positive `theta`,
full pivot coverage, and `BALLS4SCANLEVON=1`. Exact controls and smoothing retain
the body-pivot fallback. MPI shear methods continue to reject this experimental
option. The independent 2PCF pass retains its own acceptance rules.

- `CBALLS_SHEAR_PIVOT_TOL`: phase budget in radians, finite `[0,3]`, default `0.1`.
- `CBALLS_SHEAR_BIN_THETA`: internal radial-bin assignment allowance, finite `[0,1]`
  bin widths, default `0`. Zero requires complete containment inside a bin.
  Positive values allow center-based assignment when the combined cap radius
  fits within the selected fraction of the local bin width. The minimum and
  maximum separation cuts remain strict. This is not an exact translation of
  dual node's bin-slop rule.
- `nsmooth=1`: exposes finer pivot groups in the KD and PCA ball trees. It increases
  tree storage and can improve reuse; the octree already has individual-body leaves.

For a performance/accuracy trial through the Python drivers, add
`--nsmooth 1 --more-options shear-pivot-reuse` and explicitly set the two environment
controls, for example `CBALLS_SHEAR_PIVOT_TOL=3 CBALLS_SHEAR_BIN_THETA=0.5`.
These are approximate trial settings, not an accuracy guarantee. Compare raw
numerators, windows, and corrected complex multipoles against an exact reference
for the actual catalog. Poorly conditioned windows can amplify small raw errors.

Native provenance records `shear_hierarchical_reuse.enabled`,
`phase_budget_radians`, and `radial_bin_slop` from the completed run. With
`CBALLS_SHEAR_PROFILE=1`, `radial_pairs` counts actual radial-pair combinations,
`represented_pairs` counts the equivalent individual-pivot combinations, and
`partial_reductions` identifies completion above unresolved descendants. The
64 MiB per-worker scratch limit includes completion masks. Reuse walk timers
include reductions and must not be added to the reduction timer.
