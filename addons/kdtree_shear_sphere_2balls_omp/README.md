# kdtree-shear-sphere-2balls-omp

This addon computes full-sky weak-lensing shear `xi+`, `xi-`, and all four
natural Gamma-x 3PCF multipole components over an independent median KD tree.
Input positions are observer-centered 3D vectors and shear components are
defined in each sample's local east/north frame.

Every KD node stores weighted first and second spin-2 moments transported into
the tangent frame at its spherical center. The 2PCF uses a symmetric
dual-node traversal, combined node radii, same-bin containment, and the 0.585
split heuristic. The 3PCF keeps exact body pivots and scans accepted neighbor
KD nodes into radial LogMultipole rings, including exact repeated-neighbor
second-moment subtraction.

Catalog positions are normalized once before the tree is built. Tree
aggregation and the hot transport kernel then reuse those unit vectors rather
than taking another square root per node membership or accepted neighbor.
Radial rings use a bin-major layout and store only nonnegative weight modes;
transport geometry is shared across the complex products, whose two-component
updates are SIMD-friendly. Ring buffers are reused and only their active spans
are cleared.

With `BALLS4SCANLEVON=1`, 3PCF work is split into a work-estimated pivot-cell
frontier rather than fixed pivot ranges. Tasks prefilter conservative neighbor
roots, run under dynamic OpenMP scheduling, and merge task-local histograms in
spatial order under a bounded memory budget. The independent 2PCF path uses a
symmetric dual-tree frontier. Set `CBALLS_SHEAR_PROFILE=1` to print wall and
per-thread time in radial lookup, transport, ring accumulation, ring clearing,
tree walking, and reduction.

Unsmoothed repeated catalog roles share one immutable KD tree within a
correlation call. Smoothed pivot and neighbor roles keep separate trees;
nothing is cached across calls or retained after model cleanup. Profiling also
reports `binary_tree_build` and the number of `unique_trees` actually built.

Tree construction uses OpenMP tasks for independent median-split branches
containing at least 8192 bodies each. Subtrees reserve their exact node ranges
before tasks start, preserving the serial depth-first layout and each node's
moment summation order. Smaller branches are built serially. This also applies
within each MPI rank and does not change approximation controls. Node storage
is allocated for the actual leaf capacity rather than always for body leaves.

The binary traversal rejects cells before logarithmic bin lookup or bearing
calculation when cheaper range and angular-extent tests suffice. The 2PCF
reuses pair distance for its split decision and prepares its fixed acceptance
settings once per traversal. Bin edges, opening tolerances, and histogram
publication order are unchanged. These optimizations apply to the shared
binary scan used by the spherical ball-tree and MPI partners as well.

Coincident cell centers are subdivided when their bounding spheres overlap the
search range. They are not treated as zero-distance body pairs: cross catalogs
can have identical cell centers while containing valid nonzero-distance pairs.

`only-2pcf` and `only-3pcf` skip the unused statistic. `no-two-balls` or
`no-one-ball` forces body-level results; `dual-node-bin-theta` enables the looser
dual-node radial criterion. Masks and the shared shear mode-coupling edge solve
are supported. With `SMOOTHPIVOTON=1`, deterministic transport-aware smoothing
is available through `smooth-pivot`; it is off by default and
`no-smooth-pivot` disables it explicitly. The safety requirement
`2*rsmooth <= rminHist` is inherited from the spherical shear estimator.

The two-node split strategy is adapted from dual-node by Mike Jarvis under its
BSD-style license.

## KD OpenMP fast kernels

The KD OpenMP wrapper enables guarded small-histogram lookup, direct
normalization of squared spin-2 orientations, and shared positive/negative
ring products. The logarithmic shortcut is limited to double precision,
at most 32 bins, `deltaR >= 1e-8`, `rminHist >= 1e-100`, and
`rangeN <= 1e100`. Distances within `1e-10*distance` of an edge, linear bins,
and other domains retain the original bin formula.

Standalone 2PCF projects both weighted shears onto the connecting great circle,
normalizes squared bearings directly, and updates auto-pair multiplicity once.
Cross pairs retain their orientation and imaginary xi+ component. This reuses
the geodesic pair implementation already available to the spherical ball-tree.

For unsmoothed exact combined runs, `no-one-ball`, `no-two-balls`, or `theta=0`
allows the pair contribution and 3PCF rings to share each body visit. The fused
walk evaluates radial bin, bearing, and spin transport once and caches the
pivot weight and weighted shear. It preserves distinct-neighbor subtraction
and node second moments. Approximate unsmoothed combined runs keep a separate
symmetric pair traversal because pair and ring acceptance differ. Smoothed
combined runs already use pivot traversal and now share its arithmetic too.

The KD OpenMP and MPI wrappers enable these specializations. They do not change
opening criteria, smoothing defaults, normalization, or the estimator.
Roundoff and summation order can differ from
older binaries. Thread reductions retain their deterministic order.

## Reproducible performance and qualification

Run `make test-shear-sphere-kdtree-2balls` for the independent spin reference calculation,
linear/logarithmic bins, poles, leaf capacities, masks, one/two/three catalog
roles, all three exact controls, order selection, smoothing, accepted cells,
and thread determinism. Parallel construction is tested on odd populations
with leaf capacities 1, 8 and 31, masks, zero weights and repeated cross-catalog
roles. The radial helper also receives over two million
classifications against its original expression, including edge ULPs and
fallback domains.

The benchmark driver records fixtures, build identity, settings, timings, and
raw/corrected observables for separate or combined orders. Example:

```bash
python scripts/benchmark_shear_sphere.py \
  --engine kdtree-shear-sphere-2balls-omp --output /tmp/kd-exact \
  --n 8192 --geometry clustered --order both --threads 4 --exact
python scripts/benchmark_shear_sphere.py \
  --engine kdtree-shear-sphere-2balls-omp --output /tmp/kd-approx \
  --n 8192 --geometry clustered --order both --threads 4 --theta .1 \
  --reference /tmp/kd-exact/products-0.npz
```

Use `--module` to measure a preserved extension in a separate process. Keep
fixture, bins, modes, weights, mask, threads, and parameters matched for A/B
implementation comparisons. Match the reference catalog and observable
settings before using `--reference`; the NPZ alone does not identify the catalog and settings used for the calculation.
The driver reports each observable's finite mask, scaled relative L2 error,
and absolute/bin errors. Its 2% L2 plus `1e-12` absolute criterion is an example
qualification policy, not a guarantee from the tree opening bound.

A DES-Y3 Takahashi calibration (r002/bin1/region1, one million sampled objects,
seed 8675309, unit weights, local east/north shear, 20 Sofia figure-1 bins,
signal order 7, window order 14, 16 threads) measured parallel construction
at 0.64 s versus 3.24 s for the serial builder. Four interleaved timings per
implementation averaged 22.70 s versus 24.00 s for the complete correlation
with leaf capacity 1, phase tolerance 3 and radial allowance 0.5. Mean CPU
seconds increased by 4.7%; all output arrays remained bitwise identical.

Against exact sums on that sample, leaf capacity 1 gave 2.62% raw and 2.19%
corrected global relative L2 error, while corrected order +/-7 reached 8.82%.
Leaf capacity 8 with the same tolerances took 23.99 s in one trial and gave
2.46% raw and 0.166% corrected global error; its largest corrected order error
was 0.621%. Changing leaf capacity changes approximate accepted interactions.
These are sample measurements, not full-catalogue timing or error guarantees.
The source change does not loosen any approximation controls.

`theta` controls runtime opening/phase accuracy; `nsmooth` is KD leaf capacity
(default 8 in the active unit-sphere profile). Smaller leaves can change the
speed and accepted approximate interactions, so qualify each setting.
`rsmooth` is a spherical radius in arcminutes; the internal chord bound remains
`2*rsmooth <= rminHist`. Smoothing is an approximation even with exact neighbor
visits. `BALLS4SCANLEVON=1` retains the adaptive pivot-cell scheduler.
Compile-time `THETA` remains a legacy scan-level parameter, not a replacement
for runtime error qualification. Experimental `shear-pivot-reuse` is supported; see the hierarchical controls below.
Existing control defaults remain unchanged.


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
the body-pivot fallback. The octree, KD-tree and ball-tree MPI shear methods also support this
option; controls must agree across ranks. The independent 2PCF pass retains its own acceptance rules.

- `CBALLS_SHEAR_PIVOT_TOL`: phase budget in radians, finite `[0,3]`, default `3`.
- `CBALLS_SHEAR_BIN_THETA`: internal radial-bin assignment allowance, finite `[0,1]`
  bin widths, default `0.5`. Zero requires complete containment inside a bin.
  Positive values allow center-based assignment when the combined cap radius
  fits within the selected fraction of the local bin width. The minimum and
  maximum separation cuts remain strict. This is not an exact translation of
  dual node's bin-slop rule.
- `nsmooth` defaults to `8` in native cballs and the active Python drivers.
  `nsmooth=1` exposes finer pivot groups in the KD and PCA ball trees. It increases
  tree storage and can improve reuse; the octree already has individual-body leaves.

For a performance/accuracy trial through the Python drivers, add
`--more-options shear-pivot-reuse` with the shipped defaults
`CBALLS_SHEAR_PIVOT_TOL=3 CBALLS_SHEAR_BIN_THETA=0.5`.
A tighter radial allowance such as `0.05` remains an explicit choice. Both defaults can be changed with `?=` assignments in
`addons/Makefile_addons_settings` or overridden by the runtime environment.
These are approximate trial settings, not an accuracy guarantee. Compare raw
numerators, windows, and corrected complex multipoles against an exact reference
for the actual catalog. Poorly conditioned windows can amplify small raw errors.

Native run information records `shear_hierarchical_reuse.enabled`,
`phase_budget_radians`, and `radial_bin_slop` from the completed run. With
`CBALLS_SHEAR_PROFILE=1`, `radial_pairs` counts actual radial-pair combinations,
`represented_pairs` counts the equivalent individual-pivot combinations, and
`partial_reductions` identifies completion above unresolved descendants. The
64 MiB per-worker scratch limit includes completion masks. Reuse walk timers
include reductions and must not be added to the reduction timer.
