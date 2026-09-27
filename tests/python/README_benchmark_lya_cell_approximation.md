# Ly-alpha cell-geometry calibration

`benchmark_lya_cell_approximation.py` measures the anisotropic, five-dimensional
forest 3PCF, preserving the observer LOS convention, weighted numerator and
denominator, and three distinct forest IDs. It uses the same public FITS reader
as `lya_corr_all_engines.py`. It supports Cartesian DESI examples as well as the
reader's forest FITS formats, normalized NPZ, six-column `lya-ascii`, and synthetic
input. This is a native-process benchmark; Python preprocessing is excluded.

## Controls and implementation order

1. `lya3Kernel=0 lya3MuSlop=s` reuses the existing per-pivot, same-forest,
   same-radial/polar-bin segments. Only the opening cosine can change bins.
2. `lya3Kernel=3` builds a persistent binary hierarchy within each forest once
   per search, sorted by observer distance. Pivots remain individual pixels.
3. `lya3Kernel=4` also starts from pivot cells with at most
   `lya3PivotCellMax=8` pixels (configurable 1..64). A scheduling heuristic also limits the cell
   bounding radius to one eighth of a radial-bin width, avoiding large pivot
   cells that almost always need subdivision. Unresolved triple products
   select a splittable cell using its estimated radial, polar and opening-angle
   bin uncertainty; the pivot score includes both legs and LOS variation.
   Already accepted products are retained.

All three slops default to zero. `lya3MuSlop`, `lya3RadialSlop`, and
`lya3PolarSlop` each accept a finite number from 0 to 1. They are fractions of
**their own linear bin widths**: `2/lya3MuBins`, `lya3RMax/lya3RBins`, and
`pi/lya3ThetaBins`. A bounded interval may extend that fraction beyond its
representative cell-center bin. Positive slop is an explicit approximation.
These parameters are **not** bounds on relative correlation error.

Radial/polar slop requires kernel 3 or 4. Kernels 1 (independent direct reference)
and 2 (tiled direct) reject positive slop. Kernels 3/4 currently support the four
OpenMP 3PCF/combined names only; MPI selection fails with an explicit message.
The original and `lya-los-tree-*` names share the persistent backend when kernel
3/4 is requested, so their 3PCF traversal is then the same. Kernels 0..2 retain
the existing engine-specific discovery path. Pair-only methods reject active
3PCF approximation settings. The separate experimental multipole method is not
combined with these controls.

## What remains exact

Every accepted cell is from one forest. All three forest IDs must differ. Both
pivot-neighbor legs must satisfy `0 < r < lya3RMax`; slop never expands this
physical domain. With zero slop a cell product must fit one bin on every axis,
otherwise it is subdivided until reference pixel arithmetic resolves it.
Floating-point sums can differ because contributions are reassociated.

Each node stores `W=sum(w)` and `Q=sum(w*delta)`. Certified products deposit
`N=Qp*Qq*Qr` and `D=Wp*Wq*Wr` for both leg orders. Cell positions are geometric
bounding-box centers, never signed-field centroids. Cartesian boxes are intersected
with forest-aligned capsule bounds: the segment joining a node's endpoint pixels,
expanded by a radius enclosing every actual pixel in that node. Conservative
radial and directional-cone bounds include the pivot's extent. Cached observer-
direction extrema include the changing pivot LOS. Bent/noncollinear forests do
not assume an ideal ray; their capsule simply becomes wider. Unsafe exponent
ranges fall back to pixels.

Combined runs compute their 2PCF with the existing exact pixel traversal; these
slops affect only 3PCF. Forest nodes are private to one search and freed afterward.
Point pivots reject unreachable subtrees before computing their pixel geometry,
then cache actual leaf geometry and bottom-up extrema on demand. Pivot-cell
subdivisions reuse a bounded 8,192-entry cache keyed by both node IDs; collisions
only cause recomputation. Direct tiles use dedicated pixel-cache slots for the
active pivot task, so subdivision cannot repeatedly evict those pixel geometries.
The pixel record retains full-precision displacement, radius, polar angle and
bin/status fields; interval fields are reconstructed only for node traversal.
On a 64-bit double build this is 56 bytes plus an 8-byte tag per slot.
Singleton pivots reuse their dense node cache. All caches are private to a worker
and discarded at the end of the search.

Dedicated pixel storage scales with `threads * pixels * largest_actual_pivot_task`;
its size is zero when every task is a singleton. The node hash and dense node
cache add per-worker storage. Checked arithmetic and the combined
histogram+nodes+geometry-cache+frontier memory preflight guard these allocations.
A larger pivot cap may therefore increase memory and is not a speed guarantee.
Task blocks and reduction order are fixed independently of thread count.
`lya3PivotBlock` counts pivot-cell tasks for kernels 3/4 (automatic block size 8),
and pixel pivots for kernels 0..2. `accepted`/`nbbcalc` count pixel-tree visits;
pure persistent 3PCF does not make such visits, so use its separate cell counters.

## Traversal diagnostics

The native log reports `evaluations`, `leaf_evaluations`, `cache_hits`,
`pair_cache_hits`, and `pruned_nodes`. The benchmark retains them per run in
`summary.json`, together with the existing aggregation counters. `evaluations`
counts cache misses that construct node or pixel geometry; `leaf_evaluations`
counts actual pixel geometry computations. `pair_cache_hits` is the internal-node
hash subset of all cache hits. `pruned_nodes` counts point-pivot subtree rejections before
pixel-cache construction. These are work counters, not elapsed phase timings.
Compare complete-process wall time, process CPU, peak RSS, counts and raw sums
alongside the counters.

At zero slop these refinements change traversal and summation order, not the
estimator. Positive-slop bin assignments may change when tighter bounds or new
subdivision decisions accept a different cell product; recalibrate approximate
settings after updating. There is no new automatic 5% error guarantee.

## Run a calibration

From the cTreeBalls source root, with the benchmark Python environment active:

```bash
python tests/python/benchmark_lya_cell_approximation.py \
  --fits /path/to/catalogs/lya_15_xyz_raw_with_losid.fits \
  --max-forests 1000 --pixel-stride 10 \
  --threads 16 --warmups 1 --repeats 3 \
  --reference-kernel 0 \
  --case 0:.001:0:0 --case 0:.01:0:0 --case 0:.05:0:0 \
  --case 3:.01:0:0 --case 4:.01:0:0 \
  --case 4:.01:.001:.001 --case 4:.01:.01:.01 \
  --relative-floor 1e-6 --max-relative-error .05 \
  --small-signal-atol 5e-8 --require-accepted \
  --outdir results/lya-cell-calibration
```

Use a fresh output directory. A case is
`kernel:mu-slop:radial-slop:polar-slop[:pivot-cap]`. The default pivot cap is 8.
Start with `--synthetic` or a representative small FITS selection and
`--reference-kernel 1` to compare directly with the independent kernel.
Kernel 0 is the retained exact segment reference for expensive full selections.
Use `--method lya-los-tree-3pcf-omp` for that named engine, or either combined
OpenMP name to additionally check that 2PCF agrees. The combined 2PCF uses its
native default domain, recorded in each run's metadata.

The script retains all native products, metadata and commands, hashes of the
executable and selected catalog, repeat timings, CPU time, peak native RSS, cell
counters, and per-bin raw sums/error tables. `summary.json` records acceptance:

- At most 5% relative zeta error in every occupied reference bin with
  `abs(zeta_exact) > relative-floor` (limits are configurable).
- Absolute zeta error at most `small-signal-atol` in every remaining occupied
  reference bin; the default is an explicit policy choice, not a physical law.
- No positive denominator in an empty reference bin.
- Identical physical ordered-triplet count and finite outputs.

Numerator and denominator relative L2 errors and percentile correlation errors
are also reported, but do not replace the per-bin acceptance conditions. Exact
cases additionally require tight agreement of raw sums. Repeated runs must
produce identical tables. `--require-accepted` returns a failure status if no
positive-slop case passes; without it, failed calibrations remain reportable.
A setting must be calibrated again when catalog, weights, bins, or scales change.
Choose among passing cases by timing; persistent cells can be slower when few
products can be aggregated. No default approximation or speedup is assumed.

## Files and diagnostics

Implementation: `addons/lya_forest_omp/lya_triplet_exact.h` (mu-only slop),
`lya_triplet_cells.h` (persistent and pivot cells), and
`search_lya_forest_omp.c` (dispatch/reduction/output). Parameters pass through the
shared native/Cython parser. No new histogram getter is introduced: the public
5D file remains the correlation product. Its comment header and the
`lya_geometry` run-metadata object record all controls and approximation status.
`getRunMetadata()` exposes the latter through Cython.

`approximate_pairs` counts pivot/neighbor-pair combinations accepted using a
slop-relaxed bound; `ordered_triplets` counts both leg orders.
`aggregated_pairs` includes exact and approximate aggregated work;
`pivot_aggregates` counts accepted node products with more than one pivot.
These diagnostics are not a count of actually misplaced triplets.

Regression module: `tests/make_tests/test_lya_cell_approximation.py`. Run it with
pytest and `CBALLS=/absolute/path/to/cballs`. It is also selected by the capability
manifest and active release gate.
