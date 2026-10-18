# Ly-alpha separation-plane quadtree: exact mode and optional bin slop

Use the existing engine name and options:

```bash
python3 tests/python/lya_corr_all_engines.py \
  --fits catalogs/desi-lya_downloaded/combined_lya_0500.fits \
  --fits-layout desi --engine lya-rprt-quadtree-2pcf-omp \
  --statistics 2pcf --threads 16 --max-forests 100000 --pixel-stride 2 \
  --lya2-kernel 1 --output Output_rprt
```

With the default zero slops, the engine computes the unordered cross-forest estimator
`N = sum(w_i delta_i w_j delta_j)`, `D = sum(w_i w_j)`, `xi = N/D`.
It preserves the observer-centered, bisector LOS geometry and exact requested
half-open rp/rt bins. Its separation plane depends on both endpoints; there is
no absolute-coordinate 2D pixel tree or materialized quadratic pair catalog.
The complete-forest bounds index remains the broad phase.

## Changes

1. Read chi, weight and weighted delta once into a compact, immutable cache in
   forest/radial order. Certified small products stream these three hot fields
   instead of repeatedly gathering large body records. Pair ownership and
   Cartesian discovery-sphere containment must be certified before this path.
   A bin-edge ambiguity immediately uses the original per-pixel calculation.
2. Reject whole-bin aggregation attempts before the quadtree lookup when the
   projected rectangle is wider than one bin or fails the outer containment
   certificate. The remaining lookups start from the inherited bin region.
   The four reference bin checks still guard every whole-cell acceptance.
3. Reuse an ancestor Cartesian containment certificate for descendants. Reuse
   an ancestor LOS enclosure while its uncertainty is small in bin units;
   otherwise recompute the enclosure. The 1/8-bin threshold controls when to
   tighten a conservative bound. It is not a bin slop or acceptance tolerance.
   Bent forests remain valid because the ancestor encloses all actual pixels.

Raw products retain the same arithmetic at individual accepted pixel pairs.
Whole-cell summation can change floating-point addition order, so comparison
must include raw N/D sums and pair counts, with appropriate roundoff tolerance.
Fixed ordered forest-row commits preserve repeatability across thread counts.

The compact cache costs 24 bytes per retained pixel in the double profile,
plus the existing bin tree and forest workspace. This allocation is included
in the checked, per-rank memory plan and freed on normal and failure exits.
Other forest engines do not allocate this cache and retain their search paths.

## Legacy controls

The useful ideas from the legacy octree/kdtree/balltree implementations are
hierarchical rejection, cached immutable bounds and bulk certified sums.
`SMOOTHPIVOT` substitutes grouped pivots, while `BALLS4SCANLEV` builds a legacy
scan frontier. Those structures do not implement this engine's pair-dependent
rp/rt bin plane. The forest index and radial hierarchy already provide its
frontier and complete-extent bounds.

`theta`, compile-time `THETA`, `rsmooth` and `nsmooth` do not tune this engine's
bin accuracy or compact tile capacity. Keep `no-smooth-pivot`; the Python
forest driver sets it. Forest pivot frontier controls remain rejected; pair
slops are now optional, explicit controls described below. Do not change these legacy values to seek
an uncalibrated speed-up. The cache optimization itself remains exact; positive pair slop explicitly
enables approximate bin placement.

## Validation and timing

Run the engine regressions against the pixel reference and independent small
catalog oracle:

```bash
python3 -m pytest -q tests/python/test_lya_rprt_quadtree.py \
  tests/python/test_lya_pair_cells.py tests/python/test_lya_forest_bounds_index.py
```

Keep an importable copy of the old cyballs extension in a separate directory.
The benchmark now accepts the quadtree engine directly:

```bash
python3 scripts/benchmark_lya_hierarchy.py \
  --reference-module /path/to/before --candidate-module . \
  --methods lya-rprt-quadtree-2pcf-omp \
  --reference-parameters '{"lya2Kernel":1}' \
  --candidate-parameters '{"lya2Kernel":1}' \
  --catalog /path/to/retained-input.npz --threads 16 \
  --radius 200 --pair-bins 50 --warmups 1 --repeats 5 \
  --outdir benchmark_rprt
```

Both extensions must use the same Python ABI, compiler and build settings.
The NPZ needs positions, delta, weights and integer forest_ids. This benchmark
excludes FITS decoding and Python catalog registration. It includes native
MainLoop tree/cache construction and search, alternates variant order, and
retains raw arrays, timings, peak process RSS, logs and provenance. It reports
native wall time and process CPU separately. Small millisecond cases are
sensitive to scheduling; use a quiet host and inspect every timing sample.
`verbose=2` also reports compact tile and certified cached-pair counts.

Reference methods reviewed: FCFC dual-tree range rejection and traversal,
Corrfunc sorted-coordinate pruning and data locality, dual node cell-pair
recursion, lya2pcf forest geometry and weighted-delta input, and the local
legacy octree, kdtree and balltree pair paths. No external source was copied.
A radial-window scan was also tried but did not improve the tested sparse
samples consistently, so it is not part of the installed change.


## Optional dual node-style bin slop

For the command in the introduction, add `--bin-slop 0.1` to set both axes.
`--lya2-bin-slop` is an alias. These are fractions of the **linear bin width**.
With 50 bins over 200 Mpc/h, 0.1 permits up to 0.4 Mpc/h leakage beyond each
edge of an assigned bin in each axis. It is not a 10% correlation-error bound.

Use `--lya2-rp-slop` and `--lya2-rt-slop` to override either axis, for example:

```bash
python3 tests/python/lya_corr_all_engines.py \
  --fits catalogs/desi-lya_downloaded/combined_lya_0500.fits \
  --fits-layout desi --engine lya-rprt-quadtree-2pcf-omp \
  --threads 16 --max-forests 100000 --pixel-stride 1 \
  --compare-desi-dr1 --lya2-kernel 1 --project-delta \
  --bin-slop 0.1 --output Output_rprt_projected_slop01
```

All slops must be finite in `[0,1]`; both default to zero. Keep the same
kernel 1: no new kernel number is needed. For native CLI/Cython, set
`lya2RpSlop` and `lya2RtSlop`. The Python driver also maps these controls to
other eligible 3D pair-cell engines when explicitly using pair kernel 1.
It resets them to zero for unrelated estimators when reusing a native object.

The exact single-bin acceptance runs first. With positive slop, a remaining
product may instead use a target leaf selected from the **geometric enclosure
midpoint**, provided its entire rp/rt rectangle lies in that expanded bin.
An axis with zero slop requires exact endpoint bin agreement. Individual
pixel products stay exact. Same-forest exclusions, original pivot ownership,
zero-weight pair counting and the original outer rp/rt and Cartesian discovery
cutoffs remain exact. No out-of-range pairs enter through expanded outer bins.

Numerators and denominators are summed as before, but assigned bins can change.
Global sums agree up to floating-point reduction order; occupied bins need not
match. Positive slop adds a direct leaf directory of `sizeof(size_t)*rp_bins*rt_bins`
bytes (20,000 bytes for 50 x 50 bins on this build), included in preflight.
The directory avoids walking the bin tree again for each expanded-bin query.
Metadata declares nonzero
slop approximate, and marks `lya_rprt_quadtree.exact_bin_containment=false`;
`outer_cuts_exact` remains true. Verbose work logs report `slop_products` and
`approximate_pairs`. This describes requested mode; a catalog may admit no
approximate products at small slop.

`--project-delta` runs before pair searching and remains compatible. The
continuum distortion matrix currently describes exact hard-bin assignment;
the driver rejects nonzero slop with distortion computation/forward modelling
rather than presenting that response as a response for the approximate bins.

Measure errors against the same projected zero-slop input. The driver can
record this comparison using `--accuracy-reference PATH_TO_EXACT_OUTPUT`;
`--accuracy-report-only` retains a report even when a target is exceeded.
Inspect absolute differences near zero and both N/D and xi, in addition to
relative L2 and per-bin relative errors away from zero. A small slop or a
subsample calibration does not establish an accuracy bound for a full catalog.

For a native calibration on a saved projected NPZ, the existing tool accepts
this engine:

```bash
python3 tests/python/benchmark_lya_pair_cells.py \
  --catalog /path/to/projected-input.npz --method lya-rprt-quadtree-2pcf-omp \
  --threads 16 --case 0:0 --case .1:.1 --case .25:.25 \
  --outdir benchmark_rprt_slop
```

That tool reports complete native-process wall time, including ASCII input.
Use the hierarchy benchmark above with `lya2RpSlop`/`lya2RtSlop` JSON overrides
for native MainLoop-only timing; its strict exactness gate will report an
expected numerical mismatch for genuinely approximate candidates and retains
the arrays for calibration. No speed or error improvement is assumed.


## Two-ball adaptation (2026-10-16)

The search refines two forest cells simultaneously. It now chooses one or both
cells using their radii projected into bin units, following the split-selection
idea in dual node (`src/Corr2.cpp`, `include/Split.h`) and the local
`support/native_octree_pair/native_octree_pair.c`. A partner more than twice as
small is kept intact; comparable cells split together when both exceed the
0.585 scheduling threshold relative to half the permitted bin footprint.
This threshold controls work scheduling only.

For a radial cell, let `m=(chi_min+chi_max)/2`, `h=(chi_max-chi_min)/2`.
For a pair of these balls, `H=h1+h2` bounds the radial separation error.
At fixed half-angle, the two projected errors are bounded by `H*cos(alpha/2)`
and `H*sin(alpha/2)`. This is the sum-of-radii argument in Jarvis, Bernstein &
Jain (2004), section 4.1, equations 62-63 and figure 3. The logarithmic
small-error approximation in that paper is not needed for our linear bins.

The implementation evaluates interval endpoints, with outward roundoff pads,
rather than subtracting rounded centers and radii. LOS angle enclosures account
for every actual pixel, including bent forests. Whole-product acceptance still
requires the complete projected rectangle inside an exact bin, or an explicitly
requested expanded bin, and inside the original rp/rt and Cartesian discovery
cutoffs. Isotropic distance alone cannot certify these two anisotropic bins.
Signed delta sums never define a position or radius.

For dense unresolved leaves, products of two small balls (at most two pixels
each) stream from the compact cache. Exact products use `(sum w*delta)_1 *
(sum w*delta)_2` and `(sum w)_1 * (sum w)_2`. Failed certificates retain
per-pixel arithmetic. These fine blocks themselves use zero bin slop; the
existing whole-cell slop rule above them is unchanged. A sparse-sampling gate
only chooses which exact path to run. It does not reject pairs. This avoids the
cost of recursive tiny-cell visits seen in an initial prototype.

Every source product is partitioned into disjoint child products. Forest
exclusions, original-ID pivot ownership, zero-weight pair counts and ordered
forest-row reductions remain intact. Moment regrouping changes roundoff, so
before/after raw arrays should agree within floating-point tolerance; repeated
runs and thread counts retain deterministic results.

No new command-line control is required: use the same engine, `--lya2-kernel 1`,
and zero slop for exact binning. Existing `--bin-slop`, per-axis overrides, and
`--project-delta` work as before. The implementation does not map legacy
`theta`, `THETA`, `rsmooth` or `nsmooth` to an unrequested approximation.
Counters distinguish tested/accepted two-ball blocks and represented pairs,
and one-cell versus two-cell subdivisions. They count work, not speedup.

The hot bin scales are precomputed once; extra storage consists of six scalar
constants per bin tree, two 16-entry index arrays and eight four-scalar ball
records on the tile stack. No
additional per-pixel allocation is required. Memory preflight for the retained
pixel cache and bin tree is unchanged.

For a new catalog, benchmark old/new modules with
`scripts/benchmark_lya_hierarchy.py`, supplying identical quadtree parameter
JSON to both `--reference-parameters` and `--candidate-parameters`, the same
retained catalog, threads and bin settings. Use at least five interleaved
repeats and one warm-up per process. Inspect raw numerator/denominator error
and pair counts as well as native wall time. The saved session report provides
actual DESI subset measurements; performance depends on sampling and bins.


## Parallel construction (2026-10-16)

For `lya-rprt-quadtree-2pcf-omp`, the existing `--threads N` now also controls
native preparation before the pair search:

- Sort retained pixel pointers in fixed 32768-pixel chunks. Parallel co-ranked
  merge slices use all workers even at the final merge level. The original
  forest/distance/ID comparator gives the same total order as the serial sort.
  This path starts at 65536 pixels.
- Reserve each forest's exact binary-tree preorder range, then build its
  radial hierarchy, complete Cartesian/LOS bounds and weighted moments in
  parallel. Forest construction starts at 32768 pixels and two forests.
  The 2PCF-only quadtree omits an unused 3PCF pivot-frontier traversal.
- Fill the compact chi/weight/weighted-field cache in parallel, starting at
  32768 pixels.
- Construct sufficiently large rp/rt bin quadtrees with OpenMP subtree tasks,
  starting at 16384 bins, with a minimum task size of 1024 bin leaves. Each
  subtree has deterministic compact preorder slots reserved before workers
  run. Odd, rectangular and singleton dimensions retain the same node layout,
  split coordinates and leaf directory as the serial tree.

The Rp/Rt coordinates of actual pixel pairs depend on both endpoints. Their
projection and two-ball certificates run inside the existing OpenMP forest-row
search. The quadtree's physical Rp/Rt cut coordinates are initialized by the
subtree tasks. Pair-space projection remains implicit and no pair catalog is
allocated. Small preparations use the serial path to amortize thread overhead.

No floating-point reductions are parallelized within a source cell. Its child
sum order, original-ID ownership, bins and strict cuts remain unchanged.
Only disjoint trees, merge output slices, cache entries and bin subtrees are
written concurrently. The pair search retains its ordered forest-row commits.

Parallel sorting temporarily allocates one pointer per retained pixel; its
combined peak with persistent forest storage and histogram plans is checked
before allocation. Scratch is freed before the bin tree and pair workers are
allocated. The existing node-capacity bounds are unchanged. Serial forest
boundary discovery and the small broad-phase forest-bounds index remain.

Verbose output now includes forest `build_wall`, `sort_wall`, `nodes_wall`,
`sort_threads`, `node_threads`, and `sort_scratch_bytes`; a separate rp/rt line
reports `pixel_wall`, `bin_wall`, `pixel_threads`, `bin_threads`, and
`projection_threads`. These are local-rank elapsed seconds and actual OpenMP
team sizes. CPU timings can rise while elapsed construction time falls.
Metadata records construction capability, implicit projection and the requested
construction mode; actual size-dependent team choices are in the log.

For a controlled native API benchmark, append `lya-serial-build` to `options`
to force serial preparation while leaving the pair-search thread count
unchanged. This option is for validation and measurement; normal Python driver
commands use automatic parallel preparation through `--threads`. Existing
zero/default and positive per-axis bin slops and `--project-delta` are retained.
