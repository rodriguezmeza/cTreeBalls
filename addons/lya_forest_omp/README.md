# Lyman-alpha forest OpenMP addon

This addon computes weighted Lyman-alpha forest correlations from individual
forest pixels. It uses the cTreeBalls octree for exact range rejection; accepted
leaves are always evaluated as individual pixels.

## Build and methods

Set `LYAFORESTOMPON=1` in `addons/Makefile_addons_settings`, then build normally.
The available search methods are:

- `lya-2pcf-omp`: weighted 2PCF only.
- `lya-3pcf-omp`: weighted anisotropic 3PCF only.
- `lya-2pcf-3pcf-omp`: both estimators in one tree traversal.
- `lya-los-tree-2pcf-omp`: exact 3D 2PCF using octree forest discovery and
  per-forest radial trees.
- `lya-los-tree-3pcf-omp`: the same hybrid search for the exact 5D 3PCF.
- `lya-los-tree-2pcf-3pcf-omp`: both 3D estimators with shared forest discovery.
- `lya-1d-2pcf-omp`: weighted radial-only 2PCF.
- `lya-1d-tree-2pcf-omp`: exact interval-tree radial-only 2PCF.
- `lya-1d-tree-same-los-2pcf-omp`: equal-LOS average of exact within-forest
  radial 2PCFs.
- `lya-1d-tree-3pcf-omp`: exact interval-tree radial-only 3PCF.
- `lya-1d-3pcf-omp`: weighted radial-only 3PCF.
- `lya-1d-2pcf-3pcf-omp`: both radial estimators in one scan.

All methods require `DEFDIMENSION=3`, OpenMP, `usePeriodic=false`, exactly one
input file, and `infileformat=lya-ascii`.

Set `LYAFORESTMPION=1` for an MPI+OpenMP counterpart of every method above
except `lya-1d-tree-same-los-2pcf-omp` and the three `lya-los-tree-*` methods,
replacing the final `-omp` with `-mpi`.
See
[`addons/lya_forest_mpi/README.md`](../lya_forest_mpi/README.md) for the
replicated-catalog contract, example parameters, and rank-comparison tests.

## Octree Discovery and LOS Trees

The `lya-los-tree-*` methods implement a **1D search within each forest, not
a radial-only estimator**. They retain transverse distances, the existing
anisotropic bins, statistical weights, and all distinct-forest exclusions.
Their output formats match `lya-2pcf-omp` and `lya-3pcf-omp`.

![Forest discovery followed by radial queries](los_tree_schematic.png)

The schematic is illustrative, not to scale; the equations below define the
actual conservative interval used by the implementation.

1. Build the native 3D octree, a compact threaded view tagged with homogeneous
   forest IDs, and balanced radial interval trees for all forests. Each radial
   leaf holds at most eight pixels. Trees are immutable during the search.
2. For each pivot, walk the octree until one valid pixel inside the search
   sphere witnesses a forest. Record that forest once. Skip homogeneous
   branches of the pivot forest and already-discovered forests. Mixed-forest
   branches still require traversal; discovery cost is not always one test per
   forest. Thread-local generation stamps avoid clearing the full forest table
   for every pivot.
3. For a straight forest with unit sightline `n_f`, pivot observer-position
   `x_p`, and search radius `R`, use
   `a = dot(x_p, n_f)`, `b^2 = |x_p|^2 - a^2`, and
   `chi_min,max = a +/- sqrt(R^2 - b^2)`. The radial coordinate `chi` is
   measured from the observer, not from the pivot. The forest tree automatically
   clips this interval to the pixel support.
4. Query each discovered forest's radial tree, apply the exact 3D distance cut
   to each candidate pixel, and feed the original 2PCF/3PCF accumulators.

The 2PCF discovery radius is `hypot(lya2RpMax, lya2RtMax)`, which encloses the
full rectangular pair domain. It is **not** either axis maximum alone. The
3PCF uses `lya3RMax`; combined mode uses the larger radius and then applies
each estimator's own cuts. The pivot forest is entirely excluded. For 3PCF,
the two neighbor forests must also differ. Ordered pivot-side permutations
and the existing `mu = cos(opening angle)` bins are unchanged.

Interval calculations use a conservatively inflated sphere. A measured bound
on each forest's departure from its reference sightline covers non-collinear
pixels, recentering, and mixed-precision position storage; final distance/bin
tests use the original stored geometry. No `theta`-controlled approximation
or smooth-pivot aggregation is introduced. Arithmetic order differs from the
old traversal, so agreement is to floating-point rounding, not necessarily
byte-for-byte between different algorithms. Fixed 64-pivot publication blocks
preserve determinism across OpenMP thread counts for the same method/build,
while avoiding a synchronized merge and scratch clearing after every pivot.
The old 3D methods keep their original per-pivot publication order.

This accelerates neighbor discovery, not the quadratic neighbor-pair loop
of the five-dimensional 3PCF. Dense 3PCF neighborhoods can therefore remain
dominated by triangle accumulation. Index memory is linear in pixels and
native octree nodes, with `O(forests * threads)` discovery scratch. Trees are
rebuilt per correlation call and do not retain Python pointers after cleanup.
With `verbose=2`, `LOS-tree:` reports index-build CPU time and discovery/radial
visit counters. `build_CPU` excludes the preceding native octree construction.

Example using a six-column forest catalog in comoving units:

```sh
./cballs search=lya-los-tree-2pcf-3pcf-omp \
  infile=forests.txt infileformat=lya-ascii iCatalogs=1 \
  usePeriodic=false numberThreads=8 rootDir=Output_lya_los_tree \
  lya2RpMax=200 lya2RtMax=200 lya2RpBins=50 lya2RtBins=50 \
  lya3RMax=80 lya3RBins=8 lya3ThetaBins=6 lya3MuBins=8
```

For only one statistic, change the search name to `lya-los-tree-2pcf-omp`
or `lya-los-tree-3pcf-omp`. The old `lya-1d-tree-*` methods remain unchanged
and measure different, radial-only statistics.

Run `make test-lya-forest-los-tree` (requires pytest and NumPy) for independent
oracle checks, old/new 3D histogram equality, forest exclusion, near-boundary
and non-collinear geometry, combined/separate modes, and thread determinism.

## Input

The multi-engine driver `tests/python/lya_corr_all_engines.py` reads DESI DR1 delta
FITS, NPZ, or this ASCII format once and retains the catalog across engines.
Use `cyballs.cballs.set_forest_catalog(positions, delta, weights, forest_ids)`
for direct NumPy input. Equal integer IDs identify the same quasar; observer
distances and sightlines are initialized before tree construction. The driver
broadcasts the retained arrays once for MPI. See
[`tests/python/README_lya_corr_all_engines.md`](../../tests/python/README_lya_corr_all_engines.md)
for a small public DESI download and run examples.

The ASCII interchange format has six columns:

```text
x y z delta weight forest_id
```

Coordinates are observer-centered comoving Cartesian coordinates. `forest_id`
identifies the quasar sightline. The standard forest estimators exclude pairs
and triplets containing two pixels from the same forest. The explicitly named
`same-los` method is the exception: it accepts only pairs with equal IDs.

Convert one or more lya2pcf `data*.npy` files with:

```bash
python3 addons/lya_forest_omp/lya2pcf_to_cballs.py \
  --lya2pcf-source /path/to/lya2pcf \
  -o forests.txt /path/to/data*.npy
```

The converter writes `delta = dw / we`, so the weighted delta stored by lya2pcf
is preserved exactly. Unpickling requires the Python dependencies used by the
lya2pcf checkout.

## Estimators

The 2PCF accumulates `sum(w_i delta_i w_j delta_j) / sum(w_i w_j)` in
`(r_parallel, r_transverse)` bins. Configure it with `lya2RpMax`, `lya2RtMax`,
`lya2RpBins`, and `lya2RtBins`.

The 3PCF accumulates `sum(w_i delta_i w_j delta_j w_k delta_k) /
sum(w_i w_j w_k)` for three distinct forests. Its five coordinates are the two
pivot-side lengths, their two angles to the pivot line of sight, and the cosine
of their opening angle. Configure it with `lya3RMax`, `lya3RBins`,
`lya3ThetaBins`, and `lya3MuBins`.

The radial methods build a sorted index of the observer-centered distance and
ignore the angular coordinates during candidate selection and binning. The
radial 2PCF uses `|chi_j-chi_i|` with `lya2RpMax` and `lya2RpBins`. The radial
3PCF uses the two signed pivot lags `(chi_j-chi_i, chi_k-chi_i)`. Each axis has
`2*lya3RBins` bins spanning `(-lya3RMax, lya3RMax)`, so configurations on the
same and opposite sides of a pivot remain distinguishable. Pixels in the same
quasar forest are excluded from both estimators; 3PCF triplets require three
distinct forest IDs. The build remains `DEFDIMENSION=3` because that describes
catalog storage, while the radial search itself is one-dimensional.

OpenMP work is grouped into fixed radial-pivot blocks. Every block processes
its pivots in radial order, and completed block histograms are merged in block
order. This makes the floating-point result independent of the OpenMP thread
count while avoiding a serialized merge after every pivot. Configure the
compile-time block size with `LYA1D_OMP_PIVOT_BLOCK_SIZE` in
`addons/Makefile_addons_settings`; the default is 256. Smaller values expose
more parallel tasks for small or strongly imbalanced catalogs, while larger
values reduce merge synchronization for inexpensive 2PCF workloads.

The radial 3PCF tree leaf capacity is configured independently with
`LYA1D_TREE3_LEAF_SIZE` in the same file; the default is 8. Smaller leaves
increase traversal depth and exact-node resolution, while larger leaves reduce
tree overhead but leave more work to the leaf kernel. Both values appear in
`cballs options=make-info` so benchmark profiles are reproducible.

`lya-1d-tree-2pcf-omp` is an independent, exact alternative for dense radial
2PCF workloads. It builds a balanced binary interval tree over the same radial
ordering. A node pair is accumulated in bulk only when its minimum and maximum
possible separations lie in the same `lya2RpBins` bin; otherwise the wider node
is split. The tree first accumulates all unordered pairs, then subtracts pairs
within each `forest_id` using a sorted per-forest range scan. Long-double
internal sums and exact per-bin pair counts protect this subtraction against
false occupied bins. Tree tasks and histogram commits use fixed logical order,
so output is byte-identical across OpenMP thread counts. The ordinary
`lya-1d-2pcf-omp` range scan remains the preferred baseline for sparse catalogs
or very narrow radial windows.

`lya-1d-tree-same-los-2pcf-omp` builds an independent balanced interval tree
for every `forest_id` and traverses only `(tree_f, tree_f)`. It therefore never
forms a cross-forest pair. For each radial bin and forest it computes
`xi_f = sum(w_i delta_i w_j delta_j) / sum(w_i w_j)`, then publishes the
unweighted mean of `xi_f` over forests whose denominator is nonzero in that
bin. This gives each occupied sightline equal weight rather than weighting the
answer by its pixel-pair count. Forests are committed in sorted ID order, so
the output is byte-identical across OpenMP thread counts.

`lya-1d-tree-3pcf-omp` computes the same signed-lag matrix as
`lya-1d-3pcf-omp`, but avoids its quadratic neighbor-pair loop. For each pivot,
the interval tree bulk-accumulates weight and weighted-delta moments only when
a complete node lies in one signed lag bin. Their outer products represent all
ordered neighbor pairs. The kernel then removes the pivot forest, self-pairs,
and all pairs whose two neighbors share a forest ID. The latter correction uses
a sparse per-worker forest/bin table, so memory does not scale as the full
number of forests times the number of bins. Fixed pivot blocks make output
byte-identical across OpenMP thread counts. This method is usually preferable
to the range-scan 3PCF when radial windows contain many pixels.

Outputs are `<histXi2pcfFileName>_lya.txt` and
`<histZetaFileName>M_lya5d.txt`. The 3PCF output is sparse by default; add
`lya-output-empty-bins` to `options` to include zero-weight bins. Both
estimators define a zero-denominator bin to have correlation zero.

Radial outputs are `<histXi2pcfFileName>_lya1d.txt` and
`<histZetaFileName>M_lya1d.txt`. They use the same zero-denominator and sparse
3PCF policies.

The same-LOS output is `<histXi2pcfFileName>_lya1d_same_los.txt`. Its five
columns are `bin radial_separation xi sum_xi contributing_los`; consequently
`xi = sum_xi / contributing_los`, with zero used when no LOS contributes.

## Tests

The regular regression uses independent NumPy 2PCF and 3PCF oracles and checks
OpenMP determinism:

```bash
make LYAFORESTOMPON=1 test-lya-forest-omp
```

The radial regression additionally checks direct Python 2PCF and 3PCF oracles,
the scan and tree implementations against one another, same-quasar subtraction,
the same-LOS equal-forest estimator, invariance under arbitrary changes to
transverse coordinates, and one-thread versus multi-thread byte-for-byte
repeatability:

```bash
make LYAFORESTOMPON=1 test-lya-forest-1d-omp
```

An additional end-to-end test runs the original lya2pcf CPU implementation on
real `forest_class.quasar` objects, converts the same saved dictionary, and
compares every cTreeBalls 2PCF numerator, denominator, and normalized bin:

```bash
make LYAFORESTOMPON=1 \
  LYA2PCF_SOURCE=/path/to/lya2pcf \
  test-lya2pcf-reference
```

This reference checkout implements only the 2PCF computation. The independent
oracle remains the regression reference for the addon 3PCF.

Compare the native radial estimators, the radial interval tree, and the 3D
octree estimators on the same exactly collinear catalog with:

```bash
python3 examples/compare_lya_1d_3d.py --threads 4
```

The script projects the 3D numerator and denominator onto the 1D bins before
normalization, checks the dedicated and combined methods, records wall times,
and writes CSV, JSON, and PNG results.  The companion notebook is
`examples/compare_lya_1d_3d.ipynb`.  Exact numerical agreement is meaningful
only for a collinear catalog; on a general survey the 1D and 3D methods apply
different transverse-selection rules.

Per-thread 3PCF storage is approximately
`lya3RBins^2 * lya3ThetaBins^2 * lya3MuBins *
(2*sizeof(real) + sizeof(size_t))` bytes. This
addon computes raw estimators; survey distortion matrices, covariance
estimation, continuum fitting, and metal corrections remain preprocessing or
post-processing responsibilities.
