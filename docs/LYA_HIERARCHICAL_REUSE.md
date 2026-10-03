# Hierarchical reuse in the 3D Ly-alpha engines

The six `lya-{2pcf,3pcf,2pcf-3pcf}-{omp,mpi}` methods support certified
forest-cell reuse. Select `lya2Kernel=1` for pairs and `lya3Kernel=5` for
triples. Combined runs share the forest hierarchy. Kernel-0 defaults remain
unchanged. No new Makefile or environment flag is needed.

## Estimator and controls

The 2PCF uses parallel/transverse separation bins. The 3PCF measures a
five-dimensional histogram: two radial bins, two
polar bins relative to the pivot LOS, and an opening-cosine bin. The new
path reuses `W=sum(w)` and `Q=sum(w*delta)` moments and certified radial/polar
combinations. It preserves this histogram; no finite Legendre reconstruction
is substituted. The separate anisotropic-multipole method has its own
[Legendre moment reuse](LYA_MULTIPOLE_HIERARCHICAL_REUSE.md); its finite-order
mu-bin reconstruction remains a different observable.

| Parameter | Selection |
|---|---|
| `lya2Kernel=0` | Original pixel pair traversal (default) |
| `lya2Kernel=1` | Persistent pair cells and certified radial range sums; OMP/MPI |
| `lya3Kernel=0` | Per-pivot segments (default), now built from child moments |
| `lya3Kernel=1,2` | Direct reference/tiled paths |
| `lya3Kernel=3,4` | Persistent/pivot cells, now also supported by MPI |
| `lya3Kernel=5` | Adaptive radial/polar moment hierarchy and pivot cells |
| `lya3PivotCellMax=8` | Maximum pivot task population for kernels 4/5; 1–64 |
| `lya3PivotBlock=0` | Automatic task block size; explicit values supported |

All geometry slops default to zero, preserving the estimator up to floating
summation order. Existing finite `[0,1]` pair and triple slops remain explicit
approximations in their own bin-width units. Radial/polar slop supports
kernels 3/4/5. Physical cuts stay strict. Slop is not a correlation-error bound.
The new fallback applies only at zero triple slop.

`theta`, compile-time `THETA`, `nsmooth`, and `rsmooth` are not translated into
forest error controls. The separate `lyaScanLevel`/`lyaPivotRadius` features
retain their pixel-pivot scope and reject persistent kernels and MPI.

## Reused work and correctness

Immutable per-forest trees retain geometry and weight/field sums. For a pivot
cell, kernel 5 groups the source frontier by radial/polar bin and azimuth,
then combines child moments and direction enclosures. Products are deposited
once when the opening-cosine interval fits a bin and source forests are
certified disjoint. Both leg orders and every represented pivot are included.
Unresolved products refine their own branches; completed contributions remain.

Forest-ID ranges and 256-bit hashed membership masks certify disjointness.
Hash collisions only cause refinement, never false acceptance. This avoids
subtracting a large same-forest auto product from a nearly equal total.
Bounds cover actual noncollinear coordinates. Unsafe geometry/exponent
ranges and non-double/fast-math profiles descend to retained pixel arithmetic.

At zero triple slop, if pivot tasks exceed one quarter of the active pixels,
kernel 5 falls back to exact segments before allocating geometry caches.
This fixed geometric decision is independent of threads/ranks. A completed
pair-cell pass is retained in combined mode. It avoids a large slowdown seen
on sparsely sampled forests but is not a universal speed guarantee.

The pair scan reuses observer-angle enclosures. An exponential/binary search
finds contiguous radial ranges of at least eight pixels certified to occupy
one `(rp,rt)` bin. A cheap spacing check retains direct tiles on sparse regular
sampling, where moment queries cost more than their short products. Dyadic
tree moments sum those ranges without prefix subtraction. Cartesian boxes
certify the discovery cut; ambiguous ranges retain pixel tests. The original
segment kernel now combines child sums/bounds instead of rescanning every
node's complete pixel range.

## MPI, memory and provenance

Catalogs, trees and histograms remain replicated. Pair cells use cyclic forest
rows; persistent triples use cyclic fixed pivot-task blocks. Pixel fallback
retains the existing per-pivot MPI ownership. OpenMP publication is ordered
within a rank. Raw sums and integer counts are reduced before normalization;
changing ranks can change last-bit summation. Only rank zero publishes.

Startup agreement covers kernels, bins, caps and slops. Rank-local errors use
collective failure boundaries. Checked memory plans include new scratch and
reallocation overlap; the budget is not a whole-process RSS cap. Persistent
geometry caches can be much larger than segment scratch.

`getRunMetadata()['lya_hierarchy']` reports kernel-5 selection, actual
`pixel_fallback`, hierarchy counts, represented products and certified pair
ranges. Counters are global on the publishing MPI rank. Existing qualification
metadata remains unmeasured until compared against a supplied reference.

## Use and reproduce

For combined runs use `lya2Kernel=1 lya3Kernel=5` with the appropriate native
search name, catalog, domains/bins and `options=no-smooth-pivot`. Set only the
applicable kernel for pair-only or triple-only runs. Pair-only methods reject
persistent triple kernels. MPI thread counts are **per rank**.

The all-engines driver accepts `--lya2-kernel 1 --lya3-kernel 5` and
`--lya3-pivot-cell-max 8`, applying them only to eligible 3D Ly-alpha statistics.
Radial-only and other multipole families retain their existing controls.

```sh
python3 scripts/benchmark_lya_hierarchy.py \
  --reference-module /path/to/preserved-before-extension \
  --geometry clustered --forests 64 --pixels 64 --threads 8 \
  --warmups 1 --repeats 3 --outdir results/lya-hierarchy
```

The benchmark retains raw arrays, occupancy/ratio errors, settings, hashes,
metadata, wall/CPU samples and rank RSS maxima. Use `--catalog` for NPZ arrays
`positions`, `delta`, `weights`, `forest_ids`; the existing driver can prepare
FITS selections. Select MPI names with `--methods`, `--ranks` and
`--mpi-command` (excluding `-n`); that Python needs mpi4py.

To isolate improvement over previous cell kernels, add
`--reference-parameters '{"lya2Kernel":1,"lya3Kernel":4}'`. The default
reference is kernel 0. Fresh processes alternate before/after order. Timings
cover native MainLoop; complete Python `Run` wall/CPU values are also retained.
Both exclude Python loading, catalog registration and result extraction.
Summed rank RSS maxima are not simultaneous peak memory. Failed exact
comparisons cause a failing exit while retaining evidence.

```sh
make test-lya-hierarchy
make test-lya-hierarchy-mpi
mpiexec -n 2 python3 -O tests/python/test_lya_hierarchy.py --cython-mpi-worker
```

The final command requires a rebuilt extension and mpi4py. Tests cover direct
oracles, adversarial geometry/weights, forest exclusions, strict cuts, actual
aggregation, thread/rank agreement, and rank-local failures with recovery.
Measure the intended catalog and bins before choosing a performance profile.

The LOS-tree variants now retain their LOS index in the sparse fallback and
include matching MPI registrations. Their per-pivot mixed-forest hierarchy and
block discovery reuse are described in [LOS reuse](LYA_LOS_HIERARCHICAL_REUSE.md).
