# LOS-tree hierarchical reuse

The six search methods are `lya-los-tree-{2pcf,3pcf,2pcf-3pcf}-{omp,mpi}`.
The OpenMP names require `LYAFORESTOMPON=1`; MPI names require
`LYAFORESTMPION=1`. Both use the shared OpenMP implementation. MPI IDs are
210, 211, 212; `print-search-methods` lists the compiled registrations.

## Estimator and exact reuse

These 3PCFs are five-dimensional hard-bin histograms of two radii, two polar
angles relative to the pivot LOS, and the opening cosine. A finite Legendre
expansion would change their estimator. Instead, certified nodes reuse
`W = sum(w)` and `Q = sum(w * delta)`. They deposit products of these moments
only when every represented leg and opening angle belongs to the same bin.
All three forest IDs must differ. Forest-ID ranges and hashed membership
masks certify disjointness; mask collisions force descent, never acceptance.
No same-forest total is subtracted from a large all-forest total.

Use `lya2Kernel=1 lya3Kernel=5` to opt into the existing persistent pair and
triple trees plus the new LOS fallback hierarchy. Defaults remain zero.

* Pair kernel 1 reuses forest-pair geometry and certified radial range moments.
* Triple kernel 5 uses persistent forest/pivot cells where compact pivots
  amortize geometry. This path was already available to LOS OpenMP methods.
* Sparse pivot groups now retain LOS-tree discovery. After building per-forest
  segments, 32 or more segment roots form a hierarchy sorted by radial/polar
  bin and azimuth about the pivot LOS. The LOS-aligned ordering keeps polar
  rings local; acceptance still uses the original Cartesian bounds.
  Certified mixed-forest products reuse the segment sums.
  Unresolved node products descend to the retained exact segment/pixel kernel.
  Smaller frontiers, or frontiers averaging more than three segment roots
  per forest, retain the existing segment loop. Many repeated forest IDs
  across leg bins otherwise prevent enough mixed-forest aggregation to
  amortize the additional hierarchy. This selection is geometry-dependent
  and independent of threads/ranks.
* LOS pixel paths share discovery within a fixed block when its pivots fit
  within a sphere of radius at most one quarter of the search cutoff. The
  enlarged sphere includes all pivot forests; the current pivot's own forest
  is excluded separately. Each pivot applies its own exact geometry. Wide
  blocks use separate discovery. Sorting candidate forests preserves thread
  determinism. This reuse also benefits default pair-only kernel 0.

Zero slops preserve the hard-bin estimator, with changed floating-point
summation order across algorithms. Positive existing slops remain explicit
bin-leakage allowances, not correlation-error bounds. Kernel 5's sparse
fallback is used only at zero triple slop. Extreme weights or unsafe geometric
exponents retain direct refinement. `theta`, legacy `SMOOTHPIVOT` and
`BALLS4SCANLEV` do not control these certificates.

## Running

```sh
./cballs search=lya-los-tree-2pcf-3pcf-omp \
  infile=forests.txt infileformat=lya-ascii iCatalogs=1 \
  numberThreads=8 usePeriodic=false lya2Kernel=1 lya3Kernel=5

mpiexec -n 2 ./cballs search=lya-los-tree-2pcf-3pcf-mpi \
  infile=forests.txt infileformat=lya-ascii iCatalogs=1 \
  numberThreads=4 usePeriodic=false lya2Kernel=1 lya3Kernel=5
```

Specify radial/polar/opening bins for the observable being measured. Use
`lya-los-tree-2pcf-*` or `lya-los-tree-3pcf-*` for a single statistic.
The Python driver accepts `--lya2-kernel 1 --lya3-kernel 5` and the new MPI
names, including `all-mpi` / `all-tree`. Existing spatial-frontier and positive
pivot-radius controls remain unsupported with MPI or persistent kernels.

MPI distributes fixed LOS pivot blocks, forest-pair rows or persistent pivot
tasks. Catalogs, trees and histograms remain replicated. Raw histograms and
counters reduce before rank 0 normalizes/publishes. MPI calls occur outside
OpenMP work regions. Allocation and computation failures enter the existing
collective error boundary. Results are repeatable across thread counts at a
fixed rank count; different rank counts may change rounding.

## Evidence and qualification

`run-metadata.json` / `getRunMetadata()` expose:

* `lya_hierarchy`: requested kernel 5, sparse fallback, hierarchy nodes,
  visited products, certified products, represented source pairs and pair ranges.
  Counts include both persistent and sparse LOS hierarchy work across ranks.
* `lya_los_tree`: global discovery traversals, reused pivot queries and exact
  pixel-distance tests. These remain zero when a persistent path completes
  without constructing the LOS index.

These work counters do not measure correlation accuracy. Numerical
qualification remains `UNMEASURED` until a reference comparison is supplied.

```sh
make test-lya-los-hierarchy
make test-lya-los-hierarchy-mpi
python3 scripts/benchmark_lya_hierarchy.py \
  --reference-module /path/to/pre-change/module \
  --methods lya-los-tree-2pcf-omp lya-los-tree-3pcf-omp lya-los-tree-2pcf-3pcf-omp \
  --threads 16 --geometry forests --forests 128 --pixels 64 \
  --warmups 1 --repeats 3 --outdir results/los-hierarchy
```

The benchmark retains raw arrays, numerical comparisons, native MainLoop and
complete Python Run timings, RSS and extension hashes. Use identical retained
catalogs, cuts, bins and slops. To isolate this implementation from the earlier
kernel selection benefit, add `--reference-parameters '{"lya2Kernel":1,
"lya3Kernel":5}'`. A pre-change build has no LOS MPI registrations; compare
current MPI against current OpenMP or against the current MPI pixel kernels,
and label that comparison explicitly.

The retained regression checks independent triangles, exact counts, changed
weights, bent/translated/boundary geometry, forest-mask collisions, combined
versus separate outputs, real hierarchy/discovery reuse, MPI partitioning,
rank-local failure and Cython recovery. It runs under Python `-O` using explicit
checks. Small tests are not a universal accuracy or scaling guarantee.

Speedup depends on geometry, bin occupancy and forest length. Hierarchy setup
and replicated memory can outweigh reuse in difficult geometries. Benchmark
representative data before choosing kernels or rank counts; multinode scaling
is a separate qualification.
