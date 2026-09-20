# Octree two-ball OpenMP addon

`search=octree-2balls-omp` applies dual-node-style scans to the native
cTreeBalls octree. It does not construct the FCFC PCA ball tree. A compact
binary view groups each octree cell's live children while retaining the native
octree hierarchy. The 2PCF uses a dual-node traversal. The production 3PCF
uses LogMultipole pivot-neighbor scans and forms every `(r1,r2,m)` bin from
products of radial moments, with second moments removing `q == r` exactly.
The octree adapter keeps pivots at exact body positions and applies two-ball
acceptance to neighbor nodes, avoiding failed coarse-pivot scans while still
reusing each accepted neighbor moment across all 3PCF radial-bin pairs.
For large multithreaded builds, the compact view first assigns deterministic
node and packed-point ranges, then computes leaf moments and parent bounds in
a post-order OpenMP task tree for subtrees above 32768 points. Single-thread
and smaller builds compute geometry inline to retain cache locality. This
parallelizes the expensive large-tree geometry pass without a concurrent
allocator or nondeterministic node numbering.

The process keeps the two most recently used compact views in a content-keyed
cache. Repeating a 2PCF over an unchanged catalog reuses its topology, bounds,
packed points, and scalar moments; changes to positions, masks, weights, or the
scalar field invalidate the entry. The cache also distinguishes adaptive leaf
capacities. Catalog fingerprints use fixed deterministic chunks and parallelize
above 262144 bodies, avoiding a serial validation floor on cache hits. A probe
before native tree construction lets a cache hit skip both the temporary native
octree and compact-view builds. Use
`options=no-native-tree-cache` for a cold-build benchmark or when retaining the
compact views is undesirable.

For 2PCF, `dual-node-bin-slop` enables controlled approximate radial-bin
acceptance with a `radius1 + radius2 <= theta * bin_width` criterion and a
near-boundary fallback. Without that option, the full interval
`distance - radius1 - radius2` through `distance + radius1 + radius2` must stay
inside one radial bin. Such same-bin scalar aggregation preserves the binned
pair counts and field products up to floating-point summation differences.
Otherwise the larger node is split with the dual-node
`0.585` split rule. Acceptance constants are precomputed once, and a boundary
case reuses its logarithmic bin coordinate rather than evaluating a second log. The
production 3PCF body-pivot path always applies theta-sized radial acceptance
and bounds angular phase error by `theta*pi/(2*mChebyshev+1)`; its radial
acceptance does not switch to exact bins when `dual-node-bin-slop` is absent.
With edge correction enabled, the angular bound uses the highest required
window order instead. Use `options=no-two-balls,no-smooth-pivot` for the
same-engine exact body-pair and body-moment reference.

## Accuracy calibration

The modern two-ball path has these parameter roles:

| Control | Role |
| --- | --- |
| `theta` | Controls allowed radial/angular node extent in production 3PCF and approximate 2PCF with `dual-node-bin-slop`. It is not a percentage bound on the final correlation. |
| `nsmooth` | Controls compact-tree leaf capacity, not pivot smoothing. It affects speed and can change the approximate node partition, so it must be included in a calibration profile. |
| `rsmooth` / `smooth-pivot` | Pivot smoothing is unsupported by the modern octree two-ball path. Changing a smoothing radius is not an accuracy/speed tuning strategy for this path. |
| Compile-time `THETA` | Controls the legacy B4 scan-table criterion. Modern octree two-ball traversal does not consume that table; uppercase `THETA` is not its acceptance tolerance. |
| `sizeHistN`, radial limits, `mChebyshev` | Define the observable and its resolution. Keep them fixed between exact and approximate runs. Lowering multipole order is a different truncation approximation. |
| `BALLS4SCANLEV`, thread count, tree cache | Affect scheduling/build work. Keep their configuration and the cache policy fixed for timing comparisons. |

Calibrate 2PCF and 3PCF separately on the same retained positions, scalar
values, weights, mask, bins, multipoles, and normalization as the reference.
In modern mode `no-one-ball` is not a substitute for `no-two-balls`.
Do not mix a `legacy-one-ball` reference with a modern-mode approximation.

For a strict 5% test, require every requested populated 2PCF bin and every
complex 3PCF coefficient to satisfy
`abs(approximate - exact) <= 0.05 * abs(exact) + roundoff_allowance`.
Use only a floating-point-scale allowance for zero coefficients, and reject
nonfinite outputs. An RMS error below 5% does not establish this maximum-error
condition, particularly for signed coefficients near zero. This test compares
complex magnitudes; it does not separately bound the relative errors of a
nearly zero real or imaginary part.

Sweep `theta` and leaf capacities rather than assuming their effect is
monotonic. Check candidate profiles on held-out regions and densities. Use
warmups and repeated timings; `no-native-tree-cache` makes every trial include
a fresh tree build. Revalidate after changing the catalog, field, mask, bins,
multipole order, normalization, edge-correction mode, or build. A calibration
on a subset is empirical evidence for that subset, not a guaranteed bound for
an arbitrary larger catalog. A tuned exact run can be preferable when the
strict error condition leaves little useful node aggregation.

## Completed-pivot progress

The production 3PCF reports completed body pivots in both exact and approximate
mode. Set `stepState=10000 verbose=1 verbose_log=1` to print progress to the
terminal and `cballs.log` after at least 10000 additional pivots have finished:

```text
octree-2balls-omp: 3PCF progress: completed pivots 10000 / 100000 (10.0%); elapsed 2.35 s
```

The count starts at zero, ends at the number of active pivots, and excludes
masked bodies when `read-mask` is enabled. OpenMP workers publish batches of
at most 64 pivots, so intermediate counts can pass the requested interval
slightly. Counts describe completed work, not catalog indices. Lines are
serialized and flushed immediately; `verbose=0 verbose_log=0` disables them.
Elapsed time measures the pivot traversal; the percentage is a pivot fraction,
not an estimate of remaining time. At 100 percent, histogram reduction,
normalization, any requested edge correction, and file output still follow.
The independent `only-2pcf` and `dual-node-direct-triples` paths do not use this
pivot counter.

## Numerical and parallel behavior

For 3D catalogs the angular error bound is on projected tangent bearings.
Catalogs retain the observer origin and radial bins are chord distances.
`no-normalize-HistZeta,weights-norm` selects raw weighted distinct triplets;
do not subtract repeated neighbors again. See [the scalar contract](../../docs/3pcf.rst).

`TWOPCFON=1` and `TPCFON=1` enable the two correlation orders at build time.
When both are active, `only-2pcf` and `only-3pcf` select one at runtime.
`compute-HistN`, `weights-norm`, `no-normalize-HistZeta`, and
`out-m-HistZeta` have the same meanings as for the active ball-tree 2-ball estimator.

The fixed pivot frontier owns private moment and histogram buffers and is
reduced in a fixed order, so OpenMP worker count does not alter the result.
Within each body-pivot task, the 3PCF now passes an inherited neighbor frontier
to child pivots. Its retention policy is adaptive: when splitting would exceed
the candidate budget, it retains the unresolved parent cell for a descendant
to refine instead of abandoning the list and restarting at the catalog root.
The first body in a dense leaf measures the local root cost before the
remaining bodies are allowed to reuse the inherited list. This keeps frontier
traffic bounded while preserving every disjoint candidate subtree.

The independent 2PCF path uses one dense histogram per OpenMP thread, reuses
its body-pair batch buffer across dynamic frontier tasks, and combines thread
histograms with a binary tree reduction. Its approximate-mode leaf capacity
uses catalog surface density, the outer radial-bin width, and the selected
`theta` error budget, with a calibrated 75--125 percent bound around
`nsmooth`; exact mode retains the configured `nsmooth` capacity. On macOS,
logarithmic body bins use
Accelerate vForce. On Linux, `NATIVE_PAIR_VECTOR_LOG=auto` enables the batched
path when `pkg-config sleef` is available. Use
`NATIVE_PAIR_VECTOR_LOG=sleef` to require that backend, or
`NATIVE_PAIR_VECTOR_LOG=compiler` with a toolchain configured for vector libm;
otherwise Linux retains the scalar logarithm path.

`options=dual-node-direct-triples` retains the cubic triple-node traversal as a
validation oracle for moderate catalogs. The traversal is adapted from
dual-node by Mike Jarvis under its BSD license; the full notice is in
`addons/balltree_2balls_omp/dual-node_LICENSE`.

## Octree-GGG compatibility

`options=legacy-one-ball` switches this search name to the privately linked
one-ball native-octree implementation. This is a real kernel dispatch: tree loading
builds the threaded native octree and GGG scan frontier instead of the compact
binary view. Consequently `behavior-ball`/`no-one-ball`, `compute-HistN`,
masking, normalization, complex edge correction, `ggg-full-window`, and
`ggg-profile` have exactly their GGG meanings. With `SMOOTHPIVOTON=1`, pivot
smoothing is enabled by default in compatibility mode and
`options=no-smooth-pivot` disables it.

Do not combine `legacy-one-ball` with `no-two-balls`, `dual-node-bin-slop`, or
`dual-node-direct-triples`; those select features of the native two-ball
kernel. `only-2pcf` is supported. `only-3pcf` is rejected because the GGG
kernel does not yet provide a true skip-2PCF execution path. Without
`legacy-one-ball`, smooth-pivot remains unsupported and the compact two-ball
algorithm is unchanged.
