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

With `dual-node-bin-slop`, two-ball acceptance uses the dual-node controlled
`radius1 + radius2 <= theta * bin_width` criterion and its conservative
near-boundary fallback. Without that option, the full interval
`distance - radius1 - radius2` through `distance + radius1 + radius2` must stay
inside one radial bin. Otherwise the larger node is split with dual-node's
`0.585` split rule. Acceptance constants are precomputed once, and a boundary
case reuses its logarithmic bin coordinate rather than evaluating a second log. The
3PCF also bounds angular phase error by `theta*pi/(2*mChebyshev+1)`. Use
`options=no-two-balls` for exact body-pair and body-moment accumulation.

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
