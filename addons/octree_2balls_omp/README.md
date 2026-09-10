# Octree two-ball OpenMP addon

`search=octree-2balls-omp` applies dual-node-style scans to the native
cTreeBalls octree. It does not construct the FCFC PCA ball tree. A temporary
binary view groups each octree cell's live children while retaining the native
octree hierarchy. The 2PCF uses a dual-node traversal. The production 3PCF
uses LogMultipole pivot-neighbor scans and forms every `(r1,r2,m)` bin from
products of radial moments, with second moments removing `q == r` exactly.
The octree adapter keeps pivots at exact body positions and applies two-ball
acceptance to neighbor nodes, avoiding failed coarse-pivot scans while still
reusing each accepted neighbor moment across all 3PCF radial-bin pairs.

The two-ball acceptance requires the full interval
`distance - radius1 - radius2` through `distance + radius1 + radius2` to stay
inside one radial bin and satisfy the `theta`-scaled bin-slop tolerance.
Otherwise the larger node is split with dual-node's `0.585` split rule. The
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
