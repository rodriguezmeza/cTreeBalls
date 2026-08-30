# Octree two-ball OpenMP addon

`search=octree-2balls-omp` applies TreeCorr-style scans to the native
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
Otherwise the larger node is split with TreeCorr's `0.585` split rule. The
3PCF also bounds angular phase error by `theta*pi/(2*mChebyshev+1)`. Use
`options=no-two-balls` for exact body-pair and body-moment accumulation.

`TWOPCFON=1` and `TPCFON=1` enable the two correlation orders at build time.
When both are active, `only-2pcf` and `only-3pcf` select one at runtime.
`compute-HistN`, `weights-norm`, `no-normalize-HistZeta`, and
`out-m-HistZeta` have the same meanings as for `balltree-2balls-omp_3pcf`.

The fixed pivot frontier owns private moment and histogram buffers and is
reduced in a fixed order, so OpenMP worker count does not alter the result.
`options=treecorr-direct-triples` retains the cubic triple-node traversal as a
validation oracle for moderate catalogs. The traversal is adapted from
TreeCorr by Mike Jarvis under its BSD license; the full notice is in
`addons/balltree_2balls_omp/TreeCorr_LICENSE`.
