# kdtree-shear-sphere-2balls-omp

This addon computes full-sky weak-lensing shear `xi+`, `xi-`, and all four
natural Gamma-x 3PCF multipole components over an independent median KD tree.
Input positions are observer-centered 3D vectors and shear components are
defined in each sample's local east/north frame.

Every KD node stores weighted first and second spin-2 moments transported into
the tangent frame at its spherical center. The 2PCF uses a dual-node-style
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
`no-one-ball` forces body-level results; `dual-node-bin-slop` enables the looser
dual-node radial criterion. Masks and the shared shear mode-coupling edge solve
are supported. With `SMOOTHPIVOTON=1`, deterministic transport-aware smoothing
is on by default and `no-smooth-pivot` disables it. The safety requirement
`2*rsmooth <= rminHist` is inherited from the spherical shear estimator.

The two-node split strategy is adapted from dual-node by Mike Jarvis under its
BSD-style license; see `dual-node_LICENSE`.
