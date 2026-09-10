# balltree-shear-sphere-2balls-omp

This addon computes full-sky weak-lensing shear `xi+`, `xi-`, and all four
natural Gamma-x 3PCF multipole components over an FCFC-style PCA ball tree.
Input positions are observer-centered 3D vectors and shear components are
defined in each sample's local east/north frame.

Nodes are split at the median of their dominant principal axis and carry a
conservative spherical chord-radius bound. Their weighted first and second
spin-2 moments are parallel transported into the tangent frame at the node
center. The estimator shares the dual-node-style dual-node 2PCF and accepted
neighbor-node LogMultipole 3PCF traversal with the KD-tree variant, including
the combined-radius opening criterion and 0.585 split heuristic.

`only-2pcf` and `only-3pcf` skip unused work. `no-two-balls` or `no-one-ball`
forces body-level results; `dual-node-bin-slop` enables the looser dual-node
radial criterion, and `nsmooth` sets leaf capacity. Masks, shear
mode-coupling edge correction, deterministic `BALLS4SCANLEVON`, and the
`SMOOTHPIVOTON`/`no-smooth-pivot` contract match the spherical octree and KD
engines. Explicit smoothing radii must satisfy `2*rsmooth <= rminHist`.

The PCA split and initial enclosing-sphere strategy are adapted from FCFC by
Cheng Zhao under the MIT license. The two-node traversal strategy is adapted
from dual-node by Mike Jarvis under its BSD-style license; see
`dual-node_LICENSE`.
