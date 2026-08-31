
octree-balls4-omp

1. Include searching method: octree-balls4-omp

This module works only in log scale binning: set useLogHist=true

Enable OCTREEBALLS4OMPON=1. This is the renamed octree-kkk-balls4-omp
addon (method ID 69); use the new name in parameter files and Python scripts.
Its normal-tree B4 partition does not require OCTREESMOOTHINGON, BALLSON,
or BALLS4SCANLEVON. The separate triangle addon has not been renamed.

TWOPCFON=1 enables 2PCF and TPCFON=1 enables 3PCF.
Global SMOOTHPIVOTON may be 0 or 1 for other engines. This addon has no
smooth-pivot capability and rejects options=smooth-pivot in C and Cython.
The no-one-ball and no-two-balls options still control cell aggregation.

Run make SMOOTHPIVOTON=0 test-octree-balls4-no-smoothing and repeat with
SMOOTHPIVOTON=1 to exercise both profiles.

read-mask accepts an embedded mask or a cyballs in-memory mask. Add
edge-corrections,no-normalize-HistZeta for complex corrected multipoles;
weights-norm enables weighted signal/window moments. Edge runs use body
pivots in the B4 partition and exclude repeated neighbors. Uncorrected
no-normalize-HistZeta runs use this same distinct body-pivot kernel.
Only legacy normalized non-edge runs retain the cell-pivot estimator.
These are not the legacy cell-pivot 3PCF followed by a correction. Window orders reach
2*mChebyshev. Either no-one-ball or no-two-balls makes edge moments exact.
The edge walker builds geometric enclosing radii independently of theta;
theta controls only its angular approximation. Cells with excessive projected
bearing extent or crossing a radial bin boundary are opened.
Three-dimensional sky coordinates preserve the observer origin and use chord
bins. See docs/3pcf.rst for raw weights and complex-mode conventions.
Run make test-octree-balls4-edge for the numerical oracle and driver tests.

OCTREEBALLS4MPION=1 adds octree-balls4-mpi with the same options. See
addons/octree_balls4_mpi/README.md and docs/addons.rst for MPI details.
