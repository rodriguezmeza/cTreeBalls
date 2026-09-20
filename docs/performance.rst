Performance and Parallelization
===============================

cTreeBalls performance is controlled by the search method, tree acceptance
settings, histogram resolution, catalog size, and OpenMP configuration.

Search Methods
--------------

The maintained profile provides octree, median-KD-tree, and PCA-ball-tree
dual-node angular engines, full-sky shear variants, forest searches, physical
3D multipoles, and two periodic box methods. Availability is controlled by
``addons/Makefile_addons_settings`` and reported by
``options=print-search-methods``.

Do not compare method timing without first verifying that geometry, bins,
normalization, smoothing, and enabled statistics are equivalent.

Search Controls
---------------

``theta``
    Controls tree-cell acceptance and therefore the speed/accuracy tradeoff.
    Validate any value changed from the repository default.

``nsmooth``
    Usually controls leaf capacity in KD/ball trees; meaning is engine-specific.

``rsmooth``
    Smoothing radius for supported methods. ``SMOOTHPIVOTON=1`` enables the
    smoothing path by default; ``no-smooth-pivot`` disables it for one run.

``options=no-one-ball``
    Disables the one-ball acceleration path in supported searches.  It is
    useful for validation but is normally slower.

``options=no-two-balls``
    Exact body-pair limit for 2PCF, with smoothing separately disabled.
    For exact scalar 3PCF across the active methods use
    ``no-one-ball,no-two-balls,no-smooth-pivot``;
    native octree 3PCF still permits neighbor-cell acceptance with
    ``no-two-balls`` alone. ``BALLS4SCANLEVON=1`` may remain enabled.

``options=dual-node-bin-slop``
    Uses dual-node-compatible Log/Linear near-bin acceptance in two-ball 2PCF
    scans. The conservative default requires both ball-distance bounds to
    remain in the center bin. Calibrate each method to the same numerical
    error target before comparing performance; matching ``theta`` alone does
    not establish equivalent accuracy.

``octree-2balls`` tree preparation
    These engines stop native-octree preparation after production cell
    aggregation. They do not construct the threaded walk, scan-level arrays,
    or pruning products used by legacy octree searches. Their compact binary
    view owns its point order and exact centroid radii independently.

In 3D, scalar angular acceptance bounds projected bearings rather than only
chord lengths. This can open more cells near degenerate bearings. A faster
older calculation using a different angular estimator is not an equivalent
performance reference. See :doc:`3pcf`.

Histogram Cost
--------------

``sizeHistN`` determines the radial grid size.  3PCF multipole matrices scale
with two radial dimensions, so memory and output grow rapidly with this value.
``mChebyshev`` controls the number of 3PCF multipoles; use convergence tests
before increasing both settings together.


OpenMP
------

OpenMP is enabled at build time with ``OPENMPMACHINE = 1`` and
``OMPFLAG = -fopenmp``.  Set threads at runtime with ``numberThreads`` and,
where appropriate, the environment:

.. code-block:: bash

   export OMP_NUM_THREADS=8
   ./cballs numberThreads=8 rootDir=Output_threads

Benchmark physical cores first.  More threads can increase memory traffic and
may not improve small-catalog runs.

Convergence and Benchmarking
----------------------------

For a production configuration:

* compare at least two ``theta`` values;
* vary ``sizeHistN`` and ``mChebyshev`` independently;
* compare accelerated results with a direct or
  ``no-one-ball,no-two-balls,no-smooth-pivot`` scalar validation run
  on a reduced catalog;
* repeat timings and record thread affinity and hardware;
* retain the used-values file and Makefile settings.

See :doc:`benchmarks` for estimator contracts, MPI timing scopes, and the
maintained scripts in ``tests/python``. Compare saved arrays as well as timing
tables when checking a new build. Synthetic measurements are not predictions
for a full-resolution survey catalog.
