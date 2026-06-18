Performance and Parallelization
===============================

cTreeBalls performance is controlled by the search method, tree acceptance
settings, histogram resolution, catalog size, and OpenMP configuration.

Search Methods
--------------

The repository contains octree, k-d-tree, neighbor-box, and historical
tree/ball implementations.  Available methods depend on
``addons/Makefile_addons_settings``.  ``octree-ggg-omp`` is the principal
three-point path in current examples; box-oriented methods are used for
selected Cartesian 2PCF workflows.

Do not compare method timing without first verifying that geometry, bins,
normalization, smoothing, and enabled statistics are equivalent.

Search Controls
---------------

``theta``
    Controls tree-cell acceptance and therefore the speed/accuracy tradeoff.
    Validate any value changed from the repository default.

``nsmooth`` and ``rsmooth``
    Control bucket or pivot-neighborhood smoothing in supported methods.

``options=no-one-ball``
    Disables the one-ball acceleration path in supported searches.  It is
    useful for validation but is normally slower.

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
* compare accelerated results with a direct or ``no-one-ball`` validation run
  on a reduced catalog;
* repeat timings and record thread affinity and hardware;
* retain the used-values file and Makefile settings.

Benchmark scripts and notes are available under ``tests/python`` and
``tests/Readme_benchmarks.txt``.
