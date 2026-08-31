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

``nsmooth``
    Usually controls leaf capacity in KD/ball trees; meaning is engine-specific.

``rsmooth``
    Smoothing radius for supported methods, active only with ``smooth-pivot``.
    Compiling ``SMOOTHPIVOTON=1`` or setting a radius alone does not activate it.

``options=no-one-ball``
    Disables the one-ball acceleration path in supported searches.  It is
    useful for validation but is normally slower.

``options=no-two-balls``
    Exact body-level reference for the two-ball engines.
    Raw KD/legacy balltree multipoles and all BALLS4 runs reject ``smooth-pivot``.

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

For ``octree-ggg-omp`` and ``octree-ggg-mpi``, ordinary 3PCF runs compute
signal and window modes from zero through ``mChebyshev``. Edge correction
automatically extends the window through ``2*mChebyshev``. Add
``options=ggg-full-window`` to retain that extended window without applying
edge correction, for example when saving files for later correction.
``out-m-HistZeta`` alone does not request the extended window. With
``mChebyshev=7``, ordinary runs write eight window modes instead of fifteen;
signal multipoles, normalization, and self-pair subtraction are unchanged.
Use a fresh output directory when changing the window policy so files from
an earlier extended-window run are not mistaken for current results.

The signal and window share one angular recurrence. Their four 3PCF
components are accumulated directly without per-pivot temporary matrices.
Both masked-cell acceptance and deterministic chunk ordering are preserved.

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

Add ``ggg-profile`` to the options list and enable ``verbose=1`` or
``verbose_log=1`` to profile either GGG engine. The reported ``work``,
``wait``, and ``merge`` are sums of worker wall-clock intervals, not process
CPU times. ``wait`` measures waiting to enter the ordered chunk reduction;
``wall`` measures elapsed search time through completion of the OpenMP workers.
For MPI these are per-rank measurements, before the MPI histogram reduction.
The default path makes no profiling clock calls. Use enough pivots for
multiple 4096-pivot chunks per worker when measuring parallel scaling.

Convergence and Benchmarking
----------------------------

For a production configuration:

* compare at least two ``theta`` values;
* vary ``sizeHistN`` and ``mChebyshev`` independently;
* compare accelerated results with a direct or ``no-one-ball`` validation run
  on a reduced catalog;
* repeat timings and record thread affinity and hardware;
* retain the used-values file and Makefile settings.

See :doc:`benchmarks` for the current multi-code suite, estimator contracts,
MPI timing scopes, and environment. Historical scripts remain in ``tests/python``.

For the default 3D GGG profile, run the focused regressions and a bounded
synthetic masked-sky benchmark with:

.. code-block:: bash

   make test-octree-ggg-fast-path
   python tests/make_tests/test_octree_ggg_fast_path.py \
       --benchmark ggg-timing.json --nbody 1048576 --threads 16 --repeat 3

The benchmark writes timing JSON and an adjacent NPZ containing the retained
histograms. Add ``--full-window`` to measure the extended-window path or
``--no-one-ball`` to measure direct neighbor accumulation. Compare NPZ values
as well as timings when checking a new build. These synthetic measurements
are not predictions for a full-resolution HEALPix catalog.
