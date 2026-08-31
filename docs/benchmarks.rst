Benchmarks and Numerical Comparisons
====================================

A timing comparison is meaningful only when geometry, input selection,
weights, bins, estimator, normalization, approximation, and computed orders
match. The :doc:`search_methods` families deliberately measure different
quantities. See :doc:`performance` for tuning after validating correctness.

Correctness Before Timing
-------------------------

For the scalar angular family, use
``no-normalize-HistZeta,weights-norm`` and the distinct-triplet contract in
:doc:`3pcf`. All listed raw kernels remove repeated neighbors internally,
including KD-tree, legacy balltree, and BALLS4. Do not subtract them again in
Python. Reconstruct complex modes before comparing; separate sine/cosine
component matrices depend on the pivot tangent basis.

On a reduced catalog, compare ``no-one-ball,no-two-balls`` with decreasing
``theta`` and record errors per multipole/bin. Exact numerical tests cover
rotations, nonuniform weights, partial skies, undefined bearings, file versus
memory input, and OpenMP/MPI execution::

   python3 -m pytest -q tests/make_tests/test_scalar_numerical_contract.py
   CBALLS_TEST_MPI=1 mpiexec -n 2 python3 tests/make_tests/test_scalar_numerical_contract.py

Use matching C/Cython builds, with both ``SMOOTHPIVOTON=0`` and
``SMOOTHPIVOTON=1`` tested without runtime smoothing. TreeCorr is an optional
reference dependency for the pytest run. CI repeats these contracts; passing
small catalogs does not certify full-resolution performance.

General CPU Suite
-----------------

The optional local workspace ``addons/python_env/cputime_comparison`` contains
``benchmark_kappa_corr.py``, focused ENCORE scripts, and environment tools.
It may be absent from a published checkout or sdist; it is not a runtime
dependency. This directory is separate from the older script of the same
name under ``tests/python``.

When that workspace and the external source snapshots are present, run from
its directory::

   ./create_benchmark_environment
   conda activate ctreeballs-bench
   ./verify_benchmark_environment.py
   python3 benchmark_kappa_corr.py --help

The environment bootstrap installs/builds local TreeCorr, Corrfunc, FCFC,
lya2pcf, ENCORE, and cTreeBalls sources, and registers the
``Python (cTreeBalls benchmarks)`` Jupyter kernel. It rebuilds native code:
review its README and source-path overrides before running it.
ENCORE and FCFC are native executables, not Python packages.

A bounded scalar comparison, from that directory::

   python3 benchmark_kappa_corr.py \
       --scenarios sphere-counts,sphere-3pcf-convergence \
       --backends ctreeballs,treecorr \
       --sphere-methods octree-ggg-omp,kdtree-omp,balltree-omp \
       --sphere-3pcf-methods octree-ggg-omp,kdtree-omp,balltree-omp,octree-balls4-omp \
       --sizes 32,96 --threads 1,3 --nbins 4 \
       --theta-min 1 --theta-max 120 --warmups 0 --repeats 1 \
       --no-plots --fail-on-mismatch --outdir results/scalar-check

These benchmark angles are degrees and are converted to chord distances.
The shared contract is ``scalar-3pcf-logmultipole-distinct``. The historical
``--keep-repeated-neighbors`` option does not restore repeated neighbors in
the current distinct-only kernels.

Use ``--help`` for scenarios/method controls and the native executable's
``options=print-search-methods`` for compiled methods.
The parent benchmark launches selected MPI workers; do not launch the parent
under ``mpiexec``. Use ``--mpi-ranks`` and the same MPI runtime as cyballs.

Outputs and Interpretation
--------------------------

* ``metadata.json`` records command, software paths, build information, and timing policy.
* ``values.csv`` records bins and real/imaginary multipoles.
* ``comparisons.csv`` records zero-safe symmetric relative differences.
* ``speedups.csv`` groups compatible ``work_contract`` values.

Time comparable work: a combined 2PCF/3PCF engine is not a pure-2PCF baseline.
MPI timings can include launch and collective costs; inspect the recorded
timing scope. Repeat measurements and record hardware, affinity, ranks,
threads, warmups, and approximation settings. Failed numerical comparisons
must not be promoted to speedup claims.

Other References
----------------

TreeCorr supports matched scalar/shear reference workflows. Corrfunc and FCFC
serve selected 2PCF/counting scenarios, not every 3PCF or field.
lya2pcf is a CPU 2PCF reference; it is not an independent 3PCF reference.
ENCORE has separate periodic scalar and nonperiodic data/random survey
benchmarks. Use the matching estimator and executable.

For catalog workflows and plots outside the optional timing workspace, use
the kappa and forest all-engines drivers and the examples under ``examples/``.
See :doc:`user/python`, :doc:`lyman_alpha`, and :doc:`scalar_3d`.
