Benchmarks and Numerical Comparisons
====================================

A timing comparison is meaningful only when geometry, input selection,
weights, bins, estimator, normalization, approximation, and computed orders
match. Validate a reduced catalog before timing a production one.

Driver Benchmarks
-----------------

The maintained convergence, shear, and forest drivers are in
``tests/python``. They read a catalog once, retain it for all selected native
engines, report wall and process CPU time, and save relative differences::

   python3 tests/python/kappa_corr_all_engines.py --list-engines
   python3 tests/python/shear_corr_all_engines.py --list-engines
   python3 tests/python/lya_corr_all_engines.py --list-engines

For a scalar full-sky comparison::

   python3 tests/python/kappa_corr_all_engines.py \
      --fits catalogs/map.fits --engine all-omp --statistics both \
      --threads 16 --outdir results/kappa

Use ``--statistics 2pcf`` to keep pair timings free of compiled 3PCF work.
``--max-points`` performs deterministic input thinning for scaling tests; zero
retains every selected point.

Masks and Edge Correction
-------------------------

All six active scalar angular engines support masks and complex 3PCF edge
correction::

   python3 tests/python/kappa_corr_all_engines.py \
      --fits catalogs/map.fits --mask catalogs/mask.fits \
      --engine octree-2balls-omp,kdtree-2balls-omp \
      --edge-corrections --threads 16 --outdir results/masked

The full input and active count are recorded separately. Masking selects
bodies; edge correction solves the angular window system and requires 3PCF.

General CPU Suite
-----------------

The optional local workspace ``addons/python_env/cputime_comparison`` contains
the broader ``benchmark_kappa_corr.py`` runner. It compares active cTreeBalls
methods with compatible Corrfunc, FCFC, and lya2pcf workloads. This directory
is deliberately excluded from publication branches and source distributions.

From a configured local workspace::

   ./create_benchmark_environment
   conda activate ctreeballs-bench
   python3 benchmark_kappa_corr.py --help

Example scalar run::

   python3 benchmark_kappa_corr.py \
      --scenarios sphere-counts,sphere-convergence \
      --backends ctreeballs,corrfunc \
      --sphere-methods kdtree-2balls-omp,balltree-2balls-omp,octree-2balls-omp \
      --sizes 1000,10000 --threads 1,8 --repeats 5 \
      --outdir results/scalar

Corrfunc provides selected pair-count workloads, FCFC provides periodic
isotropic pair counts, and lya2pcf provides a forest 2PCF reference. Unsupported
combinations are written to ``skipped.csv`` instead of being compared as if
their estimators matched.

MPI Timing
----------

The parent benchmark launches MPI workers; do not launch the parent itself
under ``mpiexec``. Use ``--mpi-ranks`` and repeated ``--mpi-extra-arg`` values.
The C extension, ``mpi4py``, compiler wrapper, and launcher must use the same
MPI implementation. Record ranks and threads per rank because catalog/tree
memory is generally replicated.

Outputs
-------

``summary.json`` and ``timing_report.txt`` are produced by the field drivers.
The general suite writes ``timings.csv``, ``summary.csv``, ``speedups.csv``,
``comparisons.csv``, ``values.csv``, ``skipped.csv``, and ``metadata.json``.
Three-point plots include per-mode matrices and flattened radial-bin views.

Report hardware, affinity, compile flags, smoothing, approximation controls,
warmups, repeats, and timing scope with any performance claim. A faster result
with a different normalization or triplet policy is not a valid speedup.
