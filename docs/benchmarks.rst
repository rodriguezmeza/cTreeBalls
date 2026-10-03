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

Each native driver separates setup and MainLoop wall time from process CPU
time. MPI wall is the maximum across participating ranks, while CPU is summed.
Per-rank measurements and native numerical parameters are retained. The kappa
and forest drivers exclude catalog registration; the shear setup timer includes
it. MainLoop includes any requested native histogram writing; Python extraction,
plots and final cleanup are excluded. Forest ``wall_seconds`` retains its older
cleanup-inclusive scope. Compare the explicitly labeled compute columns only
after checking that requested outputs and statistics also match.

The ``mainloop_wall_s`` / ``mainloop_cpu_s`` columns exclude Python provenance
capture. The older ``compute_*`` columns time the complete Python ``Run`` call,
including that capture. JSON stores the native scope as
``native_mainloop_wall_time`` and ``native_mainloop_cpu_time``. Both scopes use
maximum-rank wall time and summed rank CPU; do not mix them in speed ratios.
See :doc:`capabilities` for the cold-process validation workflow.

The public ``tests/python/benchmark_kappa_corr.py`` is an entry-point alias for
the convergence all-engines driver and uses that driver's command-line options.
It does not depend on the private CPU suite.

For exact unsmoothed scalar references across engines, use
``no-one-ball,no-two-balls,no-smooth-pivot``.
Native octree 3PCF is not necessarily exact with ``no-two-balls`` alone.
For shear pivot reuse, record ``CBALLS_SHEAR_PIVOT_TOL``, ``theta``, leaf
capacity, bins, multipole order and window conditioning. A phase tolerance in
radians is not a 5% coefficient guarantee; validate before reporting speedups.

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

The drivers can launch MPI workers using ``--mpi-ranks`` and repeated
``--mpi-extra-arg`` values. Kappa and forest drivers also recognize an existing
MPI launch; do not nest MPI launchers. The shear driver launches workers itself.
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
