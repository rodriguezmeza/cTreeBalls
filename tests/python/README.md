# Python programs for cTreeBalls

This directory contains the maintained analysis programs, CPU-time benchmarks,
regression tests, catalogue readers and plotting utilities, together with their
README files. Run them from the repository root with a matching native and Cython
build. The `python/` directory contains only the Cython binding sources.

- [CPU-time benchmarks](README_benchmark_suite.md): Lyman-alpha, kappa and shear.
- [Regression tests](README_regression_tests.md): numerical and interface checks.
- [Cluster jobs](../slurm_scripts/README.md): SLURM launchers remain beside their
  shell configuration; Python job tools are in this directory.

The three entry points are `lya_corr_all_engines.py`,
`kappa_corr_all_engines.py` and `shear_corr_all_engines.py`.
Use `--help` for input units, correlation orders and timing controls. CPU time is
consumed process time, including OpenMP workers; MPI totals sum participating
ranks. Elapsed time is the maximum participating-rank duration. Compare the same
catalogue, estimator, geometry and numerical accuracy before comparing timings.

The optional external dual node adapter is a reference calculation, not another
native search method. Its package is credited in the main README Acknowledgments.
The executable reports the active methods with `options=print-search-methods`.

Native C tests, fixtures and shell launchers keep their locations under `tests/`.
Historical plotting examples may require their original input catalogues;
`example_scalar_catalog.py` resolves its catalogue relative to the repository.
