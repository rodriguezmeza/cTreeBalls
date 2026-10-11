# Comparison and CPU benchmark drivers

Run these scripts from the repository root with the matching C/Cython build:

```sh
python tests/python/kappa_corr_all_engines.py --help
python tests/python/shear_corr_all_engines.py --help
python tests/python/lya_corr_all_engines.py --help
```

The corresponding [kappa](README_kappa_corr_all_engines.md),
[shear](README_shear_corr_all_engines.md), and
[Ly-alpha](README_lya_corr_all_engines.md) guides document inputs and controls.
The `benchmark_lya_*.py` programs and their README files provide specialized
accuracy and CPU-time measurements. `benchmark_kappa_corr.py` is a CLI alias for
the kappa comparison driver. The DES shear batch runner lives here too.

The shared `benchmark_timing.py` records wall time, process CPU time, native
MainLoop clocks, and explicit tuning settings. Compare equal catalogs,
estimators, precision, geometry, output controls, and accuracy when comparing
performance. MPI wall time is the maximum rank duration; CPU is summed over
participating ranks.

These drivers and their FITS/catalogue and analysis helpers share this directory.
Preserve the helper files when copying a driver.
The separate portable DES workflows retain their self-contained layout.

Automated tests live in `tests/python/`. Existing `make test-*` commands are
unchanged; their rules are included from `Makefile_tests`.
