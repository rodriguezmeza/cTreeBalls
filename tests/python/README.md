# Python tests, analyses and benchmarks

Run from the repository root with a matching freshly built C/Cython profile.
The Cython binding source stays in `python/`; all Python tests and analysis
scripts are here. C test programs and shell/native launchers remain under
`tests/` and `tests/make_tests/`.

CPU benchmark entry points:

- [Scalar/convergence](README_kappa_corr_all_engines.md): `kappa_corr_all_engines.py`.
- [Shear](README_shear_corr_all_engines.md): `shear_corr_all_engines.py`.
- [Forest](README_lya_corr_all_engines.md): `lya_corr_all_engines.py`.

All three share `benchmark_timing.py`. Compare equal catalogs, estimator,
precision, geometry, output controls and accuracy before comparing times.
Complete Python Run, native MainLoop and launcher scopes are recorded separately.
MPI wall time is the maximum rank time; CPU is summed over participating ranks.

Run focused regression checks:

```sh
PYTHONPATH=. python -m pytest -q tests/python/test_benchmark_timing.py \
  tests/python/test_public_profile.py tests/python/test_capability_contracts.py
make test-make-info test-search-methods
```

[Regression guide](README_regression_tests.md) lists the wider test suite.
Plotting and catalog conversion tools retain their required input data paths;
run historical data-specific examples from `tests/` unless their CLI states
otherwise. `example_scalar_catalog.py` resolves its example input relative to
this repository and performs no work on import.
