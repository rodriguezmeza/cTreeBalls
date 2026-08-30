# cTreeBalls executable tests

This directory contains the shell `run_test*` launchers and their Python
regression programs. Catalogs, parameter files, expected outputs, and C unit
test sources remain one level above in `tests/`.

Run an individual launcher from this directory, for example:

```sh
cd tests/make_tests
./run_test_make_info
./run_test_search_methods
./run_test_openmp_determinism
```

The launchers resolve the repository and `tests/` directories from their own
locations, so they may also be invoked from the repository root:

```sh
tests/make_tests/run_test_make_info
```

The root Makefile remains the preferred way to build the required profile and
run a test together:

```sh
make test-make-info
make test-search-methods
make test-openmp-determinism
```

MPI launchers require an MPI-enabled build and `mpiexec`. Optional-profile
tests may rebuild cTreeBalls with the corresponding Make setting.
