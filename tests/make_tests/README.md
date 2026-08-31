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
make test-scalar-numerical-contract
make test-lya-forest-mpi
make test-lya-corr-all-engines
make test-octree-3pcf-3d-mpi
```

MPI launchers require an MPI-enabled build and `mpiexec`. Optional-profile
tests may rebuild cTreeBalls with the corresponding Make setting.

`test_scalar_numerical_contract.py` checks raw weighted scalar multipoles
against an independent oracle, rotations, inactive smoothing, degenerate
bearings, and C file input. With matching C/Cython binaries, run it through
pytest to include the optional TreeCorr reference. Set `CBALLS_TEST_MPI=1`
and launch it with `mpiexec -n 2` to check the MPI sibling engines. CI repeats
these checks with `SMOOTHPIVOTON=0` and `SMOOTHPIVOTON=1`.
