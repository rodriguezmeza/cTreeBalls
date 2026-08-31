# Lyman-alpha forest MPI addons

Enable `LYAFORESTMPION=1` in `addons/Makefile_addons_settings`. The
OpenMP-only siblings (`LYAFORESTOMPON`) may independently be enabled or disabled.
Build with an MPI C compiler (`MPICC`, normally `mpicc`), OpenMP, and
`DEFDIMENSION=3`. Catalog storage stays 3D even for radial searches.

| MPI method | ID | Estimator |
| --- | --- | --- |
| lya-2pcf-mpi | 185 | 3D 2PCF in parallel/transverse bins |
| lya-3pcf-mpi | 186 | 3D anisotropic five-dimensional 3PCF |
| lya-2pcf-3pcf-mpi | 187 | Both 3D estimators |
| lya-1d-2pcf-mpi | 188 | Radial-only 2PCF |
| lya-1d-3pcf-mpi | 189 | Radial-only signed-lag 3PCF |
| lya-1d-2pcf-3pcf-mpi | 190 | Both radial estimators |
| lya-1d-tree-2pcf-mpi | 191 | Exact interval-tree radial 2PCF |

## Run

From the cTreeBalls root, the included eight-pixel example is runnable with:

```sh
make -j4 LYAFORESTMPION=1 cballs
mpiexec -n 2 ./cballs addons/lya_forest_mpi/parameters.ini
```

Replace `infile` with a six-column `x y z delta weight forest_id` catalog.
All ranks must use the same catalog, parameters, and method. `numberThreads`
is the number of OpenMP threads **per rank**. Do not request more total
threads than the machine/job allocation supplies.

The three-dimensional modes use `lya2RpMax`, `lya2RtMax`, `lya2RpBins`,
`lya2RtBins` and the `lya3RMax`, `lya3RBins`, `lya3ThetaBins`, `lya3MuBins`
parameters. Radial modes ignore transverse separation; their 2PCF uses
`lya2RpMax`/`lya2RpBins`, and their 3PCF uses two signed radial lags with
`2*lya3RBins` bins per axis. Require one `lya-ascii` catalog and
`usePeriodic=false`. These estimators do not approximate pivots by smoothing.

## Parallel and numerical contract

Each rank reads the catalog and builds its own tree/index. Memory for the
catalog and multidimensional histograms is replicated, not distributed.
The shared OpenMP kernels evaluate disjoint cyclic subsets of work:

- 3D modes partition body pivots.
- Radial scans partition fixed `LYA1D_OMP_PIVOT_BLOCK_SIZE` pivot blocks.
- The radial tree partitions fixed blocks of node-pair tasks and, separately,
  same-forest subtraction tasks.

There must be enough logical blocks to occupy the requested ranks. Small
catalogs may use only one rank for a radial/tree workload.

Numerators, denominators, and exact integer counts are reduced to rank 0.
The radial tree preserves long-double accumulation and reduces all-pair and
same-forest sums **before subtraction**. Normalization is performed only after
these reductions. Pixels in the same forest are excluded from pairs, and
triplets require three different forests. Zero-denominator bins return zero.

Ordered local accumulation is repeatable across thread counts for a fixed
rank count. Changing the number of ranks changes floating-point grouping:
results agree to rounding precision, not necessarily byte for byte.

The 3D octree uses a conservative geometric cell bound, independent of the
theta-scaled opening radius. MPI calls are outside OpenMP regions and require
at least `MPI_THREAD_FUNNELED`. Rank-wide status checks cover catalog reads,
tree/search allocation and evaluation, and publication errors. Only rank 0
writes histograms and logs. Output filenames and columns match the corresponding
`-omp` method; `no-out-Hist` and `lya-output-empty-bins` retain their meaning.

## Cython

Cython accepts the same seven method names and file-based input. Launch Python
under `mpiexec`, import `mpi4py.MPI`, and have **every rank** call `Run` and
cleanup in the same order. Read the output files only on rank 0 after `Run`
returns. MPI initialized by Python is not finalized by a cyballs object.
Generic kappa in-memory arrays do not include forest IDs or the required
radial metadata and are not a replacement for `lya-ascii` input.
Use the forest-specific `set_forest_catalog(positions, delta, weights, forest_ids)`
API instead. The `python/lya_corr_all_engines.py` driver reads DESI FITS, ASCII
or NPZ once on rank 0, broadcasts NumPy arrays once, and retains them across all
selected engines. See `python/README_lya_corr_all_engines.md`.

## Tests

```sh
make test-lya-forest-omp test-lya-forest-1d-omp
make MPIEXEC='mpiexec --oversubscribe' test-lya-forest-mpi
PYTHONPATH="$PWD:$PWD/python" python3 tests/make_tests/test_lya_forest_mpi.py \
  --cballs ./cballs --mpi-command 'mpiexec -n 2' --cython
```

The Cython MPI check additionally requires `mpi4py` in that Python environment.
Use `--mpi-only` (or `LYA_MPI_TEST_ARGS=--mpi-only` with make) when the OpenMP-only
siblings were intentionally disabled. Tests cover every method, independent
2PCF/3PCF oracles, rank/thread comparisons, radial scan versus interval-tree
results, zero windows, same-quasar exclusion, and recoverable input/output errors.
