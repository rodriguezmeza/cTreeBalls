# octree-3pcf-3d-mpi

Independent MPI+OpenMP sibling of `octree-3pcf-3d-omp`, with the same
spherical-harmonic scalar 2PCF/3PCF and ENCORE-style survey estimator.
The alias `octree-ggg-3d-mpi` selects the same method (ID 192).

## Build and run

Enable `OCTREE3PCF3DMPION=1` in `addons/Makefile_addons_settings`.
`OCTREE3PCF3DOMPON` independently controls registration of the OpenMP sibling.
Build with an MPI C compiler (`MPICC`, normally `mpicc`), OpenMP, and
`DEFDIMENSION=3`. Rebuild Cython with the same make profile.

From the project root:

```sh
make -j4 OCTREE3PCF3DMPION=1 cballs cyballs-static-lib
mpiexec -n 2 ./cballs addons/octree_3pcf_3d_mpi/parameters.ini
mpiexec -n 2 ./cballs addons/octree_3pcf_3d_mpi/parameters_survey.ini
```

The included small catalogs are correctness examples, not scaling benchmarks.
Replace them with your own catalogs for production use.
On hosts requiring explicit MPI placement, add `--map-by slot --bind-to none`.
`numberThreads` is the number of OpenMP threads **per rank**.

## Estimators and inputs

The scalar path takes `x y z delta weight` and computes an auto-correlation.
ASCII input uses `multi-columns-ascii`, `columns=1,2,3,4,5`, and
`options=pos-and-convergence-weight`. FITS input, per-pixel masks, and
pivot-neighbor LOS exclusions follow the OpenMP sibling's conventions.
The exact physical cutoff remains `rangeN`, including with `Rcut/theta`.

By default only the 3PCF is computed. Select `only-2pcf-3d`,
`only-3pcf-3d`, or `compute-2pcf-3d,compute-3pcf-3d`.
`no-out-Hist` suppresses the correlation files.

For `survey-estimator-3d`, supply two catalogs, data then randoms,
with `iCatalogs=1,2`. Coordinates remain in a common Cartesian frame.
The value column is ignored; non-negative weights define
`alpha=sum(w_data)/sum(w_random)`, `N=D-alpha*R`, and `R=alpha*R`.

Every rank constructs the same owned survey work catalogs and measures its
subset of the N and R pivots. Global raw multipoles are reduced separately
before rank 0 performs the 2PCF ratio and solves the 3PCF window-coupling
matrix. No rank-local normalized correlations are averaged.
Output columns, multipole truncation, and empty/singular-window validity flags
match the OpenMP implementation. Survey mode rejects `read-mask`,
`remove-mean`, and LOS exclusions; the random catalog encodes the selection.

For an empty window, corrected values are zero and `valid=0`; the
`pivot_condition` diagnostic remains infinity because no matrix was solved.

See `../octree_3pcf_3d_omp/README.md` for the detailed estimator definitions.

## Parallel contract

- Each rank reads or receives the same complete catalog(s), and builds its own
  tree. Catalog, tree, histogram, and survey-workspace memory is replicated.
- Fixed blocks of `CB3D_OMP_PIVOT_BLOCK_SIZE` pivots (default 64) are
  assigned cyclically to ranks. OpenMP dynamically schedules local blocks,
  which are committed in block order.
- Thread-count changes are deterministic for a fixed rank count. Rank-count
  changes regroup floating-point sums and can change the last rounding bits.
- Small catalogs with fewer blocks than ranks cannot occupy all ranks.
  There is no measured speedup guarantee; sufficiently large workloads are
  needed to amortize replicated startup and MPI reduction.
- Only rank 0 writes output. Recoverable reader, tree, measurement, and output
  errors are coordinated at the corresponding collective boundaries.
- All ranks must run the same method, parameters, catalog contents, and stage
  sequence. Launch-time/parser failures before MPI initialization and process
  crashes are outside the recoverable collective contract.
- MPI calls occur outside OpenMP regions and require MPI_THREAD_FUNNELED or
  stronger support.

## Cython

Launch Python under `mpiexec` and import `mpi4py.MPI` first. Every
rank must call `Run` and cleanup in the same order. Both file input and
`set_catalog` work; for a survey, register data as catalog 0 and randoms
as catalog 1 on every rank. This is replicated input, not a rank-0-only
catalog broadcast. Read output files on rank 0 after `Run` returns.
MPI initialized by Python is not finalized by cyballs; reuse after
`MPI.Finalize()` is rejected.

The angular-kappa driver deliberately excludes this physical 3D statistic.

## Regressions

```sh
make test-octree-3pcf-3d-omp
make MPIEXEC='mpiexec --oversubscribe' test-octree-3pcf-3d-mpi
PYTHONPATH="$PWD:$PWD/python" python3 tests/make_tests/test_octree_3pcf_3d_mpi.py \
  --cballs ./cballs --mpi-command 'mpiexec -n 2' --cython
```

Cython MPI tests require `mpi4py`. Add `--mpi-only` (or
`CB3D_MPI_TEST_ARGS=--mpi-only` with make) when the OpenMP sibling is disabled.
Add `--fits` to check FITS 2PCF/3PCF and LOS exclusion when CFITSIO and Astropy
are available.
Tests include independent scalar and survey oracles, modes, multi-block
rank/thread comparisons, periodic and logarithmic bins, Rcut/theta, in-memory
masks, collective failures, repeated Cython recovery, and post-finalization
rejection. CI builds both with and without the OpenMP sibling.
