# Active-profile release gate

Run this gate from the source tree after changing source, settings, dependencies,
or the compiler. It rebuilds the native executable/static library and Cython
extension together, resolves the actual Make settings and native registry, then
runs the active numerical and interface checks. It never enables an addon to
make a test pass. A missing MPI launcher, failed check, unknown registered engine,
missing result, or changed build identity makes the gate fail.

The maintained matrix covers the current 3D, double-precision profile: core
sine/cosine octree; scalar KD/ball/octree two-ball methods; three full-sky shear
methods; KD-box and neighbor boxes; physical 3D correlations; and all forest
methods. There are 22 OpenMP and 12 MPI search names. Other dimensions, precision
profiles, or newly registered engines need an explicit matrix extension.

## Invocation

Use the Python interpreter that will load `cyballs`. It needs NumPy, Cython,
setuptools, pytest, SciPy, Astropy, healpy, matplotlib, and mpi4py. Install
`requirements/release.txt` in a fresh virtual environment; native prerequisites
are a working OpenMP/MPI compiler, external GSL, CFITSIO, zlib and pkg-config. The MPI used by mpi4py must
match the compiler wrapper and runtime used to build cTreeBalls. The host must
permit local MPI communication. Choose a **new** output directory for every run.

```sh
python3 scripts/active_release_gate.py \
  --output /absolute/path/to/new-release-results \
  --mpi-command 'mpiexec --host localhost:4 --oversubscribe --bind-to none'
```

The equivalent Make entry point is:

```sh
make active-release-gate PYTHON=python3 \
  RELEASE_GATE_ARGS="--output /absolute/path/to/new-release-results --mpi-command 'mpiexec --host localhost:4 --oversubscribe --bind-to none'"
```

The launcher string excludes `-n`; the gate supplies one and two ranks. OpenMP
runs use one and two threads; existing family scripts also exercise additional
thread counts. `--jobs` controls compilation parallelism only. There is no
skip-MPI release mode. A diagnostic environment without MPI cannot earn PASS.
This runner resolves the saved settings files and environment; do not use parent
Make command-line feature overrides as a substitute for a saved release profile.

On the validated macOS host, mpi4py was installed into an isolated dependency
directory. The MPI wheel needed `DYLD_LIBRARY_PATH` pointing to the matching
Open MPI library directory, with `-x DYLD_LIBRARY_PATH` in the launcher string.
Provide the dependency directory through an absolute `PYTHONPATH`. No global
Python environment needs to be changed by the gate. If the system Make/shell
launch drops `DYLD_LIBRARY_PATH`, set it on the Python invocation itself:

```sh
make active-release-gate \
  PYTHON='env DYLD_LIBRARY_PATH=/path/to/openmpi/lib /path/to/python3' \
  RELEASE_GATE_ARGS="--output /absolute/path/to/new-release-results --mpi-command 'mpiexec --host localhost:4 --oversubscribe --bind-to none -x DYLD_LIBRARY_PATH'"
```


## Gate evidence

`gate.json` records overall status, exact commands, return codes, timings, logs,
resolved settings, canonical registry, binary hashes, and numerical comparisons.
`build-fingerprint.json` retains the individual source SHA-256 hashes. A PASS
requires the complete matrix; individual inactive-engine pytest skips do not
substitute for the explicit registered-engine coverage check.

There are two retained runs per registered method, in `cases/`. Each contains
`result.npz`, `result.json`, `fixture.npz`, native numerical files, and
`run-metadata.json` when native output is enabled. Every OpenMP method is checked
across thread counts. Every MPI method is compared with the matching OpenMP
method and with another MPI rank/thread configuration. Independent pair,
Fourier/window, spin-2, Legendre/survey, and forest enumerations remain in use.
The core sine/cosine retained case checks finite products and repeatability; its
existing native family comparison scripts provide additional numerical coverage.

Cross-execution comparisons use relative and absolute tolerances of `3e-12`,
with identical finite/unsupported-bin masks. Each independent oracle keeps its
own pre-existing tolerance. The MPI 3D suite also retains byte-exact comparisons
for thread changes and aliases, and larger multiple-block cases ensure that
non-root ranks receive real work. Rank-count changes may regroup floating-point
sums; the gate does not promise bitwise equality across rank counts.

The three `lya-los-tree-*` methods each retain independent 3D oracle checks,
thread comparisons, and comparisons against their original `lya-*` counterparts.
The focused LOS-tree suite also covers larger multi-block forests, adversarial
geometry, unequal 2PCF/3PCF domains, wide forest IDs and rejected parameters.
The in-memory scalar suite uses the public octree two-ball engine. Unavailable
legacy compatibility comparisons appear as explicit pytest skips with reasons.
The startup suite also checks successful large-tree runs and cleanup after the
cellRadius overflow-bin repair; it no longer expects that repaired input to fail.
It remains part of the active gate. Its memory check uses live malloc bytes on
macOS (where freed allocator pages can remain resident) and RSS elsewhere.
The gate records the two-rank mpi4py vendor/library probe and Python package list.

The gate invokes the actual CLI entry points of mask, scalar edge, shear,
physical 3D, forest, I/O, and runtime scripts. It also invokes all three maintained
analysis drivers (`kappa_corr_all_engines.py`, `shear_corr_all_engines.py`, and
`lya_corr_all_engines.py`) on small generated catalogs, retaining their products.
These driver demonstrations run active OpenMP methods; MPI is exercised by the
retained per-engine workers and the native/Cython MPI regression entry points.

`io/routes.json` and its catalog directories retain successful format checks.
These include native ASCII/binary forms, Gadget, IOLIB coordinate readers,
FITS XYZ/RA-Dec/RA-Dec-radius, HEALPix, raw-double maps, and mask branches.
Earlier recorded Gadget/FITS and runtime failure fixtures run separately.

Some legacy names describe asymmetric formats: `columns-ascii-pos` writes a
headerless position table but its reader expects a native header and assigns a
constant scalar of 2. The `numpy-healpix` writer emits an XYZ/KAPPA FITS table;
its ordinary reader consumes headerless native-endian doubles, while its mask
branches use FITS. The route tests check these actual contracts explicitly and
do not silently fall back to an unrecognized format name.

## Build identity and result provenance

`make print-build-fingerprint` prints the resolved identity. The native query is
`cballs options=build-fingerprint verbose=0 verbose_log=0`; Python exposes
`cyballs.build_info()`. The identity includes resolved switches, compiler/linker
flags, source digest, compiler and MPI wrapper information, precision profile,
build-package versions, GSL/CFITSIO versions when their configuration tools are
available, architecture, and byte order. Source hashing is conservative and also
includes disabled implementations and regression inputs; packaged trees with
different archived source content can therefore have different IDs. External
library versions/paths are recorded, rather than recursively hashing an entire
system installation. Native and extension artifact SHA-256 hashes are retained
in the gate manifest separately from the source/profile identity.

The same identity is compiled into the native library and Cython translation
unit. Python checks their agreement, then compares structure sizes returned by
the **native library** against its declarations. A stale native/Cython pair is
rejected even when its structure sizes happen to match.

Successful native runs with file output write `run-metadata.json` alongside the
products. A metadata open/write/close failure fails the run and participates in
MPI failure consensus for correlation output. For Python in-memory results:

```python
balls.Run()
metadata = balls.getRunMetadata()  # independent, JSON-serializable copy
immutable = balls.run_settings['provenance']
```

The snapshot remains readable after `struct_cleanup()`. Changes to parameters or
catalogs invalidate it. The retained metadata contains the build ID and settings,
engine/name ID, estimator, bin edges (including forest axes), coordinate and
weight/mask policies, input descriptors, effective smoothing, opening controls,
precision, and rank/thread counts. In-memory catalogs additionally carry SHA-256
hashes, shapes, and dtypes for registered positions, fields, weights, masks, and
forest identifiers. Retained gate fixtures allow those inputs to be reconstructed.

`openmp_max_threads` is the runtime limit; `openmp_probe_threads` is the size of an
actual post-run OpenMP team. Neither is a claim that every worker received useful
work in every phase. `world_ranks` records the MPI communicator size and
`estimator_ranks` distinguishes distributed estimators from an OpenMP method
executed inside an MPI process. Exact per-phase utilization is outside this gate.

Always check the process result. Reusing a native output directory outside the
gate can leave products from an earlier run after a later failure. The gate
refuses existing result directories and never treats those stale files as a pass.

A passing gate validates the retained fixtures and active runtime contracts. It
is not a production-scale accuracy, performance, memory-scaling, or unrestricted
Python concurrency certification.

## Source distribution and clean installation

The public distribution is `cyballs` (project: cTreeBalls). It ships the active
kernels and build scripts, and requires external GSL/CFITSIO rather than vendored
library source trees. Both CI and the publishing workflow run the full active
gate before accepting the release, retain its evidence even on failure, check
the source archive, and exercise a fresh installation outside the checkout:

```sh
python -m build --sdist
python -m twine check --strict dist/*.tar.gz
python scripts/check_sdist.py dist/*.tar.gz
python -m venv /absolute/path/to/install-env
/absolute/path/to/install-env/bin/python -m pip install /absolute/path/to/dist/cyballs-VERSION.tar.gz
cd /absolute/path/outside/the/checkout
/absolute/path/to/install-env/bin/python /absolute/path/to/checkout/scripts/check_installed_package.py \
  --source-root /absolute/path/to/checkout --output /absolute/path/to/install-result.json
```

The installation check rejects a module loaded from the checkout, verifies the
installed distribution and public registry, and computes known weighted pair
and triplet products with all three LOS-tree modes. A successful import alone
does not count as an installation test. The source archive has its own build
fingerprint because it deliberately excludes inactive development source trees.

## Resource and accuracy acceptance

The gate also runs `test-resource-contracts`, the Python resource/lifecycle suite,
and retained exact-subsample accuracy cases. See [Resource and scientific contracts](RESOURCE_SCIENTIFIC_CONTRACTS.md) for the budget scope, geometry and timing semantics, accepted engine settings, and explicitly unaccepted legacy diagnostics.

## Capability and maintenance contracts

The engine matrix comes from `capabilities/engines.json`. Builds check generated
outputs before compiling. The gate retains its oracle plan and conservative
regression selection, checks direct runtime ownership and two-rank MPI lifecycle
and consensus, and runs cold-process phase/peak-RSS/accuracy benchmarks. See
[MAINTAINABILITY_CONTRACTS.md](MAINTAINABILITY_CONTRACTS.md) and
[ENGINE_CAPABILITIES.md](ENGINE_CAPABILITIES.md).
