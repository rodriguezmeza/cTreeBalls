# Development guidelines: engines, execution and measurements

The active engine set is resolved from the saved Makefiles and recorded in the
release checks’s `active-profile.md`, with independent numerical references.

## Declaring engines

`capabilities/engines.json` supplies stable IDs, compile conditions, aliases,
geometry/estimator help, source ownership and reference calculations and test declarations. Run
`python scripts/generate_capabilities.py` after editing it. Retain the generated
C table, legacy registration fragments, Python gate plan and
`ENGINE_CAPABILITIES.md`. Every native build checks generated outputs and rejects
stale files. The build fingerprint includes the capability JSON.

Backends and dispatch still need implementation. Registration comes only from
a capability entry. Public engines must declare an existing independent reference calculation
family and reference modules. Legacy entries explicitly declare missing public
coverage; enabling one makes the public release checks fail until a reference calculation is supplied.
Aliases resolve to canonical IDs without inflating the test matrix. Reference
numerical formulae are never generated. Optional/development entries do not imply
public accuracy support or verification of all alternate builds.

## Selecting regressions

`make affected-regressions REGRESSION_ARGS='--changed source/smooth_pivots.c'`
writes selected engines, reference modules and reasons. Use `--baseline` with a
retained build-fingerprint.json on non-Git trees, or `--base REF` for tracked
changes against a Git base. Quoted C includes are followed in reverse,
transitively, including shared addon implementation headers. Unknown paths and
top-level shared definitions select the complete active matrix. CI retains the
selection and runs the complete public release checks for selected active engines.
Optional-profile references are listed separately with required build overrides;
they are not promoted into the public release checks merely by being listed.
Selection never shortens release verification.

## Runtime ownership and file boundaries

- `runtime_context.c`: activation, guarded entry points and context lifetime.
- `memory_catalog.c`: validation and copying of embedded catalogs.
- `common_histogram.c`: common scalar allocation, destruction and normalization.
- `smooth_pivots.c`: deterministic smoothing, grouping and accumulation.
- `engine_registry.c`: lookup/help using generated capability data.
- `mpi_runtime.c`: one MPI owner and shared failure consensus.

The extracted histogram, smoothing and catalog routine bodies preserve their
arithmetic and iteration order. RNG, GSL histogram workspace and six I/O column
pointers now belong directly to the active context. Native-octree construction
counters, radius histograms and diagnostic scratch also live there. Each of
13 MPI backends has a context-local active/rank/size view. Native context tests
and Python interleaving/recovery tests protect these ownership boundaries.

Catalogue/tree pointer globals still use the serialized save/restore adapter
because OpenMP clauses name them. The current-context selector, GSL error
handler and parameter parser retain process-level
constraints. Concurrent native entry is **not** supported; Python retains the
GIL. These are verified migration steps, not a reentrancy claim. Rebuild C
consumers with the matching headers and library.

## MPI execution rules

Only `mpi_runtime.c` initializes/finalizes MPI and implements error consensus.
All 13 backend wrappers delegate to it. One process owner is tracked; MPI supplied
by mpi4py or another host is never finalized by cTreeBalls. FUNNELED support and
MPI main-thread entry are required. Calls after finalization fail. All ranks
receive the first failed rank's diagnostic when a stage fails.

MPI_COMM_WORLD and MPI_ERRORS_RETURN remain the policy. Backend task ownership,
reduction datatypes, chunks and summation order are unchanged. Two real ranks
exercise owned/borrowed lifecycles, local failure, consensus, recovery and
finalized rejection. Retained engine reference calculations check distributed science results.

## Performance, memory and accuracy

```
python scripts/benchmark_accuracy_and_timing.py --output /new/benchmark-directory --repeats 3
```

Defaults cover scalar KD, ball and octree engines on two 512-point exact
subsamples from `accuracy_acceptance.py`. `--engines` also accepts the three
validated spherical-shear engines. Approximation and supported clustered-field
smoothing use the existing fixed tolerances. No tolerance is fitted to results.

Each sample runs alone in a fresh process with a cold tree cache. Records retain
build/extension identity, fixtures, products, configuration, wall/process CPU
seconds for construction, registration, initialization, MainLoop, extraction and
cleanup, plus process peak RSS in bytes. Available native build, frontier,
traversal and reduction timers are retained. Native `*_thread` timers sum thread
work; they cannot be added to elapsed `*_wall` timers.

RSS includes Python, imports, fixtures and C allocations. Phase snapshots are
cumulative process high-water marks, **not** independent phase peaks or totals
of live allocations. Import/fixture RSS is recorded separately. Product
serialization follows the measured workload peak. Darwin bytes and Linux KiB
are normalized explicitly.

Exact/approximate arrays, pair/3PCF L2 errors and acceptance limits accompany
phase samples and median/min/max times. Exact repetitions use the stricter
exact tolerance; approximation uses the existing 2% Frobenius L2 plus absolute
floor as the acceptance criterion. Accuracy failure fails the benchmark and checks for the selected build. These
small workloads establish reproducible baselines, not large-survey scalability.
Timing variation is reported; this refactor claims no speedup. Compare matching
machines, workloads, threads and settings when judging future optimizations.

The native-octree and PCA compact-tree caches are now runtime-object attachments.
Explicit cleanup/destruction frees them; internal parameter recomputation can
reuse them within the same object. See RELEASE_WORK_SEQUENCE.md for lifecycle
and concurrency guarantees and the extracted radial/acceptance/scheduling owners.

## Incremental engineering contracts

See [engineering contracts](docs/engineering_contracts.rst) for numerical helper
owners, shared scalar parsing, the process-wide serialized native API contract,
acceptance-policy identity, platform/dependency receipts and the dedicated-node
benchmark workflow. No GIL release or native reentrancy is introduced.
