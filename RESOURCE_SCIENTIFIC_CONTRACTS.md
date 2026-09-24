# Resource and scientific contracts

The public 3D profile uses an estimator-specific common histogram plan.
Scalar Fourier 3PCF runs retain the scalar tensor set required by the legacy
kernels and exporters. Pair-only scalar runs and all shear, physical 3D,
periodic-count, and forest runs omit that set and the square Python export
buffer. The unused flattened export buffer is never allocated. Small common
vectors remain available to shared startup and pair routines. Cleanup accepts
partially allocated plans. Scalar 3PCF getters reject an uncomputed product.

For 512 radial bins and mChebyshev=31, the common Python-enabled pair plan
with smoothing compiled contains 57,456 bytes. Previously, the nine scalar
3PCF tensors alone consumed 576 MiB, excluding their pointer tables. This
saving does not include estimator-owned shear/forest products or tree storage.
`getAllocationInfo()` reports the plan and actual live pointer presence;
`getRunMetadata()['resources']` retains its size and budget.

## Checked dimensions and budget

`CBALLS_MEMORY_BUDGET_MB` is a positive integer in MiB, default 1024. Set it
before a run, identically on all MPI ranks; do not mutate process environment
while workers execute. Zero, malformed and overflowing values fail closed.
For example, `CBALLS_MEMORY_BUDGET_MB=256 python analysis.py`.

This is a **per-rank preflight limit**, applied to the complete common histogram
plan and every allocation through the common malloc/calloc, byte, vector,
matrix, tensor, and Numerical Recipes allocators. Matrix/tensor accounting
includes padding and all pointer tables. Smoothing ownership maps and
anisotropic forest global-plus-worker histogram plans have combined preflights.
C-owned in-memory catalog copies also use the checked allocator. Failure returns
an error through guarded calls; checked worker helpers return status for the
engine's existing failure consensus. It is **not a process RSS ceiling** or an
accounting system for every tree, backend, cached allocation, or third-party
library. Multiple independent allocations can exceed it in aggregate.

All common shape arithmetic uses checked size_t addition/multiplication bounded
by PTRDIFF_MAX before conversion or allocation. Legacy signed-index kernels
also impose sizeHistN <= 46340 and bounds on multipole/index products, before
any m+1 or doubled-mode expressions in parameter validation. This may reject
shapes that a future size_t-only kernel could handle safely. Allocation failure
cleanup is preserved. Numeric arrays support zero, one, and negative lower
bounds through zero; invalid/reversed bounds are rejected.

`make test-resource-contracts` exercises LONG_MIN/LONG_MAX/SIZE_MAX boundaries,
compound-shape budget failure, checked calloc/malloc failure, and indexing/free
behavior without requesting enormous allocations. The Python resource suite
checks plan selection, oversize rejection, error recovery, and lifecycle state.

## Retained accuracy envelope

The active release gate invokes `scripts/accuracy_acceptance.py`. It retains
fixtures, exact and approximate arrays, settings, build identity, explicit
metrics, and timings in `accuracy/`. These deterministic **regression acceptance
cases** do not establish a universal error guarantee or a monotonic bound on
all smaller theta values.

Two 512-point weighted subsamples are used: an isotropic sky with a signed
random scalar field and a sky of 64 close groups with a smoothly varying scalar
and shear field. Weights span 0.5 to 1.5. The setup uses six logarithmic chord
bins from .02 to 1.8, scalar modes 0..5, and two OpenMP threads. Each result is
compared to the same engine's body-exact, unsmoothed run. Exact production
results must also agree across the three tree families to 1e-11 relative L2,
plus 1e-10 absolute L2. Existing independent enumeration oracles remain in the
gate; the subsample tests add coverage at larger, nontrivial tree sizes.

| Engine | Accepted theta setting | Pair and 3PCF relative L2 limit | Smoothing case |
| --- | ---: | ---: | --- |
| kdtree-2balls-omp | 0.5 | 0.02 | literal Cartesian radius 0.0003 |
| balltree-2balls-omp | 0.2 | 0.02 | literal Cartesian radius 0.0003 |
| octree-2balls-omp | 0.025 | 0.02 | unsupported in production mode |
| octree-shear-sphere-2balls-omp | 0.05 | 0.02 | 1 arcmin |
| kdtree-shear-sphere-2balls-omp | 0.05 | 0.02 | 1 arcmin |
| balltree-shear-sphere-2balls-omp | 0.05 | 0.02 | 1 arcmin |

The criterion for each product is norm(approx-exact) <= 0.02*norm(exact)+1e-10,
using a complex Frobenius norm. It is evaluated separately for the pair vector
and the full 3PCF array. Scalar products are the normalized pair estimator and
raw complex Fourier moments. Shear products are xi-plus/xi-minus and **raw
Upsilon** multipoles. Window-corrected Gamma is excluded here: sparse-window
conditioning is a different acceptance problem and has separate exact-oracle
coverage. These tolerances are not binwise relative errors or covariance-based
significance limits. Smoothing is compared against unsmoothed exact results on
the clustered fixture, so approximation and smoothing errors are both included.

The suite also retains explicitly marked **diagnostic-only**, unaccepted runs
for octree-sincos-omp and the octree-GGG compatibility path. Their clustered
results exceeded the 2% criterion, including large pair errors and an empty-bin
NaN in compatibility mode. They are outside this validated accuracy envelope;
the gate does not silently call them accurate or raise the tolerance to pass.
Use production engines and validate representative exact subsamples before
choosing a different setting. Experimental shear-pivot-reuse, corrected survey
windows, alternate precision, and other datasets need their own calibration.

Scalar smoothing now includes the representative body's own weight in its
weighted field sum, matching every claimed neighbor. A constant-field,
nonuniform-weight regression protects this correction. The core
`octree-sincos-omp` method does not implement `only-2pcf` and now rejects that
previously ignored request; select `octree-2balls-omp` for a true pair-only run.

## Lifecycle, timing, and geometry

`state` is true only after successful MainLoop while native results remain
live. Fresh, input/startup-only, failed, invalidated, and cleaned objects report
false. A retained immutable run snapshot can outlive native arrays.

`Run()` keeps its historical numeric return for compatibility: MainLoop process
CPU seconds divided by requested threads. This is neither elapsed wall time
nor total CPU consumption. `getTimings()` and Python run provenance instead
report explicitly named MainLoop `wall_seconds` (monotonic perf_counter),
`process_cpu_seconds` (process_time), requested threads, and the legacy value.
Cached runs retain the original measurement; parameter/catalog changes
invalidate it. `getCPUTime()` is the native search-stage process CPU counter.
These are local-rank quantities; job wall time requires a rank maximum and job
CPU consumption requires a rank sum, which these getters do not perform.

Metadata distinguishes planar scalar geometry, flat-sky shear, full-sky shear,
physical 3D, forest, and periodic boxes. It records compiled dimensions and
spherical-shear/observer-frame flags. Requested and effective smoothing radius
units are explicit: full-sky shear accepts arcminutes and stores chord distance;
scalar/flat-sky methods use Cartesian catalog units. A 2D profile requires
CFITSIO and the 3D-only estimators disabled. Its geometry checks are separate
from the active public 3D/double release gate.

Legacy ball-tree pair-only runs also omit worker angular vectors, matrices and tensors, angular accumulation, reductions, and 3PCF output. The native allocation test enforces this under a 1 MiB limit. Forest provenance emits only computed axes and preflights JSON expansion before allocating it.
