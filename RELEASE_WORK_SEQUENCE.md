# Release, resource, scientific, ownership and performance contracts

This change set follows the six-stage active-profile review. No saved addon
flag or estimator default is changed. Numerical reference enumerators remain
independent of the production kernels.

## 1. Release verification

The installed-package checker resolves the intended Makefile profile (or accepts
`--expected-build gate/build-fingerprint.json`) and compares exact method/ID
mappings, including unavailable canonical names. It contains no fixed engine
count. The active gate produces `active-profile.md` and fails if any required
module declared by an active engine or regression group was not successfully
executed. The Ly-alpha pair-cell module is included. Source archives must include
the new API, shared kernels, regressions and benchmark drivers.

The installation check runs outside the checkout in a virtual environment and
actually computes three LOS-tree statistics. A venv that reuses site packages
validates the installed native artifact and import isolation; it is not a test
of resolving dependencies on a pristine operating system.

## 2. Resource policy

`CBALLS_DEFAULT_MEMORY_BUDGET_MIB` in `include/resource_contracts.h` is the single
compiled default: **65536 MiB**. `CBALLS_MEMORY_BUDGET_MB` overrides it at runtime.
Set the environment before constructing/running objects, identically on MPI
ranks. No recompilation is needed. Invalid values fail closed.

The limit checks individual common allocations and combined plans. The
anisotropic forest global/worker/scratch plan now includes catalog bodies,
common histograms and retained object caches. It is **not a total RSS cap**.

Before allocating a catalog:

```sh
python scripts/resource_plan.py --engine lya-3pcf-omp \
  --pixels 12582912 --threads 16 --parameters parameters.json
```

`parameters.json` is a JSON dictionary using the native parameter names.
The equivalent installed Python API is:

```python
from cyballs import resource_plan, resource_policy
policy = resource_policy()
plan = resource_plan('lya-3pcf-omp', 12582912, 16,
                     {'lya3RBins': 20, 'lya3ThetaBins': 10, 'lya3MuBins': 20})
```

`known_components_bytes` contains catalog bodies, the common histogram plan,
optional retained cache bytes, and (for anisotropic forest) the global-plus-
worker histogram arrays. `estimated_components_bytes` separates geometry/tree,
radial and backend workspaces whose exact size depends on the workload.
`estimated_total_bytes` is a planning scenario, **not an upper bound**. Both
sums are checked against PTRDIFF_MAX. Plans describe the active run, not a
process-wide ledger of every other live Python object. Interpreter/user arrays, input readers,
MPI/libraries, allocation overhead and data-dependent scratch/frontier/window
solves remain explicit exclusions. Supply all catalogs' combined pixel count
when estimating replicated input memory. Multiply per-rank costs by co-resident
ranks when planning a node. Measure RSS for the intended workload.

`getAllocationInfo()` adds live retained cache/result bytes. Run metadata records
compiled layout sizes, the default and cache scope. `getCacheInfo()` observes
current ownership; historical run metadata does not change after cleanup.

## 3. Scientific qualification

The original acceptance fixtures remain required. `scripts/workload_acceptance.py`
adds held-out masked skies, clustered weak signed scalar/shear fields and
straight/bent forests, using seed 7260926. Each measurement runs alone in a
fresh worker process and retains the fixture, products, explicit settings,
build identity, wall/process CPU time, peak RSS and per-observable errors.

Scalar/shear candidates use theta=0.01, including optional existing smoothing
radii (scalar 0.0003 Cartesian units; spherical shear 1 arcmin). Both raw moments
and window-corrected products are examined. Forest candidates compare direct
pixels, exact persistent pair/pivot cells, and slop=0.01 in every applicable
coordinate. Their raw numerators, denominators and normalized ratios are all
compared. These are candidates, not universally recommended configurations.

Qualification is fixed in advance: per-product L2 error <= 0.05 times reference
L2 + 1e-12, with identical finite masks. Exact aggregation uses 3e-11 + 1e-12.
Maximum bin-relative error, weak-bin absolute error, nonfinite counts and the
signal floor are retained separately. **Qualification does not mean <=5% error
in every bin.** Sparse window conditioning and cancellation can amplify small
geometric changes; a change in the finite mask rejects the candidate.

The qualification campaign may complete with `status=PASS` while candidates
are `REJECTED`: PASS means exact controls passed and all candidates were
measured, not that every candidate is accurate. Only a candidate/observable/
workload combination marked QUALIFIED may be described as qualified by this
campaign. Existing accepted settings retain their original narrowly specified
fixture envelope. Arbitrary theta/slop values, experimental multipole truncation
and pivot smoothing have no general accuracy certificate.

```sh
python scripts/workload_acceptance.py --output /new/qualification-directory
```

## 4. Ownership and result APIs

The two-entry native-octree and PCA-tree caches now belong to each
`cballs_runtime_state`; different objects do not share cached trees. An object
can reuse its content-keyed packed trees across parameter-only recomputations.
`clearCaches()` releases them without invalidating completed products.
`struct_cleanup()` clears results and caches; `struct_cleanup(clear_cache=False)`
explicitly retains that object's compact trees. Destruction releases all owners.
There are no borrowed Python/body pointers in retained compact trees.

Cython objects keep successful publishing-rank forest and physical arrays by
transferring ownership before search-local cleanup. This does not duplicate the
native arrays. `getForestResults()` and `getPhysicalResults()` return independent
NumPy copies with metadata; they work with `options=no-out-Hist`. Getters reject
an uncomputed/cleaned/wrong-family object and non-publishing MPI ranks. Copies
already returned remain usable after cleanup and cannot mutate native state.

```python
result = model.getForestResults()
num = result['arrays']['triple_numerator']
den = result['arrays']['triple_denominator']
model.struct_cleanup()  # num and den remain valid owned NumPy arrays
```

Products and shapes:

| Method family | Arrays | Shape / convention |
| --- | --- | --- |
| Anisotropic forest 2PCF | pair_numerator, pair_denominator | (RpBins, RtBins) |
| Anisotropic forest 3PCF | triple_numerator, triple_denominator | (RBins, RBins, ThetaBins, ThetaBins, MuBins) |
| Experimental anisotropic multipole | moments_numerator, moments_denominator | same first four axes, final axis LMax+1; signed raw moments, not a binwise ratio |
| Radial forest 2PCF | pair_numerator, pair_denominator | (RpBins,) |
| Radial forest 3PCF | triple_numerator, triple_denominator | (2*RBins, 2*RBins), signed lags |
| Same-LOS radial mean | correlation_sum, contributing_forests | (RpBins,); divide the sum of per-forest correlations by the forest count |
| Physical 3D | pair_numerator/denominator; triple_numerator/denominator | (B,), and (LMax+1, B, B) |
| Physical survey | data_ and random_ prefixes on the physical arrays | D-R and random raw products; corrected solutions remain a separate observable |

Floating products return float64, including conversion from long-double radial
accumulators; contributing-forest counts return uint64. Axis edges/conventions
are in metadata. Raw signed moments and survey windows must not be normalized
with an assumed universal numerator/denominator rule.

Native execution still uses a serialized legacy activation adapter, MPI process
ownership and GSL/parser state. Cython retains the GIL across calls. Concurrent
threads may serialize Python calls but **concurrent direct native entry is not
supported**; use processes/MPI for independent concurrent workloads. Moving two
cache owners is an incremental ownership improvement, not a reentrancy claim.

## 5. Implementation boundaries

Common histogram planning/allocation now resides beside its lifecycle in
`source/common_histogram.c`. Shared radial classification, pair acceptance and
work-estimated scheduling, Fourier moments/pivot traversal, and search/backend
policy have separate `dual_node_*.h` components under the
scalar two-ball owner. They retain their specialization context and arithmetic.
They are not alternative reference implementations.

`affected_regressions.py` follows transitive includes; regressions assert that
changes in each extracted component select all six scalar OpenMP/MPI engines
and their independent edge/numerical references. Unknown/shared paths still
select the full active matrix. The full release gate remains the final check.

## 6. Matching-science scaling

```sh
python scripts/benchmark_scaling.py \
  --catalog /path/to/catalogs/lya_15_xyz_raw_with_losid.fits \
  --forest-counts 256 1000 3854 --pixel-stride 10 \
  --threads 1 2 4 --prefixes lya lya-los-tree \
  --statistic 2pcf-3pcf --radius 200 \
  --warmups 1 --repeats 3 --output /new/scaling-directory
```

The existing NPZ catalog is accepted too. Forest prefixes are ascending unique
IDs; pixel selection preserves original within-forest order and uses the same
stride at every size. Retained source-row indices and SHA-256 hashes make the
selection reproducible. The baseline uses existing exact pixel/segment kernels;
optimized mode uses exact persistent pair cells and the selected `--triple-kernel`
(default 0, preserving segment traversal). Kernel 2 tiles direct work, 3 uses
persistent neighbor nodes, and 4 adds pivot cells. These are workload-dependent
candidates: larger cells or more caching do not guarantee speed. Optional `--modes
reference optimized approx --slop 0.01` measures approximation separately.

The script enforces a one-thread exact reference first and compares full products
for the same selected input. It refuses a speed claim for a rejected candidate.
It records all repetitions, explicit wall and process CPU, peak process RSS,
resource forecasts, build/compiler identity, warmups and cache policy. Every
iteration changes rootDir to force computation rather than timing cached Run().
Cold configuration startup and warm repeated execution are distinguished;
forest indexes currently rebuild per call. Workers run sequentially to avoid
interference from this benchmark itself. Other host activity remains possible.

`scaling.json`, `size-scaling.png` and `strong-scaling.png` retain the evidence.
A speedup is an observation for that workload/host/configuration, not a claim
that every exact optimization or thread count is faster. The existing scalar/
shear `benchmark_contracts.py` phase and accuracy benchmark remains available.
