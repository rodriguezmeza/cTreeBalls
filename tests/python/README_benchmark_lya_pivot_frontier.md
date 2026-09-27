# Forest-local pivot smoothing and spatial scan tasks

The six 3D OpenMP Ly-alpha hard-bin methods accept two independent acceleration
controls inspired by legacy `SMOOTHPIVOT` and `BALLS4SCANLEV`. They are available
through native parameters, parameter files and `cyballs.set(...)`.

| Parameter | Default | Meaning |
|---|---:|---|
| `lyaScanLevel` | 0 | 0 keeps catalog/group order; 1..20 partitions pivots at an octree depth into spatial tasks. |
| `lyaPivotRadius` | 0 | Positive values replace nearby pivots **from the same forest** by representative geometry. Units match the Cartesian catalog. |
| `lyaPivotMax` | 8 | Maximum original pixels represented by one pivot; valid 1..1024. |
| `lya3PivotBlock` | 0 | 3PCF block-size override, 1..4096; 0 uses the existing automatic block size. It caps 3PCF spatial tasks too; pair-only tasks retain a cap of 64. |

Defaults preserve the exact estimator. A scan level with radius zero changes
work ordering, not pixel geometry, binning, or forest exclusions. Raw sums can
differ at rounding level from catalog order; fixed task boundaries and ordered
publication keep results deterministic across thread counts for a fixed build
and settings. Task execution is dynamically scheduled. Each task still walks
the neighbor tree separately for its representative pivots; the frontier is
not a dual-tree cell-pair kernel and does not guarantee a speedup.

## Smoothing contract

Active, valid pivots are sorted by forest ID, observer distance and original row
ID. Starting at the first unused pixel, a group includes consecutive pixels
from that forest within `lyaPivotRadius` of the first pixel, up to
`lyaPivotMax`. Measured 3D distance is used, including bent forests. The first
pixel's position, distance and LOS become the representative geometry. Catalog
positions, fields, weights, masks and update flags are never changed. All
neighbor pixels remain individual original pixels.

For a group G, the code stores `W_G = sum(w_i)` and `Q_G = sum(w_i * delta_i)`.
For each accepted neighbor pair j,k from two different other forests, it adds
`Q_G * (w_j delta_j) * (w_k delta_k)` to the 3PCF numerator and
`W_G * w_j * w_k` to its denominator, depositing both ordered-side permutations.
Geometric triplet counts increase by `2 * |G|`. Thus weights and forest
exclusions are preserved; geometry is approximated.

For 2PCF, original row-ID ownership remains `Id(i) < Id(j)`. If a group straddles
neighbor j's ID, only eligible members contribute to its W, Q and multiplicity.
This handles shuffled catalogs without lost or duplicate original pairs.
Because only pivot geometry is replaced, positive-radius pair results can
depend on input row ordering. Zero-radius results have the usual rounding-only
ordering dependence.

**Positive radius can change bin membership and hard-cutoff membership.** It is
not an error tolerance and does not imply a 5% correlation bound. Cancellation
can amplify even small geometric changes. Use exact output to calibrate every
catalog selection, binning and radius. `lyaPivotMax=1` gives original geometry,
although a positive requested radius is conservatively marked approximate in
metadata. This mode is separate from the certified-cell geometry bounds and
slops of the existing pair/3PCF cell kernels.

`smooth-pivot`, `rsmooth`, `theta`, `SMOOTHPIVOTON` and `BALLS4SCANLEVON` do not
enable these Ly-alpha controls implicitly. No legacy compile flag is required.

## Supported combinations

- Original `lya-2pcf-omp`, `lya-3pcf-omp`, `lya-2pcf-3pcf-omp` and their
  `lya-los-tree-*` counterparts are supported.
- Pair-only uses `lya2Kernel=0`; 3PCF uses `lya3Kernel=0`, 1 or 2.
- Combined runs may use `lya2Kernel=1`: the existing pair-cell pass retains its
  own geometry controls, while the new frontier/smoothing applies to 3PCF.
- MPI, radial-only, anisotropic multipole and persistent 3PCF kernels 3/4 reject
  active frontier/smoothing controls with an explicit error. Those kernels keep
  their existing algorithms and parameters.

The frontier, group members, sort and map buffers use checked O(N) storage.
Their conservative allocation plan is added to the histogram/worker plan and
checked against `CBALLS_MEMORY_BUDGET_MB` before allocation; this is not a
whole-process RSS cap. Group sums must remain finite. Before smoothing, a
conservative full-catalog multiplicity bound is checked against signed INTEGER;
a huge sparse catalog can fail this bound even if few tuples would survive.
Reduce the selection or disable smoothing in that case. All storage is released
on success or error. Repeated in-memory Cython calls are regression-tested.

## Reproducible calibration

Run from the source root, using Python with NumPy; FITS input also requires the
dependencies of the shared `lya_fits.py` reader. Matplotlib generates figures.
Use an otherwise idle host and a **new** output directory:

```bash
python tests/python/benchmark_lya_pivot_frontier.py \
  --fits /path/to/catalogs/lya_15_xyz_raw_with_losid.fits \
  --max-forests 100 --pixel-stride 10 \
  --threads 16 --warmups 1 --repeats 3 \
  --scan-level 2 --radii 0 4 6 12 --pivot-max 8 \
  --outdir results/lya-pivot-calibration
```

The defaults test all three original methods, 160-unit maxima, 50x50 pair bins,
and 4x4x4x4x4 triplet bins. Use `--methods` to choose a subset or LOS methods;
use `--rp-max`, `--rt-max`, `--r3-max`, `--rp-bins`, `--rt-bins`, `--r3-bins`,
`--theta-bins`, and `--mu-bins` to match the production analysis exactly.
`--catalog saved.npz`, `--ascii six-columns.txt` and `--synthetic` are alternatives
to FITS. Forest/stride selection applies only to FITS. Catalog distance units
are passed through, with no implicit Mpc/h conversion.

Start with a small representative selection, then repeat accepted settings on
the intended workload (for example `--max-forests 1000 --pixel-stride 10`).
Do not promote a smoothing radius that fails the small-selection acceptance
test. Radius zero always remains in the sweep as the exact spatial-task case.
Restart Python/Jupyter after rebuilding the Cython extension so its native
parameter layout matches the new executable and headers.

`--baseline-cballs /path/to/before/cballs` measures an earlier executable;
otherwise the exact reference uses the current executable with both new
controls disabled. Build both with the same compiler and profile. Optional
`--include-pair-cells` compares the existing exact pair-cell backend as well.
`--kernel3` and `--reference-kernel3` select kernels 0..2 independently; leave
both at 0 to isolate the new frontier/smoothing cost.

Runs execute serially, interleaved by case per repeat. Wall time includes fresh
native process startup, catalog input, trees/frontier, search and output;
Python input preprocessing is excluded. CPU time is total child user+system
time, never divided by thread count. Peak RSS is the child process high-water
mark from `wait4` when available. Warmups are retained but excluded from medians.
Small pair workloads can be dominated by startup/output, so inspect ranges.

Outputs include the input catalog, binary/catalog hashes, full commands and
logs, native metadata, dense histograms (including empty bins), per-bin errors,
`summary.json`, and timing/speed/error/RSS PNGs. Repeat products must be identical.
The printed frontier counters give active and representative counts, tasks,
actual maximum group radius and build CPU time.

## Acceptance rules

For occupied reference bins with `abs(correlation) > --relative-floor` (default
1e-6), report the maximum and 95th-percentile pointwise relative error. The
default maximum permitted is 0.05. Smaller-signal occupied bins must satisfy
`--small-signal-atol` (default 5e-8). Candidate output must be finite and retain
reference occupancy; newly occupied or lost bins are reported and rejected.
Approximate acceptance also requires at least one eligible bin. Numerator and
denominator relative L2 errors are reported separately.

Exact candidates additionally require raw numerator/denominator agreement
(`rtol=3e-11`, `atol=1e-9`) and identical integer pair/triplet counts. A failed
exact candidate gives exit status 1. Positive-radius failures are recorded but
do not abort the sweep. `--require-accepted-smoothing` gives exit status 2 unless
every selected method has an accepted positive-radius case that actually merges
pixels (a radius smaller than every pixel spacing does not qualify). Acceptance applies
only to the saved workload and stated thresholds, not to another catalog.

Independent regression tests enumerate each original contribution at the
declared representative geometry, without using products of group sums:

```bash
PYTHONPATH="$PWD" python -m pytest -q tests/make_tests/test_lya_pivot_frontier.py
python scripts/generate_capabilities.py --check
```

The test suite covers original and LOS methods, multiple threads and triplet
kernels, shuffled ownership, bent/antipodal forests, zero weights, cancellation,
cutoffs, unsupported controls, memory preflight, Cython reuse and strict
benchmark acceptance. Existing duplicate-position rejection remains in force.
