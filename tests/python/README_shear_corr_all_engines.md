# Full-sky shear all-engines driver

All three drivers share `benchmark_timing.py`. Wall time uses a monotonic clock;
CPU time includes the process's worker threads. MPI wall time takes the maximum
rank measurement; CPU time sums ranks. Total wall time is the maximum *per-rank
setup-plus-compute total*, not a sum of separate phase maxima. Invalid clock
measurements, duplicate ranks and inconsistent timing scopes are rejected.
Native `MainLoop`, complete Python `Run`, and launcher times stay separate.


The timing table separates `mainloop_wall_s` / `mainloop_cpu_s` (native
`MainLoop`, excluding Python provenance capture) from `compute_*` (the complete
Python `Run` call). JSON stores these as `native_mainloop_wall_time` and
`native_mainloop_cpu_time`. MPI uses maximum rank wall time and summed rank CPU
for both scopes; use the same scope, output settings and accuracy for comparisons.
Build the matching current Cython extension before running this driver.

`tests/python/shear_corr_all_engines.py` runs one retained spin-2 catalog through
the six full-sky dual node methods enabled in the public Makefile profile:

- `octree-shear-sphere-2balls-omp`
- `kdtree-shear-sphere-2balls-omp`
- `balltree-shear-sphere-2balls-omp`
- `octree-shear-sphere-2balls-mpi`
- `kdtree-shear-sphere-2balls-mpi`
- `balltree-shear-sphere-2balls-mpi`

All methods use observer-centered unit vectors, local east/north shear
components, great-circle parallel transport, and the same cTreeBalls
2PCF/3PCF output contract. Disabled addons are excluded from this public profile.

## Build

Build the active profile before running the driver:

```sh
make -j4 all
```

Use `--list-engines` to check the executable and Python extension selected by
the current environment.

The native executable is the authoritative source for its compiled methods,
runtime options, and resolved numerical backends:

```sh
./cballs options=print-search-methods
./cballs options=print-options
./cballs options=make-info
```

`print-options` now includes the scalar PCA ball-tree cache, sparse-frontier,
parallel-build, leaf-policy, and phase-profiling controls. `make-info` reports
the resolved SLEEF state, actual vector-log backend, active ball-tree policies,
and spherical shear input, accuracy, and MPI policies. The help distinguishes
Python FITS-reader arguments from native options; DES input does not add a
native `infileformat` value. The driver runs all three probes once and stores their return codes,
full text, parsed method/option names, and Make settings in
`summary.json["ctreeballs_runtime"]`.

## Examples

Compare every compiled compatible method on a HEALPix map:

```sh
python3 tests/python/shear_corr_all_engines.py \
  --fits catalogs/shear.fits --geometry sphere \
  --engine all --statistics both \
  --sep-units arcmin --min-sep 3 --max-sep 120 \
  --threads 4 --nsmooth 16 --no-smooth-pivot \
  --outdir Output_shear_all
```

`--threads` is the OpenMP thread count. `--engine all` only selects methods
found in both the executable and extension. Use `--engine all-omp` or `--engine all-mpi` to select one family. MPI shear
addons are enabled and shipped; `--mpi-ranks` sets the rank count and
`--threads` sets threads per rank. The driver launches MPI workers itself.
Do not wrap the parent driver in `mpiexec`.

Use separate `--statistics 2pcf` and `--statistics 3pcf` runs for performance
measurements. `both` exercises the combined production path and is useful for
cross-engine output checks.

For a first large-map test, add `--max-points 100000`. This deterministic
selection is useful for timing and regression checks but is not the full
survey estimator.

With `SMOOTHPIVOTON=1`, capable paths enable pivot smoothing by default;
`--no-smooth-pivot` disables it. A literal `--smooth-radius` is specified in
arcminutes and must obey `2*rsmooth <= min-sep`. Use `--tree-theta 0 --exact-tree --no-smooth-pivot` for an unsmoothed
body reference accepted by the numerical qualification API. `--nsmooth` controls the native
leaf/smoothing capacity; `--more-options NAME[,NAME...]` adds native options.

OpenMP native octree, KD-tree, and ball-tree 3PCF support opt-in `--more-options shear-pivot-reuse`
with `--no-smooth-pivot`, `BALLS4SCANLEVON=1`, positive `theta` and full pivot
coverage. All MPI shear methods reject this option. MPI shear also rejects
`legacy-one-ball`. `CBALLS_SHEAR_PIVOT_TOL` sets a phase
budget in radians in `[0,3]`, default `0.1`, not a percentage error guarantee.
With reuse active the budget controls 3PCF acceptance instead of the magnitude
of `theta`. Zero budget only disables reuse, not ordinary neighbor acceptance.
Use `no-one-ball` for an exact native-octree 3PCF. Calibrate against exact
results for the same catalog before choosing a production budget.

FITS input must provide `GAMMA1` and `GAMMA2`, or selectors supplied with
`--gamma1-field` and `--gamma2-field`. `--mask` is applied before catalog
registration. NPZ input uses `positions`, `gamma1`, `gamma2`, optional
`weights`, and optional `geometry`.

## Sparse DES/Takahashi FITS tables

`--fits-format auto` (the default) detects an x/y/z catalog table; use
`--fits-format desy3` or `healpix` to select the layout explicitly. DES input
requires Astropy and one table with scalar float32/float64 `x`, `y`, `z`, `gamma1`, and
`gamma2` columns (case-insensitive). `kappa` and other columns may be present
but are not used for shear. Multiple coordinate tables are rejected as
ambiguous. Directions come from the stored rows and are normalized; they are
not reconstructed from NSIDE or padded into a dense map.

The sign policy matches the DES workflow:

| `G2CONV` header | Input treatment |
| --- | --- |
| `T17RAW` | Keep gamma1, negate gamma2 once. |
| `LOCALEN` or `EASTNORTH` | Keep both components. |
| Missing or unknown | Require `--des-shear-convention takahashi` or `local-east-north`. |

The result is `gamma1+i*gamma2` in the local east/north basis. The table
already encodes its footprint and uses unit pixel weights. `WTSUM`, `G1MEAN`
and `G2MEAN` are retained as provenance, never used to divide or recenter the
fields. Extra masks, custom weights/field selectors, flat geometry and
`--conjugate-input-shear` are rejected for this format. REALIZ/TOMOBIN/REGION
are checked against standard filenames, and NPIXCAT against the row count.

```sh
python3 tests/python/shear_corr_all_engines.py \
  --fits catalogs/DESY3_Takahashi_r002_bin1_region1.fits \
  --fits-format desy3 --engine all-omp --statistics both \
  --max-points 4096 --sampling-seed 8675309 --fits-chunk-rows 65536 \
  --binning sofia-fig1 --multipoles 7 --phi-bins 32 \
  --tree-theta 0 --exact-tree --no-smooth-pivot --threads 4 \
  --outdir results/desy3-shear
```

`--max-points 0` retains all valid rows. Otherwise SplitMix64 bottom-k row
selection is repeatable across chunk sizes; invalid directions and nonfinite
shear are removed before selection. `--fits-chunk-pixels` is an alias for
`--fits-chunk-rows`. Selection memory is proportional to the chunk plus the
retained sample; thinning changes the catalog and is not a full-survey result.

`--binning sofia-fig1` reproduces the supplied workflow's 20 logarithmic chord
bins whose first/last centers correspond to 8/200 arcmin. The distinct
`paper-8-200-edges` preset uses 8/200 arcmin as outer edges. Presets override
min/max separation, unit and bin-count arguments, and reject linear bins.
`custom` preserves explicit settings. Angular edge/center arrays are saved.

For deterministic realization/bin/region jobs, including array-task selection:

```sh
python3 tests/python/run_desy3_shear_catalogs.py \
  --catalog-root catalogs --realizations 2 --tomobins 1:4 --regions 1:4 --list
python3 tests/python/run_desy3_shear_catalogs.py \
  --catalog-root catalogs --realizations 2 --tomobins 1:4 --regions 1:4 \
  --task-index 0 --outdir results/desy3-batch -- \
  --engine all-omp --max-points 4096 --threads 4 --no-smooth-pivot
```

Ranges are inclusive; task indices are zero-based in realization/bin/region
order. Missing or duplicate identities fail before launch. The batch runner
uses a fresh process per job and writes a completion marker only after all
requested engines succeed. `--resume` checks input size/mtime, arguments,
script and extension hashes, and retained output hashes. Full input content
hashing is not part of this resume policy.

## Numerical qualification and retained validity masks

Every run publishes a qualification packet and prints `UNMEASURED` until
compared with a matching exact-reference packet. `[done]` only confirms
execution. To qualify, first retain an identical catalog/binning/multipole run
with `--tree-theta 0 --exact-tree --no-smooth-pivot`, then pass its output root
as `--qualification-reference` on the candidate run. The tolerances are
`--qualification-rtol` (default 0.02) and `--qualification-atol` (1e-10).
Acceptance is per-observable L2 error; pointwise deviations are diagnostics,
not a guarantee that every bin satisfies the relative tolerance. Even a
reference run initially reports UNMEASURED. Different catalogs or NSIDEs need
separate evidence.

Histogram packets additionally retain `pair_valid`, raw `upsilon`, `window`,
`normalized=upsilon/window_monopole`, native corrected `gamma`, window condition
numbers and relative solve residuals. `gamma_valid` requires a positive real
window monopole, condition number <=1e10 and residual <=1e-8. This mask is a
linear-solve diagnostic, not approximation qualification. Raw/corrected arrays
are retained unchanged; built-in plots mask invalid bins, and other consumers
should explicitly apply the masks. DES
packets use the supplied workflow's `desy3-shear-v1` schema and retain catalog
identity/sign metadata in the companion JSON.

## Outputs and timing

Each engine writes a compressed histogram file beneath `--outdir`.
`summary.json` contains the numerical comparison, native runtime-help snapshot,
and full timing metadata, native parameters, leaf capacity and requested reuse budget;
`timing_report.txt` presents setup, compute wall, process CPU, rank count, and
threads per rank. MPI compute wall is the slowest rank and MPI CPU is summed
over ranks. Launcher startup and serialized catalog transfer are recorded but
excluded from the native compute timing.

The driver plots `xi+`, `xi-`, corrected Gamma-x components, and flattened
3PCF multipole views. `--no-plots` disables figures. Set
`CBALLS_SHEAR_PROFILE=1` to include native per-thread geometry, radial lookup,
ring accumulation, scheduler, and MPI-reduction diagnostics in the console or
MPI worker log.

For scaling, repeat a fixed catalog and statistic at each thread count with
separate output directories. Compare dual node methods at matched numerical
error, not merely equal opening parameters. Private benchmark environments
are not bundled in this public checkout.

## Optional external dual-node API

The all-engines driver times native cTreeBalls engines. Separate external
comparisons can use [dual_node_compat.py](dual_node_compat.py), which lazily
loads the installed external library and translates `bin_theta` and
`angle_theta` without changing their values. It is optional; no external
package is required for the native driver. For example:

```python
from dual_node_compat import catalog, correlation

# Use identical catalogs, units, weights, bins and estimator conventions.
ref = correlation("GG", min_sep=0.01, max_sep=1.0, nbins=10,
                  bin_theta=0.0, angle_theta=0.0)
# ref.process(catalog(...), num_threads=4)
```

`bin_theta`/`angle_theta` are geometric tolerances, not relative-error bounds.
The native `theta` parameter retains its existing meaning.


### Hierarchical 3PCF reuse

All three OpenMP spherical two-ball shear engines support `shear-pivot-reuse`.
Partial multipole rings retain their original spherical acceptance frames.
An unresolved neighbor marks only the radial bins that its enclosing distance
interval can intersect. Once both legs of a radial pair are complete, that pair
is accumulated at the current pivot cell. Descendants inherit a completion mask
and never accumulate that pair again. Mixed resolved/unresolved pairs wait until
both legs are complete; diagonal self-neighbor subtraction is retained.

The option is off by default and requires `no-smooth-pivot`, positive `theta`,
full pivot coverage, and `BALLS4SCANLEVON=1`. Exact controls and smoothing retain
the body-pivot fallback. MPI shear methods continue to reject this experimental
option. The independent 2PCF pass retains its own acceptance rules.

- `CBALLS_SHEAR_PIVOT_TOL`: phase budget in radians, finite `[0,3]`, default `0.1`.
- `CBALLS_SHEAR_BIN_THETA`: internal radial-bin assignment allowance, finite `[0,1]`
  bin widths, default `0`. Zero requires complete containment inside a bin.
  Positive values allow center-based assignment when the combined cap radius
  fits within the selected fraction of the local bin width. The minimum and
  maximum separation cuts remain strict. This is not an exact translation of
  dual node's bin-slop rule.
- `nsmooth=1`: exposes finer pivot groups in the KD and PCA ball trees. It increases
  tree storage and can improve reuse; the octree already has individual-body leaves.

For a performance/accuracy trial through the Python drivers, add
`--nsmooth 1 --more-options shear-pivot-reuse` and explicitly set the two environment
controls, for example `CBALLS_SHEAR_PIVOT_TOL=3 CBALLS_SHEAR_BIN_THETA=0.5`.
These are approximate trial settings, not an accuracy guarantee. Compare raw
numerators, windows, and corrected complex multipoles against an exact reference
for the actual catalog. Poorly conditioned windows can amplify small raw errors.

Native provenance records `shear_hierarchical_reuse.enabled`,
`phase_budget_radians`, and `radial_bin_slop` from the completed run. With
`CBALLS_SHEAR_PROFILE=1`, `radial_pairs` counts actual radial-pair combinations,
`represented_pairs` counts the equivalent individual-pivot combinations, and
`partial_reductions` identifies completion above unresolved descendants. The
64 MiB per-worker scratch limit includes completion masks. Reuse walk timers
include reductions and must not be added to the reduction timer.
