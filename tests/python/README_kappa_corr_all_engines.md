# Convergence all-engines driver

All three drivers share `benchmark_timing.py`. Wall time uses a monotonic clock;
CPU time includes the process's worker threads. MPI wall time takes the maximum
rank measurement; CPU time sums ranks. Total wall time is the maximum *per-rank
setup-plus-compute total*, not a sum of separate phase maxima. Invalid clock
measurements, duplicate ranks and inconsistent timing scopes are rejected.
Native `MainLoop`, complete Python `Run`, and launcher times stay separate.


The timing table separates `mainloop_wall_s` / `mainloop_cpu_s` (native
`MainLoop`, excluding the collection of run information in Python) from `compute_*` (the complete
Python `Run` call). JSON stores these as `native_mainloop_wall_time` and
`native_mainloop_cpu_time`. MPI uses maximum rank wall time and summed rank CPU
for both scopes; use the same scope, output settings and accuracy for comparisons.
Build the matching current Cython extension before running this driver.

The runnable drivers and their guides live in `tests/python`; `python/` is
reserved for the Cython binding sources. Run commands from the repository root.

`timing_report.txt` separates setup and compute wall time from process CPU.
MPI wall time is the maximum across participating ranks and CPU time is their
sum. `summary.json` retains each rank's times and native parameters. Catalog
registration, Python extraction/plotting, and cleanup are outside these timers;
the native MainLoop includes its requested output. Compare identical statistics,
bins, masks and calibrated accuracy, not just identical `theta` values.

`kappa_corr_all_engines.py` reads a convergence catalog once, registers its
NumPy arrays with `cyballs`, and runs the active scalar angular engines. The
current profile exposes:

- `kdtree-2balls-omp` and `kdtree-2balls-mpi`
- `balltree-2balls-omp` and `balltree-2balls-mpi`
- `octree-2balls-omp` and `octree-2balls-mpi`

These are native cTreeBalls engines. Their production traversal uses dual-node
acceptance; `--dual-node-bin-theta` enables the looser bin-aware policy, and
`--more-options no-one-ball,no-two-balls --no-smooth-pivot` requests an exact body-level
validation walk. For native octree 3PCF, `no-two-balls` alone is not an exact
reference. `--more-options legacy-one-ball` selects the privately linked compatibility
kernel, whose smoothing behavior differs from the default octree traversal.

## Examples

Run all compiled OpenMP engines on a FITS map:

```sh
python3 tests/python/kappa_corr_all_engines.py \
  --fits Tests/catalogs/allskymap_nres12r081_zs9_mag.fits \
  --engine all-omp --threads 16 --outdir Output_kappa
```

Use a mask and complex scalar 3PCF edge correction:

```sh
python3 tests/python/kappa_corr_all_engines.py \
  --fits Tests/catalogs/allskymap_nres12r081_zs9_mag.fits \
  --mask Tests/catalogs/mask_octant.fits \
  --engine octree-2balls-omp,kdtree-2balls-omp \
  --edge-corrections --threads 16 --outdir Output_kappa_masked
```

Run an MPI engine through the driver launcher:

```sh
python3 tests/python/kappa_corr_all_engines.py \
  --catalog-npz catalog.npz --engine octree-2balls-mpi \
  --mpi-ranks 2 --threads 4 --outdir Output_kappa_mpi
```

Use repeated `--mpi-extra-arg` options when the local launcher needs them.
`mpi4py`, `mpiexec`, and the extension must use the same MPI implementation.

## Angular patches

Use `--patch` to apply the same sky cut to the shared Python catalog before
any engine receives it. The existing spelling `--more-options patch` is
equivalent. Specifying bounds alone does not enable the filter, so old
commands without `patch` retain their previous full-catalog behavior.

```sh
python3 tests/python/kappa_corr_all_engines.py \
  --fits map.fits \
  --patch --phiL 73.344728 --phiR 106.677236 \
  --thetaL 73.334472 --thetaR 106.665490 \
  --engine octree-2balls-omp --statistics 3pcf \
  --no-smooth-pivot --threads 16 --outdir Output_kappa_patch
```

All four bounds are in **degrees**, independently of `--theta-scale`.
`phi` is longitude; `theta` is **colatitude**, measured from the north
pole, not declination (`declination = 90 - theta`). Selection exactly uses
the native patch inequalities:

```text
phiL < phi < phiR
thetaL < theta < thetaR
```

Edges are excluded. Bounds must satisfy `0 <= phiL < phiR <= 360` and
`0 <= thetaL < thetaR <= 180`; longitude wraparound is not supported.
The existing defaults (`0,90` for each interval) select the northern
first-longitude quadrant only when `--patch` is enabled.

For FITS input the order is resolution conversion, invalid-pixel/mask
selection, patch selection, `--max-points` thinning, and mean subtraction.
The patch intersects `--mask`; it does not replace it. After `--nside-down`,
the cut uses the resulting pixel centers. Native-resolution RING and NESTED
maps remain memory-mapped and are filtered in chunks. Use
`--no-center-field` to preserve input kappa values instead of subtracting
the mean of the retained sample.

NPZ and synthetic catalogs also support the patch. Positions are not
translated or rotated; weights and masks remain aligned with selected rows.
In the Python API, pass `patch=AngularPatch(phiL, phiR, thetaL, thetaR)`
to a catalog loader, or use `RunConfig(patch=True, phi_left=...,
phi_right=..., theta_left=..., theta_right=...)` with `run_engine_suite`.
For already-centered API catalogs, retain `metadata["centered"] = True`
to recenter the selected active sample. For FITS thinning, provide the patch
to the loader so selection precedes sampling.

Rank 0 filters before saving `--save-catalog-npz` or broadcasting to MPI
ranks. The terminal reports before/after counts, and the catalog metadata in
`summary.json` records the patch and selected population.
`patch-with-all` is a different, pivot-only selection and cannot be combined
with this shared-catalog filter.


## Input and results

`--fits` accepts an implicit HEALPix scalar map. `--catalog-npz` accepts an NPZ
file containing `positions`, `kappa`, and optional `weights` and `mask` arrays.
`--max-points` performs deterministic thinning after mask selection, which is
useful for smoke tests but changes the measured catalog. `--nside-down 0`
retains native resolution.

All active scalar engines support masks and edge correction. Masking controls
which bodies enter the estimator. Edge correction is a separate 3PCF window
deconvolution and cannot be combined with `--statistics 2pcf`.

Each engine writes into its own subdirectory. `summary.json` records input
metadata, compile settings, wall/CPU timing, and pairwise numerical
comparisons. Plotting produces ordinary multipole views and flattened radial
matrices; use `--no-plots` or `--no-flatten-plots` to disable them.

The catalog is read once and retained between runs, but every engine still
builds its own C-owned tree and histograms.

## Optional external dual-node API

The all-engines driver times native cTreeBalls engines. Separate external
comparisons can use [dual_node_compat.py](../python/dual_node_compat.py), which lazily
loads the installed external library and translates `bin_theta` and
`angle_theta` without changing their values. It is optional; no external
package is required for the native driver. For example:

```python
from dual_node_compat import catalog, correlation

# Use identical catalogs, units, weights, bins and estimator conventions.
ref = correlation("KK", min_sep=0.01, max_sep=1.0, nbins=10,
                  bin_theta=0.0, angle_theta=0.0)
# ref.process(catalog(...), num_threads=4)
```

`bin_theta`/`angle_theta` are geometric tolerances, not relative-error bounds.
The native `theta` parameter retains its existing meaning.

## Optional hierarchical scalar 3PCF

For the scalar `octree-2balls`, `kdtree-2balls`, and `balltree-2balls` engines
(OpenMP or MPI), append `--more-options scalar-pivot-reuse --no-smooth-pivot`.
Select these engines explicitly; other engine families reject this option.
`CBALLS_SCALAR_PIVOT_TOL` (default 0 radians) and
`CBALLS_SCALAR_BIN_THETA` (default 0 bin widths) are runtime environment
variables; MPI ranks must receive matching values. Benchmark and qualify
the actual scalar products before choosing looser values.
The qualification packet retains the native hierarchy controls and counters.
See [the algorithm, limits, and benchmark instructions](../../docs/SCALAR_HIERARCHICAL_REUSE.md).

### Explicit separation limits and plotted coordinates

The existing `--theta-min`, `--theta-max`, and `--theta-scale` names remain
supported. Their aliases are `--min-sep`, `--max-sep`, and `--sep-units`.
Units may be `degree`, `arcmin`, `radian`, or native distances (`au`, also
spelled `chord`). `au` is the historical native-distance spelling, not an
astronomical-unit conversion. Defaults remain unchanged: about 7.35028–217.67984
arcmin, expressed as 0.12250467913471644–3.6279974066581295 degrees, and 20 bins.
For example, `--min-sep 2.5 --max-sep 250 --sep-units arcmin --nbins 20`
uses 2.5 and 250 arcmin as outer edges for both statistics and both 3PCF legs.
Angular values are converted to sphere chords before correlation; plotted
coordinates are converted back to the requested units.

2PCF figures show those explicit edges and use logarithmic axes for logarithmic
bins. 3PCF heat maps show actual radial edges in the selected units on both axes.
Flattened 3PCF plots retain an index axis; `plots/ENGINE_radial_pairs.csv` maps
each index to its two actual radial centers. Plot titles state the measured
range. No correlation coefficients or bin assignments are changed by plotting.
NPZ products now retain `radial_edges`, and angular runs also retain
`chord_edges`, `theta_edges_arcmin` and `theta_arcmin`. Native-distance runs
need not be spherical, so angular metadata is not invented for them.
`summary.json` and per-engine JSON record `separation_binning`.

A measurement made with different edges must be recomputed; changing axis limits
cannot reconstruct missing separations. No C recompilation is needed for these
Python changes.

Benchmark controls and timing (October 2026)

The shared benchmark_timing.py records explicit CBALLS_SHEAR_*,
CBALLS_SCALAR_*, CBALLS_LYA_* and thread-placement environment overrides.
An omitted environment variable is not assumed to be a compiled default;
read native run information for the effective settings and build identity.

compute_wall_time/compute_cpu_time cover the complete Python Run call.
native_mainloop_wall_time/native_mainloop_cpu_time cover the native C loop.
MPI CPU sums participating ranks; wall takes their maximum, including the
maximum of each rank's setup+compute total rather than the sum of separate
phase maxima. CPU is never divided by thread or rank count. The CPU/wall
ratio estimates average CPU use, not elapsed-time speedup. Read timing_scope
for setup, output, cleanup and external-adapter boundaries.

--tree-theta defaults to 1 (approximate traversal) in angular drivers. Make THETA
is a distinct legacy compile-time criterion, also defaulting to 1. Opt-in
shear-pivot-reuse uses CBALLS_SHEAR_PIVOT_TOL=3 and CBALLS_SHEAR_BIN_THETA=0.5.
--nsmooth defaults to 8 in the active kappa and shear drivers, matching native cballs.
Scalar reuse uses phase tolerance 0.1 and bin tolerance 0. Environment values
override these compiled reuse controls, not theta. These are speed/accuracy
choices, not a universal percentage-error bound. Lyα keeps certified cell
aggregation with zero bin slops; positive slops remain explicit choices.
CBALLS_SHEAR_TASKS and CBALLS_SHEAR_TASK_MEMORY_MB also apply to MPI frontiers.
--verbose 1 enables available live work counters. CBALLS_SHEAR_PROFILE adds
instrumentation whenever present, including value 0; unset it to disable.

Shear retains full count/window modes by default. --more-options
no-count-multipoles selects N0 normalization; edge-corrections takes
precedence and restores higher window modes. Scalar dual-node runs already
skip higher count modes without edge correction; ggg-full-window is the
octree compatibility diagnostic override. NMultipolesON is a legacy build
switch, not a substitute for these runtime controls.

See [the approximation defaults and 5% L2 criterion](../../docs/APPROXIMATION_DEFAULTS.md) for fast defaults, measured penalties, exact references, and the limits of finite-multipole reconstruction.

Measured approximation penalties

Pass --accuracy-reference EXACT_OUTPUT_ROOT --accuracy-report-only to compare
against saved exact results for the identical input, bins and multipoles. The
reference root must contain ENGINE/qualification for each native engine.
Each native run writes ENGINE/qualification/approximation-penalty.json and
prints raw, normalized and corrected relative L2 percentages where available,
including each shear component. Without a reference the penalty is UNMEASURED.
The default comparison target is 5% (--qualification-rtol 0.05); it does not
change traversal controls. Report-only retains results above the target;
without it a failed requested numerical qualification raises an error.
Relative L2 uses no denominator floor: zero reference norms are undefined and
absolute errors are retained. A small pooled L2 does not bound every bin.


## Parameters used by each run

Native runs write `OUTDIR/ENGINE/parameters-usedvalues` before counting starts
and refresh it after MainLoop completes. The `run_stage` comment distinguishes
startup records from completed calculations; failed searches retain their startup
record. dual node saves its record after its calculation completes.
Native runs retain the usual `name = value` entries from cballs, updated with
the final effective values. Python arguments, angular units, catalogue sizes,
environment overrides, compiled reuse defaults and active reuse controls are
included as `%` comments. This also preserves the shear record after its
temporary native directory is removed. The recorded native `rootDir` describes
that directory; `infile=m0` denotes data supplied from Python, not a disk file.

MPI publishes this file only on global rank zero. dual node adapters save their
actual correlation-object bin/slop settings under the same filename, clearly
identified as dual node settings rather than native cballs parameters. Parameter
output is outside the reported correlation compute time. Explicit argument and
environment overrides are recorded; omitted values use the effective defaults.
