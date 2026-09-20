# Full-sky shear all-engines driver

`tests/python/shear_corr_all_engines.py` runs one retained spin-2 catalog through
the three full-sky dual node addons enabled in the public Makefile profile:

- `octree-shear-sphere-2balls-omp`
- `kdtree-shear-sphere-2balls-omp`
- `balltree-shear-sphere-2balls-omp`

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
the resolved SLEEF state, actual vector-log backend, and active ball-tree
policies. The driver runs all three probes once and stores their return codes,
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
found in both the executable and extension. The driver also understands MPI
worker launch for separately configured builds, but no MPI shear addon is
enabled or shipped in this public profile.

Use separate `--statistics 2pcf` and `--statistics 3pcf` runs for performance
measurements. `both` exercises the combined production path and is useful for
cross-engine output checks.

For a first large-map test, add `--max-points 100000`. This deterministic
selection is useful for timing and regression checks but is not the full
survey estimator.

With `SMOOTHPIVOTON=1`, capable paths enable pivot smoothing by default;
`--no-smooth-pivot` disables it. A literal `--smooth-radius` is specified in
arcminutes and must obey `2*rsmooth <= min-sep`. Use `--exact-tree` together with
`--no-smooth-pivot` for exact body traversal. `--nsmooth` controls the native
leaf/smoothing capacity; `--more-options NAME[,NAME...]` adds native options.

Native octree and ball-tree 3PCF support opt-in `--more-options shear-pivot-reuse`
with `--no-smooth-pivot`, `BALLS4SCANLEVON=1`, positive `theta` and full pivot
coverage. KD-tree rejects this option. `CBALLS_SHEAR_PIVOT_TOL` sets a phase
budget in radians in `[0,3]`, default `0.1`, not a percentage error guarantee.
With reuse active the budget controls 3PCF acceptance instead of the magnitude
of `theta`. Zero budget only disables reuse, not ordinary neighbor acceptance.
Use `no-one-ball` for an exact native-octree 3PCF. Calibrate against exact
results for the same catalog before choosing a production budget.

FITS input must provide `GAMMA1` and `GAMMA2`, or selectors supplied with
`--gamma1-field` and `--gamma2-field`. `--mask` is applied before catalog
registration. NPZ input uses `positions`, `gamma1`, `gamma2`, optional
`weights`, and optional `geometry`.

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
