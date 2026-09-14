# Full-sky shear all-engines driver

`shear_corr_all_engines.py` reads one spin-2 catalog and compares the three
full-sky dual-node add-ons enabled by the public build profile:

- `octree-shear-sphere-2balls-omp`
- `kdtree-shear-sphere-2balls-omp`
- `balltree-shear-sphere-2balls-omp`

All three use observer-centered unit vectors, local east/north shear
components, great-circle parallel transport, and the same cTreeBalls
2PCF/3PCF output contract.

## Runtime inventory

The executable is the authoritative source for compiled methods, recognized
options, numerical backends, and build settings:

```sh
./cballs options=print-search-methods
./cballs options=print-options
./cballs options=make-info
```

The driver runs these probes once and records their return codes, full output,
parsed method and option names, and Make settings in
`summary.json["ctreeballs_runtime"]`. Use `--list-engines` for a compact view
of the three shear methods.

## Example

```sh
python3 tests/python/shear_corr_all_engines.py \
  --fits run2/Takahasi/allskymaps_fits_shear/allskymap_nres12r081_zs12_mag_shear.fits \
  --geometry sphere --engine all --statistics both \
  --sep-units arcmin --min-sep 3 --max-sep 120 \
  --threads 16 --outdir Output_shear
```

For a first large-map test, add `--max-points 100000`. This deterministic
selection is useful for timing and regression checks but is not the full
survey estimator.

Use `--statistics 2pcf` and `--statistics 3pcf` in separate performance runs;
`both` exercises the combined production path. With `SMOOTHPIVOTON=1`, capable
methods enable pivot smoothing by default and `--no-smooth-pivot` disables it.
A literal `--smooth-radius` is specified in arcminutes and must satisfy
`2*rsmooth <= min-sep`. Use `--exact-tree` for body-level traversal and
`--more-options NAME[,NAME...]` for additional native controls.

FITS input must provide `GAMMA1` and `GAMMA2`, or selectors supplied with
`--gamma1-field` and `--gamma2-field`. The optional mask is applied before
catalog registration. NPZ input uses `positions`, `gamma1`, `gamma2`, optional
`weights`, and optional `geometry`.

## Outputs

Each method writes a compressed histogram file beneath `--outdir`.
`summary.json` contains numerical comparisons, the runtime inventory, and full
timing metadata. `timing_report.txt` lists setup, compute wall, process CPU,
and native-reported CPU time. The selected order option isolates the requested
statistic.

The driver plots `xi+`, `xi-`, corrected Gamma-x components, and flattened
3PCF multipole views. `--no-plots` disables figures. Set
`CBALLS_SHEAR_PROFILE=1` to include native per-thread transport geometry,
radial lookup, ring accumulation, scheduler, and reduction diagnostics.
