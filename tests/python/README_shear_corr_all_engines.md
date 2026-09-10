# Full-sky shear all-engines driver

`shear_corr_all_engines.py` reads one spin-2 catalog and compares the active
full-sky dual-node implementations:

- `octree-shear-sphere-2balls-omp`
- `kdtree-shear-sphere-2balls-omp`
- `balltree-shear-sphere-2balls-omp`

All three use observer-centered unit vectors, local east/north shear
components, great-circle parallel transport, and the same 2PCF/3PCF output
contract.

## Example

```sh
python3 tests/python/shear_corr_all_engines.py \
  --fits run2/Takahasi/allskymaps_fits_shear/allskymap_nres12r081_zs12_mag_shear.fits \
  --geometry sphere --engine all --statistics both \
  --sep-units arcmin --min-sep 3 --max-sep 120 \
  --threads 16 --output Output_shear
```

For a first large-map test, add `--max-points 100000`. This deterministic
selection is useful for timing and regression checks but is not the full
survey estimator.

Use `--statistics 2pcf` for the dedicated pair kernels or `--statistics 3pcf`
to skip pair work. With `SMOOTHPIVOTON=1`, capable compatibility paths enable
pivot smoothing by default; `--option no-smooth-pivot` disables it. A literal
`--smooth-radius` is specified in arcminutes and must obey
`2*rsmooth <= min-sep`.

FITS input must provide `GAMMA1` and `GAMMA2` fields, or selectors supplied by
`--g1-field` and `--g2-field`. `--mask` is applied before catalog registration.
NPZ input uses `positions`, `gamma1`, `gamma2`, and optional `weights`; use
`--geometry sphere` for three-dimensional observer-centered vectors.

Each engine writes its native files beneath `--output`. `summary.json` and
`timing_report.txt` contain setup, compute wall time, and process CPU time.
The driver plots `xi+`, `xi-`, natural 3PCF components, and flattened 3PCF
multipole views. `--no-plots` disables figures.
