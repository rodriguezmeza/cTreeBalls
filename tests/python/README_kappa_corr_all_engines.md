# Convergence all-engines driver

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
acceptance; `--dual-node-bin-slop` enables the looser bin-aware policy, and
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
