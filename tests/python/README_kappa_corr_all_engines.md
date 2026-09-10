# Convergence all-engines driver

`kappa_corr_all_engines.py` reads a convergence catalog once, registers its
NumPy arrays with `cyballs`, and runs the active scalar angular engines. The
current profile exposes:

- `kdtree-2balls-omp` and `kdtree-2balls-mpi`
- `balltree-2balls-omp` and `balltree-2balls-mpi`
- `octree-2balls-omp` and `octree-2balls-mpi`

These are native cTreeBalls engines. Their production traversal uses dual-node
acceptance; `--dual-node-bin-slop` enables the looser bin-aware policy, and
`--option no-two-balls` requests an exact body-level validation walk where
supported.

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
  --catalog catalog.npz --engine octree-2balls-mpi \
  --mpi-ranks 2 --threads 4 --outdir Output_kappa_mpi
```

Use repeated `--mpi-extra-arg` options when the local launcher needs them.
`mpi4py`, `mpiexec`, and the extension must use the same MPI implementation.

## Input and results

`--fits` accepts an implicit HEALPix scalar map. `--catalog` accepts an NPZ
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
