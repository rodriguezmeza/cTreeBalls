# Lyman-alpha forest OpenMP addon

This addon computes weighted Lyman-alpha forest correlations from individual
forest pixels. It uses the cTreeBalls octree for exact range rejection; accepted
leaves are always evaluated as individual pixels.

## Build and methods

Set `LYAFORESTOMPON=1` in `addons/Makefile_addons_settings`, then build normally.
The available search methods are:

- `lya-2pcf-omp`: weighted 2PCF only.
- `lya-3pcf-omp`: weighted anisotropic 3PCF only.
- `lya-2pcf-3pcf-omp`: both estimators in one tree traversal.
- `lya-1d-2pcf-omp`: weighted radial-only 2PCF.
- `lya-1d-tree-2pcf-omp`: exact interval-tree radial-only 2PCF.
- `lya-1d-3pcf-omp`: weighted radial-only 3PCF.
- `lya-1d-2pcf-3pcf-omp`: both radial estimators in one scan.

All methods require `DEFDIMENSION=3`, OpenMP, `usePeriodic=false`, exactly one
input file, and `infileformat=lya-ascii`.

## Input

The ASCII interchange format has six columns:

```text
x y z delta weight forest_id
```

Coordinates are observer-centered comoving Cartesian coordinates. `forest_id`
identifies the quasar sightline. Pairs and triplets containing two pixels from
the same forest are excluded.

Convert one or more lya2pcf `data*.npy` files with:

```bash
python3 addons/lya_forest_omp/lya2pcf_to_cballs.py \
  --lya2pcf-source /path/to/lya2pcf \
  -o forests.txt /path/to/data*.npy
```

The converter writes `delta = dw / we`, so the weighted delta stored by lya2pcf
is preserved exactly. Unpickling requires the Python dependencies used by the
lya2pcf checkout.

## Estimators

The 2PCF accumulates `sum(w_i delta_i w_j delta_j) / sum(w_i w_j)` in
`(r_parallel, r_transverse)` bins. Configure it with `lya2RpMax`, `lya2RtMax`,
`lya2RpBins`, and `lya2RtBins`.

The 3PCF accumulates `sum(w_i delta_i w_j delta_j w_k delta_k) /
sum(w_i w_j w_k)` for three distinct forests. Its five coordinates are the two
pivot-side lengths, their two angles to the pivot line of sight, and the cosine
of their opening angle. Configure it with `lya3RMax`, `lya3RBins`,
`lya3ThetaBins`, and `lya3MuBins`.

The radial methods build a sorted index of the observer-centered distance and
ignore the angular coordinates during candidate selection and binning. The
radial 2PCF uses `|chi_j-chi_i|` with `lya2RpMax` and `lya2RpBins`. The radial
3PCF uses the two signed pivot lags `(chi_j-chi_i, chi_k-chi_i)`. Each axis has
`2*lya3RBins` bins spanning `(-lya3RMax, lya3RMax)`, so configurations on the
same and opposite sides of a pivot remain distinguishable. Pixels in the same
quasar forest are excluded from both estimators; 3PCF triplets require three
distinct forest IDs. The build remains `DEFDIMENSION=3` because that describes
catalog storage, while the radial search itself is one-dimensional.

OpenMP work is grouped into fixed radial-pivot blocks. Every block processes
its pivots in radial order, and completed block histograms are merged in block
order. This makes the floating-point result independent of the OpenMP thread
count while avoiding a serialized merge after every pivot. Configure the
compile-time block size with `LYA1D_OMP_PIVOT_BLOCK_SIZE` in
`addons/Makefile_addons_settings`; the default is 256. Smaller values expose
more parallel tasks for small or strongly imbalanced catalogs, while larger
values reduce merge synchronization for inexpensive 2PCF workloads.

`lya-1d-tree-2pcf-omp` is an independent, exact alternative for dense radial
2PCF workloads. It builds a balanced binary interval tree over the same radial
ordering. A node pair is accumulated in bulk only when its minimum and maximum
possible separations lie in the same `lya2RpBins` bin; otherwise the wider node
is split. The tree first accumulates all unordered pairs, then subtracts pairs
within each `forest_id` using a sorted per-forest range scan. Long-double
internal sums and exact per-bin pair counts protect this subtraction against
false occupied bins. Tree tasks and histogram commits use fixed logical order,
so output is byte-identical across OpenMP thread counts. The ordinary
`lya-1d-2pcf-omp` range scan remains the preferred baseline for sparse catalogs
or very narrow radial windows.

Outputs are `<histXi2pcfFileName>_lya.txt` and
`<histZetaFileName>M_lya5d.txt`. The 3PCF output is sparse by default; add
`lya-output-empty-bins` to `options` to include zero-weight bins. Both
estimators define a zero-denominator bin to have correlation zero.

Radial outputs are `<histXi2pcfFileName>_lya1d.txt` and
`<histZetaFileName>M_lya1d.txt`. They use the same zero-denominator and sparse
3PCF policies.

## Tests

The regular regression uses independent NumPy 2PCF and 3PCF oracles and checks
OpenMP determinism:

```bash
make LYAFORESTOMPON=1 test-lya-forest-omp
```

The radial regression additionally checks a direct Python oracle, the scan and
tree implementations against one another, same-quasar subtraction, invariance
under arbitrary changes to transverse coordinates, and one-thread versus
multi-thread byte-for-byte repeatability:

```bash
make LYAFORESTOMPON=1 test-lya-forest-1d-omp
```

An additional end-to-end test runs the original lya2pcf CPU implementation on
real `forest_class.quasar` objects, converts the same saved dictionary, and
compares every cTreeBalls 2PCF numerator, denominator, and normalized bin:

```bash
make LYAFORESTOMPON=1 \
  LYA2PCF_SOURCE=/path/to/lya2pcf \
  test-lya2pcf-reference
```

This reference checkout implements only the 2PCF computation. The independent
oracle remains the regression reference for the addon 3PCF.

Compare the native radial estimators, the radial interval tree, and the 3D
octree estimators on the same exactly collinear catalog with:

```bash
python3 examples/compare_lya_1d_3d.py --threads 4
```

The script projects the 3D numerator and denominator onto the 1D bins before
normalization, checks the dedicated and combined methods, records wall times,
and writes CSV, JSON, and PNG results.  The companion notebook is
`examples/compare_lya_1d_3d.ipynb`.  Exact numerical agreement is meaningful
only for a collinear catalog; on a general survey the 1D and 3D methods apply
different transverse-selection rules.

Per-thread 3PCF storage is approximately
`lya3RBins^2 * lya3ThetaBins^2 * lya3MuBins *
(2*sizeof(real) + sizeof(size_t))` bytes. This
addon computes raw estimators; survey distortion matrices, covariance
estimation, continuum fitting, and metal corrections remain preprocessing or
post-processing responsibilities.
