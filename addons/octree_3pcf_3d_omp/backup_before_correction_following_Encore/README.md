# `octree-3pcf-3d-omp`

cBalls addon for scalar-field 2PCF and 3PCF multipoles in three dimensions.

## Purpose

This mode keeps the existing projected/2D cBalls estimators unchanged and adds a
separate 3D estimator for catalogs with columns

```text
x  y  z  delta  weight
```

where `delta` is stored internally as `Kappa(p)` and `weight` as `Weight(p)`.
The estimator is a no-random, no-edge-correction scalar-field estimator. It is
appropriate for periodic boxes, homogeneous mocks, and algorithm validation. For
survey geometry/mask correction, a random-catalog estimator must be added later.

## Estimator

For each primary `i`, the code accumulates shell-binned spherical harmonic
coefficients

```text
a_lm^i(b) = sum_{j in shell b} w_j delta_j Y_lm^*(rhat_ij)
```

and then forms the Legendre multipoles by the addition theorem. Same-shell
self-pairs `j=k` are subtracted.

The output convention is normalized so that a constant scalar field gives
`zeta_0 = 1` and `zeta_l>0 ~ 0` for an isotropic catalog.

The optional 2PCF is accumulated in the same directed pivot-neighbor loop:

```text
xi(b) = sum_(i,j in b) w_i delta_i w_j delta_j
      / sum_(i,j in b) w_i w_j
```

The harmonic implementation uses a Cartesian recurrence for
`P_l^m(z) exp(i m phi)`. Normalization factors are prepared once per worker, so
the neighbor loop no longer calls `atan2`, trigonometric functions, or `lgamma`.

## Search method

The method string is

```text
searchMethod = octree-3pcf-3d-omp
```

The alias

```text
searchMethod = octree-ggg-3d-omp
```

is also accepted.

The code reuses the cTreeBalls octree traversal and OpenMP parallel loop pattern
from `addons/octree_ggg_omp`, but the angular estimator is exact-in-angle: cells
are used for pruning only, and accepted leaves are individual bodies.

`rangeN` is always the physical radial cutoff for this exact estimator. The
legacy `Rcut/theta` option may change the global cutoff used by other search
methods, but it does not reduce either body acceptance or cell pruning here;
therefore enabling it cannot remove pairs from the outer radial bins.

## ASCII input

For a plain 5-column ASCII catalog:

```text
infile       = catalog_xyz_delta_w.txt
infileformat = multi-columns-ascii
columns      = 1,2,3,4,5
options      = pos-and-convergence-weight,KKKCorrelation
```

Column order is `x y z delta weight`.

Existing cBalls ASCII formats are preserved. The old `columns-ascii-all` path can
also be used when the file has the native cBalls header and includes `weight` and
`mask`.

## FITS input

For FITS tables, use the existing cfitsio path with five columns:

```text
infile       = catalog_xyz_delta_w.fits
infileformat = fits
columns      = 1,2,3,4,5
options      = with-weight,KKKCorrelation
```

Column order is again `x y z delta weight`. Without `with-weight`, the FITS
reader keeps the previous behavior and sets all weights to 1.

FITS input may provide a sixth integer `LOS_ID` column. When any of
`exclude-same-los`, `exclude-los`, or `exclude-pivot-los` is selected, six
columns are required and pivot-neighbor pairs with equal IDs are skipped:

```text
columns = 1,2,3,4,5,6
options = with-weight,exclude-same-los,KKKCorrelation
```

The LOS rule removes `LOS_i=LOS_j` and `LOS_i=LOS_k` terms. It does not remove
`LOS_j=LOS_k` terms within the harmonic product. Inputs without a LOS column
receive unique IDs, so the exclusion is then a well-defined no-op.

## Computation modes

With no mode option the historical behavior is retained and only the 3PCF is
computed. The available mode options are:

```text
compute-2pcf-3d,compute-3pcf-3d  # compute both
only-2pcf-3d                     # 2PCF only
only-3pcf-3d                     # 3PCF only
```

## Main bin/multipole parameters

```text
rangeN        = 200.0
rminHist      = 0.0
sizeHistN     = 20
useLogHist    = false
mChebyshev    = 4      # interpreted as lmax in this 3D mode
numberThreads = 16
```

## Output

The output file is written as

```text
<histZetaMFileName>_3d.txt
```

with columns

```text
ell  bin1  bin2  r1_center  r2_center  zeta  numerator  denominator
```

When the 2PCF is enabled, it is written to

```text
<histXi2pcfFileName>_3d.txt
```

with columns

```text
bin  r_center  xi  numerator  denominator
```
