# `octree-3pcf-3d-omp`

cBalls addon for scalar-field 2PCF and 3PCF multipoles in three dimensions.

The independent MPI sibling is described in
[octree-3pcf-3d-mpi](../octree_3pcf_3d_mpi/README.md). Both use the same
measurement kernel and fixed, ordered pivot-block reduction.

## Purpose

This mode keeps the existing projected/2D cBalls estimators unchanged and adds a
separate 3D estimator for catalogs with columns

```text
x  y  z  delta  weight
```

where `delta` is stored internally as `Kappa(p)` and `weight` as `Weight(p)`.
The default estimator is a no-random scalar-field estimator. It is appropriate
for periodic boxes, homogeneous mocks, and algorithm validation. An opt-in
ENCORE-style point-process survey estimator is also available; it combines a
data and random catalog, measures the survey window, and applies the isotropic
3PCF edge-correction matrix.

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

## Survey estimator and edge correction

Enable the survey path with:

    infile       = data_xyz_w.txt,random_xyz_w.txt
    infileformat = multi-columns-ascii,multi-columns-ascii
    iCatalogs    = 1,2
    columns      = 1,2,3,4,5
    options      = pos-and-convergence-weight,survey-estimator-3d,
                   compute-2pcf-3d,compute-3pcf-3d

Both files use 'x y z value weight'. In survey mode 'value' is ignored and the
catalog is treated as a weighted point process. Weights must be finite and
non-negative. cBalls computes

    alpha = sum_i w_i^D / sum_j w_j^R

and constructs two owned work catalogs:

    N: D - alpha R
    R: alpha R

The two input catalogs are kept in one common Cartesian coordinate frame.
Separately translating or centering the data and random files would change the
cross terms and is therefore explicitly avoided.

For the 2PCF the output is the Landy-Szalay/ENCORE ratio

    xi = (D - alpha R)^2 / (alpha R)^2.

For the 3PCF, cBalls measures signed numerator multipoles 'N_l' and positive
random multipoles 'R_l', forms 'f_l=R_l/R_0', and solves independently for each
distinct radial-bin pair:

    sum_lpp M[l,lpp] zeta_lpp = N_l/R_0

    M[l,lpp] = sum_lp f_lp sqrt((2l+1)(2lp+1)(2lpp+1))
                Wigner3j(l,lp,lpp;0,0,0)^2 / (4 pi).

As in ENCORE, 'mChebyshev' is the measured maximum multipole. By default the
edge-corrected file emits only 'ell <= mChebyshev-1', since the highest measured
mode supplies the extra window-coupling order. The diagnostic option
'survey-keep-top-multipole' also emits that last, truncation-sensitive mode.

Survey mode requires exactly two distinct catalogs and is incompatible with
'read-mask', 'remove-mean', and the LOS-exclusion options. Encode the angular
and radial selection directly in the random catalog. Radial shell pairs use
'bin1 < bin2', matching ENCORE's 3PCF estimator.

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

Survey outputs are written separately, leaving legacy output names untouched:

    <histXi2pcfFileName>_3d_survey.txt
    <histZetaMFileName>_3d_survey.txt

The survey 3PCF file contains both the standard Legendre coefficient and the
ENCORE-basis coefficient, along with raw 'N_l', raw 'R_l', 'R_0', a pivot-based
conditioning diagnostic, and a validity flag. Empty-window or singular radial
pairs follow the project zero-denominator policy: corrected values are written
as zero with 'valid=0'.

## cTreeBalls versus ENCORE example

The runnable example

    examples/compare_octree_3pcf_3d_encore.py

generates weighted data and random catalogs in a non-periodic survey window,
runs cTreeBalls on the two catalogs, builds a native CPU ENCORE executable with
matching radial bins and multipole order, and compares:

- the raw `(D-alpha R)^2` and `(alpha R)^2` pair counts;
- the raw ENCORE-basis `N_l` and `R_l` multipoles;
- the edge-corrected 2PCF and 3PCF.

Run the in-memory Cython example with:

    python examples/compare_octree_3pcf_3d_encore.py \
        --backend cython --encore-source /path/to/encore

Use `--backend cli --cballs /path/to/cballs` to exercise the standalone C
executable instead. The script writes comparison CSV files, a JSON summary,
logs, and two plots. The companion notebook is

    examples/compare_octree_3pcf_3d_encore.ipynb

Set `ENCORE_SOURCE` before launching Jupyter if ENCORE is not at the default
path used by the example.

## Regression tests

Run:

    tests/make_tests/run_test_octree_3pcf_3d_omp

The launcher executes the legacy scalar-field oracle and a separate direct
survey oracle. The latter brute-forces '(D-alpha R)^2', '(D-alpha R)^3', random
multipoles, the Gaunt coupling matrix, and the corrected solution, then compares
one-thread and multi-thread outputs.

When the legacy 2PCF is enabled, it is written to

```text
<histXi2pcfFileName>_3d.txt
```

with columns

```text
bin  r_center  xi  numerator  denominator
```
