Compare Lyman-alpha correlations from different catalogs
======================================================

compare_lya_catalogs.py reads the saved native histograms produced by
lya_corr_all_engines.py. It does not reread FITS or recompute the correlations.
Requirements: Python 3.10+, NumPy, and Matplotlib (plots only). cyballs, MPI,
Astropy and the original input catalogs are not needed for this comparison.

From run2, compare two runs made with the same engine and bin settings:

  python3 tests/python/compare_lya_catalogs.py Output_DESI Output_Rodrigo \
      --engine lya-2pcf-3pcf-omp \
      --labels DESI Rodrigo --reference DESI \
      --output Output_DESI_vs_Rodrigo

Add further input directories and labels to compare three or more catalogs.
--reference is a LABEL, and defaults to the first input. --statistics both
(the default) compares all saved product families shared by every input.
Use --statistics 2pcf or --statistics 3pcf to restrict the order.

For example, first compute both orders for the two FITS catalogs:

  python3 tests/python/lya_corr_all_engines.py \
      --fits catalogs/desi-lya_downloaded/combined_lya_0500.fits \
      --fits-layout desi --engine lya-2pcf-3pcf-omp --statistics both \
      --threads 16 --output Output_DESI

  python3 tests/python/lya_corr_all_engines.py \
      --fits Rodrigo/deltas-NoDLAm-NoNoise-True-dz02_LYA.fits \
      --fits-layout radec-r --engine lya-2pcf-3pcf-omp --statistics both \
      --threads 16 --output Output_Rodrigo

Use the same --rp-max/--rt-max/--r3-max, --rp-bins/--rt-bins/--r3-bins,
--theta-bins and --mu-bins in both computations. The examples use the driver's
defaults. Their catalog size, projection, blinding, noise, redshift, weights
and cosmology still determine the interpretation of the comparison.

If the runs used different engines, pass the engine directories themselves:

  python3 tests/python/compare_lya_catalogs.py \
      Output_DESI/lya-rprt-quadtree-2pcf-omp \
      Output_Rodrigo/lya-los-tree-2pcf-omp \
      --labels DESI Rodrigo --statistics 2pcf --output Output_comparison

Native histogram FILES are also accepted. Keep summary.json with its original
run when possible: it provides physical grid and catalog metadata. A root
containing multiple engines requires --engine or explicit engine directories.

Output
------

For each shared product family:

  FAMILY_bins.png       correlation overlay, absolute differences, percent
  FAMILY_maps.png       2D correlation/difference/percentage maps, when applicable
  FAMILY_radial.png     pooled radial curves and differences, when applicable
  FAMILY__*.csv         all native bins for each candidate/reference pair
  FAMILY_map__*.csv     the common-support map values, when applicable
  FAMILY_radial__*.csv  the common-support radial curves, when applicable
  comparison_summary.json  inputs, selections, grid notes, coverage, and metrics

--formats png pdf produces both image formats; svg is also supported.
--no-plots writes only CSV and JSON and does not import Matplotlib.
Percent panels use a symmetric log scale by default: it preserves signs and
zero and keeps small differences visible alongside large outliers. The linear
region is +/-1 percentage point (--relative-linthresh). Use
--relative-plot-scale linear for a fully linear percentage axis/colorbar.
This only changes the plot scale; data and reported metrics are unchanged.
Native-bin curves are displayed in lexicographic index order (or physical
separation for a 1D product). CSV includes every bin index and coordinate.
For large tables --max-plot-bins (default 6000) samples the native-bin PLOT
only; all data remain in the CSV and metrics. A very large/nonuniform detached
2D grid may omit its image; native-bin curves and CSV remain available.

Relative differences and missing bins
------------------------------------

  difference         = candidate - reference
  relative_difference = difference / reference
  relative_percent   = 100 * relative_difference

The denominator is the SIGNED reference correlation. A negative reference
therefore changes the sign relative to the absolute difference. Native bins
are matched by indices and checked against physical coordinates and available
grid metadata. Different grids are rejected; they are never interpolated.
Different estimator families are compared separately: anisotropic 3D,
radial-only, same-LOS, and scalar multipoles are not interchangeable.

Absent rows and zero-denominator bins are undefined, exported as NaN and
excluded from differences; they are not replaced by zero correlation.
Coverage metrics separately count bins occupied only in one catalog.

Percent differences are undefined where |reference| <= --relative-floor
(default 1e-12). --relative-scale-floor 0.001 additionally excludes reference
values below 0.1% of the largest absolute reference in that view. Absolute
differences remain available. Thresholds are recomputed separately for native
bins, maps, and radial projections. The JSON records their effective values.
The reported relative L2 percentage uses ALL common occupied bins, without a
per-bin floor: 100 * ||candidate-reference||_2 / ||reference||_2. It is undefined
for zero reference norm and does not bound the largest individual difference.

Projection and slices
---------------------

For projected maps and radial curves, native bins must be occupied in EVERY
catalog in the comparison. Each catalog then pools its own raw numerator and
denominator on that same support. Correlations are never averaged without
their weights. Catalogs with different pair/triplet weights can consequently
have different projected values even if some native correlations agree.

2PCF maps retain (r_parallel, r_perp). Radial curves pool over r_perp.
Anisotropic 3PCF maps pool the angular dimensions into (r1, r2), and radial
curves additionally pool r2. Radial-only 3PCF maps retain signed lags.
Scalar multipole maps/radial curves use ell=0; native CSV retains all ell.
If the all-catalog intersection is empty, projections are omitted and native
pairwise comparisons are still exported.

To inspect a particular 3PCF angular slice instead of pooling every angle:

  python3 tests/python/compare_lya_catalogs.py Output_DESI Output_Rodrigo \
      --engine lya-2pcf-3pcf-omp --labels DESI Rodrigo \
      --product 3d_3pcf --select t1=1 --select t2=2 --select bmu=0 \
      --output Output_3pcf_slice

Indices are those printed in the histogram (normally zero based for forest
engines). --product can select 3d_2pcf, 3d_3pcf, 1d_2pcf, 1d_3pcf,
1d_same_los_2pcf, or multipole_3pcf. A selection applies before projection
and comparisons and must be valid for every selected product.

Scientific interpretation
-------------------------

These are relative DIFFERENCES between catalog estimates, not estimated
statistical errors or detection significances. No Poisson error is inferred
from pair/triplet counts. Assessing statistical consistency requires suitable
covariance, including cross-covariance for overlapping samples. This script
does not infer that covariance, apply a distortion matrix, unblind data, or
remove differences in the input selections. It compares the saved estimators
as supplied. The JSON retains available catalog metadata for review.
