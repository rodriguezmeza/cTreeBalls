Weighted continuum projection for cTreeBalls Ly-alpha estimators
==============================================================

For every forest q, using the final correlation weights w and the true
observed log wavelength L=log(lambda), define

  x_i = L_i - sum(w*L)/sum(w)
  delta_projected_i = delta_i - sum(w*delta)/sum(w)
                     - x_i * sum(w*x*delta)/sum(w*x*x).

This is the projector eta of arXiv:2507.00129 equation (2.7), and agrees with
quasar.fill_dw in the user's lya2pcf_2026-03-14/forest_class.py. The reference
returns w*delta_projected; cTreeBalls receives delta_projected and w separately.
Every native engine applies the weights in its usual estimator, exactly once.

Apply the correction once before building/searching trees. All engines in one
driver invocation reuse the same corrected catalog: octree, LOS-tree,
quadtree, radial, same-LOS, multipole, OMP and MPI variants. No binning or
forest-exclusion rule is changed. The correction affects both the 2PCF and
3PCF. It is not an inverse distortion-matrix correction.

Python driver
-------------
Append this flag to the existing run command:

  --project-delta

The FITS path already supported mean+slope removal; it now uses the shared
projector and records per-forest constraint checks. DESI, eBOSS/PICCA and
RA-DEC-R formats supply real wavelengths. Cartesian-only inputs without
wavelength metadata cannot define this slope and are rejected when projection
is requested. Distance is never substituted for log wavelength.

NPZ caches containing positions, delta, weights, forest_ids and log_lambda now
also support the same flag, for example:

  python3 tests/python/lya_corr_all_engines.py \
      --catalog Output_multipoles_calibration/catalog.npz --project-delta \
      --engine lya-los-tree-2pcf-3pcf-omp --statistics both --threads 16 \
      --r3-max 80 --r3-bins 6 --theta-bins 6 --mu-bins 6 \
      --lya3-kernel 5 --output Output_lya_projected

Projection is explicitly requested; the default preserves the supplied delta
field. Existing calibrated runs therefore retain their current preprocessing.
If a cache is already projected with identical pixels and weights, applying
the projector again is mathematically idempotent (up to roundoff), but it is
not necessary. Changing selection or weights defines a different projector.

The existing FITS convention is preserved: projection fits the retained
pixels after cuts/stride. To match an upstream fit on a full forest, project
with all its valid pixels before later thinning. This scope is recorded in
the catalog metadata. The optional redshift-weight factor is applied before
projection. To additionally reproduce the reference reader's weighting, use
--redshift-weight-exponent 1.9 --weight-z-ref 2.25 if appropriate to the input;
the projection flag alone does not introduce that factor.

Direct cyballs API
------------------
The rebuilt extension supports:

  balls.set_forest_catalog(positions, delta, weights, forest_ids,
                           project_delta=True, log_lambda=log_lambda)

It projects before registering the catalog, so subsequent searches use the
corrected field. The four-argument call keeps the previous behavior. Exact
integer identities are retained, including IDs larger than 2**53 and forests
split into nonadjacent chunks. Original NumPy inputs are not modified.

A reusable projector is also exported:

  from cyballs import project_lya_delta
  delta_p, diagnostics = project_lya_delta(delta, weights, forest_ids,
      log_lambda, return_diagnostics=True)

Register delta_p with project_delta=False (the default), then run any engine.
The compiled API and FITS/NPZ helpers share the same source implementation.
The native six-column lya-ascii format has no wavelength field; prepare its
delta column beforehand or use the wavelength-aware Python interface.

Numerical behavior and provenance
---------------------------------
Constant-wavelength or single-positive-pixel forests use mean removal only.
Zero-weight pixels stay unchanged and cannot enter a weighted estimator.
Invalid shapes, nonfinite values, negative weights and floating-point forest
IDs are rejected. Coordinate/weight scaling stabilizes weighted sums without
altering the projector. No dense per-forest projection matrix is allocated.

summary.json and --save-catalog archives retain project_delta, fitting-domain
scope, and projection_constraints. The latter records the maximum absolute
weighted mean and normalized first-moment residual across forests. The first
moment is normalized by sqrt(sum(w)*sum(w*x*x)), in delta units. These are
numerical constraints, not statistical uncertainties or an accuracy test for
the scientific model.

An NPZ explicitly reprojected on its retained pixels drops any earlier full
fitting-domain mapping. Computing a 2PCF distortion matrix then uses this same
retained domain and the same weights; projection order must be 1. The existing
distortion matrix models the 2PCF response only, not a full 3PCF response.
