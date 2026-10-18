RA-DEC-R inputs and comparison with DESI DR1
==========================================

The driver supports an explicit distance convention:

  --radec-r-distance-source input     use supplied rcdis (default)
  --radec-r-distance-source redshift  use observer distance chi(zobs)

The input mode requires rcdis to be the absolute radial distance from the
observer, in --radec-r-distance-unit (default Mpc/h) unless TUNIT is present.
A radial offset is not a Cartesian translation: it changes transverse pair
separations. Relative distances within a box or redshift slice cannot simply
be multiplied by the sky unit vector for this estimator.

The redshift mode uses FlatLambdaCDM with --omega-m (default .315), --h
(default .674), and Tcmb0=0, matching the existing DESI reader's fiducial
distance convention. The resulting positions are in Mpc/h. It uses zobs,
checks consistency with wobs/1215.67-1, and does not use the rcdis values to
place pixels. This is an explicit analysis choice, not a determination of
the original catalog's cosmology or coordinate convention.

--project-delta removes each retained forest's weighted mean and its slope
against log wavelength, using the same final pixel weights as the estimator.
It operates after cuts/stride, on the joined tid forest. The original FITS
is never modified. An upstream unweighted zero mean does not imply zero
mean with the wdelta weights used by this estimator.

For the supplied Rodrigo catalog, a diagnostic comparison can be run with:

  python3 tests/python/lya_corr_all_engines.py \
      --fits Rodrigo/deltas-NoDLAm-NoNoise-True-dz02_LYA.fits \
      --fits-layout radec-r --radec-r-distance-source redshift \
      --project-delta --engine lya-los-tree-2pcf-omp \
      --threads 16 --max-forests 10000 --pixel-stride 2 \
      --compare-desi-dr1 --lya2-kernel 1 \
      --output Output_Rodrigo_redshift_projected

--outdir is accepted as an alias for --output. The driver prints the retained
weighted and unweighted delta means and stores field_diagnostics in the
catalog metadata (and any --save-catalog NPZ). The distance source and
fiducial cosmology are also recorded. Missing rcdis TUNIT is identified as
an assumption in the per-extension metadata.

Findings for the original 10,000-forest, stride-2 run (2,256,093 pixels):
  weighted mean delta             -0.031042782371839946
  unweighted mean delta           -1.5559593011595908e-6
  squared weighted mean            9.636543373854169e-4
  pair-weighted xi at 160<r<200     8.428602664977159e-4

A nonzero weighted mean can add a disconnected contribution to the product
estimator. The global mean squared is only a magnitude estimate, because
the contributing pixels/weights differ between separation bins. At r=200
it would contribute about 38.55 to r^2 xi, close to the large rise in the
original plot. This is not a prescription to subtract a constant from xi.

An example stored point has zobs=2.185036936 and rcdis=776.254. The fiducial
observer distance is 5572.443 Mpc, or 3755.827 Mpc/h. Neither interpretation
of the unlabelled rcdis as Mpc or as Mpc/h agrees with that observer distance.
The catalog producer's coordinate and weight definitions remain unknown.

A controlled 512-forest run at the original supplied positions found:
  original delta:     weighted mean -0.0294633; xi(120<r<200) = 7.93204e-4
  projected delta:    weighted mean 3.84e-18; xi(120<r<200) = -1.53421e-6
Both engines (LOS-tree and r_parallel/r_transverse quadtree) agree on the
unprojected 512-forest correlation to maximum absolute difference 6.94e-18.

DESI uses projected fluctuations and corresponding forward-model distortion
matrices (Adame et al., section 3.1, equations 7-9; section 4.3):
https://arxiv.org/html/2404.03001v3#S3.SS1

--compute-distortion computes a forward-model response. It does not remove
an input mean or perform the pixel projection; use --project-delta when
that is the desired field definition. Changing weights or projection can
change the estimator even for a true-continuum mock. This particular mock's
upstream chunk normalization, selections, redshift distribution and weight
definitions are not established to match DESI DR1. These options therefore
do not guarantee agreement with Figure 4 or a statistically valid goodness
of fit. The displayed published error bars are not uncertainties on this
mock/subset. No original output or original data file is overwritten.
