Production use of lya-anisotropic-multipole-3pcf-omp
=================================================

The regular lya_corr_all_engines.py driver now accepts this engine explicitly.
It is excluded from the automatic all/all-tree/all-multipole groups to keep
their existing estimator selection and computation cost. Its controls are:

  --lya3-lmax 8                   raw anisotropic Legendre order, 0..32
  --lya3-multipole-kernel 0        0 automatic hierarchy; 1 prefix; 2 forced groups
  --lya3-mu-mode exact            moments plus additional exact opening-angle bins
  --lya3-mu-mode finite-multipole  moments plus finite-L reconstruction only

The defaults are Lmax=8, moment kernel=0, and mu-mode=exact. This reproduces
the mode that passed the supplied calibration. No native rebuild is required;
these controls expose existing compiled parameters. --multipole-lmax is a
separate control for the scalar octree-3pcf-3d engines, not this engine.

From run2:

  python3 tests/python/lya_corr_all_engines.py \
      --fits catalogs/desi-lya_downloaded/combined_lya_0500.fits \
      --fits-layout desi --engine lya-anisotropic-multipole-3pcf-omp \
      --statistics 3pcf --threads 16 --max-forests 100000 --pixel-stride 8 \
      --r3-max 80 --r3-bins 6 --theta-bins 6 --mu-bins 6 \
      --lya3-lmax 8 --lya3-multipole-kernel 0 --lya3-mu-mode exact \
      --output Output_multipoles_production

Do not supply --lya3-kernel 5 for a multipole-only run. That control belongs
to hard-bin forest engines. The frontend rejects this combination before
reading the catalog. When running a hard-bin engine alongside this engine,
--lya3-kernel applies only to the hard-bin engine, and the separate
--lya3-multipole-kernel controls the anisotropic moments.

The native 3PCF geometry retains radial bins r1,r2 and both polar angles
theta1,theta2 relative to the pivot LOS. Here mu is the cosine of the opening
angle between the two triangle legs. It is distinct from the two-point
separation-to-LOS cosine used for DESI 2PCF wedges.

Outputs in ENGINE/ under the selected output root
------------------------------------------------

  histZetaM_lya_multipoles.txt
    Exact signed raw numerator and denominator Legendre moments, to floating
    precision. Columns b1 b2 t1 t2 ell numerator denominator. Higher denominator
    moments can be negative or zero. No N_ell/D_ell correlation is manufactured.

  histZetaM_lya5d.txt
    Additional exact angular-bin correlation/numerator/denominator, only with
    --lya3-mu-mode exact. It has the same 13-column contract as the exact
    3D forest estimators. It participates in ordinary exact-engine comparisons
    and radial projection plots. Empty bins follow the usual native convention.

  histZetaM_lya5d_multipole.txt
    APPROXIMATE finite-L reconstruction in the 13-column grid. Signed sums are
    retained, and denominator<=0 produces NaN correlation. This is never included
    in the exact-product comparison pipeline, even when exact bins also exist.

  anisotropic_raw_moments.png
    Raw signed numerator/denominator maps for selected ell (0,1,2,Lmax as
    available), summed over the two polar-angle bins, without forming ratios.

  anisotropic_finite_reconstruction.png
    Approximate reconstructed values and signed sums by native bin index.
    At most 6000 bins are drawn; native tables always retain all bins.

summary.json records engine parameters, native provenance, the three product
contracts, counts of undefined/negative-denominator reconstructed bins, and
the extra plot paths. Exact mode also reports reconstruction differences
against the exact companion on the same input. These diagnostics do not
qualify a different catalog or larger selection. --no-plots omits images.

The supplied calibration
------------------------

Output_multipoles_calibration/summary.json describes 1000 DESI forests,
73538 pixels, stride 8, Rmax=80, six radial/polar/mu bins, Lmax=8, kernel=0,
and mu-mode=exact. The exact companion passed the direct-reference test:
maximum eligible-bin relative correlation error approximately 3.15e-12.

The finite-L reconstruction DID NOT pass the 5% criterion. Its numerator
relative L2 error was about 40.7%, its 95th-percentile eligible-bin relative
correlation error about 568.6%, and its maximum about 57449.6%. Acceptance of
the exact companion is not acceptance of this approximate product.

Use exact mode for the binned estimator supported by this calibration. If
your intended observables are the raw moments themselves, finite-multipole
mode avoids the additional exact-bin computation; the reconstructed file
still remains an unaccepted diagnostic for this L=8 calibration. Exact mode
does additional work and is not a promise of moment-only execution speed.

To compare exact bins in one production invocation, explicitly select both:

  --engine lya-anisotropic-multipole-3pcf-omp lya-los-tree-3pcf-omp \
  --statistics 3pcf --lya3-mu-mode exact --lya3-kernel 5 --fail-on-mismatch

The native result-packet --qualification-reference workflow currently rejects
anisotropic moment packets as exact references, so this frontend rejects that
combination early. Use the calibration script for raw moment/reconstruction
studies or a hard-bin engine alongside this engine for exact-bin comparisons.
Zero engine comparisons is explicitly identified as no independent accuracy
certificate. Statistical uncertainty, model accuracy and observational
systematics are outside these numerical comparisons.

Multipole maps and radial slices (2026-10-16)
-------------------------------------------
New runs create multipole maps and radial slices automatically unless
--no-plots is supplied. All available orders are shown by default, paginated
at nine orders per figure. To select orders, add for example:

  --plot-multipole-orders 0 1 2 3 4

The scalar octree-3pcf-3d-omp/mpi engines also get all-order maps and radial
slices, preserving their native normalization and independent shell pairs.

To plot existing results WITHOUT loading FITS/NPZ or running the search:

  python3 tests/python/lya_corr_all_engines.py \
      --plot-existing Output_multipoles_production

This reads the saved summary.json and histogram files and writes to
Output_multipoles_production/multipole_plots/. Optionally choose another
destination with --output Output_multipole_plots. Existing native histograms
and the original summary.json are preserved; plot_summary.json records the
source files, formulas, support counts, orders and generated paths. Run this
command as a single Python process. No native rebuild is needed.

New anisotropic files under the engine subdirectory:

  anisotropic_multipoles_maps.png / anisotropic_multipoles_slices.png
    Q_l(r1,r2) = (2l+1) sum_(theta1,theta2) N_l /
                          sum_(theta1,theta2) D_0.
    These are exact sampled, window-weighted Legendre moments. The survey
    angular window still mixes orders, so Q_l is not a window-corrected
    correlation multipole. Never divide by a signed D_l. Q_0 is the same
    angular-summed numerator/denominator correlation as the standard radial
    image. Empty D_0 cells remain NaN (gray in maps, gaps in lines).

  anisotropic_window_multipoles_maps.png / ..._slices.png
    W_l(r1,r2) = (2l+1) sum_(theta1,theta2) D_l / sum_(theta1,theta2) D_0.
    These show angular sampling and window structure. W_0=1 in occupied cells;
    higher moments may be signed. They are not a correction applied to Q_l.

  anisotropic_multipoles.csv
    All computed orders, regardless of figure selection. Columns:
    b1,b2,ell,r1_mpc_h,r2_mpc_h,N_ell,D_ell,D_0,Q_ell,window_ell.

Only when exact angular bins are available:

  anisotropic_binned_multipoles_maps.png / ..._slices.png
    First pool N_mu and D_mu over the polar-angle cells, form zeta_mu=N_mu/D_mu,
    then integrate its piecewise-constant Legendre projection:

      zeta_l^binned = (2l+1)/2 sum_mu zeta_mu integral_(bin) P_l(mu) dmu.

    The Legendre polynomial is integrated over each bin, not evaluated only
    at its center. This is a projection of the exact-bin input, with finite
    mu-bin resolution. It is distinct from the native finite-L reconstruction.
    All mu bins in a radial cell must have positive pooled weight; otherwise
    every coefficient in that cell remains NaN. The angular window is not
    deconvolved within each mu bin or within the retained polar bins.
    Orders are limited to min(lya3_lmax,mu_bins-1). This limit does not certify
    convergence: increase mu resolution and compare if high orders matter.
    For six saved mu bins, these figures contain at most l=0..5; Q_l and W_l
    still contain all saved raw orders, including l=8.

  anisotropic_binned_multipoles.csv
    b1,b2,ell,r1_mpc_h,r2_mpc_h,zeta_ell_binned,mu_coverage.
    Coverage is the supported fraction of [-1,1], without filling missing bins.

Each slice panel shows up to three occupied fixed-r1 rows (first, middle,
last), plus the r1=r2 diagonal when it is populated. There are no inferred
statistical error bars; these tables do not supply a multipole covariance.
