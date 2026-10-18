Lyman-alpha Forests
===================

The forest addons use pixel positions, flux contrast ``delta``, statistical
weights, and integer forest/quasar IDs. Standard estimators exclude same-forest
pairs and require three different forests for triplets. A separately named
same-LOS estimator deliberately measures within-forest pairs. The dedicated
addons produce anisotropic or radial estimators; the all-engines driver can
additionally run the active scalar Legendre-multipole octree paths.

Choose a Method
---------------

Enable ``LYAFORESTOMPON=1`` and/or ``LYAFORESTMPION=1``. In each name below,
replace ``omp`` with ``mpi`` for its distributed sibling:

* ``lya-2pcf-omp``, ``lya-3pcf-omp``, ``lya-2pcf-3pcf-omp``:
  anisotropic physical-3D estimators.
* ``lya-1d-2pcf-omp``, ``lya-1d-3pcf-omp``, ``lya-1d-2pcf-3pcf-omp``:
  signed radial-lag estimators ignoring transverse distance.
* ``lya-1d-tree-2pcf-omp`` and ``lya-1d-tree-3pcf-omp``: exact interval-tree
  radial estimators with same-forest subtraction.
* ``lya-1d-tree-same-los-2pcf-omp``: one exact interval tree per forest;
  within-forest radial 2PCFs are normalized separately and averaged with equal
  LOS weight in each occupied bin. This method currently has no MPI sibling.
* ``lya-los-tree-2pcf-omp``, ``lya-los-tree-3pcf-omp`` and
  ``lya-los-tree-2pcf-3pcf-omp``: exact anisotropic 3D statistics with an octree
  for forest discovery and per-LOS radial trees for neighbor ranges. Unlike
  radial-only methods, transverse distance matters. MPI siblings use fixed
  pivot blocks, replicated catalogs and collective raw-sum reductions.
* ``octree-3pcf-3d-omp`` and ``octree-3pcf-3d-mpi``: scalar Legendre 3PCF
  multipoles with ``exclude-all-same-los`` enforcing three distinct forests.

All modes still require ``DEFDIMENSION=3`` and ``usePeriodic=false``.
"1D" describes the search coordinate, not the body-storage dimension.
Smooth-pivot is unsupported. Radial and 3D outputs measure different statistics
and should not be treated as interchangeable reference values.
``LYA1D_OMP_PIVOT_BLOCK_SIZE`` controls deterministic radial work blocks and
``LYA1D_TREE3_LEAF_SIZE`` controls the exact 3PCF interval-tree leaf capacity;
both are reported by ``options=make-info``.

Inputs and Bins
---------------

``infileformat=lya-ascii`` uses six columns ``x y z delta weight forest_id``.
Coordinates are observer-centered comoving distances. Forest IDs must remain
integers, including IDs larger than ``2**53``. The in-memory API is
``set_forest_catalog(positions, delta, weights, forest_ids)``; generic scalar
registration alone does not provide forest metadata.

The 3D pair grid uses ``lya2RpMax``, ``lya2RtMax``, ``lya2RpBins``, and
``lya2RtBins``. The 3D triplet grid also uses ``lya3RMax``, ``lya3RBins``,
``lya3ThetaBins``, and ``lya3MuBins``. Radial pairs use only the parallel
range/bins; radial triples have two signed-lag axes with ``2*lya3RBins`` bins
per axis. Weighted numerators and denominators are reduced before division;
zero-denominator bins return zero.

The same-LOS method is a distinct statistical family. Its output is
``histXi2pcf_lya1d_same_los.txt`` with columns ``bin radial_separation xi
sum_xi contributing_los``. Thus the published value is the equal-LOS mean
``sum_xi / contributing_los`` in each bin, rather than the pair-weighted
cross-forest correlation produced by the other radial 2PCF methods.

Examples
--------

From a checkout with the MPI addon enabled::

   mpiexec -n 2 ./cballs addons/lya_forest_mpi/parameters.ini

The all-engines driver reads DESI or eBOSS/PICCA delta FITS, NPZ, or six-column ASCII once,
broadcasts arrays once for MPI, and retains registered catalogs between engines::

   python3 tests/python/lya_corr_all_engines.py --list-engines
   bash examples/download_desi_lya_example.sh /tmp/desi-lya
   python3 tests/python/lya_corr_all_engines.py \
       --fits /tmp/desi-lya/delta-1019.fits.gz --max-forests 6 \
       --pixel-stride 30 --engine all-omp --threads 2

Use ``--engine all-tree --statistics both`` to select tree-based forest engines,
including radial interval trees, the same-LOS estimator and anisotropic LOS-tree
methods. Comparisons remain separated by estimator family.
Use ``--engine all-multipole --statistics 3pcf`` for the octree OpenMP/MPI
estimators. Multipole and five-dimensional products are reported as
different families and are never compared bin by bin.

The default driver statistic is 2PCF. Add ``--statistics both`` for 3PCF
and ``--engine all --mpi-ranks 2`` for an MPI-inclusive suite. Subsampling,
distance cosmology, and DESI blinding metadata matter scientifically; these
commands are examples, not a production DESI analysis.

See :download:`driver README <../tests/python/README_lya_corr_all_engines.md>`,
:download:`MPI addon README <../addons/lya_forest_mpi/README.md>`,
and ``examples/lya_corr_all_engines.ipynb``.
The script/notebook ``examples/compare_lya_1d_3d`` compares radial and 3D
workflows while keeping their estimator definitions distinct.

RA-DEC-R vector forest FITS
---------------------------

The driver also auto-detects ``SPECTRA_*`` tables with per-row sky coordinates
``ra, dec``, target identifier ``tid``, and pixel vectors
``rcdis, delta, wdelta, wobs, zobs``. Select ``--fits-layout radec-r`` to
require this format. It uses the supplied radius, with
``x = R cos(dec) cos(ra)``, ``y = R cos(dec) sin(ra)``, ``z = R sin(dec)``.
``delta`` is the field and ``wdelta`` the statistical weight. NaN padding and
invalid/nonpositive radii or weights are removed using one joint pixel mask.

FITS column units take precedence. For unlabelled columns, the fallbacks are
degrees and Mpc/h; select ``--radec-r-angle-unit rad`` or
``--radec-r-distance-unit Mpc`` explicitly when needed. Mpc radii are
multiplied by ``h`` (``--h``, default 0.674), giving native Mpc/h coordinates.
The supplied Rodrigo file has no unit labels, so confirm its radial units;
the loader preserves ``rcdis`` by default and records this assumption. It
does not replace ``rcdis`` using a cosmological distance-redshift relation.

Rows with the same ``tid`` across redshift extensions are joined into one
forest; ``cid`` identifies a chunk and does not create a new sightline.
``--max-forests`` counts complete target IDs. Redshift cuts use ``zobs``;
stride is applied after joining and wavelength sorting. True ``wobs`` values
provide log wavelength for optional projection, redshift weighting and the
continuum response. The response can retain all valid chunks even when
correlation cuts/stride retain only part of a forest. Floating ``tid`` exports
must be safely representable integer values; native integer IDs retain 64 bits.

For example, from ``run2`` (using the Mpc/h fallback)::

   python3 tests/python/lya_corr_all_engines.py \
       --fits Rodrigo/deltas-NoDLAm-NoNoise-True-dz02_LYA.fits \
       --fits-layout radec-r --radec-r-distance-unit Mpc/h \
       --max-forests 100 --pixel-stride 1 \
       --engine lya-los-tree-2pcf-omp --lya2-kernel 1 --threads 4

Live counting progress
----------------------

Native Ly-alpha methods report live progress with ``verbose=1`` (terminal) or
``verbose_log=1`` (``cballs.log``). The all-engines Python driver enables terminal
progress by default; use ``--verbose 0`` to silence native output or
``--verbose-log 1`` to record progress in the native log. No profiling variable
is required.

Each line identifies the engine, counting phase, rank, completed/total tasks,
percentage, and elapsed wall seconds. Tasks are coarse pivot blocks, forest-pair
rows, cell blocks, or radial tree/correction blocks, as indicated by the phase.
They have unequal costs: the percentage is completed work units, not a runtime
estimate or a count of individual pairs/triples. Updates are flushed about every
one percent, or after five seconds at the next completed task. A long individual
task can therefore leave the display unchanged for longer than five seconds.

Combined estimators may show separate pair and triple phases. Reaching 100%
finishes that counting phase; later phases, MPI reductions and file output may
still be pending. MPI follows the existing rank-zero logging policy and shows
only rank zero's assigned tasks, not a global fraction. Logging does not change
the histogram arithmetic, task order, or scientific approximations.

Validation
----------

The driver computes all pairwise comparisons within each compatible estimator
family, plots correlation functions and relative errors, and can fail on
specified tolerances. ``--fits-layout eboss`` selects eBOSS/PICCA delta HDUs;
these are prepared forest deltas, not raw spectrograph exposures. Optional
``--lya2pcf-source`` selects an external CPU 2PCF reference implementation.
Distortion matrices act on a compatible model vector, not on the measured
correlation as an inverse correction. Wedges, covariance and distortion
analysis options are documented in the driver README; no 3PCF external
reference is implied. Native and external timing scopes are recorded separately.

Run ``make test-lya-forest-omp test-lya-forest-1d-omp`` and
``make test-lya-forest-mpi`` for enabled profiles.
The MPI implementation replicates input/tree memory and partitions pivot or
task blocks. Thread-count changes are deterministic at fixed rank count;
changing ranks may change rounding. Native and Cython tests cover independent
reference calculations, scan/tree agreement, exclusions, strict in-memory multipoles, and
recoverable failures.

LOS Hierarchical Reuse
----------------------

Select ``lya2Kernel=1 lya3Kernel=5`` (Python: ``--lya2-kernel 1
--lya3-kernel 5``) for certified pair-range and adaptive triple moments.
The LOS sparse fallback retains forest discovery and builds a hierarchy
above per-forest radial/polar segments. Compact fixed pivot blocks also
share a conservative discovery frontier, including with the default pair
kernel. Every pivot retains its own cuts and forest exclusions.

Zero slops preserve the hard-bin estimator to summation rounding. These
methods do not replace the five-dimensional histogram by truncated
multipoles. MPI counterparts support the same kernels. See the
:download:`LOS implementation and validation guide <LYA_LOS_HIERARCHICAL_REUSE.md>`
for definitions, work counters, benchmark commands and limits.

Anisotropic Legendre moment acceleration
----------------------------------------

``lya-anisotropic-multipole-3pcf-omp`` now reuses forest-block harmonic
moments with an occupancy-based automatic choice at ``lya3Kernel=0``.
For this engine, kernel 1 preserves the prefix reference and kernel 2 forces
four-forest groups. Raw moments retain exact radial/polar bins; finite-order
reconstruction of hard opening-angle bins remains approximate. See
``docs/LYA_MULTIPOLE_HIERARCHICAL_REUSE.md`` and the addon README.


Exact angular-bin companion
---------------------------

For the anisotropic multipole OpenMP engine, ``lya3MuMode=1`` additionally
computes exact hard opening-angle bins using shared neighbor discovery and
certified hierarchy products with pixel fallback. Their ``_lya5d`` output is
independent of Lmax; ``_lya5d_multipole`` remains approximate. Mode 0 preserves
the existing behavior. The exact companion requires extra memory and runtime.
See :download:`the reconstruction guide <LYA_MULTIPOLE_RECONSTRUCTION.md>`.

Spatial index over complete forests
-----------------------------------

For ``lya2Kernel=1``, the pair search first uses a balanced kdtree over whole
forest bounds. Each node encloses all of its forests' actual pixel positions,
observer directions and radial extents. Conservative distance and parallel/
transverse separation bounds reject groups of distant forests together.
Median splits use bounding-box centres only for partitioning; a forest centre
is never used as a substitute for its spatial extent.

Retained forest pairs use the existing radial-segment hierarchy, moment sums,
pixel fallbacks and exact binning. Candidate bitmaps restore the previous
forest-row order, preserving summation order without sorting every candidate
list. The index requires at least 32 forests and a strict double build;
smaller or unsupported profiles retain the flat forest-row traversal.
OpenMP and MPI pair ownership are unchanged. Combined engines use this index
for their pair pass. The separate distortion-matrix routine is unchanged.

The normal log reports ``Ly-alpha 2PCF forest index`` with the forest/node
counts, visited/pruned nodes, candidate forest pairs, per-thread bitmap bytes
and build CPU time. Work counters are local to each MPI rank. Additional
index/bitmap memory is checked with the existing native memory budget.
Run ``make test-lya-forest-index`` for full-extent, bent-forest, strict-boundary,
candidate-order and independent-estimator regressions. Set
``CTREEBALLS_MPI_COMMAND`` to a launcher with an explicit allocation to include
the two-rank comparison. Tree speed depends on forest distribution and cutoff;
compare native timings on the intended catalog.

Quadtree in the two-dimensional separation plane
------------------------------------------------

``lya-rprt-quadtree-2pcf-omp`` is an alternative 2PCF search with exact
default binning and optional bin slop. It uses
a **quadtree** (four quadrants in two dimensions) over the requested
``(r_parallel, r_transverse)`` bin plane. These are pair coordinates, not
absolute coordinates of individual pixels. If the observer distances are
``chi_a, chi_b`` and the angle between their sightlines is ``theta_ab``, then::

   r_parallel   = abs(chi_a - chi_b) * cos(theta_ab / 2)
   r_transverse = (chi_a + chi_b) * sin(theta_ab / 2)

Moving either endpoint changes the separation coordinates. A single global
pixel tree in this plane would therefore be invalid. Building a catalog of
every pixel pair first would instead require quadratic memory and pair work.

The engine builds the bin-plane quadtree once and projects **products of
forest segments** into it lazily. Each projected rectangle bounds every
pixel pair in that product, including the true angular extent of bent
sightlines. The quadtree locates the smallest bin region containing the
whole rectangle. A leaf is one physical histogram bin. If the rectangle fits
inside a leaf and inside the strict discovery cutoff, the entire product
contributes at once::

   numerator   += sum(w * delta)_a * sum(w * delta)_b
   denominator += sum(w)_a * sum(w)_b

Unresolved products split along both radial hierarchies, yielding up to four
child products. At small ambiguous cells the original exact pixel arithmetic
decides the bin in zero-slop mode. Products in that mode are not assigned
by a centre or by a representative pixel. Forest exclusions, original pixel-row ownership, zero-weight pair
counts, and strict upper cutoffs remain intact. Exact bin membership does not
imply bitwise-identical weighted sums to a different traversal; use the
normal numerical comparison tolerance.

Complete-forest bounds use the existing top-level kdtree to reject distant
groups. This engine **does not construct the global Cartesian pixel octree**.
Its added bin-tree storage is at most ``2 * rp_bins * rt_bins - 1`` nodes;
forest hierarchy storage remains linear in the retained pixel count. It never
materializes the full pixel-pair cloud. All allocations are included in the
normal memory preflight. Arbitrary rectangular and non-power-of-two grids,
including a singleton axis, are supported.

Select it directly in the all-engines driver::

   python3 tests/python/lya_corr_all_engines.py \
       --fits Rodrigo/deltas-NoDLAm-NoNoise-True-dz02_LYA.fits \
       --fits-layout radec-r --radec-r-distance-unit Mpc/h \
       --engine lya-rprt-quadtree-2pcf-omp --threads 4 \
       --max-forests 100 --verbose 2

The radius-unit assumption above must match the input product. DESI and the
other existing input layouts work through the same shared reader. The new
method has an intrinsic effective ``lya2Kernel=1``; no kernel flag is required.
It is OpenMP only, compiled with ``LYAFORESTOMPON=1``. Optional
``lya2RpSlop``/``lya2RtSlop`` in ``[0,1]`` permit expanded-bin assignment,
while the physical cuts remain exact. The Python driver exposes
``--bin-slop`` for both axes, and ``--lya2-rp-slop``/``--lya2-rt-slop``
for independent overrides. Defaults remain zero. These are fractions of bin
widths, not correlation-error bounds; calibrate against zero slop. Pivot
frontier/smoothing remain rejected. In precision/fast-math profiles outside
the strict double interval certificates it retains exact pixel evaluation.
The existing 2PCF output and model comparison formats remain unchanged.
Continuum projection is compatible with bin slop. The separate exact hard-bin
distortion response requires zero pair slops. See
``docs/LYA_RPRT_OPTIMIZATION.md`` for the full slop contract and calibration.

``verbose=2`` prints bin-node visits, single-bin products, represented pairs
and remaining pixel tests. ``make test-lya-rprt-quadtree`` checks independent
pixel sums, difficult geometry, boundaries, rectangular grids, whole-product
acceptance, memory recovery and thread consistency. A speed improvement is
data- and bin-dependent; narrow bins and broad segment images can require
many subdivisions. Measure this engine against both original pair kernels.

Continuum distortion matrix for the 2PCF
----------------------------------------

``cyballs.lya_distortion_matrix`` computes the continuum-projection response
in native C with OpenMP. It implements equations (2.5)--(2.7) of
`Lya2pcf (arXiv:2507.00129) <https://arxiv.org/abs/2507.00129>`_:
``xi_measured[A] = sum_B D[A,B] xi_unprojected[B]``. The matrix acts on a
theoretical model; it is not inverted to correct the observed correlation.

The forest driver accepts ``--compute-distortion`` and saves
``distortion/distortion.npz``, ``distortion.fits`` and ``distortion.json``.
With ``--model-correlation model.npy`` it also uses this matrix in the existing
2PCF model comparison. Compute it once per catalog, weights and grid; it is
shared by the pixel, LOS-tree and combined 2PCF/3PCF search engines. MPI runs
compute it on rank zero with OpenMP.

For example, from the directory containing the driver::

   python lya_corr_all_engines.py --fits delta-*.fits.gz \
       --engine lya-los-tree-2pcf-omp --threads 4 \
       --compute-distortion --output Output_with_distortion

The default removes each forest's weighted mean and log-wavelength slope.
``--distortion-projection-order 0`` removes only the mean. The matrix depends
on positions, wavelength and weights, not delta. It does not reproject delta
unless ``--project-delta`` is explicitly supplied. The declared projector must
match the continuum treatment of the input product. Survey-wide mean-flux
subtraction, model redshift evolution, metals, and a 3PCF distortion operator
are outside this calculation.

DESI/eBOSS FITS readers preserve true log wavelength. For input projected
upstream, ``--compute-distortion`` retains the full valid forest and original
``WEIGHT`` for the fitting operator: cuts, thinning, or correlation-weight
evolution affect only the measured-pair weights. With ``--project-delta``,
the retained pixels are explicitly refitted using the evolved weights and
the response uses that fitting domain. NPZ caches saved during computation
preserve both the selected pixels and any full fitting domain. A previously
thinned cache without that domain must be regenerated from FITS. An NPZ may
also supply an aligned ``log_lambda`` array; distance is never substituted
for wavelength. ASCII/Cartesian input without wavelength supports mean-only
computation or the lower-level API with explicit wavelength metadata.

The Python API supports separate data and model grids::

   from cyballs import lya_distortion_matrix
   response = lya_distortion_matrix(
       positions, weights, forest_ids, log_lambda,
       projection_weights=fit_weights, projection_order=1,
       rp_max=200, rt_max=200, rp_bins=50, rt_bins=50,
       model_rp_max=300, model_rt_max=300,
       model_rp_bins=75, model_rt_bins=75,
       number_threads=4, max_memory_mb=1024)
   predicted = response['matrix'] @ model_vector

Supply all fitting-domain pixels. Zero data weights remove pixels from the
measured pairs without removing their projection contribution. Both grids
are folded observer-bisector coordinates in Mpc/h, flattened with parallel
bin major and transverse bin minor. Model correlations are assumed zero
outside the model grid. Increase model support to study truncation; never
renormalize rows to sum to one. With complete model support, mean removal
annihilates a constant model, so row sums should approach zero.

The return value contains ``matrix``, unnormalised ``numerator``, row
``weights``, ``valid_rows``, bin centres, grid definitions and work statistics.
Empty data bins have NaN matrix rows and zero numerator/weight. FITS exports
these as ``DM``, ``RAW_DM``, ``WDM`` and ``VALID``. These are compatible with
``lya_analysis``; the driver model-comparison path uses identical data/model
grids. Rectangular grids are available through the API.

The exact low-rank expansion replaces four pixel-index loops by two pair
geometry passes and sparse binned moment products. Fixed primary-forest tasks
reduce in a fixed order, giving identical outputs across thread counts.
No pixel/forest subsampling or angular approximation is used. Singleton or
constant-wavelength forests fall back to mean-only projection and are counted.
Dense storage grows as the product of the data and model bin counts, with
one accumulator per worker. ``--distortion-max-mib`` (default 512), the API's
``max_memory_mb``, and ``CBALLS_MEMORY_BUDGET_MB`` bound the planned workspace.
Use fewer bins or threads for large matrices. A strict double-precision 3D
Ly-alpha build is required.

The stateless C interface is declared in ``include/lya_distortion.h``;
``source/lya_distortion.c`` owns the implementation. Its output is the raw
numerator and row weights, with caller-owned arrays. The Python binding
normalises valid rows. Scientific regression tests in
``tests/python/test_lya_distortion.py`` compare against a literal four-index
oracle and direct covariance projection, including distinct fitting/data
weights, full-domain selection, rectangular grids, and thread reproducibility.
