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
  radial-only methods, transverse distance matters. These have no MPI siblings.
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
workflows while keeping their estimator contracts distinct.

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
oracles, scan/tree agreement, exclusions, strict in-memory multipoles, and
recoverable failures.
