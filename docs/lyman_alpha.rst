Lyman-alpha Forests
===================

The forest addons use pixel positions, flux contrast ``delta``, statistical
weights, and integer forest/quasar IDs. They exclude same-forest pairs and
require three different forests for triplets. They do not use scalar-kappa
Fourier multipoles.

Choose a Method
---------------

Enable ``LYAFORESTOMPON=1`` and/or ``LYAFORESTMPION=1``. In each name below,
replace ``omp`` with ``mpi`` for its distributed sibling:

* ``lya-2pcf-omp``, ``lya-3pcf-omp``, ``lya-2pcf-3pcf-omp``:
  anisotropic physical-3D estimators.
* ``lya-1d-2pcf-omp``, ``lya-1d-3pcf-omp``, ``lya-1d-2pcf-3pcf-omp``:
  signed radial-lag estimators ignoring transverse distance.
* ``lya-1d-tree-2pcf-omp``: exact interval-tree radial 2PCF with same-forest subtraction.

All modes still require ``DEFDIMENSION=3`` and ``usePeriodic=false``.
"1D" describes the search coordinate, not the body-storage dimension.
Smooth-pivot is unsupported. Radial and 3D outputs measure different statistics
and should not be treated as interchangeable reference values.

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

Examples
--------

From a checkout with the MPI addon enabled::

   mpiexec -n 2 ./cballs addons/lya_forest_mpi/parameters.ini

The all-engines driver reads DESI delta FITS, NPZ, or six-column ASCII once,
broadcasts arrays once for MPI, and retains registered catalogs between engines::

   python3 python/lya_corr_all_engines.py --list-engines
   bash examples/download_desi_lya_example.sh /tmp/desi-lya
   python3 python/lya_corr_all_engines.py \
       --fits /tmp/desi-lya/delta-1019.fits.gz --max-forests 6 \
       --pixel-stride 30 --engine all-omp --threads 2

The default driver statistic is 2PCF. Add ``--statistics both`` for 3PCF
and ``--engine all --mpi-ranks 2`` for an MPI-inclusive suite. Subsampling,
distance cosmology, and DESI blinding metadata matter scientifically; these
commands are examples, not a production DESI analysis.

See :download:`driver README <../python/README_lya_corr_all_engines.md>`,
:download:`MPI addon README <../addons/lya_forest_mpi/README.md>`,
and ``examples/lya_corr_all_engines.ipynb``.
The script/notebook ``examples/compare_lya_1d_3d`` compares radial and 3D
workflows while keeping their estimator contracts distinct.

Validation
----------

Run ``make test-lya-forest-omp test-lya-forest-1d-omp`` and
``make test-lya-forest-mpi`` for enabled profiles.
The MPI implementation replicates input/tree memory and partitions pivot or
task blocks. Thread-count changes are deterministic at fixed rank count;
changing ranks may change rounding. Native and Cython tests cover independent
oracles, scan/tree agreement, exclusions, and recoverable failures.
