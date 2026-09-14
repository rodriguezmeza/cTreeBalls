Search Methods
==============

The executable is the authority for the current build profile::

   ./cballs options=make-info
   ./cballs options=print-options
   ./cballs options=print-search-methods

The current maintained profile exposes the core method and only the enabled
addon methods below. An addon source directory does not make it available
unless its Makefile switch is 1.

Scalar Angular Methods
----------------------

``octree-sincos-omp`` is the core octree method. It computes scalar 2PCF and
sine/cosine 3PCF multipoles with OpenMP and is not controlled by an addon
switch.

.. list-table::
   :header-rows: 1
   :widths: 28 30 42

   * - Runtime name
     - Build switch
     - Contract
   * - ``kdtree-2balls-omp`` / ``kdtree-2balls-mpi``
     - ``KDTREE2BALLSOMPON`` / ``KDTREE2BALLSMPION``
     - Median KD-tree dual-node 2PCF and body-pivot LogMultipole 3PCF.
   * - ``balltree-2balls-omp`` / ``balltree-2balls-mpi``
     - ``BALLTREE2BALLSOMPON`` / ``BALLTREE2BALLSMPION``
     - PCA ball-tree dual-node 2PCF and body-pivot LogMultipole 3PCF.
   * - ``octree-2balls-omp`` / ``octree-2balls-mpi``
     - ``OCTREE2BALLSOMPON`` / ``OCTREE2BALLSMPION``
     - Native-octree dual-node 2PCF and LogMultipole 3PCF.

For the addon methods, ``TWOPCFON`` and ``TPCFON`` compile the two correlation orders.
``only-2pcf`` and ``only-3pcf`` select work at runtime. The default dual-node
acceptance requires the complete pair-distance interval to remain inside one
radial bin. ``dual-node-bin-slop`` enables the less conservative Log/Linear
bin-position policy, and ``no-two-balls`` requests exact body pairs.

All six methods accept masks and complex scalar 3PCF edge correction. Use::

   options=read-mask,edge-corrections,no-normalize-HistZeta

``weights-norm`` applies catalog weights to signal and window moments. Empty or
singular correction systems publish finite zero. The KD and PCA ball-tree
methods support the compiled smooth-pivot default; ``no-smooth-pivot`` disables
it. Native octree dual-node mode does not use smooth pivots.

Full-Sky Shear Methods
----------------------

.. list-table::
   :header-rows: 1

   * - Runtime name
     - Build switch
     - Tree
   * - ``octree-shear-sphere-2balls-omp``
     - ``OCTREESHEARSPHERE2BALLSOMPON``
     - Native octree
   * - ``kdtree-shear-sphere-2balls-omp``
     - ``KDTREESHEARSPHERE2BALLSOMPON``
     - Median KD tree
   * - ``balltree-shear-sphere-2balls-omp``
     - ``BALLTREESHEARSPHERE2BALLSOMPON``
     - PCA ball tree

These consume observer-centered three-dimensional vectors and
``gamma1+i*gamma2`` in each point's local east/north basis. They normalize
positions to the unit sphere, bin chord distance, parallel transport spin-2
fields along great circles, and compute xi+/xi- plus natural 3PCF multipoles.
See :doc:`shear`.

With ``BALLS4SCANLEVON=1``, ``balltree-shear-sphere-2balls-omp`` schedules its
body-pivot 3PCF through an adaptive PCA-cell frontier.  It repeatedly splits
the task with the largest estimated pivot-by-neighbor work, prefilters a
256-node neighbor multipole frontier by conservative spherical-ball overlap,
and executes the resulting tasks dynamically.  Each task owns a private
histogram, capped collectively at 256 MiB, and task histograms are merged in
spatial order after the parallel region for deterministic results.  The
``only-2pcf`` path continues to use the symmetric dual-tree scan.

Lyman-alpha Forest Methods
--------------------------

``LYAFORESTOMPON=1`` enables nine OpenMP names. ``LYAFORESTMPION=1`` enables
the eight MPI counterparts; the same-LOS 2PCF is OpenMP only.

The families are:

* ``lya-2pcf-*``, ``lya-3pcf-*``, and ``lya-2pcf-3pcf-*`` for anisotropic 3D
  forest statistics;
* ``lya-1d-2pcf-*``, ``lya-1d-3pcf-*``, and
  ``lya-1d-2pcf-3pcf-*`` for radial scans;
* ``lya-1d-tree-2pcf-*`` and ``lya-1d-tree-3pcf-*`` for interval-tree radial
  scans;
* ``lya-1d-tree-same-los-2pcf-omp`` for an equal-forest average of within-LOS
  pairs.

Input is ``x y z delta weight forest_id`` or an in-memory
``set_forest_catalog`` call. General forest pairs exclude equal IDs and forest
triplets require three distinct IDs. See :doc:`lyman_alpha`.

Physical 3D and Box Methods
---------------------------

``octree-3pcf-3d-omp`` and ``octree-3pcf-3d-mpi`` are enabled by
``OCTREE3PCF3DOMPON`` and ``OCTREE3PCF3DMPION``. They provide physical-3D
Legendre multipoles and the data/random survey-window estimator described in
:doc:`scalar_3d`.

``kdtree-box-omp`` and ``neighbor-boxes-omp`` are enabled for periodic
Cartesian 2PCF workloads. They are not angular convergence or forest
estimators.

MPI Rules
---------

MPI engines require one MPI implementation shared by the compiler wrappers,
runtime launcher, C extension, and ``mpi4py``. ``numberThreads`` is per rank.
Every rank enters the same run and cleanup sequence; rank 0 publishes output.
The all-engines drivers under ``tests/python`` load a catalog once on rank 0,
broadcast it, and register one retained copy per rank.
