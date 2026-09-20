Add-ons
*******

Add-ons are selected in ``addons/Makefile_addons_settings``. This page covers
the maintained profile whose switches are 1. Use
``options=print-search-methods`` after compilation for the exact runtime list.

Angular Scalar Trees
--------------------

The active KD-tree, PCA-ball-tree, and native-octree families each provide an
OpenMP and MPI method:

* ``kdtree-2balls-omp`` and ``kdtree-2balls-mpi``;
* ``balltree-2balls-omp`` and ``balltree-2balls-mpi``;
* ``octree-2balls-omp`` and ``octree-2balls-mpi``.

They compute 2PCF and LogMultipole 3PCF, selected with ``only-2pcf`` or
``only-3pcf``. Dual-node acceptance is conservative by default;
``dual-node-bin-slop`` enables bin-aware approximation and ``no-two-balls``
selects exact body pairs. Exact unsmoothed 3PCF uses
``no-one-ball,no-smooth-pivot``; native octree 3PCF is not exact with
``no-two-balls`` alone. ``BALLS4SCANLEVON=1`` supplies balanced task frontiers.

All six support ``read-mask`` and complex scalar edge correction. Enable
``TPCFON=1`` and use::

   options=KKKCorrelation,edge-corrections,no-normalize-HistZeta

Add ``weights-norm`` for weighted signal and window moments. The correction
uses signal modes through :math:`M` and window modes through :math:`2M` in a
complex Toeplitz solve. Empty or singular scalar systems publish NaN with
validity and conditioning diagnostics. Masking is
selection; it does not itself request edge correction.

The KD and PCA ball-tree families use smooth pivots by default when
``SMOOTHPIVOTON=1``. ``no-smooth-pivot`` disables this. The native-octree
dual-node path retains exact body pivots.

Full-Sky Spin-2 Trees
---------------------

The active shear methods are:

* ``octree-shear-sphere-2balls-omp``;
* ``kdtree-shear-sphere-2balls-omp``;
* ``balltree-shear-sphere-2balls-omp``.

They share unit-sphere normalization, great-circle parallel transport, local
east/north shear input, 2PCF normalization, natural 3PCF multipoles, masks,
edge correction, order selectors, and deterministic OpenMP reductions. Their
tree partition and enclosing-cell geometry differ. See :doc:`shear`.

Lyman-alpha Forests
-------------------

``LYAFORESTOMPON=1`` and ``LYAFORESTMPION=1`` provide anisotropic 3D, radial
scan, and radial interval-tree 2PCF/3PCF estimators. The same-LOS equal-forest
2PCF is OpenMP only. Forest-aware input preserves integer IDs and enforces the
documented pair/triplet exclusions. See :doc:`lyman_alpha` and the addon
READMEs under ``addons/lya_forest_omp`` and ``addons/lya_forest_mpi``.

Physical 3D Multipoles
----------------------

``OCTREE3PCF3DOMPON=1`` and ``OCTREE3PCF3DMPION=1`` enable
``octree-3pcf-3d-omp`` and ``octree-3pcf-3d-mpi``. They compute Legendre
multipoles and can form data-minus-random and random multipoles for the survey
window solve. See :doc:`scalar_3d`.

Input and Utility Add-ons
-------------------------

The maintained profile also enables ``GADGETIOON``, ``CLASSLIBON``,
``PXDON``, ``IOLIBON``, and ``CFITSIOON``. CFITSIO is discovered as an
external dependency. ``kdtree-box-omp`` and ``neighbor-boxes-omp``
provide periodic Cartesian 2PCF methods.

Parallel Rules
--------------

MPI methods require the same MPI implementation in ``MPICC``, ``mpiexec``,
the C extension, and ``mpi4py``. ``numberThreads`` is per rank. Every rank
enters matching run and cleanup collectives; rank 0 writes output. OpenMP and
MPI methods reduce raw sums before normalization.

Validation
----------

Run the focused launchers from ``tests/make_tests`` and the driver contracts::

   make test-search-methods test-two-ball-edge
   python3 -m pytest -q tests/make_tests/test_kappa_corr_all_engines.py
   python3 -m pytest -q tests/make_tests/test_shear_corr_all_engines.py
   python3 -m pytest -q tests/make_tests/test_lya_corr_all_engines.py

Tests requiring MPI or large external catalogs remain opt-in.
