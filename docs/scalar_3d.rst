Physical 3D Scalar and Survey Estimators
========================================

``octree-3pcf-3d-omp`` and ``octree-3pcf-3d-mpi`` measure physical
three-dimensional scalar correlations using spherical-harmonic moments and
Legendre multipoles. They are not the observer-tangent Fourier estimator
described in :doc:`3pcf`, and are excluded from the angular-kappa driver.

Enable ``OCTREE3PCF3DOMPON=1`` and/or ``OCTREE3PCF3DMPION=1`` with
``DEFDIMENSION=3``. The ``octree-ggg-3d-omp``/``mpi`` names are aliases.
The default computes 3PCF; select ``only-2pcf-3d``, ``only-3pcf-3d``,
or ``compute-2pcf-3d,compute-3pcf-3d`` explicitly when comparing timings.

Scalar Input
------------

Use ``x y z delta weight`` with ``multi-columns-ascii``,
``columns=1,2,3,4,5`` and ``options=pos-and-convergence-weight``,
or register the corresponding arrays with ``set_catalog``.
Distances use the physical coordinate units. ``rangeN`` remains the exact
outer cutoff even when ``Rcut/theta`` is selected.

Survey Input and Correction
---------------------------

Select ``options=survey-estimator-3d`` and provide data then random catalogs
with ``iCatalogs=1,2``. Both use the same Cartesian frame. Nonnegative
weights define ``alpha=sum(w_data)/sum(w_random)``, the signed field
``N=D-alpha*R``, and the scaled random field ``R=alpha*R``.
The scalar value column is ignored in this mode.

The implementation accumulates N and R multipoles separately, forms the 2PCF
ratio, and solves the 3PCF survey-window coupling matrix after global
reduction. It does not average rank-local normalized correlations.
The random catalog represents the selection; survey mode rejects
``read-mask``, ``remove-mean``, and line-of-sight exclusions.

Empty windows return zero with ``valid=0``; ``pivot_condition`` remains
infinite when no matrix was solved. Inspect validity and conditioning rather
than interpreting every zero as a measured zero correlation.

Examples and Reference
----------------------

With the MPI profile enabled, the bundled tiny examples run from the root::

   mpiexec -n 2 ./cballs addons/octree_3pcf_3d_mpi/parameters.ini
   mpiexec -n 2 ./cballs addons/octree_3pcf_3d_mpi/parameters_survey.ini

Use ``examples/compare_octree_3pcf_3d_encore.py --help`` and the matching
``.ipynb`` for ENCORE comparisons. The general and focused CPU benchmarks are
described in :doc:`benchmarks`.

Detailed formulas, outputs, and options are in the
:download:`OpenMP README <../addons/octree_3pcf_3d_omp/README.md>` and
:download:`MPI README <../addons/octree_3pcf_3d_mpi/README.md>`.

MPI catalog/tree memory is replicated and threads are per rank. In direct
Python MPI use, register full data/random arrays as catalog slots 0/1 on
every rank and enter the same lifecycle. Only rank 0 writes outputs.
Changing rank count may change floating-point grouping.

Run ``make test-octree-3pcf-3d-omp`` and
``make test-octree-3pcf-3d-mpi`` for enabled profiles. Their regressions
cover scalar/survey oracles, cutoff correctness, multiple ranks and threads,
and input/output failure recovery.
