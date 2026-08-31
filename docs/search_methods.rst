Search Methods
==============

Start with the field and geometry, then choose an enabled engine. This guide
covers the maintained families; the executable is the authority for a particular
build::

   ./cballs options=make-info
   ./cballs options=print-options
   ./cballs options=print-search-methods

The first command reports compiled settings, the second registered options,
and the third method names, availability, and usage. A source directory alone
does not make a method available. See :doc:`build_profiles` before changing flags.

Scalar Angular Engines
----------------------

These measure counts/convergence using scalar Fourier multipoles, not shear
or the physical-3D Legendre estimator. In 3D preserve the observer origin and
use chord separations. See :doc:`3pcf` for the precise raw estimator.

.. list-table::
   :header-rows: 1
   :widths: 27 30 43

   * - Runtime name
     - Build switch
     - Search and controls
   * - ``octree-ggg-omp`` / ``octree-ggg-mpi``
     - ``OCTREEGGGOMPON`` / ``OCTREEGGGMPION``
     - Pivot-neighbor octree moments; log bins; ``no-one-ball`` for exact neighbors.
   * - ``kdtree-omp``
     - ``KDTREEOMPON``
     - KD-tree; exact by default, ``behavior-ball`` enables aggregation with log bins.
   * - ``balltree-omp`` / ``balltree-mpi``
     - ``BALLTREEOMPON`` / ``BALLTREEMPION``
     - FCFC-style ball tree; exact by default; ``behavior-ball`` enables aggregation with log bins.
   * - ``octree-balls4-omp`` / ``octree-balls4-mpi``
     - ``OCTREEBALLS4OMPON`` / ``OCTREEBALLS4MPION``
     - Native B4 partition; 3D, nonperiodic, log bins; raw mode uses distinct body-pivot moments.
   * - ``balltree-2balls-omp`` / ``balltree-2balls-mpi``
     - ``BALLTREE2BALLSOMPON`` / ``BALLTREE2BALLSMPION``
     - Dual-node 2PCF and explicit triple-node 3PCF; exact triples can have cubic cost.
   * - ``octree-2balls-omp`` / ``octree-2balls-mpi``
     - ``OCTREE2BALLSOMPON`` / ``OCTREE2BALLSMPION``
     - Native-octree view; dual-node pairs and production LogMultipole 3PCF.
   * - ``balltree-2balls-omp_3pcf`` / ``balltree-2balls-mpi_3pcf``
     - ``BALLTREE2BALLSOMP3PCFON`` / ``BALLTREE2BALLSMPI3PCFON``
     - 3PCF-only LogMultipole scans with inherited coarse-pivot moments.

Set the relevant switch to ``1``. ``TWOPCFON`` and ``TPCFON`` control compiled
correlation orders. Where supported, ``only-2pcf`` and ``only-3pcf`` select
runtime work; they are not universal options for every historical engine.
The ``_3pcf`` addons do not provide 2PCF. Use ``no-two-balls`` for the exact
two-ball reference and reduce ``theta`` to check convergence of aggregation.

``nsmooth`` generally sets leaf capacity in KD/ball trees; it does not imply
pivot smoothing. ``SMOOTHPIVOTON=1`` only compiles support. Supported engines
require an explicit ``smooth-pivot`` option. Raw KD/legacy balltree multipoles
reject smoothing, and BALLS4 rejects it in every mode.

Masks and Windows
`````````````````

GGG, octree two-ball, and BALLS4 OpenMP/MPI engines support ``read-mask``.
Masking removes invalid pivots and neighbors; it does not perform edge correction.

Complex scalar edge corrections are available in GGG, BALLS4, and all listed
two-ball engines. Use ``edge-corrections,no-normalize-HistZeta``;
``weights-norm`` selects statistical weights. GGG also needs
``NMultipolesON=1`` and ``NONORMHISTON=1``. The two-ball and BALLS4 correction
paths do not need these extra switches. Windows extend through twice the
signal multipole order. Empty or singular correction systems return zero.

The ``python/kappa_corr_all_engines.py`` driver filters engine groups by
compiled availability, masks, and edge support. Explicit incompatible requests
fail. It is an angular scalar driver, not a generic frontend for every family.

Other Field Families
--------------------

.. list-table::
   :header-rows: 1

   * - Family
     - Build switches
     - Guide
   * - ``octree-shear-omp``
     - ``OCTREESHEAROMPON=1``
     - :doc:`shear`: flat-sky spin-2 natural components and window correction.
   * - Seven ``lya-...-omp`` and seven ``lya-...-mpi`` modes
     - ``LYAFORESTOMPON=1`` / ``LYAFORESTMPION=1``
     - :doc:`lyman_alpha`: forest IDs, radial-only and anisotropic 3D estimators.
   * - ``octree-3pcf-3d-omp`` / ``octree-3pcf-3d-mpi``
     - ``OCTREE3PCF3DOMPON=1`` / ``OCTREE3PCF3DMPION=1``
     - :doc:`scalar_3d`: physical 3D scalar and data/random survey estimators.
   * - Box and historical scalar engines
     - Profile-dependent
     - Use executable help and the corresponding addon README; do not infer angular contracts.

Aliases ``octree-ggg-3d-omp`` and ``octree-ggg-3d-mpi`` select the physical
3D family, not the angular GGG engines.

MPI Rules
---------

MPI engines require an MPI compiler/runtime and OpenMP support. Thread count
is **per rank**. Catalogs and trees are generally replicated, so MPI does not
automatically reduce per-rank memory requirements. Only rank 0 publishes output.

Every Python rank must enter the same run and cleanup stages. Import
``mpi4py.MPI`` first and use the same MPI implementation as the native build.
The all-engines drivers broadcast catalogs; calling ``set_catalog`` directly
requires complete input on every rank. Do not read rank-0-only getters elsewhere.

Fixed task reductions support repeatability, but not every family promises
bitwise equality when changing MPI rank count. Compare numerical tolerances,
record rank/thread counts, and follow the addon-specific parallel contract.
