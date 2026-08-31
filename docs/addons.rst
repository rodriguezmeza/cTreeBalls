
Addons
======

Add-ons extend cBalls with optional search methods, catalog formats, and
supporting libraries without making every core file permanently depend on those
features. They are enabled at compile time through ``Makefile_settings`` and the
add-on make fragments under ``addons/``.

For method selection and build switches, start with :doc:`search_methods`.
The scalar angular contract is in :doc:`3pcf`; forests and physical-3D
surveys have separate guides, :doc:`lyman_alpha` and :doc:`scalar_3d`.

Include files
----------------------------------

Core headers and modules contain guarded ``//B socket:`` blocks such as:

.. code-block:: c

   #ifdef ADDONS
   #include "protodefs_include.h"
   #endif

An add-on uses these include sockets to add prototypes, struct fields,
parameter definitions, switch cases, and implementation hooks. Keep add-on
specific declarations in the add-on include files rather than editing unrelated
core code.

Common socket targets include:

``cmdline_defs_include.h``
    Adds new command-line parameters.

``protodefs_include.h``
    Adds prototypes that must be visible across C modules.

``globaldefs_include_*.h`` and ``global_data_include_*.h``
    Adds globals or run-state fields required by the add-on.

``cballs_include_*.h``, ``startrun_include_*.h``, ``cballsio_include_*.h``
    Adds search-method dispatch, startup, and I/O cases.

Source files
----------------------------------

Each add-on should keep its implementation and build fragment in its own
subdirectory, for example ``addons/cfitsio``, ``addons/gadget_io``,
``addons/kdtree_box_omp``, or ``addons/octree_2balls_omp``. The add-on
``Readme.txt`` and ``To_add_in_Makefile_addons`` files explain which fragments
must be included in the top-level add-on build.

When adding a new search add-on, wire it through the same conceptual path as a
core search method:

1. Define the command-line search string and integer method tag.
2. Add parameter parsing or validation hooks if the method needs extra options.
3. Add the method's prototype through the prototype socket.
4. Add the ``EvalHist`` switch case through the ``cballs`` socket.
5. Make sure any method-specific allocations are released by an existing or new
   cleanup hook.

When adding a new I/O add-on, add the format string conversion, the parser or
writer implementation, and the ``InputData``/``OutputData`` switch case through
the I/O sockets.

Octree BALLS4
-------------

The ``octree-balls4-omp`` search is enabled with ``OCTREEBALLS4OMPON=1``.
Its source is ``addons/octree_balls4_omp/search_octree_balls4_omp.c``.
It replaces the name ``octree-kkk-balls4-omp`` while retaining method ID 69;
update parameter files, Python engine names, and make overrides accordingly.
The separate triangle addon retains its existing name.

This method uses the normal tree loader and its dedicated B4 pivot partition,
with ``OCTREESMOOTHINGON=0``, ``BALLSON=0``, and ``BALLS4SCANLEVON=0`` supported.
It requires 3D, logarithmic radial bins, and nonperiodic geometry.
``TWOPCFON=1`` enables 2PCF and ``TPCFON=1`` enables 3PCF.

Smooth-pivot support has been removed from this addon. Global
``SMOOTHPIVOTON`` may be 0 or 1 for other methods without changing this
method's result. An explicit ``options=smooth-pivot`` returns an error in C
and Cython. Ordinary cell aggregation remains controlled by ``no-one-ball``
and ``no-two-balls``.

Run ``make SMOOTHPIVOTON=0 test-octree-balls4-no-smoothing`` and repeat with
``SMOOTHPIVOTON=1``. The regression covers C/Cython, exact and accelerated
searches, multiple threads, masks, and recoverable rejection of smoothing.

Mask and edge corrections
~~~~~~~~~~~~~~~~~~~~~~~~~

``read-mask`` filters pivots and neighbors. Mixed-mask cells are opened,
and count-based density normalization uses the number of valid bodies.
One ``columns-ascii-all`` or ``binary-all`` catalog can carry its own mask;
``cyballs.set_catalog(..., mask=mask)`` accepts the same selection in memory.
The FITS driver can read a separate mask and attach it without rereading the map.

Add ``edge-corrections,no-normalize-HistZeta`` for complex angular window
deconvolution. The corrected path uses body pivots from the B4 partition and
excludes repeated neighbors, including diagonal-bin self terms. It accumulates
signal modes through ``mChebyshev`` and window modes through twice that order.
``weights-norm`` enables weighted signal/window moments; otherwise the 3PCF
window uses unit weights. Neighbor-cell aggregation preserves weight/field
sums and squared sums; either ``no-one-ball`` or ``no-two-balls`` makes this
edge path exact. Empty or singular systems return zero. Cython exposes
``getHistZetaM_EE_complex(m)`` with one-based multipole indices.

Both edge runs and uncorrected ``no-normalize-HistZeta`` runs now use this
distinct-neighbor body-pivot estimator. It is not the legacy cell-pivot 3PCF
followed by a correction. Only the legacy normalized non-edge path retains
the original cell-pivot contract.
The edge path does not require ``NMultipolesON`` or ``NONORMHISTON``.

MPI variant
~~~~~~~~~~~

``OCTREEBALLS4MPION=1`` builds ``octree-balls4-mpi`` (method ID 184),
independently of ``OCTREEBALLS4OMPON``. All ranks construct the same native
tree. Ordinary runs distribute disjoint B4 frontier nodes and reduce raw
sums before normalization; edge runs distribute fixed body-pivot blocks and
reduce them in order, independently of rank/thread count. Only rank 0 publishes
results or writes histograms. ``numberThreads`` is the number of threads per
rank; ``SMOOTHPIVOT`` remains unsupported.

Run ``make test-octree-balls4-edge`` for C/Cython numerical tests and
``make test-octree-balls4-mpi`` for the two-rank comparison.
The latter accepts ``MPIEXEC='mpiexec --oversubscribe'`` when needed.

Two-Ball Angular Edge Corrections
----------------------------------

The shared TreeCorr-style kernel supports ``edge-corrections`` in
``balltree-2balls-omp``, ``balltree-2balls-mpi``,
``balltree-2balls-omp_3pcf``, ``balltree-2balls-mpi_3pcf``,
``octree-2balls-omp``, and ``octree-2balls-mpi``. Enable ``TPCFON=1`` and use
``options=KKKCorrelation,edge-corrections,no-normalize-HistZeta``.
Add ``weights-norm`` for weighted signal and window moments. Native octree
methods can also use ``read-mask``; masking alone does not deconvolve the window.

For each pair of radial bins the kernel accumulates distinct-triplet signal
moments :math:`S_\ell` through :math:`M=\mathtt{mChebyshev}` and window moments
:math:`W_\ell` through :math:`2M`. With :math:`u_i=w_i` when weighted, otherwise
:math:`u_i=1`, their common angular phase is
:math:`\exp[i\ell(\phi_{ij}-\phi_{ik})]`. The signal multiplies
:math:`u_i u_j u_k` by :math:`\kappa_i\kappa_j\kappa_k`; the window does not.
The corrected multipoles solve

.. math::

   \sum_{n=-M}^{M} \frac{W_{\ell-n}}{W_0}\,\zeta_n
   = \frac{S_\ell}{W_0}, \qquad -M \leq \ell \leq M.

Both positive and negative modes, including imaginary parts, enter this
complex Toeplitz system. Negative modes are complex conjugates. This is the
truncated angular-window correction used for scalar fields, not a replacement
for a data-minus-random galaxy estimator or a spherical-harmonic survey matrix.
It retains the engine's existing angular-coordinate convention in 2D and 3D.
Higher true signal modes can leak across any finite truncation; check stability
with increasing ``mChebyshev`` for narrow or complicated masks.
Narrow angular windows can also be ill-conditioned without being singular,
amplifying roundoff and sampling noise. Check convergence with binning and
traversal tolerances as well as multipole order.

LogMultipole searches accumulate signal and window in the same neighbor scan.
The window second-moment subtraction removes ``j==k`` for every order, and
completed window moments are inherited when a pivot splits. Angular acceptance
uses the highest required window mode. The direct triple-node traversal uses
the same window convention. MPI reduces the unnormalized task histograms in
fixed task order and solves on rank 0. Failures are propagated to every rank.

Empty and numerically singular systems produce zero corrected modes, with
counts reported at verbosity 2. No inverse is fabricated for an unsupported
radial-bin pair. Ordinary 2PCF pair-weight normalization is unchanged, and
``only-2pcf`` with edge correction is rejected. ``no-two-balls`` selects the
exact validation path; without edge correction the existing estimator remains.
No ``NMultipolesON`` or ``NONORMHISTON`` build setting is needed for these six
engines' window accumulation.

Output files are ``histZetaM_EE_N.txt`` and ``histZetaM_EE_Im_N.txt``. Raw window
files are ``histZetaM_window_Re_N.txt`` and ``histZetaM_window_Im_N.txt`` with
``N=1`` for order zero. ``no-out-Hist`` suppresses files but not the correction.
Use ``balls.getHistZetaM_EE_complex(N)`` before ``struct_cleanup()`` in Cython.
The in-memory kappa driver enables this with ``--type edge_effects`` and retains
both complex components in its NPZ output.

Run ``make test-two-ball-edge`` for active-engine independent-oracle tests.
The test script also accepts ``--mpi-command 'mpiexec -n 2'``. A 2D build can
be compared to TreeCorr with ``--dimension 2 --treecorr``. The Cython and driver
regressions are in ``tests/make_tests/test_two_ball_edge_cython.py``.

Masked Octree Two-Ball Searches
----------------------------------

``octree-2balls-omp`` and ``octree-2balls-mpi`` support masked scalar
autocorrelation with ``options=read-mask``. A mask value of 1 means valid;
0 excludes a body as both a pivot and a neighbor. Mixed native cells are
opened while constructing the binary search tree, whose centers, radii,
field sums, and squared moments are then computed from valid bodies only.
Binary-tree storage scales with the unmasked count; the input catalog and
native octree still contain all input bodies.

Both 2PCF and 3PCF support the mask, in exact and aggregate traversal.
``weights-norm`` uses statistical weights of valid bodies only. Normalized
3PCF uses distinct valid triplets; ``no-normalize-HistZeta`` returns their
raw sums. One or two valid bodies produce zero 3PCF; an entirely empty mask
returns a recoverable error. Masking is selection, not survey-window edge
correction. In particular, ``and-CF`` is not a masked-survey random-catalog
estimator.

The existing FITS data-plus-mask input remains supported. Alternatively,
a single ``columns-ascii-all`` catalog can carry its own mask:

.. code-block:: text

   # nbody NDIM Lx Ly Lz
   # 4 3 2 2 2
   -0.4 -0.3 -0.2  0.10 1.0 1
    0.4  0.2  0.3  0.20 0.8 1
    0.1 -0.2  0.4 -0.15 1.2 0
   -0.2  0.4 -0.1  0.05 0.9 1

The columns are ``x y z kappa weight mask`` in a 3D build. Run, for example:

.. code-block:: bash

   ./cballs search=octree-2balls-omp in=masked.txt \
       infmt=columns-ascii-all rootDir=Output-mask \
       rangeN=1.5 rminHist=0.01 sizeHistN=8 mChebyshev=3 \
       options=read-mask,weights-norm,compute-HistN,out-m-HistZeta

   mpiexec -n 2 ./cballs search=octree-2balls-mpi in=masked.txt \
       infmt=columns-ascii-all numberThreads=4 rootDir=Output-mask-mpi \
       rangeN=1.5 rminHist=0.01 sizeHistN=8 mChebyshev=3 \
       options=read-mask,weights-norm,compute-HistN,out-m-HistZeta

For an in-memory catalog, no companion file is needed:

.. code-block:: python

   from cyballs import cballs

   balls = cballs()
   balls.set({"searchMethod": "octree-2balls-omp",
              "options": "read-mask,weights-norm,compute-HistN",
              "rangeN": 0.1, "rminHist": 0.001, "sizeHistN": 16})
   balls.set_catalog(positions, kappa=kappa, weights=weights, mask=valid)
   try:
       balls.Run(level=["MainLoop"])
       xi = balls.getHistXi2pcf().copy()
       zeta = balls.getHistZetaMsincos(1, 1).copy()
   finally:
       balls.struct_cleanup()

Supplying ``mask=`` stores the mask; include ``read-mask`` to activate it.
This mode follows the existing masked-autocorrelation input convention,
not the unmasked cross-catalog ``iCatalogs=1,2`` convention.

With ``TWOPCFON=1`` and ``TPCFON=1``, run the regressions using:

.. code-block:: bash

   make test-octree-2balls-mask
   python tests/make_tests/test_octree_2balls_mask.py --cython
   python tests/make_tests/test_octree_2balls_mask.py --cballs ./cballs --mpi
   python tests/make_tests/test_octree_2balls_mask.py --cballs ./cballs --mpi --fits

The Python-only check must import a freshly rebuilt extension. The existing
``test-octree-2balls-omp`` and ``test-octree-2balls-mpi`` runners also include
the masked CLI regressions.
The optional ``--fits`` check additionally requires CFITSIO and ``healpy``.

Lyman-alpha Forest MPI
----------------------

``LYAFORESTMPION=1`` enables seven MPI+OpenMP searches independently of
``LYAFORESTOMPON``:

* ``lya-2pcf-mpi``, ``lya-3pcf-mpi``, and ``lya-2pcf-3pcf-mpi``;
* ``lya-1d-2pcf-mpi``, ``lya-1d-3pcf-mpi``, and ``lya-1d-2pcf-3pcf-mpi``;
* ``lya-1d-tree-2pcf-mpi``.

They share the OpenMP estimators, bins, and output columns. Require a single
``lya-ascii`` catalog with columns ``x y z delta weight forest_id``,
``DEFDIMENSION=3``, and ``usePeriodic=false``. Radial searches deliberately
ignore transverse separation. Correlated pixels must belong to different forests.

From the project root:

.. code-block:: bash

   make -j4 LYAFORESTMPION=1 cballs
   mpiexec -n 2 ./cballs addons/lya_forest_mpi/parameters.ini
   make MPIEXEC='mpiexec --oversubscribe' test-lya-forest-mpi

``numberThreads`` is per rank. Catalog/tree and histogram memory is replicated
on each rank; pivots or task blocks are partitioned cyclically. Raw sums are
reduced before normalization, including the radial tree's same-forest
subtraction. Only rank 0 writes output. Results are deterministic across thread
counts for a fixed rank count; changing rank count can change rounding.

Under Cython, every MPI rank must call ``Run`` and cleanup in the same order.
Import ``mpi4py.MPI`` before running and read output files only on rank 0.
These modes require forest-aware input, so they are deliberately excluded from
the scalar-kappa driver. See ``addons/lya_forest_mpi/README.md`` for details.

The forest-aware driver ``python/lya_corr_all_engines.py`` supports all fourteen
OpenMP/MPI forest engines. It reads DESI DR1 delta FITS, NPZ or six-column ASCII
once, broadcasts arrays once for MPI, and uses
``set_forest_catalog(positions, delta, weights, forest_ids)``. Original integer
IDs remain intact; compact internal IDs preserve same-quasar exclusions.

.. code-block:: bash

   bash examples/download_desi_lya_example.sh /tmp/desi-lya
   python python/lya_corr_all_engines.py --fits /tmp/desi-lya/delta-1019.fits.gz \
       --max-forests 6 --pixel-stride 30 --engine all-omp --threads 2

Add ``--statistics both`` to include 3PCF and ``--engine all --mpi-ranks 2``
for MPI comparisons. The default is 2PCF only. Plots and CSV differences compare
equivalent estimators; radial and anisotropic 3D products are kept separate.
DESI blinding labels, the chosen distance cosmology and subsampling are recorded.
See ``python/README_lya_corr_all_engines.md`` and
``examples/lya_corr_all_engines.ipynb`` for details and limitations.

3D Scalar MPI
-------------

``OCTREE3PCF3DMPION=1`` enables ``octree-3pcf-3d-mpi``, with alias
``octree-ggg-3d-mpi``. It shares the spherical-harmonic kernel and scalar
and ENCORE-style survey estimators with ``octree-3pcf-3d-omp``.
``OCTREE3PCF3DOMPON`` may independently be disabled.

From the project root:

.. code-block:: bash

   make -j4 OCTREE3PCF3DMPION=1 cballs
   mpiexec -n 2 ./cballs addons/octree_3pcf_3d_mpi/parameters.ini
   mpiexec -n 2 ./cballs addons/octree_3pcf_3d_mpi/parameters_survey.ini
   make MPIEXEC='mpiexec --oversubscribe' test-octree-3pcf-3d-mpi

Use ``DEFDIMENSION=3``, OpenMP, and an MPI C compiler. Catalog/tree memory
is replicated and ``numberThreads`` is per rank. Fixed pivot blocks are
partitioned across ranks and locally committed in order. Raw sums, including
the data-minus-random and random multipoles, are reduced before rank 0
normalizes or solves the survey edge-correction matrix.

Select ``only-2pcf-3d``, ``only-3pcf-3d``, or
``compute-2pcf-3d,compute-3pcf-3d``. Survey input uses two catalogs and
``options=survey-estimator-3d``, ``iCatalogs=1,2``.
The output formats and zero/invalid-window policy match the OpenMP sibling.
Only rank 0 writes results. Thread-count changes are deterministic at a fixed
rank count; changing the rank count can change floating-point rounding.

Cython supports file and in-memory catalogs. Under ``mpiexec``, import
``mpi4py.MPI`` first, provide complete catalogs on every rank, and call
``Run`` and cleanup in the same order. This physical 3D estimator is not an
angular-kappa engine. See ``addons/octree_3pcf_3d_mpi/README.md``.
