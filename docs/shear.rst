Weak-Lensing Shear
==================

The active shear add-ons measure 2PCFs and the four natural shear 3PCF
components in the Porth x-projection, :math:`\Gamma^\times`. They keep
unit-sphere positions, parallel transport every neighbor shear into the pivot
tangent plane, and use dual-node tree traversal. Enable the octree, KD-tree,
and ball-tree implementations when building C and Cython:

.. code-block:: bash

   make OCTREESHEARSPHERE2BALLSOMPON=1 \
      KDTREESHEARSPHERE2BALLSOMPON=1 \
      BALLTREESHEARSPHERE2BALLSOMPON=1 PYTHON=python3 all

The octree two-ball addon uses its native dual-node behavior for pair work and
its spherical pivot-ring frontier for 3PCF work. Both paths support masks,
edge correction, normalization, and deterministic ``BALLS4SCANLEV`` reduction.

Scientific Convention
---------------------

The estimator uses a right-handed Cartesian tangent plane. Angles increase
counterclockwise from the positive x axis, and input shear is
``gamma = gamma1 + i gamma2``. For pair direction ``phi`` the code accumulates

.. math::

   \xi_+ = \frac{\sum_{ij} w_i w_j\,\gamma_i\gamma_j^*}
                  {\sum_{ij} w_i w_j}, \qquad
   \xi_- = \frac{\sum_{ij} w_i w_j\,\gamma_i\gamma_j
                         e^{-4\mathrm{i}\phi_{ij}}}
                  {\sum_{ij} w_i w_j}.

For the 3PCF, ring sums ``G_n`` generate the raw natural-component
multipoles ``Upsilon``. The implementation is specifically the Porth
x-projection with :math:`\zeta_1=(\varphi_1+\varphi_2)/2`,
:math:`\zeta_2=\varphi_1`, and :math:`\zeta_3=\varphi_2`. It does not return
the centroid projection. Convert projections before comparing with software
or measurements that use the centroid convention.

Same-neighbor diagonal products are subtracted when both triangle legs use the
same catalog, and window multipoles ``N_n`` are used to solve the finite-survey
mode-coupling system for ``Gamma_0^x`` through ``Gamma_3^x``.
``mChebyshev`` is the maximum positive Fourier order for this add-on, so the
returned multipole axis has ``2*mChebyshev + 1`` entries.

Catalog selection supports tomography. One input catalog computes
``Z1,Z1,Z1``; two compute ``Z1,Z2,Z2``; three or more compute ``Z1,Z2,Z3``
from the first three entries of ``iCatalogs``. The 2PCF uses ``Z1,Z2``. The
core parser requires one ``iCatalogs`` entry per input catalog even though this
estimator ignores entries after the third. A single common translation is
applied before the trees are built, so cross-catalog separations are preserved.

Full-Sky Geometry
-----------------

The spherical engine requires a 3D build and nonzero observer-centered input
vectors. It normalizes positions to the unit sphere. Input
``gamma1 + i*gamma2`` is expressed in the right-handed local east/north basis.
For a pivot :math:`\boldsymbol p` and neighbor :math:`\boldsymbol q`, the pair
bearing is the normalized projection of :math:`\boldsymbol q` into the pivot
tangent plane. The neighbor shear is parallel transported along the shortest
great circle into that same basis before the pair and ring products are formed.

Spherical radial bins use chord distance
:math:`r=|\boldsymbol p-\boldsymbol q|=2\sin(\alpha/2)`, so ``rangeN`` cannot
exceed 2. Angular command-line limits are converted to this coordinate by the
Python driver; ``--sep-units`` accepts ``arcmin``, ``degree``, and ``radian``.
The 3PCF uses a separate pivot tangent plane at every sky
position, matching the local-plane full-sky construction used by the projected
scalar estimator.

Spherical tree cells carry basis-aware first and second shear moments. Child
moments are parallel transported into the tangent basis at each cell center,
and accepted-node bounds include the accumulated transport error.

The two-ball addon scans pairs of octree nodes. For node radii :math:`s_1`
and :math:`s_2` and center chord separation :math:`d`, logarithmic-bin
acceptance requires :math:`(s_1+s_2)/d` to fit the ``theta``-scaled logarithmic
bin width; linear bins use the corresponding absolute-width test. Radial
bounds, spin phase, and transport-error tests must also pass. Failed pairs
split the larger node and may split both comparable nodes using dual-node's
empirical ``0.585`` rule. ``no-two-balls`` and ``no-one-ball`` force the exact
body-pair limit. Combined runs use this pair scan for 2PCF and retain the
spherical pivot-ring scan for the compatible Gamma-x LogMultipole 3PCF.

With ``SMOOTHPIVOTON=1``, smoothing is available with explicit ``smooth-pivot``. Each claimed body's
weighted shear is parallel transported into the representative pivot's tangent
basis before the group moment is accumulated. An explicit ``rsmooth`` is in
arcmin and is converted to unit-sphere chord distance without applying the
tree-opening ``THETA``. It must satisfy ``2*rsmooth <= rminHist`` so no smooth
group can contain a pair from the measured domain. Automatically selected
radii already use tree-coordinate units and are capped at the same limit.
Smoothed ``only-2pcf`` runs use the body-pivot traversal because the symmetric
dual tree has no unique pivot side. Add ``options=no-smooth-pivot`` to retain
raw pivots and the faster dual-tree 2PCF path.

``kdtree-shear-sphere-2balls-omp`` applies the same spin-2 estimator and
acceptance criteria to an independently built median KD tree. Its 2PCF uses
the dual-node two-node split rule; its 3PCF scans accepted KD nodes from each
pivot and uses their transported first and second shear moments. ``nsmooth``
sets the KD leaf capacity. Masks, edge correction, ``only-2pcf``,
``only-3pcf``, smooth pivots, and deterministic OpenMP block reductions match
the spherical octree engines.

``balltree-shear-sphere-2balls-omp`` replaces the median axis-aligned KD
partition with FCFC's principal-axis median partition and conservative
spherical enclosing balls. It shares the KD engine's dual-node and accepted
neighbor-node traversal, so masks, normalization, edge correction, runtime
order selection, smoothing, and deterministic reduction have the same
behavior. ``nsmooth`` sets the ball-tree leaf capacity.

Pivot-Cell Reuse
----------------

The octree, KD-tree and PCA ball-tree OpenMP and MPI methods support
``options=shear-pivot-reuse,no-smooth-pivot,only-3pcf`` with
``BALLS4SCANLEVON=1``. They inherit sparse unresolved neighbor lists, keep
accepted rings in their original frames, and transport directly to the final
pivot frame. Safe completed cells use aggregate pivot moments. The separate
2PCF path is unchanged.

``CBALLS_SHEAR_PIVOT_TOL`` is a finite phase budget in radians in ``[0,3]``
(default ``3``), not a relative-error bound on a correlation coefficient.
``CBALLS_SHEAR_BIN_THETA`` allows internal bin-assignment slop in ``[0,1]``
bin widths (default ``0.5``); zero restores full radial-bin containment.
The minimum and maximum radial cutoffs remain strict. Spherical transport
bounds still apply. Zero phase budget disables reuse; the zero opening
tolerance also disables ordinary approximate neighbor acceptance.
With reuse active, positive ``theta`` enables the path but its magnitude does
not replace this phase budget. Zero budget only disables reuse. Exact
validation uses ``theta=0`` and ``no-one-ball,no-two-balls,no-smooth-pivot``. Smoothing, partial pivot
coverage and nonpositive ``theta`` disable reuse; compatibility mode rejects it.

The defaults live in ``addons/Makefile_addons_settings`` as conditional Make
assignments::

    THETA ?= 1
    CBALLS_SHEAR_PIVOT_TOL ?= 3
    CBALLS_SHEAR_BIN_THETA ?= 0.5

For example, ``make THETA=0.8 CBALLS_SHEAR_PIVOT_TOL=0.1
CBALLS_SHEAR_BIN_THETA=0`` overrides them before compilation. These values
are recorded by ``make-info``. Environment variables of the same shear names
override the compiled values at runtime. Unset them to use the compiled
defaults. Runtime ``theta`` (also default 1; ``--tree-theta`` in the Python
drivers) is separate from uppercase ``THETA``, the legacy scan-table control.
Native ``nsmooth`` and the active shear/kappa Python drivers default to 8;
explicit ``--nsmooth`` values override this. This changes binary leaf capacity,
not the shear octree's individual-body leaves, and does not enable smoothing.
None of these settings enables ``shear-pivot-reuse`` automatically. Kappa
engines do not read the shear environment controls; scalar reuse has its own
``CBALLS_SCALAR_*`` controls.

Calibrate every tree/catalog/bin/multipole combination against its own exact
reference. Cancellations and window conditioning can amplify errors, and sparse
catalogs may not benefit. ``CBALLS_SHEAR_PROFILE=1`` prints per-thread phase
times and ``SHEAR_REUSE`` counters. The ball-tree builder fuses deterministic
center/PCA statistics and can temporarily cache source tangent frames (at most
256 MiB, freed before traversal); ``no-balltree-shear-member-cache`` tests the
equivalent uncached path. These construction changes do not loosen accuracy.

The native octree two-ball engines cache normalized cell centers and tighten
spherical radii during cell-moment preparation. Each radius encloses the union
of its live children's bounds by the triangle inequality, with outward
roundoff padding, and is intersected with the existing cube bound. This reduces
unnecessary descent on sparse spherical footprints. The cached geometry is
used by dual-node pair work and pivot reuse in both OpenMP and MPI builds.
Ordinary neighbor walks retain their existing radius acceptance rule; their
cell-to-pivot transport also uses the cached centers. It is rebuilt for each search, including masked catalogs.

This does not change ``theta``, the phase budget or bin slop. Tighter valid
bounds can admit different approximate cells, so approximate products need not
match an earlier build bit for bit; retain an exact-reference calibration for
the chosen controls. Exact body sums and strict outer radial cutoffs retain
their existing behavior. The cache adds 32 bytes per native cell in a 3D double build
with shear enabled, and none per body. Rebuild the executable, static library
and Python extension together after this cell-layout change.

MPI pivot reuse and task scheduling
-----------------------------------

All three ``shear-sphere-2balls-mpi`` engines support ``shear-pivot-reuse``
with ``BALLS4SCANLEVON=1``. Whole pivot subtrees, including inherited rings and
completed-bin ownership, stay within a task on one rank. Each rank reduces
its local task histograms before the global reduction; rank zero alone
normalizes and solves the window system. The phase/bin controls and effective
reuse eligibility must agree across ranks. Invalid controls and local reuse
failures are agreed collectively before histogram reduction, including idle
ranks. No MPI calls occur within the OpenMP traversal.

The MPI pivot target is at least eight tasks per worker (OpenMP threads times
MPI ranks), with a minimum of 64 for binary trees and 256 for the native octree.
It is limited by tree leaves and histogram memory rather than a 256-task
ceiling. ``CBALLS_SHEAR_TASKS`` sets an explicit global target (0 selects auto).
``CBALLS_SHEAR_TASK_MEMORY_MB`` sets the task-histogram budget per rank in MiB
(default 256); one histogram is the unavoidable minimum. Catalogs, trees,
neighbor lists, output histograms and per-thread reuse scratch are additional
memory. A native octree split can exceed the requested target by up to six
tasks, while still respecting the rank histogram capacity.

Longest estimated tasks are assigned first to the least-loaded rank with
remaining histogram capacity; each rank schedules its tasks dynamically among
OpenMP threads. Zero-work ties favor ranks with fewer tasks. The profile shows
``frontier_target``, ``frontier_tasks``, ``local_tasks``, ``histogram_bytes`` and
``rank_task_limit`` per rank. The frozen ``shear_frontier`` run metadata
records the effective target, global/local task counts and local histogram bytes. Ordinary OpenMP runs retain their existing task
target unless the new environment controls are supplied.

Changing the task partition can change where approximate reuse begins. For
close numerical comparisons between layouts, set the same explicit
``CBALLS_SHEAR_TASKS`` and a sufficient memory budget, with identical trees,
bins, phase/bin controls and catalogs. MPI summation still introduces roundoff.
More tasks may improve rank balance but can reduce reuse; benchmark both time
and accuracy. A fixed-core MPI run need not beat OpenMP: catalogs and trees
remain replicated, and launch, construction, reduction and driver I/O add cost.
Use a new output directory when changing environment controls during tuning.

On Perlmutter, benchmark CPU nodes with the same total core budget first,
then add nodes. The Python MPI launcher creates an uncompressed temporary catalog packet
to avoid compression and repeated decompression on every rank;
``TMPDIR`` must be a shared scratch directory visible at the same path on every
node. Launch the batch Python driver once with ``--mpiexec srun`` and set
``--mpi-ranks`` and ``--threads`` (threads per rank). Do not run one parent
Python driver on every rank. Use a consistent MPI/mpi4py/native-extension
software stack and export the same controls to every node. This is a CPU
kernel; it does not use GPUs.

Live progress
-------------

Spherical shear frontier searches print completed 3PCF tasks, percentage, and
elapsed traversal time when ``verbose >= 1``. This includes the ordinary
body-pivot path and ``shear-pivot-reuse`` with ``BALLS4SCANLEVON=1``. Each
MPI rank reports its own local tasks. Task costs differ: the percentage is
completed tasks, not a runtime fraction or a count of galaxy triples. The
100% line marks traversal completion; histogram reduction and final output
still follow. Pair-only searches do not print 3PCF progress.

``verbose=0`` suppresses terminal progress; ``verbose_log >= 1`` independently
routes it to the native log. ``CBALLS_SHEAR_PROFILE=1`` remains optional and
adds detailed phase/reuse statistics at completion. Progress uses a single
synchronized update per frontier task, outside the neighbor and multipole
loops, and does not change histogram accumulation or merge order.

``tests/python/run_desy3_shear_catalogs.py`` forwards each child's output live
to the terminal while retaining the same bytes in its per-catalog ``run.log``.
The MPI launcher also forwards worker output live while enforcing its timeout.
Interrupting the driver terminates and reaps its child. New code/extension
fingerprints require a new output directory for previously completed runs.

Python Example
--------------

.. code-block:: python

   import numpy as np
   from cyballs import cballs

   positions = np.zeros((64, 3))
   angle = np.linspace(0.0, 2.0*np.pi, 64, endpoint=False)
   positions[:, 0] = 0.8*np.cos(angle)
   positions[:, 1] = 0.8*np.sin(angle)
   gamma = 0.03*np.exp(2j*angle)

   model = cballs()
   model.set({
       "searchMethod": "octree-shear-sphere-2balls-omp",
       "iCatalogs": "1",
       "usePeriodic": "false",
       "useLogHist": "false",
       "rminHist": 0.01,
       "rangeN": 2.0,
       "sizeHistN": 8,
       "sizeHistPhi": 16,
       "mChebyshev": 4,
       "lengthBox": 2.0,
       "numberThreads": 4,
       "rootDir": "Output_shear",
       "options": "no-out-Hist",
   })
   model.set_catalog(positions, gamma1=gamma.real, gamma2=gamma.imag)
   model.Run(level=["MainLoop"])

   xi_plus = model.getShearXiPlus()
   xi_minus = model.getShearXiMinus()
   gamma_m = model.getShearGammaXMultipoles()
   gamma_phi = model.getShearGammaX()
   model.clean_all()

The array shapes are ``(B,)`` for each 2PCF, ``(4, 2*nmax+1, B, B)`` for
corrected x-projection multipoles, and ``(4, P, B, B)`` for angular
x-projection natural components. The old ``getShearGammaMultipoles()`` and
``getShearGamma()`` names remain compatibility aliases. Read all getters before
``clean_all()``.

Validation
----------

Run one full-sky catalog through every active shear engine with:

.. code-block:: bash

   python3 tests/python/shear_corr_all_engines.py \
      --geometry sphere --fits shear.fits --engine all --statistics both \
      --min-sep 1 --max-sep 100 --sep-units arcmin --max-points 4096 \
      --outdir Output_shear_comparison

The driver reads the FITS map once, retains unit vectors and local east/north
spin-2 components, and passes the same arrays to all three native engines. See
``tests/python/README_shear_corr_all_engines.md`` for the Takahashi example, mask
handling, spin convention, output arrays, dual-node acceptance controls,
exact-validation overrides, and the ``timing_report.txt``
wall/process-CPU report.

With
``no-smooth-pivot``, the native ``only-2pcf`` path uses a deterministic
symmetric dual-octree traversal. Smoothed runs use representative body pivots.
Spherical cells carry basis-aware shear moments transported to the cell-center
tangent frame, and are accepted only when radial-bin, angular-phase, and
transport-error bounds all pass. Add ``--exact-tree`` to select the body-level
reference calculation. The driver reports pairwise relative differences between matching
native output definitions.

Run the independent reference calculation and OpenMP determinism test with:

.. code-block:: bash

   tests/make_tests/run_test_shear_sphere_2balls_octree_omp
   tests/make_tests/run_test_shear_sphere_2balls_kdtree_omp
   tests/make_tests/run_test_shear_sphere_2balls_balltree_omp

The derivation implemented by the add-on follows the ring-multipole and
finite-window equations in Lucas-Porth et al., *The three-point correlation
function of cosmic shear* (arXiv:2309.08601). The paper is a scientific
reference; cTreeBalls runtime and build behavior are defined by this source
tree.


Optional count/window multipoles
--------------------------------

Shear 3PCF retains its full window and edge-corrected output by default.
Add ``options=no-count-multipoles`` (Python driver:
``--more-options no-count-multipoles``) to compute only the triplet-weight
monopole N0. The output then contains Upsilon/N0, which remains convolved with
the survey window; it is not a finite-window solution. This mode skips higher
weight-ring moments, their radial-pair contractions and the coupling solve.
It also reduces ring workspace and the phase order needed by tree acceptance.
Consequently approximate tree results can differ beyond the final normalization;
use the exact body limit to isolate normalization differences.

An explicit ``edge-corrections`` option always takes precedence and restores
the full window and solve, even alongside ``no-count-multipoles``. Shear does
not require the scalar-only ``no-normalize-HistZeta`` companion option.
Two-point statistics are unchanged. Scalar production two-ball engines already
omit higher count modes without ``edge-corrections`` and always retain their
normalization denominator; this new token adds no further scalar optimization.

``getShearWindowMultipoles()`` returns shape (1,B,B) in N0-only mode and the
usual (4*nmax+1,B,B) otherwise. ``getShearWindowMultipoleOrders()`` returns the
actual window axis. Existing Gamma getters retain their names but return the
selected normalization; ``getResults()`` uses ``gamma_normalized`` instead of
``gamma_corrected`` when correction is off. Native run information, text headers,
driver JSON and plots explicitly identify the selected normalization. Skipped
window orders are omitted, never represented as measured zeros. No matrix
condition/residual diagnostics are published when no solve was performed.

Measured approximation penalty
------------------------------

The comparison drivers accept ``--accuracy-reference EXACT_OUTPUT_ROOT`` and
``--accuracy-report-only``. They print relative L2 percentages independently
for raw moments, windows, corrected results and all four shear components,
and save ``ENGINE/qualification/approximation-penalty.json``. Without an exact
reference they report ``UNMEASURED``. The 5% comparison target does not alter
the selected controls; report-only preserves results that exceed it. See
:download:`the calibration and settings guide <APPROXIMATION_DEFAULTS.md>`.
