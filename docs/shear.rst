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

The octree two-ball addon keeps its native dual-node behavior by default. Add
``options=legacy-one-ball`` to use its privately linked one-node compatibility
kernel. This is a complete compatibility dispatch, including combined 2PCF/3PCF traversal,
spherical transport, smoothing, masks, edge correction, normalization,
output, and deterministic ``BALLS4SCANLEV`` reduction. The compatibility
kernel is linked privately when the ordinary addon is disabled.

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

With ``SMOOTHPIVOTON=1``, smoothing is enabled by default. Each claimed body's
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
acceptance contracts to an independently built median KD tree. Its 2PCF uses
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
contract. ``nsmooth`` sets the ball-tree leaf capacity.

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
oracle. The driver reports pairwise relative differences between matching
native output contracts.

Run the independent oracle and OpenMP determinism test with:

.. code-block:: bash

   tests/make_tests/run_test_shear_sphere_2balls_octree_omp
   tests/make_tests/run_test_shear_sphere_2balls_kdtree_omp
   tests/make_tests/run_test_shear_sphere_2balls_balltree_omp

The derivation implemented by the add-on follows the ring-multipole and
finite-window equations in Lucas-Porth et al., *The three-point correlation
function of cosmic shear* (arXiv:2309.08601). The paper is a scientific
reference; cTreeBalls runtime and build behavior are defined by this source
tree.
