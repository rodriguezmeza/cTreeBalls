Weak-Lensing Shear
==================

The optional ``octree-shear-omp`` add-on measures flat-sky shear 2PCFs and the
four natural shear 3PCF components in the Porth x-projection,
:math:`\Gamma^\times`. Enable it when building C and Cython:

.. code-block:: bash

   make OCTREESHEAROMPON=1 clean
   make OCTREESHEAROMPON=1 PYTHON=python3 all

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

Geometry and Binning
--------------------

The estimator supports a native 2D build or a 3D build whose third coordinate
is constant. It rejects periodic geometry, non-flat catalogs, and selected
catalogs lying on different constant-coordinate planes. ``rminHist``,
``rangeN``, ``sizeHistN``, and ``useLogHist`` control radial bins;
``sizeHistPhi`` controls the angular reconstruction and must satisfy the core
parameter checks.

Catalog selection supports tomography. One input catalog computes
``Z1,Z1,Z1``; two compute ``Z1,Z2,Z2``; three or more compute ``Z1,Z2,Z3``
from the first three entries of ``iCatalogs``. The 2PCF uses ``Z1,Z2``. The
core parser requires one ``iCatalogs`` entry per input catalog even though this
estimator ignores entries after the third. A single common translation is
applied before the trees are built, so cross-catalog separations are preserved.

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
       "searchMethod": "octree-shear-omp",
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

Run the independent oracle and OpenMP determinism test with:

.. code-block:: bash

   tests/make_tests/run_test_shear_octree_omp

The derivation implemented by the add-on follows the ring-multipole and
finite-window equations in Lucas-Porth et al., *The three-point correlation
function of cosmic shear* (arXiv:2309.08601). The paper is a scientific
reference; cTreeBalls runtime and build behavior are defined by this source
tree.
