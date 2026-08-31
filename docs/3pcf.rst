
3-point correlation functions
=============================

**cBalls** can also compute 3-point correlations (3pcf) when counts (N) or a scalar field are involved (like convergence :math:`\kappa`). In particular to compute the correlation KKK the numerical code use the harmonic base.

Harmonic base
-------------

Scalar angular contract
~~~~~~~~~~~~~~~~~~~~~~~

The scalar Fourier engines (GGG, KD-tree, balltree, BALLS4, and the two-ball
engines, including their MPI variants) use observer-relative tangent angles
in three dimensions. Positions must retain their original observer origin.
Both in-memory catalogs and the ASCII readers preserve this frame; a masked
sky catalog must not be translated to its bounding-box center.

Radial bins use Euclidean separation, which is chord distance for unit-sphere
catalogs. At pivot :math:`i`, project both neighbor displacements onto the
plane perpendicular to :math:`\hat{x}_i`. The oriented angle is

.. math::

   \phi_{ijk} = \operatorname{atan2}
   \left[\hat{x}_i\cdot(\hat{t}_{ij}\times\hat{t}_{ik}),
   \hat{t}_{ij}\cdot\hat{t}_{ik}\right].

With ``options=no-normalize-HistZeta,weights-norm``, raw mode :math:`m` is

.. math::

   Z_m(a,b) = \sum_{i,j,k\;\mathrm{distinct}}
      w_i\kappa_i w_j\kappa_j w_k\kappa_k
      \mathbf{1}_a(r_{ij})\mathbf{1}_b(r_{ik})\exp(i m\phi_{ijk}).

The two legs are ordered. Repeated-neighbor terms (:math:`j=k`) are removed
inside the native search, including when a cell represents several neighbors.
No post-processing subtraction or division by the catalog size is needed.
The four C/Cython component matrices reconstruct this mode as
``coscos + sinsin + 1j*(sincos - cossin)``; Cython indexes order :math:`m`
with argument ``m+1``. This matches TreeCorr's positive LogMultipole modes.
Individual component matrices depend on the chosen tangent basis; the
reconstructed complex modes do not.

A pivot at the observer, or a coincident/radial/antipodal neighbor leg, has no
defined bearing and does not contribute to angular 3PCF. Ordinary 2PCF pair
counts still include such pairs when their distances fall in a radial bin.
The NDIM=2 engines retain their planar angular convention. The separate
``octree-3pcf-3d-omp``/MPI spherical-harmonic survey estimator is not this
angular Fourier estimator.

Smoothing and approximations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``SMOOTHPIVOTON=1`` only compiles support. Smoothing is inactive unless
``options=smooth-pivot`` is explicitly requested; setting ``rsmooth`` alone
does not enable it. Raw distinct-multipole mode in KD-tree and legacy balltree
rejects runtime smoothing because grouped pivots do not implement this raw
weighted estimator. Other engines retain their documented option restrictions.

Cell aggregation is approximate. Its angular acceptance tests bound projected
bearings, and inherited two-ball moments are transported between pivot tangent
bases. Check production accuracy against ``no-one-ball,no-two-balls`` and
decreasing ``theta`` for the desired angular range and maximum multipole.

Rebuild both C and Cython after updating. Higher modes and raw KD-tree/balltree
or BALLS4 outputs from older implementations are not numerically interchangeable
with this contract and should be recomputed.

Regression checks
~~~~~~~~~~~~~~~~~

From the source root, with matching C and Cython builds::

   python3 -m pytest -q tests/make_tests/test_scalar_numerical_contract.py
   CBALLS_TEST_MPI=1 mpiexec -n 2 python3 tests/make_tests/test_scalar_numerical_contract.py

The tests compare independent ordered-triplet sums, TreeCorr, full/partial-sky
rotations, nonuniform weights, degenerate bearings, inactive smoothing,
OpenMP threads, MPI ranks, and both ASCII input formats. CI runs both compiled
smoothing profiles. These small-catalog checks validate numerical contracts;
they are not a full-sky performance certification.
