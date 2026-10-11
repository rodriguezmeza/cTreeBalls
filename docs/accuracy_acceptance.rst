Bounded numerical acceptance
============================

Generated from ``scripts/accuracy_cases.py`` by ``generate_accuracy_docs.py``.
The same definitions drive the numerical acceptance runner. These are
retained regression settings, not the fast runtime defaults or a universal
error guarantee. Every listed route is a required acceptance case.

.. list-table:: Retained acceptance cases
   :header-rows: 1
   :widths: 46 10 14 30

   * - Engine / route
     - Theta
     - Relative L2 limit
     - Explicit smoothing case
   * - kdtree-2balls-omp
     - 0.5
     - 0.02
     - 0.0003 Cartesian units
   * - balltree-2balls-omp
     - 0.2
     - 0.02
     - 0.0003 Cartesian units
   * - octree-2balls-omp
     - 0.025
     - 0.02
     - unsupported
   * - octree-sincos-omp (core estimator)
     - 0.025
     - 0.02
     - 0.0003 arcmin
   * - octree-2balls-omp (legacy-one-ball)
     - 0.025
     - 0.02
     - 0.0003 arcmin
   * - octree-shear-sphere-2balls-omp
     - 0.05
     - 0.02
     - 1 arcmin
   * - kdtree-shear-sphere-2balls-omp
     - 0.05
     - 0.02
     - 1 arcmin
   * - balltree-shear-sphere-2balls-omp
     - 0.05
     - 0.02
     - 1 arcmin

The fixtures have 512 weighted points, six logarithmic chord bins from
0.02 to 1.8, orders 0 through 5, and two OpenMP threads. Each pair and
3PCF product must satisfy L2(error) <= 1e-10 + rtol * L2(reference).
Smoothing is explicitly enabled only for the clustered fixture.

Modern scalar products are raw distinct-neighbor Fourier moments;
the core estimator retains its separate normalized neighbor-mean definition.
Shear acceptance measures raw Upsilon, not corrected Gamma. Window
conditioning, arbitrary data, different opening/smoothing settings and
experimental reuse need separate qualification. The default theta is one;
successful execution alone remains UNMEASURED.
