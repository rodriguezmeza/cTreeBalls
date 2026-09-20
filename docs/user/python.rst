Python and In-memory Catalogs
=============================

Import ``cballs`` from ``cyballs`` built for the current interpreter and
feature profile. See :doc:`../installation` and :doc:`../api`.

Scalar Example
--------------

Run from the checkout root with the GGG, 2PCF, and 3PCF features enabled:

.. code-block:: python

   import numpy as np
   from cyballs import cballs

   rng = np.random.default_rng(123)
   xyz = rng.normal(size=(128, 3))
   xyz /= np.linalg.norm(xyz, axis=1)[:, None]
   kappa = rng.normal(size=len(xyz))
   weights = rng.uniform(0.5, 1.5, len(xyz))

   model = cballs()
   try:
       model.set(searchMethod="octree-2balls-omp", rootDir="Output_memory",
                 rminHist=0.05, rangeN=1.0, sizeHistN=6, mChebyshev=3,
                 useLogHist=True, usePeriodic=False, numberThreads=2,
                 options="no-normalize-HistZeta,weights-norm,no-one-ball")
       model.set_catalog(xyz, kappa=kappa, weights=weights)
       model.Run()
       radius = model.getrBins()
       xi = model.getHistXi2pcf()
       cc, ss, sc, cs = [model.getHistZetaMsincos(2, t) for t in range(1, 5)]
       zeta_m1 = cc + ss + 1j*(sc-cs)
   finally:
       model.clean_all()

``m=1`` in a getter denotes the monopole; argument ``2`` above denotes
physical order one. Getters return NumPy copies but require a live, completed
compatible run. See :doc:`../3pcf` for normalization and angular conventions.

Registration and Ownership
--------------------------

``set_catalog`` accepts ``positions`` shaped ``(N, compiled_ndim)`` and optional
one-dimensional ``kappa``, ``weights``, ``mask``, ``gamma1`` and ``gamma2``.
At least three bodies and finite values are required. Omitted scalar values and
weights default to one. Supply both shear components or neither.
A mask must be boolean or 0/1; select ``read-mask`` on a mask-capable engine
to apply it. Supplying a mask array does not itself request edge correction.

Registration retains prepared contiguous NumPy arrays. They are copied into
C-owned body storage when a run starts, and the input is not modified by C.
This avoids rereading catalogs but is not a zero-copy native API. Do not mutate
registered arrays between runs unless the changed catalog is intentional.

``struct_cleanup()``
    Free active C run state, keeping staged parameters and registered catalogs.

``clear_catalogs()``
    Free active C state and discard registered catalog references.

``clean_all()``
    Free C state and discard both staged parameters and registered catalogs.

To reuse one catalog, call ``struct_cleanup`` between methods, change parameters,
then call ``Run`` again. Use a different ``rootDir`` for each engine.
``set_catalog(..., catalog=0)`` is zero-based; multiple catalogs must use
contiguous slots. C ``iCatalogs="1,2"`` is one-based.

All-engines Drivers
-------------------

``tests/python/kappa_corr_all_engines.py`` reads scalar FITS/NPZ input once,
filters active native engines by mask/edge capabilities, and retains the
catalog across runs. For example::

   python3 tests/python/kappa_corr_all_engines.py --list-engines
   python3 tests/python/kappa_corr_all_engines.py --fits map.fits --mask mask.fits \
       --engine all-omp --edge-corrections --threads 4 --outdir Output_kappa
   python3 tests/python/kappa_corr_all_engines.py --fits map.fits \
       --engine octree-2balls-omp,kdtree-2balls-omp \
       --threads 4 --outdir Output_compare

The driver finds the executable independently of the current directory. An
explicit ``--cballs`` argument has priority, followed by
``CTREEBALLS_CBALLS``/``CBALLS``, the source-tree executable, and ``cballs`` on
``PATH``. This avoids hard-coded ``DEFAULT_CBALLS`` edits in external scripts.

Normal plots include 2PCF curves and one radial-bin heat map per complex 3PCF
mode. Flattened plots place each radial-bin matrix on one row-major bin index,
with one panel per mode; use ``--no-flatten-plots`` or ``--no-plots`` to disable
them.

These are angular scalar products, not shear or physical-3D survey products.
Read the :download:`kappa driver README <../../tests/python/README_kappa_corr_all_engines.md>`
for current angle controls, normalization, MPI, and output formats.
For large native-resolution HEALPix comparisons, add ``--statistics 2pcf`` to
avoid hidden 3PCF work and ``--max-points N`` for deterministic bottom-hash
sampling. Map and mask columns are scanned in bounded chunks, and masked-out
pixels are discarded before the Cartesian catalog is built. A thinned run
measures the sampled catalog; use ``--max-points 0`` for the full retained map.
All compared arrays use the native cTreeBalls normalization contract.
``no-normalize-HistZeta,weights-norm`` selects weighted raw sums. With
``--edge-corrections``, active scalar engines build window modes through order
``2M`` and apply the per-bin scalar mode-coupling solve.

The all-engine drivers write ``timing_report.txt`` and structured
``summary.json["timings"]`` entries with setup, compute, and total wall/process
CPU seconds. Process CPU can exceed wall time under OpenMP. MPI CPU is summed
across participating ranks, and wall time is their maximum; per-rank measurements
are also retained. Check each driver's timing scope before comparing results.

For forests, use ``set_forest_catalog(positions, delta, weights, forest_ids)``
or ``tests/python/lya_corr_all_engines.py``. Forest IDs must remain integer arrays.
See :doc:`../lyman_alpha`. For shear and physical-3D catalogs see
:doc:`../shear` and :doc:`../scalar_3d`.

MPI and Notebooks
-----------------

Import ``mpi4py.MPI`` before starting MPI runs. All ranks must use the same
engine, parameters, complete catalogs, and lifecycle sequence. The supplied
drivers perform broadcasting; direct registration does not. Read native
published histograms only on rank 0.

Use ``examples/cyballs_in_memory_catalog.py`` and the corresponding notebook
for a runnable catalog example. Restart the kernel after rebuilding cyballs:
a Python process does not reload an already loaded native library in place.
