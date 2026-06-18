Python Wrapper
==============

The compiled extension is imported as ``cyballs`` and exposes the ``cballs``
class.

Minimal Run
-----------

.. code-block:: python

   from cyballs import cballs

   model = cballs()
   model.set({
       "nbody": 4096,
       "sizeHistN": 12,
       "mChebyshev": 3,
       "rootDir": "Output_python",
       "numberThreads": 1,
       "verbose": 0,
       "verbose_log": 0,
   })

   cpu_time = model.Run()
   radius = model.getrBins().copy()
   xi = model.getHistXi2pcf().copy()
   model.clean_all()

Parameters
----------

``set`` accepts a mapping or keyword arguments:

.. code-block:: python

   model.set({"searchMethod": "octree-ggg-omp"})
   model.set(rootDir="Output_python", numberThreads=4)

The wrapper rejects names that the compiled C parser does not consume.  This
helps expose spelling errors and mismatches between Python code and active
build options.

Catalog Workflow
----------------

.. code-block:: python

   model = cballs()
   model.set(
       infile="map.fits",
       infileformat="fits-healpix",
       iCatalogs="1",
       rminHist=0.00213811,
       rangeN=0.0633205,
       sizeHistN=20,
       options="compute-HistN",
       rootDir="Output_map",
       numberThreads=8,
       verbose=0,
       verbose_log=0,
   )
   model.Run()

   radius = model.getrBins().copy()
   xi = model.getHistXi2pcf().copy()
   counts = model.getHistNN().copy()
   model.clean_all()

The ``fits-healpix`` format requires the corresponding CFITSIO/I/O add-ons in
the compiled library.

Cleanup
-------

Copy required NumPy arrays before ``clean_all`` and call cleanup after every
completed or failed workflow.  See :doc:`../python` for the method reference
and available getters.
