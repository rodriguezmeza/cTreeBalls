Python Interface Reference
==========================

The current extension module is ``cyballs`` and the wrapped class is
``cballs``.  It calls the same C lifecycle as the command-line executable.

Basic Lifecycle
---------------

.. code-block:: python

   from cyballs import cballs

   model = cballs()
   model.set(
       nbody=4096,
       sizeHistN=12,
       mChebyshev=3,
       rootDir="Output_python",
       numberThreads=1,
       verbose=0,
       verbose_log=0,
   )
   cpu_time = model.Run()
   radius = model.getrBins()
   xi = model.getHistXi2pcf()
   model.clean_all()

``set`` accepts either one mapping or keyword arguments.  ``Run`` expands the
requested stage dependencies, converts the Python parameters into the C parser
representation, initializes the search, and leaves result arrays available to
the getters.  ``clean_all`` releases C allocations and clears the parameter
dictionary.

Run Levels
----------

``Run(level=[...])`` accepts the last C stage required by the caller.  Its
default, ``["MainLoop"]``, executes the full calculation while keeping arrays
available to Python.  The common order is:

1. ``input`` parses the staged Python parameters.
2. ``StartRun_Common`` validates settings and allocates run state.
3. ``PrintParameterFile`` records the used values.
4. ``SetNumberThreads`` applies OpenMP thread control.
5. ``MainLoop`` builds trees, performs searches, and evaluates histograms.

Available Getters
-----------------

``getrBins``
    Radial-bin values.

``getHistNN`` and ``getHistCF``
    Neighbor-count and count-derived correlation histograms.

``getHistXi2pcf``
    Two-point correlation function.

``getHistZetaMsincos(m, type)``
    3PCF multipole matrix.  ``type`` selects cosine, sine, sine-cosine, or
    cosine-sine components.

``getHistZetaM_EE(m)``
    Edge-corrected 3PCF multipole matrix when that build path is enabled.

``getNThreads``, ``getnMultipoles``, ``getTheta``, ``getrsmooth``,
``getCPUTime``, ``getVersion``, ``getRootDir``, and ``getNBody`` expose
selected scalar state.

Memory Ownership
----------------

The C core owns catalogs, trees, and histogram buffers after initialization.
Call ``clean_all`` after copying required NumPy results.  For repeated runs in a
long-lived process, do not retain references to arrays that depend on C memory
after cleanup.

See Also
--------

See :doc:`user/python` for a task-oriented workflow and
:doc:`code_structure` for Cython implementation details.

