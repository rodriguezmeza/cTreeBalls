Troubleshooting and Common Errors
=================================

``python: command not found``
------------------------------

The Makefiles default to ``PYTHON ?= python``.  Build with:

.. code-block:: bash

   make clean
   make PYTHON=python3 all

``No module named Cython`` or ``numpy``
-----------------------------------------

Install the wrapper build dependencies:

.. code-block:: bash

   python3 -m pip install --user numpy Cython

``No module named cyballs``
---------------------------

Confirm the extension was built for the Python interpreter being used:

.. code-block:: bash

   make clean
   make PYTHON=python3 all
   python3 -c "from cyballs import cballs; print(cballs)"

Older notes may mention ``cballys``; the current source module is ``cyballs``.

OpenMP Link or Runtime Errors
-----------------------------

The default GNU build uses ``-fopenmp`` and links ``-lgomp``.  Verify that
``CC`` and ``OMPFLAG`` in ``Makefile_machine`` match the installed compiler.
If OpenMP is unavailable, set ``OPENMPMACHINE = 0``, rebuild cleanly, and omit
``numberThreads`` expectations from benchmarks.

GSL or CFITSIO Headers/Libraries Not Found
------------------------------------------

The repository defaults use bundled GSL and CFITSIO paths.  If switching to
system libraries, set ``GSLINTERNAL = 0`` or ``CFITSIOLIBON = 0`` and update
the include/library paths in ``Makefile_machine``.  Always run ``make clean``
after changing these switches.

Parameter Was Not Read
----------------------

Command-line tokens must use ``name=value`` with no spaces.  Parameter files
should use full parameter names.  Check the compiled ``./cballs --help``
because enabled add-ons determine which names are available.

Empty or Unexpected Histograms
------------------------------

Check:

* input format and column order;
* catalog units versus ``rminHist`` and ``rangeN``;
* ``iCatalogs`` selection;
* search method compatibility with the geometry;
* requested options and compile-time ``TWOPCFON``/``TPCFON`` settings;
* the generated used-values file.

Segmentation Fault in a Python Loop
-----------------------------------

Copy NumPy results before cleanup and call ``clean_all`` once per completed
run.  If repeated in-process runs remain unstable, isolate calculations in
separate processes and report the minimal parameter set, platform, compiler,
and commit hash to the issue tracker.
