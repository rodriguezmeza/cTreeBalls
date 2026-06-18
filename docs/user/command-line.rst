Command-Line Usage
==================

The executable accepts ``name=value`` tokens.  Do not put spaces around ``=``
on the command line.

.. code-block:: bash

   ./cballs searchMethod=octree-ggg-omp nbody=4096 rootDir=Output

Syntax
------

Run entirely from command-line parameters:

.. code-block:: text

   ./cballs name=value name=value ...

or pass a parameter file as the first argument:

.. code-block:: text

   ./cballs path/to/parameters_file

Parameter files use ``name = value`` and may contain explanatory comments.
Use full names in files even when the CLI offers a short alias.

Discovery Commands
------------------

.. code-block:: bash

   ./cballs --help
   ./cballs --clue

The available defaults depend on ``Makefile_settings`` and enabled add-ons, so
the help output from the compiled executable takes precedence over copied
examples.

Common Parameters
-----------------

.. list-table::
   :header-rows: 1
   :widths: 23 17 60

   * - Parameter
     - Typical default
     - Purpose
   * - ``searchMethod``
     - ``octree-sincos-omp`` for the unit-sphere defaults
     - Select the tree, octree, k-d-tree, or neighbor-box implementation.
   * - ``infile``
     - empty
     - Input catalog or comma-separated catalogs.
   * - ``infileformat``
     - ``columns-ascii``
     - Input format or comma-separated format list.
   * - ``iCatalogs``
     - ``1``
     - Catalog labels used for auto- or cross-correlations.
   * - ``rootDir``
     - ``Output``
     - Directory for histograms, parameters, and logs.
   * - ``sizeHistN``
     - ``20`` in the unit-sphere defaults
     - Number of radial bins.
   * - ``rminHist`` / ``rangeN``
     - build dependent
     - Minimum and maximum radial separations.
   * - ``mChebyshev``
     - ``7``
     - Highest 3PCF multipole index; ``mChebyshev + 1`` multipoles are used.
   * - ``numberThreads``
     - ``16`` in OpenMP builds
     - Requested OpenMP thread count.
   * - ``verbose`` / ``verbose_log``
     - ``0``
     - Console and log verbosity.
   * - ``options``
     - empty
     - Comma-separated behavior flags such as ``compute-HistN`` and ``and-CF``.

The complete historical parameter descriptions remain in :doc:`../params` and
``tests/In/parameters_explained``.

Reproducible Runs
-----------------

For production work, save the command or parameter file together with the
generated ``*-usedvalues`` record, the commit hash, compiler version, active
Makefile settings, and input-catalog metadata.
