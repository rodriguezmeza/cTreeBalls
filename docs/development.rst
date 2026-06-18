Developer Guide
===============

Repository Layout
-----------------

.. list-table::
   :header-rows: 1
   :widths: 28 72

   * - Path
     - Purpose
   * - ``main/``
     - Command-line entry point.
   * - ``source/``
     - Startup, I/O, tree construction, searches, histograms, and cleanup.
   * - ``include/``
     - Parameter structures, run state, data structures, and prototypes.
   * - ``addons/``
     - Search engines, catalog formats, parser/Cython hooks, and bundled libraries.
   * - ``python/``
     - ``cyballs`` Cython source and declarations.
   * - ``tests/``
     - Parameter files, catalogs, scripts, references, and benchmarks.
   * - ``docs/``
     - Sphinx sources and generated documentation.

For the detailed runtime flow and module ownership map, see
:doc:`code_structure`.

Validation Checklist
--------------------

Before proposing a code or documentation change:

.. code-block:: bash

   make clean
   make PYTHON=python3 all
   ./cballs nbody=4096 sizeHistN=12 mChebyshev=3 \
      rootDir=Output_check numberThreads=1 verbose=0 verbose_log=0
   python3 -c "from cyballs import cballs; print(cballs)"
   python3 -m sphinx -E -a -W --keep-going -b html docs docs/_build/html

The fuller regression workflow is described in ``tests/Readme_first.txt`` and
the scripts under ``tests/scripts``.

Adding Parameters
-----------------

Coordinate changes across:

* ``include/cmdline_defs.h`` for defaults, help text, and aliases;
* ``include/cmdline_data.h`` for persistent command state;
* parameter reading, checking, and printing in ``source/startrun.c``;
* parser/add-on declarations used by ``cyballs``;
* documentation and compact validation examples.

Every new user parameter must appear in the generated used-values file.

Adding Formats or Search Methods
--------------------------------

Input/output format dispatch is implemented in ``source/cballsio.c`` and
extended by add-on sockets.  Search method names are mapped during startup and
dispatched by ``EvalHist`` in ``source/cballs.c``.  Keep optional modules under
``addons/`` and document their compile-time switch, runtime name, supported
geometry, outputs, and cleanup ownership.

Documentation Style
-------------------

* begin with the user task or scientific concept;
* state the working directory for commands;
* use current names such as ``cyballs``;
* distinguish smoke-test settings from science settings;
* preserve old workflows when reorganizing navigation;
* build Sphinx with warnings treated as errors.
