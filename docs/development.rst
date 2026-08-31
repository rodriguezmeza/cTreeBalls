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

   make -j4 cballs cyballs-static-lib
   CBALLS_STATIC_LIBRARY_READY=1 python3 setup.py build_ext --inplace --force
   ./cballs nbody=4096 sizeHistN=12 mChebyshev=3 \
      rootDir=Output_check numberThreads=1 verbose=0 verbose_log=0
   python3 -c "from cyballs import cballs; print(cballs)"
   python3 -m pytest -q tests/make_tests/test_scalar_numerical_contract.py
   python3 -m sphinx -E -a -n -W --keep-going -b html docs docs/_build/html

Use the same build profile in both commands. Documentation-only edits require
only the Sphinx command, not a C rebuild. The API pages are maintained without
importing a compiled extension.

The maintained launchers and their Python programs are in ``tests/make_tests``.
``tests/scripts`` retains older catalog-specific workflows. See
:doc:`benchmarks` and ``tests/make_tests/README.md`` for numerical contracts.
Scalar CI checks both smoothing profiles, MPI ranks, and TreeCorr references;
the main test workflow also covers C/Cython, packaging, and sanitizer jobs.
The documentation workflow builds HTML with warnings treated as errors.

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

Shared Numerical Helpers
------------------------

``include/angular_contracts.h`` owns observer-frame selection, tangent phases,
and projected angular acceptance. Readers, startup, trees, and searches must
agree on this frame. Two-ball inherited moments must be transported when the
pivot tangent basis changes.

``include/scalar_moments.h`` implements raw legacy KD/balltree moments and
neighbor self-term subtraction. Aggregates need a sum of squared weighted
fields, not the square of their sum. Document estimator changes in :doc:`3pcf`
and update independent tests instead of merely copying production formulas.
