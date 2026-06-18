Quickstart
==========

This page validates a source build with a small synthetic catalog.  It avoids
the large external Takahashi data files used by the full tutorial.

Build From a Source Checkout
----------------------------

.. code-block:: bash

   git clone https://github.com/rodriguezmeza/cTreeBalls.git
   cd cTreeBalls
   python3 -m pip install --user numpy Cython
   make clean
   make PYTHON=python3 all

Inspect Runtime Help
--------------------

.. code-block:: bash

   ./cballs --help
   ./cballs --clue

``--help`` lists compiled defaults.  ``--clue`` prints a compact reminder of
the command-line syntax.

Run a Compact CLI Calculation
-----------------------------

.. code-block:: bash

   ./cballs nbody=4096 sizeHistN=12 mChebyshev=3 \
      rootDir=Output_quick numberThreads=1 verbose=0 verbose_log=0

With the repository defaults, cTreeBalls generates a random scalar field on
the unit sphere and evaluates enabled two- and three-point statistics.  The
reduced ``nbody`` and histogram sizes are intended for a smoke test, not a
science analysis.

Inspect the Results
-------------------

.. code-block:: bash

   find Output_quick -maxdepth 2 -type f | sort
   head Output_quick/histXi2pcf.txt
   cat Output_quick/parameters_null-cballs-usedvalues

Output names can vary with compile-time add-ons and runtime options.  The
used-values file is the authoritative record of the parameters applied.

Run With a Parameter File
-------------------------

The repository includes an explained template:

.. code-block:: bash

   less tests/In/parameters_explained
   cd tests
   ../cballs ./In/parameters_explained

That example expects a catalog named in the parameter file.  Update ``infile``
before running it, or continue with the synthetic example above.

Verify the Python Extension
---------------------------

.. code-block:: bash

   python3 -c "from cyballs import cballs; print(cballs)"

For a full wrapper calculation, continue with
:doc:`tutorials/python-wrapper`.

Next Steps
----------

* :doc:`user/inputs` explains catalog formats and multiple-catalog runs.
* :doc:`user/command-line` describes parameter files and common controls.
* :doc:`user/outputs` identifies histogram and provenance outputs.
* :doc:`performance` explains search methods, bins, and OpenMP settings.
* :doc:`tutorials/index` contains longer CLI, Python, 2PCF, and 3PCF workflows.
