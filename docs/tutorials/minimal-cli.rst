Minimal Command-Line Tutorial
=============================

This workflow creates a small synthetic unit-sphere catalog and computes the
statistics enabled by the current build.

1. Build cTreeBalls
-------------------

From the repository root:

.. code-block:: bash

   make clean
   make PYTHON=python3 all

2. Run a Reduced Calculation
----------------------------

.. code-block:: bash

   ./cballs nbody=4096 sizeHistN=12 mChebyshev=3 \
      rootDir=Output_tutorial numberThreads=1 \
      verbose=1 verbose_log=1

The generated catalog avoids external downloads.  The low body count and
histogram resolution keep this tutorial compact.

3. Inspect Provenance
---------------------

.. code-block:: bash

   cat Output_tutorial/parameters_null-cballs-usedvalues

Confirm ``nbody``, ``sizeHistN``, ``mChebyshev``, ``rootDir``, and
``numberThreads`` before interpreting the histograms.

4. Inspect the 2PCF
-------------------

.. code-block:: bash

   head Output_tutorial/rBins.txt
   head Output_tutorial/histXi2pcf.txt

The first file contains radial coordinates and the second contains the 2PCF.
Units follow the generated unit-sphere geometry.

5. Identify 3PCF Products
-------------------------

.. code-block:: bash

   find Output_tutorial -maxdepth 1 -type f -name 'histZeta*' | sort

3PCF filenames depend on the compiled engine and requested options.  Each
multipole matrix corresponds to the radial grid recorded by ``rBins.txt``.

6. Scale Up Carefully
---------------------

Increase one of ``nbody``, ``sizeHistN``, ``mChebyshev``, or
``numberThreads`` at a time.  Follow :doc:`../performance` before selecting
science settings.
