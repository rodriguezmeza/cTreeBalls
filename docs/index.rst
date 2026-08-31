cTreeBalls
==========

**cTreeBalls** measures two- and three-point correlation functions from point
catalogs, scalar fields, shear, and Lyman-alpha forests using tree- and
ball-based search methods. Supported estimators depend on the build. The
project provides:

* the compiled C command-line executable ``cballs``;
* the static library ``libcballs.a``;
* the Cython extension ``cyballs``;
* OpenMP and MPI addon engines, test catalogs, plotting scripts, and benchmarks.

This guide describes the 1.1.0 source tree, including the repaired scalar
numerical contract. Published packages may lag behind local changes.
Begin with :doc:`search_methods` to select the correct field and geometry.

Basic Usage
-----------

Install a published Python release from PyPI:

.. code-block:: bash

   python3 -m pip install cTreeBalls

Then import the compiled wrapper as ``cyballs``:

.. code-block:: python

   from cyballs import cballs

For a no-checkout notebook workflow, open the standalone Colab example:
`cTreeBalls minimal Colab notebook <https://colab.research.google.com/github/rodriguezmeza/cTreeBalls/blob/main/examples/cTreeBalls_minimal_colab.ipynb>`_.

Use a source checkout for the ``cballs`` executable, current unpublished
changes, test catalogs, or development files:

.. code-block:: bash

   git clone https://github.com/rodriguezmeza/cTreeBalls.git
   cd cTreeBalls
   python3 -m pip install numpy Cython setuptools wheel
   make -j4 cballs cyballs-static-lib
   CBALLS_STATIC_LIBRARY_READY=1 python3 setup.py build_ext --inplace --force

Configure one matching profile first; see :doc:`installation`.

Run a compact synthetic-catalog calculation:

.. code-block:: bash

   ./cballs nbody=4096 sizeHistN=12 mChebyshev=3 \
      rootDir=Output_quick numberThreads=1 verbose=0 verbose_log=0

Call the same C core through Python:

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
   model.Run()
   radius = model.getrBins()
   xi = model.getHistXi2pcf()
   model.clean_all()

How to Use This Guide
---------------------

Read :doc:`installation` first for the PyPI and source-build paths.  Then use
:doc:`quickstart` or the :doc:`tutorials/colab-python-wrapper` notebook for a
small smoke test.  For production runs, consult :doc:`user/command-line`,
:doc:`user/inputs`, :doc:`user/outputs`, and :doc:`performance`.  Existing
detailed material on parameters, formats, pre/post-processing, 2PCF, 3PCF, and
add-ons remains available under Tutorials and Reference.

.. toctree::
   :maxdepth: 2
   :caption: User Guide

   overview
   installation
   quickstart
   user/command-line
   user/configuration
   user/inputs
   user/outputs
   user/python
   search_methods
   shear
   lyman_alpha
   scalar_3d
   performance
   benchmarks

.. toctree::
   :maxdepth: 2
   :caption: Tutorials

   tutorials/index

.. toctree::
   :maxdepth: 2
   :caption: Reference

   params
   catalog_files
   api
   python
   code_structure
   build_profiles
   addons
   troubleshooting
   development
   citing
