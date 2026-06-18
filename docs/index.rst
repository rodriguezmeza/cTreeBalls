cTreeBalls
==========

**cTreeBalls** measures two- and three-point correlation functions from point
catalogs and scalar fields using tree- and ball-based search methods.  The
project provides:

* the compiled C command-line executable ``cballs``;
* the static library ``libcballs.a``;
* the Cython extension ``cyballs``;
* test catalogs, plotting scripts, and benchmark workflows.

The guide is organized like the companion 3ptWL projects: begin with the
overview, installation, and quickstart; use the task-oriented pages for real
runs; then consult the tutorials and developer reference.

Basic Usage
-----------

Build from a source checkout:

.. code-block:: bash

   git clone https://github.com/rodriguezmeza/cTreeBalls.git
   cd cTreeBalls
   make clean
   make PYTHON=python3 all

Run a compact synthetic-catalog calculation:

.. code-block:: bash

   ./cballs nbody=4096 sizeHistN=12 mChebyshev=3 \
      rootDir=Output_quick numberThreads=1 verbose=0 verbose_log=0

Or call the same C core through Python:

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

Read :doc:`installation` and :doc:`quickstart` first.  For production runs,
consult :doc:`user/command-line`, :doc:`user/inputs`, :doc:`user/outputs`, and
:doc:`performance`.  Existing detailed material on parameters, formats,
pre/post-processing, 2PCF, 3PCF, and add-ons remains available under Tutorials
and Reference.

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
   performance

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
   addons
   troubleshooting
   development
   citing
