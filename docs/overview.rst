Overview
========

cTreeBalls: Correlation Functions with Tree/Ball Methods
========================================================

Project Team
------------

:Author: Mario A. Rodriguez-Meza
:Major contributors: Alejandro Aviles, Eladio Moreno, Gustavo Niz, and
   Axel Romero Tisnado
:Repository: `Source code and issue tracking
   <https://github.com/rodriguezmeza/cTreeBalls>`_
:Documentation: `cTreeBalls on Read the Docs
   <https://ctreeballs.readthedocs.io/en/latest/>`_

Related Projects
----------------

* `3ptWL-cov <https://github.com/sadirs/3ptWL-cov>`_: computes Gaussian
  weak-lensing three-point covariance terms in a harmonic basis on the sphere.
* `3ptWL-mod <https://github.com/sadirs/3ptWL-mod>`_: models weak-lensing
  three-point correlation functions using perturbative, effective-field-theory,
  and halo-model methods.

Introduction
------------

cTreeBalls is a C code for computing two-point and three-point correlation
functions with tree- and ball-based neighbor searches.  It supports number
counts and scalar fields such as weak-lensing convergence, Cartesian and
unit-sphere geometries, multiple catalog formats, OpenMP parallelism, and a
Cython interface for Python workflows.

The three related projects cover complementary stages of a weak-lensing
analysis: cTreeBalls measures correlation functions from catalogs or maps,
3ptWL-mod predicts the three-point signal, and 3ptWL-cov computes its Gaussian
covariance contribution.

Installing and Getting Started
------------------------------

Clone and build the executable, static library, and Python extension:

.. code-block:: bash

   git clone https://github.com/rodriguezmeza/cTreeBalls.git
   cd cTreeBalls
   make clean
   make PYTHON=python3 all

Then run the compact synthetic example from the repository root:

.. code-block:: bash

   ./cballs nbody=4096 sizeHistN=12 mChebyshev=3 \
      rootDir=Output_quick numberThreads=1 verbose=0 verbose_log=0

The run writes radial bins, 2PCF data, 3PCF multipole files when enabled, a
used-values parameter file, and optional logs beneath ``Output_quick``.  See
:doc:`quickstart` for verification steps and :doc:`user/outputs` for the file
layout.

Configuration
-------------

Build-time switches live in ``Makefile_settings``,
``Makefile_machine``, and ``addons/Makefile_addons_settings``.  Runtime
parameters may be supplied as ``name=value`` command-line tokens or through a
parameter file.  See :doc:`user/configuration` and :doc:`user/command-line`.

Python
------

The Python extension is imported as:

.. code-block:: python

   from cyballs import cballs

It wraps the same initialization, search, histogram, and cleanup stages used by
the executable.  See :doc:`user/python` for a complete lifecycle example.

Documentation Builds
--------------------

Build the HTML guide from the repository root:

.. code-block:: bash

   python3 -m pip install -r docs/requirements.txt
   make -C docs html

The generated site is written to ``docs/_build/html/index.html``.

License
-------

cTreeBalls is open source and distributed under the MIT license.  See
:doc:`citing` for citation and acknowledgement guidance.

