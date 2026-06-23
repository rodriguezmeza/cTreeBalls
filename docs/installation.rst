Installation
============

The Python interface is distributed on PyPI as ``cTreeBalls`` and imported as
``cyballs``.  The installation builds a native extension from source using the
repository's bundled GSL and CFITSIO code.

Prerequisites
-------------

The source build requires a C compiler, ``make``, ``ar``, and Python
development headers.  On Debian or Ubuntu:

.. code-block:: bash

   sudo apt-get update
   sudo apt-get install build-essential python3-dev zlib1g-dev

Install the Python Interface
----------------------------

.. code-block:: bash

   python3 -m pip install cTreeBalls

Verify the compiled extension:

.. code-block:: bash

   python3 -c "from cyballs import cballs; print(cballs)"

The distribution name and import name intentionally differ: pip installs
``cTreeBalls``, while Python programs import the compiled extension as
``cyballs``.  The first installation can take several minutes because the
bundled native libraries are compiled locally.

Build All Interfaces from a Checkout
------------------------------------

Clone the repository when you also need the command-line executable, static
library, test catalogs, or an editable development checkout:

.. code-block:: bash

   git clone https://github.com/rodriguezmeza/cTreeBalls.git
   cd cTreeBalls
   make clean
   make PYTHON=python3 all

``make all`` builds ``cballs`` and ``libcballs.a`` and invokes the Python build
for ``cyballs``.  If your system exposes Python only as ``python3``, retain the
``PYTHON=python3`` override.

Machine-Specific Settings
-------------------------

Edit ``Makefile_machine`` only when the compiler, optimization flags, OpenMP
flags, or external-library paths differ from the defaults.  The principal
variables are ``CC``, ``PYTHON``, ``OPTFLAG``, ``OMPFLAG``, and ``NAGBODYDIR``.

Build Features
--------------

User-facing switches in ``Makefile_settings`` include:

.. list-table::
   :header-rows: 1
   :widths: 24 16 60

   * - Setting
     - Default
     - Purpose
   * - ``DEFDIMENSION``
     - ``3``
     - Select two- or three-dimensional coordinates.
   * - ``USEGSL``
     - ``1``
     - Enable GSL-backed functionality.
   * - ``GSLINTERNAL``
     - ``1``
     - Use the bundled GSL sources instead of a system installation.
   * - ``OPENMPMACHINE``
     - ``1``
     - Compile OpenMP search paths and runtime thread control.
   * - ``ADDONSON``
     - ``1``
     - Enable optional search methods, formats, parser support, and Cython hooks.

Changing a switch requires a clean rebuild:

.. code-block:: bash

   make clean
   make PYTHON=python3 all

Verify the Checkout Build
-------------------------

From the repository root:

.. code-block:: bash

   ./cballs --help
   ./cballs nbody=4096 sizeHistN=12 mChebyshev=3 \
      rootDir=Output_check numberThreads=1 verbose=0 verbose_log=0
   python3 -c "from cyballs import cballs; print(cballs)"

If the executable runs but the Python import fails, see :doc:`troubleshooting`.

Build the Documentation
-----------------------

.. code-block:: bash

   python3 -m pip install -r docs/requirements.txt
   make -C docs html
