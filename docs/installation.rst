Installation
============

Prerequisites
-------------

A source build requires a C compiler, ``make``, and Python development tools
when the ``cyballs`` extension is needed.  The default configuration enables
OpenMP, bundled GSL sources, the add-on system, and bundled CFITSIO support.

On Debian or Ubuntu, a typical preparation is:

.. code-block:: bash

   sudo apt-get update
   sudo apt-get install build-essential python3-dev python3-pip
   python3 -m pip install --user numpy Cython

Clone and Build
---------------

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

Verify the Installation
-----------------------

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

