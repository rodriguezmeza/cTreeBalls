Build Configuration
===================

cTreeBalls uses three layers of build configuration:

``Makefile_settings``
    Core feature switches such as dimension, GSL, OpenMP, and add-ons.

``Makefile_machine``
    Compiler, optimization, Python executable, OpenMP flag, and external
    library paths.

``addons/Makefile_addons_settings``
    Optional search methods, I/O libraries, Cython declarations, and
    experimental features.

Core Settings
-------------

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   * - Setting
     - Default
     - Meaning
   * - ``DEFDIMENSION``
     - ``3``
     - Compile for two- or three-dimensional positions.
   * - ``USEGSL``
     - ``1``
     - Enable GSL-dependent routines.
   * - ``GSLINTERNAL``
     - ``1``
     - Build bundled GSL sources; set to ``0`` for a configured system GSL.
   * - ``OPENMPMACHINE``
     - ``1``
     - Enable OpenMP compilation and ``numberThreads``.
   * - ``ADDONSON``
     - ``1``
     - Include optional engines, formats, parser support, and Python hooks.

Important Add-on Settings
-------------------------

``OCTREEGGGOMPON``, ``KDTREEOMPON``, ``KDTREEBOXOMPON``, and
``NEIGHBORBOXESOMPON`` select search implementations.  ``CLASSLIBON`` and
``PXDON`` support the Cython wrapper.  ``IOLIBON`` and ``CFITSIOON`` enable
additional catalog formats.  ``TWOPCFON`` and ``TPCFON`` control compilation
of two- and three-point statistics.

Several lower sections of ``addons/Makefile_addons_settings`` are explicitly
marked experimental or development-only.  Do not enable them for production
without validating against the repository test suite.

Machine Variables
-----------------

``CC``
    C compiler, normally ``gcc``.

``PYTHON``
    Python executable used for the extension build.  Override with
    ``PYTHON=python3`` at make time.

``OPTFLAG``
    Optimization flags, currently ``-O3 -ffast-math``.

``OMPFLAG``
    OpenMP compiler flag, currently ``-fopenmp``.

``NAGBODYDIR``
    Base path used when external GSL or CFITSIO libraries are selected.

After Changing Settings
-----------------------

.. code-block:: bash

   make clean
   make PYTHON=python3 all

Record the three configuration files with benchmark and production outputs.
