Build Configuration
===================

Three files select the native and Python profile:

* ``Makefile_settings``: dimensions, GSL, OpenMP, and addon support.
* ``Makefile_machine``: toolchain, precision, optimization, and library discovery.
* ``addons/Makefile_addons_settings``: engines, I/O, Cython hooks, and optional features.

See :doc:`../build_profiles` for dependency and ABI details.
Use :doc:`../search_methods` to select the corresponding engine switches.

Important Controls
------------------

``DEFDIMENSION``
    Body coordinate dimension, 2 or 3. Radial forest searches still require 3.

``TWOPCFON`` / ``TPCFON``
    Compiled scalar correlation orders. Runtime selectors are engine-specific.

``SINGLEPON``
    ``0`` retains double storage/computation; ``1`` selects mixed precision:
    compact storage for selected data with double computation/accumulation.
    This is not a blanket float32 histogram build. Validate numerical errors.

``SMOOTHPIVOTON``
    Compile pivot-smoothing capability. It does not activate smoothing without
    ``options=smooth-pivot``. Unsupported engines reject the runtime request.

``CLASSLIBON`` / ``PXDON``
    Parser and wrapper hooks required by cyballs, together with ``ADDONSON=1``
    and ``USEGSL=1``.

``CC`` / ``MPICC``
    C compiler and MPI wrapper. Enabled MPI addons select MPI-aware compilation.
    The selected compiler and its OpenMP runtime must agree.

``OPTFLAG``
    The source profile uses ``-O3 -fno-fast-math``. Do not assume historical
    ``-ffast-math`` instructions preserve the numerical contract.

``PYTHON``
    Interpreter used by Make's Python targets. Select the same environment
    for setup.py, pip, tests, and notebooks.

Inspect and Rebuild
-------------------

``./cballs options=make-info`` reports the executable's compiled settings.
``make --no-print-directory -s print-cyballs-build-env`` reports the profile
that a new wrapper build will consume. They answer different questions.

Rebuild both C and Cython after changing a profile; use the matching commands
in :doc:`../installation`. A command-line make override alone is not carried
into an independent Python build. Retain the profile and actual module path
alongside scientific outputs.
