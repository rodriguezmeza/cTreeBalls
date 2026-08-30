Build Profiles
==============

cTreeBalls is compiled from a Makefile profile. The active profile controls the
C executable, ``libcballs.a``, and the ``cyballs`` Cython extension. The C and
Cython builds must use the same feature macros.

Profile Layers
--------------

``Makefile_settings``
    User-facing core switches: dimension, GSL, OpenMP, and add-ons.

``Makefile_machine``
    Compiler, Python executable, optimization flags, OpenMP flags, integer
    precision, and external GSL/CFITSIO discovery.

``addons/Makefile_addons_settings``
    Optional search methods, I/O formats, Cython/PXD support, and experimental
    features.

Current Default Profile
-----------------------

The current default source profile is:

.. list-table::
   :header-rows: 1

   * - Setting
     - Value
     - Meaning
   * - ``DEFDIMENSION``
     - ``3``
     - Compile with ``THREEDIMCODE``.
   * - ``USEGSL``
     - ``1``
     - Enable GSL-backed paths.
   * - ``GSLINTERNAL``
     - ``1``
     - Build the bundled/internal GSL sources.
   * - ``CFITSIOON``
     - ``1``
     - Enable FITS/HEALPix I/O support.
   * - ``CFITSIOLIBON``
     - ``1``
     - Build/use the bundled/internal CFITSIO library.
   * - ``OPENMPMACHINE``
     - ``1``
     - Enable OpenMP and ``OPENMPCODE``.
   * - ``ADDONSON``
     - ``1``
     - Include add-on engines, formats, CLASS parser, and PXD hooks.
   * - ``LONGINTON``
     - ``1``
     - Compile ``INTEGER`` as ``long``.

Native Libraries
----------------

When ``GSLINTERNAL = 0``, GSL is discovered with ``gsl-config`` unless
``GSL_INCLUDE`` and ``GSL_LIB`` are provided manually.

When ``CFITSIOON = 1`` and ``CFITSIOLIBON = 0``, CFITSIO is discovered with
``pkg-config cfitsio`` unless ``CFITSIO_INCLUDE`` and ``CFITSIO_LIB`` are
provided manually.

Inspecting The Active Profile
-----------------------------

Use this command before debugging a build or Python ABI mismatch:

.. code-block:: bash

   make --no-print-directory -s print-cyballs-build-env

It prints the flags consumed by ``setup.py`` when building ``cyballs``.

Cython ABI Coupling
-------------------

The Python extension is not independent of the C profile. ``setup.py`` queries
the Makefile profile, generates ``python/ccyballs.pxd`` from
``python/ccyballs.pxd.in``, builds ``libcballs.a``, and compiles ``cyballs``
with the same macro set.

The wrapper requires:

* ``ADDONSON = 1``
* ``CLASSLIBON = 1``
* ``PXDON = 1``
* ``USEGSL = 1``

The generated PXD must match the C structs for active choices such as
``LONGINT``, ``THREEDIMCODE``, ``NMultipoles``, and enabled add-ons. At import
time, ``cyballs`` checks the C and Cython sizes of ``cmdline_data`` and
``global_data``.

Safe Rebuild Workflow
---------------------

After changing any build profile setting:

.. code-block:: bash

   make clean
   make PYTHON=python3 all
   python3 -c "from cyballs import cballs; print(cballs().abi_sizes())"

For documentation-only validation:

.. code-block:: bash

   python3 -m sphinx -M html docs /tmp/cTreeBalls-sphinx-build -W --keep-going
