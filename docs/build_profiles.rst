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

Typical Source Profile
----------------------

The source selects the following core features; local edits and packaging
profiles can differ. Inspect the executable rather than inferring its build
from this table:

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
   * - ``SINGLEPON``
     - ``0``
     - Double storage/computation; ``1`` selects mixed storage with double arithmetic.
   * - ``SMOOTHPIVOTON``
     - ``0``
     - Capability only; when compiled, ``smooth-pivot`` remains runtime opt-in.

For MPI addons, the build selects the configured ``MPICC`` wrapper.
The runtime must provide ``MPI_THREAD_FUNNELED`` or stronger support.
See :doc:`search_methods` for independent OMP/MPI addon switches.

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
   ./cballs options=make-info

The Make query prints flags for a new wrapper build. ``make-info`` prints
the settings embedded in the existing executable; stale binaries may differ.

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

   make -j4 cballs cyballs-static-lib
   CBALLS_STATIC_LIBRARY_READY=1 python3 setup.py build_ext --inplace --force
   python3 -c "from cyballs import cballs; print(cballs().abi_sizes())"

Use exactly the same profile for both steps. A one-off make override does not
automatically reach setup.py's Make queries; use profile files or the consistent
environment workflow in :doc:`installation`. Keep ready-library builds limited
to freshly built matching archives. Restart Python after rebuilding.

For documentation-only validation:

.. code-block:: bash

   python3 -m sphinx -M html docs /tmp/cTreeBalls-sphinx-build -W --keep-going
