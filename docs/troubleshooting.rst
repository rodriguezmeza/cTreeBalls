Troubleshooting
===============

Missing or Wrong cyballs
------------------------

Check the interpreter and imported file::

   python3 -c "import sys, cyballs; print(sys.executable); print(cyballs.__file__)"

Build from the source root using :doc:`installation`. A stale extension in
another directory can shadow the new build. Rebuild for the correct Python
minor version and restart notebooks after replacing the library.
Never manually edit generated ``ccyballs.pxd`` declarations to repair an ABI mismatch.

Missing mpi.h or MPI Runtime Errors
-----------------------------------

An enabled MPI addon needs an MPI compiler wrapper, not only MPI link flags.
Make selects the configured ``MPICC`` for MPI profiles. Ensure ``mpicc`` is
on PATH and uses a compiler compatible with OpenMP and the rest of the build.
Direct ``CC=...`` overrides must remain MPI-aware.

Python MPI jobs need mpi4py linked against the same implementation. Every
rank must enter the same method and cleanup sequence. Rank-0-only startup
or getters on worker ranks can fail or hang. Only the supplied drivers
automatically broadcast their catalog input.

GSL and CFITSIO Discovery
-------------------------

Bundled builds use ``GSLINTERNAL=1`` and ``CFITSIOLIBON=1``.
For external libraries, verify::

   gsl-config --cflags --libs
   pkg-config --cflags --libs cfitsio

Point ``PKG_CONFIG_PATH`` at the installation's pkgconfig directory when
needed. Make derives CFITSIO rpaths from the library search paths.
An import error such as ``symbol not found: _ffclos`` indicates unresolved
CFITSIO linkage or a stale/mismatched extension. Rebuild the static library
and wrapper with the same dependency profile; changing only the executable
does not repair an already built extension.

OpenMP Failures
---------------

The compiler, OpenMP flags, and runtime must agree. A macOS shared-memory
initialization failure in a subprocess can be an environment restriction,
not a failed numerical comparison. Report the exact runtime error.
Do not mix GNU and LLVM OpenMP unintentionally.

An OpenMP-named engine may require OpenMP at build time. Disabling
``OPENMPMACHINE`` is not a universal fallback for all addons.

C and Python Results Differ
---------------------------

Check the actual imported extension, ``cballs options=make-info``, compiled
method availability, selected catalogs, weights, masks, bins, and options.
Use a fresh output directory to avoid reading old files.

Scalar angular 3D positions must keep the observer origin. Current file and
in-memory paths preserve it. Old results made after catalog recentering or
with the previous angular formula are not interchangeable with corrected modes.
Use chord distances, not large-angle radians; see :doc:`3pcf`.

Raw KD-tree, legacy balltree, and BALLS4 now exclude repeated neighbors inside
the native kernel. Do not subtract those terms again in benchmark scripts.
The numerical test in :doc:`benchmarks` checks both input paths independently.

Smoothing Changed an Unsmoothed Run
-----------------------------------

With current code, ``SMOOTHPIVOTON=1`` only compiles capability;
``options=smooth-pivot`` must be explicitly requested. ``rsmooth`` alone
does not activate smoothing. Rebuild both interfaces if old behavior persists.
Raw KD/legacy balltree multipoles reject smoothing; BALLS4 has no smooth-pivot mode.

Zero or Unexpected Results
--------------------------

Check units, bin bounds, catalog selection, field type, compiled statistics,
and whether the requested output was actually computed. Coincident/radial/
antipodal legs have undefined angular bearings and are excluded from angular
3PCF while ordinary pairs still count.

Empty/singular angular edge systems return zero. Physical-3D survey products
also expose validity and conditioning diagnostics. These zeros are not
automatically evidence of a measured null signal. Narrow windows can amplify
noise even when a solver succeeds.

Getter Errors or Repeated-run Failures
--------------------------------------

Read getters only after the appropriate successful stage and before cleanup.
Use ``try/finally`` and ``struct_cleanup`` between retained-catalog runs or
``clean_all`` at completion. Getters return independent NumPy copies.
Report a minimal reproducer with profile, method, interpreter/module path,
catalog shape, and parameters if a crash persists.
