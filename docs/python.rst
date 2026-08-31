Python Interface Reference
==========================

The module is ``cyballs`` and the class is ``cballs``. The complete callable
reference is :doc:`api`; :doc:`user/python` provides an in-memory workflow.

Lifecycle
---------

``set_default(**overrides)`` restores compiled C defaults and stages overrides.
``set`` accepts one mapping or keyword arguments. Unknown names are rejected
when ``Run`` passes them through the compiled parser.

``Run(level=["MainLoop"])`` executes initialization and computation while
keeping results alive. Read getters before requesting ``EndRun`` or cleanup.
Arrays returned by histogram getters are NumPy copies.

Use ``struct_cleanup`` to release C state while retaining registered catalogs
for another engine. Use ``clean_all`` when finished; it also clears parameters
and catalog references. Cleanup belongs in a ``finally`` block.

Catalog and Result APIs
-----------------------

* ``set_catalog`` registers scalar, weighted, masked, or shear catalogs.
* ``set_forest_catalog`` registers forest pixels with integer quasar IDs.
* ``clear_catalogs`` drops registrations.
* ``getrBins``, ``getHistNN`` and ``getHistXi2pcf`` expose supported pair products.
* ``getHistZetaMsincos(m, type)`` exposes scalar angular components.
* ``getHistZetaM_EE_complex(m)`` returns both parts of complex corrected modes.

Scalar multipole arguments are one-based: ``m=1`` means order zero.
The raw reconstruction is ``coscos+sinsin + 1j*(sincos-cossin)``.
Individual components depend on the tangent basis; use the complex combination
for rotation and reference comparisons. See :doc:`3pcf`.

Registration prepares/retains NumPy arrays and each run copies them to C-owned
body storage. It avoids repeated I/O, not all allocation or copying.
``getNBody`` preserves the compiled INTEGER width as a Python integer.

Build Coupling
--------------

Only the root ``setup.py`` builds the extension. It generates
``python/ccyballs.pxd`` from the template and selected feature flags.
Do not manually set NDIM or add profile-specific fields to the generated file.
After changes, rebuild C and Cython together and restart Python.
See :doc:`build_profiles`.

For MPI, use matching mpi4py/native MPI libraries, identical rank lifecycle
calls, and rank-0 result access. See :doc:`search_methods`.
