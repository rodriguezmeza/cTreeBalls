Environment Controls
====================

Print the controls recognized by this source version and the environment
inherited by the executable:

.. code-block:: bash

   ./cballs options=print-environment
   ./cballs options=print-env

The report groups cTreeBalls memory limits, profiling, scalar/shear reuse and
task controls, external OpenMP and numerical-library settings, launcher
context, and Python helper settings. It displays quoted values, ``<unset>``
for absent variables, units, accepted ranges where applicable, and defaults.
Compiled reuse defaults come from this executable's build settings. External
library defaults depend on the linked runtime; their environment values may
be overridden by program parameters.

The sole-argument report exits successfully before catalog allocation or MPI
initialization, including when ``CBALLS_MEMORY_BUDGET_MB`` is invalid. It
neither changes definitions nor creates output files. Additional arguments
use normal checked startup. In a parameter file or Cython ``options`` setting,
the same report runs during startup. Append ``no-stop`` to continue:

.. code-block:: bash

   ./cballs options=print-environment,no-stop nbody=4096 numberThreads=4 rootDir=Output

Clear or temporarily override a variable
----------------------------------------

In bash or zsh, remove a definition from the current shell:

.. code-block:: bash

   unset CBALLS_MEMORY_BUDGET_MB
   unset CBALLS_SHEAR_PROFILE
   ./cballs options=print-environment

Unsetting the memory variable restores the compiled default of 65536 MiB
(64 GiB) per MPI rank for checked allocation/histogram plans. This is not a
whole-process memory cap; Python arrays, native trees, and other allocations
also use RAM. Increasing a budget does not supply physical memory.

``CBALLS_SHEAR_PROFILE`` uses presence semantics: setting it to ``0`` or
``""`` still enables profiling. An empty memory definition is invalid.
Use ``unset NAME`` rather than ``NAME=""`` to remove a definition.

To remove overrides for one child command while preserving the shell:

.. code-block:: bash

   env -u CBALLS_MEMORY_BUDGET_MB -u CBALLS_SHEAR_PROFILE ./cballs options=print-env
   env -u CBALLS_MEMORY_BUDGET_MB python3 tests/python/lya_corr_all_engines.py --synthetic

The report prints complete, copyable ``unset`` commands for cTreeBalls-owned
controls, plus separate optional commands for OpenMP, numerical libraries,
and frontend helpers. Unsetting a reuse tolerance restores its compiled
fallback; remove the corresponding ``scalar-pivot-reuse`` or
``shear-pivot-reuse`` option token if you also want to stop requesting reuse.
The shear-specific controls do not tune Ly-alpha engines. Ly-alpha kernel
and bin-slop controls are parameters; general search parameters such as
``theta``, ``nsmooth``, and ``rsmooth`` have engine-specific support.
See ``options=make-info`` and
``options=print-search-methods`` for compiled settings and engine support.

For an existing Python process:

.. code-block:: python

   import os
   os.environ.pop("CBALLS_MEMORY_BUDGET_MB", None)

Change native controls before ``Run()``; do not mutate environment variables
while a calculation is running. OpenMP/BLAS libraries may consume settings
at import or initialization, so set those before launching Python and restart
an existing interpreter when needed. ``numberThreads`` (``--threads`` in
the drivers) explicitly requests the native OpenMP team size and overrides
the initial ``OMP_NUM_THREADS`` value, subject to runtime limits.

A subprocess cannot unset a variable in its parent shell. To prevent a
definition from returning in future sessions, remove its ``export`` from
the relevant shell startup file, scheduler script, or environment activation
configuration. Keep MPI/scheduler-provided variables intact inside jobs;
they appear in the diagnostic section but are excluded from reset commands.
