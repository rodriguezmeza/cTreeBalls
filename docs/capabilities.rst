Capabilities and Validation
===========================

The public catalogue, ``capabilities/engines.json``, declares the 34 active
search names, stable IDs, compile conditions, aliases, estimator descriptions
and regression ownership. Inactive addon declarations are not distributed.
Required shared internal kernels do not register extra search methods.

Run ``python3 scripts/generate_capabilities.py`` after editing the catalogue.
It generates the C registry, registration fragments, Python gate declarations
and ``ENGINE_CAPABILITIES.md``. Native builds reject stale generated files.
``options=print-search-methods`` uses this registry and displays only compiled
entries. ``options=make-info`` reports the resolved build, while
``options=print-options`` describes runtime controls.

Validation Commands
-------------------

Run from the repository root with matching compiler, MPI and Python environments::

   make check-capabilities
   make test-search-methods test-make-info test-benchmark-drivers
   make test-resource-contracts test-runtime-context
   python3 scripts/affected_regressions.py --changed source/mpi_runtime.c --output affected.json
   python3 scripts/active_release_gate.py --output /new/release-check --jobs 4

The gate requires real MPI execution for compiled MPI methods. Its retained
fixtures are numerical/runtime regressions, not guarantees for every survey,
approximation setting or production catalog size. Affected-regression selection
never replaces complete release verification.

Timing and Accuracy
-------------------

The three all-engines drivers live in ``tests/python``. Their
``mainloop_wall_s`` and ``mainloop_cpu_s`` columns exclude Python provenance
capture. ``compute_*`` times the complete ``Run`` call. Each driver documents
its setup, catalog-registration and cleanup scopes. MPI uses maximum rank wall
time and summed process CPU, with per-rank measurements retained.

For cold-process phase, memory and accuracy records::

   python3 scripts/benchmark_contracts.py --output /new/benchmark --threads 2 --repeats 3

These fixed small workloads retain exact references and approximation errors.
L2 acceptance limits are not per-bin relative-error guarantees. Peak RSS is a
cumulative process high-water mark, not independent phase allocations. See
:doc:`benchmarks` for comparison requirements.
