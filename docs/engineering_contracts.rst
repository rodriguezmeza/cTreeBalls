Engineering contracts and acceptance ownership
==============================================

Numerical implementation boundaries
-----------------------------------

The two largest active implementation units delegate numerical work to private
headers in their own support directories:

* ``support/scalar_octree/scalar_harmonic_accumulation.h`` owns the angular
  recurrence and self-pair-subtracted outer-product accumulation.
* ``support/scalar_octree/scalar_edge_correction.h`` owns the complex harmonic
  window matrix and its pivoted solve, including existing validity states.
* ``support/shear_kernel/shear_complex_arithmetic.h`` owns scalar/vector complex
  arithmetic and accumulation.

These are internal implementation headers, included after the engine's types.
They preserve expression order, reduction order and compile-time feature guards.
They do not create a public ABI or replace the distinct standard-Legendre survey
solver. Existing independent scalar/shear oracles, exact/thread/MPI comparisons,
range tests and sanitizer runs remain required after edits. Reverse-include
regression selection follows these new headers; extraction does not reduce the
release matrix.

Configuration ownership
-----------------------

``general_libs/clib.c`` owns complete-token scalar conversion. The standalone
``GetbParam`` adapter now uses its strict Boolean conversion, like file and
Python routes. ``addons/class_lib/parser.c`` owns a shared unique-key lookup for
integer, double, Boolean and string scalar readers. Missing parameters preserve
existing defaults; duplicate or invalid scalar values do not mark the entry as
read or alter its destination. Parameter-file duplicate diagnostics and the
intentional CLI-over-file override remain the merger's responsibility.

The native ``test-configuration-readers`` target covers missing keys, duplicates,
unchanged destinations, valid conversions and rejected standalone Boolean suffixes.
The existing configuration suite covers CLI, parameter files, Python and recovery.
Addon-specific parameters and list grammars remain separate incremental units;
this is not a change to every configuration language.

Native top-level concurrency contract
-------------------------------------

**Serialize all public native entry in each process, across every model.** This
includes construction/activation, configuration and catalog registration, runs,
result access, cache operations, cleanup and destruction. There is no public
process-wide native mutex and separate contexts do not imply reentrancy.

The active-context selector and catalogue/tree globals still use save/restore;
parameter-parser state and GSL error-handler replacement are process-wide.
The functions named ``guarded`` protect error/allocation recovery, not simultaneous
entry. OpenMP worker parallelism is internal to one serialized top-level call.
Do not call the public API from these workers or reenter it from a callback.

The Cython entry points retain the conventional CPython GIL during native calls.
This does not promise an atomic multi-method Python transaction: NumPy conversion
and other Python work can yield. Applications sharing models or contexts between
threads must put one process-wide application lock around the complete operation,
including cleanup and mutation. Parallel independent jobs should use separate
processes created with ``spawn``; do not fork an initialized MPI/OpenMP runtime.
Free-threaded CPython and concurrent subinterpreters are outside the support
matrix. No GIL behavior is changed by this repair.

MPI additionally requires FUNNELED support and MPI main-thread entry. Every rank
participates in matching run order on a healthy communicator; an application lock
cannot repair inconsistent collective order. Cleanup is noncollective, but still
obeys the process serialization rule. Borrowed MPI is never finalized by a model.
See :doc:`performance` for parallel execution guidance.

Acceptance policy and source identity
-------------------------------------

``scripts/acceptance_metadata.py`` owns release-component obligations, fixed
accuracy thresholds and the canonical policy identity. ``accuracy_cases.py``
retains engine-specific measured settings. Release components, scalar/shear
acceptance, survey campaigns and performance measurements embed the same policy
record. Collectors reject missing/stale policy IDs and mixed source identities.
Thresholds are never selected from observed errors. Packaging pins Make queries
and static-library builds to the actual Python build interpreter, even when a
conflicting ``PYTHON`` variable or another interpreter is present on ``PATH``.

The source identity covers packaged numerical, parser, test and build inputs,
plus documentation, requirements, CI configuration and examples/notebooks.
Runtime timings and generated binaries are not source inputs. A software PASS
and a default-mode science qualification are distinct fields. Section 5.3's
rejected default science cases are not promoted by an engineering refactor.

Platform and dependency coverage
--------------------------------

``capabilities/support_matrix.json`` is the declared CI coverage matrix. Each
Linux x86-64 and macOS ARM64 platform has CPython 3.10, 3.11, 3.12 and 3.13 rows
with currently resolvable compatible dependencies. Each also has a CPython 3.10
minimum-dependency row: setuptools 68.2.2, wheel 0.41.3, Cython 0.29.36 and NumPy
2.0.2 at build time, followed by NumPy 1.22.4 at runtime. The other rows build and
run with NumPy >=2.0 and Cython >=0.29.36,<4. SciPy is installed in every runtime row, with a compatible version resolved
by pip. Actual resolved versions are retained.

Every row builds a wheel, installs it in a fresh environment outside the checkout,
and exercises installed scientific/lifecycle oracles. Native dependencies are
external GSL, CFITSIO, zlib and MPI; macOS uses a GNU OpenMP compiler consistently
through the MPI wrapper. The full public MPI gate and sanitizer jobs remain
separate mandatory checks. Hosted CI requires every matrix receipt and rejects
missing rows; a local single-platform PASS is not matrix certification.

The matrix describes standard-GIL CPython and the public double-precision profile.
Windows, free-threaded Python, other architectures and multi-node execution are
not certified by these jobs. A configured row is a target until its current
receipt succeeds. Installability beyond these rows is not an accuracy claim.

The runner labels follow `GitHub's runner reference
<https://docs.github.com/en/actions/reference/runners/github-hosted-runners>`_.
Build/runtime NumPy separation follows `NumPy's downstream ABI guidance
<https://numpy.org/doc/2.0/dev/depending_on_numpy.html>`_. Current receipts, rather
than these external documents, determine this project's actual compatibility.

Repeated dedicated-node measurement
-----------------------------------

Use ``scripts/benchmark_dedicated.py`` only after the combined software receipt
and bounded-accuracy acceptance pass for the exact source/build/policy. It times
the already accepted settings, not the unqualified fast defaults. Each requested
configuration must be present in the retained accuracy evidence.

The driver discards at least one cold-process warmup per configuration, shuffles
measured configurations with a retained seed, runs workers sequentially and
requires at least five measured repetitions. It retains every array, fixture,
phase timer, peak RSS and numerical comparison. Summaries give sample counts,
median, quartiles, median absolute deviation and extrema; no outlier deletion or
universal speed ranking is applied. Every measured sample must pass its fixed
exact-reference criterion. The existing 512-point workloads are calibration
workloads, not large-survey scalability measurements.

A cooperative per-user host lock prevents overlapping instances of this driver.
It does not detect unrelated workloads. ``--exclusive-node --allocation-id ID``
records an operator attestation; use an actual exclusive scheduler allocation.
``tests/slurm_scripts/benchmark_dedicated.slurm`` requests one exclusive node.
Set the site account/partition, CPUs per task and matching environment explicitly.
Do not benchmark a login node.

Local workflow validation is explicitly marked ``LOCAL_SMOKE_ONLY``::

   python scripts/benchmark_dedicated.py --smoke --repeats 2 --threads 2 \
     --engines kdtree-2balls-omp \
     --release /path/to/release-readiness.json \
     --accuracy /path/to/gate/accuracy/accuracy.json \
     --output /new/local-smoke

A dedicated run uses the same prerequisites and workload, without ``--smoke``::

   sbatch --account=YOUR_ACCOUNT --partition=YOUR_PARTITION --cpus-per-task=8 \
     tests/slurm_scripts/benchmark_dedicated.slurm \
     /path/to/release-readiness.json /path/to/gate/accuracy/accuracy.json \
     /new/dedicated-results

A different build on the dedicated node needs fresh matching correctness evidence
there before measurement. A local smoke never supplies dedicated-node evidence.
