Code structure
==============

This page is a developer-oriented map of the cTreeBalls source tree. It is
intended to answer two questions before editing the code: where does a run move
through the program, and which files own each piece of state.

Runtime flow
------------

The command-line executable follows this high-level sequence:

1. ``main/main.c`` creates the run state objects and calls ``StartRun``.
2. ``source/startrun.c`` reads parameters, checks them, initializes random
   generators, prepares output paths, and allocates the run-owned arrays.
3. ``source/cballs.c`` enters ``MainLoop`` and dispatches to ``EvalHist``.
4. ``source/treeload.c`` builds one tree per selected input catalog with
   ``MakeTree``.
5. ``source/search.c`` walks the tree, accepts or rejects node pairs, and
   accumulates thread-local histogram contributions.
6. ``source/cballsutils.c`` supplies shared search utilities, histogram
   allocation and normalization helpers, geometry transforms, and acceptance
   predicates.
7. ``source/cballs.c`` writes evaluated histograms through its ``PrintHist*``
   helpers.
8. ``source/cballsio.c`` closes files and releases memory in ``EndRun``.

The Cython wrapper in ``python/cyballs.pyx`` calls the same C stages, but it
splits them into a dependency-checked ``Run(level=[...])`` workflow so Python
callers can compute and then read arrays back without going through the command
line.

Module map
----------

.. list-table::
   :header-rows: 1
   :widths: 28 72

   * - Path
     - Responsibility
   * - ``main/main.c``
     - CLI entry point. Creates ``struct cmdline_data`` and
       ``struct global_data``, records start times, reads the parameter file
       argument when ``GETPARAM`` is disabled, and runs ``StartRun``,
       ``MainLoop``, then ``EndRun``.
   * - ``include/globaldefs.h``
     - Central include file. Pulls in common libraries, declares global tables
       such as ``bodytable``, ``roottable``, and ``nodetablescanlev``, selects
       GSL or Numerical Recipes support, and exposes project-wide error and
       debug macros.
   * - ``include/cmdline_defs.h``
     - Parameter registry used by the ``getparam`` interface. Add command-line
       parameters here first, then wire them through ``startrun.c``.
   * - ``include/cmdline_data.h``
     - Parsed user options and derived command state. This is the user-facing
       configuration object passed through almost every major routine.
   * - ``include/global_data.h``
     - Run state that is derived, allocated, or accumulated during execution:
       file paths, histogram arrays, tree statistics, timing counters, input
       catalog metadata, and memory ownership flags.
   * - ``include/datastruc_defs.h``
     - Body, cell, and node layout. ``body`` and ``cell`` both embed a common
       ``node`` header so tree traversal can work on ``nodeptr`` values.
   * - ``include/datastruc_hist.h``
     - Per-thread and global histogram workspaces used by the search code.
   * - ``include/protodefs.h``
     - Shared function prototypes and add-on include hook for APIs that need to
       be visible across C modules.
   * - ``source/startrun.c``
     - Startup lifecycle: command-line and parameter-file parsing,
       ``StartRun_Common``, validation, data loading or test-data generation,
       output setup, and memory allocation.
   * - ``source/cballsio.c``
     - Input/output formats, output directory and filename setup, writing
       catalogs, and all teardown routines called from C and Cython.
   * - ``source/testdata.c``
     - Synthetic catalog generation for simple cubic and unit-sphere runs.
   * - ``source/treeload.c``
     - Octree construction, root-cell setup, body insertion, cell properties,
       tree threading, scan-level tables, and tree cleanup.
   * - ``source/search.c``
     - The default ``octree-sincos-omp`` search. Owns the OpenMP pivot loop,
       tree walk, node summation, and merging of thread-local histograms.
   * - ``source/cballsutils.c``
     - Shared search utilities: histogram initialization, body property
       reduction, acceptance/rejection checks, periodic wrapping, coordinate
       conversion, and histogram diagnostics.
   * - ``source/cballsutils_with-RADECPOS.c``
     - RA/DEC-position variant of ``cballsutils.c``. Keep it in sync with
       ``cballsutils.c`` when enabling or maintaining that build path.
   * - ``source/cballs.c``
     - Main computation driver and histogram output dispatch. This is where
       search-method selection reaches the concrete search implementation.
   * - ``python/cyballs.pyx``
     - Python class wrapper around the C lifecycle. Converts a Python parameter
       dictionary into CLASS-style ``file_content``, calls C modules, exposes
       histogram getters, and releases C allocations.
   * - ``python/ccyballs.pxd``
     - Cython declarations for the C structs and functions that Python needs.
       Only expose fields here when the wrapper really consumes them.
   * - ``addons/``
     - Optional search methods, input formats, and bundled libraries. Add-ons
       connect through the ``//B socket:`` include points in core files.

Core state model
----------------

Most C routines receive two pointers:

``struct cmdline_data *cmd``
    Values requested by the user or read from a parameter file. Treat this as
    configuration. It contains items such as ``searchMethod``, histogram sizes,
    input/output names, verbosity, and option strings.

``struct global_data *gd``
    Values computed during the run. Treat this as mutable execution state. It
    owns histogram buffers, input catalog metadata, output paths, timing
    counters, tree statistics, and allocation flags used during cleanup.

The large catalog and tree tables are declared as project globals in
``include/globaldefs.h``:

``bodytable[MAXITEMS]``
    Per-catalog body arrays.

``roottable[MAXITEMS]``
    One tree root per loaded catalog.

``nodetablescanlev[MAXITEMS]`` and ``nodetablescanlev_root[MAXITEMS]``
    Flattened node tables used by scan-level traversal and diagnostics.

The body/cell hierarchy is defined in ``include/datastruc_defs.h``. A ``body``
contains a ``node`` named ``bodynode`` and a body id. A ``cell`` contains a
``node`` named ``cellnode``, a size, a linked ``more`` pointer, and ``NSUB``
subcells. The macros in that header, for example ``Pos(x)``, ``Mass(x)``,
``Radius(x)``, and ``Subp(x)``, are the intended field access style throughout
the C code.

Tree and search lifecycle
-------------------------

``MakeTree`` in ``source/treeload.c`` starts with ``newtree`` and a root cell,
then loads bodies into subcells using ``loadbody`` and ``subindex``. After
insertion, ``hackcellprop`` computes cell-level properties used by the search,
``threadtree`` links traversal order, and ``scanLevel`` prepares node tables
for level-based operations.

``searchcalc_normal_sincos`` in ``source/search.c`` is the default search path
for ``searchMethod=tree-omp-sincos``/``octree-sincos-omp`` style runs. It:

1. Initializes global and per-thread histogram workspaces.
2. Runs an OpenMP loop over pivot bodies from the selected pivot catalog.
3. Walks the scanning catalog tree with ``normal_walktree_sincos``.
4. Uses acceptance/rejection helpers from ``source/cballsutils.c`` to decide
   whether a body/cell or cell/cell pair contributes.
5. Accumulates 2PCF and, when enabled at compile time, 3PCF Chebyshev/sine/cosine
   components.
6. Merges thread-local histograms back into ``gd``.

Histogram ownership
-------------------

Histogram buffers are allocated during startup and search initialization, then
freed through the ``EndRun_FreeMemory_*`` routines in ``source/cballsio.c``.
The important ownership rule is: if a new allocation is attached to ``cmd`` or
``gd``, it needs a matching allocation flag or an existing cleanup stage that
cannot miss it. This is especially important for the Python wrapper, because
``cballs.struct_cleanup`` calls the granular C cleanup routines instead of only
exiting the process.

Common histogram outputs include:

``histNN``
    Neighbor count histogram.

``histCF``
    Correlation-function style output derived from counts.

``histXi2pcf``
    Two-point correlation function output.

``histZetaM*`` and ``mhistZetaM*``
    Three-point/Chebyshev multipole outputs, split by sine/cosine component and
    multipole index when the corresponding compile-time options are enabled.

Adding or changing parameters
-----------------------------

To add a user-visible parameter:

1. Add the default, description, and optional short name in
   ``include/cmdline_defs.h``.
2. Add fields to ``struct cmdline_data`` in ``include/cmdline_data.h`` if the
   value needs to persist.
3. Read the parameter in ``ReadParametersCmdline`` and ``ReadParameterFile`` in
   ``source/startrun.c``.
4. Print it from ``PrintParameterFile`` so generated parameter snapshots remain
   reproducible.
5. Validate derived constraints in ``CheckParameters`` when invalid
   combinations would otherwise fail later.
6. If the parameter is exposed to Python, mirror it in ``python/ccyballs.pxd``
   and set or document it in ``python/cyballs.pyx``.

Adding input/output formats
---------------------------

Input format names are converted to integer tags by
``infilefmt_string_to_int`` in ``source/cballsio.c``. Add a new tag in the same
style as the existing formats, implement the parser near the other
``inputdata_*`` functions, and add the switch case in ``InputData`` or
``InputData_all_in_one``.

Output formats follow the matching pattern through
``outfilefmt_string_to_int``, ``outputdata``, and the ``outputdata_*`` helpers.
Keep filename construction in ``setFilesDirs``/``setFilesDirs_log`` so command
line, tests, and Python runs all agree on output paths.

Adding search methods
---------------------

Search-method selection happens in two places:

1. ``search_method_string_to_int`` in ``source/startrun.c`` maps user strings
   to integer constants.
2. ``EvalHist`` in ``source/cballs.c`` switches on ``gd->searchMethod_int`` and
   calls the concrete implementation.

New search methods usually need prototypes in ``include/protodefs.h``, compile
settings in the Makefiles or add-on make fragments, and output handling in
``source/cballs.c`` only if they produce different histogram families.

Add-on sockets
--------------

Core files contain ``//B socket:`` blocks guarded by ``#ifdef ADDONS``. These
blocks include generated or configured files such as ``*_include.h`` so optional
modules can add fields, prototypes, switch cases, and setup logic without
editing every core file manually. When adding an add-on, keep its ownership
localized under ``addons/<name>/`` and document which socket files it expects to
be copied or generated.

Python wrapper notes
--------------------

The wrapper class ``cballs`` in ``python/cyballs.pyx`` keeps a Python parameter
dictionary in ``self._pars``. ``Run`` converts that dictionary into
``file_content``, calls ``input_read_from_file``, then advances through
``StartRun_Common``, ``PrintParameterFile``, ``SetNumberThreads``, and
``MainLoop`` according to the requested dependency level.

Getter methods such as ``getrBins``, ``getHistNN``, ``getHistXi2pcf``, and
``getHistZetaMsincos`` copy C arrays into NumPy arrays. They assume the relevant
run stage has already completed and that ``sizeHistN`` still matches the
allocated histogram dimensions.

When changing C struct fields consumed by Python, update both
``python/ccyballs.pxd`` and the getter or setup code in ``python/cyballs.pyx``.
The wrapper intentionally exposes only a subset of the C structs; avoid adding
fields to the ``.pxd`` unless Python needs direct access.

Generated and vendored code
---------------------------

The repository includes generated Sphinx output under ``docs/_build`` and large
vendored libraries under ``addons/``. Prefer editing the source ``.rst`` files,
core ``source/`` and ``include/`` code, or add-on source fragments rather than
generated HTML/LaTeX or third-party library internals.
