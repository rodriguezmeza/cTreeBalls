Installation
============

The distribution name is ``cTreeBalls``; import it as ``cyballs``.
These docs describe this source checkout. A published package or GitHub branch
may not yet contain local changes. A source push does not publish a PyPI release.

Prerequisites
-------------

Use an isolated Python virtual or Conda environment. A source build requires
a C compiler, make, ar, Python development headers, and zlib. OpenMP engines
require compatible compiler/runtime support. Any enabled MPI addon additionally
requires an MPI C compiler and runtime, normally ``mpicc`` and ``mpiexec``.

For an MPI-enabled source build on Debian/Ubuntu::

   sudo apt-get update
   sudo apt-get install build-essential python3-dev zlib1g-dev libopenmpi-dev openmpi-bin

On macOS, configure a compiler with working OpenMP support and an MPI wrapper
using a compatible toolchain. Do not mix unrelated MPI or OpenMP installations.
The deployment target is propagated through native and extension builds.

Python Installation
-------------------

Install a published release::

   python3 -m pip install cTreeBalls

Install a source checkout into the active environment::

   python3 -m pip install .

Install the version on the repository's default branch::

   python3 -m pip install "git+https://github.com/rodriguezmeza/cTreeBalls.git"

Pin a branch or commit explicitly for reproducible VCS installations. Pip builds
the extension using that source's selected profile; it does not install the
checkout's standalone ``cballs`` executable or external benchmark environments.

Verify the interpreter and module location::

   python3 -c "import cyballs; print(cyballs.__file__); print(cyballs.cballs().abi_sizes())"

Matching C and Cython Checkout Build
------------------------------------

For the executable, static library, examples, and a local extension::

   git clone https://github.com/rodriguezmeza/cTreeBalls.git
   cd cTreeBalls
   python3 -m pip install numpy Cython setuptools wheel
   make -j4 cballs cyballs-static-lib
   CBALLS_STATIC_LIBRARY_READY=1 python3 setup.py build_ext --inplace --force

Select features in the three Make configuration files before both commands.
The ready flag skips rebuilding the static library and is safe only when that
library has just been built with exactly the same profile. Restart existing
Python processes and notebook kernels after replacing a compiled extension.

``make PYTHON=python3 all`` is also supported, but includes the Makefile's
Python installation steps. Use the explicit workflow above when you want a
local checkout build without installing into the environment.

Overrides and Rebuilds
----------------------

Command-line overrides used only for make do not automatically reach a separate
setup.py invocation. Prefer editing the profile files. To use environment
overrides consistently in both commands, for example::

   export MAKEFLAGS=-e
   export SMOOTHPIVOTON=0
   make -j4 cballs cyballs-static-lib
   CBALLS_STATIC_LIBRARY_READY=1 python3 setup.py build_ext --inplace --force
   unset MAKEFLAGS SMOOTHPIVOTON

Environment precedence applies to all Make variables under ``-e``; use a clean,
intentional environment. See :doc:`build_profiles` for the ABI requirements.
Do not manually edit generated PXD declarations or copy an extension between
different Python minor versions.

Native Dependencies
-------------------

``GSLINTERNAL=1`` and ``CFITSIOLIBON=1`` select bundled sources.
``USEGSL=1`` is required by the current wrapper. FITS support additionally
requires ``CFITSIOON=1``.

For external GSL use ``GSLINTERNAL=0`` and ``gsl-config``, or explicit
``GSL_INCLUDE``/``GSL_LIB``. For external CFITSIO use ``CFITSIOLIBON=0``
and ``pkg-config cfitsio``. Set ``PKG_CONFIG_PATH`` to the directory containing
``cfitsio.pc`` for a custom installation. Make derives include, library, and
runtime-search flags; do not hardcode another user's absolute library path.

Verify and Continue
-------------------

From the checkout root::

   ./cballs options=make-info
   ./cballs options=print-search-methods

See :doc:`quickstart` for a small calculation, :doc:`search_methods` for
field/geometry selection, and :doc:`troubleshooting` for build/runtime errors.

Documentation Only
------------------

The pinned Sphinx toolchain requires Python 3.11 or newer. Documentation
builds do not import cyballs or compile C::

   python3 -m pip install -r docs/requirements.txt
   python3 -m sphinx -E -a -n -W --keep-going -b html docs docs/_build/html

Open ``docs/_build/html/index.html``. See ``docs/README.rst`` for other builders.
