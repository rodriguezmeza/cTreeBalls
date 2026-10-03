# Source distribution contract

The public `cyballs` source archive uses external GSL and CFITSIO. Its staged
`Makefile_settings` contains `GSLINTERNAL = 0` and its staged
`addons/Makefile_addons_settings` contains `CFITSIOLIBON = 0`. The `sdist` command
normalizes those two settings before archiving. It does not edit the development
checkout or change engine selection, precision, estimator options or other
settings. Replacements detach files even if setuptools stages them as hard links.
The bundled source trees remain excluded from the public archive.

Build prerequisites are make/ar, a compatible C compiler and OpenMP runtime,
Python development headers, zlib, external GSL and CFITSIO, and MPI when the
packaged profile enables MPI engines. Python build isolation does not install
these operating-system dependencies. On Debian/Ubuntu the native prerequisites
include `build-essential python3-dev pkg-config libgsl-dev libcfitsio-dev
zlib1g-dev libopenmpi-dev openmpi-bin`. On macOS use a consistent OpenMP/MPI C
toolchain and installed external libraries.

Make discovers GSL with `gsl-config` (or `GSL_CONFIG`) and CFITSIO with
`pkg-config cfitsio` (`PKG_CONFIG` and `CFITSIO_PKG` are configurable).
`GSL_INCLUDE`/`GSL_LIB` and `CFITSIO_INCLUDE`/`CFITSIO_LIB` support custom prefixes.
The Cython build consumes the resolved Make discovery flags too, including
library and runtime-search flags. Failed discovery produces an actionable error;
it does not silently assume a bundled source directory or `/usr/local` install.

## Build and verify the artifact

From a checkout with its own native/Python build prerequisites installed:

```sh
python -m pip install build twine
python -m build --sdist
python -m twine check --strict dist/*.tar.gz
python scripts/check_sdist.py dist/cyballs-VERSION.tar.gz
python scripts/verify_sdist_install.py dist/cyballs-VERSION.tar.gz \
  --output /absolute/path/outside-the-checkout/artifact-check --jobs 4
```

The checker inspects the actual archive without executing its Makefiles. It
requires the external defaults, generated/public build inputs, recent numerical
fixes and their regression tests, while rejecting bundled libraries, unsafe
paths and binary build artifacts. This structural check is only the first stage.

The verification script requires a fresh output directory. It extracts the
archive twice, builds the standalone program with the archive defaults in one
copy, and builds a wheel with PEP 517 Python build isolation in the other. It
installs that wheel into a fresh virtual environment and runs outside both
source copies, with user-site and checkout `PYTHONPATH` imports disabled. No
`GSLINTERNAL`/`CFITSIOLIBON` overrides are injected by the verifier.

Verification checks the resolved/native profile, external dependency settings,
installed source identity and complete engine registry, and C/Cython ABI through
object construction. It runs independent scalar complex-window oracles and known
weighted LOS pair/triplet totals. Logs, commands, exit codes, archive/wheel hashes,
build identity and numerical results are retained in `verification.json` and
adjacent files. A failing subprocess makes verification fail. Checks remain
active under optimized Python; numerical-oracle workers run with assertions on.

Use the build prerequisites selected by `pyproject.toml`; build isolation can
select a different NumPy/Cython version from the development environment. The
resulting wheel is specific to its Python ABI/platform and external native
libraries. This workflow does not make a portable manylinux/macOS wheel and does
not publish anything. The publishing workflow continues to publish the checked
source archive only.

Developer builds with bundled libraries remain supported when those source trees
are present. Public source archives cannot enable those missing trees without
first supplying them. Use the archive profile as the expected profile for an
installed-artifact check, rather than a developer checkout's bundled defaults.
