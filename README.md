<p align="left">
    <img src="https://github.com/rodriguezmeza/cTreeBalls/blob/main/addons/cBalls_04.png" width="300" height="300">
</p>

# cTreeBalls: Correlation functions computation with Tree/Balls methods

Author: Mario A. Rodriguez-Meza

Major contributors:

- [Alejandro Aviles](https://github.com/alejandroaviles)
- [Eladio Moreno](https://github.com/Eladio-Moreno)
- [Gustavo Niz](https://github.com/gnizq64)
- [Sadi Ramirez](https://github.com/sadirs)
- Axel Romero Tisnado
- Sofia Samario


## Introduction

**cTreeBalls 1.1.0** provides the C executable `cballs`, the static library
`libcballs.a`, and the Python extension `cyballs`. Optional addons measure:

- Counts and convergence: scalar angular 2PCF/3PCF on flat or spherical catalogs.
- Shear: flat-sky spin-2 correlations and angular window correction.
- Lyman-alpha forests: anisotropic 3D and radial-only 2PCF/3PCF with forest-ID exclusions.
- Physical 3D scalar fields and ENCORE-style data/random survey estimators.

OpenMP and MPI availability depends on the selected build. Not every engine
accepts every geometry, field, mask, or normalization.

## Documentation

- [Sphinx guide](https://ctreeballs.readthedocs.io/en/latest/)
- [Installation](docs/installation.rst) and [build profiles](docs/build_profiles.rst)
- [Search-method guide](docs/search_methods.rst)
- [Scalar numerical contract](docs/3pcf.rst)
- [Python API](docs/api.rst) and [in-memory catalogs](docs/user/python.rst)
- [Benchmarking](docs/benchmarks.rst) and [regressions](tests/make_tests/README.md)

These source docs describe this checkout. Published packages and hosted docs
may lag behind local changes; uploading source does not publish a PyPI release.

## Build and Install

Work inside a Python virtual or Conda environment. Native builds need a C
compiler, OpenMP support for the OpenMP engines, make, ar, Python development
headers, and zlib. Enabled MPI addons additionally require an MPI C compiler
and runtime, normally `mpicc` and `mpiexec`.

```sh
git clone https://github.com/rodriguezmeza/cTreeBalls.git
cd cTreeBalls
python3 -m pip install numpy Cython setuptools wheel
make -j4 cballs cyballs-static-lib
CBALLS_STATIC_LIBRARY_READY=1 python3 setup.py build_ext --inplace --force
```

The last command is a checkout build, not a package installation. Use
`CBALLS_STATIC_LIBRARY_READY=1` only immediately after building the static
library with the same profile. To install this checkout in the active environment:

```sh
python3 -m pip install .
```

To install a published PyPI release, use `python3 -m pip install cTreeBalls`.
The import name is always `cyballs`. A pip installation does not install the
checkout's `cballs` command, examples, or benchmark source trees.

Configure `Makefile_settings`, `Makefile_machine`, and
`addons/Makefile_addons_settings` before building. Bundled GSL and CFITSIO
are selected with `GSLINTERNAL=1` and `CFITSIOLIBON=1`; system-library
discovery uses `gsl-config` and `pkg-config cfitsio`. Keep C and Cython
flags identical and rebuild both after changing a profile. Do not hand-edit
the generated `python/ccyballs.pxd`.

## Discover and Run

From the checkout root:

```sh
./cballs --help
./cballs options=make-info
./cballs options=print-options
./cballs options=print-search-methods
./cballs search=octree-ggg-omp nbody=512 sizeHistN=6 mChebyshev=3 \
    numberThreads=2 rootDir=Output_quick verbose=0 verbose_log=0
```

`make-info` reports the settings embedded in the executable, not merely the
current Makefiles. `print-search-methods` lists registered methods and their
requirements. Save the used-values file and build information with each run.

For scalar angular 3PCF in 3D, use observer-centered positions and Euclidean
chord bins. On a unit sphere, a separation angle `alpha` has chord length
`2*sin(alpha/2)`. **Do not recenter a sky catalog.**

The shared raw benchmark contract uses
`options=no-normalize-HistZeta,weights-norm` and ordered distinct triplets.
Repeated neighbors are removed inside the search. See the
[3PCF contract](docs/3pcf.rst) for complex-mode conventions and normalization.
Older affected multipoles must be recomputed.

`SMOOTHPIVOTON=1` only compiles support: `options=smooth-pivot` is required
to activate it on supported methods. Setting `rsmooth` alone does not.
Raw KD-tree/legacy balltree multipoles reject smoothing; BALLS4 has no
smooth-pivot support. Validate approximate tree acceptance against exact
small-catalog runs before interpreting speedups.

## Python and Notebooks

Register a catalog already in memory:

```python
import numpy as np
from cyballs import cballs

rng = np.random.default_rng(123)
xyz = rng.normal(size=(128, 3))
xyz /= np.linalg.norm(xyz, axis=1)[:, None]
model = cballs()
try:
    model.set(searchMethod="octree-ggg-omp", rootDir="Output_memory",
              rangeN=1.0, rminHist=0.05, sizeHistN=6, mChebyshev=3,
              numberThreads=2, useLogHist=True,
              options="no-normalize-HistZeta,weights-norm,no-one-ball")
    model.set_catalog(xyz, kappa=rng.normal(size=len(xyz)))
    model.Run()
    xi = model.getHistXi2pcf()
    cc, ss, sc, cs = [model.getHistZetaMsincos(2, t) for t in range(1, 5)]
    zeta_m1 = cc + ss + 1j*(sc-cs)
finally:
    model.clean_all()
```

Getters return NumPy copies. Registration retains prepared NumPy arrays;
each run copies them into C-owned storage. This avoids repeated file reads,
not all memory copies. Use `struct_cleanup()` between engines to retain
registered catalogs, and `clean_all()` when finished.

Ready-to-use examples:

- [In-memory Python example](examples/cyballs_in_memory_catalog.py) and
  [notebook](examples/cyballs_in_memory_catalog.ipynb)
- [Kappa all-engines driver](python/README_kappa_corr_all_engines.md):
  FITS/NPZ input, masks, complex edge corrections, OpenMP/MPI
- [Forest all-engines driver](python/README_lya_corr_all_engines.md) and
  [notebook](examples/lya_corr_all_engines.ipynb)
- [Physical 3D/ENCORE comparison](examples/compare_octree_3pcf_3d_encore.py)
  and [notebook](examples/compare_octree_3pcf_3d_encore.ipynb)

MPI Python runs require `mpi4py` linked to the same MPI implementation.
Every rank must enter the same run and cleanup sequence; read published
histograms on rank 0. The all-engines drivers manage catalog broadcasting.

## Tests and Documentation Builds

The documentation toolchain requires Python 3.11 or newer.

```sh
make test-make-info test-search-methods
python3 -m pytest -q tests/make_tests/test_scalar_numerical_contract.py
CBALLS_TEST_MPI=1 mpiexec -n 2 python3 tests/make_tests/test_scalar_numerical_contract.py
python3 -m pip install -r docs/requirements.txt
python3 -m sphinx -E -a -n -W --keep-going -b html docs docs/_build/html
```

Scalar numerical tests require a matching scalar C/Cython profile; MPI checks
also require the MPI siblings and mpi4py. Repeat with both smoothing build
settings. The main CI workflow covers profile-specific C/Cython, packaging,
MPI, and sanitizer jobs. Small correctness tests are not full-sky timing claims.

Open `docs/_build/html/index.html` after the Sphinx build. Additional plotting
notebooks are maintained in
[CBalls_plots](https://github.com/joar-cafe/CBalls_plots/tree/main/benchmarks).

## License

**cBalls** is written by Mario A. Rodriguez-Meza, is open source and distributed under the [MIT license](LICENSE). If you use this program in research work that results in publications, please cite the following paper:

Abraham Arvizu et al., [arXiv:2408.16847](https://arxiv.org/abs/2408.16847)

## Acknowledgements

cBalls use/is based on the following codes or projects:
-   [Zeno](https://home.ifa.hawaii.edu/users/barnes/zeno/index.html)
-   [Gadget-2](https://wwwmpa.mpa-garching.mpg.de/gadget/)
-   [CUTE](https://github.com/damonge/CUTE)
-   [Numerical Recipes](https://numerical.recipes/)
-   [GSL](https://www.gnu.org/software/gsl/)
-   [CLASS](https://github.com/lesgourg/class_public)
-   [CFITSIO](https://heasarc.gsfc.nasa.gov/fitsio/fitsio.html)
-   [HEALPix](https://healpix.sourceforge.io/)
-   [TreeCorr](https://github.com/rmjarvis/treecorr)
-   [FCFC](https://github.com/cheng-zhao/FCFC)

Also author acknowledges for helpful discussion and testing to the following people:

- Abraham Arvizu
- Alejandro Aviles
- Juan Carlos Hidalgo
- Eladio Moreno
- Gustavo Niz
- Axel Romero Tisnado
- Sofia Samario
