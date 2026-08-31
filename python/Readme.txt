cyballs Python interface

Build from the cTreeBalls root, not this python directory. setup.py is at
the repository root; python/ccyballs.pxd is generated and must not be edited.

In the intended Python environment, after selecting one consistent Make profile:

    python3 -m pip install numpy Cython setuptools wheel
    make -j4 cballs cyballs-static-lib
    CBALLS_STATIC_LIBRARY_READY=1 python3 setup.py build_ext --inplace --force

The ready flag is safe only for the static library just built with the same
profile. For installation into the environment instead, run:

    python3 -m pip install .

Import with: from cyballs import cballs
The distribution name is cTreeBalls. Keep -fPIC and matching OpenMP/MPI
runtimes. GSL is supported and USEGSL=1 is required for this wrapper;
GSLINTERNAL=1 uses the bundled sources. Rebuild for every Python ABI or
feature-profile change, then restart existing notebook kernels.

set_catalog registers in-memory scalar/shear arrays; set_forest_catalog also
carries integer forest IDs. Prepared NumPy arrays are retained and copied into
C-owned body storage at each run. struct_cleanup frees C state but keeps the
registered arrays; clean_all also clears the arrays and parameters.
Histogram getters return copies and require the appropriate completed run.

See ../docs/api.rst and ../docs/user/python.rst.
See README_kappa_corr_all_engines.md and README_lya_corr_all_engines.md for
one-read multi-engine drivers, masks, MPI, plotting, and output conventions.
Examples and notebooks are in ../examples.
