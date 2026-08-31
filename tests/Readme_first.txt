
Current regression entry points

Use tests/make_tests/README.md and docs/development.rst. The maintained
launchers and associated Python tests are in tests/make_tests.
With matching C and Cython binaries, run from the repository root:

    make test-make-info test-search-methods
    python3 -m pytest -q tests/make_tests/test_scalar_numerical_contract.py

Repeat scalar numerical tests with both SMOOTHPIVOTON build settings.
That switch compiles capability; it does not activate smoothing without
options=smooth-pivot. Preserve your own build setting after testing.

Historical catalog/plot workflow

The following notes require external catalogs and reference outputs.
Older multipoles may predate the corrected angular/normalization contract;
regenerate affected references before interpreting their differences.

1. To test cBalls:

In addons/Makefile_addons_settings do:

SMOOTHPIVOTON = 0
TPCFON = 1

and in cTreeBalls directory:

$ make clean; make all

$ cd tests

Check in the run_all_tests script that definition of CBALLS variable points to where cBalls exe is located. Then,

$ ./scripts/run_all_tests

You will get several pdf plots showing horizontal lines indicating comparisons versus reference results in Outputs_to_compare_with folder are OK.

Or test individual cases:

$ time ../cballs ./In/parameters_to-test_nside256_octree-ggg-omp

to plot results (vs reference results):

$ python python/compare_xi2pcf_curves.py --scale loglog --xscale radian --plot-mul-theta --ref b --file-a Output/histXi2pcf.txt --file-b Outputs_to_compare_with/Output_nside256_octree-ggg-omp/histXi2pcf.txt --outdir ./

and for 3pcf

$ python python/compare_xi3pcf_flatten_curves.py --file-a Output/ --file-b Outputs_to_compare_with/Output_nside256_octree-ggg-omp/ --ref b --bin-min 40 --bin-max 150

2. To test vs cute_box:

$ more test_cute_box/Readme.txt

to read the readme files in there.

For a local wrapper build, use the source root and matching profile:

$ make -j4 cballs cyballs-static-lib
$ CBALLS_STATIC_LIBRARY_READY=1 python3 setup.py build_ext --inplace --force
$ python3 -c "from cyballs import cballs; print(cballs)"
