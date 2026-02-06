
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

$ cd test_cute_box

Read the readme files.




