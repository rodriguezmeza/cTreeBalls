
In addons/Makefile_addons_settings do:

SMOOTHPIVOTON = 1
TPCFON = 1
PRUNEON = 1
THETA = 1.25

and in cTreeBalls directory:

$ make clean; make all

Then in directory "tests"

$ time python python/kappa_corr.py --fits catalogs/Takahashi/allskymap_nres12r081_zs9_mag.fits --threads 16

$ python python/compare_xi2pcf_curves.py --outdir ./ --scale loglog --xscale arcmin --plot-mul-theta --file-a Output/histXi2pcf.txt --file-b Outputs_to_compare_with/Output_nside4096_octree-ggg-omp_NMultipoles_NONORMHIST/histXi2pcf.txt --ref b

$ python python/compare_xi3pcf_flatten_curves.py --file-a Output --bin-min 100 --bin-max 400 --scale semilogy --file-b Outputs_to_compare_with/Output_nside4096_octree-ggg-omp_NMultipoles_NONORMHIST/ --ref b


