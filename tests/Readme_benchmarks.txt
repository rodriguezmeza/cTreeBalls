
First, download a Takahashi's realization (http://cosmo.phys.hirosaki-u.ac.jp/takahasi/allsky_raytracing/)

$ cd tests/catalogs/Takahashi
$ wget http://cosmo.phys.hirosaki-u.ac.jp/takahasi/allsky_raytracing/sub2/nres12/allskymap_nres12r081.zs9.mag.dat

In addons/Makefile_addons_settings do:

SMOOTHPIVOTON = 1
TPCFON = 0

and in cTreeBalls directory:

$ make clean; make all

Then in directory "tests"

$ time python3 python/benchmark_kappa_corr.py --kappa catalogs/Takahashi/allskymap_nres12r081_zs9_mag.fits --nsides "128,256,512,1024,2048" --threads 16 --nbins 20 --repeats 10 --outdir cputime_2pcf

Now, do in addons/Makefile_addons_settings:

TPCFON = 1

and in cTreeBalls directory:

$ make clean; make all

and in tests folder:

$ time python3 python/benchmark_kappa_corr.py --kappa catalogs/Takahashi/allskymap_nres12r081_zs9_mag.fits --nsides "128,256,512,1024,2048" --threads 16 --nbins 20 --repeats 10 --outdir cputime_3pcf


Note: comparisons against other correlation codes can be found in tests/comparisons.


