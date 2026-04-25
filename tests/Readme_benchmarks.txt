
First, download a Takahashi's realization (http://cosmo.phys.hirosaki-u.ac.jp/takahasi/allsky_raytracing/)

$ cd tests/catalogs/Takahashi
$ wget http://cosmo.phys.hirosaki-u.ac.jp/takahasi/allsky_raytracing/sub2/nres12/allskymap_nres12r081.zs9.mag.dat

Transform this Takahashi's file to HEALPix format. First edit file takahasi_to_fits.py to adapt to the name of file you just downloaded, and then execute:

$ python takahasi_to_fits.py

In addons/Makefile_addons_settings do:

SMOOTHPIVOTON = 1
TPCFON = 0

and in cTreeBalls directory:

$ make clean; make all

Then in directory "tests"

$ time python3 python/benchmark_kappa_corr.py --kappa catalogs/Takahashi/allskymap_nres12r081_zs9_mag.fits --nsides "128,256,512,1024,2048" --threads 16 --nbins 20 --repeats 5 --outdir cputime_2pcf

Now, do in addons/Makefile_addons_settings:

TPCFON = 1

and in cTreeBalls directory:

$ make clean; make all

and in tests folder:

$ time python3 python/benchmark_kappa_corr.py --kappa catalogs/Takahashi/allskymap_nres12r081_zs9_mag.fits --nsides "128,256,512,1024,2048" --threads 16 --nbins 20 --repeats 5 --outdir cputime_3pcf


Note: in some machines "make all" does not install cyballs. Try:

$ python3.xx setup.py build
$ python3.xx setup.py install

and if last line does not work use:

$ python3.xx setup.py install --user

In general, we recommend the use of a python environment.
