

First download convergence maps from DESY3 (DES Y3 mass map reconstruction — https://arxiv.org/abs/2105.13539) webpage:

https://des.ncsa.illinois.edu/releases/y3a2/Y3massmaps

$ wget -c --tries=5 --timeout=30 --read-timeout=30 \
       --retry-connrefused --show-progress \
       -O KS_tomo1.fits https://desdr-server.ncsa.illinois.edu/despublic/y3a2_files/massmaps/KS_tomo1.fits

$ wget -c --tries=5 --timeout=30 --read-timeout=30 \
       --retry-connrefused --show-progress \
       -O glimpse_mask.fits https://desdr-server.ncsa.illinois.edu/despublic/y3a2_files/massmaps/glimpse_mask.fits

$ time python python/kappa_2pcf_with-plotting.py --fits catalogs_working/DES/KS_tomo2.fits --mask catalogs_working/DES/glimpse_mask.fits --outdir Output_DES --threads 16 --y-mul-theta --xscale arcmin --scale loglog

To see how kappa_2pcf_with-plotting.py execute:

$ python python/kappa_2pcf_with-plotting.py

or

$ python python/kappa_2pcf_with-plotting.py --help

If a Takahashi realization has been downloaded, then you can apply DES mask:

$ time python python/kappa_2pcf_with-plotting.py --fits catalogs_working/Takahashi/allskymap_nres12r081_zs9_mag.fits --mask catalogs_working/DES/glimpse_mask.fits --outdir Output_Takahashi --threads 16 --y-mul-theta --xscale arcmin --scale loglog --nside-down 1024

And then compare versus DES KS_tomo2 result:

$ python python/compare_xi2pcf_curves.py --outdir compare --xscale arcmin --scale loglog --file-a Output_DES/wtheta_KS_tomo2_cBalls.txt --file-b Output_Takahashi/wtheta_allskymap_nres12r081_zs9_mag_cBalls.txt --ref b --plot-mul-theta

Python script compare_xi2pcf_curves.py has help as kappa_2pcf_with-plotting.py:

$ python python/compare_xi2pcf_curves.py --help


Note: when kappa_2pcf_with-plotting.py is run it loads cballys module. If TPCFON = 1 (active) in addons/Makefile_addons_settings the 3pcf will be computed along with 2pcf. If you are only interested in 2pcf then set TPCFON = 0 and recompile cBalls again (make clean; make all). You will notice a speed up in the computation and in the folder results you will see only files related to 2pcf.



