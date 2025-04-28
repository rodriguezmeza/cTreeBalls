
Takahasi et al

Full-sky Gravitational Lensing Mock Catalogs

http://cosmo.phys.hirosaki-u.ac.jp/takahasi/allsky_raytracing/

Download a realization (redshift 0.5078, zs9, r081, size 3.1 GB):

$ cd tests/catalogs/Takahasi
$ wget http://cosmo.phys.hirosaki-u.ac.jp/takahasi/allsky_raytracing/sub2/nres12/allskymap_nres12r081.zs9.mag.dat

Post processing files:

Extract only kappa:

Edit takahasi_to_fits.py to have the input file name as: allskymap_nres12r081.zs9.mag.dat

$ python takahasi_to_fits.py

This code above is from Tahahasi's web page. There is also f90 version code.

Now, experiment this convergence catalog with cBalls:
(we assume have installed Healpix)

$ cballs in=allskymap_nres12r081_zs9_mag.fits infmt=fits options=header-info,stop
$ map2gif -inp allskymap_nres12r081_zs9_mag.fits -out allskymap_nres12r081_zs9_mag.gif
$ ud_grade
$ mv outmap.fits allskymap_nres08r081_zs9_mag.fits
$ map2gif -inp allskymap_nres11r081_zs9_mag.fits -out allskymap_nres11r081_zs9_mag.gif


Down grade above ($ ud_grade) was done from Nside 1024 to 256 therefore the final name is:

allskymap_nres08r081_zs9_mag.fits

where nres08 means Nside 256. You may try another Nside's, like 1024, for instance.

Now, go tests directory

$ cd ../..




This catalog is not available due to GitHub space limitations...


Catalog of a Takahasi realization zs9 and r081 to side = 1024

full_sky_whole_XYZK__zs9_r081_nside1024.txt

and converted to binary format with:

cballs in=./Takahasi/full_sky_whole_XYZK__zs9_r081_nside1024.txt o=full_sky_whole_XYZK__zs9_r081_nside1024.bin ofmt=binary options=stop



Note: cballs add ".txt" extension to the file name. We removed it.

