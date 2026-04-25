
To process a Takahashi realization do:

$ time python python/kappa_corr.py --fits catalogs/Takahashi/allskymap_nres12r081_zs9_mag.fits --threads 16 --nside-down 1024

And the plot results:

$ python python/compare_xi2pcf_curves.py --outdir ./ --scale loglog --xscale arcmin --plot-mul-theta --file-a Output/histXi2pcf.txt

$ python python/compare_xi3pcf_flatten_curves.py --file-a Output --bin-min 1 --bin-max 400 --scale semilogy

Note: we were supposed to download first a Takahashi realization and that it was converted to HEALPix format. Go to catalogs/Takahashi and read the Readme file there.
