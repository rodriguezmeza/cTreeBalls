

- Individual computations

$ time python3 des_kappa_corr.py --fits DES/KS_tomo2.fits --outdir simulaciones/des/cballs --engine cballs --threads=16 --nbins 20 --y-mul-theta --scale loglog --xscale arcmin --mask DES/glimpse_mask.fits

$ time python3 des_kappa_corr.py --fits DES/KS_tomo2.fits --outdir simulaciones/des/treecorr --engine treecorr --threads=16 --nbins 20 --y-mul-theta --scale loglog --xscale arcmin --mask DES/glimpse_mask.fits

$ time python3 des_kappa_corr.py --fits DES/KS_tomo2.fits --outdir simulaciones/des/corrfunc --engine corrfunc --threads=16 --nbins 20 --y-mul-theta --scale loglog --xscale arcmin --mask DES/glimpse_mask.fits

- Comparison

$ python compare_des_curves.py --file-a simulaciones/des/treecorr/wtheta_KS_tomo2_treecorr.txt --label-a "TreeCorr" --file-b simulaciones/des/corrfunc/wtheta_KS_tomo2_corrfunc.txt --label-b "Corrfunc" --file-c simulaciones/des/cballs/wtheta_KS_tomo2_cballs.txt --label-c "cBalls" --ref a --theta-min 10 --theta-max 190 --scale loglog --outdir simulaciones/des --out "des_vs_treecorr" --plot-mul-theta


Comparison vs octree-ggg-omp old results


$ time python python/kappa_corr.py --fits Takahasi/allskymaps_fits/allskymap_nres10r081_zs9_mag.fits --outdir Output --xscale arcmin

$ python python/compare_xi2pcf_curves.py --file-a Output/histXi2pcf.txt --file-b Output/histXi2pcf.txt --file-c Outputs_to_compare_with/Output_nside1024_zs9r081_octree-ggg-omp_NMultipoles_NONORMHIST/histXi2pcf.txt --plot-mul-theta --outdir ./



Use after 24 line:

# perfect comparison vs Output_nside1024_zs9r081_octree-ggg-omp_NMultipoles_NONORMHIST
#   Using  16 threads...
#   Searching cputime= 230.66754774999998  sec.
#   cputime: 4m4.830s
#        Balls.set({'infile':'./Takahasi/full_sky_whole_XYZK__zs9_r081_nside1024.txt'})
#        Balls.set({'infileformat':'columns-ascii'})
#


Without any approximation, only PRUNEON = 1 (THETA=1, theta=1)

$ time python python/kappa_corr.py --fits Takahasi/allskymaps_fits/allskymap_nres12r081_zs9_mag.fits --outdir Output --threads 16

$ python python/compare_xi2pcf_curves.py --file-a Output/histXi2pcf.txt --label-a cBalls --file-b TreeCorr/Output_TC_kk/histXi2pcf.txt --label-b TreeCorr --file-c Outputs_to_compare_with/Output_nside4096_zs9r081_octree-ggg-omp_NMultipoles_NONORMHIST/histXi2pcf.txt --label-c cBalls-good --ref c --plot-mul-theta --outdir ./

Output -> Output_nside4096_zs9r081_octree-ggg-omp_NMultipoles_NONORMHIST_pruned_Xi2pcf



$ time python python/kappa_corr.py --fits Takahasi/allskymaps_fits/allskymap_nres12r081_zs9_mag.fits --threads 16


Last ones:

time python python/kappa_corr.py --fits Takahasi/allskymaps_fits/allskymap_nres10r081_zs9_mag.fits --threads 16

python python/compare_xi2pcf_curves.py --outdir ./ --scale loglog --xscale arcmin --plot-mul-theta --file-a Output/histXi2pcf.txt --file-b Outputs_to_compare_with_v1.0.1/Output_nside4096_octree-ggg-omp_NMultipoles_NONORMHIST_v1.0.1_pruned_theta-1p00_THETA-1p25/histXi2pcf.txt

python python/compare_xi3pcf_flatten_curves.py --file-a Output --bin-min 100 --bin-max 400 --scale semilogy --file-b Outputs_to_compare_with_v1.0.1/Output_nside4096_octree-ggg-omp_NMultipoles_NONORMHIST_v1.0.1_pruned_theta-1p05_THETA-1p00/


