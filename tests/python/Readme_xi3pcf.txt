cTreeBalls scalar 3PCF notes
===========================

Build the active profile with ``TPCFON=1`` and run, from the checkout root:

    python3 tests/python/kappa_corr_all_engines.py \
        --fits tests/catalogs/allskymap_nres12r081_zs9_mag.fits \
        --engine octree-2balls-omp,kdtree-2balls-omp,balltree-2balls-omp \
        --statistics 3pcf --threads 16 --outdir Output_3pcf

The driver uses the same in-memory catalog for every selected engine and
writes radial-bin and flattened 3PCF comparison plots. For a strict numerical
check, use exact body traversal on a small catalog before enabling approximate
dual-node acceptance.
