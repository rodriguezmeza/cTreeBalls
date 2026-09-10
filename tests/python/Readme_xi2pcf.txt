cTreeBalls 2PCF comparison notes
================================

Use the active-profile scalar driver from any working directory:

    python3 tests/python/kappa_corr_all_engines.py \
        --fits DES/KS_tomo2.fits \
        --mask DES/glimpse_mask.fits \
        --engine octree-2balls-omp,kdtree-2balls-omp,balltree-2balls-omp \
        --threads 16 --outdir simulations/des

The driver reads the catalog once, passes it to each native engine in memory,
writes timing information, and produces ordinary and flattened comparison
plots. Use ``--list-engines`` to inspect the methods registered by the active
build profile.
