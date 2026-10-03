#ifndef _cmdline_defs_lya_forest_omp_h
#define _cmdline_defs_lya_forest_omp_h

    "lya2RpMax=200.0",       ";Maximum absolute line-of-sight separation for the Lyman-alpha 2PCF",
    "lya2RtMax=200.0",       ";Maximum transverse separation for the Lyman-alpha 2PCF",
    "lya2RpBins=50",         ";Number of line-of-sight bins for the Lyman-alpha 2PCF",
    "lya2RtBins=50",         ";Number of transverse bins for the Lyman-alpha 2PCF",
    "lya3RMax=80.0",         ";Maximum side length for the anisotropic Lyman-alpha 3PCF",
    "lya3RBins=20",          ";Number of side-length bins for the anisotropic Lyman-alpha 3PCF",
    "lya3ThetaBins=10",      ";Number of [0,pi] line-of-sight angle bins for the Lyman-alpha 3PCF",
    "lya3MuBins=20",         ";Number of [-1,1] opening-angle cosine bins for the Lyman-alpha 3PCF",

    "lya3Kernel=0",
    ";3PCF kernel: 0 segments, 1 reference, 2 tiled, 3 persistent nodes, 4 pivot cells, 5 hierarchical radial moments; multipole method: 0 auto, 1 legacy prefix, 2 forced hierarchy",

    "lya3LMax=8",
    ";Maximum anisotropic Legendre order (0..32), multipole method only",

    "lya3MuMode=0", ";Multipole mu output: 0 finite-L reconstruction only; 1 also compute exact hard bins with shared discovery",

    "lya3PivotBlock=0", ";3PCF pivot block: 0 automatic, or 1..4096; independent of thread count",

    "lya3MuSlop=0", ";Allowed opening-cosine bin leakage in bin-width units (0 exact, 0..1)",

    "lya3RadialSlop=0", ";Allowed radial bin leakage in bin-width units; kernels 3/4/5 only",

    "lya3PolarSlop=0", ";Allowed polar-angle bin leakage in bin-width units; kernels 3/4/5 only",

    "lya3PivotCellMax=8", ";Maximum pixels per pivot-cell task for kernels 4/5 (1..64)",

    "lya2Kernel=0", ";2PCF: 0 reference pixel traversal; 1 persistent forest cell pairs",
    "lya2RpSlop=0", ";Parallel bin leakage in bin-width units; lya2Kernel=1 only",
    "lya2RtSlop=0", ";Transverse bin leakage in bin-width units; lya2Kernel=1 only",

    "lyaScanLevel=0", ";Exact octree pivot frontier level (0 catalog order, 1..20 spatial tasks)",
    "lyaPivotRadius=0", ";Opt-in approximate same-forest pivot smoothing radius in catalog units (0 off)",
    "lyaPivotMax=8", ";Maximum original pixels per smoothed pivot representative (1..1024)",

#endif
