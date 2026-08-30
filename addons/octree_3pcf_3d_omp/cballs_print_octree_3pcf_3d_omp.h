// Use inside PrintEvalHist switch in source/cballs.c

#ifndef _cballs_print_octree_3pcf_3d_omp_h
#define _cballs_print_octree_3pcf_3d_omp_h

case 166:                   // search=octree-3pcf-3d-omp
    /* The 3D estimator writes its own multipole table inside
     * searchcalc_octree_3pcf_3d_omp(). Avoid printing the legacy 2D
     * sin/cos/Chebyshev histograms for this mode. */
    PrintHistrBins(cmd, gd);
    break;

#endif // !_cballs_print_octree_3pcf_3d_omp_h
