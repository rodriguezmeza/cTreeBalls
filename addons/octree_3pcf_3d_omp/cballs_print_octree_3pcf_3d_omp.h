// Use inside PrintEvalHist switch in source/cballs.c

#ifndef _cballs_print_octree_3pcf_3d_omp_h
#define _cballs_print_octree_3pcf_3d_omp_h

#ifdef OCTREE3PCF3DOMP
case 166:
#endif
#ifdef OCTREE3PCF3DMPI
case 192:
#endif
    /* The 3D estimator writes its own multipole table inside
     * searchcalc_octree_3pcf_3d_omp(). Avoid printing the legacy 2D
     * sin/cos/Chebyshev histograms for this mode. */
    if (PrintHistrBins(cmd, gd) == FAILURE) return FAILURE;
    break;

#endif // !_cballs_print_octree_3pcf_3d_omp_h
