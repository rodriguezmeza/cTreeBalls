// Use:
//#include "cballs_print_octree_kkk_omp.h"

#ifndef _cballs_print_octree_kkk_omp_h
#define _cballs_print_octree_kkk_omp_h

//#define OCTREEKKKOMP     61

case 61:
    verb_print(cmd->verbose,
               "\n\tprintEvalHist: printing octree-kkk-omp method...\n\n");
    PrintHistrBins(cmd, gd);

switch(correlation_int) {
    case KKCORRELATION:
        if (cballs_opt_compute_histn(cmd)) PrintHistNN(cmd, gd);
            PrintHistXi2pcf(cmd, gd);
        break;
    case NNCORRELATION:
        if (cballs_opt_compute_histn(cmd)) PrintHistNN(cmd, gd);
        PrintHistN2pcf(cmd, gd);
        break;
    case NNEstimator:
        if (cballs_opt_compute_histn(cmd)) PrintHistNN(cmd, gd);
        PrintHistN2pcf(cmd, gd);
        break;
}
    break;

#endif	// ! _cballs_print_octree_kkk_omp_h
