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

        //B This case is done in search_octree_kkk_omp routine
//    case KKKCORRELATION:
//        break;
        //E

    case KKCORRELATION:
        if (scanopt(cmd->options, "compute-HistN")) PrintHistNN(cmd, gd);
            PrintHistXi2pcf(cmd, gd);
        break;
    case NNCORRELATION:
        if (scanopt(cmd->options, "compute-HistN")) PrintHistNN(cmd, gd);
        PrintHistN2pcf(cmd, gd);
        break;
    case NNEstimator:
        if (scanopt(cmd->options, "compute-HistN")) PrintHistNN(cmd, gd);
        PrintHistN2pcf(cmd, gd);
        break;
}
    break;

#endif	// ! _cballs_print_octree_kkk_omp_h
