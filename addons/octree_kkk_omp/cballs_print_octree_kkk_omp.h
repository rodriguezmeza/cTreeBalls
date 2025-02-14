// Use:
//#include "cballs_print_octree_kkk_omp.h"

#ifndef _cballs_print_octree_kkk_omp_h
#define _cballs_print_octree_kkk_omp_h

//#define OCTREEKKKOMPMETHOD     61

case 61:
    verb_print(cmd->verbose,
               "\n\tprintEvalHist: printing octree-kkk-omp method...\n\n");
    PrintHistrBins(cmd, gd);

switch(correlation_int) {
    case KKKCORRELATION:
//B (computeTPCF)
#ifdef NONORMHIST
        if (scanopt(cmd->options, "no-normalize-HistZeta"))
            PrintHistZetaM_sincos(cmd, gd);
        else
            PrintHistZetaM_sincos_normalized(cmd, gd);
#else
        PrintHistZetaM_sincos(cmd, gd);
#endif
#ifdef NMultipoles
        PrintHistZetaM_sincos_N(cmd, gd);
#endif
        if (scanopt(cmd->options, "out-m-HistZeta")) {
#ifdef NONORMHIST
            if (scanopt(cmd->options, "no-normalize-HistZeta"))
                PrintHistZetaMm_sincos(cmd, gd);
            else
                PrintHistZetaMm_sincos_normalized(cmd, gd);
#else
            PrintHistZetaMm_sincos(cmd, gd);
#endif
#ifdef NMultipoles
            PrintHistZetaMm_sincos_N(cmd, gd);
#endif
        }
        if (scanopt(cmd->options, "out-HistZetaG")) {
            PrintHistZetaGm_sincos(cmd, gd);
        }
//E (computeTPCF)
        break;
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

#endif	// ! _cballs_print_octree_kkk_omp_h
