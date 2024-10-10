// Use:
//#include "cballs_print_octree_kkk_omp.h"

#ifndef _cballs_print_octree_kkk_omp_h
#define _cballs_print_octree_kkk_omp_h

//#define OCTREEKKKOMPMETHOD     61

case 61:
    verb_print(cmd->verbose,
               "\n\tprintEvalHist: printing tc method (octree-kkk-omp)\n\n");
    PrintHistrBins(cmd, gd);
#ifdef NONORMHIST
    if (scanopt(cmd->options, "normalize-HistZeta"))
        PrintHistZetaM_sincos_normalized(cmd, gd);
    else
        PrintHistZetaM_sincos(cmd, gd);
#else
    PrintHistZetaM_sincos(cmd, gd);
#endif
#ifdef NMultipoles
    PrintHistZetaM_sincos_N(cmd, gd);
#endif
    if (scanopt(cmd->options, "out-m-HistZeta")) {
#ifdef NONORMHIST
        if (scanopt(cmd->options, "normalize-HistZeta"))
            PrintHistZetaMm_sincos_normalized(cmd, gd);
        else
            PrintHistZetaMm_sincos(cmd, gd);
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
    break;

#endif	// ! _cballs_print_octree_kkk_omp_h
