// Use:
//#include "cballs_octree_sincos_omp.h"

// included in addons_include/source/cballs/cballs_include_03.h

#ifndef _cballs_print_octree_sincos_omp_h
#define _cballs_print_octree_sincos_omp_h

//#define OCTREESINCOSOMP         74

        case 74:
            verb_print(cmd->verbose,
                    "\n\tevalHist: printing octree-sincos-omp method\n\n");
            if (scanopt(cmd->options, "compute-HistN"))
                PrintHistNN(cmd, gd);
            PrintHistrBins(cmd, gd);
            PrintHistXi2pcf(cmd, gd);
            if (cmd->computeTPCF) {
                PrintHistZetaM_sincos(cmd, gd);
                if (scanopt(cmd->options, "out-m-HistZeta"))
                    PrintHistZetaMm_sincos(cmd, gd);
                if (scanopt(cmd->options, "out-HistZetaG")) {
                    PrintHistZetaMZetaGm_sincos(cmd, gd);
                }
            }
            break;

#endif	// ! _cballs_print_octree_sincos_omp_h
