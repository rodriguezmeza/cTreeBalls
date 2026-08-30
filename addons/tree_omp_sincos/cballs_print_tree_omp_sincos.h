// Use:
//#include "cballs_print_tree_omp_sincos.h"

// included in addons_include/source/cballs/cballs_include_03.h

#ifndef _cballs_print_tree_omp_sincos_h
#define _cballs_print_tree_omp_sincos_h

//#define TREEOMPSINCOS             75

        case 75:
            verb_print(cmd->verbose,
                    "\n\tevalHist: printing octree-sincos-omp method\n\n");
            if (cballs_opt_compute_histn(cmd))
                PrintHistNN(cmd, gd);
            PrintHistrBins(cmd, gd);
            PrintHistXi2pcf(cmd, gd);
#ifdef TPCF
                PrintHistZetaM_sincos(cmd, gd);
                if (cballs_opt_out_m_histzeta(cmd))
                    PrintHistZetaMm_sincos(cmd, gd);
                if (cballs_opt_out_histzetag(cmd)) {
                    PrintHistZetaMZetaGm_sincos(cmd, gd);
                }
#endif
            break;

#endif	// ! _cballs_print_tree_omp_sincos_h
