// Use:
//#include "cballs_kdtree_omp.h"

#ifndef _cballs_print_kdtree_omp_h
#define _cballs_print_kdtree_omp_h

//#define KDTREEOMP         59

        case 59:
            verb_print(cmd->verbose,
                       "\n\tevalHist: printing kdtree-omp method\n\n");
            if (scanopt(cmd->options, "compute-HistN")) PrintHistNN(cmd, gd);
                PrintHistrBins(cmd, gd);
            PrintHistXi2pcf(cmd, gd);
#ifdef NMultipoles
            PrintHistXi2pcf_N(cmd, gd);
#endif
            if (cmd->computeTPCF) {
                PrintHistZetaM_sincos(cmd, gd);
#ifdef NMultipoles
                PrintHistZetaM_sincos_N(cmd, gd);
#endif
                if (scanopt(cmd->options, "out-m-HistZeta")) {
                    PrintHistZetaMm_sincos(cmd, gd);
#ifdef NMultipoles
                    PrintHistZetaMm_sincos_N(cmd, gd);
#endif
                }
                if (scanopt(cmd->options, "out-HistZetaG")) {
                    PrintHistZetaMZetaGm_sincos(cmd, gd);
                }
            }
            break;

#endif	// ! _cballs_print_kdtree_omp_h
