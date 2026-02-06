// Use:
//#include "cballs_kdtree_omp.h"

// NMultipoles has been switched off for kdtree_omp
//  NMultipoles -> NMultipoles_kdtree

#ifndef _cballs_print_kdtree_omp_h
#define _cballs_print_kdtree_omp_h

//#define KDTREEOMP         59

        case 59:
            verb_print(cmd->verbose,
                       "\n\tevalHist: printing kdtree-omp method\n\n");
            if (scanopt(cmd->options, "compute-HistN")) PrintHistNN(cmd, gd);
                PrintHistrBins(cmd, gd);
            PrintHistXi2pcf(cmd, gd);

#ifdef TPCF
                PrintHistZetaM_sincos(cmd, gd);
                if (scanopt(cmd->options, "out-m-HistZeta")) {
                    PrintHistZetaMm_sincos(cmd, gd);
                }
                if (scanopt(cmd->options, "out-HistZetaG")) {
                    PrintHistZetaMZetaGm_sincos(cmd, gd);
                }
#endif
            break;

#endif	// ! _cballs_print_kdtree_omp_h
