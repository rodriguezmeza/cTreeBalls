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
            if (cballs_opt_compute_histn(cmd)) PRINT_OR_FAIL(PrintHistNN(cmd, gd));
                PRINT_OR_FAIL(PrintHistrBins(cmd, gd));
            PRINT_OR_FAIL(PrintHistXi2pcf(cmd, gd));

#ifdef TPCF
            PRINT_OR_FAIL(PrintHistZetaM_sincos(cmd, gd));
                if (cballs_opt_out_m_histzeta(cmd)) {
                    PRINT_OR_FAIL(PrintHistZetaMm_sincos(cmd, gd));
                }
                if (cballs_opt_out_histzetag(cmd)) {
                    PRINT_OR_FAIL(PrintHistZetaMZetaGm_sincos(cmd, gd));
                }
#endif
            break;

#endif	// ! _cballs_print_kdtree_omp_h
