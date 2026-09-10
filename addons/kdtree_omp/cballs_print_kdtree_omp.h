// Use:
//#include "cballs_kdtree_omp.h"

// NMultipoles has been switched off for kdtree_omp
//  NMultipoles -> NMultipoles_kdtree

#ifndef _cballs_print_kdtree_omp_h
#define _cballs_print_kdtree_omp_h

//#define KDTREEOMP         59

#ifdef KDTREEOMP
        case KDTREEOMPMETHOD:
#endif
#ifdef KDTREEMPI
        case KDTREEMPIMETHOD:
#endif
            verb_print(cmd->verbose,
                       "\n\tevalHist: printing %s method\n\n",
                       cmd->searchMethod);
            if (!cballs_opt_only_3pcf(cmd) && cballs_opt_compute_histn(cmd))
                PRINT_OR_FAIL(PrintHistNN(cmd, gd));
                PRINT_OR_FAIL(PrintHistrBins(cmd, gd));
            if (!cballs_opt_only_3pcf(cmd))
                PRINT_OR_FAIL(PrintHistXi2pcf(cmd, gd));

#ifdef TPCF
            if (!cballs_opt_only_2pcf(cmd))
                PRINT_OR_FAIL(PrintHistZetaM_sincos(cmd, gd));
                if (!cballs_opt_only_2pcf(cmd) && cballs_opt_out_m_histzeta(cmd)) {
                    PRINT_OR_FAIL(PrintHistZetaMm_sincos(cmd, gd));
                }
                if (!cballs_opt_only_2pcf(cmd) && cballs_opt_out_histzetag(cmd)) {
                    PRINT_OR_FAIL(PrintHistZetaMZetaGm_sincos(cmd, gd));
                }
#endif
            break;

#endif	// ! _cballs_print_kdtree_omp_h
