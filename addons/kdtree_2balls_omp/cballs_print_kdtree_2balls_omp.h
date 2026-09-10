#ifndef _cballs_print_kdtree_2balls_omp_h
#define _cballs_print_kdtree_2balls_omp_h

#ifdef KDTREE2BALLSOMP
        case KDTREE2BALLSOMPMETHOD:
#endif
#ifdef KDTREE2BALLSMPI
        case KDTREE2BALLSMPIMETHOD:
#endif
            verb_print(cmd->verbose,
                       "\n\tprintEvalHist: printing %s method\n\n",
                       cmd->searchMethod);
#ifdef TWOPCF
            if (!cballs_opt_only_3pcf(cmd)) {
                if (cballs_opt_compute_histn(cmd))
                    PRINT_OR_FAIL(PrintHistNN(cmd, gd));
                PRINT_OR_FAIL(PrintHistXi2pcf(cmd, gd));
            }
#endif
            PRINT_OR_FAIL(PrintHistrBins(cmd, gd));
#ifdef THREEPCFCONVERGENCE
            if (!cballs_opt_only_2pcf(cmd)) {
                PRINT_OR_FAIL(PrintHistZetaM_sincos(cmd, gd));
                if (cballs_opt_out_m_histzeta(cmd))
                    PRINT_OR_FAIL(PrintHistZetaMm_sincos(cmd, gd));
            }
#endif
            break;

#endif /* !_cballs_print_kdtree_2balls_omp_h */
