#ifndef _cballs_print_balltree_2balls_omp_3pcf_h
#define _cballs_print_balltree_2balls_omp_3pcf_h

        case BALLTREE2BALLSOMP3PCFMETHOD:
            verb_print(cmd->verbose,
                       "\n\tevalHist: printing TreeCorr-style LogMultipole 3PCF\n\n");
            PRINT_OR_FAIL(PrintHistrBins(cmd, gd));
            PRINT_OR_FAIL(PrintHistZetaM_sincos(cmd, gd));
            if (cballs_opt_out_m_histzeta(cmd))
                PRINT_OR_FAIL(PrintHistZetaMm_sincos(cmd, gd));
            break;

#endif /* !_cballs_print_balltree_2balls_omp_3pcf_h */
