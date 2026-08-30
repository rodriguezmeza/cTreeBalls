// Use:
//#include "cballs_print_balls_omp_0357.h"

#ifndef _cballs_print_balls_omp_0357_h
#define _cballs_print_balls_omp_0357_h

//#define BALLSOMPMETHOD         79

#ifdef BALLS0357
        case 79:
            verb_print(cmd->verbose,
                "\n\tprintEvalHist: printing  balls tree-omp method\n\n");
            if (cballs_opt_compute_histn(cmd))
                PRINT_OR_FAIL(PrintHistNN(cmd, gd));
            PRINT_OR_FAIL(PrintHistrBins(cmd, gd));
            PRINT_OR_FAIL(PrintHistXi2pcf(cmd, gd));
if (gd->computeTPCF) {
#ifdef SINCOS
    PRINT_OR_FAIL(PrintHistZetaM_sincos(cmd, gd));
    if (cballs_opt_out_m_histzeta(cmd))
        PRINT_OR_FAIL(PrintHistZetaMm_sincos(cmd, gd));
    if (cballs_opt_out_histzetag(cmd)) {
        PRINT_OR_FAIL(PrintHistZetaMZetaGm_sincos(cmd, gd));
    }
#else // This segment is obsolete... delete
    PRINT_OR_FAIL(PrintHistZetaM(cmd, gd));
    if (cballs_opt_out_m_histzeta(cmd))
        PRINT_OR_FAIL(PrintHistZeta_theta2_fix(cmd, gd));
#endif // ! SINCOS
}
            break;

#endif // ! BALLS0357


#endif	// ! _cballs_print_balls_omp_0357_h
