// Use:
//#include "cballs_print_balls_omp.h"

#ifndef _cballs_print_balls_omp_h
#define _cballs_print_balls_omp_h

//#define BALLSOMPMETHOD         46

#ifdef BALLS
        case 46:
            verb_print(cmd->verbose,
                "\n\tprintEvalHist: printing  balls tree-omp method\n\n");
            if (cballs_opt_compute_histn(cmd)) PrintHistNN(cmd, gd);
            PrintHistrBins(cmd, gd);
            PrintHistXi2pcf(cmd, gd);
#ifdef TPCF
#ifdef SINCOS
    PrintHistZetaM_sincos(cmd, gd);
    if (cballs_opt_out_m_histzeta(cmd))
        PrintHistZetaMm_sincos(cmd, gd);
    if (cballs_opt_out_histzetag(cmd)) {
        PrintHistZetaMZetaGm_sincos(cmd, gd);
    }
#else // This segment is obsolete... delete
    PrintHistZetaM(cmd, gd);
    if (cballs_opt_out_m_histzeta(cmd))
        PrintHistZeta_theta2_fix(cmd, gd);
#endif // ! SINCOS
#endif
            break;

#endif // ! BALLS


#endif	// ! _cballs_print_balls_omp_h
