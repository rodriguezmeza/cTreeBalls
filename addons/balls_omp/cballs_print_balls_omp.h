// Use:
//#include "cballs_print_balls_omp.h"

#ifndef _cballs_print_balls_omp_h
#define _cballs_print_balls_omp_h

//#define BALLSOMPMETHOD         46

#ifdef BALLS
        case 46:
            verb_print(cmd->verbose,
                "\n\tprintEvalHist: printing  balls tree-omp method\n\n");
            if (scanopt(cmd->options, "compute-HistN")) PrintHistNN(cmd, gd);
            PrintHistrBins(cmd, gd);
            PrintHistXi2pcf(cmd, gd);
if (cmd->computeTPCF) {
#ifdef SINCOS
    PrintHistZetaM_sincos(cmd, gd);
    if (scanopt(cmd->options, "out-m-HistZeta"))
        PrintHistZetaMm_sincos(cmd, gd);
    if (scanopt(cmd->options, "out-HistZetaG")) {
        PrintHistZetaMZetaGm_sincos(cmd, gd);
    }
#else // This segment is obsolete... delete
    PrintHistZetaM(cmd, gd);
    if (scanopt(cmd->options, "out-m-HistZeta"))
        PrintHistZeta_theta2_fix(cmd, gd);
#endif // ! SINCOS
}
            break;

#endif // ! BALLS


#endif	// ! _cballs_print_balls_omp_h
