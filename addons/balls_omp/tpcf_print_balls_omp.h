// Use:
//NOLSST:
//#include "tpcf_print_balls_omp.h"

#ifndef _tpcf_print_balls_omp_h
#define _tpcf_print_balls_omp_h

//#define BALLSOMPMETHOD         46

#ifdef BALLS
        case 46:
            verb_print(cmd->verbose,
                "\n\tprintEvalHist: printing  balls tree-omp method\n\n");
            if (scanopt(cmd->options, "compute-HistN")) printHistN(cmd, gd);
            printHistrBins(cmd, gd);
            printHistXi2pcf(cmd, gd);
//#ifdef TPCF
if (cmd->computeTPCF) {
#ifdef SINCOS
    printHistZetaM_sincos(cmd, gd);
    if (scanopt(cmd->options, "out-m-HistZeta"))
        printHistZeta_sincos(cmd, gd);
#else
    printHistZetaM(cmd, gd);
    if (scanopt(cmd->options, "out-m-HistZeta"))
        printHistZeta(cmd, gd);
#endif // ! SINCOS
//#endif // ! TPCF
}
            break;

#endif // ! BALLS


#endif	// ! _tpcf_print_balls_omp_h
