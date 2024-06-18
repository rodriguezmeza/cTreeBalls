// Use:
//NOLSST:
//#include "tpcf_print_tree_3pcf_direct_omp.h"

#ifndef _tpcf_print_tree_3pcf_direct_omp_h
#define _tpcf_print_tree_3pcf_direct_omp_h

//#define TREE3PCFBFOMPMETHOD     45

case 45:
    verb_print(cmd->verbose,
"\n\tprintEvalHist: printing normal tree method (tree-3pcf-direct-omp)\n\n");
    if (scanopt(cmd->options, "compute-HistN")) printHistN(cmd, gd);
    printHistrBins(cmd, gd);
    printHistXi2pcf(cmd, gd);
if (cmd->computeTPCF) {
    printHistZetaM(cmd, gd);
    if (scanopt(cmd->options, "out-m-HistZeta"))
        printHistZeta(cmd, gd);
}
    break;

#endif	// ! _tpcf_print_tree_3pcf_direct_omp_h
