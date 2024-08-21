// Use:
//#include "tpcf_print_direct_method.h"

#ifndef _tpcf_print_direct_method_h
#define _tpcf_print_direct_method_h

//#define DIRECTSIMPLESINCOS      19

case 19:
    verb_print(cmd->verbose, 
    "\n\tprintEvalHist: printing direct method simple (base sincos)\n\n");
    if (scanopt(cmd->options, "compute-HistN")) PrintHistNN(cmd, gd);
    PrintHistrBins(cmd, gd);
    PrintHistXi2pcf(cmd, gd);
    if (cmd->computeTPCF) {
        PrintHistZetaM_sincos(cmd, gd);
        if (scanopt(cmd->options, "out-m-HistZeta"))
            PrintHistZetaMm_sincos(cmd, gd);
        if (scanopt(cmd->options, "out-HistZetaG")) {
            PrintHistZetaMZetaGm_sincos(cmd, gd);
        }
    }
    break;

#endif	// ! _tpcf_print_direct_method_h
