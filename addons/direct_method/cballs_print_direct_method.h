// Use:
//#include "cballs_print_direct_method.h"

#ifndef _cballs_print_direct_method_h
#define _cballs_print_direct_method_h

//#define DIRECTSINCOS      19

case 19:
    verb_print(cmd->verbose, 
    "\n\tprintEvalHist: printing direct method (base sincos)\n\n");
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

#endif	// ! _cballs_print_direct_method_h
