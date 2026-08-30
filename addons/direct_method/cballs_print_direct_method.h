// Use:
//#include "cballs_print_direct_method.h"

#ifndef _cballs_print_direct_method_h
#define _cballs_print_direct_method_h

//#define DIRECTSINCOS      19

case 19:
    verb_print(cmd->verbose, 
    "\n\tprintEvalHist: printing direct method (base sincos)\n\n");
    if (cballs_opt_compute_histn(cmd)) PrintHistNN(cmd, gd);
    PrintHistrBins(cmd, gd);
    PrintHistXi2pcf(cmd, gd);
#ifdef TPCF
        PrintHistZetaM_sincos(cmd, gd);
        if (cballs_opt_out_m_histzeta(cmd))
            PrintHistZetaMm_sincos(cmd, gd);
        if (cballs_opt_out_histzetag(cmd)) {
            PrintHistZetaMZetaGm_sincos(cmd, gd);
        }
#endif
    break;

#endif	// ! _cballs_print_direct_method_h
