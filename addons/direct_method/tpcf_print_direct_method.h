// Use:
//NOLSST:
//#include "tpcf_print_direct_method.h"

#ifndef _tpcf_print_direct_method_h
#define _tpcf_print_direct_method_h

//#define DIRECTSIMPLE            12
//#define DIRECTSIMPLESINCOS      19

case 12:
    verb_print(cmd->verbose,
               "\n\tprintEvalHist: printing direct method simple\n\n");
    if (scanopt(cmd->options, "compute-HistN")) printHistN(cmd, gd);
    printHistXi2pcf(cmd, gd);
#ifdef TPCF
    printHistZetaM(cmd, gd);
    printHistZeta(cmd, gd);
#endif
    break;

case 19:
    verb_print(cmd->verbose, 
    "\n\tprintEvalHist: printing direct method simple (base sincos)\n\n");
    if (scanopt(cmd->options, "compute-HistN")) printHistN(cmd, gd);
    printHistrBins(cmd, gd);
    printHistXi2pcf(cmd, gd);
#ifdef TPCF
    printHistZetaM_sincos(cmd, gd);
    if (scanopt(cmd->options, "out-m-HistZeta"))
    printHistZeta_sincos(cmd, gd);
#endif
    break;

#endif	// ! _tpcf_print_direct_method_h
