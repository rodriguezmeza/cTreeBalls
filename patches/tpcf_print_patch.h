// Use:
//NOLSST:
//#include "tpcf_print_patch.h"

#ifndef _tpcf_print_patch_h
#define _tpcf_print_patch_h


        case DIRECTSIMPLE:
            verb_print(cmd.verbose, "\n\tprintEvalHist: printing direct method simple\n\n");
            if (scanopt(cmd.options, "compute-HistN")) printHistN();
            printHistXi2pcf();
#ifdef TPCF
            printHistZetaM();
            printHistZeta();
#endif
            break;

#endif	// ! _tpcf_print_patch_h
