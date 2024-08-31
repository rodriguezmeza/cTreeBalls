// Use:
//#include "cballs_direct_method.h"

#ifndef _cballs_direct_method_h
#define _cballs_direct_method_h

//#define DIRECTSIMPLESINCOS      19

case 19:
    verb_print(cmd->verbose,
               "\n\tEvalHist: simple direct method (base sincos)\n\n");
    searchcalc_direct_simple_sincos(cmd, gd, bodytable, gd->nbodyTable, 1,
            gd->nbodyTable, gd->iCatalogs[0], gd->iCatalogs[1]);
    break;

#endif	// ! _cballs_direct_method_h
