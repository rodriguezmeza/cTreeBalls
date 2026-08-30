// Use:
//#include "cballs_direct_method_simple.h"
//
//  it is included in:
//      addons/addons_include/source/cballs/cballs_include_02.h
//

#ifndef _cballs_direct_method_simple_h
#define _cballs_direct_method_simple_h

//
// see:
//  addons/addons_include/source/startrun/startrun_include_11.h
//  for the tag numbering
//

//#define DIRECTSIMPLE      67

case 67:
    verb_print(cmd->verbose,
               "\n\tEvalHist: simple direct method (base sincos)\n\n");
    searchcalc_direct_simple_sincos(cmd, gd, bodytable, gd->nbodyTable, 1,
            gd->nbodyTable, gd->iCatalogs[0], gd->iCatalogs[1]);
    break;

#endif	// ! _cballs_direct_method_simple_h
