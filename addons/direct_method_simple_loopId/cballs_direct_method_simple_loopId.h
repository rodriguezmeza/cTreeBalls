// Use:
//#include "cballs_direct_method_simple_loopId.h"
//
//  it is included in:
//      addons/addons_include/source/cballs/cballs_include_02.h
//

#ifndef _cballs_direct_method_simple_loopId_h
#define _cballs_direct_method_simple_loopId_h

//
// see:
//  addons/addons_include/source/startrun/startrun_include_11.h
//  for the tag numbering
//

//#define DIRECTSIMPLE      78

case 78:
    verb_print(cmd->verbose,
               "\n\tEvalHist: simple direct method loopId (base sincos)\n\n");
    searchcalc_direct_simple_sincos_loopId(cmd, gd, bodytable, gd->nbodyTable, 1,
            gd->nbodyTable, gd->iCatalogs[0], gd->iCatalogs[1]);
    break;

#endif	// ! _cballs_direct_method_simple_loopId_h
