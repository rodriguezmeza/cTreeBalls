// Use:
//#include "cballs_kdtree_cute_box.h"
//
//  it is included in:
//      addons/addons_include/source/cballs/cballs_include_02.h
//

#ifndef _cballs_kdtree_cute_box_h
#define _cballs_kdtree_cute_box_h

//
// see:
//  addons/addons_include/source/startrun/startrun_include_11.h
//  for the tag numbering
//

//#define KDTREECUTEBOX        48

    case 48:
        verb_print(cmd->verbose,
                   "\n\tevalHist: with kdtree-cute-box method\n\n");
        searchcalc_kdtree_cute_box(cmd, gd,
                                bodytable, gd->nbodyTable, 1, gd->nbodyTable,
                                gd->iCatalogs[0], gd->iCatalogs[1]);
        break;

#endif	// ! _cballs_kdtree_cute_box_h
