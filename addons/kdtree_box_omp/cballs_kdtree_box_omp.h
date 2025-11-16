// Use:
//#include "cballs_kdtree_box_omp.h"

// included in addons_include/source/cballs/cballs_include_02.h

#ifndef _cballs_kdtree_box_omp_h
#define _cballs_kdtree_box_omp_h

//#define KDTREEBOXOMP         71

    case 71:
        verb_print(cmd->verbose,
                   "\n\tevalHist: with kdtree-box-omp method\n\n");
        searchcalc_kdtree_box_omp(cmd, gd,
                              bodytable, gd->nbodyTable,
                              1, gd->nbodyTable,
                              gd->iCatalogs[0], gd->iCatalogs[1]);
        break;

#endif	// ! _cballs_kdtree_box_omp_h
