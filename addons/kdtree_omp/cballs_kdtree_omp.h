// Use:
//#include "cballs_kdtree_omp.h"

#ifndef _cballs_kdtree_omp_h
#define _cballs_kdtree_omp_h

//#define KDTREEOMP         59

    case 59:
        verb_print(cmd->verbose,
                   "\n\tevalHist: with kdtree-omp method\n\n");
        searchcalc_kdtree_omp(cmd, gd,
                              bodytable, gd->nbodyTable,
                              1, gd->nbodyTable,
                              gd->iCatalogs[0], gd->iCatalogs[1]);
        break;

#endif	// ! _cballs_kdtree_omp_h
