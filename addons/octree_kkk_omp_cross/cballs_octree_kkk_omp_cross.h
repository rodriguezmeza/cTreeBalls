// Use:
//#include "cballs_octree_ggg_omp.h"

//
// included in file: addons/addons_include/source/cballs/cballs_include02.h
//

#ifndef _cballs_octree_ggg_omp_h
#define _cballs_octree_ggg_omp_h

//
// look search routines tag numbers in file:
//              addons/addons_include/source/startrun/startrun_include_11.h
//

//#define OCTREEGGGOMP     66

case 66:                   // search=octree-ggg-omp
    verb_print(cmd->verbose,
    "\n\tevalHist: with octree-ggg-omp method\n\n");
    for (ifile=0; ifile<gd->ninfiles; ifile++) {
        DO_BODY(p,bodytable[ifile],bodytable[ifile]+gd->nbodyTable[ifile]) {
            Update(p) = TRUE;
        }
        if (MakeTree(cmd, gd, bodytable[ifile], gd->nbodyTable[ifile], ifile)
            == FAILURE)
            return FAILURE;
    }
    switch(correlation_int) {
        case GGGCORRELATION:
            searchcalc_octree_ggg_omp(cmd, gd, bodytable, gd->nbodyTable, 1,
                                      gd->nbodyTable,
                                      gd->iCatalogs[0], gd->iCatalogs[1]);
            break;
        default:
            searchcalc_octree_ggg_omp(cmd, gd, bodytable, gd->nbodyTable, 1,
                                     gd->nbodyTable,
                                     gd->iCatalogs[0], gd->iCatalogs[1]);
            break;
    }
    break;

#endif	// ! _cballs_octree_ggg_omp_h
