// Use:
//#include "cballs_octree_ggg_omp_triangles.h"

//
// included in file: addons/addons_include/source/cballs/cballs_include02.h
//

#ifndef _cballs_octree_ggg_omp_triangles_h
#define _cballs_octree_ggg_omp_triangles_h

//
// look search routines tag numbers in file:
//              addons/addons_include/source/startrun/startrun_include_11.h
//

//#define OCTREEGGGOMPTRIANGLES     68

case 68:                   // search=octree-ggg-omp-triangles
    verb_print(cmd->verbose,
    "\n\tevalHist: with octree-ggg-omp-triangles method\n\n");
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
            searchcalc_octree_ggg_omp_triangles(cmd, gd, bodytable,
                                                gd->nbodyTable, 1,
                                      gd->nbodyTable,
                                      gd->iCatalogs[0], gd->iCatalogs[1]);
            break;
        default:
            searchcalc_octree_ggg_omp_triangles(cmd, gd, bodytable,
                                                gd->nbodyTable, 1,
                                     gd->nbodyTable,
                                     gd->iCatalogs[0], gd->iCatalogs[1]);
            break;
    }
    break;

#endif	// ! _cballs_octree_ggg_omp_triangles_h
