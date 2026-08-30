// Use:
//#include "cballs_octree_kkk_balls4_omp_triangles.h"

#ifndef _cballs_octree_kkk_balls4_omp_triangles_h
#define _cballs_octree_kkk_balls4_omp_triangles_h

//#define OCTREEKKKBALLS4OMPTRIANGLES     70

case 70:                   // search=octree-kkk-omp-triangles
    verb_print(cmd->verbose,
    "\n\tevalHist: with octree method (octree-kkk-balls4-omp-triangles)\n\n");
    for (ifile=0; ifile<gd->ninfiles; ifile++) {
        DO_BODY(p,bodytable[ifile],bodytable[ifile]+gd->nbodyTable[ifile])
        Update(p) = TRUE;
        if (MakeTree(cmd, gd, bodytable[ifile], gd->nbodyTable[ifile], ifile)
            == FAILURE)
            return FAILURE;
    }
    searchcalc_octree_kkk_balls4_omp_triangles(cmd, gd,
                                               bodytable, gd->nbodyTable, 1,
                              gd->nbodyTable,
                              gd->iCatalogs[0], gd->iCatalogs[1]);
    break;

#endif	// ! _cballs_octree_kkk_balls4_omp_triangles_h
