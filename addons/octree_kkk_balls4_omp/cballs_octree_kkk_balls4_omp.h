// Use:
//#include "cballs_octree_kkk_balls4_omp.h"

#ifndef _cballs_octree_kkk_balls4_omp_h
#define _cballs_octree_kkk_balls4_omp_h

//#define OCTREEKKKBALLS4OMP     69

case OCTREEKKKBALLS4OMPMETHOD: // search=octree-kkk-balls4-omp
    verb_print(cmd->verbose,
    "\n\tevalHist: with octree method (octree-kkk-balls4-omp)\n\n");
    for (ifile=0; ifile<gd->ninfiles; ifile++) {
        DO_BODY(p,bodytable[ifile],bodytable[ifile]+gd->nbodyTable[ifile])
        Update(p) = TRUE;
        if (MakeTree(cmd, gd, bodytable[ifile], gd->nbodyTable[ifile], ifile)
            == FAILURE)
            return FAILURE;
    }
    if (searchcalc_octree_kkk_balls4_omp(
            cmd, gd, bodytable, gd->nbodyTable, 1, gd->nbodyTable,
            gd->iCatalogs[0], gd->iCatalogs[1]) == FAILURE)
        return FAILURE;
    break;

#endif	// ! _cballs_octree_kkk_balls4_omp_h
