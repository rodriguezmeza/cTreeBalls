// Use:
//#include "cballs_octree_balls4_omp.h"

#ifndef _cballs_octree_balls4_omp_h
#define _cballs_octree_balls4_omp_h

//#define OCTREEBALLS4OMP     69

case OCTREEBALLS4OMPMETHOD: // search=octree-balls4-omp
    verb_print(cmd->verbose,
    "\n\tevalHist: with octree method (octree-balls4-omp)\n\n");
    for (ifile=0; ifile<gd->ninfiles; ifile++) {
        DO_BODY(p,bodytable[ifile],bodytable[ifile]+gd->nbodyTable[ifile])
        Update(p) = TRUE;
        if (MakeTree(cmd, gd, bodytable[ifile], gd->nbodyTable[ifile], ifile)
            == FAILURE)
            return FAILURE;
    }
    if (searchcalc_octree_balls4_omp(
            cmd, gd, bodytable, gd->nbodyTable, 1, gd->nbodyTable,
            gd->iCatalogs[0], gd->iCatalogs[1]) == FAILURE)
        return FAILURE;
    break;

#endif	// ! _cballs_octree_balls4_omp_h
