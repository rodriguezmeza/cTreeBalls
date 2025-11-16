// Use:
//#include "cballs_octree_box_omp.h"

// included in addons_include/source/cballs/cballs_include_02.h

#ifndef _cballs_octree_box_omp_h
#define _cballs_octree_box_omp_h

//#define OCTREEBOXOMP         72

    case 72:
        verb_print(cmd->verbose,
                   "\n\tevalHist: with octree-box-omp method\n\n");
        for (ifile=0; ifile<gd->ninfiles; ifile++) {
            DO_BODY(p,bodytable[ifile],bodytable[ifile]+gd->nbodyTable[ifile])
            Update(p) = TRUE;
            MakeTree(cmd, gd, bodytable[ifile], gd->nbodyTable[ifile], ifile);
        }
        searchcalc_octree_box_omp(cmd, gd,
                              bodytable, gd->nbodyTable,
                              1, gd->nbodyTable,
                              gd->iCatalogs[0], gd->iCatalogs[1]);
        break;

#endif	// ! _cballs_kdtree_box_omp_h
