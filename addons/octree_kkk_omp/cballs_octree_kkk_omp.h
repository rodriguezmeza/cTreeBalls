// Use:
//#include "cballs_octree_kkk_omp.h"

#ifndef _cballs_octree_kkk_omp_h
#define _cballs_octree_kkk_omp_h

//#define OCTREEKKKOMP     61

case 61:                   // search=octree-kkk-omp
    verb_print(cmd->verbose,
    "\n\tevalHist: with octree-kkk-omp method\n\n");
    for (ifile=0; ifile<gd->ninfiles; ifile++) {
        DO_BODY(p,bodytable[ifile],bodytable[ifile]+gd->nbodyTable[ifile]) {
            Update(p) = TRUE;
        }
        MakeTree(cmd, gd, bodytable[ifile], gd->nbodyTable[ifile], ifile);
    }
    switch(correlation_int) {
        case KKKCORRELATION:
            searchcalc_octree_kkk_omp(cmd, gd, bodytable, gd->nbodyTable, 1,
                                      gd->nbodyTable,
                                      gd->iCatalogs[0], gd->iCatalogs[1]);
            break;
        case KKCORRELATION:
            searchcalc_octree_kk_omp(cmd, gd, bodytable, gd->nbodyTable, 1,
                                     gd->nbodyTable,
                                     gd->iCatalogs[0], gd->iCatalogs[1]);
            break;
        case NNCORRELATION:
            searchcalc_octree_nn_omp(cmd, gd, bodytable, gd->nbodyTable, 1,
                                     gd->nbodyTable,
                                     gd->iCatalogs[0], gd->iCatalogs[1]);
            break;
        case NNEstimator:
            nncorrelation_octree_nn_omp(cmd, gd, bodytable, gd->nbodyTable, 1,
                                     gd->nbodyTable,
                                     gd->iCatalogs[0], gd->iCatalogs[1]);
            break;
        default:
            searchcalc_octree_kkk_omp(cmd, gd, bodytable, gd->nbodyTable, 1,
                                     gd->nbodyTable,
                                     gd->iCatalogs[0], gd->iCatalogs[1]);
            break;
    }
    break;

#endif	// ! _cballs_octree_kkk_omp_h
