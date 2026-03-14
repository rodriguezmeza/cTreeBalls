// ============================================================================*/
//        1          2          3          4        ^ 5          6          7

// Use:
//#include "cballs_octree_g1g2g3_omp.h"

//
// included in file: addons/addons_include/source/cballs/cballs_include_02.h
//

#ifndef _cballs_octree_g1g2g3_omp_h
#define _cballs_octree_g1g2g3_omp_h

//
// look search routines tag numbers in file:
//              addons/addons_include/source/startrun/startrun_include_11.h
//

//#define OCTREEG1G2G3OMP     74

case 74:                   // search=octree-ggg-omp
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                        "\n\tevalHist: with octree-g1g2g3-omp method\n\n",
                           routineName);

    if (scanopt(cmd->options, "read-mask")) {
        ifile=0;
        DO_BODY(p,bodytable[ifile],bodytable[ifile]
                +gd->nbodyTable[ifile]) {
            Update(p) = TRUE;
        }
        MakeTree(cmd, gd, bodytable[ifile], gd->nbodyTable[ifile], ifile);
    } else {
        for (ifile=0; ifile<gd->ninfiles; ifile++) {
            DO_BODY(p,bodytable[ifile],bodytable[ifile]
                    +gd->nbodyTable[ifile]) {
                Update(p) = TRUE;
            }
            MakeTree(cmd, gd, bodytable[ifile],
                     gd->nbodyTable[ifile], ifile);
        }
    }

    switch(correlation_int) {
        case GGGCORRELATION:
            if (scanopt(cmd->options, "read-mask")) {
                searchcalc_octree_g1g2g3_omp(cmd, gd,
                                             bodytable, gd->nbodyTable, 1,
                                             gd->nbodyTable,
                                             gd->iCatalogs[0],
                                             gd->iCatalogs[1],
                                             gd->iCatalogs[2]);
            } else {
                searchcalc_octree_g1g2g3_omp(cmd, gd,
                                             bodytable, gd->nbodyTable, 1,
                                             gd->nbodyTable,
                                             gd->iCatalogs[0],
                                             gd->iCatalogs[1],
                                             gd->iCatalogs[2]);
            }
            break;
        default:
            if (scanopt(cmd->options, "read-mask")) {
                searchcalc_octree_g1g2g3_omp(cmd, gd,
                                             bodytable, gd->nbodyTable, 1,
                                             gd->nbodyTable,
                                             gd->iCatalogs[0],
                                             gd->iCatalogs[1],
                                             gd->iCatalogs[2]);
            } else {
                searchcalc_octree_g1g2g3_omp(cmd, gd,
                                             bodytable, gd->nbodyTable, 1,
                                             gd->nbodyTable,
                                             gd->iCatalogs[0],
                                             gd->iCatalogs[1],
                                             gd->iCatalogs[2]);
            }
            break;
    }

break;

#endif	// ! _cballs_octree_g1g2g3_omp_h
