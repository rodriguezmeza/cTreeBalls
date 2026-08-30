// Use:
//#include "cballs_octree_ggg.h"
//
// included in file: addons/addons_include/source/cballs/cballs_include02.h
//

#ifndef _cballs_octree_ggg_h
#define _cballs_octree_ggg_h

//
// look search routines tag numbers in file:
//              addons/addons_include/source/startrun/startrun_include_11.h
//

//#define OCTREEGGG     75

case 75:                   // search=octree-ggg
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                           "\n\tevalHist: with octree-ggg method\n\n",
                           routineName);

if (cballs_opt_read_mask(cmd)) {
    ifile=0;
    DO_BODY(p,bodytable[ifile],bodytable[ifile]+gd->nbodyTable[ifile]) {
        Update(p) = TRUE;
    }
    if (MakeTree(cmd, gd, bodytable[ifile], gd->nbodyTable[ifile], ifile)
        == FAILURE)
        return FAILURE;
} else {
    for (ifile=0; ifile<gd->ninfiles; ifile++) {
        DO_BODY(p,bodytable[ifile],bodytable[ifile]+gd->nbodyTable[ifile]) {
            Update(p) = TRUE;
        }
        if (MakeTree(cmd, gd, bodytable[ifile], gd->nbodyTable[ifile], ifile)
            == FAILURE)
            return FAILURE;
    }
}

    switch(correlation_int) {
        case GGGCORRELATION:
            if (cballs_opt_read_mask(cmd)) {
                searchcalc_octree_ggg(cmd, gd, bodytable, gd->nbodyTable, 1,
                                      gd->nbodyTable,
                                      gd->iCatalogs[0], gd->iCatalogs[0]);
            } else {
                searchcalc_octree_ggg(cmd, gd, bodytable, gd->nbodyTable, 1,
                                      gd->nbodyTable,
                                      gd->iCatalogs[0], gd->iCatalogs[1]);
            }
            break;
        default:
            if (cballs_opt_read_mask(cmd)) {
                searchcalc_octree_ggg(cmd, gd, bodytable, gd->nbodyTable, 1,
                                      gd->nbodyTable,
                                      gd->iCatalogs[0], gd->iCatalogs[0]);
            } else {
                searchcalc_octree_ggg(cmd, gd, bodytable, gd->nbodyTable, 1,
                                      gd->nbodyTable,
                                      gd->iCatalogs[0], gd->iCatalogs[1]);
            }
            break;
    }
    break;

#endif	// ! _cballs_octree_ggg_h
