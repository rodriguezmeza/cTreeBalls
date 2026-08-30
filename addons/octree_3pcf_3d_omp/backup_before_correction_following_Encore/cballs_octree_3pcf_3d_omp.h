// Use:
//#include "cballs_octree_3pcf_3d_omp.h"

#ifndef _cballs_octree_3pcf_3d_omp_h
#define _cballs_octree_3pcf_3d_omp_h

case 166:                   // search=octree-3pcf-3d-omp
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                           "\n\tevalHist: with octree-3pcf-3d-omp method\n\n");

#if NDIM != 3
    cBALLS_FAIL(cmd,
                "octree-3pcf-3d-omp requires DEFDIMENSION=3 in Makefile_settings");
#endif

    if (cballs_opt_read_mask(cmd)) {
        ifile = 0;
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

    if (cballs_opt_read_mask(cmd)) {
        if (searchcalc_octree_3pcf_3d_omp(
                cmd, gd, bodytable, gd->nbodyTable, 1, gd->nbodyTable,
                gd->iCatalogs[0], gd->iCatalogs[0]) == FAILURE)
            return FAILURE;
    } else {
        if (searchcalc_octree_3pcf_3d_omp(
                cmd, gd, bodytable, gd->nbodyTable, 1, gd->nbodyTable,
                gd->iCatalogs[0], gd->iCatalogs[1]) == FAILURE)
            return FAILURE;
    }
    break;

#endif // !_cballs_octree_3pcf_3d_omp_h
