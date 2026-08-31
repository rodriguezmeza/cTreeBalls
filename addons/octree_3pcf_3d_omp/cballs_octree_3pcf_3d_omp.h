// Use:
//#include "cballs_octree_3pcf_3d_omp.h"

#ifndef _cballs_octree_3pcf_3d_omp_h
#define _cballs_octree_3pcf_3d_omp_h

#ifdef OCTREE3PCF3DOMP
case 166:
#endif
#ifdef OCTREE3PCF3DMPI
case 192:
#endif
    {
    int survey_mode = scanopt(cmd->options, "survey-estimator-3d")
                   || scanopt(cmd->options, "encore-survey-estimator")
                   || scanopt(cmd->options, "survey-edge-correction");
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                           "\n\tevalHist: with 3D scalar octree method\n\n");

#if NDIM != 3
    cBALLS_FAIL(cmd,
                "3D scalar octree requires DEFDIMENSION=3 in Makefile_settings");
#endif

    if (survey_mode) {
        int dmr_cat;
        int random_cat;
        REAL random_scale;

        if (cballs_opt_read_mask(cmd))
            cBALLS_FAIL(cmd,
                        "3D scalar octree survey mode does not use read-mask; "
                        "encode the selection in the data and random catalogs");
        if (gd->ninfiles != 2 || gd->iCatalogs[0] == gd->iCatalogs[1])
            cBALLS_FAIL(cmd,
                        "3D scalar octree survey mode requires exactly two "
                        "distinct catalogs: data,random");

        if (cb3d_parallel_consensus(cmd, cb3d_prepare_survey_catalogs(
                cmd, gd, bodytable, gd->nbodyTable,
                gd->iCatalogs[0], gd->iCatalogs[1],
                &dmr_cat, &random_cat, &random_scale),
                "3D catalog/tree preparation") == FAILURE)
            return FAILURE;
        if (cb3d_parallel_consensus(cmd, MakeTree(cmd, gd, bodytable[dmr_cat],
                     gd->nbodyTable[dmr_cat], dmr_cat),
                "3D catalog/tree preparation") == FAILURE)
            return FAILURE;
        if (cb3d_parallel_consensus(cmd, MakeTree(cmd, gd, bodytable[random_cat],
                     gd->nbodyTable[random_cat], random_cat),
                "3D catalog/tree preparation") == FAILURE)
            return FAILURE;
        if (searchcalc_octree_3pcf_3d_omp_survey(
                cmd, gd, bodytable, gd->nbodyTable,
                dmr_cat, random_cat, random_scale) == FAILURE)
            return FAILURE;
        break;
    }

    if (cballs_opt_read_mask(cmd)) {
        ifile = 0;
        DO_BODY(p,bodytable[ifile],bodytable[ifile]+gd->nbodyTable[ifile]) {
            Update(p) = TRUE;
        }
        if (cb3d_parallel_consensus(cmd, MakeTree(cmd, gd, bodytable[ifile], gd->nbodyTable[ifile], ifile),
                "3D catalog/tree preparation")
            == FAILURE)
            return FAILURE;
    } else {
        for (ifile=0; ifile<gd->ninfiles; ifile++) {
            DO_BODY(p,bodytable[ifile],bodytable[ifile]+gd->nbodyTable[ifile]) {
                Update(p) = TRUE;
            }
            if (cb3d_parallel_consensus(cmd, MakeTree(cmd, gd, bodytable[ifile], gd->nbodyTable[ifile], ifile),
                "3D catalog/tree preparation")
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
    }

#endif // !_cballs_octree_3pcf_3d_omp_h
