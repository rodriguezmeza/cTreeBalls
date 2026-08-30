// Use:
//#include "cballs_octree_3pcf_3d_omp.h"

#ifndef _cballs_octree_3pcf_3d_omp_h
#define _cballs_octree_3pcf_3d_omp_h

case 166:                   // search=octree-3pcf-3d-omp
    {
    int survey_mode = scanopt(cmd->options, "survey-estimator-3d")
                   || scanopt(cmd->options, "encore-survey-estimator")
                   || scanopt(cmd->options, "survey-edge-correction");
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                           "\n\tevalHist: with octree-3pcf-3d-omp method\n\n");

#if NDIM != 3
    cBALLS_FAIL(cmd,
                "octree-3pcf-3d-omp requires DEFDIMENSION=3 in Makefile_settings");
#endif

    if (survey_mode) {
        int dmr_cat;
        int random_cat;
        REAL random_scale;

        if (cballs_opt_read_mask(cmd))
            cBALLS_FAIL(cmd,
                        "octree-3pcf-3d-omp survey mode does not use read-mask; "
                        "encode the selection in the data and random catalogs");
        if (gd->ninfiles != 2 || gd->iCatalogs[0] == gd->iCatalogs[1])
            cBALLS_FAIL(cmd,
                        "octree-3pcf-3d-omp survey mode requires exactly two "
                        "distinct catalogs: data,random");

        if (cb3d_prepare_survey_catalogs(
                cmd, gd, bodytable, gd->nbodyTable,
                gd->iCatalogs[0], gd->iCatalogs[1],
                &dmr_cat, &random_cat, &random_scale) == FAILURE)
            return FAILURE;
        if (MakeTree(cmd, gd, bodytable[dmr_cat],
                     gd->nbodyTable[dmr_cat], dmr_cat) == FAILURE)
            return FAILURE;
        if (MakeTree(cmd, gd, bodytable[random_cat],
                     gd->nbodyTable[random_cat], random_cat) == FAILURE)
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
    }

#endif // !_cballs_octree_3pcf_3d_omp_h
