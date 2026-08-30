#ifndef _cballs_lya_forest_omp_h
#define _cballs_lya_forest_omp_h

case 169:
case 170:
case 171:
    if (gd->ninfiles != 1) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s requires exactly one flattened Lyman-alpha catalog",
                 cmd->searchMethod);
        return FAILURE;
    }
    ifile = gd->iCatalogs[0];
    DO_BODY(p, bodytable[ifile],
            bodytable[ifile] + gd->nbodyTable[ifile]) {
        Update(p) = TRUE;
    }
    if (MakeTree(cmd, gd, bodytable[ifile], gd->nbodyTable[ifile], ifile)
        == FAILURE)
        return FAILURE;
    if (searchcalc_lya_forest_omp(
            cmd, gd, bodytable, gd->nbodyTable, 1, gd->nbodyTable, ifile,
            gd->searchMethod_int != 170,
            gd->searchMethod_int != 169) == FAILURE)
        return FAILURE;
    break;

case LYAFOREST1D2PCFMETHOD:
case LYAFOREST1D3PCFMETHOD:
case LYAFOREST1D2PCF3PCFMETHOD:
    if (gd->ninfiles != 1) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s requires exactly one flattened Lyman-alpha catalog",
                 cmd->searchMethod);
        return FAILURE;
    }
    ifile = gd->iCatalogs[0];
    DO_BODY(p, bodytable[ifile],
            bodytable[ifile] + gd->nbodyTable[ifile]) {
        Update(p) = TRUE;
    }
    if (searchcalc_lya_forest_1d_omp(
            cmd, gd, bodytable[ifile], gd->nbodyTable[ifile],
            gd->searchMethod_int != LYAFOREST1D3PCFMETHOD,
            gd->searchMethod_int != LYAFOREST1D2PCFMETHOD) == FAILURE)
        return FAILURE;
    break;

case LYAFOREST1DTREE2PCFMETHOD:
    if (gd->ninfiles != 1) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s requires exactly one flattened Lyman-alpha catalog",
                 cmd->searchMethod);
        return FAILURE;
    }
    ifile = gd->iCatalogs[0];
    DO_BODY(p, bodytable[ifile],
            bodytable[ifile] + gd->nbodyTable[ifile]) {
        Update(p) = TRUE;
    }
    if (searchcalc_lya_forest_1d_tree_omp(
            cmd, gd, bodytable[ifile], gd->nbodyTable[ifile]) == FAILURE)
        return FAILURE;
    break;

#endif
