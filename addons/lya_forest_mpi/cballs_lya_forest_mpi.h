case 185:
case 186:
case 187:
case 188:
case 189:
case 190:
case 191: {
    int lya_status = SUCCESS;
    const int kind = lya_forest_method_kind(cmd->searchMethod);
    if (gd->ninfiles != 1) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s requires exactly one flattened Ly-alpha catalog", cmd->searchMethod);
        lya_status = FAILURE;
    }
    if (lya_forest_mpi_consensus(cmd, lya_status, "Ly-alpha MPI catalog") == FAILURE)
        return FAILURE;
    ifile = gd->iCatalogs[0];
    DO_BODY(p, bodytable[ifile], bodytable[ifile] + gd->nbodyTable[ifile])
        Update(p) = TRUE;
    if (kind < 3) {
        lya_status = MakeTree(cmd, gd, bodytable[ifile], gd->nbodyTable[ifile], ifile);
        if (lya_forest_mpi_consensus(cmd, lya_status, "Ly-alpha MPI tree") == FAILURE)
            return FAILURE;
        lya_status = searchcalc_lya_forest_omp(cmd, gd, bodytable, gd->nbodyTable,
                            1, gd->nbodyTable, ifile, kind != 1, kind != 0);
    } else if (kind < 6) {
        lya_status = searchcalc_lya_forest_1d_omp(cmd, gd, bodytable[ifile],
                            gd->nbodyTable[ifile], kind != 4, kind != 3);
    } else {
        lya_status = searchcalc_lya_forest_1d_tree_omp(cmd, gd, bodytable[ifile],
                                                    gd->nbodyTable[ifile]);
    }
    if (lya_status == FAILURE) return FAILURE;
    break;
}
