#ifndef _cballs_octree_ggg_mpi_h
#define _cballs_octree_ggg_mpi_h

case 172: {                 /* search=octree-ggg-mpi */
    int tree_status = SUCCESS;
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                           "\n\tevalHist: with octree-ggg-mpi method\n\n",
                           routineName);
    if (cballs_opt_read_mask(cmd)) {
        ifile = 0;
        DO_BODY(p, bodytable[ifile],
                bodytable[ifile] + gd->nbodyTable[ifile]) Update(p) = TRUE;
        tree_status = MakeTree(cmd, gd, bodytable[ifile],
                               gd->nbodyTable[ifile], ifile);
    } else {
        for (ifile = 0; ifile < gd->ninfiles; ifile++) {
            DO_BODY(p, bodytable[ifile],
                    bodytable[ifile] + gd->nbodyTable[ifile]) Update(p) = TRUE;
            if (tree_status == SUCCESS)
                tree_status = MakeTree(cmd, gd, bodytable[ifile],
                                       gd->nbodyTable[ifile], ifile);
        }
    }
    tree_status = fcfc_octree_ggg_mpi_consensus(
        cmd, tree_status, "MPI octree construction");
    if (tree_status == FAILURE) return FAILURE;
    if (cballs_opt_read_mask(cmd)) {
        if (searchcalc_octree_ggg_mpi(
                cmd, gd, bodytable, gd->nbodyTable, 1, gd->nbodyTable,
                gd->iCatalogs[0], gd->iCatalogs[0]) == FAILURE)
            return FAILURE;
    } else {
        if (searchcalc_octree_ggg_mpi(
                cmd, gd, bodytable, gd->nbodyTable, 1, gd->nbodyTable,
                gd->iCatalogs[0], gd->iCatalogs[1]) == FAILURE)
            return FAILURE;
    }
    break;
}

#endif
