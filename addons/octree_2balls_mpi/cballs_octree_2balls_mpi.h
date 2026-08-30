#ifndef _cballs_octree_2balls_mpi_h
#define _cballs_octree_2balls_mpi_h

case OCTREE2BALLSMPIMETHOD: {
    int tree_status = SUCCESS;

    verb_print(cmd->verbose,
               "\n\tevalHist: with distributed octree two-ball 2PCF "
               "and LogMultipole 3PCF\n\n");
    if (cballs_opt_read_mask(cmd)) {
        ifile = gd->iCatalogs[0];
        DO_BODY(p, bodytable[ifile],
                bodytable[ifile] + gd->nbodyTable[ifile])
            Update(p) = TRUE;
        tree_status = MakeTree(cmd, gd, bodytable[ifile],
                               gd->nbodyTable[ifile], ifile);
    } else {
        for (ifile = 0; ifile < gd->ninfiles; ifile++) {
            DO_BODY(p, bodytable[ifile],
                    bodytable[ifile] + gd->nbodyTable[ifile])
                Update(p) = TRUE;
            if (tree_status == SUCCESS)
                tree_status = MakeTree(cmd, gd, bodytable[ifile],
                                       gd->nbodyTable[ifile], ifile);
        }
    }
    tree_status = fcfc_octree_2balls_mpi_consensus(
        cmd, tree_status, "MPI native-octree construction");
    if (tree_status == FAILURE) return FAILURE;

    if (cballs_opt_read_mask(cmd)) {
        if (searchcalc_octree_2balls_mpi(
                cmd, gd, bodytable, gd->nbodyTable, 1, gd->nbodyTable,
                ifile, ifile) == FAILURE)
            return FAILURE;
    } else if (searchcalc_octree_2balls_mpi(
            cmd, gd, bodytable, gd->nbodyTable, 1, gd->nbodyTable,
            gd->iCatalogs[0], gd->iCatalogs[1]) == FAILURE) {
        return FAILURE;
    }
    break;
}

#endif /* !_cballs_octree_2balls_mpi_h */
