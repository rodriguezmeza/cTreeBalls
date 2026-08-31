case OCTREEBALLS4MPIMETHOD: {
    int tree_status = SUCCESS;
    for (ifile = 0; ifile < gd->ninfiles; ifile++) {
        DO_BODY(p, bodytable[ifile], bodytable[ifile] + gd->nbodyTable[ifile])
            Update(p) = TRUE;
        if (tree_status == SUCCESS)
            tree_status = MakeTree(cmd, gd, bodytable[ifile],
                                   gd->nbodyTable[ifile], ifile);
    }
    if (fcfc_octree_balls4_mpi_consensus(cmd, tree_status,
            "BALLS4 MPI tree construction") == FAILURE) return FAILURE;
    if (searchcalc_octree_balls4_omp(cmd, gd, bodytable, gd->nbodyTable,
            1, gd->nbodyTable, gd->iCatalogs[0], gd->iCatalogs[1]) == FAILURE)
        return FAILURE;
    break;
}
