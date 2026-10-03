#ifndef _cballs_octree_shear_sphere_2balls_mpi_h
#define _cballs_octree_shear_sphere_2balls_mpi_h

case OCTREESHEARSPHERE2BALLSMPIMETHOD:
{
    int tree_status = SUCCESS;
    int shear_cat1;
    int shear_cat2;
    int shear_cat3;

    verb_print_normal_info(
        cmd->verbose, cmd->verbose_log, gd->outlog,
        "\n\tevalHist: with octree-shear-sphere-2balls-mpi method\n\n");
    cballs_mpi_local_begin();
    for (ifile = 0; ifile < gd->ninfiles; ifile++) {
        DO_BODY(p, bodytable[ifile],
                bodytable[ifile] + gd->nbodyTable[ifile])
            Update(p) = TRUE;
        if (tree_status == SUCCESS)
            tree_status = MakeTree(cmd, gd, bodytable[ifile],
                                   gd->nbodyTable[ifile], ifile);
    }
    cballs_mpi_local_end();
    if (fcfc_octree_shear_sphere_2balls_mpi_consensus(
            cmd, tree_status, "MPI shear native-octree construction") == FAILURE)
        return FAILURE;
    shear_cat1 = gd->iCatalogs[0];
    shear_cat2 = gd->ninfiles >= 2 ? gd->iCatalogs[1] : shear_cat1;
    shear_cat3 = gd->ninfiles >= 3 ? gd->iCatalogs[2] : shear_cat2;
    if (searchcalc_octree_shear_sphere_2balls_mpi(
            cmd, gd, bodytable, gd->nbodyTable, 1, gd->nbodyTable,
            shear_cat1, shear_cat2, shear_cat3) == FAILURE)
        return FAILURE;
    break;
}

#endif
