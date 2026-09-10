#ifndef CTREEBALLS_CBALLS_KDTREE_SHEAR_SPHERE_2BALLS_OMP_H
#define CTREEBALLS_CBALLS_KDTREE_SHEAR_SPHERE_2BALLS_OMP_H

case KDTREESHEARSPHERE2BALLSOMPMETHOD:
{
    int shear_cat1;
    int shear_cat2;
    int shear_cat3;

    verb_print_normal_info(
        cmd->verbose, cmd->verbose_log, gd->outlog,
        "\n\tevalHist: with kdtree-shear-sphere-2balls-omp method\n\n");
    for (ifile = 0; ifile < gd->ninfiles; ifile++)
        DO_BODY(p, bodytable[ifile],
                bodytable[ifile] + gd->nbodyTable[ifile])
            Update(p) = TRUE;
    shear_cat1 = gd->iCatalogs[0];
    shear_cat2 = gd->ninfiles >= 2 ? gd->iCatalogs[1] : shear_cat1;
    shear_cat3 = gd->ninfiles >= 3 ? gd->iCatalogs[2] : shear_cat2;
    if (searchcalc_kdtree_shear_sphere_2balls_omp(
            cmd, gd, bodytable, gd->nbodyTable, 1, gd->nbodyTable,
            shear_cat1, shear_cat2, shear_cat3) == FAILURE)
        return FAILURE;
    break;
}

#endif
