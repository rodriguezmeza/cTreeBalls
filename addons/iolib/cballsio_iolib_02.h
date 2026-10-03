
#ifndef _cballsio_iolib_02_h
#define _cballsio_iolib_02_h

local bool iolib_preserve_common_catalog_frame(
        struct cmdline_data *cmd, struct global_data *gd)
{
    bool preserve = cballs_observer_frame(cmd);

    (void)cmd;
    (void)gd;

#if defined(OCTREE3PCF3DOMP) || defined(OCTREE3PCF3DMPI)
    preserve |= cb3d_is_method(cmd->searchMethod)
        && (scanopt(cmd->options, "survey-estimator-3d")
            || scanopt(cmd->options, "encore-survey-estimator")
            || scanopt(cmd->options, "survey-edge-correction"));
#endif
#ifdef BALLTREE2BALLSOMP
    preserve |= gd->searchMethod_int == BALLTREE2BALLSMETHOD
        && gd->iCatalogs[0] != gd->iCatalogs[1];
#endif
#ifdef BALLTREE2BALLSMPI
    preserve |= gd->searchMethod_int == BALLTREE2BALLSMPIMETHOD
        && gd->iCatalogs[0] != gd->iCatalogs[1];
#endif
#ifdef KDTREE2BALLSOMP
    preserve |= gd->searchMethod_int == KDTREE2BALLSOMPMETHOD
        && gd->iCatalogs[0] != gd->iCatalogs[1];
#endif
#ifdef KDTREE2BALLSMPI
    preserve |= gd->searchMethod_int == KDTREE2BALLSMPIMETHOD
        && gd->iCatalogs[0] != gd->iCatalogs[1];
#endif
#ifdef OCTREE2BALLSOMP
    preserve |= gd->searchMethod_int == OCTREE2BALLSMETHOD
        && gd->iCatalogs[0] != gd->iCatalogs[1];
#endif
#ifdef BALLTREE2BALLSOMP3PCF
    preserve |= gd->searchMethod_int == BALLTREE2BALLSOMP3PCFMETHOD
        && gd->iCatalogs[0] != gd->iCatalogs[1];
#endif
#ifdef BALLTREE2BALLSMPI3PCF
    preserve |= gd->searchMethod_int == BALLTREE2BALLSMPI3PCFMETHOD
        && gd->iCatalogs[0] != gd->iCatalogs[1];
#endif
    return preserve;
}

// Layouts 2 and 3 share the checked native ASCII parser and cleanup.
local int inputdata_ascii_pos(struct cmdline_data *cmd, struct global_data *gd,
                               string filename, int ifile)
{
    return inputdata_ascii_columns(cmd, gd, filename, ifile, 2);
}

#if NDIM == 3
local int inputdata_ascii_2d_to_3d(struct cmdline_data *cmd, struct global_data *gd,
                                   string filename, int ifile)
{
    return inputdata_ascii_columns(cmd, gd, filename, ifile, 3);
}
#endif

/* Selected-column buffers are private to the active runtime. */
local void iolib_free_columns(void)
{
    for (int j = 0; j < 6; ++j) {
        free(cballs_runtime_current()->io_columns[j]);
        cballs_runtime_current()->io_columns[j] = NULL;
    }
}

// infileformat: multi-columns-ascii
local int inputdata_ascii_mcolumns(struct cmdline_data* cmd, struct  global_data* gd,
                               string filename, int ifile)
{
    string routineName = "inputdata_ascii_mcolumns";
    stream instr;
    int ndim;
    bodyptr p;
    real mass=1;
    real weight=1;
    int convergence_weight =
        scanopt(cmd->options, "pos-and-convergence-weight");
    const bool preserve_common_survey_frame =
        iolib_preserve_common_catalog_frame(cmd, gd);

    gd->input_comment = "Multi-column position input file";

    verb_print_warning(cmd->verbose, "\nWarning!! %s\n\t%s... \t%s\n",
                       "need to use \"columns\" option to give positions columns correctly",
                       "not doing so, may result in wrong results ",
                       "(and even get segmentation fault)");
//
// Upgrade to cross-correlations: kappa-gamma
//
    if ( (scanopt(cmd->options, "pos-and-convergence"))
            && (scanopt(cmd->options, "pos-and-shear")))
        cBALLS_FAIL(cmd, "inputdata_ascii_mcolumns: mutually excluyent options: %s",
                    "'pos-and-convergence' and 'pos-and-shear'");
    if (convergence_weight && scanopt(cmd->options, "pos-and-shear"))
        cBALLS_FAIL(cmd, "inputdata_ascii_mcolumns: mutually exclusive options: %s",
                    "'pos-and-convergence-weight' and 'pos-and-shear'");
    if (convergence_weight && scanopt(cmd->options, "pos-and-convergence"))
        cBALLS_FAIL(cmd, "inputdata_ascii_mcolumns: mutually exclusive options: %s",
                    "'pos-and-convergence-weight' and 'pos-and-convergence'");
#if NDIM != 3
    if (convergence_weight)
        cBALLS_FAIL(cmd,
                    "inputdata_ascii_mcolumns: pos-and-convergence-weight requires NDIM=3");
#endif

    int vnpoint;

#if NDIM == 3
    if (scanopt(cmd->options, "only-pos")) {
        if (InputData_3c(cmd, filename,
                         gd->columns[0], gd->columns[1], gd->columns[2], &vnpoint) == FAILURE)
            return FAILURE;
    } else {
        if (convergence_weight) {
            if (InputData_5c(cmd, filename,
                             gd->columns[0], gd->columns[1], gd->columns[2],
                             gd->columns[3], gd->columns[4],
                             &vnpoint) == FAILURE)
                return FAILURE;
        } else if (scanopt(cmd->options, "pos-and-convergence")) {
            if (InputData_4c(cmd, filename,
                             gd->columns[0], gd->columns[1], gd->columns[2],
                             gd->columns[3],
                             &vnpoint) == FAILURE)
                return FAILURE;
        } else {
            if (scanopt(cmd->options, "pos-and-shear")) {
                if (InputData_5c(cmd, filename,
                                 gd->columns[0], gd->columns[1], gd->columns[2],
                                 gd->columns[3], gd->columns[4],
                                 &vnpoint) == FAILURE)
                    return FAILURE;
            } else
                cBALLS_FAIL(cmd, "\n\t%s: options need one of: \n\t%s, %s, %s or %s\n",
                            routineName,
                            "only-pos", "pos-and-convergence",
                            "pos-and-convergence-weight", "pos-and-shear");
        }
    }
#else
    if (scanopt(cmd->options, "only-pos")) {
        if (InputData_2c(cmd, filename, gd->columns[0], gd->columns[1], &vnpoint) == FAILURE)
            return FAILURE;
    } else {
        if (scanopt(cmd->options, "pos-and-convergence")) {
            if (InputData_3c(cmd, filename,
                             gd->columns[0], gd->columns[1],
                             gd->columns[2],
                             &vnpoint) == FAILURE)
                return FAILURE;
        } else {
            if (scanopt(cmd->options, "pos-and-shear")) {
                if (InputData_4c(cmd, filename,
                                 gd->columns[0], gd->columns[1],
                                 gd->columns[2], gd->columns[3],
                                 &vnpoint) == FAILURE)
                    return FAILURE;
            }
        }
    }
#endif

    ndim = NDIM;
    cmd->nbody = vnpoint;
    gd->nbodyTable[ifile] = cmd->nbody;

    if (cballs_calloc_checked((void **)&bodytable[ifile], (size_t)cmd->nbody,
                               sizeof(body), "ASCII column catalog",
                               cmd->error_message, _ERRORMSGSIZE_) == FAILURE) {
        iolib_free_columns();
        return FAILURE;
    }
    gd->bodytable_allocated = TRUE;

    DO_BODY(p, bodytable[ifile], bodytable[ifile]+cmd->nbody) {
        Type(p) = BODY;
        Mass(p) = mass;
        Weight(p) = weight;
        Kappa(p) = 2.0;
        Mask(p) = MASK_NODE_VALID;
        Id(p) = p-bodytable[ifile]+1;
#if defined(OCTREE3PCF3DOMP) || defined(OCTREE3PCF3DMPI)
        Octree3pcf3dLosId(p) = Id(p);
#endif
#if NDIM == 3
        Pos(p)[0] = inout_xval[Id(p)-1];
        Pos(p)[1] = inout_yval[Id(p)-1];
        Pos(p)[2] = inout_zval[Id(p)-1];
        if (scanopt(cmd->options, "pos-and-convergence"))
            Kappa(p) = inout_uval[Id(p)-1];
        if (convergence_weight) {
            Kappa(p) = inout_uval[Id(p)-1];
            Weight(p) = inout_vval[Id(p)-1];
        }

#ifdef THREEPCFSHEAR
        if (scanopt(cmd->options, "pos-and-shear")) {
            Gamma1(p) = inout_uval[Id(p)-1];
            Gamma2(p) = inout_vval[Id(p)-1];
        }
#endif

#else // ! NDIM
        Pos(p)[0] = inout_xval[Id(p)-1];
        Pos(p)[1] = inout_yval[Id(p)-1];
        if (scanopt(cmd->options, "pos-and-convergence"))
            Kappa(p) = inout_zval[Id(p)-1];

#ifdef THREEPCFSHEAR
        if (scanopt(cmd->options, "pos-and-shear")) {
            Gamma1(p) = inout_zval[Id(p)-1];
            Gamma2(p) = inout_uval[Id(p)-1];
        }
#endif

#endif // ! NDIM
    }

    iolib_free_columns();

    if (cballs_input_validate_bodies(cmd, filename, bodytable[ifile],
                                     gd->nbodyTable[ifile]) == FAILURE)
        return FAILURE;

//B
//B Set (0,0,...) as the center of the box
// By now it is only working with boxes centered at (0,0,...)
    cellptr root;                                   // Set it up a temporal root
    int root_status = FAILURE;
    root = (cellptr) allocate(1 * sizeof(body));

    if (FindRootCenter(cmd, gd, bodytable[ifile], gd->nbodyTable[ifile],
                       ifile, root) == FAILURE)
        goto cleanup_root;
    if (!preserve_common_survey_frame
        && centerBodies(bodytable[ifile], gd->nbodyTable[ifile], ifile, root)
           == FAILURE)
        goto cleanup_root;
    if (FindRootCenter(cmd, gd, bodytable[ifile], gd->nbodyTable[ifile],
                       ifile, root) == FAILURE)
        goto cleanup_root;
    if (!preserve_common_survey_frame)
        CLRV(Pos(root));
//E
    gd->rSizeTable[ifile] = 1.0;
    if (expandbox(cmd, gd, bodytable[ifile], gd->nbodyTable[ifile],
                  ifile, root) == FAILURE)
        goto cleanup_root;
    root_status = SUCCESS;

cleanup_root:
    free(root);
    if (root_status == FAILURE)
        return FAILURE;
#if NDIM == 3
    real xmin, ymin, zmin;
    real xmax, ymax, zmax;
    xmin = Pos(bodytable[ifile])[0];
    ymin = Pos(bodytable[ifile])[1];
    zmin = Pos(bodytable[ifile])[2];
    xmax = Pos(bodytable[ifile])[0];
    ymax = Pos(bodytable[ifile])[1];
    zmax = Pos(bodytable[ifile])[2];
    DO_BODY(p, bodytable[ifile], bodytable[ifile]+cmd->nbody) {
        xmin = MIN(xmin,Pos(p)[0]);
        ymin = MIN(ymin,Pos(p)[1]);
        zmin = MIN(zmin,Pos(p)[2]);
        xmax = MAX(xmax,Pos(p)[0]);
        ymax = MAX(ymax,Pos(p)[1]);
        zmax = MAX(zmax,Pos(p)[2]);
    }
    gd->Box[0] = xmax-xmin;
    gd->Box[1] = ymax-ymin; gd->Box[2] = zmax-zmin;
#else
    real xmin, ymin;
    real xmax, ymax;
    xmin = Pos(bodytable[ifile])[0];
    ymin = Pos(bodytable[ifile])[1];
    xmax = Pos(bodytable[ifile])[0];
    ymax = Pos(bodytable[ifile])[1];
    DO_BODY(p, bodytable[ifile], bodytable[ifile]+cmd->nbody) {
        xmin = MIN(xmin,Pos(p)[0]);
        ymin = MIN(ymin,Pos(p)[1]);
        xmax = MAX(xmax,Pos(p)[0]);
        ymax = MAX(ymax,Pos(p)[1]);
    }
    gd->Box[0] = xmax-xmin;
    gd->Box[1] = ymax-ymin;
#endif
    verb_print_q(2, cmd->verbose, 
                 "\tinputdata_ascii_mcolumns: nbody and ndim: %" INTEGER_FMT " %d...\n",
                 cmd->nbody, ndim);
    verb_print_q(2, cmd->verbose,
                 "\tinputdata_ascii_mcolumns: Lx, Ly, Lz: %g %g %g...\n",
               gd->Box[0], gd->Box[1], gd->Box[2]);

    real kavg=0.0;
    DO_BODY(p, bodytable[ifile], bodytable[ifile]+cmd->nbody) {
        Type(p) = BODY;
        Mass(p) = mass;
        if (!convergence_weight)
            Weight(p) = weight;
        Mask(p) = MASK_NODE_VALID;
        if (scanopt(cmd->options, "kappa-constant")) {
            if (scanopt(cmd->options, "kappa-constant-one"))
                Kappa(p) = 1.0;
            else
                Kappa(p) = 2.0;
        }
        kavg += Kappa(p);
    }
#if defined(OCTREE3PCF3DOMP) || defined(OCTREE3PCF3DMPI)
    gd->octree3pcf3d_los_ids[ifile] = TRUE;
#endif
    verb_print_q(2, cmd->verbose,
               "inputdata_ascii_mcolumns: average of kappa (%" INTEGER_FMT
               " particles) = %le\n",
               cmd->nbody, kavg/((real)cmd->nbody) );

    //B If needed locate particles with same position.
    //  This is a slow process, use if it necessary...
    if (scanopt(cmd->options, "check-eq-pos")) {
        bodyptr q;
        real dist2;
        vector distv;
        INTEGER pqequals=0;
        bool flag=0;
        DO_BODY(p, bodytable[ifile], bodytable[ifile]+cmd->nbody-1) {
            DO_BODY(q, p+1, bodytable[ifile]+cmd->nbody) {
                if (p != q) {
                    DOTPSUBV(dist2, distv, Pos(p), Pos(q));
                    if (dist2 == 0.0) {
                        flag=1;
                        pqequals++;
                        verb_print(cmd->verbose,
                        "inputdata_ascii_mcolumns: (p,q) = %ld %ld :: %g %g %g :: %g %g %g\n",
                        Id(p), Id(q),
                        Pos(p)[0], Pos(p)[1], Pos(p)[2], Pos(q)[0], Pos(q)[1], Pos(q)[2]);
                    }
                }
            }
        }
        verb_print(cmd->verbose,
        "inputdata_ascii_mcolumns: Total equal pairs: %ld\n",pqequals);
        if (flag)
            cBALLS_FAIL(cmd,
                        "inputdata_ascii_mcolumns: at least two bodies have same position\n");
    }
    //E

    return SUCCESS;
}

// Routine to read ra-dec ascii files
//  by default:
//      ra is column 1
//      dec is column 2
//      kappa is column 3
local int inputdata_ascii_ra_dec(struct cmdline_data* cmd, struct  global_data* gd,
                               string filename, int ifile)
{
    string routineName = "inputdata_ascii_ra_dec";
    stream instr;
    int ndim = NDIM;
    bodyptr p;
    real mass=1;
    real weight=1;
    real phi, theta;
    const bool preserve_common_survey_frame =
        iolib_preserve_common_catalog_frame(cmd, gd);

    gd->input_comment = "ra-dec-kappa input file";

    int vnpoint;

#if NDIM == 2
    cBALLS_FAIL(cmd, "\n%s only works in 3D", routineName);
#endif

    if (InputData_3c(cmd, filename,
                     gd->columns[0], gd->columns[1], gd->columns[2], &vnpoint) == FAILURE)
        return FAILURE;

    cmd->nbody = vnpoint;
    gd->nbodyTable[ifile] = cmd->nbody;

    if (cballs_calloc_checked((void **)&bodytable[ifile], (size_t)cmd->nbody,
                               sizeof(body), "ASCII column catalog",
                               cmd->error_message, _ERRORMSGSIZE_) == FAILURE) {
        iolib_free_columns();
        return FAILURE;
    }
    gd->bodytable_allocated = TRUE;

    DO_BODY(p, bodytable[ifile], bodytable[ifile]+cmd->nbody) {
        Id(p) = p-bodytable[ifile]+1;
        Weight(p) = 1.0;
        Mask(p) = MASK_NODE_VALID;
        phi = inout_xval[Id(p)-1];
        theta = inout_yval[Id(p)-1];
        if (scanopt(cmd->options, "in-degrees")) {
            phi = radian(phi);
            theta = radian(theta);
        }
        coordinate_transformation(cmd, gd, theta, phi, Pos(p));
        Kappa(p) = inout_zval[Id(p)-1];
    }

    iolib_free_columns();
    if (cballs_input_validate_bodies(cmd, filename, bodytable[ifile],
                                     gd->nbodyTable[ifile]) == FAILURE)
        return FAILURE;

//B Set (0,0,...) as the center of the box
// By now it is only working with boxes centered at (0,0,...)
    cellptr root;                                   // Set it up a temporal root
    int root_status = FAILURE;
    root = (cellptr) allocate(1 * sizeof(body));

    if (FindRootCenter(cmd, gd, bodytable[ifile], gd->nbodyTable[ifile],
                       ifile, root) == FAILURE)
        goto cleanup_root;
    if (!preserve_common_survey_frame
        && centerBodies(bodytable[ifile], gd->nbodyTable[ifile], ifile, root)
           == FAILURE)
        goto cleanup_root;
    if (FindRootCenter(cmd, gd, bodytable[ifile], gd->nbodyTable[ifile],
                       ifile, root) == FAILURE)
        goto cleanup_root;
    if (!preserve_common_survey_frame)
        CLRV(Pos(root));
//E
    gd->rSizeTable[ifile] = 1.0;
    if (expandbox(cmd, gd, bodytable[ifile], gd->nbodyTable[ifile],
                  ifile, root) == FAILURE)
        goto cleanup_root;
    root_status = SUCCESS;

cleanup_root:
    free(root);
    if (root_status == FAILURE)
        return FAILURE;

    real xmin, ymin, zmin;
    real xmax, ymax, zmax;
    xmin = Pos(bodytable[ifile])[0];
    ymin = Pos(bodytable[ifile])[1];
    zmin = Pos(bodytable[ifile])[2];
    xmax = Pos(bodytable[ifile])[0];
    ymax = Pos(bodytable[ifile])[1];
    zmax = Pos(bodytable[ifile])[2];
    DO_BODY(p, bodytable[ifile], bodytable[ifile]+cmd->nbody) {
        xmin = MIN(xmin,Pos(p)[0]);
        ymin = MIN(ymin,Pos(p)[1]);
        zmin = MIN(zmin,Pos(p)[2]);
        xmax = MAX(xmax,Pos(p)[0]);
        ymax = MAX(ymax,Pos(p)[1]);
        zmax = MAX(zmax,Pos(p)[2]);
    }
    gd->Box[0] = xmax-xmin;
    gd->Box[1] = ymax-ymin; gd->Box[2] = zmax-zmin;
    verb_print_q(2, cmd->verbose,
                 "\tinputdata_ascii_mcolumns: nbody and ndim: %" INTEGER_FMT " %d...\n",
                 cmd->nbody, ndim);
    verb_print_q(2, cmd->verbose,
                 "\tinputdata_ascii_mcolumns: Lx, Ly, Lz: %g %g %g...\n",
               gd->Box[0], gd->Box[1], gd->Box[2]);

    real kavg=0.0;
    DO_BODY(p, bodytable[ifile], bodytable[ifile]+cmd->nbody) {
        Type(p) = BODY;
        Mass(p) = mass;
        Weight(p) = weight;
        if (scanopt(cmd->options, "kappa-constant")) {
            if (scanopt(cmd->options, "kappa-constant-one"))
                Kappa(p) = 1.0;
            else
                Kappa(p) = 2.0;
        }
        kavg += Kappa(p);
    }
    verb_print_q(2, cmd->verbose,
               "inputdata_ascii_mcolumns: average of kappa (%ld particles) = %le\n",
               cmd->nbody, kavg/((real)cmd->nbody) );

    //B If needed locate particles with same position.
    //  This is a slow process, use if it necessary...
    if (scanopt(cmd->options, "check-eq-pos")) {
        bodyptr q;
        real dist2;
        vector distv;
        INTEGER pqequals=0;
        bool flag=0;
        DO_BODY(p, bodytable[ifile], bodytable[ifile]+cmd->nbody-1) {
            DO_BODY(q, p+1, bodytable[ifile]+cmd->nbody) {
                if (p != q) {
                    DOTPSUBV(dist2, distv, Pos(p), Pos(q));
                    if (dist2 == 0.0) {
                        flag=1;
                        pqequals++;
                        verb_print(cmd->verbose,
                                   "%s: (p,q) = %ld %ld :: %g %g %g :: %g %g %g\n",
                                   routineName, Id(p), Id(q),
                                   Pos(p)[0], Pos(p)[1], Pos(p)[2],
                                   Pos(q)[0], Pos(q)[1], Pos(q)[2]);
                    }
                }
            }
        }
        verb_print(cmd->verbose,
        "%s: Total equal pairs: %ld\n",
                   routineName, pqequals);
        if (flag)
            cBALLS_FAIL(cmd, "inputdata_ascii_mcolumns: at least two bodies have same position\n");
    }
    //E

    return SUCCESS;
}

#endif	// ! _cballsio_iolib_02_h
