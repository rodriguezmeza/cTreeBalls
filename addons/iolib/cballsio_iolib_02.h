
#ifndef _cballsio_iolib_02_h
#define _cballsio_iolib_02_h

// infileformat: columns-ascii-pos
local int inputdata_ascii_pos(struct cmdline_data* cmd, struct  global_data* gd,
                               string filename, int ifile)
{
    stream instr;
    int ndim;
    bodyptr p;
    char gato[1], firstline[20];
    real mass=1;
    real weight=1;

    gd->input_comment = "Column position input file";

    instr = stropen(filename, "r");

    fgets(firstline,200,instr);
    fscanf(instr,"%1s",gato);
    in_int_long(instr, &cmd->nbody);
    if (cmd->nbody < 1)
        error("inputdata: nbody = %d is absurd\n", cmd->nbody);
    in_int(instr, &ndim);
    if (ndim != NDIM)
        error("inputdata: ndim = %d; expected %d\n", ndim, NDIM);

    gd->nbodyTable[ifile] = cmd->nbody;

// Check the center of the box!!!
#if NDIM == 3
    real Lx, Ly, Lz;
#ifdef SINGLEP
    in_real_double(instr, &Lx);
    in_real_double(instr, &Ly);
    in_real_double(instr, &Lz);
#else
    in_real(instr, &Lx);
    in_real(instr, &Ly);
    in_real(instr, &Lz);
#endif
    gd->Box[0] = Lx;
    gd->Box[1] = Ly;
    gd->Box[2] = Lz;
#else
    real Lx, Ly;
    in_real(instr, &Lx);
    in_real(instr, &Ly);
    gd->Box[0] = Lx;
    gd->Box[1] = Ly;
#endif

    verb_print(cmd->verbose, "\tInput: nbody and ndim: %d %d...\n", cmd->nbody, ndim);
    verb_print(cmd->verbose, "\tInput: Lx, Ly, Lz: %g %g %g...\n",
               gd->Box[0], gd->Box[1], gd->Box[2]);

    bodytable[ifile] = (bodyptr) allocate(cmd->nbody * sizeof(body));

    DO_BODY(p, bodytable[ifile], bodytable[ifile]+cmd->nbody) {
        in_vector(instr, Pos(p));
    }

    fclose(instr);

    real kavg=0.0;
    DO_BODY(p, bodytable[ifile], bodytable[ifile]+cmd->nbody) {
        Type(p) = BODY;
        Mass(p) = mass;
        Weight(p) = weight;
        Kappa(p) = 2.0;
        Id(p) = p-bodytable[ifile]+1;
        kavg += Kappa(p);
    }
    verb_print(cmd->verbose,
               "inputdata_ascii_pos: average of kappa (%ld particles) = %le\n",
               cmd->nbody, kavg/((real)cmd->nbody) );

    //B If needed locate particles with same position.
    //  This is a slow process, use if it necessary...
    if (scanopt(cmd->options, "check-eq-pos")) {
    bodyptr q;
    real dist2;
    vector distv;
    bool flag=0;
    DO_BODY(p, bodytable[ifile], bodytable[ifile]+cmd->nbody-1)
        DO_BODY(q, p+1, bodytable[ifile]+cmd->nbody)
            if (p != q) {
            DOTPSUBV(dist2, distv, Pos(p), Pos(q));
                if (dist2 == 0.0) {
                    flag=1;
                }
            }
    if (flag) error("inputdata_ascii: at least two bodies have same position\n");
    }
    //E

    return SUCCESS;
}

#if NDIM == 3
// infileformat: columns-ascii-2d-to-3d
local int inputdata_ascii_2d_to_3d(struct cmdline_data* cmd, struct  global_data* gd,
                                    string filename, int ifile)
{
    stream instr;
    int ndim;
    bodyptr p;
    char gato[1], firstline[20];
    real mass=1;
    real weight=1;

    gd->input_comment = "Column form input file (2d-to-3d)";

    instr = stropen(filename, "r");
    fgets(firstline,200,instr);
    fscanf(instr,"%1s",gato);
    in_int_long(instr, &cmd->nbody);
    if (cmd->nbody < 1)
        error("inputdata: nbody = %d is absurd\n", cmd->nbody);
    in_int(instr, &ndim);
    if (ndim != 2)
        error("inputdata: ndim = %d; expected 2\n", ndim);

    gd->nbodyTable[ifile] = cmd->nbody;

// Check the center of the box!!!
    real Lx, Ly;
#ifdef SINGLEP
    in_real_double(instr, &Lx);
    in_real_double(instr, &Ly);
#else
    in_real(instr, &Lx);
    in_real(instr, &Ly);
#endif
    gd->Box[0] = Lx;
    gd->Box[1] = Ly;
// Added this line to set lbox in z direction. Check if it es necessary
    gd->Box[2] = Ly;

    verb_print(cmd->verbose, "\tInput: nbody and ndim: %d %d...\n",
               cmd->nbody, ndim);
    bodytable[ifile] = (bodyptr) allocate(cmd->nbody * sizeof(body));

    DO_BODY(p, bodytable[ifile], bodytable[ifile]+cmd->nbody) {
        in_real(instr, &Pos(p)[0]);
        in_real(instr, &Pos(p)[1]);
        in_real(instr, &Kappa(p));
        if (scanopt(cmd->options, "kappa-constant"))
            Kappa(p) = 2.0;
    }

    fclose(instr);

//B Find MIN and MAX
    real theta, phi;
    real theta_min, theta_max;
    real phi_min, phi_max;
    p = bodytable[ifile];
    theta_max = theta_min = Pos(p)[0];
    phi_max = phi_min  = Pos(p)[1];

    DO_BODY(p, bodytable[ifile], bodytable[ifile]+cmd->nbody) {
        theta = Pos(p)[0];
        phi = Pos(p)[1];
        theta_min = MIN(theta_min,theta);
        theta_max = MAX(theta_max,theta);
        phi_min = MIN(phi_min,phi);
        phi_max = MAX(phi_max,phi);
    }
    verb_print(cmd->verbose, "\n\tinputdata_AA: min and max of theta = %f %f\n",
               theta_min, theta_max);
    verb_print(cmd->verbose, "\tinputdata_AA: min and max of phi = %f %f\n",
               phi_min, phi_max);
//E

//    real ra, dec;                                   // phi, theta from pix2ang ::
                                                    //  column 2, column 1,
                                                    //  respectively
/*    if (scanopt(cmd->options, "arfken")) {
        DO_BODY(p, bodytable[ifile], bodytable[ifile]+cmd->nbody) {
            ra = Pos(p)[0];
            dec = Pos(p)[1];
            Pos(p)[0] = rsin(dec)*rcos(ra);
            Pos(p)[1] = rsin(dec)*rsin(ra);
            Pos(p)[2] = rcos(dec);
        }
    } else {
        DO_BODY(p, bodytable[ifile], bodytable[ifile]+cmd->nbody) {
            ra = Pos(p)[0];
            dec = Pos(p)[1];
            Pos(p)[0] = rcos(dec)*rcos(ra);
            Pos(p)[1] = rcos(dec)*rsin(ra);
            Pos(p)[2] = rsin(dec);
        }
    } */

//    real theta, phi;
    DO_BODY(p, bodytable[ifile], bodytable[ifile]+cmd->nbody) {
        theta = Pos(p)[0];
        phi = Pos(p)[1];
        spherical_to_cartesians(cmd, gd, theta, phi, Pos(p));
    }

    DO_BODY(p, bodytable[ifile], bodytable[ifile]+cmd->nbody) {
        Type(p) = BODY;
        Mass(p) = mass;
        Weight(p) = weight;
        Id(p) = p-bodytable[ifile]+1;
    }

    return SUCCESS;
}
#endif // ! NDIM == 3

// infileformat: multi-columns-ascii
local int inputdata_ascii_mcolumns(struct cmdline_data* cmd, struct  global_data* gd,
                               string filename, int ifile)
{
    stream instr;
    int ndim;
    bodyptr p;
    real mass=1;
    real weight=1;

    gd->input_comment = "Multi-column position input file";
//
// Upgrade to cross-correlations: kappa-gamma
//
    if ( (scanopt(cmd->options, "pos-and-convergence"))
            && (scanopt(cmd->options, "pos-and-shear")))
        error("inputdata_ascii_mcolumns: mutually excluyent options: %s",
              "'pos-and-convergence' and 'pos-and-shear'");

    int vnpoint;
#if NDIM == 3
    if (scanopt(cmd->options, "only-pos")) {
        InputData_3c(filename, gd->columns[0], gd->columns[1], gd->columns[2], &vnpoint);
    } else {
        if (scanopt(cmd->options, "pos-and-convergence")) {
            InputData_4c(filename,
                         gd->columns[0], gd->columns[1], gd->columns[2],
                         gd->columns[3],
                         &vnpoint);
        } else {
            if (scanopt(cmd->options, "pos-and-shear")) {
                InputData_5c(filename,
                             gd->columns[0], gd->columns[1], gd->columns[2],
                             gd->columns[3], gd->columns[4],
                             &vnpoint);
            }
        }
    }
#else
    if (scanopt(cmd->options, "only-pos")) {
        InputData_2c(filename, gd->columns[0], gd->columns[1], &vnpoint);
    } else {
        if (scanopt(cmd->options, "pos-and-convergence")) {
            InputData_3c(filename,
                         gd->columns[0], gd->columns[1],
                         gd->columns[2],
                         &vnpoint);
        } else {
            if (scanopt(cmd->options, "pos-and-shear")) {
                InputData_4c(filename,
                             gd->columns[0], gd->columns[1],
                             gd->columns[2], gd->columns[3],
                             &vnpoint);
            }
        }
    }
#endif

    ndim = NDIM;
    cmd->nbody = vnpoint;
    gd->nbodyTable[ifile] = cmd->nbody;

    bodytable[ifile] = (bodyptr) allocate(cmd->nbody * sizeof(body));

    DO_BODY(p, bodytable[ifile], bodytable[ifile]+cmd->nbody) {
        Id(p) = p-bodytable[ifile]+1;
#if NDIM == 3
        Pos(p)[0] = inout_xval[Id(p)-1];
        Pos(p)[1] = inout_yval[Id(p)-1];
        Pos(p)[2] = inout_zval[Id(p)-1];
        if (scanopt(cmd->options, "pos-and-convergence"))
            Kappa(p) = inout_uval[Id(p)-1];
        if (scanopt(cmd->options, "pos-and-shear")) {
            Gamma1(p) = inout_uval[Id(p)-1];
            Gamma2(p) = inout_vval[Id(p)-1];
        }
#else
        Pos(p)[0] = inout_xval[Id(p)-1];
        Pos(p)[1] = inout_yval[Id(p)-1];
        if (scanopt(cmd->options, "pos-and-convergence"))
            Kappa(p) = inout_zval[Id(p)-1];
        if (scanopt(cmd->options, "pos-and-shear")) {
            Gamma1(p) = inout_zval[Id(p)-1];
            Gamma2(p) = inout_uval[Id(p)-1];
        }
#endif
    }

#if NDIM == 3
    free(inout_xval);
    free(inout_yval);
    free(inout_zval);
    free(inout_uval);
    if (scanopt(cmd->options, "pos-and-shear"))
        free(inout_vval);
#else
    free(inout_xval);
    free(inout_yval);
    free(inout_zval);
    if (scanopt(cmd->options, "pos-and-shear"))
        free(inout_uval);
#endif

//B
//B Set (0,0,...) as the center of the box
// By now it is only working with boxes centered at (0,0,...)
    cellptr root;                                   // Set it up a temporal root
    root = (cellptr) allocate(1 * sizeof(body));

    findRootCenter(cmd, gd, bodytable[ifile], gd->nbodyTable[ifile], ifile, root);
    centerBodies(bodytable[ifile], gd->nbodyTable[ifile], ifile, root);
    findRootCenter(cmd, gd, bodytable[ifile], gd->nbodyTable[ifile], ifile, root);
    CLRV(Pos(root));
//E
    gd->rSizeTable[ifile] = 1.0;
    expandbox(cmd, gd, bodytable[ifile], gd->nbodyTable[ifile], ifile, root);
    free(root);
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
                 "\tinputdata_ascii_mcolumns: nbody and ndim: %d %d...\n",
                 cmd->nbody, ndim);
    verb_print_q(2, cmd->verbose,
                 "\tinputdata_ascii_mcolumns: Lx, Ly, Lz: %g %g %g...\n",
               gd->Box[0], gd->Box[1], gd->Box[2]);

    real kavg=0.0;
    DO_BODY(p, bodytable[ifile], bodytable[ifile]+cmd->nbody) {
        Type(p) = BODY;
        Mass(p) = mass;
        Weight(p) = weight;
        if (scanopt(cmd->options, "kappa-constant"))
            Kappa(p) = 2.0;
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
        error("inputdata_ascii_mcolumns: at least two bodies have same position\n");
    }
    //E

    return SUCCESS;
}


#endif	// ! _cballsio_iolib_02_h
