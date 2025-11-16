
//=============================================================================
//        1          2          3          4        ^ 5          6          7

#ifndef _cballsio_cfitsio_02_h
#define _cballsio_cfitsio_02_h


// input in: addons/source/cballsio/cballsio_include_11a.h

local int inputdata_cfitsio(struct cmdline_data* cmd, struct  global_data* gd,
                               string filename, int ifile)
{
    fitsfile *fptr;
    char card[FLEN_CARD];
    int status = 0, nkeys, ii;                      // MUST initialize status

    gd->input_comment = "fits input file";

    //B Was CFITSIO compiled with the -D_REENTRANT flag?  1 = yes, 0 = no.
    verb_print(cmd->verbose,
               "\tinputdata_cfitsio: -D_REENTRANT flag: %d...\n",
               fits_is_reentrant());
    //E

    verb_print(cmd->verbose,
               "\tinputdata_cfitsio: opening fits file: %s...\n",
               filename);
    if (scanopt(cmd->options, "fits-type-file"))
        fits_open_file(&fptr, filename, READONLY, &status);
    else
        fits_open_data(&fptr, filename, READONLY, &status);
    if (status) {                                   // print any error messages
        verb_print(cmd->verbose,
                   "\tinputdata_cfitsio: open status: %d...\n\n",
                   status);
        fits_report_error(stderr, status);
        if (status==104) exit(1);
    }

    fits_get_hdrspace(fptr, &nkeys, NULL, &status);
    if (status) {                                   // print any error messages
        verb_print(cmd->verbose,
                   "\tinputdata_cfitsio: get_hdrspace status: %d...\n\n",
                   status);
        fits_report_error(stderr, status);
    }
    fits_get_num_rows(fptr, &cmd->nbody, &status);
    if (status) {                                   // print any error messages
        verb_print(cmd->verbose,
                   "\tinputdata_cfitsio: get_num_rows status: %d...\n\n",
                   status);
        fits_report_error(stderr, status);
    }

    if (scanopt(cmd->options, "header-info")){
        verb_print(cmd->verbose,"\nHeader information:\n\n");
        for (ii = 1; ii <= nkeys; ii++) {
            fits_read_record(fptr, ii, card, &status); // read keyword
            printf("%s\n", card);
        }
        printf("END\n\n");                          // terminate listing
                                                    //  with END
        fits_get_num_rows(fptr, &cmd->nbody, &status);
        int ncols;
        fits_get_num_cols(fptr, &ncols, &status);
        verb_print(cmd->verbose,
            "\tinputdata_cfitsio: nbody = %d... and number of columns = %d\n\n",
            cmd->nbody, ncols);
        if (cmd->nbody < 1)
            error("tinputdata_cfitsio:: nbody = %d is absurd\n", cmd->nbody);

        int typecode;
        long repeat;
        long width;
        verb_print(cmd->verbose,
                   "Columns info details:\n");
        for (ii = 1; ii <= ncols; ii++) {
            fits_get_coltype(fptr, ii, &typecode,
                             &repeat, &width, &status);
            switch(typecode) {
                case TLONG:
                    verb_print(cmd->verbose,
                               "%d: typecode, repeat, width = %d %s %ld %ld\n",
                               ii, typecode, "TLONG", repeat, width);
                    break;
                case TFLOAT:
                    verb_print(cmd->verbose,
                               "%d: typecode, repeat, width = %d %s %ld %ld\n",
                               ii, typecode, "TFLOAT", repeat, width);
                    break;
                case TDOUBLE:
                    verb_print(cmd->verbose,
                               "%d: typecode, repeat, width = %d %s %ld %ld\n",
                               ii, typecode, "TDOUBLE", repeat, width);
                    break;
            }
        }
        verb_print(cmd->verbose,"\n");
        verb_print(cmd->verbose,"end of header information.\n\n");
    } else { // ! header-info
        for (ii = 1; ii <= nkeys; ii++) {
            fits_read_record(fptr, ii, card, &status); // read keyword
        }
        fits_get_num_rows(fptr, &cmd->nbody, &status);
        int ncols;
        fits_get_num_cols(fptr, &ncols, &status);
        verb_print(cmd->verbose,
            "\tinputdata_cfitsio: nbody = %d... and number of columns = %d\n",
            cmd->nbody, ncols);
        if (cmd->nbody < 1)
            error("tinputdata_cfitsio:: nbody = %d is absurd\n", cmd->nbody);
    }

    if (scanopt(cmd->options, "stop-fits")) {
        if (strnull(cmd->outfile)) {
            fits_close_file(fptr, &status);
            if (status) {                           // print error message
                verb_print(cmd->verbose,
                           "\tinputdata_cfitsio: status: %d...\n\n",
                           status);
                fits_report_error(stderr, status);
            }
            exit(1);
        }
    }

    inputdata_cfitsio_xyz(cmd, gd, filename, ifile, fptr);

    //B once data has been read, close fits file
    verb_print(cmd->verbose,
               "\tinputdata_cfitsio: closing fits file: %s...\n",
               filename);
    fits_close_file(fptr, &status);
    if (status) {                                   // print error message
        verb_print(cmd->verbose,
                   "\tinputdata_cfitsio: status: %d...\n\n",
                   status);
        fits_report_error(stderr, status);
    }
    //E

    return SUCCESS;
}

// Routine to read ra-dec healpix fits files as given by TC
local int inputdata_cfitsio_radec_tc(struct cmdline_data* cmd,
                                             struct  global_data* gd,
                                             string filename, int ifile)
{
    fitsfile *fptr;
    char card[FLEN_CARD];
    int status = 0, nkeys, ii;                      // MUST init status

    gd->input_comment = "fits-healpix-radec-tc input file";

    //B Was CFITSIO compiled with the -D_REENTRANT flag?  1 = yes, 0 = no.
    verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                          "\tinputdata_cfitsio: -D_REENTRANT flag: %d...\n",
                          fits_is_reentrant());
    //E

    verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                        "\tinputdata_cfitsio: opening fits file: %s...\n",
                        filename);

    if (scanopt(cmd->options, "fits-type-file"))
        fits_open_file(&fptr, filename, READONLY, &status);
    else
        fits_open_data(&fptr, filename, READONLY, &status);
    if (status) {                                   // print error message
        verb_print(cmd->verbose,
                   "\tinputdata_cfitsio: open status: %d...\n\n",
                   status);
        fits_report_error(stderr, status);
        if (status==104) exit(1);
    }

    fits_get_hdrspace(fptr, &nkeys, NULL, &status);
    if (status) {                                   // print error message
        verb_print(cmd->verbose,
                   "\tinputdata_cfitsio: get_hdrspace status: %d...\n\n",
                   status);
        fits_report_error(stderr, status);
    }
    fits_get_num_rows(fptr, &cmd->nbody, &status);
    if (status) {                                   // print error message
        verb_print(cmd->verbose,
                   "\tinputdata_cfitsio: get_num_rows status: %d...\n\n",
                   status);
        fits_report_error(stderr, status);
    }

    if (scanopt(cmd->options, "header-info")){
        verb_print(cmd->verbose,"\nHeader information:\n\n");
        for (ii = 1; ii <= nkeys; ii++) {
            fits_read_record(fptr, ii, card, &status); // read keyword
            printf("%s\n", card);
        }
        printf("END\n\n");                          // terminate listing
                                                    //  with END
        fits_get_num_rows(fptr, &cmd->nbody, &status);
        int ncols;
        fits_get_num_cols(fptr, &ncols, &status);
        verb_print(cmd->verbose,
            "\tinputdata_cfitsio: nbody = %d... and number of columns = %d\n\n",
            cmd->nbody, ncols);
        if (cmd->nbody < 1)
            error("tinputdata_cfitsio:: nbody = %d is absurd\n", cmd->nbody);

        int typecode;
        long repeat;
        long width;
        verb_print(cmd->verbose,
                   "Columns info details:\n");
        for (ii = 1; ii <= ncols; ii++) {
            fits_get_coltype(fptr, ii, &typecode,
                             &repeat, &width, &status);
            switch(typecode) {
                case TLONG:
                    verb_print(cmd->verbose,
                               "%d: typecode, repeat, width = %d %s %ld %ld\n",
                               ii, typecode, "TLONG", repeat, width);
                    break;
                case TFLOAT:
                    verb_print(cmd->verbose,
                               "%d: typecode, repeat, width = %d %s %ld %ld\n",
                               ii, typecode, "TFLOAT", repeat, width);
                    break;
                case TDOUBLE:
                    verb_print(cmd->verbose,
                               "%d: typecode, repeat, width = %d %s %ld %ld\n",
                               ii, typecode, "TDOUBLE", repeat, width);
                    break;
            }
        }
        verb_print(cmd->verbose,"\n");
        verb_print(cmd->verbose,"end of header information.\n\n");
    } else { // ! header-info
        for (ii = 1; ii <= nkeys; ii++) {
            fits_read_record(fptr, ii, card, &status); // read keyword
        }
        fits_get_num_rows(fptr, &cmd->nbody, &status);
        int ncols;
        fits_get_num_cols(fptr, &ncols, &status);
        verb_print(cmd->verbose,
            "\tinputdata_cfitsio: nbody = %d... and number of columns = %d\n",
            cmd->nbody, ncols);
        if (cmd->nbody < 1)
            error("tinputdata_cfitsio:: nbody = %d is absurd\n", cmd->nbody);
    }

    if (scanopt(cmd->options, "stop-fits")) {
        if (strnull(cmd->outfile)) {
            fits_close_file(fptr, &status);
            if (status) {                           // print any error messages
                verb_print(cmd->verbose,
                           "\tinputdata_cfitsio: status: %d...\n\n",
                           status);
                fits_report_error(stderr, status);
            }
            exit(1);
        }
    }

    inputdata_cfitsio_ra_dec(cmd, gd, filename, ifile, fptr);

    //B once data has been read, close fits file
    verb_print(cmd->verbose,
               "\tinputdata_cfitsio: closing fits file: %s...\n",
               filename);
    fits_close_file(fptr, &status);
    if (status) {                                   // print any error messages
        verb_print(cmd->verbose,
                   "\tinputdata_cfitsio: status: %d...\n\n",
                   status);
        fits_report_error(stderr, status);
    }
    //E

    return SUCCESS;
}

// Routine to read ra-dec fits files
//  by default:
//      ra is column 2 and it is float type
//      dec is column 3 and it is float type
//      kappa is column 8 and it is float type
//  Check header with:
//      options=header-info,stop,ra-dec
local int inputdata_cfitsio_ra_dec(struct cmdline_data* cmd,
                                   struct  global_data* gd,
                                   string filename, int ifile, fitsfile *fptr)
{
    string routineName = "searchcalc_octree_ggg_omp";
    bodyptr p;
    real mass=1;
    real weight=1;

    int datatype;
    int colnum;
    LONGLONG firstrow;
    LONGLONG firstelem;
    LONGLONG nelements;
    float nulval;
    float *arrayKappa;
    int anynul;
    int status = 0;

    datatype = 42;                                  // TFLOAT
    colnum = 8;
    firstrow = 1;
    firstelem = 1;
    nelements = cmd->nbody;
    arrayKappa = (float*) allocate(cmd->nbody * sizeof(float));
    fits_read_col(fptr, datatype, colnum, firstrow, firstelem,
                  nelements, &nulval, arrayKappa, &anynul, &status);
    bodytable[ifile] = (bodyptr) allocate(cmd->nbody * sizeof(body));
    DO_BODY(p, bodytable[ifile], bodytable[ifile]+cmd->nbody) {
        Kappa(p) = arrayKappa[p-bodytable[ifile]];
        if (scanopt(cmd->options, "kappa-constant")) {
            if (scanopt(cmd->options, "kappa-constant-one"))
                Kappa(p) = 1.0;
            else
                Kappa(p) = 2.0;
        }
    }

    free(arrayKappa);

#if NDIM == 3
        real ra, dec;
        float *arrayRA;
        float *arrayDEC;
        arrayRA = (float*) allocate(cmd->nbody * sizeof(float));
        arrayDEC = (float*) allocate(cmd->nbody * sizeof(float));
        colnum = 2;
        fits_read_col(fptr, datatype, colnum, firstrow, firstelem,
                      nelements, &nulval, arrayRA, &anynul, &status);
        if (status) {                               // print any error messages
            verb_print(cmd->verbose,
                       "\tinputdata_cfitsio: status: %d...\n\n",
                       status);
            fits_report_error(stderr, status);
        }
        colnum = 3;
        fits_read_col(fptr, datatype, colnum, firstrow, firstelem,
                      nelements, &nulval, arrayDEC, &anynul, &status);
        if (status) {                               // print any error messages
            verb_print(cmd->verbose,
                       "\tinputdata_cfitsio: status: %d...\n\n",
                       status);
            fits_report_error(stderr, status);
        }
        if (scanopt(cmd->options, "no-arfken")) {
            DO_BODY(p, bodytable[ifile], bodytable[ifile]+cmd->nbody) {
                ra = arrayRA[p-bodytable[ifile]];
                dec = arrayDEC[p-bodytable[ifile]];
                Pos(p)[0] = rcos(dec)*rcos(ra);
                Pos(p)[1] = rcos(dec)*rsin(ra);
                Pos(p)[2] = rsin(dec);
            }
        } else {
            DO_BODY(p, bodytable[ifile], bodytable[ifile]+cmd->nbody) {
                ra = arrayRA[p-bodytable[ifile]];
                dec = arrayDEC[p-bodytable[ifile]];
                Pos(p)[0] = rsin(dec)*rcos(ra);
                Pos(p)[1] = rsin(dec)*rsin(ra);
                Pos(p)[2] = rcos(dec);
            }
        }

    free(arrayRA);
    free(arrayDEC);

    gd->nbodyTable[ifile] = cmd->nbody;

    real kavg=0.0;
    DO_BODY(p, bodytable[ifile], bodytable[ifile]+gd->nbodyTable[ifile]) {
        Type(p) = BODY;
        Mass(p) = mass;
        Weight(p) = weight;
        Id(p) = p-bodytable[ifile]+1;
        kavg += Kappa(p);
    }
    verb_print(cmd->verbose,
               "\tinputdata_cfitsio: average of kappa (%ld particles) = %le\n",
               cmd->nbody, kavg/((real)gd->nbodyTable[ifile]) );

    //B If needed locate particles with same position.
    //  This is a slow process, use if it is necessary...
    if (scanopt(cmd->options, "check-eq-pos")) {
        bodyptr q;
        real dist2;
        vector distv;
        bool flag=0;
        int k;
        verb_print(cmd->verbose,
                   "inputdata_cfitsio: checking eq-pos...");
        DO_BODY(p, bodytable[ifile], bodytable[ifile]+cmd->nbody-1)
            DO_BODY(q, p+1, bodytable[ifile]+cmd->nbody)
                if (p != q) {
                    DOTPSUBV(dist2, distv, Pos(p), Pos(q));
                    if (dist2 == 0.0) {
                        verb_print(cmd->verbose,
                        "\nIds: %ld and %ld have the same position\n",
                                   Id(p),Id(q));
                        DO_COORD(k)
                            verb_print(cmd->verbose,"Pos[k]: %le %le\n",
                                       Pos(p)[k],Pos(q)[k]);
                        flag=1;
                    }
                }
        verb_print(cmd->verbose,
                   "inputdata_cfitsio: done.\n");
        if (flag)
        error("inputdata_cfitsio: at least two bodies have same position\n");
    }
    //E
#else
// two dimension (ra,dec) is comming soon!
//#error `DEFDIMENSION` is not set to 3. Do so in Makefile_settings
    verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
        "\n%s: `DEFDIMENSION` is not set to 3. Do so in Makefile_settings\n",
        routineName);

#endif
}

local int inputdata_cfitsio_xyz(struct cmdline_data* cmd,
                                struct  global_data* gd,
                                string filename, int ifile, fitsfile *fptr)
{
    bodyptr p;
    real mass=1;
    real weight=1;

    int datatype;
    int colnum;
    LONGLONG firstrow;
    LONGLONG firstelem;
    LONGLONG nelements;
    double nulval;
    double *arrayKappa;
    int anynul;
    int status = 0;

    datatype = 82;                                  // double
    firstrow = 1;
    firstelem = 1;
    nelements = cmd->nbody;
    arrayKappa = (double*) allocate(cmd->nbody * sizeof(double));

    colnum = gd->columns[3];                        // kappa
    verb_print(cmd->verbose,
               "\tinputdata_cfitsio: colnum: %d...\n",
               colnum);
    fits_read_col(fptr, datatype, colnum, firstrow, firstelem,
                  nelements, &nulval, arrayKappa, &anynul, &status);

    bodytable[ifile] = (bodyptr) allocate(cmd->nbody * sizeof(body));
    DO_BODY(p, bodytable[ifile], bodytable[ifile]+cmd->nbody) {
        Kappa(p) = arrayKappa[p-bodytable[ifile]];
        if (scanopt(cmd->options, "kappa-constant")) {
            if (scanopt(cmd->options, "kappa-constant-one"))
                Kappa(p) = 1.0;
            else
                Kappa(p) = 2.0;
        }
    }
    free(arrayKappa);

#if NDIM == 3
    double *arrayX;
    double *arrayY;
    double *arrayZ;
    arrayX = (double*) allocate(cmd->nbody * sizeof(double));
    arrayY = (double*) allocate(cmd->nbody * sizeof(double));
    arrayZ = (double*) allocate(cmd->nbody * sizeof(double));

    colnum = gd->columns[0];
    verb_print(cmd->verbose,
               "\tinputdata_cfitsio: colnum X: %d...\n",
               colnum);
    fits_read_col(fptr, datatype, colnum, firstrow, firstelem,
                  nelements, &nulval, arrayX, &anynul, &status);
    if (status) {                                   // print any error messages
        verb_print(cmd->verbose,
                   "\tinputdata_cfitsio: status: %d...\n\n",
                   status);
        fits_report_error(stderr, status);
    }
    
    colnum = gd->columns[1];
    verb_print(cmd->verbose,
               "\tinputdata_cfitsio: colnum Y: %d...\n",
               colnum);
    fits_read_col(fptr, datatype, colnum, firstrow, firstelem,
                  nelements, &nulval, arrayY, &anynul, &status);
    if (status) {                                   // print any error messages
        verb_print(cmd->verbose,
                   "\tinputdata_cfitsio: status: %d...\n\n",
                   status);
        fits_report_error(stderr, status);
    }

    colnum = gd->columns[2];
    verb_print(cmd->verbose,
               "\tinputdata_cfitsio: colnum Z: %d...\n\n",
               colnum);
    fits_read_col(fptr, datatype, colnum, firstrow, firstelem,
                  nelements, &nulval, arrayZ, &anynul, &status);
    if (status) {                                   // print any error messages
        verb_print(cmd->verbose,
                   "\tinputdata_cfitsio: status: %d...\n\n",
                   status);
        fits_report_error(stderr, status);
    }
    DO_BODY(p, bodytable[ifile], bodytable[ifile]+cmd->nbody) {
        Pos(p)[0] = arrayX[p-bodytable[ifile]];
        Pos(p)[1] = arrayY[p-bodytable[ifile]];
        Pos(p)[2] = arrayZ[p-bodytable[ifile]];
    }

    free(arrayX);
    free(arrayY);
    free(arrayZ);

    gd->nbodyTable[ifile] = cmd->nbody;

    real kavg=0.0;
    DO_BODY(p, bodytable[ifile], bodytable[ifile]+gd->nbodyTable[ifile]) {
        Type(p) = BODY;
        Mass(p) = mass;
        Weight(p) = weight;
        Id(p) = p-bodytable[ifile]+1;
        kavg += Kappa(p);
    }
    verb_print(cmd->verbose,
               "\tinputdata_cfitsio: average of kappa (%ld particles) = %le\n",
               cmd->nbody, kavg/((real)gd->nbodyTable[ifile]) );

    //B If needed locate particles with same position.
    //  This is a very slow process, use if it necessary...
    if (scanopt(cmd->options, "check-eq-pos")) {
        bodyptr q;
        real dist2;
        vector distv;
        bool flag=0;
        int k;
        verb_print(cmd->verbose,
                   "inputdata_cfitsio: checking eq-pos...");
        DO_BODY(p, bodytable[ifile], bodytable[ifile]+cmd->nbody-1)
            DO_BODY(q, p+1, bodytable[ifile]+cmd->nbody)
                if (p != q) {
                    DOTPSUBV(dist2, distv, Pos(p), Pos(q));
                    if (dist2 == 0.0) {
                        verb_print(cmd->verbose,
                                   "\nIds: %ld and %ld have the same position\n",
                                   Id(p),Id(q));
                        DO_COORD(k)
                            verb_print(cmd->verbose,"Pos[k]: %le %le\n",
                                       Pos(p)[k],Pos(q)[k]);
                        flag=1;
                    }
                }
        verb_print(cmd->verbose,
                   "inputdata_cfitsio: done.\n");
        if (flag)
            error("inputdata_cfitsio: %s\n",
                  "at least two bodies have same position");
    }
    //E
#else
    error("inputdata_cfitsio_xyz: %s\n\n",
        "this routine works only in 3D. exiting...");
//#error `DEFDIMENSION` is not set to 3. Do so in Makefile_settings
#endif
}

// Routine to read fits-healpix files
//  so far only RING scheme
local int inputdata_cfitsio_healpix(struct cmdline_data* cmd,
                                    struct  global_data* gd,
                                    string filename, int ifile)
{
    string routineName = "inputdata_cfitsio_healpix";
    fitsfile *fptr;
    char card[FLEN_CARD];
    int status = 0, nkeys, ii;                      // MUST initialize status

    gd->input_comment = "fits-healpix input file";

    //B Was CFITSIO compiled with the -D_REENTRANT flag?  1 = yes, 0 = no.
    verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                          "\t%s: -D_REENTRANT flag: %d...\n",
                          routineName,fits_is_reentrant());
    //E

    verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                          "\t%s: opening fits file: %s...\n",
                          routineName,filename);
    if (scanopt(cmd->options, "fits-type-file"))
        fits_open_file(&fptr, filename, READONLY, &status);
    else
        fits_open_data(&fptr, filename, READONLY, &status);
    if (status) {                                   // print any error messages
        verb_print(cmd->verbose,
                   "\t%s: open status: %d...\n\n",
                   routineName,status);
        fits_report_error(stderr, status);
        if (status==104) exit(1);
    }

    fits_get_hdrspace(fptr, &nkeys, NULL, &status);
    if (status) {                                   // print any error messages
        verb_print(cmd->verbose,
                   "\t%s: get_hdrspace status: %d...\n\n",
                   routineName,status);
        fits_report_error(stderr, status);
    }
    fits_get_num_rows(fptr, &cmd->nbody, &status);
    if (status) {                                   // print any error messages
        verb_print(cmd->verbose,
                   "\tinputdata_cfitsio: get_num_rows status: %d...\n\n",
                   status);
        fits_report_error(stderr, status);
    }

    if (scanopt(cmd->options, "header-info")){
        verb_print(cmd->verbose,"\nHeader information:\n\n");
        for (ii = 1; ii <= nkeys; ii++) {
            fits_read_record(fptr, ii, card, &status); // read keyword
            printf("%s\n", card);
        }
        printf("END\n\n");                          // terminate listing
                                                    //  with END
        fits_get_num_rows(fptr, &cmd->nbody, &status);
        int ncols;
        fits_get_num_cols(fptr, &ncols, &status);
        verb_print(cmd->verbose,
            "\tinputdata_cfitsio: nbody (nrows) = %d... %s = %d\n\n",
            cmd->nbody, "and number of columns", ncols);
        if (cmd->nbody < 1)
            error("tinputdata_cfitsio:: nbody = %d is absurd\n", cmd->nbody);

        int typecode;
        long repeat;
        long width;
        verb_print(cmd->verbose,
                   "Columns info details:\n");
        for (ii = 1; ii <= ncols; ii++) {
            fits_get_coltype(fptr, ii, &typecode,
                             &repeat, &width, &status);
            switch(typecode) {
                case TLONG:
                    verb_print(cmd->verbose,
                               "%d: typecode, repeat, width = %d %s %ld %ld\n",
                               ii, typecode, "TLONG", repeat, width);
                    break;
                case TFLOAT:
                    verb_print(cmd->verbose,
                               "%d: typecode, repeat, width = %d %s %ld %ld\n",
                               ii, typecode, "TFLOAT", repeat, width);
                    break;
                case TDOUBLE:
                    verb_print(cmd->verbose,
                               "%d: typecode, repeat, width = %d %s %ld %ld\n",
                               ii, typecode, "TDOUBLE", repeat, width);
                    break;
            }
        }
        verb_print(cmd->verbose,"\n");
        verb_print(cmd->verbose,"end of header information.\n\n");
    } else { // ! header-info
        for (ii = 1; ii <= nkeys; ii++) {
            fits_read_record(fptr, ii, card, &status); // read keyword
        }
        fits_get_num_rows(fptr, &cmd->nbody, &status);
        int ncols;
        fits_get_num_cols(fptr, &ncols, &status);
        verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                    "\t%s: nbody(nrows) = %d... %s = %d\n",
                    routineName, cmd->nbody, "and number of columns", ncols);
        if (cmd->nbody < 1)
            error("\tinputdata_cfitsio:: nbody = %d is absurd\n", cmd->nbody);
    }

    if (scanopt(cmd->options, "stop-fits")) {
        if (strnull(cmd->outfile)) {
            fits_close_file(fptr, &status);
            if (status) {                           // print any error messages
                verb_print(cmd->verbose,
                           "\tinputdata_cfitsio: status: %d...\n\n",
                           status);
                fits_report_error(stderr, status);
            }
            exit(1);
        }
    }

    if (scanopt(cmd->options, "read-mask")) {
        if (ifile == 0) {                           // needs mask file be second... (?)
            inputdata_cfitsio_healpix_map(cmd, gd, filename, ifile, fptr);
        } else {
            if (ifile != 1)
                error("\t%s: read-mask ifile = %d is absurd\n",
                      routineName, ifile);
            inputdata_cfitsio_healpix_map_mask(cmd, gd, filename, ifile, fptr);
        }
    } else {
        if (scanopt(cmd->options, "mask-inside")) {
            inputdata_cfitsio_healpix_map_mask_inside(cmd, gd,
                                                      filename, ifile, fptr);
        } else {
            inputdata_cfitsio_healpix_map(cmd, gd, filename, ifile, fptr);
        }
    }

    //B once data has been read, close fits file
    verb_print(cmd->verbose,
               "\t%s: closing fits file: %s...\n",
               routineName, filename);
    fits_close_file(fptr, &status);
    if (status) {                                   // print any error messages
        verb_print(cmd->verbose,
                   "\t%s: status: %d...\n\n",
                   routineName, status);
        fits_report_error(stderr, status);
    }
    //E

    return SUCCESS;
}

// reading healpix map
#if THREEDIMCODE
//  3D only... not anymore
local int inputdata_cfitsio_healpix_map(struct cmdline_data* cmd,
                                        struct  global_data* gd,
                                        string filename,
                                        int ifile, fitsfile *fptr)
{
    string routineName = "inputdata_cfitsio_healpix_map";
    bodyptr p;
    real mass=1;
    real weight=1;

    long ipix;
    double theta;
    double phi;
    double thetamin, thetamax;
    double phimin, phimax;

    float *map;
    long npixel, nside;
    char order1[10];
    char order2[10];
    char coord[10];

#if THREEDIMCODE
    verb_print(cmd->verbose, "\nWorking 3D map...\n");
#else
    error("\nOnly 3D is implemented so far... exiting...\n\n")
#endif

    npixel = get_fits_size(filename, &nside, order1);
    verb_print(cmd->verbose,
        "filename, ifile, nside, npixel, order1:");
    verb_print(cmd->verbose,
        "%s %d %ld %ld %s\n", filename, ifile, nside, npixel, order1);
    map = read_healpix_map(filename, &nside, coord, order2);
    verb_print(cmd->verbose,
               "\nAllocated %g MByte for temporal map pixel (%ld) storage.\n",
               npixel*sizeof(float)/(1024.0*1024.0),
               npixel);

    if (strcmp(order1,order2)!=0) {printf("Error: Bad ordering\n");}
    verb_print(cmd->verbose,
               "%s: nbody = %d...\n",
               routineName, npixel);
    if (npixel < 1)
        error("%s: npixel = %d is absurd\n", routineName, npixel);

    bodyptr bodytabtmp;
    bodytabtmp = (bodyptr) allocate(npixel * sizeof(body));
    verb_print(cmd->verbose,
               "\nAllocated %g MByte for temporal particle (%ld) storage.\n",
               npixel*sizeof(body)/(1024.0*1024.0),
               npixel);

    //B save optional RA-DEC to a file... 3-columns: RA, DEC, Kappa
    char namebuf[256];
    stream outstr;
    if (scanopt(cmd->options, "save-ra-dec")&&!strnull(cmd->outfile)) {
        sprintf(namebuf, gd->fpfnameOutputFileName);
        outstr = stropen(namebuf, "w!");
    }
    //E

    cmd->nbody=npixel;

    real xmin, ymin, zmin;
    real xmax, ymax, zmax;
    xmin=0., ymin=0., zmin=0.;
    xmax=0., ymax=0., zmax=0.;

    INTEGER i;
    INTEGER iselect = 0;
    for(ipix=0;ipix<npixel;ipix++) {                // RING loop order
        p = bodytabtmp+ipix;
        Update(p) = FALSE;
        Mask(p) = 1;                                // initialize body's Mask
        pix2ang_ring(nside, ipix, &theta, &phi);
        if (scanopt(cmd->options, "patch")) {
            if (cmd->thetaL < theta && theta < cmd->thetaR) {
                if (cmd->phiL < phi && phi < cmd->phiR) {
                    iselect++;
                    coordinate_transformation(cmd, gd, theta, phi, Pos(p));
                    if (!scanopt(cmd->options, "kappa-constant"))
                        Kappa(p) = map[ipix];
                    else {
                        Kappa(p) = 2.0;
                        if (scanopt(cmd->options, "kappa-constant-one"))
                            Kappa(p) = 1.0;
                    }
                    //B save optional RA-DEC to a file...
                    if (scanopt(cmd->options, "save-ra-dec")
                        &&!strnull(cmd->outfile)) {
                        out_real_mar(outstr, phi);
                        out_real_mar(outstr, theta);
                        out_real(outstr, Kappa(p));
                    }
                    //E
                    Type(p) = BODY;
                    Mass(p) = mass;
                    Weight(p) = weight;
                    Id(p) = p-bodytabtmp+iselect;
                    Update(p) = TRUE;
                    xmin = Pos(p)[0];
                    ymin = Pos(p)[1];
                    zmin = Pos(p)[2];
                    xmax = Pos(p)[0];
                    ymax = Pos(p)[1];
                    zmax = Pos(p)[2];
                }
            }
        } else { // ! all
            iselect++;
            coordinate_transformation(cmd, gd, theta, phi, Pos(p));
            if (!scanopt(cmd->options, "kappa-constant"))
                Kappa(p) = map[ipix];
            else {
                Kappa(p) = 2.0;
                if (scanopt(cmd->options, "kappa-constant-one"))
                    Kappa(p) = 1.0;
            }
            //B save optional RA-DEC to a file...
            if (scanopt(cmd->options, "save-ra-dec")
                &&!strnull(cmd->outfile)) {
                out_real_mar(outstr, phi);
                out_real_mar(outstr, theta);
                out_real(outstr, Kappa(p));
            }
            //E
            Type(p) = BODY;
            Mass(p) = mass;
            Weight(p) = weight;
            Id(p) = p-bodytabtmp+iselect;
            Update(p) = TRUE;
            xmin = Pos(p)[0];
            ymin = Pos(p)[1];
            zmin = Pos(p)[2];
            xmax = Pos(p)[0];
            ymax = Pos(p)[1];
            zmax = Pos(p)[2];
            
            //B correction 2025-05-03 :: look for edge-effects
            // activate a flag for this catalog, that it is using patch-with-all
            //  and use it in EvalHist routine...
    #if defined(NMultipoles) && defined(NONORMHIST)
            if (scanopt(cmd->options, "patch-with-all")) {
                UpdatePivot(p) = TRUE;               
                if (cmd->thetaL < theta && theta < cmd->thetaR) {
                    if (cmd->phiL < phi && phi < cmd->phiR) {
                        UpdatePivot(p) = TRUE;
                        gd->pivotCount += 1;
                    } else {
                        UpdatePivot(p) = FALSE;
                    }
                } else {
                    UpdatePivot(p) = FALSE;
                }
            }
    #endif
            //E

        } // ! all
    } // ! end loop ipix

    //B save optional RA-DEC to a file...
    if (scanopt(cmd->options, "save-ra-dec")&&!strnull(cmd->outfile)
                            &&scanopt(cmd->options, "stop")) {
        fclose(outstr);
        exit(1);            // need to free all mem allocation...
    }
    //E

    //B correction 2025-05-03 :: look for edge-effects
#if defined(NMultipoles) && defined(NONORMHIST)
    if (scanopt(cmd->options, "patch-with-all")) {
        verb_print(cmd->verbose,
                   "\n%s: total number of pixels to be pivots: %ld\n",
                   routineName, gd->pivotCount);
    }
#endif
    //E

    free(map);
    verb_print(cmd->verbose,
               "\nFreed %g MByte for temporal map pixel (%ld) storage.\n",
               npixel*sizeof(float)/(1024.0*1024.0),
               npixel);

    bodyptr q;
    if (scanopt(cmd->options, "patch"))
        cmd->nbody = iselect;

    gd->nbodyTable[ifile] = cmd->nbody;
    bodytable[ifile] = (bodyptr) allocate(cmd->nbody * sizeof(body));

    real kavg = 0;
    INTEGER ij=0;
    for(ipix=0;ipix<npixel;ipix++){
        q = bodytabtmp+ipix;
        if(Update(q)) {
            p = bodytable[ifile]+ij;
            Pos(p)[0] = Pos(q)[0];
            Pos(p)[1] = Pos(q)[1];
            Pos(p)[2] = Pos(q)[2];
            Kappa(p) = Kappa(q);
            Type(p) = Type(q);
            Mass(p) = mass;
            Weight(p) = weight;
            Id(p) = p-bodytable[ifile]+ipix;
            Mask(p) = Mask(q);
            xmin = MIN(xmin,Pos(p)[0]);
            ymin = MIN(ymin,Pos(p)[1]);
            zmin = MIN(zmin,Pos(p)[2]);
            xmax = MAX(xmax,Pos(p)[0]);
            ymax = MAX(ymax,Pos(p)[1]);
            zmax = MAX(zmax,Pos(p)[2]);
            ij++;
            kavg += Kappa(p);
        }
    }

    char file[180] = "outputmap.fits" ;
    char fileforce[180] ;
    if (scanopt(cmd->options, "plot-map-gif")) {
        //B add ifile tag...
        sprintf(fileforce, "!%s/%s",                // leading !
                gd->outputDir, file);               //  to allow overwrite
        //E
        verb_print(cmd->verbose,
                   "\t%s: %s %s...\n",
                   routineName, "\n\t\tsaving map to a fits file:", fileforce);
        map = (float *)malloc(npixel*sizeof(float));
        for(ipix=0;ipix<npixel;ipix++){             // all pixels in the sphere are filled
            q = bodytabtmp+ipix;
            if(Update(q)) {
                map[ipix] = Kappa(q);
            } else {
                map[ipix] = 0.0;
            }
        }
        //B write healpix map in RING order: "0"
        //      and
        //      "G = Galactic, E = ecliptic, C = celestial = equatorial"
        write_healpix_map(map, nside, fileforce, 0, "C");  // choose COORDSYS...
        fprintf(stdout,"\t\tfile written\n");
        free(map);
    }

    verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                          "\n\t%s: min and max of x = %f %f\n",
                          routineName, xmin, xmax);
    verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                          "\t%s: min and max of y = %f %f\n",
                          routineName, ymin, ymax);
    verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                          "\t%s: min and max of z = %f %f\n",
                          routineName, zmin, zmax);

    free(bodytabtmp);
    verb_print(cmd->verbose,
            "\nFreed %g MByte for temporal particle (%ld) storage.\n",
            npixel*sizeof(body)/(1024.0*1024.0),npixel);

    if (scanopt(cmd->options, "all"))
        verb_print(cmd->verbose,
            "\n\t%s: selected read points and nbody: %ld %ld\n",
                   routineName, iselect, cmd->nbody);
    else
        verb_print(cmd->verbose,
            "\n\t%s: selected read points = %ld\n", routineName, iselect);

    verb_print(cmd->verbose,
            "%s: average of kappa (%ld particles) = %le\n",
               routineName, cmd->nbody, kavg/((real)cmd->nbody) );

    return SUCCESS;
}
#else // !  THREEDIMCODE

local int inputdata_cfitsio_healpix_map(struct cmdline_data* cmd,
                                        struct  global_data* gd,
                                        string filename,
                                        int ifile, fitsfile *fptr)
{
    string routineName = "inputdata_cfitsio_healpix_map";
    bodyptr p;
    real mass=1;
    real weight=1;

    long ipix;
    double theta;
    double phi;
    double thetamin, thetamax;
    double phimin, phimax;

    float *map;
    long npixel, nside;
    char order1[10];
    char order2[10];
    char coord[10];

    verb_print(cmd->verbose, "\nWorking 2D map...\n");

    npixel = get_fits_size(filename, &nside, order1);
    verb_print(cmd->verbose,
        "filename, ifile, nside, npixel, order1:");
    verb_print(cmd->verbose,
        "%s %d %ld %ld %s\n", filename, ifile, nside, npixel, order1);
    map = read_healpix_map(filename, &nside, coord, order2);
    verb_print(cmd->verbose,
               "\nAllocated %g MByte for temporal map pixel (%ld) storage.\n",
               npixel*sizeof(float)*INMB, npixel);

    if (strcmp(order1,order2)!=0) {printf("Error: Bad ordering\n");}
    verb_print(cmd->verbose,
               "%s: nbody = %d...\n",
               routineName, npixel);
    if (npixel < 1)
        error("%s: npixel = %d is absurd\n", routineName, npixel);

    bodyptr bodytabtmp;
    bodytabtmp = (bodyptr) allocate(npixel * sizeof(body));
    verb_print(cmd->verbose,
               "\nAllocated %g MByte for temporal particle (%ld) storage.\n",
               npixel*sizeof(body)*INMB, npixel);

    //B save optional RA-DEC to a file... 3-columns: RA, DEC, Kappa
    char namebuf[256];
    stream outstr;
    if (scanopt(cmd->options, "save-ra-dec")&&!strnull(cmd->outfile)) {
        sprintf(namebuf, gd->fpfnameOutputFileName);
        outstr = stropen(namebuf, "w!");
    }
    //E

    cmd->nbody=npixel;

    real xmin, ymin;
    real xmax, ymax;
    xmin=0., ymin=0.;
    xmax=0., ymax=0.;

    INTEGER i;
    INTEGER iselect = 0;
    for(ipix=0;ipix<npixel;ipix++) {                // RING loop order
        p = bodytabtmp+ipix;
        Update(p) = FALSE;
        Mask(p) = 1;                                // initialize body's Mask
        pix2ang_ring(nside, ipix, &theta, &phi);
        if (scanopt(cmd->options, "patch")) {
            if (cmd->thetaL < theta && theta < cmd->thetaR) {
                if (cmd->phiL < phi && phi < cmd->phiR) {
                    iselect++;
//                    coordinate_transformation(cmd, gd, theta, phi, Pos(p));
                    Pos(p)[0] = phi;
                    Pos(p)[1] = theta;
                    if (!scanopt(cmd->options, "kappa-constant"))
                        Kappa(p) = map[ipix];
                    else {
                        Kappa(p) = 2.0;
                        if (scanopt(cmd->options, "kappa-constant-one"))
                            Kappa(p) = 1.0;
                    }
                    //B save optional RA-DEC to a file...
                    if (scanopt(cmd->options, "save-ra-dec")
                        &&!strnull(cmd->outfile)) {
                        out_real_mar(outstr, Pos(p)[0]);
                        out_real_mar(outstr, Pos(p)[1]);
                        out_real(outstr, Kappa(p));
                    }
                    //E
                    Type(p) = BODY;
                    Mass(p) = mass;
                    Weight(p) = weight;
                    Id(p) = p-bodytabtmp+iselect;
                    Update(p) = TRUE;
                    xmin = Pos(p)[0];
                    ymin = Pos(p)[1];
                    xmax = Pos(p)[0];
                    ymax = Pos(p)[1];
                }
            }
        } else { // ! all
            iselect++;
//            coordinate_transformation(cmd, gd, theta, phi, Pos(p));
            Pos(p)[0] = phi;
            Pos(p)[1] = theta;
            if (!scanopt(cmd->options, "kappa-constant"))
                Kappa(p) = map[ipix];
            else {
                Kappa(p) = 2.0;
                if (scanopt(cmd->options, "kappa-constant-one"))
                    Kappa(p) = 1.0;
            }
            //B save optional RA-DEC to a file...
            if (scanopt(cmd->options, "save-ra-dec")
                &&!strnull(cmd->outfile)) {
                out_real_mar(outstr, Pos(p)[0]);
                out_real_mar(outstr, Pos(p)[1]);
                out_real(outstr, Kappa(p));
            }
            //E
            Type(p) = BODY;
            Mass(p) = mass;
            Weight(p) = weight;
            Id(p) = p-bodytabtmp+iselect;
            Update(p) = TRUE;
            xmin = Pos(p)[0];
            ymin = Pos(p)[1];
            xmax = Pos(p)[0];
            ymax = Pos(p)[1];
            //B correction 2025-05-03 :: look for edge-effects
            // activate a flag for this catalog, that it is using patch-with-all
            //  and use it in EvalHist routine...
    #if defined(NMultipoles) && defined(NONORMHIST)
            if (scanopt(cmd->options, "patch-with-all")) {
                UpdatePivot(p) = TRUE;
                if (cmd->thetaL < theta && theta < cmd->thetaR) {
                    if (cmd->phiL < phi && phi < cmd->phiR) {
                        UpdatePivot(p) = TRUE;
                        gd->pivotCount += 1;
                    } else {
                        UpdatePivot(p) = FALSE;
                    }
                } else {
                    UpdatePivot(p) = FALSE;
                }
            }
    #endif
            //E
        } // ! all
    } // ! end loop ipix

    //B save optional RA-DEC to a file...
    if (scanopt(cmd->options, "save-ra-dec")&&!strnull(cmd->outfile)
                            &&scanopt(cmd->options, "stop")) {
        fclose(outstr);
        exit(1);            // need to free all mem allocation...
    }
    //E

    //B correction 2025-05-03 :: look for edge-effects
#if defined(NMultipoles) && defined(NONORMHIST)
    if (scanopt(cmd->options, "patch-with-all")) {
        verb_print(cmd->verbose,
                   "\n%s: total number of pixels to be pivots: %ld\n",
                   routineName, gd->pivotCount);
    }
#endif
    //E

    free(map);
    verb_print(cmd->verbose,
               "\nFreed %g MByte for temporal map pixel (%ld) storage.\n",
               npixel*sizeof(float)*INMB, npixel);

    bodyptr q;
    if (scanopt(cmd->options, "patch"))
        cmd->nbody = iselect;

    gd->nbodyTable[ifile] = cmd->nbody;
    bodytable[ifile] = (bodyptr) allocate(cmd->nbody * sizeof(body));

    real kavg = 0;
    INTEGER ij=0;
    for(ipix=0;ipix<npixel;ipix++){
        q = bodytabtmp+ipix;
        if(Update(q)) {
            p = bodytable[ifile]+ij;
            Pos(p)[0] = Pos(q)[0];
            Pos(p)[1] = Pos(q)[1];
            Kappa(p) = Kappa(q);
            Type(p) = Type(q);
            Mass(p) = mass;
            Weight(p) = weight;
            Id(p) = p-bodytable[ifile]+ipix;
            Mask(p) = Mask(q);
            xmin = MIN(xmin,Pos(p)[0]);
            ymin = MIN(ymin,Pos(p)[1]);
            xmax = MAX(xmax,Pos(p)[0]);
            ymax = MAX(ymax,Pos(p)[1]);
            ij++;
            kavg += Kappa(p);
        }
    }

    char file[180] = "outputmap.fits" ;
    char fileforce[180] ;
    if (scanopt(cmd->options, "plot-map-gif")) {
        //B add ifile tag...
        sprintf(fileforce, "!%s/%s",                // leading !
                gd->outputDir, file);               //  to allow overwrite
        //E
        verb_print(cmd->verbose,
                   "\t%s: %s %s...\n",
                   routineName, "\n\t\tsaving map to a fits file:", fileforce);
        map = (float *)malloc(npixel*sizeof(float));
        for(ipix=0;ipix<npixel;ipix++){             // all pixels in
                                                    //  the sphere are filled
            q = bodytabtmp+ipix;
            if(Update(q)) {
                map[ipix] = Kappa(q);
            } else {
                map[ipix] = 0.0;
            }
        }
        //B write healpix map in RING order: "0"
        //      and
        //      "G = Galactic, E = ecliptic, C = celestial = equatorial"
        write_healpix_map(map, nside, fileforce, 0, "C");  // choose COORDSYS...
        fprintf(stdout,"\t\tfile written\n");
        free(map);
    }

    verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                          "\n\t%s: min and max of x = %f %f\n",
                          routineName, xmin, xmax);
    verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                          "\t%s: min and max of y = %f %f\n",
                          routineName, ymin, ymax);

    free(bodytabtmp);
    verb_print(cmd->verbose,
            "\nFreed %g MByte for temporal particle (%ld) storage.\n",
            npixel*sizeof(body)*INMB,npixel);

    if (scanopt(cmd->options, "all"))
        verb_print(cmd->verbose,
            "\n\t%s: selected read points and nbody: %ld %ld\n",
                   routineName, iselect, cmd->nbody);
    else
        verb_print(cmd->verbose,
            "\n\t%s: selected read points = %ld\n", routineName, iselect);

    verb_print(cmd->verbose,
            "%s: average of kappa (%ld particles) = %le\n",
               routineName, cmd->nbody, kavg/((real)cmd->nbody) );

    return SUCCESS;
}

#endif // !  THREEDIMCODE



//#ifdef MASKED
// reading healpix mask map
//  3D only... and full-sky only
local int inputdata_cfitsio_healpix_map_mask(struct cmdline_data* cmd,
                                        struct  global_data* gd,
                                        string filename,
                                        int ifile, fitsfile *fptr)
{
    string routine_name = "inputdata_cfitsio_healpix_map_mask";
    bodyptr p;
    vector q;
    int k;

    long ipix;
    double theta;
    double phi;

    float *map;
    long npixel, nside;
    char order1[10];
    char order2[10];
    char coord[10];

#if THREEDIMCODE
    verb_print(cmd->verbose, "\nWorking 3D map...\n");
#else
    error("\nOnly 3D is implemented so far... exiting...\n\n");
#endif

    npixel = get_fits_size(filename, &nside, order1);
    verb_print(cmd->verbose,
        "filename, ifile, nside, npixel, order1:");
    verb_print(cmd->verbose,
        "%s %d %ld %ld %s\n", filename, ifile, nside, npixel, order1);
    map = read_healpix_map(filename, &nside, coord, order2);
    verb_print(cmd->verbose,
               "\nAllocated %g MByte for temporal map pixel (%ld) storage.\n",
               npixel*sizeof(float)/(1024.0*1024.0),
               npixel);

    if (strcmp(order1,order2)!=0) {printf("Error: Bad ordering\n");}
        verb_print(cmd->verbose, "%s: nbody = %d...\n",
               routine_name, npixel);
    if (npixel < 1)
        error("%s: npixel = %d is absurd\n", routine_name, npixel);

    if (npixel != gd->nbodyTable[gd->iCatalogs[0]])
        error("%s: npixel = %ld is not equal to npixel in cat 0: %ld\n",
              routine_name, npixel, gd->nbodyTable[gd->iCatalogs[0]]);

    INTEGER iselect = 0;
    for(ipix=0;ipix<npixel;ipix++) {                // RING loop order
        p = bodytable[gd->iCatalogs[0]]+ipix;
        pix2ang_ring(nside, ipix, &theta, &phi);
        coordinate_transformation(cmd, gd, theta, phi, q);
        DO_COORD(k) {
            if (Pos(p)[k] != q[k])
                error("%s: mask position is not equal: %g %g\n",
                      routine_name, Pos(p)[k], q[k]);
        }

        Mask(p) = map[ipix];
        if (Mask(p) == 0) {
            iselect++;
        }
    } // ! end loop ipix

    verb_print(cmd->verbose,
        "\n\t%s: masked pixels = %ld\n",
               routine_name, iselect);
    verb_print(cmd->verbose,
        "\t%s: unmasked pixels = %ld\n",
               routine_name, gd->nbodyTable[gd->iCatalogs[0]]-iselect);
    
    return SUCCESS;
}
//#endif



// reading healpix map
//  3D only
local int inputdata_cfitsio_healpix_map_mask_inside(struct cmdline_data* cmd,
                                        struct  global_data* gd,
                                        string filename,
                                        int ifile, fitsfile *fptr)
{
    string routineName = "inputdata_cfitsio_healpix_map_mask_inside";
    bodyptr p;
    real mass=1;
    real weight=1;

    long ipix;
    double theta;
    double phi;
    double thetamin, thetamax;
    double phimin, phimax;

    float *map;
    long npixel, nside;
    char order1[10];
    char order2[10];
    char coord[10];

#if THREEDIMCODE
    verb_print(cmd->verbose, "\nWorking 3D map...\n");
#else
    error("\nOnly 3D is implemented so far... exiting...\n\n");
#endif

    npixel = get_fits_size(filename, &nside, order1);
    verb_print(cmd->verbose,
        "filename, ifile, nside, npixel, order1:");
    verb_print(cmd->verbose,
        "%s %d %ld %ld %s\n", filename, ifile, nside, npixel, order1);
    map = read_healpix_map(filename, &nside, coord, order2);
    verb_print(cmd->verbose,
               "\nAllocated %g MByte for temporal map pixel (%ld) storage.\n",
               npixel*sizeof(float)*INMB,
               npixel);

    if (strcmp(order1,order2)!=0) {printf("Error: Bad ordering\n");}
    verb_print(cmd->verbose,
               "%s: nbody = %d...\n",
               routineName, npixel);
    if (npixel < 1)
        error("%s: npixel = %d is absurd\n", routineName, npixel);

    bodyptr bodytabtmp;
    bodytabtmp = (bodyptr) allocate(npixel * sizeof(body));
    verb_print(cmd->verbose,
               "\nAllocated %g MByte for temporal particle (%ld) storage.\n",
               npixel*sizeof(body)*INMB,
               npixel);

    cmd->nbody=npixel;

    real xmin, ymin, zmin;
    real xmax, ymax, zmax;
    xmin=0., ymin=0., zmin=0.;
    xmax=0., ymax=0., zmax=0.;

    INTEGER i;
    INTEGER iselect = 0;
    for(ipix=0;ipix<npixel;ipix++) {                // RING loop order
        p = bodytabtmp+ipix;
        Update(p) = FALSE;
        Mask(p) = 1;                                // initialize body's Mask
        pix2ang_ring(nside, ipix, &theta, &phi);
        if (scanopt(cmd->options, "patch")) {
            if (cmd->thetaL < theta && theta < cmd->thetaR) {
                if (cmd->phiL < phi && phi < cmd->phiR) {
                    iselect++;
                    coordinate_transformation(cmd, gd, theta, phi, Pos(p));
                    if (!scanopt(cmd->options, "kappa-constant"))
                        Kappa(p) = map[ipix];
                    else {
                        Kappa(p) = 2.0;
                        if (scanopt(cmd->options, "kappa-constant-one"))
                            Kappa(p) = 1.0;
                    }
                    Type(p) = BODY;
                    Mass(p) = mass;
                    Weight(p) = weight;
                    Id(p) = p-bodytabtmp+iselect;
                    Update(p) = TRUE;
                    xmin = Pos(p)[0];
                    ymin = Pos(p)[1];
                    zmin = Pos(p)[2];
                    xmax = Pos(p)[0];
                    ymax = Pos(p)[1];
                    zmax = Pos(p)[2];
                }
            }
        } else { // ! all
            coordinate_transformation(cmd, gd, theta, phi, Pos(p));
            Kappa(p) = map[ipix];
            if (Kappa(p)!=0.) {                     // set unmasked pixel
                if (scanopt(cmd->options, "kappa-constant")) {
                    Kappa(p) = 2.0;
                    if (scanopt(cmd->options, "kappa-constant-one"))
                        Kappa(p) = 1.0;
                }
                Mask(p) = 1;
                iselect++;
                Type(p) = BODY;
                Mass(p) = mass;
                Weight(p) = weight;
                Id(p) = p-bodytabtmp+iselect;
                Update(p) = TRUE;
                xmin = Pos(p)[0];
                ymin = Pos(p)[1];
                zmin = Pos(p)[2];
                xmax = Pos(p)[0];
                ymax = Pos(p)[1];
                zmax = Pos(p)[2];

                //B correction 2025-05-03 :: look for edge-effects
                // activate a flag for this catalog, that it is using patch-with-all
                //  and use it in EvalHist routine...
#if defined(NMultipoles) && defined(NONORMHIST)
                if (scanopt(cmd->options, "patch-with-all")) {
                    UpdatePivot(p) = TRUE;
                    if (cmd->thetaL < theta && theta < cmd->thetaR) {
                        if (cmd->phiL < phi && phi < cmd->phiR) {
                            UpdatePivot(p) = TRUE;
                            gd->pivotCount += 1;
                        } else {
                            UpdatePivot(p) = FALSE;
                        }
                    } else {
                        UpdatePivot(p) = FALSE;
                    }
                }
#endif
                //E
            }
        } // ! all
    } // ! end loop ipix

    verb_print(cmd->verbose,
        "\n\t%s: masked pixels = %ld\n",
               routineName, cmd->nbody-iselect);
    verb_print(cmd->verbose,
        "\t%s: unmasked pixels = %ld\n",
               routineName, iselect);

    //B correction 2025-05-03 :: look for edge-effects
#if defined(NMultipoles) && defined(NONORMHIST)
    if (scanopt(cmd->options, "patch-with-all")) {
        verb_print(cmd->verbose,
                   "\n%s: total number of pixels to be pivots: %ld\n",
                   routineName, gd->pivotCount);
    }
#endif
    //E

    free(map);
    verb_print(cmd->verbose,
               "\nFreed %g MByte for temporal map pixel (%ld) storage.\n",
               npixel*sizeof(float)/(1024.0*1024.0),
               npixel);


    bodyptr q;
    cmd->nbody = iselect;
    gd->nbodyTable[ifile] = cmd->nbody;
    bodytable[ifile] = (bodyptr) allocate(cmd->nbody * sizeof(body));

    real kavg = 0;
    INTEGER ij=0;
    for(ipix=0;ipix<npixel;ipix++){
        q = bodytabtmp+ipix;
        if(Update(q)) {
            p = bodytable[ifile]+ij;
            Pos(p)[0] = Pos(q)[0];
            Pos(p)[1] = Pos(q)[1];
            Pos(p)[2] = Pos(q)[2];
            Kappa(p) = Kappa(q);
            Type(p) = Type(q);
            Mass(p) = mass;
            Weight(p) = weight;
            Id(p) = p-bodytable[ifile]+ipix;
            Mask(p) = Mask(q);
            xmin = MIN(xmin,Pos(p)[0]);
            ymin = MIN(ymin,Pos(p)[1]);
            zmin = MIN(zmin,Pos(p)[2]);
            xmax = MAX(xmax,Pos(p)[0]);
            ymax = MAX(ymax,Pos(p)[1]);
            zmax = MAX(zmax,Pos(p)[2]);
            ij++;
            kavg += Kappa(p);
        }
    }

    char file[180] = "outputmap.fits" ;
    char fileforce[180] ;
    if (scanopt(cmd->options, "plot-map-gif")) {
        sprintf(fileforce, "!%s/%s",                // leading !
                gd->outputDir, file);               //  to allow overwrite
        verb_print(cmd->verbose,
                   "\t%s: %s %s...\n",
                   routineName, "\n\t\tsaving map to a fits file:", fileforce);
        map = (float *)malloc(npixel*sizeof(float));
        for(ipix=0;ipix<npixel;ipix++){
            q = bodytabtmp+ipix;
            if(Update(q)) {
                map[ipix] = Kappa(q);
            } else {
                map[ipix] = 0.0;
            }
        }
        //B write healpix map in RING order "0"
        //      and
        //      "G = Galactic, E = ecliptic, C = celestial = equatorial"
        write_healpix_map(map, nside, fileforce, 0, "C");
        fprintf(stdout,"\t\tfile written\n");
        free(map);
    }

    verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                          "\n\t%s: min and max of x = %f %f\n",
                          routineName, xmin, xmax);
    verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                          "\t%s: min and max of y = %f %f\n",
                          routineName, ymin, ymax);
    verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                          "\t%s: min and max of z = %f %f\n",
                          routineName, zmin, zmax);

    free(bodytabtmp);
    verb_print(cmd->verbose,
            "\nFreed %g MByte for temporal particle (%ld) storage.\n",
            npixel*sizeof(body)/(1024.0*1024.0),npixel);

    if (scanopt(cmd->options, "all"))
        verb_print(cmd->verbose,
            "\n\t%s: selected read points and nbody: %ld %ld\n",
                   routineName, iselect, cmd->nbody);
    else
        verb_print(cmd->verbose,
            "\n\t%s: selected read points = %ld\n", routineName, iselect);

    verb_print(cmd->verbose,
            "%s: average of kappa (%ld particles) = %le\n",
               routineName, cmd->nbody, kavg/((real)cmd->nbody) );

    return SUCCESS;
}

// Routine to read numpy-healpix files
//  so far only RING scheme
local int inputdata_numpy_healpix(struct cmdline_data* cmd,
                                    struct  global_data* gd,
                                    string filename, int ifile)
{
    string routineName = "inputdata_numpy_healpix";
    char card[FLEN_CARD];
    int status = 0, nkeys, ii;                      // MUST initialize status
    int ncols;
    char buf[200];
    string headerScript = "python header_script.py";

    gd->input_comment = "numpy-healpix input file";

    verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                          "\t%s: opening numpy-healpix file: %s...\n",
                          routineName,filename);
    if (scanopt(cmd->options, "header-info")){
        verb_print(cmd->verbose,"\nHeader information:");
        sprintf(buf,"%s",headerScript);
        verb_print(cmd->verbose,
                   "\n\t%s: header processing: executing %s...\n",
                   routineName, headerScript);
        system(buf);
        verb_print(cmd->verbose, "\tdone.\n");
        if (scanopt(cmd->options, "stop-numpy")) {
            if (strnull(cmd->outfile)) {
                exit(1);
            }
        }
    }

    if (scanopt(cmd->options, "read-mask")) {
        if (ifile == 0) {
            inputdata_numpy_healpix_map(cmd, gd, filename, ifile);
        } else {
            if (ifile != 1)
                error("\t%s: read-mask ifile = %d is absurd\n",
                      routineName, ifile+1);
            inputdata_numpy_healpix_map_mask(cmd, gd, filename, ifile);
        }
    } else {
        if (scanopt(cmd->options, "mask-inside")) {
            inputdata_numpy_healpix_map_mask_inside(cmd, gd,
                                                      filename, ifile);
        } else {
            inputdata_numpy_healpix_map(cmd, gd, filename, ifile);
        }
    }

    return SUCCESS;
}

// reading numpy-healpix map
//  3D only...
local int inputdata_numpy_healpix_map(struct cmdline_data* cmd,
                                        struct  global_data* gd,
                                        string filename,
                                        int ifile)
{
    string routineName = "inputdata_numpy_healpix_map";
    bodyptr p;
    real mass=1;
    real weight=1;

    long j;
    long ipix;
    double theta;
    double phi;
    double thetamin, thetamax;
    double phimin, phimax;

//    float *map;
    double *map;
    long npixel, nside;
    string order = "RING";
    FILE *fp;

#if THREEDIMCODE
    verb_print(cmd->verbose, "\nWorking 3D map...\n");
#else
    error("\nOnly 3D is implemented so far... exiting...\n\n");
#endif

    // info runing header_script.py and set it in cmd->nbody
    npixel = cmd->nbody;            // must be in the form 12*nside^2
    nside = npix2nside(npixel);
    verb_print(cmd->verbose,
        "filename, ifile, nside, npixel, order:");
    verb_print(cmd->verbose,
        "%s %d %ld %ld %s\n", filename, ifile, nside, npixel, order);

    fp = stropen(filename, "rb");       // open numpy file
//    map=(float *)malloc(sizeof(float)*npixel);
    map=(double *)malloc(sizeof(double)*npixel);
//    verb_print(cmd->verbose,
//               "\nAllocated %g MByte for temporal map pixel (%ld) storage.\n",
//               npixel*sizeof(float)/(1024.0*1024.0),
//               npixel);
    verb_print(cmd->verbose,
               "\nAllocated %g MByte for temporal map pixel (%ld) storage.\n",
               npixel*sizeof(double)*INMB,
               npixel);
    for(j=0;j<npixel;j++){
//        fread(&map[j], sizeof(float), 1, fp);
        fread(&map[j], sizeof(double), 1, fp);
    }
    fclose(fp);             // close numpy file.

    verb_print(cmd->verbose,
               "%s: nbody = %ld...\n",
               routineName, npixel);
    if (npixel < 1)
        error("%s: npixel = %ld is absurd\n", routineName, npixel);

    bodyptr bodytabtmp;
    bodytabtmp = (bodyptr) allocate(npixel * sizeof(body));
    verb_print(cmd->verbose,
               "\nAllocated %g MByte for temporal particle (%ld) storage.\n",
               npixel*sizeof(body)*INMB,
               npixel);

//    cmd->nbody=npixel;

    real xmin, ymin, zmin;
    real xmax, ymax, zmax;
    xmin=0., ymin=0., zmin=0.;
    xmax=0., ymax=0., zmax=0.;

    INTEGER i;
    INTEGER iselect = 0;
    for(ipix=0;ipix<npixel;ipix++) {                // RING loop order
        p = bodytabtmp+ipix;
        Update(p) = FALSE;
        Mask(p) = 1;                                // initialize body's Mask
        pix2ang_ring(nside, ipix, &theta, &phi);
        if (scanopt(cmd->options, "patch")) {
            if (cmd->thetaL < theta && theta < cmd->thetaR) {
                if (cmd->phiL < phi && phi < cmd->phiR) {
                    iselect++;
                    coordinate_transformation(cmd, gd, theta, phi, Pos(p));
                    if (!scanopt(cmd->options, "kappa-constant"))
                        Kappa(p) = map[ipix];
                    else {
                        Kappa(p) = 2.0;
                        if (scanopt(cmd->options, "kappa-constant-one"))
                            Kappa(p) = 1.0;
                    }
                    Type(p) = BODY;
                    Mass(p) = mass;
                    Weight(p) = weight;
                    Id(p) = p-bodytabtmp+iselect;
                    Update(p) = TRUE;
                    xmin = Pos(p)[0];
                    ymin = Pos(p)[1];
                    zmin = Pos(p)[2];
                    xmax = Pos(p)[0];
                    ymax = Pos(p)[1];
                    zmax = Pos(p)[2];
                }
            }
        } else { // ! all
            iselect++;
            coordinate_transformation(cmd, gd, theta, phi, Pos(p));
            if (!scanopt(cmd->options, "kappa-constant"))
                Kappa(p) = map[ipix];
            else {
                Kappa(p) = 2.0;
                if (scanopt(cmd->options, "kappa-constant-one"))
                    Kappa(p) = 1.0;
            }
            Type(p) = BODY;
            Mass(p) = mass;
            Weight(p) = weight;
            Id(p) = p-bodytabtmp+iselect;
            Update(p) = TRUE;
            xmin = Pos(p)[0];
            ymin = Pos(p)[1];
            zmin = Pos(p)[2];
            xmax = Pos(p)[0];
            ymax = Pos(p)[1];
            zmax = Pos(p)[2];
            
            //B correction 2025-05-03 :: look for edge-effects
            // activate a flag for this catalog, that it is using patch-with-all
            //  and use it in EvalHist routine...
    #if defined(NMultipoles) && defined(NONORMHIST)
            if (scanopt(cmd->options, "patch-with-all")) {
                UpdatePivot(p) = TRUE;
                if (cmd->thetaL < theta && theta < cmd->thetaR) {
                    if (cmd->phiL < phi && phi < cmd->phiR) {
                        UpdatePivot(p) = TRUE;
                        gd->pivotCount += 1;
                    } else {
                        UpdatePivot(p) = FALSE;
                    }
                } else {
                    UpdatePivot(p) = FALSE;
                }
            }
    #endif
            //E

        } // ! all
    } // ! end loop ipix

    //B correction 2025-05-03 :: look for edge-effects
#if defined(NMultipoles) && defined(NONORMHIST)
    if (scanopt(cmd->options, "patch-with-all")) {
        verb_print(cmd->verbose,
                   "\n%s: total number of pixels to be pivots: %ld\n",
                   routineName, gd->pivotCount);
    }
#endif
    //E

    free(map);
    verb_print(cmd->verbose,
               "\nFreed %g MByte for temporal map pixel (%ld) storage.\n",
               npixel*sizeof(float)/(1024.0*1024.0),
               npixel);


    bodyptr q;
    cmd->nbody = iselect;

    gd->nbodyTable[ifile] = cmd->nbody;
    bodytable[ifile] = (bodyptr) allocate(cmd->nbody * sizeof(body));

    real kavg = 0;
    INTEGER ij=0;
    for(ipix=0;ipix<npixel;ipix++){
        q = bodytabtmp+ipix;
        if(Update(q)) {
            p = bodytable[ifile]+ij;
            Pos(p)[0] = Pos(q)[0];
            Pos(p)[1] = Pos(q)[1];
            Pos(p)[2] = Pos(q)[2];
            Kappa(p) = Kappa(q);
            Type(p) = Type(q);
            Mass(p) = mass;
            Weight(p) = weight;
            Id(p) = p-bodytable[ifile]+ipix;
            Mask(p) = Mask(q);
            xmin = MIN(xmin,Pos(p)[0]);
            ymin = MIN(ymin,Pos(p)[1]);
            zmin = MIN(zmin,Pos(p)[2]);
            xmax = MAX(xmax,Pos(p)[0]);
            ymax = MAX(ymax,Pos(p)[1]);
            zmax = MAX(zmax,Pos(p)[2]);
            ij++;
            kavg += Kappa(p);
        }
    }

    float *mapout;
    char file[180] = "outputmap.fits" ;
    char fileforce[180] ;
    if (scanopt(cmd->options, "plot-map-gif")) {
        sprintf(fileforce, "!%s/%s",                // leading !
                gd->outputDir, file);               //  to allow overwrite
        verb_print(cmd->verbose,
                   "\t%s: %s %s...\n",
                   routineName, "\n\t\tsaving map to a fits file:", fileforce);
        mapout = (float *)malloc(npixel*sizeof(float));
        for(ipix=0;ipix<npixel;ipix++){
            q = bodytabtmp+ipix;
            if(Update(q)) {
                mapout[ipix] = Kappa(q);
            } else {
                mapout[ipix] = 0.0;
            }
        }
        //B write healpix map in RING order "0"
        //      and
        //      "G = Galactic, E = ecliptic, C = celestial = equatorial"
        write_healpix_map(mapout, nside, fileforce, 0, "C");
        fprintf(stdout,"\t\tfile written\n");
        free(map);
    }

    verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                          "\n\t%s: min and max of x = %f %f\n",
                          routineName, xmin, xmax);
    verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                          "\t%s: min and max of y = %f %f\n",
                          routineName, ymin, ymax);
    verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                          "\t%s: min and max of z = %f %f\n",
                          routineName, zmin, zmax);

    free(bodytabtmp);
    verb_print(cmd->verbose,
            "\nFreed %g MByte for temporal particle (%ld) storage.\n",
            npixel*sizeof(body)*INMB,npixel);

    if (scanopt(cmd->options, "all"))
        verb_print(cmd->verbose,
            "\n\t%s: selected read points and nbody: %ld %ld\n",
                   routineName, iselect, cmd->nbody);
    else
        verb_print(cmd->verbose,
            "\n\t%s: selected read points = %ld\n", routineName, iselect);

    verb_print(cmd->verbose,
            "%s: average of kappa (%ld particles) = %le\n",
               routineName, cmd->nbody, kavg/((real)cmd->nbody) );

    return SUCCESS;
}

// reading numpy-healpix mask map
//  3D only... and full-sky only
//          working but as a fits-healpix map,
//              then, the mask file must be in healpix fmt
local int inputdata_numpy_healpix_map_mask(struct cmdline_data* cmd,
                                        struct  global_data* gd,
                                        string filename,
                                        int ifile)
{
    string routine_name = "inputdata_cfitsio_healpix_map_mask";
    bodyptr p;
    vector q;
    int k;

    long ipix;
    double theta;
    double phi;

    float *map;
    long npixel, nside;
    char order1[10];
    char order2[10];
    char coord[10];

#if THREEDIMCODE
    verb_print(cmd->verbose, "\nWorking 3D map...\n");
#else
    error("\nOnly 3D is implemented so far... exiting...\n\n");
#endif

    npixel = get_fits_size(filename, &nside, order1);
    verb_print(cmd->verbose,
        "filename, ifile, nside, npixel, order1:");
    verb_print(cmd->verbose,
        "%s %d %ld %ld %s\n", filename, ifile, nside, npixel, order1);
    map = read_healpix_map(filename, &nside, coord, order2);
    verb_print(cmd->verbose,
               "\nAllocated %g MByte for temporal map pixel (%ld) storage.\n",
               npixel*sizeof(float)*INMB,
               npixel);

    if (strcmp(order1,order2)!=0) {printf("Error: Bad ordering\n");}
        verb_print(cmd->verbose, "%s: nbody = %d...\n",
               routine_name, npixel);
    if (npixel < 1)
        error("%s: npixel = %d is absurd\n", routine_name, npixel);

    if (npixel != gd->nbodyTable[gd->iCatalogs[0]])
        error("%s: npixel = %ld is not equal to npixel in cat 0: %ld\n",
              routine_name, npixel, gd->nbodyTable[gd->iCatalogs[0]]);

    INTEGER iselect = 0;
    for(ipix=0;ipix<npixel;ipix++) {                // RING loop order
        p = bodytable[gd->iCatalogs[0]]+ipix;
        pix2ang_ring(nside, ipix, &theta, &phi);
        coordinate_transformation(cmd, gd, theta, phi, q);
        DO_COORD(k) {
            if (Pos(p)[k] != q[k])
                error("%s: mask position is not equal: %g %g\n",
                      routine_name, Pos(p)[k], q[k]);
        }

        Mask(p) = map[ipix];
        if (Mask(p) == 0) {
            iselect++;
        }
    } // ! end loop ipix

    verb_print(cmd->verbose,
        "\n\t%s: masked pixels = %ld\n",
               routine_name, iselect);
    verb_print(cmd->verbose,
        "\t%s: unmasked pixels = %ld\n",
               routine_name, gd->nbodyTable[gd->iCatalogs[0]]-iselect);
    
    return SUCCESS;
}

// reading numpy-healpix map
//  3D only
//      ... not working...
local int inputdata_numpy_healpix_map_mask_inside(struct cmdline_data* cmd,
                                        struct  global_data* gd,
                                        string filename,
                                        int ifile)
{
    string routineName = "inputdata_cfitsio_healpix_map_mask_inside";
    bodyptr p;
    real mass=1;
    real weight=1;

    long ipix;
    double theta;
    double phi;
    double thetamin, thetamax;
    double phimin, phimax;

    float *map;
    long npixel, nside;
    char order1[10];
    char order2[10];
    char coord[10];

#if THREEDIMCODE
    verb_print(cmd->verbose, "\nWorking 3D map...\n");
#else
    error("\nOnly 3D is implemented so far... exiting...\n\n");
#endif

    npixel = get_fits_size(filename, &nside, order1);
    verb_print(cmd->verbose,
        "filename, ifile, nside, npixel, order1:");
    verb_print(cmd->verbose,
        "%s %d %ld %ld %s\n", filename, ifile, nside, npixel, order1);
    map = read_healpix_map(filename, &nside, coord, order2);
    verb_print(cmd->verbose,
               "\nAllocated %g MByte for temporal map pixel (%ld) storage.\n",
               npixel*sizeof(float)*INMB,
               npixel);

    if (strcmp(order1,order2)!=0) {printf("Error: Bad ordering\n");}
    verb_print(cmd->verbose,
               "%s: nbody = %d...\n",
               routineName, npixel);
    if (npixel < 1)
        error("%s: npixel = %d is absurd\n", routineName, npixel);

    bodyptr bodytabtmp;
    bodytabtmp = (bodyptr) allocate(npixel * sizeof(body));
    verb_print(cmd->verbose,
               "\nAllocated %g MByte for temporal particle (%ld) storage.\n",
               npixel*sizeof(body)*INMB,
               npixel);

    cmd->nbody=npixel;

    real xmin, ymin, zmin;
    real xmax, ymax, zmax;
    xmin=0., ymin=0., zmin=0.;
    xmax=0., ymax=0., zmax=0.;

    INTEGER i;
    INTEGER iselect = 0;
    for(ipix=0;ipix<npixel;ipix++) {                // RING loop order
        p = bodytabtmp+ipix;
        Update(p) = FALSE;
        pix2ang_ring(nside, ipix, &theta, &phi);
        if (scanopt(cmd->options, "patch")) {
            if (cmd->thetaL < theta && theta < cmd->thetaR) {
                if (cmd->phiL < phi && phi < cmd->phiR) {
                    iselect++;
                    coordinate_transformation(cmd, gd, theta, phi, Pos(p));
                    if (!scanopt(cmd->options, "kappa-constant"))
                        Kappa(p) = map[ipix];
                    else {
                        Kappa(p) = 2.0;
                        if (scanopt(cmd->options, "kappa-constant-one"))
                            Kappa(p) = 1.0;
                    }
                    Type(p) = BODY;
                    Mass(p) = mass;
                    Weight(p) = weight;
                    Id(p) = p-bodytabtmp+iselect;
                    Update(p) = TRUE;
                    xmin = Pos(p)[0];
                    ymin = Pos(p)[1];
                    zmin = Pos(p)[2];
                    xmax = Pos(p)[0];
                    ymax = Pos(p)[1];
                    zmax = Pos(p)[2];
                }
            }
            //E
        } else { // ! all
            coordinate_transformation(cmd, gd, theta, phi, Pos(p));
            Kappa(p) = map[ipix];
            if (Kappa(p)!=0.) {                      // set unmasked pixel
                Mask(p) = 1;
                iselect++;
                Type(p) = BODY;
                Mass(p) = mass;
                Weight(p) = weight;
                Id(p) = p-bodytabtmp+iselect;
                Update(p) = TRUE;
                xmin = Pos(p)[0];
                ymin = Pos(p)[1];
                zmin = Pos(p)[2];
                xmax = Pos(p)[0];
                ymax = Pos(p)[1];
                zmax = Pos(p)[2];
                //E
                
                //B correction 2025-05-03 :: look for edge-effects
                // activate a flag for this catalog, that it is using patch-with-all
                //  and use it in EvalHist routine...
#if defined(NMultipoles) && defined(NONORMHIST)
                if (scanopt(cmd->options, "patch-with-all")) {
                    UpdatePivot(p) = TRUE;
                    if (cmd->thetaL < theta && theta < cmd->thetaR) {
                        if (cmd->phiL < phi && phi < cmd->phiR) {
                            UpdatePivot(p) = TRUE;
                            gd->pivotCount += 1;
                        } else {
                            UpdatePivot(p) = FALSE;
                        }
                    } else {
                        UpdatePivot(p) = FALSE;
                    }
                }
#endif
                //E

                if (scanopt(cmd->options, "kappa-constant")) {
                    Kappa(p) = 2.0;
                    if (scanopt(cmd->options, "kappa-constant-one"))
                        Kappa(p) = 1.0;
                }
            }
        } // ! all
    } // ! end loop ipix

    verb_print(cmd->verbose,
        "\n\t%s: masked pixels = %ld\n",
               routineName, cmd->nbody-iselect);
    verb_print(cmd->verbose,
        "\t%s: unmasked pixels = %ld\n",
               routineName, iselect);

    //B correction 2025-05-03 :: look for edge-effects
#if defined(NMultipoles) && defined(NONORMHIST)
    if (scanopt(cmd->options, "patch-with-all")) {
        verb_print(cmd->verbose,
                   "\n%s: total number of pixels to be pivots: %ld\n",
                   routineName, gd->pivotCount);
    }
#endif
    //E

    free(map);
    verb_print(cmd->verbose,
               "\nFreed %g MByte for temporal map pixel (%ld) storage.\n",
               npixel*sizeof(float)*INMB,
               npixel);


    bodyptr q;
    cmd->nbody = iselect;

    gd->nbodyTable[ifile] = cmd->nbody;
    bodytable[ifile] = (bodyptr) allocate(cmd->nbody * sizeof(body));

    real kavg = 0;
    INTEGER ij=0;
    for(ipix=0;ipix<npixel;ipix++){
        q = bodytabtmp+ipix;
        if(Update(q)) {
            p = bodytable[ifile]+ij;
            Pos(p)[0] = Pos(q)[0];
            Pos(p)[1] = Pos(q)[1];
            Pos(p)[2] = Pos(q)[2];
            Kappa(p) = Kappa(q);
            Type(p) = Type(q);
            Mass(p) = mass;
            Weight(p) = weight;
            Id(p) = p-bodytable[ifile]+ipix;
            Mask(p) = Mask(q);
            xmin = MIN(xmin,Pos(p)[0]);
            ymin = MIN(ymin,Pos(p)[1]);
            zmin = MIN(zmin,Pos(p)[2]);
            xmax = MAX(xmax,Pos(p)[0]);
            ymax = MAX(ymax,Pos(p)[1]);
            zmax = MAX(zmax,Pos(p)[2]);
            ij++;
            kavg += Kappa(p);
        }
    }

    char file[180] = "outputmap.fits" ;
    char fileforce[180] ;
    if (scanopt(cmd->options, "plot-map-gif")) {
        sprintf(fileforce, "!%s/%s",                // leading !
                gd->outputDir, file);               //  to allow overwrite
        verb_print(cmd->verbose,
                   "\t%s: %s %s...\n",
                   routineName, "\n\t\tsaving map to a fits file:", fileforce);
        map = (float *)malloc(npixel*sizeof(float));
        for(ipix=0;ipix<npixel;ipix++){
            q = bodytabtmp+ipix;
            if(Update(q)) {
                map[ipix] = Kappa(q);
            } else {
                map[ipix] = 0.0;
            }
        }
        //B write healpix map in RING order "0"
        //      and
        //      "G = Galactic, E = ecliptic, C = celestial = equatorial"
        write_healpix_map(map, nside, fileforce, 0, "C");
        fprintf(stdout,"\t\tfile written\n");
        free(map);
    }

    verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                          "\n\t%s: min and max of x = %f %f\n",
                          routineName, xmin, xmax);
    verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                          "\t%s: min and max of y = %f %f\n",
                          routineName, ymin, ymax);
    verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                          "\t%s: min and max of z = %f %f\n",
                          routineName, zmin, zmax);

    free(bodytabtmp);
    verb_print(cmd->verbose,
            "\nFreed %g MByte for temporal particle (%ld) storage.\n",
            npixel*sizeof(body)*INMB,npixel);

    if (scanopt(cmd->options, "all"))
        verb_print(cmd->verbose,
            "\n\t%s: selected read points and nbody: %ld %ld\n",
                   routineName, iselect, cmd->nbody);
    else
        verb_print(cmd->verbose,
            "\n\t%s: selected read points = %ld\n", routineName, iselect);

    verb_print(cmd->verbose,
            "%s: average of kappa (%ld particles) = %le\n",
               routineName, cmd->nbody, kavg/((real)cmd->nbody) );

    return SUCCESS;
}

#endif	// ! _cballsio_cfitsio_02_h
