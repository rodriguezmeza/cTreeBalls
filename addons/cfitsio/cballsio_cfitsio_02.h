
//=============================================================================
//        1          2          3          4        ^ 5          6          7

#ifndef _cballsio_cfitsio_02_h
#define _cballsio_cfitsio_02_h

#include "fitsio.h"

local int inputdata_cfitsio_ra_dec(struct cmdline_data* cmd, struct  global_data* gd,
                            string filename, int ifile, fitsfile *fptr);
local int inputdata_cfitsio_xyz(struct cmdline_data* cmd, struct  global_data* gd,
                            string filename, int ifile, fitsfile *fptr);

local int inputdata_cfitsio(struct cmdline_data* cmd, struct  global_data* gd,
                               string filename, int ifile)
{
    fitsfile *fptr;
    char card[FLEN_CARD];
    int status = 0, nkeys, ii;                      // MUST initialize status

    gd->input_comment = "fits input file";

    //B Was CFITSIO compiled with the -D_REENTRANT flag?  1 = yes, 0 = no.
    verb_print(cmd->verbose,
               "\tinputdata_cfitsio: -D_REENTRANT flag: %d...\n\n",
               fits_is_reentrant());
    //E

    if (scanopt(cmd->options, "fits-type-file"))
        fits_open_file(&fptr, filename, READONLY, &status);
    else
        fits_open_data(&fptr, filename, READONLY, &status);
    if (status) {                               // print any error messages
        verb_print(cmd->verbose,
                   "\tinputdata_cfitsio: open status: %d...\n\n",
                   status);
        fits_report_error(stderr, status);
    }

    fits_get_hdrspace(fptr, &nkeys, NULL, &status);
    if (status) {                               // print any error messages
        verb_print(cmd->verbose,
                   "\tinputdata_cfitsio: get_hdrspace status: %d...\n\n",
                   status);
        fits_report_error(stderr, status);
    }
    fits_get_num_rows(fptr, &cmd->nbody, &status);
    if (status) {                               // print any error messages
        verb_print(cmd->verbose,
                   "\tinputdata_cfitsio: get_num_rows status: %d...\n\n",
                   status);
        fits_report_error(stderr, status);
    }

    if (scanopt(cmd->options, "fits-header-info")){
        for (ii = 1; ii <= nkeys; ii++) {
            fits_read_record(fptr, ii, card, &status); // read keyword
            printf("%s\n", card);
        }
        printf("END\n\n");                          // terminate listing with END */
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
    } else { // ! fits-header-info
        for (ii = 1; ii <= nkeys; ii++) {
            fits_read_record(fptr, ii, card, &status); // read keyword
        }
        fits_get_num_rows(fptr, &cmd->nbody, &status);
        int ncols;
        fits_get_num_cols(fptr, &ncols, &status);
        verb_print(cmd->verbose,
                   "\tinputdata_cfitsio: nbody = %d... and number of columns = %d\n\n",
                   cmd->nbody, ncols);
        if (cmd->nbody < 1)
            error("tinputdata_cfitsio:: nbody = %d is absurd\n", cmd->nbody);
    }

    if (scanopt(cmd->options, "stop")) {
        fits_close_file(fptr, &status);
        if (status) {                               // print any error messages
            verb_print(cmd->verbose,
                       "\tinputdata_cfitsio: status: %d...\n\n",
                       status);
            fits_report_error(stderr, status);
        }
        exit(1);
    }

    if (scanopt(cmd->options, "ra-dec")) {
        inputdata_cfitsio_ra_dec(cmd, gd,
                                 filename, ifile, fptr);
    } else { // ! ra-dec
        inputdata_cfitsio_xyz(cmd, gd,
                                 filename, ifile, fptr);
    }

    return SUCCESS;
}

local int inputdata_cfitsio_ra_dec(struct cmdline_data* cmd, struct  global_data* gd,
                               string filename, int ifile, fitsfile *fptr)
{
    bodyptr p;
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

    datatype = 42;
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
        if (status) {                                   // print any error messages
            verb_print(cmd->verbose,
                       "\tinputdata_cfitsio: status: %d...\n\n",
                       status);
            fits_report_error(stderr, status);
        }
        colnum = 3;
        fits_read_col(fptr, datatype, colnum, firstrow, firstelem,
                      nelements, &nulval, arrayDEC, &anynul, &status);
        if (status) {                                   // print any error messages
            verb_print(cmd->verbose,
                       "\tinputdata_cfitsio: status: %d...\n\n",
                       status);
            fits_report_error(stderr, status);
        }
        if (scanopt(cmd->options, "arfken")) {
            DO_BODY(p, bodytable[ifile], bodytable[ifile]+cmd->nbody) {
                ra = arrayRA[p-bodytable[ifile]];
                dec = arrayDEC[p-bodytable[ifile]];
                Pos(p)[0] = rsin(dec)*rcos(ra);
                Pos(p)[1] = rsin(dec)*rsin(ra);
                Pos(p)[2] = rcos(dec);
            }
        } else {
            DO_BODY(p, bodytable[ifile], bodytable[ifile]+cmd->nbody) {
                ra = arrayRA[p-bodytable[ifile]];
                dec = arrayDEC[p-bodytable[ifile]];
                Pos(p)[0] = rcos(dec)*rcos(ra);
                Pos(p)[1] = rcos(dec)*rsin(ra);
                Pos(p)[2] = rsin(dec);
            }
        }

    free(arrayRA);
    free(arrayDEC);

    DO_BODY(p, bodytable[ifile], bodytable[ifile]+cmd->nbody) {
        Type(p) = BODY;
        Weight(p) = weight;
        Id(p) = p-bodytable[ifile]+1;
    }

    gd->nbodyTable[ifile] = cmd->nbody;

    real kavg=0.0;
    DO_BODY(p, bodytable[ifile], bodytable[ifile]+cmd->nbody) {
        Type(p) = BODY;
        Weight(p) = weight;
        Id(p) = p-bodytable[ifile]+1;
        kavg += Kappa(p);
    }
    verb_print(cmd->verbose,
               "inputdata_cfitsio: average of kappa (%ld particles) = %le\n",
               cmd->nbody, kavg/((real)cmd->nbody) );

    //B If needed locate particles with same position.
    //  This is a slow process, use if it necessary...
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
        if (flag) error("inputdata_cfitsio: at least two bodies have same position\n");
    }
    //E
#endif
}

local int inputdata_cfitsio_xyz(struct cmdline_data* cmd, struct  global_data* gd,
                               string filename, int ifile, fitsfile *fptr)
{
    bodyptr p;
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

    datatype = 82;                                      // double
    firstrow = 1;
    firstelem = 1;
    nelements = cmd->nbody;
    arrayKappa = (double*) allocate(cmd->nbody * sizeof(double));

    colnum = gd->columns[3];                            // kappa
    verb_print(cmd->verbose,
               "\tinputdata_cfitsio: colnum: %d...\n",
               colnum);
    fits_read_col(fptr, datatype, colnum, firstrow, firstelem,
                  nelements, &nulval, arrayKappa, &anynul, &status);

    bodytable[ifile] = (bodyptr) allocate(cmd->nbody * sizeof(body));
    DO_BODY(p, bodytable[ifile], bodytable[ifile]+cmd->nbody) {
        Kappa(p) = arrayKappa[p-bodytable[ifile]];
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

    DO_BODY(p, bodytable[ifile], bodytable[ifile]+cmd->nbody) {
        Type(p) = BODY;
        Weight(p) = weight;
        Id(p) = p-bodytable[ifile]+1;
    }

    gd->nbodyTable[ifile] = cmd->nbody;

    real kavg=0.0;
    DO_BODY(p, bodytable[ifile], bodytable[ifile]+cmd->nbody) {
        Type(p) = BODY;
        Weight(p) = weight;
        Id(p) = p-bodytable[ifile]+1;
        kavg += Kappa(p);
    }
    verb_print(cmd->verbose,
               "inputdata_cfitsio: average of kappa (%ld particles) = %le\n",
               cmd->nbody, kavg/((real)cmd->nbody) );

    //B If needed locate particles with same position.
    //  This is a slow process, use if it necessary...
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
        if (flag) error("inputdata_cfitsio: at least two bodies have same position\n");
    }
    //E
#endif
}

#endif	// ! _cballsio_cfitsio_02_h
