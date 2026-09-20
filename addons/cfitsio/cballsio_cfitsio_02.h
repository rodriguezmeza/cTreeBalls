
//=============================================================================
//        1          2          3          4        ^ 5          6          7

#ifndef _cballsio_cfitsio_02_h
#define _cballsio_cfitsio_02_h

#include <ctype.h>

/* Preserve the first read error while always closing with a fresh status. */
static int fits_failure(struct cmdline_data *cmd, int status, const char *where);

local int cfitsio_close_input(struct cmdline_data *cmd, fitsfile *fptr, int rc)
{
    int close_status = 0;
    if (fptr != NULL && fits_close_file(fptr, &close_status) && rc == SUCCESS)
        return fits_failure(cmd, close_status, "fits_close_file (input)");
    return rc;
}

/* All table readers share header validation and ownership of the open handle.
   Never consume metadata after CFITSIO reports a failure. */
local int cfitsio_open_table(struct cmdline_data *cmd, string filename,
                            fitsfile **fptr)
{
    int status = 0, nkeys = 0, ncols = 0;
    char card[FLEN_CARD];
    *fptr = NULL;
    if (scanopt(cmd->options, "fits-type-file"))
        fits_open_file(fptr, filename, READONLY, &status);
    else
        fits_open_data(fptr, filename, READONLY, &status);
    if (status) goto fail;
    fits_get_hdrspace(*fptr, &nkeys, NULL, &status);
    if (status) goto fail;
    fits_get_num_rows(*fptr, &cmd->nbody, &status);
    if (status) goto fail;
    fits_get_num_cols(*fptr, &ncols, &status);
    if (status) goto fail;
    if (cmd->nbody < 1) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "FITS input '%s': table must contain at least one row", filename);
        cfitsio_close_input(cmd, *fptr, FAILURE);
        *fptr = NULL;
        return FAILURE;
    }
    if (scanopt(cmd->options, "header-info")) {
        for (int ii = 1; ii <= nkeys; ii++) {
            fits_read_record(*fptr, ii, card, &status);
            if (status) goto fail;
            printf("%s\n", card);
        }
        printf("END\n\n");
        for (int ii = 1; ii <= ncols; ii++) {
            int typecode;
            long repeat, width;
            fits_get_coltype(*fptr, ii, &typecode, &repeat, &width, &status);
            if (status) goto fail;
            verb_print(cmd->verbose, "%d: typecode, repeat, width = %d %ld %ld\n",
                       ii, typecode, repeat, width);
        }
    }
    verb_print(cmd->verbose, "FITS input: nbody = %ld, columns = %d\n",
               (long)cmd->nbody, ncols);
    return SUCCESS;
fail:
    fits_failure(cmd, status, filename);
    cfitsio_close_input(cmd, *fptr, FAILURE);
    *fptr = NULL;
    return FAILURE;
}

local int cfitsio_column_error(struct cmdline_data *cmd, int status,
                               int column, string filename)
{
    fits_report_error(stderr, status);
    snprintf(cmd->error_message, _ERRORMSGSIZE_,
             "FITS input '%s': column %d failed with status=%d",
             filename, column, status);
    return FAILURE;
}

// input in: addons/source/cballsio/cballsio_include_11a.h

local int inputdata_cfitsio(struct cmdline_data* cmd,
                               struct global_data* gd, string filename, int ifile)
{
    fitsfile *fptr = NULL;
    int rc;
    gd->input_comment = "fits input file";
    if (cfitsio_open_table(cmd, filename, &fptr) == FAILURE)
        return FAILURE;
    if (scanopt(cmd->options, "stop-fits") && strnull(cmd->outfile)) {
        rc = cfitsio_close_input(cmd, fptr, SUCCESS);
        if (rc == FAILURE) return FAILURE;
        gd->inputHeaderFlag = TRUE;
        gd->stopflag = TRUE;
        return SUCCESS;
    }
    rc = inputdata_cfitsio_xyz(cmd, gd, filename, ifile, fptr);
    return cfitsio_close_input(cmd, fptr, rc);
}

// Routine to read ra-dec-field fits files
//  use columns=1,2,3,4 -> to adjust ra-dec-field-weight information
local int inputdata_cfitsio_radec_field(struct cmdline_data* cmd,
                               struct global_data* gd, string filename, int ifile)
{
    fitsfile *fptr = NULL;
    int rc;
    gd->input_comment = "fits-radec-field input file";
    if (cfitsio_open_table(cmd, filename, &fptr) == FAILURE)
        return FAILURE;
    if (scanopt(cmd->options, "stop-fits") && strnull(cmd->outfile)) {
        rc = cfitsio_close_input(cmd, fptr, SUCCESS);
        if (rc == FAILURE) return FAILURE;
        gd->inputHeaderFlag = TRUE;
        gd->stopflag = TRUE;
        return SUCCESS;
    }
    rc = inputdata_cfitsio_ra_dec(cmd, gd, filename, ifile, fptr);
    return cfitsio_close_input(cmd, fptr, rc);
}

// Routine to read ra-dec-field fits files
//  columns order is field,ra,dec,weight (angles in radians)
//  Check header with:
//      options=header-info,stop-fits
//  and if needed add to options 'with-weight'
local int inputdata_cfitsio_ra_dec(struct cmdline_data* cmd,
                                   struct  global_data* gd,
                                   string filename, int ifile, fitsfile *fptr)
{
    string routineName = "inputdata_cfitsio_ra_dec";
    bodyptr p;
    double *arrayKappa = NULL;
    double *arrayRA = NULL;
    double *arrayDEC = NULL;
    double *arrayWEIGHT = NULL;

    real mass=1;
    real weight=1;

    int datatype;
    int colnum;
    LONGLONG firstrow;
    LONGLONG firstelem;
    LONGLONG nelements;
    double nulval = 0.0;
    int anynul;
    int status = 0;

    verb_print_debug(1, "\n%s: columns: %d %d %d %d\n",
                     routineName,
                     gd->columns[0], gd->columns[1],
                     gd->columns[2], gd->columns[3]);

    datatype = TDOUBLE;
    colnum = gd->columns[0];
    firstrow = 1;
    firstelem = 1;
    nelements = cmd->nbody;
    arrayKappa = (double*) allocate(cmd->nbody * sizeof(double));
    fits_read_col(fptr, datatype, colnum, firstrow, firstelem,
                  nelements, &nulval, arrayKappa, &anynul, &status);
    if (status) goto fits_read_fail;
    bodytable[ifile] = (bodyptr) allocate(cmd->nbody * sizeof(body));
    gd->bodytable_allocated = TRUE;
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
    arrayKappa = NULL;

#if NDIM == 3
    real ra, dec;
    arrayRA = (double*) allocate(cmd->nbody * sizeof(double));
    arrayDEC = (double*) allocate(cmd->nbody * sizeof(double));

    //E
    colnum = gd->columns[1];
    fits_read_col(fptr, datatype, colnum, firstrow, firstelem,
                  nelements, &nulval, arrayRA, &anynul, &status);
    if (status) goto fits_read_fail;

    colnum = gd->columns[2];
    fits_read_col(fptr, datatype, colnum, firstrow, firstelem,
                  nelements, &nulval, arrayDEC, &anynul, &status);
    if (status) goto fits_read_fail;

    if (scanopt(cmd->options, "with-weight")) {
        arrayWEIGHT = (double*) allocate(cmd->nbody * sizeof(double));
        colnum = gd->columns[3];
        fits_read_col(fptr, datatype, colnum, firstrow, firstelem,
                      nelements, &nulval, arrayWEIGHT, &anynul, &status);
        if (status) goto fits_read_fail;
    }

    if (scanopt(cmd->options, "no-arfken")) {
        DO_BODY(p, bodytable[ifile], bodytable[ifile]+cmd->nbody) {
            ra = arrayRA[p-bodytable[ifile]];
            dec = arrayDEC[p-bodytable[ifile]];
            Pos(p)[0] = rcos(dec)*rcos(ra);
            Pos(p)[1] = rcos(dec)*rsin(ra);
            Pos(p)[2] = rsin(dec);
            if (scanopt(cmd->options, "with-weight")) {
                Weight(p) = arrayWEIGHT[p-bodytable[ifile]];
            } else {
                Weight(p) = weight;
            }
        }
    } else {
        DO_BODY(p, bodytable[ifile], bodytable[ifile]+cmd->nbody) {
            ra = arrayRA[p-bodytable[ifile]];
            dec = arrayDEC[p-bodytable[ifile]];
            Pos(p)[0] = rsin(dec)*rcos(ra);
            Pos(p)[1] = rsin(dec)*rsin(ra);
            Pos(p)[2] = rcos(dec);
            if (scanopt(cmd->options, "with-weight")) {
                Weight(p) = arrayWEIGHT[p-bodytable[ifile]];
            } else {
                Weight(p) = weight;
            }
        }
    }

    free(arrayWEIGHT);
    arrayWEIGHT = NULL;
    free(arrayDEC);
    arrayDEC = NULL;
    free(arrayRA);
    arrayRA = NULL;

    gd->nbodyTable[ifile] = cmd->nbody;

    real kavg=0.0;
    DO_BODY(p, bodytable[ifile], bodytable[ifile]+gd->nbodyTable[ifile]) {
        Type(p) = BODY;
        Mass(p) = mass;
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
            cBALLS_FAIL(cmd,
                "inputdata_cfitsio: at least two bodies have same position\n");
    }
    //E
#else
    verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
        "\n%s: `DEFDIMENSION` is not set to 3. Do so in Makefile_settings\n",
        routineName);

#endif

    return SUCCESS;

fits_read_fail:
    free(arrayKappa);
    free(arrayRA);
    free(arrayDEC);
    free(arrayWEIGHT);
    return cfitsio_column_error(cmd, status, colnum, filename);
}

//B RADECR_FIELD
// Routine to read ra-dec-field fits files
//  use columns=1,2,3,4 -> to adjust ra-dec-r-field-weight information
local int inputdata_cfitsio_radecr_field(struct cmdline_data* cmd,
                               struct global_data* gd, string filename, int ifile)
{
    fitsfile *fptr = NULL;
    int rc;
    gd->input_comment = "fits-radecr-field input file";
    if (cfitsio_open_table(cmd, filename, &fptr) == FAILURE)
        return FAILURE;
    if (scanopt(cmd->options, "stop-fits") && strnull(cmd->outfile)) {
        rc = cfitsio_close_input(cmd, fptr, SUCCESS);
        if (rc == FAILURE) return FAILURE;
        gd->inputHeaderFlag = TRUE;
        gd->stopflag = TRUE;
        return SUCCESS;
    }
    rc = inputdata_cfitsio_ra_dec_r(cmd, gd, filename, ifile, fptr);
    return cfitsio_close_input(cmd, fptr, rc);
}


// Routine to read ra-dec-r-field fits files
//  columns order is ra,dec,r,field,weight
//  Check header with:
//      options=header-info,stop-fits
//  and if needed add to options 'with-weight'
local int inputdata_cfitsio_ra_dec_r(struct cmdline_data* cmd,
                                     struct  global_data* gd,
                                     string filename, int ifile,
                                     fitsfile *fptr)
{
    string routineName = "inputdata_cfitsio_ra_dec_r";
    bodyptr p;
    double *arrayKappa = NULL;
    double *arrayRA = NULL;
    double *arrayDEC = NULL;
    double *arrayR = NULL;
    double *arrayWEIGHT = NULL;
    bodyptr bodytabtmp = NULL;

    real mass=1;
    real weight=1;

    int datatype;
    int colnum;
    LONGLONG firstrow;
    LONGLONG firstelem;
    LONGLONG nelements;
    double nulval = 0.0;
    int anynul;
    int status = 0;
    INTEGER nrows;

    int typecode;
    long repeat;
    long width;

    verb_print_debug(1,
                     "\n%s: columns (ra-dec-r-field-weight): %d %d %d %d %d\n",
                     routineName,
                     gd->columns[0], gd->columns[1],
                     gd->columns[2], gd->columns[3], gd->columns[4]);

    //B arrayKappa
    datatype = 82;                                  // TDOUBLE
    colnum = gd->columns[3];
    verb_print(cmd->verbose, "Column info details (Kappa):\n");
    fits_get_coltype(fptr, colnum, &typecode, &repeat, &width, &status);
    if (status) goto fits_read_fail;
    switch(typecode) {
        case TLONG:
            verb_print(cmd->verbose,
                       "%d: typecode, repeat, width = %d %s %ld %ld\n",
                       colnum, typecode, "TLONG", repeat, width);
            break;
        case TFLOAT:
            verb_print(cmd->verbose,
                       "%d: typecode, repeat, width = %d %s %ld %ld\n",
                       colnum, typecode, "TFLOAT", repeat, width);
            break;
        case TDOUBLE:
            verb_print(cmd->verbose,
                       "%d: typecode, repeat, width = %d %s %ld %ld\n",
                       colnum, typecode, "TDOUBLE", repeat, width);
            break;
    }
    firstrow = 1;
    firstelem = 1;
    nelements = cmd->nbody*repeat;
    INTEGER nelementsKappa = nelements;

    verb_print_debug(1, "\n%s: rows, nelements: %ld %ld\n",
                     routineName, cmd->nbody, nelements);
    cmd->nbody = nelements;
    bodytabtmp = (bodyptr) allocate(cmd->nbody * sizeof(body));

    arrayKappa = (double*) allocate(nelements * sizeof(double));
    fits_read_col(fptr, datatype, colnum, firstrow, firstelem,
                  nelements, &nulval, arrayKappa, &anynul, &status);
    if (status) goto fits_read_fail;

    INTEGER nonvalid=0;
    INTEGER valid=0;
    INTEGER ij=1;
    INTEGER ip=1;                                   // to track row
    DO_BODY(p, bodytabtmp, bodytabtmp+cmd->nbody) {
        Update(p) = FALSE;
        Mask(p) = TRUE;                             // initialize body's Mask
        if (arrayKappa[p-bodytabtmp] < SMALLESTDOUBLE) {
            nonvalid++;
        } else {
            valid++;
            Update(p) = TRUE;
        }
        Kappa(p) = arrayKappa[p-bodytabtmp];
        if (scanopt(cmd->options, "kappa-constant")) {
            if (scanopt(cmd->options, "kappa-constant-one"))
                Kappa(p) = 1.0;
            else
                Kappa(p) = 2.0;
        }
        if (ij%repeat == 0) {
            Id(p) = ip++;
        }
        ij++;
    }
    nrows = ip-1;
    verb_print_debug(1,
                "\n%s: valid, nonvalid, nelements, nrows: %ld %ld %ld %ld\n",
                routineName, valid, nonvalid, nelements, nrows);
    //B not needed this line... check
    nelements = cmd->nbody;
    //E
    int repeatKappa = repeat;
#ifdef DIRECTMETHODSIMPLELOOPID
#include "iolib_direct_method_simple_loopId.h"
#endif
    //E arrayKappa

#if NDIM == 3
    real ra = 0.0;
    real dec = 0.0;

    //B arrayRA
    colnum = gd->columns[0];
    verb_print(cmd->verbose, "Column info details (RA):\n");
    fits_get_coltype(fptr, colnum, &typecode, &repeat, &width, &status);
    if (status) goto fits_read_fail;
    nelements = nrows*repeat;
    arrayRA = (double*) allocate(nelements * sizeof(double));
    switch(typecode) {
        case TLONG:
            verb_print(cmd->verbose,
                       "%d: typecode, repeat, width = %d %s %ld %ld\n",
                       colnum, typecode, "TLONG", repeat, width);
            break;
        case TFLOAT:
            verb_print(cmd->verbose,
                       "%d: typecode, repeat, width = %d %s %ld %ld\n",
                       colnum, typecode, "TFLOAT", repeat, width);
            break;
        case TDOUBLE:
            verb_print(cmd->verbose,
                       "%d: typecode, repeat, width = %d %s %ld %ld\n",
                       colnum, typecode, "TDOUBLE", repeat, width);
            break;
    }
    fits_read_col(fptr, datatype, colnum, firstrow, firstelem,
                  nelements, &nulval, arrayRA, &anynul, &status);
    if (status) goto fits_read_fail;

    //E

    //B arrayDEC
    colnum = gd->columns[1];
    verb_print(cmd->verbose, "Column info details (DEC):\n");
    fits_get_coltype(fptr, colnum, &typecode, &repeat, &width, &status);
    if (status) goto fits_read_fail;
    nelements = nrows*repeat;
    arrayDEC = (double*) allocate(nelements * sizeof(double));
    switch(typecode) {
        case TLONG:
            verb_print(cmd->verbose,
                       "%d: typecode, repeat, width = %d %s %ld %ld\n",
                       colnum, typecode, "TLONG", repeat, width);
            break;
        case TFLOAT:
            verb_print(cmd->verbose,
                       "%d: typecode, repeat, width = %d %s %ld %ld\n",
                       colnum, typecode, "TFLOAT", repeat, width);
            break;
        case TDOUBLE:
            verb_print(cmd->verbose,
                       "%d: typecode, repeat, width = %d %s %ld %ld\n",
                       colnum, typecode, "TDOUBLE", repeat, width);
            break;
    }
    fits_read_col(fptr, datatype, colnum, firstrow, firstelem,
                  nelements, &nulval, arrayDEC, &anynul, &status);
    if (status) goto fits_read_fail;

    //E

    //B arrayR
    colnum = gd->columns[2];
    verb_print(cmd->verbose, "Column info details (R):\n");
    fits_get_coltype(fptr, colnum, &typecode, &repeat, &width, &status);
    if (status) goto fits_read_fail;
    nelements = nrows*repeat;
    arrayR = (double*) allocate(nelements * sizeof(double));
    switch(typecode) {
        case TLONG:
            verb_print(cmd->verbose,
                       "%d: typecode, repeat, width = %d %s %ld %ld\n",
                       colnum, typecode, "TLONG", repeat, width);
            break;
        case TFLOAT:
            verb_print(cmd->verbose,
                       "%d: typecode, repeat, width = %d %s %ld %ld\n",
                       colnum, typecode, "TFLOAT", repeat, width);
            break;
        case TDOUBLE:
            verb_print(cmd->verbose,
                       "%d: typecode, repeat, width = %d %s %ld %ld\n",
                       colnum, typecode, "TDOUBLE", repeat, width);
            break;
    }
    fits_read_col(fptr, datatype, colnum, firstrow, firstelem,
                  nelements, &nulval, arrayR, &anynul, &status);
    if (status) goto fits_read_fail;


    INTEGER nonvalidR=0;
    INTEGER validR=0;
    ij=1;
    real r;
    real rmin=BIGGESTDOUBLE, rmax=0.;
    DO_BODY(p, bodytabtmp, bodytabtmp+cmd->nbody) {
        r = arrayR[p-bodytabtmp];
        if (Update(p)==TRUE && r < SMALLESTRMPC) {
            nonvalidR++;
        } else {
            validR++;
        }
        if (Update(p)) {
            rmax = MAX(rmax,r);
            rmin = MIN(rmin,r);
        }
        ij++;
    }
    verb_print_debug(1,
                "\n%s: validR, nonvalidR, nelements, nrows: %ld %ld %ld %ld\n",
                routineName, validR, nonvalidR, nelements, nrows);
    verb_print(cmd->verbose, "\t%s: min and max of r = %g %g\n",
               routineName, rmin, rmax);

    //E arrayR

    //B arrayWEIGHT
    if (scanopt(cmd->options, "with-weight")) {
        colnum = gd->columns[4];
        verb_print(cmd->verbose, "Column info details (weight):\n");
        fits_get_coltype(fptr, colnum, &typecode, &repeat, &width, &status);
        if (status) goto fits_read_fail;
        nelements = nrows*repeat;
        arrayWEIGHT = (double*) allocate(nelements * sizeof(double));
        switch(typecode) {
            case TLONG:
                verb_print(cmd->verbose,
                           "%d: typecode, repeat, width = %d %s %ld %ld\n",
                           colnum, typecode, "TLONG", repeat, width);
                break;
            case TFLOAT:
                verb_print(cmd->verbose,
                           "%d: typecode, repeat, width = %d %s %ld %ld\n",
                           colnum, typecode, "TFLOAT", repeat, width);
                break;
            case TDOUBLE:
                verb_print(cmd->verbose,
                           "%d: typecode, repeat, width = %d %s %ld %ld\n",
                           colnum, typecode, "TDOUBLE", repeat, width);
                break;
        }
        fits_read_col(fptr, datatype, colnum, firstrow, firstelem,
                      nelements, &nulval, arrayWEIGHT, &anynul, &status);
        if (status) goto fits_read_fail;
    }

    int repeatWeight = repeat;
    if (scanopt(cmd->options, "with-weight") && repeatKappa != repeatWeight) {
        status = BAD_TFORM;
        goto fits_read_fail;
    }
    //E

    ij=1;
    ip=1;                                           // to track row
    real ramin=0., ramax=0.;
    real decmin=0., decmax=0.;
    rmin=BIGGESTDOUBLE, rmax=0.;
    if (scanopt(cmd->options, "no-arfken")) {
        DO_BODY(p, bodytabtmp, bodytabtmp+cmd->nbody) {
            if (ij%repeatKappa == 0) {
                ra = arrayRA[ip-1] * 60.0/RADTOARCMIN;
                dec = arrayDEC[ip-1] * 60.0/RADTOARCMIN;
                ramin = MIN(ramin,ra);
                decmin = MIN(decmin,dec);
                ramax = MAX(ramax,ra);
                decmax = MAX(decmax,dec);
                ip++;
            }
            r = arrayR[p-bodytabtmp];
            if(Update(p)) {
                rmax = MAX(rmax,r);
                rmin = MIN(rmin,r);
            }
            Pos(p)[0] = r*rcos(dec)*rcos(ra);
            Pos(p)[1] = r*rcos(dec)*rsin(ra);
            Pos(p)[2] = r*rsin(dec);
            if (scanopt(cmd->options, "with-weight")) {
                Weight(p) = arrayWEIGHT[p-bodytabtmp];
            } else {
                Weight(p) = weight;
            }
            ij++;
        }
    } else {
        DO_BODY(p, bodytabtmp, bodytabtmp+cmd->nbody) {
            if (ij%repeatKappa == 0) {
                ra = arrayRA[ip-1] * 60.0/RADTOARCMIN;
                dec = arrayDEC[ip-1] * 60.0/RADTOARCMIN;
                ramin = MIN(ramin,ra);
                decmin = MIN(decmin,dec);
                ramax = MAX(ramax,ra);
                decmax = MAX(decmax,dec);
                ip++;
            }
            Id(p) = ip-1;
            r = arrayR[p-bodytabtmp];
            //B this segment will give too much bodies with equal positions...
            //      randomize positions a little bit here instead in treeload
            if (scanopt(cmd->options, "ra-dec-only"))
                r = 1.0;
            //E
            if(Update(p)) {
                rmax = MAX(rmax,r);
                rmin = MIN(rmin,r);
            }
            Pos(p)[0] = r*rsin(dec)*rcos(ra);
            Pos(p)[1] = r*rsin(dec)*rsin(ra);
            Pos(p)[2] = r*rcos(dec);
            if (scanopt(cmd->options, "with-weight")) {
                Weight(p) = arrayWEIGHT[p-bodytabtmp];
            } else {
                Weight(p) = weight;
            }
            ij++;
        }
    }
    verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,"\n");
    verb_print(cmd->verbose, "\n\t%s: min and max of ra = %f %f\n",
               routineName, ramin, ramax);
    verb_print(cmd->verbose, "\t%s: min and max of dec = %f %f\n",
               routineName, decmin, decmax);
    verb_print(cmd->verbose, "\t%s: min and max of r = %f %f\n",
               routineName, rmin, rmax);

    free(arrayWEIGHT);
    arrayWEIGHT = NULL;
    free(arrayR);
    arrayR = NULL;
    free(arrayDEC);
    arrayDEC = NULL;
    free(arrayRA);
    arrayRA = NULL;
    free(arrayKappa);
    arrayKappa = NULL;

    cmd->nbody = valid;
    gd->nbodyTable[ifile] = cmd->nbody;
    bodyptr q;
    bodytable[ifile] = (bodyptr) allocate(cmd->nbody * sizeof(body));
    gd->bodytable_allocated = TRUE;
    verb_print(cmd->verbose_log, "\n%s: Created bodies = %ld",
               routineName, cmd->nbody);

    real kavg=0.0;
    ij = 0;
    ip=1;                                   // to track row
    real xmin, xmax;
    real ymin, ymax;
    real zmin, zmax;
    xmin=0., ymin=0.0, zmin=0.;
    xmax=0., ymax=0., zmax=0.;
    verb_print(cmd->verbose,
               "\t%s: columns 5 : %d\n",
               routineName, gd->columns[5]);
    DO_BODY(q, bodytabtmp, bodytabtmp+nelementsKappa) {
        if(Update(q)) {
            p = bodytable[ifile]+ij;
            Pos(p)[0] = Pos(q)[0];
            Pos(p)[1] = Pos(q)[1];
            Pos(p)[2] = Pos(q)[2];
            Kappa(p) = Kappa(q);

            Type(p) = Type(q);
            Mass(p) = mass;
            Weight(p) = Weight(q);
            Mask(p) = Mask(q);
            Id(p) = ij+1;
            xmin = MIN(xmin,Pos(p)[0]);
            ymin = MIN(ymin,Pos(p)[1]);
            zmin = MIN(zmin,Pos(p)[2]);
            xmax = MAX(xmax,Pos(p)[0]);
            ymax = MAX(ymax,Pos(p)[1]);
            zmax = MAX(zmax,Pos(p)[2]);
            ij++;
            kavg += Kappa(p);
            Update(p) = Update(q);
            ip++;
        }
    }
    verb_print(cmd->verbose,
               "\t%s: selected kappas = %ld\n",
               routineName, ij);
    verb_print(cmd->verbose,
               "\t%s: average of kappa (%ld particles) = %le\n",
               routineName, cmd->nbody, kavg/((real)gd->nbodyTable[ifile]) );

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
            cBALLS_FAIL(cmd, "%s: at least two bodies have same position\n",
                        routineName);
    }
    //E

    verb_print(cmd->verbose, "\n\t%s: min and max of x = %f %f\n",
               routineName, xmin, xmax);
    verb_print(cmd->verbose, "\t%s: min and max of y = %f %f\n",
               routineName, ymin, ymax);
    verb_print(cmd->verbose, "\t%s: min and max of z = %f %f\n",
               routineName, zmin, zmax);

    verb_print(cmd->verbose,
               "\n\t%s: selected read points (out of nbody=%ld) = %ld\n",
               routineName, cmd->nbody, valid);

    free(bodytabtmp);

#else
#error `DEFDIMENSION` is not set to 3. Do so in Makefile_settings
//    verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
//        "\n%s: `DEFDIMENSION` is not set to 3. Do so in Makefile_settings\n",
//        routineName);
#endif

    return SUCCESS;

fits_read_fail:
    free(arrayKappa);
    free(arrayRA);
    free(arrayDEC);
    free(arrayR);
    free(arrayWEIGHT);
    free(bodytabtmp);
    return cfitsio_column_error(cmd, status, colnum, filename);
}
//E RADECR_FIELD

#if (defined(OCTREE3PCF3DOMP) || defined(OCTREE3PCF3DMPI)) && NDIM == 3
local int cb3d_cfitsio_count_columns_tokens(string columns)
{
    int count = 0;
    int in_token = FALSE;
    const char *cursor;

    if (columns == NULL) return 0;
    for (cursor = columns; *cursor != '\0'; cursor++) {
        if (*cursor == ',' || isspace((unsigned char)*cursor)) {
            in_token = FALSE;
        } else if (!in_token) {
            count++;
            in_token = TRUE;
        }
    }
    return count;
}

local int cb3d_cfitsio_exclude_same_los_option(struct cmdline_data* cmd)
{
    return scanopt(cmd->options, "exclude-same-los")
        || scanopt(cmd->options, "exclude-los")
        || scanopt(cmd->options, "exclude-pivot-los")
        || scanopt(cmd->options, "exclude-all-same-los")
        || scanopt(cmd->options, "lya-distinct-forests");
}
#endif

local int inputdata_cfitsio_xyz(struct cmdline_data* cmd,
                                struct  global_data* gd,
                                string filename, int ifile, fitsfile *fptr)
{
    bodyptr p;
    double *arrayKappa = NULL;
    double *arrayX = NULL;
    double *arrayY = NULL;
    double *arrayZ = NULL;

    real mass=1;
    real weight=1;

    int datatype;
    int colnum;
    LONGLONG firstrow;
    LONGLONG firstelem;
    LONGLONG nelements;
    double nulval = 0.0;
    int anynul;
    int status = 0;
#if (defined(OCTREE3PCF3DOMP) || defined(OCTREE3PCF3DMPI)) && NDIM == 3
    double *arrayWeight = NULL;
    LONGLONG *arrayLosId = NULL;
    int read_weight = scanopt(cmd->options, "with-weight");
    int read_los_id = cb3d_cfitsio_count_columns_tokens(cmd->columns) >= 6;
    int invalid_los_id = FALSE;

    if (cb3d_cfitsio_exclude_same_los_option(cmd) && !read_los_id)
        cBALLS_FAIL(cmd,
                    "inputdata_cfitsio_xyz: LOS exclusion requires six FITS columns: x,y,z,delta,weight,los_id");
#endif

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
    if (status) goto fits_read_fail;

    bodytable[ifile] = (bodyptr) allocate(cmd->nbody * sizeof(body));
    gd->bodytable_allocated = TRUE;
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
    arrayKappa = NULL;

#if NDIM == 3
    arrayX = (double*) allocate(cmd->nbody * sizeof(double));
    arrayY = (double*) allocate(cmd->nbody * sizeof(double));
    arrayZ = (double*) allocate(cmd->nbody * sizeof(double));

    colnum = gd->columns[0];
    verb_print(cmd->verbose,
               "\tinputdata_cfitsio: colnum X: %d...\n",
               colnum);
    fits_read_col(fptr, datatype, colnum, firstrow, firstelem,
                  nelements, &nulval, arrayX, &anynul, &status);
    if (status) goto fits_read_fail;

    
    colnum = gd->columns[1];
    verb_print(cmd->verbose,
               "\tinputdata_cfitsio: colnum Y: %d...\n",
               colnum);
    fits_read_col(fptr, datatype, colnum, firstrow, firstelem,
                  nelements, &nulval, arrayY, &anynul, &status);
    if (status) goto fits_read_fail;


    colnum = gd->columns[2];
    verb_print(cmd->verbose,
               "\tinputdata_cfitsio: colnum Z: %d...\n\n",
               colnum);
    fits_read_col(fptr, datatype, colnum, firstrow, firstelem,
                  nelements, &nulval, arrayZ, &anynul, &status);
    if (status) goto fits_read_fail;


#if defined(OCTREE3PCF3DOMP) || defined(OCTREE3PCF3DMPI)
    if (read_weight) {
        colnum = gd->columns[4];
        arrayWeight = (double*) allocate(cmd->nbody * sizeof(double));
        fits_read_col(fptr, datatype, colnum, firstrow, firstelem,
                      nelements, &nulval, arrayWeight, &anynul, &status);
        if (status) goto fits_read_fail;
    }
    if (read_los_id) {
        LONGLONG nulval_los = 0;
        colnum = gd->columns[5];
        arrayLosId = (LONGLONG*) allocate(cmd->nbody * sizeof(LONGLONG));
        fits_read_col(fptr, TLONGLONG, colnum, firstrow, firstelem,
                      nelements, &nulval_los, arrayLosId, &anynul, &status);
        if (status) goto fits_read_fail;
    }

#endif

    DO_BODY(p, bodytable[ifile], bodytable[ifile]+cmd->nbody) {
        size_t index = (size_t)(p-bodytable[ifile]);
        Pos(p)[0] = arrayX[index];
        Pos(p)[1] = arrayY[index];
        Pos(p)[2] = arrayZ[index];
#if defined(OCTREE3PCF3DOMP) || defined(OCTREE3PCF3DMPI)
        Weight(p) = read_weight ? arrayWeight[index] : weight;
        if (arrayLosId != NULL) {
            Octree3pcf3dLosId(p) = (INTEGER)arrayLosId[index];
            if ((LONGLONG)Octree3pcf3dLosId(p) != arrayLosId[index])
                invalid_los_id = TRUE;
        } else {
            Octree3pcf3dLosId(p) = (INTEGER)index + 1;
        }
        Mask(p) = MASK_NODE_VALID;
#endif
    }

    free(arrayX);
    arrayX = NULL;
    free(arrayY);
    arrayY = NULL;
    free(arrayZ);
    arrayZ = NULL;
#if defined(OCTREE3PCF3DOMP) || defined(OCTREE3PCF3DMPI)
    free(arrayWeight);
    free(arrayLosId);
    if (invalid_los_id)
        cBALLS_FAIL(cmd,
                    "inputdata_cfitsio_xyz: LOS_ID cannot be represented by INTEGER");
    gd->octree3pcf3d_los_ids[ifile] = TRUE;
#endif

    gd->nbodyTable[ifile] = cmd->nbody;

    real kavg=0.0;
    DO_BODY(p, bodytable[ifile], bodytable[ifile]+gd->nbodyTable[ifile]) {
        Type(p) = BODY;
        Mass(p) = mass;
#if !defined(OCTREE3PCF3DOMP) && !defined(OCTREE3PCF3DMPI)
        Weight(p) = weight;
#endif
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
            cBALLS_FAIL(cmd, "inputdata_cfitsio: %s\n",
                        "at least two bodies have same position");
    }
    //E
#else
    cBALLS_FAIL(cmd, "inputdata_cfitsio_xyz: %s\n\n",
                "this routine works only in 3D. exiting...");
#endif

    return SUCCESS;

fits_read_fail:
    free(arrayKappa);
    free(arrayX);
    free(arrayY);
    free(arrayZ);
#if (defined(OCTREE3PCF3DOMP) || defined(OCTREE3PCF3DMPI)) && NDIM == 3
    free(arrayWeight);
    free(arrayLosId);
#endif
    return cfitsio_column_error(cmd, status, colnum, filename);
}

// Routine to read fits-healpix files
//  so far only RING scheme
local int inputdata_cfitsio_healpix(struct cmdline_data* cmd,
                                   struct global_data* gd, string filename, int ifile)
{
    string routineName = "inputdata_cfitsio_healpix";
    fitsfile *fptr = NULL;
    int rc;
    gd->input_comment = "fits-healpix input file";
    if (cfitsio_open_table(cmd, filename, &fptr) == FAILURE)
        return FAILURE;
    if (scanopt(cmd->options, "stop-fits") && strnull(cmd->outfile)) {
        rc = cfitsio_close_input(cmd, fptr, SUCCESS);
        if (rc == FAILURE) return FAILURE;
        gd->inputHeaderFlag = TRUE;
        gd->stopflag = TRUE;
        return SUCCESS;
    }
    if (scanopt(cmd->options, "read-mask")) {
        if (ifile == 0) {
            rc = inputdata_cfitsio_healpix_map(cmd, gd, filename, ifile, fptr);
        } else {
            if (ifile != 1) {
                cfitsio_close_input(cmd, fptr, FAILURE);
                cBALLS_FAIL(cmd, "\t%s: read-mask ifile = %d is absurd\n",
                            routineName, ifile);
            }

            rc = inputdata_cfitsio_healpix_map_mask(cmd, gd, filename, ifile, fptr);
        }
    } else {
        if (scanopt(cmd->options, "mask-inside")) {
            rc = inputdata_cfitsio_healpix_map_mask_inside(cmd, gd,
                                                           filename, ifile, fptr);
        } else {
            rc = inputdata_cfitsio_healpix_map(cmd, gd, filename, ifile, fptr);
        }
    }

    return cfitsio_close_input(cmd, fptr, rc);
}

local int cfitsio_healpix_map_to_ring(struct cmdline_data *cmd,
                                      const char *routine_name,
                                      long nside,
                                      const char *size_ordering,
                                      const char *map_ordering,
                                      float **map)
{
    int hp_status;

    if (strcmp(size_ordering, map_ordering) != 0) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: inconsistent HEALPix ordering '%s' and '%s'",
                 routine_name, size_ordering, map_ordering);
        return FAILURE;
    }

    hp_status = healpix_map_to_ring_status(map, nside, map_ordering);
    if (hp_status != 0) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: cannot normalize HEALPix ordering '%s' to RING "
                 "for nside=%ld status=%d",
                 routine_name, map_ordering, nside, hp_status);
        return FAILURE;
    }

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

    int hp_status;
    int rc = SUCCESS;
    
#if THREEDIMCODE
    verb_print(cmd->verbose, "\nWorking 3D map...\n");
#else
    cBALLS_FAIL(cmd, "\nOnly 3D is implemented so far... exiting...\n\n");
#endif

    hp_status = get_fits_size_status(filename, &nside, order1, &npixel);
    if (hp_status != 0) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: get_fits_size failed for '%s' status=%d",
                 routineName, filename, hp_status);
        return FAILURE;
    }
    
    verb_print(cmd->verbose,
        "filename, ifile, nside, npixel, order1:");
    
    //B
    map = NULL;
    hp_status = read_healpix_map_status(filename, &nside, coord, order2, &map);
    if (hp_status != 0 || map == NULL) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: read_healpix_map failed for '%s' status=%d",
                 routineName, filename, hp_status);
        return FAILURE;
    }
    if (cfitsio_healpix_map_to_ring(cmd, routineName, nside,
                                    order1, order2, &map) == FAILURE) {
        free(map);
        return FAILURE;
    }
    
    verb_print(cmd->verbose,
        "%s %d %ld %ld %s\n", filename, ifile, nside, npixel, order1);
    //E
    
    verb_print(cmd->verbose,
               "\nAllocated %g MByte for temporal map pixel (%ld) storage.\n",
               npixel*sizeof(float)/(1024.0*1024.0),
               npixel);

    verb_print(cmd->verbose,
               "%s: nbody = %ld...\n",
               routineName, npixel);
    if (npixel < 1)
        cBALLS_FAIL(cmd, "%s: npixel = %ld is absurd\n", routineName, npixel);

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
        OPEN_OUTPUT_OR_FAIL(outstr, gd->fpfnameOutputFileName, "w!");
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
                        if (out_real_mar_checked(outstr, phi, routineName,
                                                 gd->fpfnameOutputFileName,
                                    cmd->error_message, _ERRORMSGSIZE_) == FAILURE) {
                            if (outstr != NULL) fclose(outstr);
                            return FAILURE;
                        }
                        if (out_real_mar_checked(outstr, theta, routineName,
                                                 gd->fpfnameOutputFileName,
                                    cmd->error_message, _ERRORMSGSIZE_) == FAILURE) {
                            if (outstr != NULL) fclose(outstr);
                            return FAILURE;
                        }
                        if (out_real_checked(outstr, Kappa(p), routineName,
                                             gd->fpfnameOutputFileName,
                                    cmd->error_message, _ERRORMSGSIZE_) == FAILURE) {
                            if (outstr != NULL) fclose(outstr);
                            return FAILURE;
                        }
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
                //B
                if (out_real_mar_checked(outstr, phi, routineName, gd->fpfnameOutputFileName,
                    cmd->error_message, _ERRORMSGSIZE_) == FAILURE) {
                    if (outstr != NULL) fclose(outstr);
                        return FAILURE;
                }
                if (out_real_mar_checked(outstr, theta, routineName, gd->fpfnameOutputFileName,
                    cmd->error_message, _ERRORMSGSIZE_) == FAILURE) {
                    if (outstr != NULL) fclose(outstr);
                        return FAILURE;
                }
                if (out_real_checked(outstr, Kappa(p), routineName, gd->fpfnameOutputFileName,
                    cmd->error_message, _ERRORMSGSIZE_) == FAILURE) {
                    if (outstr != NULL) fclose(outstr);
                        return FAILURE;
                }
                //E
            }
            //E
            Type(p) = BODY;
            Mass(p) = mass;
            Weight(p) = weight;
            Id(p) = p-bodytabtmp+iselect;
            Update(p) = TRUE;
            Mask(p) = TRUE;                             // initialize body's Mask
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
        CLOSE_OUTPUT_OR_FAIL(outstr, gd->fpfnameOutputFileName);
        gd->stopflag = TRUE;
        rc = SUCCESS;
        goto cleanup;
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
    map = NULL;
    verb_print(cmd->verbose,
               "\nFreed %g MByte for temporal map pixel (%ld) storage.\n",
               npixel*sizeof(float)/(1024.0*1024.0),
               npixel);

    bodyptr q;
    if (scanopt(cmd->options, "patch"))
        cmd->nbody = iselect;

    gd->nbodyTable[ifile] = cmd->nbody;
    bodytable[ifile] = (bodyptr) allocate(cmd->nbody * sizeof(body));
    gd->bodytable_allocated = TRUE;

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
            Id(p) = p-bodytable[ifile]+ipix;        // define a counter and use it...
            Update(p) = Update(q);
#if defined(NMultipoles) && defined(NONORMHIST)
            UpdatePivot(p) = UpdatePivot(q);
#endif
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
        // leading ! to allow overwrite
        if (format_checked(fileforce, sizeof(fileforce),
            "fileforce", "!%s/%s", cmd->rootDir, file) != 0)
            return FAILURE;
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
        hp_status = write_healpix_map_status(map, nside, fileforce, 0, "C");
        if (hp_status != 0) {
            free(map);
            map = NULL;
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     "%s: write_healpix_map failed for '%s' status=%d",
                     routineName, fileforce, hp_status);
            return FAILURE;
        }
        
        fprintf(stdout,"\t\tfile written\n");
        free(map);
        map = NULL;
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
    bodytabtmp = NULL;
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

    rc = SUCCESS;

cleanup:
    if (map != NULL) free(map);
    if (bodytabtmp != NULL) free(bodytabtmp);

    return rc;

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

    int hp_status;

    verb_print(cmd->verbose, "\nWorking 2D map...\n");

    hp_status = get_fits_size_status(filename, &nside, order1, &npixel);
    if (hp_status != 0) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: get_fits_size failed for '%s' status=%d",
                 routineName, filename, hp_status);
        return FAILURE;
    }
    
    verb_print(cmd->verbose,
        "filename, ifile, nside, npixel, order1:");
    verb_print(cmd->verbose,
        "%s %d %ld %ld %s\n", filename, ifile, nside, npixel, order1);
    hp_status = read_healpix_map_status(filename, &nside, coord, order2, &map);
    if (hp_status != 0 || map == NULL) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: read_healpix_map failed for '%s' status=%d",
                 routineName, filename, hp_status);
        return FAILURE;
    }
    if (cfitsio_healpix_map_to_ring(cmd, routineName, nside,
                                    order1, order2, &map) == FAILURE) {
        free(map);
        return FAILURE;
    }

    verb_print(cmd->verbose,
               "\nAllocated %g MByte for temporal map pixel (%ld) storage.\n",
               npixel*sizeof(float)*INMB, npixel);

    verb_print(cmd->verbose,
               "%s: nbody = %ld...\n",
               routineName, npixel);
    if (npixel < 1)
        cBALLS_FAIL(cmd, "%s: npixel = %ld is absurd\n", routineName, npixel);

    bodyptr bodytabtmp;
    bodytabtmp = (bodyptr) allocate(npixel * sizeof(body));
    verb_print(cmd->verbose,
               "\nAllocated %g MByte for temporal particle (%ld) storage.\n",
               npixel*sizeof(body)*INMB, npixel);

    //B save optional RA-DEC to a file... 3-columns: RA, DEC, Kappa
    char namebuf[256];
    stream outstr;
    if (scanopt(cmd->options, "save-ra-dec")&&!strnull(cmd->outfile)) {
        OPEN_OUTPUT_OR_FAIL(outstr, gd->fpfnameOutputFileName, "w!");
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
                    //B tranform coordinate if necesarry...
                    Pos(p)[0] = phi;
                    Pos(p)[1] = theta;
                    //E
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
                        //B
                        if (out_real_mar_checked(outstr, Pos(p)[0], routineName, gd->fpfnameOutputFileName,
                            cmd->error_message, _ERRORMSGSIZE_) == FAILURE) {
                            if (outstr != NULL) fclose(outstr);
                                return FAILURE;
                        }
                        if (out_real_mar_checked(outstr, Pos(p)[1], routineName, gd->fpfnameOutputFileName,
                            cmd->error_message, _ERRORMSGSIZE_) == FAILURE) {
                            if (outstr != NULL) fclose(outstr);
                                return FAILURE;
                        }
                        if (out_real_checked(outstr, Kappa(p), routineName,
                                             gd->fpfnameOutputFileName,
                            cmd->error_message, _ERRORMSGSIZE_) == FAILURE) {
                            if (outstr != NULL) fclose(outstr);
                                return FAILURE;
                        }
                        //E
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
            //B tranform coordinate if necesarry...
            Pos(p)[0] = phi;
            Pos(p)[1] = theta;
            //E
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
                //B
                if (out_real_mar_checked(outstr, Pos(p)[0], routineName,
                                         gd->fpfnameOutputFileName,
                    cmd->error_message, _ERRORMSGSIZE_) == FAILURE) {
                    if (outstr != NULL) fclose(outstr);
                        return FAILURE;
                }
                if (out_real_mar_checked(outstr, Pos(p)[1], routineName, gd->fpfnameOutputFileName,
                    cmd->error_message, _ERRORMSGSIZE_) == FAILURE) {
                    if (outstr != NULL) fclose(outstr);
                        return FAILURE;
                }
                if (out_real_checked(outstr, Kappa(p), routineName, gd->fpfnameOutputFileName,
                    cmd->error_message, _ERRORMSGSIZE_) == FAILURE) {
                    if (outstr != NULL) fclose(outstr);
                        return FAILURE;
                }
                //E
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
        CLOSE_OUTPUT_OR_FAIL(outstr, gd->fpfnameOutputFileName);
        gd->stopflag = TRUE;
        return SUCCESS;
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
    gd->bodytable_allocated = TRUE;

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
        // leading ! to allow overwrite
        if (format_checked(fileforce, sizeof(fileforce),
            "fileforce", "!%s/%s", cmd->rootDir, file) != 0)
            return FAILURE;
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
        hp_status = write_healpix_map_status(map, nside, fileforce, 0, "C");
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

//    free(map);

    return SUCCESS;
}

#endif // !  THREEDIMCODE


// reading healpix mask map
//  3D only... and full-sky only
local int inputdata_cfitsio_healpix_map_mask(struct cmdline_data* cmd,
                                        struct  global_data* gd,
                                        string filename,
                                        int ifile, fitsfile *fptr)
{
    string routineName = "inputdata_cfitsio_healpix_map_mask";
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

    int hp_status;

#if THREEDIMCODE
    verb_print(cmd->verbose, "\nWorking 3D map...\n");
#else
    cBALLS_FAIL(cmd, "\nOnly 3D is implemented so far... exiting...\n\n");
#endif

    hp_status = get_fits_size_status(filename, &nside, order1, &npixel);
    if (hp_status != 0) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: get_fits_size failed for '%s' status=%d",
                 routineName, filename, hp_status);
        return FAILURE;
    }
    
    verb_print(cmd->verbose,
        "filename, ifile, nside, npixel, order1:");
    verb_print(cmd->verbose,
        "%s %d %ld %ld %s\n", filename, ifile, nside, npixel, order1);
    hp_status = read_healpix_map_status(filename, &nside, coord, order2, &map);
    if (hp_status != 0 || map == NULL) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: read_healpix_map failed for '%s' status=%d",
                 routineName, filename, hp_status);
        return FAILURE;
    }
    if (cfitsio_healpix_map_to_ring(cmd, routineName, nside,
                                    order1, order2, &map) == FAILURE) {
        free(map);
        return FAILURE;
    }
    verb_print(cmd->verbose,
               "\nAllocated %g MByte for temporal map pixel (%ld) storage.\n",
               npixel*sizeof(float)/(1024.0*1024.0),
               npixel);

    verb_print(cmd->verbose, "%s: nbody = %ld...\n",
               routineName, npixel);

    //B
    if (npixel < 1) {
        free(map);
        cBALLS_FAIL(cmd, "%s: npixel = %ld is absurd\n", routineName, npixel);
    }

    if (npixel != gd->nbodyTable[gd->iCatalogs[0]]) {
        free(map);
        cBALLS_FAIL(cmd,
            "%s: npixel = %ld is not equal to npixel in cat 0: %ld\n",
            routineName, npixel, gd->nbodyTable[gd->iCatalogs[0]]);
    }
    //E

    INTEGER iselect = 0;
    for(ipix=0;ipix<npixel;ipix++) {                // RING loop order
        p = bodytable[gd->iCatalogs[0]]+ipix;
        pix2ang_ring(nside, ipix, &theta, &phi);
        coordinate_transformation(cmd, gd, theta, phi, q);
        DO_COORD(k) {
            if (Pos(p)[k] != q[k]) {
                free(map);
                cBALLS_FAIL(cmd, "%s: mask position is not equal: %g %g\n",
                            routineName, Pos(p)[k], q[k]);
            }
        }

        Mask(p) = map[ipix];
        if (Mask(p) == 0) {
            iselect++;
        }
    } // ! end loop ipix

    verb_print(cmd->verbose,
        "\n\t%s: masked pixels = %ld\n",
               routineName, iselect);
    verb_print(cmd->verbose,
        "\t%s: unmasked pixels = %ld\n",
               routineName, gd->nbodyTable[gd->iCatalogs[0]]-iselect);
    
    free(map);

    return SUCCESS;
}


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

//    float *map;
    long npixel, nside;
    char order1[10];
    char order2[10];
    char coord[10];
    
    int hp_status;
    
    //B
    float *map = NULL;
    float *plot_map = NULL;
    bodyptr bodytabtmp = NULL;
    int bodytable_owned = FALSE;
    //E

#if THREEDIMCODE
    verb_print(cmd->verbose, "\nWorking 3D map...\n");
#else
    cBALLS_FAIL(cmd, "\nOnly 3D is implemented so far... exiting...\n\n");
#endif

    hp_status = get_fits_size_status(filename, &nside, order1, &npixel);
    if (hp_status != 0) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: get_fits_size failed for '%s' status=%d",
                 routineName, filename, hp_status);
        return FAILURE;
    }
    
    verb_print(cmd->verbose,
        "filename, ifile, nside, npixel, order1:");
    verb_print(cmd->verbose,
        "%s %d %ld %ld %s\n", filename, ifile, nside, npixel, order1);
    hp_status = read_healpix_map_status(filename, &nside, coord, order2, &map);
    if (hp_status != 0 || map == NULL) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: read_healpix_map failed for '%s' status=%d",
                 routineName, filename, hp_status);
        goto fail;
    }
    if (cfitsio_healpix_map_to_ring(cmd, routineName, nside,
                                    order1, order2, &map) == FAILURE)
        goto fail;
    //B
    if (npixel < 1) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: npixel = %ld is absurd", routineName, npixel);
        goto fail;
    }
    //E
    verb_print(cmd->verbose,
               "\nAllocated %g MByte for temporal map pixel (%ld) storage.\n",
               npixel*sizeof(float)*INMB,
               npixel);

    verb_print(cmd->verbose,
               "%s: nbody = %ld...\n",
               routineName, npixel);

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
    map = NULL;
    verb_print(cmd->verbose,
               "\nFreed %g MByte for temporal map pixel (%ld) storage.\n",
               npixel*sizeof(float)/(1024.0*1024.0),
               npixel);


    bodyptr q;
    cmd->nbody = iselect;
    gd->nbodyTable[ifile] = cmd->nbody;
    //B
    if (iselect < 1) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: mask-inside selected no pixels", routineName);
        goto fail;
    }
    //E
    bodytable[ifile] = (bodyptr) allocate(cmd->nbody * sizeof(body));
    gd->bodytable_allocated = TRUE;
    bodytable_owned = TRUE;

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
        //B
        // leading ! to allow overwrite
        if (format_checked(fileforce, sizeof(fileforce),
                           "fileforce", "!%s/%s",cmd->rootDir, file) != 0) {
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     "%s: output FITS path too long", routineName);
            goto fail;
        }
        //E
        verb_print(cmd->verbose,
                   "\t%s: %s %s...\n",
                   routineName, "\n\t\tsaving map to a fits file:", fileforce);
        //B
        plot_map = (float *) malloc(npixel * sizeof(float));
        if (plot_map == NULL) {
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     "%s: cannot allocate plot map", routineName);
            goto fail;
        }
        //E
        for(ipix=0;ipix<npixel;ipix++){
            q = bodytabtmp+ipix;
            if(Update(q)) {
                plot_map[ipix] = Kappa(q);
            } else {
                plot_map[ipix] = 0.0;
            }
        }
        //B write healpix map in RING order "0"
        //      and
        //      "G = Galactic, E = ecliptic, C = celestial = equatorial"
        hp_status = write_healpix_map_status(plot_map, nside, fileforce, 0, "C");
        //B
        if (hp_status != 0) {
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     "%s: write_healpix_map failed for '%s' status=%d",
                     routineName, fileforce, hp_status);
            goto fail;
        }
        //E
        
        fprintf(stdout,"\t\tfile written\n");
        free(plot_map);
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

fail:
    if (plot_map != NULL)
        free(plot_map);
    if (map != NULL)
        free(map);
    if (bodytabtmp != NULL)
        free(bodytabtmp);
    if (bodytable_owned && bodytable[ifile] != NULL) {
        free(bodytable[ifile]);
        bodytable[ifile] = NULL;
        gd->nbodyTable[ifile] = 0;
        cmd->nbody = 0;
    }
    return FAILURE;
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

    int rc = SUCCESS;

    gd->input_comment = "numpy-healpix input file";

    verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                          "\t%s: opening numpy-healpix file: %s...\n",
                          routineName,filename);

    if (scanopt(cmd->options, "header-info")){
        verb_print(cmd->verbose,"\nHeader information:");
        if (format_checked(buf, sizeof(buf),
            "buf", "%s",headerScript) != 0)
            return FAILURE;
        verb_print(cmd->verbose,
                   "\n\t%s: header processing: executing %s...\n",
                   routineName, headerScript);
        if (cballs_system_checked(cmd, routineName, buf) == FAILURE)
            return FAILURE;

        verb_print(cmd->verbose, "\tdone.\n");

        if (scanopt(cmd->options, "stop-numpy")) {
            if (strnull(cmd->outfile)) {
                gd->inputHeaderFlag = TRUE;
                gd->stopflag = TRUE;
                return SUCCESS;
            }
        }
    }

    if (scanopt(cmd->options, "read-mask")) {
        if (ifile == 0) {
            rc = inputdata_numpy_healpix_map(cmd, gd, filename, ifile);
        } else {
            if (ifile != 1) {
                cBALLS_FAIL(cmd, "\t%s: read-mask ifile = %d is absurd\n",
                            routineName, ifile+1);
            }

            rc = inputdata_numpy_healpix_map_mask(cmd, gd, filename, ifile);
        }
    } else {
        if (scanopt(cmd->options, "mask-inside")) {
            rc = inputdata_numpy_healpix_map_mask_inside(cmd, gd,
                                                         filename, ifile);
        } else {
            rc = inputdata_numpy_healpix_map(cmd, gd, filename, ifile);
        }
    }

    if (rc == FAILURE) {
        return FAILURE;
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

    double *map;
    long npixel, nside;
    string order = "RING";
    FILE *fp;
    
    int hp_status;

#if THREEDIMCODE
    verb_print(cmd->verbose, "\nWorking 3D map...\n");
#else
    cBALLS_FAIL(cmd, "\nOnly 3D is implemented so far... exiting...\n\n");
#endif

    // info runing header_script.py and set it in cmd->nbody
    npixel = cmd->nbody;            // must be in the form 12*nside^2
    nside = npix2nside(npixel);
    verb_print(cmd->verbose,
        "filename, ifile, nside, npixel, order:");
    verb_print(cmd->verbose,
        "%s %d %ld %ld %s\n", filename, ifile, nside, npixel, order);

    fp = stropen(filename, "rb");       // open numpy file
    map=(double *)malloc(sizeof(double)*npixel);
    verb_print(cmd->verbose,
               "\nAllocated %g MByte for temporal map pixel (%ld) storage.\n",
               npixel*sizeof(double)*INMB,
               npixel);
    for(j=0;j<npixel;j++){
        fread(&map[j], sizeof(double), 1, fp);
    }
    fclose(fp);             // close numpy file.

    verb_print(cmd->verbose,
               "%s: nbody = %ld...\n",
               routineName, npixel);
    if (npixel < 1)
        cBALLS_FAIL(cmd, "%s: npixel = %ld is absurd\n", routineName, npixel);

    bodyptr bodytabtmp;
    bodytabtmp = (bodyptr) allocate(npixel * sizeof(body));
    verb_print(cmd->verbose,
               "\nAllocated %g MByte for temporal particle (%ld) storage.\n",
               npixel*sizeof(body)*INMB,
               npixel);

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
    gd->bodytable_allocated = TRUE;

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
#if defined(NMultipoles) && defined(NONORMHIST)
            UpdatePivot(p) = UpdatePivot(q);
#endif
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
        // leading ! to allow overwrite
        if (format_checked(fileforce, sizeof(fileforce),
            "fileforce", "!%s/%s",cmd->rootDir, file) != 0)
            return FAILURE;
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
        hp_status = write_healpix_map_status(mapout, nside, fileforce, 0, "C");
        if (hp_status != 0) {
            free(mapout);
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     "%s: write_healpix_map failed for '%s' status=%d",
                     routineName, fileforce, hp_status);
            return FAILURE;
        }
        
        fprintf(stdout,"\t\tfile written\n");
        free(mapout);
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
    string routineName = "inputdata_cfitsio_healpix_map_mask";
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
    
    int hp_status;

#if THREEDIMCODE
    verb_print(cmd->verbose, "\nWorking 3D map...\n");
#else
    cBALLS_FAIL(cmd, "\nOnly 3D is implemented so far... exiting...\n\n");
#endif

    hp_status = get_fits_size_status(filename, &nside, order1, &npixel);
    if (hp_status != 0) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: get_fits_size failed for '%s' status=%d",
                 routineName, filename, hp_status);
        return FAILURE;
    }
    
    verb_print(cmd->verbose,
        "filename, ifile, nside, npixel, order1:");
    verb_print(cmd->verbose,
        "%s %d %ld %ld %s\n", filename, ifile, nside, npixel, order1);
    hp_status = read_healpix_map_status(filename, &nside, coord, order2, &map);
    if (hp_status != 0 || map == NULL) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: read_healpix_map failed for '%s' status=%d",
                 routineName, filename, hp_status);
        return FAILURE;
    }
    if (cfitsio_healpix_map_to_ring(cmd, routineName, nside,
                                    order1, order2, &map) == FAILURE) {
        free(map);
        return FAILURE;
    }
    verb_print(cmd->verbose,
               "\nAllocated %g MByte for temporal map pixel (%ld) storage.\n",
               npixel*sizeof(float)*INMB,
               npixel);

    verb_print(cmd->verbose, "%s: nbody = %ld...\n",
               routineName, npixel);
    if (npixel < 1) {
        free(map);
        cBALLS_FAIL(cmd, "%s: npixel = %ld is absurd\n", routineName, npixel);
    }

    if (npixel != gd->nbodyTable[gd->iCatalogs[0]]) {
        free(map);
        cBALLS_FAIL(cmd,
            "%s: npixel = %ld is not equal to npixel in cat 0: %ld\n",
            routineName, npixel, gd->nbodyTable[gd->iCatalogs[0]]);
    }

    INTEGER iselect = 0;
    for(ipix=0;ipix<npixel;ipix++) {                // RING loop order
        p = bodytable[gd->iCatalogs[0]]+ipix;
        pix2ang_ring(nside, ipix, &theta, &phi);
        coordinate_transformation(cmd, gd, theta, phi, q);
        DO_COORD(k) {
            if (Pos(p)[k] != q[k]) {
                free(map);
                cBALLS_FAIL(cmd, "%s: mask position is not equal: %g %g\n",
                            routineName, Pos(p)[k], q[k]);
            }
        }

        Mask(p) = map[ipix];
        if (Mask(p) == 0) {
            iselect++;
        }
    } // ! end loop ipix

    verb_print(cmd->verbose,
        "\n\t%s: masked pixels = %ld\n",
               routineName, iselect);
    verb_print(cmd->verbose,
        "\t%s: unmasked pixels = %ld\n",
               routineName, gd->nbodyTable[gd->iCatalogs[0]]-iselect);
    
    free(map);

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

    float *map = NULL;
    float *plot_map = NULL;
    bodyptr bodytabtmp = NULL;
    bodyptr selected_bodies = NULL;
    long npixel, nside;
    char order1[10];
    char order2[10];
    char coord[10];
    
    int hp_status;
    int status = FAILURE;

#if THREEDIMCODE
    verb_print(cmd->verbose, "\nWorking 3D map...\n");
#else
    cBALLS_FAIL(cmd, "\nOnly 3D is implemented so far... exiting...\n\n");
#endif

    hp_status = get_fits_size_status(filename, &nside, order1, &npixel);
    if (hp_status != 0) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: get_fits_size failed for '%s' status=%d",
                 routineName, filename, hp_status);
        goto cleanup;
    }
    
    verb_print(cmd->verbose,
        "filename, ifile, nside, npixel, order1:");
    verb_print(cmd->verbose,
        "%s %d %ld %ld %s\n", filename, ifile, nside, npixel, order1);
    hp_status = read_healpix_map_status(filename, &nside, coord, order2, &map);
    if (hp_status != 0 || map == NULL) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: read_healpix_map failed for '%s' status=%d",
                 routineName, filename, hp_status);
        goto cleanup;
    }
    if (cfitsio_healpix_map_to_ring(cmd, routineName, nside,
                                    order1, order2, &map) == FAILURE)
        goto cleanup;
    verb_print(cmd->verbose,
               "\nAllocated %g MByte for temporal map pixel (%ld) storage.\n",
               npixel*sizeof(float)*INMB,
               npixel);

    verb_print(cmd->verbose,
               "%s: nbody = %ld...\n",
               routineName, npixel);
    if (npixel < 1) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: npixel = %ld is absurd", routineName, npixel);
        goto cleanup;
    }

    bodytabtmp = (bodyptr) allocate(npixel * sizeof(body));
    verb_print(cmd->verbose,
               "\nAllocated %g MByte for temporal particle (%ld) storage.\n",
               npixel*sizeof(body)*INMB,
               npixel);

    real xmin, ymin, zmin;
    real xmax, ymax, zmax;
    xmin=0., ymin=0., zmin=0.;
    xmax=0., ymax=0., zmax=0.;

    INTEGER i;
    INTEGER iselect = 0;
    for(ipix=0;ipix<npixel;ipix++) {                // RING loop order
        p = bodytabtmp+ipix;
        Update(p) = FALSE;
        Mask(p) = MASK_NODE_MASKED;
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
                    Mask(p) = MASK_NODE_VALID;
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
               routineName, npixel-iselect);
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

    bodyptr q;
    if (iselect < 1) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: no pixels selected from '%s'", routineName, filename);
        goto cleanup;
    }
    selected_bodies = (bodyptr) allocate(iselect * sizeof(body));

    real kavg = 0;
    INTEGER ij=0;
    for(ipix=0;ipix<npixel;ipix++){
        q = bodytabtmp+ipix;
        if(Update(q)) {
            p = selected_bodies+ij;
            Pos(p)[0] = Pos(q)[0];
            Pos(p)[1] = Pos(q)[1];
            Pos(p)[2] = Pos(q)[2];
            Kappa(p) = Kappa(q);
            Type(p) = Type(q);
            Mass(p) = mass;
            Weight(p) = weight;
            Id(p) = ij + 1;
            Mask(p) = Mask(q);
            Update(p) = TRUE;
#if defined(NMultipoles) && defined(NONORMHIST)
            if (scanopt(cmd->options, "patch-with-all"))
                UpdatePivot(p) = UpdatePivot(q);
#endif
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
        // leading ! to allow overwrite
        if (format_checked(fileforce, sizeof(fileforce),
            "fileforce", "!%s/%s",cmd->rootDir, file) != 0)
            goto cleanup;
        verb_print(cmd->verbose,
                   "\t%s: %s %s...\n",
                   routineName, "\n\t\tsaving map to a fits file:", fileforce);
        plot_map = (float *)malloc(npixel*sizeof(float));
        if (plot_map == NULL) {
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     "%s: unable to allocate output map", routineName);
            goto cleanup;
        }
        for(ipix=0;ipix<npixel;ipix++){
            q = bodytabtmp+ipix;
            if(Update(q)) {
                plot_map[ipix] = Kappa(q);
            } else {
                plot_map[ipix] = 0.0;
            }
        }
        //B write healpix map in RING order "0"
        //      and
        //      "G = Galactic, E = ecliptic, C = celestial = equatorial"
        hp_status = write_healpix_map_status(plot_map, nside,
                                             fileforce, 0, "C");
        if (hp_status != 0) {
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     "%s: write_healpix_map failed for '%s' status=%d",
                     routineName, fileforce, hp_status);
            goto cleanup;
        }
        
        fprintf(stdout,"\t\tfile written\n");
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

    if (scanopt(cmd->options, "all"))
        verb_print(cmd->verbose,
            "\n\t%s: selected read points and nbody: %ld %ld\n",
                   routineName, iselect, cmd->nbody);
    else
        verb_print(cmd->verbose,
            "\n\t%s: selected read points = %ld\n", routineName, iselect);

    verb_print(cmd->verbose,
            "%s: average of kappa (%ld particles) = %le\n",
               routineName, iselect, kavg/((real)iselect) );

    cmd->nbody = iselect;
    gd->nbodyTable[ifile] = iselect;
    bodytable[ifile] = selected_bodies;
    gd->bodytable_allocated = TRUE;
    selected_bodies = NULL;
    status = SUCCESS;

cleanup:
    free(plot_map);
    free(map);
    free(bodytabtmp);
    free(selected_bodies);
    return status;
}

#endif	// ! _cballsio_cfitsio_02_h
