
//=============================================================================
//        1          2          3          4        ^ 5          6          7

// included in: addons/addons_include/source/cballsio/cballsio_indlude_11b.h

#ifndef _cballsio_cfitsio_07_h
#define _cballsio_cfitsio_07_h

#include "fitsio.h"


local int outputdata_cfitsio(struct cmdline_data* cmd, struct  global_data* gd,
                             bodyptr bodytab, INTEGER nbody)
{
    char namebuf[256];

    gd->output_comment = "fits output file";

//    sprintf(namebuf, "%s", gd->fpfnameOutputFileName);
    if (format_checked(namebuf, sizeof(namebuf),
        "namebuf", "%s", gd->fpfnameOutputFileName) != 0)
        return FAILURE;
    verb_print(cmd->verbose,
               "\toutputdata_cfitsio: opening fits file: %s...\n",
               namebuf);

    if (scanopt(cmd->options, "kappa")) {
        return writebintable_kappa(cmd, gd, bodytab, nbody, namebuf);
    } else { // ! kappa
        verb_print(cmd->verbose,
        "\toutputdata_cfitsio: only convergence (kappa) is %s\n\n",
                   "implemented. No file is saved...");
    }

    return SUCCESS;
}

static int fits_failure(struct cmdline_data *cmd, int status, const char *where);

// write a fits file for convergence (kappa)
//  only 3D
int writebintable_kappa(struct cmdline_data* cmd,
                         struct  global_data* gd,
                         bodyptr bodytab, INTEGER nbody,
                         string filename)
{
    bodyptr p;
    int hdutype;
    long firstrow, firstelem;
    long nrows;

    int tfields   = 4;

    char extname[] = "Kappa_Binary";            // extension name

    //B define the name, datatype, and physical units for the 4 columns
    char *ttype[] = { "X", "Y", "Z", "KAPPA" };
    char *tform[] = { "1D","1D","1D","1D"    };
    char *tunit[] = { " ",  " ",  " ",  " "      };
    //E

    fitsfile *fptr = NULL;
    int status = 0;
    int rc = SUCCESS;
    double *arrayX = NULL;
    double *arrayY = NULL;
    double *arrayZ = NULL;
    double *kappa = NULL;

    status=0;

    nrows = nbody;
    arrayX = (double*) allocate(nbody * sizeof(double));
    arrayY = (double*) allocate(nbody * sizeof(double));
    arrayZ = (double*) allocate(nbody * sizeof(double));
    kappa = (double*) allocate(nbody * sizeof(double));
    DO_BODY(p, bodytab, bodytab+nbody) {
        arrayX[p-bodytab] = Pos(p)[0];
        arrayY[p-bodytab] = Pos(p)[1];
        arrayZ[p-bodytab] = Pos(p)[2];
        kappa[p-bodytab] = Kappa(p);
    }

    long naxes[] = {0,0};

    if (fits_create_file(&fptr, filename, &status)) {
        rc = fits_failure(cmd, status, "fits_create_file");
        goto cleanup;
    }

    fits_create_img(fptr, SHORT_IMG, 0, naxes, &status);

    if (fits_create_tbl(fptr, BINARY_TBL, nrows, tfields, ttype, tform,
                        tunit, extname, &status)) {
        rc = fits_failure(cmd, status, "fits_create_tbl");
        goto cleanup;
    }
    

    firstrow  = 1;
    firstelem = 1;

    fits_write_col(fptr, TDOUBLE, 1, firstrow, firstelem, nrows, arrayX,
                   &status);
    fits_write_col(fptr, TDOUBLE, 2, firstrow, firstelem, nrows, arrayY,
                   &status);
    fits_write_col(fptr, TDOUBLE, 3, firstrow, firstelem, nrows, arrayZ,
                   &status);
    fits_write_col(fptr, TDOUBLE, 4, firstrow, firstelem, nrows, kappa,
                   &status);
    
    if (fits_close_file(fptr, &status)) {
        rc = fits_failure(cmd, status, "fits_close_file");
        fptr = NULL;
        goto cleanup;
    }
    fptr = NULL;
    

cleanup:
    if (fptr != NULL) {
        int close_status = 0;
        if (fits_close_file(fptr, &close_status) && rc == SUCCESS)
            rc = fits_failure(cmd, close_status, "fits_close_file");
    }

    free(kappa);
    free(arrayZ);
    free(arrayY);
    free(arrayX);

    return rc;
}


//B printerror() above is replaced by this routine
/*
static int fits_failure(struct cmdline_data *cmd, int status, const char *where)
{
    if (!status) return SUCCESS;

    fits_report_error(stderr, status);
#ifdef CLASSLIB
    snprintf(cmd->error_message, _ERRORMSGSIZE_,
             "%s: CFITSIO status=%d", where, status);
    return FAILURE;
#else
    error("%s: CFITSIO status=%d\n", where, status);
    return FAILURE;
#endif
}
*/
static int fits_failure(struct cmdline_data *cmd, int status, const char *where)
{
    if (!status) return SUCCESS;

    fits_report_error(stderr, status);

    if (cmd != NULL) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: CFITSIO status=%d", where, status);
    } else {
        fprintf(stderr, "%s: CFITSIO status=%d\n", where, status);
    }

    return FAILURE;
}
//E

local int outputdata_numpy_healpix(struct cmdline_data* cmd, struct  global_data* gd,
                             bodyptr bodytab, INTEGER nbody)
{
    char namebuf[256];

    gd->output_comment = "numpy-healpix output file";

//    sprintf(namebuf, "%s", gd->fpfnameOutputFileName);
    if (format_checked(namebuf, sizeof(namebuf),
        "namebuf", "%s", gd->fpfnameOutputFileName) != 0)
        return FAILURE;
    verb_print(cmd->verbose,
               "\toutputdata_cfitsio: opening fits file: %s...\n",
               namebuf);

    if (scanopt(cmd->options, "kappa")) {
        return writebintable_kappa(cmd, gd, bodytab, nbody, namebuf);
    } else { // ! kappa
        verb_print(cmd->verbose,
        "\toutputdata_cfitsio: only convergence (kappa) is %s\n\n",
                   "implemented. No file is saved...");
    }

    return SUCCESS;
}


#endif	// ! _cballsio_cfitsio_07_h
