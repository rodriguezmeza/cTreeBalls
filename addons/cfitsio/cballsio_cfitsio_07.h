
//=============================================================================
//        1          2          3          4        ^ 5          6          7

#ifndef _cballsio_cfitsio_07_h
#define _cballsio_cfitsio_07_h

#include "fitsio.h"


local int outputdata_cfitsio(struct cmdline_data* cmd, struct  global_data* gd,
                             bodyptr bodytab, INTEGER nbody)
{
    char namebuf[256];

    gd->output_comment = "fits output file";

    sprintf(namebuf, "%s", gd->fpfnameOutputFileName);
    verb_print(cmd->verbose,
               "\toutputdata_cfitsio: opening fits file: %s...\n",
               namebuf);

    if (scanopt(cmd->options, "kappa")) {
        writebintable_kappa(cmd, gd, bodytab, nbody, namebuf);
    } else { // ! kappa
        verb_print(cmd->verbose,
        "\toutputdata_cfitsio: only convergence (kappa) is %s\n\n",
                   "implemented. No file is saved...");
    }

    return SUCCESS;
}

// write a fits file for convergence (kappa)
//  only 3D
void writebintable_kappa(struct cmdline_data* cmd,
                         struct  global_data* gd,
                         bodyptr bodytab, INTEGER nbody,
                         string filename)
{
    bodyptr p;
    fitsfile *fptr;

    int status=0;
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

    status=0;

    nrows = nbody;
    double *arrayX;
    double *arrayY;
    double *arrayZ;
    double *kappa;
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

    if (fits_create_file(&fptr, filename, &status))
         printerror(status);

    fits_create_img(fptr, SHORT_IMG, 0, naxes, &status);

    if (fits_create_tbl(fptr, BINARY_TBL, nrows, tfields, ttype, tform,
                        tunit, extname, &status));
        printerror(status);

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
    if (fits_close_file(fptr, &status))
        printerror(status);

    free(kappa);
    free(arrayZ);
    free(arrayY);
    free(arrayX);

    return;
}

void printerror( int status)
{
    if (status) {
       fits_report_error(stderr, status);
       exit(status);
    }

    return;
}

#endif	// ! _cballsio_cfitsio_07_h
