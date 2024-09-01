
//=============================================================================
//        1          2          3          4        ^ 5          6          7

#ifndef _cballsio_cfitsio_07_h
#define _cballsio_cfitsio_07_h

#include "fitsio.h"

local int outputdata_cfitsio(struct cmdline_data* cmd, struct  global_data* gd,
                             bodyptr bodytab, INTEGER nbody)
{
    fitsfile *fptr;
    char card[FLEN_CARD];
    int status = 0, nkeys, ii;                      // MUST initialize status

    gd->output_comment = "fits output file";

    //B Was CFITSIO compiled with the -D_REENTRANT flag?  1 = yes, 0 = no.
    verb_print(cmd->verbose,
               "\toutputdata_cfitsio: -D_REENTRANT flag: %d...\n\n",
               fits_is_reentrant());
    //E

    if (scanopt(cmd->options, "fits-type-file"))
        fits_open_file(&fptr, cmd->outfile, READONLY, &status);
    else
        fits_open_data(&fptr, cmd->outfile, READONLY, &status);
    if (status) {                               // print any error messages
        verb_print(cmd->verbose,
                   "\toutputdata_cfitsio: open status: %d...\n\n",
                   status);
        fits_report_error(stderr, status);
    }

    fits_get_hdrspace(fptr, &nkeys, NULL, &status);
    if (status) {                               // print any error messages
        verb_print(cmd->verbose,
                   "\toutputdata_cfitsio: get_hdrspace status: %d...\n\n",
                   status);
        fits_report_error(stderr, status);
    }
    fits_get_num_rows(fptr, &cmd->nbody, &status);
    if (status) {                               // print any error messages
        verb_print(cmd->verbose,
                   "\toutputdata_cfitsio: get_num_rows status: %d...\n\n",
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
                   "\toutputdata_cfitsio: nbody = %d... and number of columns = %d\n\n",
                   cmd->nbody, ncols);
        if (cmd->nbody < 1)
            error("toutputdata_cfitsio:: nbody = %d is absurd\n", cmd->nbody);

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
                   "\toutputdata_cfitsio: nbody = %d... and number of columns = %d\n\n",
                   cmd->nbody, ncols);
        if (cmd->nbody < 1)
            error("toutputdata_cfitsio:: nbody = %d is absurd\n", cmd->nbody);
    }

    if (scanopt(cmd->options, "stop")) {
        fits_close_file(fptr, &status);
        if (status) {                               // print any error messages
            verb_print(cmd->verbose,
                       "\toutputdata_cfitsio: status: %d...\n\n",
                       status);
            fits_report_error(stderr, status);
        }
        exit(1);
    }

    if (scanopt(cmd->options, "ra-dec")) {
//        outputdata_cfitsio_ra_dec(cmd, gd,
//                                 cmd->outfile, 1, fptr);
    } else { // ! ra-dec
//        outputdata_cfitsio_xyz(cmd, gd,
//                              cmd->outfile, 1, fptr);
    }

    return SUCCESS;
}


#endif	// ! _cballsio_cfitsio_07_h
