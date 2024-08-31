// Use:
//#include "startrun_iolib_06.h"

#ifndef _startrun_iolib_06_h
#define _startrun_iolib_06_h

scaniOption(cmd, gd, cmd->columns, gd->columns,
            &nitems, ndummy, 2, "columns");

// pos, kappa, gamma1, gamma2, weight, seven places maximum
// use in multi-columns-ascii and cfitsio
// default is 3D: three pos and one convergence (kappa)

//B Adapt to have pos: first three positions
//              kappa: fourth position
//              gamma1: fifth position
//              gamma1: sixth position
if (nitems==1) {
    gd->columns[1] = 2;
    gd->columns[2] = 3;
    gd->columns[3] = 4;
    gd->columns[4] = 5;
    gd->columns[5] = 6;
    gd->columns[6] = 7;
}

if (nitems==2) {
    gd->columns[2] = 3;
    gd->columns[3] = 4;
    gd->columns[4] = 5;
    gd->columns[5] = 6;
    gd->columns[6] = 7;
}

if (nitems==3) {
    gd->columns[3] = 4;
    gd->columns[4] = 5;
    gd->columns[5] = 6;
    gd->columns[6] = 7;
}

if (nitems==4) {
    gd->columns[4] = 5;
    gd->columns[5] = 6;
    gd->columns[6] = 7;
}

int ii;
if (cmd->verbose_log>=3)
    for (ii=0; ii< nitems; ii++)
        verb_log_print(cmd->verbose_log, gd->outlog,
                   "\tscaniOptions: columns[%d]: %d...\n\n",
                   ii, gd->columns[ii]);


#endif	// ! _startrun_iolib_06_h
