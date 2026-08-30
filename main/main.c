/*==============================================================================
 NAME: main.c				[cTreeBalls]
 Written by: Mario A. Rodriguez-Meza
 Starting date: april 2023
 Purpose: Main routine
 Language: C
 Major revision:
 ===============================================================================
 Use: cballs --help (or -h)
 Input: 	Command line parameters, Parameters file, data catalogs
 Output: several histograms containing 2pcf, 3pcf,...
 Units:
 History:
 Acknowledgements: Abraham Arvizu, Alejandro Aviles, Juan Carlos Hidalgo,
                    Eladio Moreno, Gustavo Niz, Sadi Ramirez,
                    Axel Romero-Tisnado and Sofia Samario
 Comments and notes: (palimsesto)... coding based on references below.
 References:    Zeno project, NEMO project, Gadget, COLA, CLASS,
                NR, GSL, Rapaport's book, cute, cfitsio, Healpix...
 github: https://github.com/rodriguezmeza/cTreeBalls.git
 Publication: cite: JCAP12(2024)049 (ArXiv ePrint: 2408.16847)
 ==============================================================================*/
//        1          2          3          4        ^ 5          6          7

#define global

#include <stdio.h>
#include "globaldefs.h"
#include "cmdline_defs.h"

/*
 Main routine:
 
 This routine is in charge of main computational flow
    as it is explained below in the comments.

 Arguments:
    * `argc`: Input: int
    * `argv`: Input: string array
 Return (the error status):
    int SUCCESS or FAILURE
 */
int main(int argc, string argv[])
{
    struct cmdline_data cmd;                        // share command parameters
    struct global_data gd;                          // share global parameters
    int run_status = FAILURE;

    memset(&cmd, 0, sizeof(cmd));
    memset(&gd, 0, sizeof(gd));

    gd.cpuinit = CPUTIME;                           // init register of cpu time
    gd.cpurealinit = rcpu_time();                   // init register of real time

#ifdef GETPARAM
    InitParam(argv, defv);                          // init command parameters
                                                    //  structure
#else
    if(argc < 2) {
        verb_print(1, "Parameters are missing.\n");
        verb_print(1, "Call with <ParameterFile>\n");
        goto finish;
    }
    if (format_checked(cmd.ParameterFile, sizeof(cmd.ParameterFile),
        "cmd.ParameterFile", "%s", argv[1]) != 0)
        goto finish;
    printf("\n -> Parameter file is %s\n",
           cmd.ParameterFile);
#endif

    //B get parameters, init global structure, and do other useful
    //      process, like param check, read data points to analyze
    if (StartRun(&cmd, &gd, argv[0], HEAD1, HEAD2, HEAD3) == FAILURE) {
        if (cmd.error_message[0] != '\0')
            fprintf(stderr, "\nError in cballs: %s\n", cmd.error_message);
        EndRun_FreeMemory(&cmd, &gd);
        goto finish;
    }
    //E

    if (MainLoop(&cmd, &gd) == FAILURE) {           // make tree and search data
        if (cmd.error_message[0] != '\0')
            fprintf(stderr, "\nError in cballs: %s\n", cmd.error_message);
        EndRun(&cmd, &gd);                          // close streams and free mem
        goto finish;
    }

    if (EndRun(&cmd, &gd) == FAILURE)               // close streams and free mem
        goto finish;

    run_status = SUCCESS;

finish:
#ifdef CBALLS_MPI_ENABLED
    if (cballs_mpi_finalize(&cmd) == FAILURE) {
        if (cmd.error_message[0] != '\0')
            fprintf(stderr, "\nError finalizing MPI: %s\n",
                    cmd.error_message);
        run_status = FAILURE;
    }
#endif
    return run_status;
}

