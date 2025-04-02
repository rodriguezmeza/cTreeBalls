/*==============================================================================
 NAME: main.c				[cTreeBalls]
 Written by: M.A. Rodriguez-Meza
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
 Comments and notes: (palimsesto)... coding based on references below.
 References:    Zeno project, NEMO project, Gadget, COLA, CLASS,
                NR, GSL, Rapaport's book, cute, cfitsio...
 github: https://github.com/rodriguezmeza/cTreeBalls.git
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

    gd.cpuinit = CPUTIME;                           // init register of cpu time
    gd.cpurealinit = rcpu_time();                   // init register of real time

#ifdef GETPARAM
    InitParam(argv, defv);                          // init command parameters
                                                    //  structure
#else
    if(argc < 2) {
        verb_print(1, "Parameters are missing.\n");
        verb_print(1, "Call with <ParameterFile>\n");
        endrun_mpi(ThisTask, 0);
    }
    strcpy(cmd.ParameterFile, argv[1]);
    printf("\n -> Parameter file is %s\n",
           cmd.ParameterFile);
#endif

    StartRun(&cmd, &gd, argv[0],                    // get parameters and
             HEAD1, HEAD2, HEAD3);                  //  init global structure
                                                    //  and do other useful 
                                                    //  process, like param check,
                                                    //  read data points to analyze
	MainLoop(&cmd, &gd);                            // make tree and search data
	EndRun(&cmd, &gd);                              // close streams and free mem

    return SUCCESS;
}

