/*==============================================================================
 NAME: main.c				[cTreeBalls]
 Written by: M.A. Rodriguez-Meza
 Starting date: april 2023
 Purpose: Main routine
 Language: C
 Major revision:
 ===============================================================================
 //       1          2          3          4          5          6          7

 Use: cballs --help
 Input: 	Command line parameters, Parameters file, data catalogs
 Output: several histograms containg 2pcf, 3pcf,...
 Units:
 History:
 Comments and notes: (palimsesto)... coding based on references below.
 References:    Barnes' Treecode, NEMO project, Gadget, COLA, CLASS,
                NR, GSL, Rapaport's book
 github: https://github.com/rodriguezmeza/cTreeBalls.git
 ==============================================================================*/

#define global

#include "globaldefs.h"
#include "cmdline_defs.h"

int main(int argc, string argv[])
{
    struct cmdline_data cmd;
    struct global_data gd;

#ifdef MPICODE
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
    MPI_Comm_size(MPI_COMM_WORLD, &NTask);
    MPI_Get_processor_name(processor_name, &name_len);

    for(PTask = 0; NTask > (1 << PTask); PTask++);
#endif

    gd.cpuinit = CPUTIME;
    gd.cpurealinit = rcpu_time();

#ifdef GETPARAM
    InitParam(argv, defv);
#else
    if(argc < 2) {
        verb_print(1, "Parameters are missing.\n");
        verb_print(1, "Call with <ParameterFile> [<RestartFlag>]\n");
        endrun_mpi(ThisTask, 0);
    }
    strcpy(cmd.ParameterFile, argv[1]);
    printf("\n -> Parameter file is %s\n", cmd.ParameterFile);
#endif

//ADDONS: Setting this line here makes mpirun -np 2 tpcf run in Cosma
//#include "tpcf_01.h"

    StartRun(&cmd, &gd, argv[0], HEAD1, HEAD2, HEAD3);
	MainLoop(&cmd, &gd);
	EndRun(&cmd, &gd);

#ifdef MPICODE
    MPI_Finalize();
#endif

    return 0;
}

