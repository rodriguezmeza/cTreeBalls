/*==============================================================================
 NAME: main.c				[tpcf]
 Written by: M.A. Rodriguez-Meza and Alejandro Aviles
 Starting date: april 2023
 Purpose: Main routine
 Language: C
 Major revision:
 ================================================================================
 
 Use: tpcf -help
 Input: 	Command line parameters, Parameters file
 Output: ...
 Units:
 History:
 Comments and notes:
 References:    Barnes' Treecode, NEMO project, Gadget, COLA, CLASS,
                NR, GSL, Rapaport's book
 ==============================================================================*/

#define global

#include "globaldefs.h"
#include "cmdline_defs.h"

int main(int argc, string argv[])
{
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

    StartRun(argv[0], HEAD1, HEAD2, HEAD3);
	MainLoop();

	EndRun();

#ifdef MPICODE
    MPI_Finalize();
#endif

    return 0;
}

