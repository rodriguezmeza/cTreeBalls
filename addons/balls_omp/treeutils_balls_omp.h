// Use:
//ADDONS:
//#include "treeutils_balls_omp.h"

#ifndef _treeutils_balls_omp_h
#define _treeutils_balls_omp_h

//#define BALLSOMPMETHOD         46

if (gd->searchMethod_int == 46) {
#ifdef TREENODEALLBODIES
    verb_print(cmd->verbose, "scanning %g bodies per thread\n",(real)nbody/(real)nthreads);
#else
    verb_print(cmd->verbose, "scanning %g nodes per thread\n",(real)gd->nnodescanlevTable[cat1]/(real)nthreads);
#endif
}

#endif	// ! _treeutils_balls_omp_h
