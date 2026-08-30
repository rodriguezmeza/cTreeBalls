// Use:
//#include "cballsutils_balls_omp_0357.h"

#ifndef _cballsutils_balls_omp_0357_h
#define _cballsutils_balls_omp_0357_h

//#define BALLSOMPMETHOD         79

if (gd->searchMethod_int == 79) {
#ifdef TREENODEALLBODIES
    verb_print(cmd->verbose, "scanning %g bodies per thread\n",(real)nbody/(real)nthreads);
#else
    verb_print(cmd->verbose, "scanning %g nodes per thread\n",(real)gd->nnodescanlevTable[cat1]/(real)nthreads);
#endif
}

#endif	// ! _cballsutils_balls_omp_0357_h
