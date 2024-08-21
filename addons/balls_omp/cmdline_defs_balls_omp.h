// Use:
//#include "cmdline_defs_balls_omp.h"

#ifndef _cmdline_defs_balls_omp_h
#define _cmdline_defs_balls_omp_h

//
/*
#ifndef BALLSSINCOSEXEC
#ifdef BALLSEXEC
    "searchMethod=balls-omp",             ";Searching method to use", ":search",
#else
    "searchMethod=tree-omp",             ";Searching method to use", ":search",
#endif
#else
    "searchMethod=balls-omp",             ";Searching method to use", ":search",
#endif
*/
//
    "scanLevel=6",                     ";Scan level to start the search (look at tdepth value, will be the maximum for this parameter)", ":scl",
// Root nodes:
    "scanLevelRoot=0",                 ";Scan level of root cells to start the search (look at tdepth value, will be the maximum for this parameter)", ":sclroot",
    "scanLevelMin=-0",                 ";Scan level of size cells to stop the search. Integer negative values (look at tdepth value, will be tdepth-1+scanLevelMin+1)", ":sclmin",
    "ntosave=1000",                     ";Number of found bodies to save; use in combination with 'bodyfound', balls4' method", ":ntsav",
//    "rsmooth=",                         ";Radius of the pivot smoothing neighbourhood", ":rsm",


#endif	// ! _cmdline_defs_balls_omp_h
