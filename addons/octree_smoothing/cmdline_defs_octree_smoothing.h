// Use:
//#include "cmdline_octree_smoothing.h"

#ifndef _cmdline_octree_smoothing_h
#define _cmdline_octree_smoothing_h

//B
    "scanLevel=6",              ";Scan level to start the search (look at tdepth value, will be the maximum for this parameter)", ":scl",
    "scanLevelRoot=0",          ";Scan level of root cells to start the search (look at tdepth value, will be the maximum for this parameter)", ":sclroot",
    "scanLevelMin=0",           ";Scan level of size cells to stop the search. Integer negative values (look at tdepth value, will be tdepth-1+scanLevelMin+1)", ":sclmin",
//E

#endif	// ! _cmdline_octree_smoothing_h
