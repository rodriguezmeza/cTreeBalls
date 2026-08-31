// Use:
//#include "protodefs_octree_balls4_omp.h"

#ifndef _protodefs_octree_balls4_omp_h
#define _protodefs_octree_balls4_omp_h

#define OCTREEBALLS4OMPMETHOD 69
global int balls4_edge_search(struct cmdline_data *, struct global_data *,
                              bodyptr *, INTEGER *, INTEGER, INTEGER *, int, int);

global int searchcalc_octree_balls4_omp(struct cmdline_data* cmd,
                                     struct  global_data* gd,
                                     bodyptr *btable, INTEGER *nbody,
                                     INTEGER ipmin, INTEGER *ipmax,
                                     int cat1, int cat2);

#endif	// ! _protodefs_octree_balls4_omp_h
