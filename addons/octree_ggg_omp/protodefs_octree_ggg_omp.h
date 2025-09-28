// Use:
//#include "protodefs_octree_ggg_omp.h"

#ifndef _protodefs_octree_ggg_omp_h
#define _protodefs_octree_gg_omp_h

global int searchcalc_octree_ggg_omp(struct cmdline_data* cmd,
                                     struct  global_data* gd,
                                     bodyptr *btable, INTEGER *nbody,
                                     INTEGER ipmin, INTEGER *ipmax,
                                     int cat1, int cat2);

#endif	// ! _protodefs_octree_ggg_omp_h
