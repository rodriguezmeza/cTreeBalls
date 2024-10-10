// Use:
//#include "protodefs_octree_kkk_omp.h"

#ifndef _protodefs_octree_kkk_omp_h
#define _protodefs_octree_kkk_omp_h

global int searchcalc_octree_kkk_omp(struct cmdline_data* cmd,
                                     struct  global_data* gd,
                                     bodyptr *btable, INTEGER *nbody,
                                     INTEGER ipmin, INTEGER *ipmax,
                                     int cat1, int cat2);

#endif	// ! _protodefs_octree_kkk_omp_h
