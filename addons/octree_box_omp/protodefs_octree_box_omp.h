// Use:
//#include "protodefs_octree_box_omp.h"

// included in addons_include/include/protodefs_include.h

#ifndef _protodefs_octree_box_omp_h
#define _protodefs_octree_box_omp_h

global int searchcalc_octree_box_omp(struct cmdline_data* cmd,
                                 struct  global_data* gd,
                                 bodyptr *btable,
                                 INTEGER *nbody,
                                 INTEGER ipmin, INTEGER *ipmax,
                                 int cat1, int cat2);

#endif	// ! _protodefs_octree_box_omp_h
