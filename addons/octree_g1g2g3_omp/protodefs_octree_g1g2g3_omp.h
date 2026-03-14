// Use:
//#include "protodefs_octree_g1g2g3_omp.h"

//
// included in file: addons/addons_include/include/protodefs_include.h
//

#ifndef _protodefs_octree_g1g2g3_omp_h
#define _protodefs_octree_g1g2g3_omp_h

global int searchcalc_octree_g1g2g3_omp(struct cmdline_data* cmd,
                                     struct  global_data* gd,
                                     bodyptr *btable, INTEGER *nbody,
                                     INTEGER ipmin, INTEGER *ipmax,
                                     int cat1, int cat2, int cat3);

#endif	// ! _protodefs_octree_g1g2g3_omp_h
