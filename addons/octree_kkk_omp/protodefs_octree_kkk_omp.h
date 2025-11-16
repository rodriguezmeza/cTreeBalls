// Use:
//#include "protodefs_octree_kkk_omp.h"

#ifndef _protodefs_octree_kkk_omp_h
#define _protodefs_octree_kkk_omp_h

#if THREEDIMCODE
global int searchcalc_octree_kkk_omp(struct cmdline_data* cmd,
                                     struct  global_data* gd,
                                     bodyptr *btable, INTEGER *nbody,
                                     INTEGER ipmin, INTEGER *ipmax,
                                     int cat1, int cat2);
global int searchcalc_octree_kk_omp(struct cmdline_data* cmd,
                                     struct  global_data* gd,
                                     bodyptr *btable, INTEGER *nbody,
                                     INTEGER ipmin, INTEGER *ipmax,
                                     int cat1, int cat2);
global int searchcalc_octree_nn_omp(struct cmdline_data* cmd,
                                     struct  global_data* gd,
                                     bodyptr *btable, INTEGER *nbody,
                                     INTEGER ipmin, INTEGER *ipmax,
                                     int cat1, int cat2);
global int nncorrelation_octree_nn_omp(struct cmdline_data* cmd,
                                       struct  global_data* gd,
                                       bodyptr *btable, INTEGER *nbody,
                                       INTEGER ipmin, INTEGER *ipmax,
                                       int cat1, int cat2);
global int searchcalc_octree_kk_balls4_omp(struct cmdline_data* cmd,
                                             struct  global_data* gd,
                                             bodyptr *btable, INTEGER *nbody,
                                             INTEGER ipmin, INTEGER *ipmax,
                                           int cat1, int cat2);
#endif // ! THREEDIMCODE


#endif	// ! _protodefs_octree_kkk_omp_h
