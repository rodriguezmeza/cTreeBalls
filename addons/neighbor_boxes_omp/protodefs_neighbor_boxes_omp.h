// Use:
//#include "protodefs_neighbor_boxes_omp.h"

// included in addons_include/include/protodefs_include.h

#ifndef _protodefs_neighbor_boxes_omp_h
#define _protodefs_neighbor_boxes_omp_h

global int searchcalc_neighbor_boxes_omp(struct cmdline_data* cmd,
                                 struct  global_data* gd,
                                 bodyptr *btable,
                                 INTEGER *nbody, INTEGER ipmin, INTEGER *ipmax,
                                 int cat1, int cat2);

#endif	// ! _protodefs_neighbor_boxes_omp_h
