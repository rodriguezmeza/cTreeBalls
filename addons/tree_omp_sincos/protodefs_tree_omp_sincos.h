// Use:
//#include "protodefs_tree_omp_sincos.h"

// included in addons_include/include/protodefs_include.h

#ifndef _protodefs_tree_omp_sincos_h
#define _protodefs_tree_omp_sincos_h

global int searchcalc_tree_omp_sincos(struct cmdline_data* cmd,
                                    struct  global_data* gd,
                                    bodyptr *btable,
                                    INTEGER *nbody,
                                    INTEGER ipmin, INTEGER *ipmax,
                                    int cat1, int cat2);

#endif	// ! _protodefs_tree_omp_sincos_h
