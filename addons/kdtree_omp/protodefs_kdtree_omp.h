// Use:
//#include "protodefs_kdtree_omp.h"

#ifndef _protodefs_kdtree_omp_h
#define _protodefs_kdtree_omp_h

global int searchcalc_kdtree_omp(struct cmdline_data* cmd,
                                 struct  global_data* gd,
                                 bodyptr *btable,
                                 INTEGER *nbody, INTEGER ipmin, INTEGER *ipmax,
                                 int cat1, int cat2);

#ifdef NMultipoles
global int search_init_sincos_omp_N(struct cmdline_data* cmd,
                                  struct  global_data* gd,
                                  gdhistNptr_sincos_omp hist);
global int search_free_sincos_omp_N(struct cmdline_data* cmd,
                                    struct  global_data* gd,
                                  gdhistNptr_sincos_omp hist);
global int computeBodyProperties_sincos_N(struct cmdline_data* cmd,
                                            struct  global_data* gd,
                                            bodyptr, int,
                                            gdhistNptr_sincos_omp);
global int search_init_gd_hist_sincos_N(struct  cmdline_data* cmd,
                                        struct  global_data* gd);
#endif

#endif	// ! _protodefs_kdtree_omp_h
