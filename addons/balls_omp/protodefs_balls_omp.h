// Use:
//#include "protodefs_balls_omp.h"

#ifndef _protodefs_balls_omp_h
#define _protodefs_balls_omp_h

#ifdef BALLS
global void searchcalc_balls_omp(struct cmdline_data* cmd, struct  global_data* gd,
                                 bodyptr *, INTEGER *, INTEGER, INTEGER *, int, int);
global int search_init_balls_omp(struct cmdline_data* cmd, struct  global_data* gd,
                                 gdhistptr_omp_balls hist, int);
global int search_init_balls_omp_cc(struct cmdline_data* cmd,
                                    struct  global_data* gd,
                                    gdhistptr_omp_balls hist);
global int computeBodyProperties_balls_omp(struct cmdline_data* cmd, 
                                           struct  global_data* gd,
                                           bodyptr, int, gdhistptr_omp_balls);
global int computeBodyProperties_balls_omp_cc(struct cmdline_data* cmd, 
                                              struct  global_data* gd,
                                              bodyptr, int, gdhistptr_omp_balls);
global int computeBodyProperties_balls_omp_cc_sincos(struct cmdline_data* cmd, 
                                                     struct  global_data* gd,
                                                     bodyptr, int, gdhistptr_sincos_omp);
global int search_free_balls_omp(struct cmdline_data* cmd,
                                 struct  global_data* gd,
                                 gdhistptr_omp_balls hist);
global int search_free_balls_omp_sincos(struct cmdline_data* cmd, 
                                        struct  global_data* gd,
                                        gdhistptr_omp_balls hist);
global int search_free_balls_omp_cc(struct cmdline_data* cmd,
                                    struct  global_data* gd,
                                    gdhistptr_omp_balls hist);
global bool nodes_condition(struct cmdline_data* cmd, struct  global_data* gd,
                            nodeptr p, nodeptr q, real *dr1, vector dr);
global bool nodes_condition5(struct cmdline_data* cmd, struct  global_data* gd,
                             nodeptr p, nodeptr q);
global bool nodes_set_bin(struct cmdline_data* cmd, struct  global_data* gd,
                          nodeptr p, nodeptr q, int *n, real *dr1, vector dr);
#ifdef SINGLEP
global bool nodes_set_bin5(struct cmdline_data* cmd, struct  global_data* gd,
                           nodeptr, nodeptr, int *, float *, float *);
#else
global bool nodes_set_bin5(struct cmdline_data* cmd, struct  global_data* gd,
                           nodeptr p, nodeptr q, int *n, real *dr1, vector dr);
#endif

//B BALLS4
global int search_init_omp_balls6(struct cmdline_data* cmd, struct  global_data* gd,
                                  gdhistptr_omp_balls6, int);
global int search_init_sincos_omp_balls6(struct cmdline_data* cmd, 
                                         struct  global_data* gd,
                                         gdhistptr_sincos_omp_balls6);
global int search_init_omp_balls6_cc(struct cmdline_data* cmd, 
                                     struct  global_data* gd,
                                     gdhistptr_omp_balls6);
global int search_free_omp_balls6(struct cmdline_data* cmd, struct  global_data* gd,
                                  gdhistptr_omp_balls6);
global int search_free_omp_balls6_cc(struct cmdline_data* cmd, 
                                     struct  global_data* gd,
                                     gdhistptr_omp_balls6);
global int search_free_sincos_omp_balls6(struct cmdline_data* cmd, 
                                         struct  global_data* gd,
                                         gdhistptr_sincos_omp_balls6);

global int computeBodyProperties_omp_balls6(struct cmdline_data* cmd, 
                                            struct  global_data* gd,
                                            bodyptr, int, gdhistptr_omp_balls6);
global int computeBodyProperties_omp_balls6_cc(struct cmdline_data* cmd, 
                                               struct  global_data* gd,
                                               bodyptr, INTEGER,
                                               gdhistptr_omp_balls6);
global int computeBodyProperties_sincos_omp_balls6_cc(struct cmdline_data* cmd, 
                                                      struct  global_data* gd,
                                                      bodyptr, INTEGER,
                                               gdhistptr_sincos_omp_balls6);
//E BALLS4

global int search_compute_HistN_balls(struct cmdline_data* cmd, 
                                      struct  global_data* gd,
                                      int nbody);

#endif

#endif	// ! _protodefs_balls_omp_h
