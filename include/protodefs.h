/*==============================================================================
 HEADER: protodefs.h				[cTreeBalls]
 Written by: Mario A. Rodriguez-Meza
 Starting date: april 2023
 Purpose: Definitions of global prototypes
 Language: C
 Use: '#include "protodefs.h"
 Major revisions:
 ==============================================================================*/
//        1          2          3          4          5          6          7

#ifndef _protodefs_h
#define _protodefs_h

int output(void);

int MainLoop(void);
int StartRun(string, string, string, string);
int StartOutput(void);
global void checkstop(void);
int EndRun(void);
global void savestate(string);
global int restorestate(string);
global int startrun_memoryAllocation(void);

int testdata(void);

int inputdata(void);

// I/O directories:
global void setFilesDirs_log(void);
global void setFilesDirs(void);

//B 3PCF section
// Search methods:
global int searchcalc_direct_omp(bodyptr btab, int nbody, INTEGER ipmin, INTEGER ipmax);


//B Tree:
global void maketree(bodyptr btab, int nbody);
global void maketreenodes(bodyptr ntab, int nnode);

global void searchcalc_normal_omp(bodyptr btab, int nbody, INTEGER ipmin, INTEGER ipmax);
global void searchcalc_normal_omp_sincos(bodyptr btab, int nbody, INTEGER ipmin, INTEGER ipmax);

//E

//B Tree utilities
global void doBoxWrapping(void);
global bool reject_cell(nodeptr, nodeptr, real);
global bool reject_cell_balls(nodeptr, nodeptr, real *, vector);
global bool reject_bodycell(nodeptr, nodeptr);
global bool reject_cellcell(nodeptr, nodeptr);
global bool accept_body(bodyptr, nodeptr, real *, vector);
global int compute_cosphi(real dr1, vector dr, real *cosphi, gdhist hist);

#ifdef OPENMPCODE
global int search_init_omp(gdhistptr_omp hist);
global int search_init_sincos_omp(gdhistptr_sincos_omp hist);
global int search_free_omp(gdhistptr_omp hist);
global int search_free_sincos_omp(gdhistptr_sincos_omp hist);
global int computeBodyProperties_omp(bodyptr p, int nbody, gdhistptr_omp hist);
global int computeBodyProperties_sincos_omp(bodyptr p, int nbody, gdhistptr_sincos_omp hist);

global int search_init_omp_3pcfbf(gdhistptr_omp_3pcfbf hist);
global int search_free_omp_3pcfbf(gdhistptr_omp_3pcfbf hist);
global int computeBodyProperties_omp_3pcfbf(bodyptr p, int nbody, gdhist_omp_3pcfbf hist);
global int search_compute_HistN_3pcfbf(int nbody);
global void searchcalc_normal_3pcf_direct_omp(bodyptr btab, int nbody, INTEGER ipmin, INTEGER ipmax);
//B BALLS4
global int search_init_omp_balls6(gdhistptr_omp_balls6 hist);
global int search_init_omp_balls6_cc(gdhistptr_omp_balls6 hist);
global int search_free_omp_balls6(gdhistptr_omp_balls6 hist);
global int search_free_omp_balls6_cc(gdhistptr_omp_balls6 hist);

global int computeBodyProperties_omp_balls6(bodyptr p, int nbody, gdhistptr_omp_balls6 hist);
global int computeBodyProperties_omp_balls6_cc(bodyptr p, int nbody, gdhistptr_omp_balls6 hist);
//E BALLS4

#endif

#ifdef BALLS
global void searchcalc_balls_omp(bodyptr btab, int nbody, INTEGER ipmin, INTEGER ipmax);
global int search_init_balls_omp(gdhistptr_omp_balls hist);
global int search_init_balls_omp_cc(gdhistptr_omp_balls hist);
global int computeBodyProperties_balls_omp(bodyptr p, int nbody, gdhistptr_omp_balls hist);
global int computeBodyProperties_balls_omp_cc(bodyptr p, int nbody, gdhistptr_omp_balls hist);
global int computeBodyProperties_balls_omp_cc_sincos(bodyptr p, int nbody, gdhistptr_sincos_omp hist);
global int search_free_balls_omp(gdhistptr_omp_balls hist);
global int search_free_balls_omp_sincos(gdhistptr_omp_balls hist);
global int search_free_balls_omp_cc(gdhistptr_omp_balls hist);
global bool nodes_condition(nodeptr p, nodeptr q, real *dr1, vector dr);
global bool nodes_condition5(nodeptr p, nodeptr q);
global bool nodes_set_bin(nodeptr p, nodeptr q, int *n, real *dr1, vector dr);
global bool nodes_set_bin5(nodeptr p, nodeptr q, int *n, real *dr1, vector dr);
#endif

global int search_init_gd_hist(void);
global int search_init_gd_hist_sincos(void);
global int search_compute_HistN(int nbody);
global int search_compute_HistN_balls(int nbody);
//E




//B Other utilities
global int ThreadCount(INTEGER);
global int spherical_to_cartesians(real, real, vector);
global int spherical_periodic_condition(real *, real *, real *, real *);
//E



#endif // ! _protodefs_h
