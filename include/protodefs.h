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

int inputdata(string filename, int);

// I/O directories:
global void setFilesDirs_log(void);
global void setFilesDirs(void);

//B 3PCF section
// Search methods:
global int searchcalc_direct_omp(bodyptr, int, INTEGER, INTEGER);


//B Tree:
global int maketree(bodyptr btab, INTEGER nbody, int ifile);
//E

global void searchcalc_normal_omp_sincos(bodyptr *,
                                         INTEGER *, INTEGER, INTEGER *, int, int);

//E

//B Tree utilities
global void doBoxWrapping(void);
global bool reject_cell(nodeptr, nodeptr, real);
global bool reject_cell_balls(nodeptr, nodeptr, real *, vector);
global bool reject_bodycell(nodeptr, nodeptr);
global bool reject_cellcell(nodeptr, nodeptr);

global bool reject_balls(nodeptr p, nodeptr q, real *drpq, vector dr);
global bool nodes_condition_balls(nodeptr p, nodeptr q, real *dr1, vector dr);

#ifdef SINGLEP
global bool accept_body(bodyptr, nodeptr, float *, float *);
#else
global bool accept_body(bodyptr, nodeptr, real *, vector);
#endif
global int compute_cosphi(real dr1, vector dr, real *cosphi, gdhist hist);

#ifdef OPENMPCODE
global int search_init_sincos_omp(gdhistptr_sincos_omp hist);
global int search_free_sincos_omp(gdhistptr_sincos_omp hist);
global int computeBodyProperties_sincos_omp(bodyptr, int, gdhistptr_sincos_omp);
#endif


global int search_init_gd_hist(void);
global int search_init_gd_hist_sincos(void);
global int search_compute_HistN(int nbody);
//E


//B Other utilities
global int ThreadCount(INTEGER, int);
global int spherical_to_cartesians(real, real, vector);
global int spherical_periodic_condition(real *, real *, real *, real *);
//E


#ifdef ADDONS
// If you have an addon that need global protodefinitions
//  go to this file and add the addon item.
#include "protodefs_include.h"
#endif

#endif // ! _protodefs_h
