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
//B 2023.11.29
global void maketree(bodyptr btab, int nbody, int ifile);
//E
//B Creating tree for nodes. In (ADDON)treenodeload.c ... move to protodefs_01.h
//global void maketreenodes(bodyptr, INTEGER);
//E

global void searchcalc_normal_omp_sincos(bodyptr, int, INTEGER, INTEGER);

//E

//B Tree utilities
global void doBoxWrapping(void);
global bool reject_cell(nodeptr, nodeptr, real);
global bool reject_cell_balls(nodeptr, nodeptr, real *, vector);
global bool reject_bodycell(nodeptr, nodeptr);
global bool reject_cellcell(nodeptr, nodeptr);
//B 2023.11.22
global bool reject_balls(nodeptr p, nodeptr q, real *drpq, vector dr);
global bool nodes_condition_balls(nodeptr p, nodeptr q, real *dr1, vector dr);
//E
global bool accept_body(bodyptr, nodeptr, real *, vector);
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


#ifdef ADDONS
#ifdef ADDONSDEVELOP
#include "protodefs_01.h"
#endif

#ifdef ADDONSDEVELOP
#include "protodefs_02.h"
#endif
#endif

//B Other utilities
global int ThreadCount(INTEGER, int);
global int spherical_to_cartesians(real, real, vector);
global int spherical_periodic_condition(real *, real *, real *, real *);
//E

#ifdef ADDONS
#ifdef DIRECTMETHOD
#include "protodefs_direct_method.h"
#endif
#ifdef TREEOMP
#include "protodefs_tree_normal_omp.h"
#endif
#ifdef TREE3PCFDIRECTOMP
#include "protodefs_tree_3pcf_direct_omp.h"
#endif
#ifdef BALLS
#include "protodefs_balls_omp.h"
#endif
#endif

#endif // ! _protodefs_h
