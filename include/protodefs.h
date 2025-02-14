/*==============================================================================
 HEADER: protodefs.h				[cTreeBalls]
 Written by: Mario A. Rodriguez-Meza
 Starting date: april 2023
 Purpose: Definitions of global prototypes
 Language: C
 Use: '#include "protodefs.h"
 Major revisions:
 ==============================================================================*/
//        1          2          3          4        ^ 5          6          7

#ifndef _protodefs_h
#define _protodefs_h

#ifdef __cplusplus
extern "C" {
#endif


int OutputData(struct cmdline_data* cmd, struct  global_data* gd,
           bodyptr *btable, INTEGER *nbody, int ifile);

int MainLoop(struct cmdline_data* cmd, struct  global_data* gd);
int StartRun(struct cmdline_data* cmd, struct  global_data* gd,
             string, string, string, string);

int StartRun_Common(struct cmdline_data*, struct  global_data*);
int PrintParameterFile(struct cmdline_data *, char *);

//B If uncommented there will be a warning in the setup.py process
//#ifdef OPENMPCODE
int SetNumberThreads(struct cmdline_data *cmd);
//#endif
//E

int StartOutput(struct cmdline_data *, struct  global_data*);
//global void checkstop(void);
int EndRun(struct cmdline_data* cmd, struct  global_data* gd);
//global void savestate(string);
//global int restorestate(string);
global int startrun_memoryAllocation(struct cmdline_data* cmd, struct  global_data* gd);

int TestData(struct cmdline_data* cmd, struct  global_data* gd);


//B I/O
global int infilefmt_string_to_int(string, int *);
int InputData(struct cmdline_data* cmd,
              struct  global_data* gd, string filename, int);

//B I/O directories:
global void setFilesDirs_log(struct cmdline_data*, struct  global_data* gd);
global void setFilesDirs(struct cmdline_data*, struct  global_data* gd);
//E
//E


int EvalHist(struct cmdline_data* cmd, struct  global_data* gd);

//B 3PCF section
// Search methods:
global int searchcalc_direct_omp(bodyptr, int, INTEGER, INTEGER);


//B Tree:
global int MakeTree(struct cmdline_data* cmd, struct  global_data* gd,
                    bodyptr btab, INTEGER nbody, int ifile);
//E

global int searchcalc_normal_sincos(struct cmdline_data* cmd,
                                        struct  global_data* gd,
                                        bodyptr *,
                                         INTEGER *, INTEGER, INTEGER *, int, int);

//E

//B treeload utilities
global int expandbox(struct  cmdline_data*, struct  global_data*,
                    bodyptr, int, int, cellptr);
global int findRootCenter(struct  cmdline_data* cmd, struct  global_data* gd,
                         bodyptr, int, int, cellptr);
global int centerBodies(bodyptr, int, int, cellptr);
//global cellptr makecell(struct  cmdline_data* cmd, struct  global_data* gd, int);
//E

//B cBalls utilities
global int doBoxWrapping(struct cmdline_data* cmd, struct  global_data* gd);
global bool reject_cell(struct cmdline_data* cmd, struct  global_data* gd,
                        nodeptr, nodeptr, real);
global bool reject_cell_balls(struct cmdline_data* cmd, struct  global_data* gd,
                              nodeptr, nodeptr, real *, vector);
global bool reject_bodycell(struct cmdline_data* cmd, struct  global_data* gd,
                            nodeptr, nodeptr);
global bool reject_cellcell(struct cmdline_data* cmd, struct  global_data* gd,
                            nodeptr, nodeptr);

global bool reject_balls(struct cmdline_data* cmd, struct  global_data* gd,
                         nodeptr p, nodeptr q, real *drpq, vector dr);
global bool nodes_condition_balls(struct cmdline_data* cmd, struct  global_data* gd,
                                  nodeptr p, nodeptr q, real *dr1, vector dr);

#ifdef SINGLEP
global bool accept_body(struct cmdline_data* cmd, struct  global_data* gd,
                        bodyptr, nodeptr, float *, float *);
#else
global bool accept_body(struct cmdline_data* cmd, struct  global_data* gd,
                        bodyptr, nodeptr, real *, vector);
#endif
//B Obsolete
//global int compute_cosphi(real dr1, vector dr, real *cosphi, gdhist hist);

//#ifdef OPENMPCODE
global int search_init_sincos_omp(struct cmdline_data* cmd, 
                                  struct  global_data* gd,
                                  gdhistptr_sincos_omp hist);
global int search_free_sincos_omp(struct cmdline_data* cmd, struct  global_data* gd,
                                  gdhistptr_sincos_omp hist);
global int computeBodyProperties_sincos(struct cmdline_data* cmd,
                                            struct  global_data* gd,
                                            bodyptr, int,
                                            gdhistptr_sincos_omp);
//#endif


global int search_init_gd_hist(struct cmdline_data* cmd, struct  global_data* gd);
global int search_init_gd_hist_sincos(struct cmdline_data* cmd, struct  global_data* gd);
global int search_compute_HistN(struct cmdline_data* cmd, struct  global_data* gd,
                                int nbody);
//E


//B Other utilities
global int ThreadCount(struct cmdline_data* cmd, struct  global_data* gd, INTEGER, int);
global int spherical_to_cartesians(struct cmdline_data* cmd, struct  global_data* gd,
                                   real, real, vector);
global int spherical_periodic_condition(real *, real *, real *, real *);


global int statHistogram(struct cmdline_data* cmd, struct  global_data* gd);
//E


#ifdef ADDONS
// If you have an addon that need global protodefinitions
//  go to this file and add the addon item.
#include "protodefs_include.h"
#endif


#ifdef __cplusplus
}
#endif

#endif // ! _protodefs_h
