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

//
// lines where there is a "//B socket:" string are places to include module files
//  that can be found in addons/addons_include folder
//

#ifndef _protodefs_h
#define _protodefs_h

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif


int OutputData(struct cmdline_data* cmd, struct  global_data* gd,
           bodyptr *btable, INTEGER *nbody, int ifile);

int MainLoop(struct cmdline_data* cmd, struct  global_data* gd);
int StartRun(struct cmdline_data* cmd, struct  global_data* gd,
             string, string, string, string);

int StartRun_Common(struct cmdline_data*, struct  global_data*);
int PrintParameterFile(struct cmdline_data *, struct  global_data*, char *);

int SetNumberThreads(struct cmdline_data *cmd);

int cballs_start_run_common_guarded(struct cmdline_data *cmd,
                                    struct global_data *gd);
int cballs_print_parameter_file_guarded(struct cmdline_data *cmd,
                                        struct global_data *gd,
                                        char *filename);
int cballs_set_number_threads_guarded(struct cmdline_data *cmd);
int cballs_main_loop_guarded(struct cmdline_data *cmd,
                             struct global_data *gd);
int cballs_end_run_guarded(struct cmdline_data *cmd,
                           struct global_data *gd);
int cballs_end_run_free_memory_guarded(struct cmdline_data *cmd,
                                       struct global_data *gd);
int cballs_compiled_ndim(void);
int cballs_max_memory_catalogs(void);
int cballs_search_method_id(const char *method);
int cballs_set_memory_forest_ids(struct cmdline_data *cmd,
                                struct global_data *gd, int ifile,
                                const int64_t *forest_ids, size_t nbody);
int cballs_load_memory_catalog(struct cmdline_data *cmd,
                               struct global_data *gd,
                               int ifile,
                               const double *positions,
                               size_t nbody,
                               int ndim,
                               const double *kappa,
                               const double *weights,
                               const unsigned char *mask,
                               const double *gamma1,
                               const double *gamma2);

int StartOutput(struct cmdline_data *, struct  global_data*);
int EndRun(struct cmdline_data* cmd, struct  global_data* gd);
global int startrun_memoryAllocation(struct cmdline_data* cmd, struct  global_data* gd);
global int EndRun_FreeMemory(struct cmdline_data* cmd,
                             struct  global_data* gd);
global int freeTree(struct  cmdline_data* cmd, struct  global_data* gd);
global int EndRun_FreeMemory_cmd(struct cmdline_data* cmd,
                             struct  global_data* gd);
global int EndRun_FreeMemory_gd(struct cmdline_data* cmd,
                             struct  global_data* gd);
global int EndRun_FreeMemory_gd_2(struct cmdline_data* cmd,
                             struct  global_data* gd);
global int EndRun_FreeMemory_histograms(struct cmdline_data* cmd,
                             struct  global_data* gd);
global int EndRun_FreeMemory_tree(struct cmdline_data* cmd,
                             struct  global_data* gd);
global int EndRun_FreeMemory_bodytable(struct cmdline_data* cmd,
                             struct  global_data* gd);

int TestData(struct cmdline_data* cmd, struct  global_data* gd);


//B I/O
global int infilefmt_string_to_int(string, int *);
int InputData(struct cmdline_data* cmd,
              struct  global_data* gd, string filename, int);
global int InputData_all_in_one(struct cmdline_data* cmd,
                                struct  global_data* gd);

//B I/O directories:
global int setFilesDirs_log(struct cmdline_data*, struct  global_data* gd);
global int setFilesDirs(struct cmdline_data*, struct  global_data* gd);
//E

//B routine to compute edge corrections using two saved histZetaM histograms
global int computeEdgeCorrections(struct cmdline_data* cmd,
                                  struct  global_data* gd);
global int matrixClm(struct cmdline_data* cmd, struct  global_data* gd,
                     double ***mat3, double ***mat4,
                     int n1, int n2, double ***mat5);
//E


int EvalHist(struct cmdline_data* cmd, struct  global_data* gd);

//B 3PCF section
// Search methods:

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
global int FindRootCenter(struct  cmdline_data* cmd, struct  global_data* gd,
                         bodyptr, int, int, cellptr);
global int centerBodies(bodyptr, int, int, cellptr);
//E

//B cBalls utilities
global int doBoxWrapping(struct cmdline_data* cmd, struct  global_data* gd);
global bool reject_cell(struct cmdline_data* cmd, struct  global_data* gd,
                        nodeptr, nodeptr, real);
global bool reject_cell_balls(struct cmdline_data* cmd, struct  global_data* gd,
                              nodeptr, nodeptr, real *, compute_vector);
global bool reject_bodycell(struct cmdline_data* cmd, struct  global_data* gd,
                            nodeptr, nodeptr);
global bool reject_cellcell(struct cmdline_data* cmd, struct  global_data* gd,
                            nodeptr, nodeptr);

global bool reject_balls(struct cmdline_data* cmd, struct  global_data* gd,
                         nodeptr p, nodeptr q, real *drpq, compute_vector dr);

global bool accept_body(struct cmdline_data* cmd, struct  global_data* gd,
                        bodyptr, nodeptr, real *, compute_vector);

#ifdef SMOOTHPIVOT
global int prepare_smooth_pivots(struct cmdline_data* cmd,
                                 struct global_data* gd,
                                 bodyptr *btable, INTEGER *nbody,
                                 INTEGER ipmin, INTEGER *ipmax,
                                 int cat1, int cat2);
#endif

global int search_init_sincos_omp(struct cmdline_data* cmd,
                                  struct  global_data* gd,
                                  gdhistptr_sincos_omp hist);
global int search_free_sincos_omp(struct cmdline_data* cmd, struct  global_data* gd,
                                  gdhistptr_sincos_omp hist);
global int computeBodyProperties_sincos(struct cmdline_data* cmd,
                                            struct  global_data* gd,
                                            bodyptr, int,
                                            gdhistptr_sincos_omp);

global int search_init_gd_hist(struct cmdline_data* cmd, struct  global_data* gd);
global int search_init_gd_hist_sincos(struct cmdline_data* cmd, struct  global_data* gd);
global int search_compute_HistN(struct cmdline_data* cmd, struct  global_data* gd,
                                int nbody);
//E


//B Other utilities
global int ThreadCount(struct cmdline_data* cmd,
                       struct  global_data* gd, INTEGER, int);
global int coordinate_string_to_int(struct cmdline_data* cmd,
                                    struct  global_data* gd);
global int coordinate_transformation(struct cmdline_data* cmd,
                                   struct  global_data* gd,
                                   real, real, vector);
global int spherical_periodic_condition(real *, real *, real *, real *);

global int statHistogram(struct cmdline_data* cmd, struct  global_data* gd);
//E

//B socket:
#ifdef ADDONS
// If you have an addon that need global proto definitions
//  go to this file and add the addon item.
#include "protodefs_include.h"
#endif
//E

// Application Binary Interface (ABI) definitions
size_t sizeof_cmdline_data(void);
size_t sizeof_global_data(void);


global int cballs_system_checked(struct cmdline_data *cmd,
                                 string routineName,
                                 string command);
global int cballs_stream_close_checked(struct cmdline_data *cmd,
                                       string routineName,
                                       stream *stream_ptr,
                                       string filename);
global real cballs_normalize_or_zero(real numerator, real denominator);
#ifdef __cplusplus
}
#endif

#endif // ! _protodefs_h
