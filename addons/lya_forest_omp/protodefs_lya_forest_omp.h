#ifndef _protodefs_lya_forest_omp_h
#define _protodefs_lya_forest_omp_h

#include "lya_forest_defs.h"

#define LYAFOREST1D2PCFMETHOD 180
#define LYAFOREST1D3PCFMETHOD 181
#define LYAFOREST1D2PCF3PCFMETHOD 182
#define LYAFOREST1DTREE2PCFMETHOD 183

global int inputdata_lya_ascii(struct cmdline_data *cmd,
                               struct global_data *gd,
                               string filename, int ifile);
global int searchcalc_lya_forest_omp(struct cmdline_data *cmd,
                                     struct global_data *gd,
                                     bodyptr *btable, INTEGER *nbody,
                                     INTEGER ipmin, INTEGER *ipmax,
                                     int cat, int compute_2pcf,
                                     int compute_3pcf);
global int searchcalc_lya_forest_1d_omp(struct cmdline_data *cmd,
                                        struct global_data *gd,
                                        bodyptr table, INTEGER nbody,
                                        int compute_2pcf,
                                        int compute_3pcf);
global int searchcalc_lya_forest_1d_tree_omp(struct cmdline_data *cmd,
                                             struct global_data *gd,
                                             bodyptr table, INTEGER nbody);

#endif
