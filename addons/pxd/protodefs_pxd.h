// Use:
//#include "protodefs_pxd.h"

#ifndef _protodefs_pxd_h
#define _protodefs_pxd_h

//B PXD functions

//B flags
int get_tree_allocated(struct global_data* gd, short *value);
int get_allocated_2(struct global_data* gd, short *value);
int get_bodytable_allocated(struct global_data* gd, short *value);
int get_histograms_allocated(struct global_data* gd, short *value);
int get_gd_allocated(struct global_data* gd, short *value);
int get_cmd_allocated(struct global_data* gd, short *value);
//E

//B parameters section
int get_nthreads(struct  cmdline_data* cmd,
              int *value);
//B version 1.0.1
int get_nmultipoles(struct  cmdline_data* cmd,
              int *value);
int get_theta(struct  cmdline_data* cmd,
              real *theta);
int get_rsmooth(struct  global_data* gd, real *value);
int get_cputime(struct  global_data* gd, real *cputime);
int get_sizeHistN(struct  cmdline_data* cmd, int *sizeHistN);
int get_version(struct  cmdline_data* cmd, char *param);
int get_rootDir(struct  cmdline_data* cmd, char *value);
//E parameters section

//B global parameter section
int get_nbody(struct  cmdline_data* cmd, struct  global_data* gd, int *value);
int get_computeTPCF(struct  cmdline_data* cmd, struct  global_data* gd, short *value);
//E global parameter section

//B histograms section
int get_rBins(struct  cmdline_data* cmd,
              struct  global_data* gd);
int get_HistNN(struct  cmdline_data* cmd, struct  global_data* gd);
int get_HistCF(struct  cmdline_data* cmd, struct  global_data* gd);
int get_HistXi2pcf(struct  cmdline_data* cmd, struct  global_data* gd);
//B cross
int get_HistXi2pcf12(struct  cmdline_data* cmd, struct  global_data* gd);
int get_HistXi2pcf13(struct  cmdline_data* cmd, struct  global_data* gd);
//E
int get_HistZetaMsincos(struct  cmdline_data* cmd,
                         struct  global_data* gd,
                         int m, int type, ErrorMsg errmsg);
// (EE) edge_effects
int get_HistZetaM_EE(struct  cmdline_data* cmd,
                         struct  global_data* gd,
                         int m, ErrorMsg errmsg);
//E histograms section

//E PXD functions

#endif	// ! _protodefs_pxd_h
