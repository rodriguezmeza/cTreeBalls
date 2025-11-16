// Use:
//#include "protodefs_pxd.h"

#ifndef _protodefs_pxd_h
#define _protodefs_pxd_h

//B PXD functions

//B parameters section
int get_theta(struct  cmdline_data* cmd,
              real *theta);
int get_sizeHistN(struct  cmdline_data* cmd, int *sizeHistN);
int get_version(struct  cmdline_data* cmd, char *param);
//E parameters section

//B histograms section
int get_rBins(struct  cmdline_data* cmd,
              struct  global_data* gd);
int get_HistNN(struct  cmdline_data* cmd, struct  global_data* gd);
int get_HistCF(struct  cmdline_data* cmd, struct  global_data* gd);
int get_HistXi2pcf(struct  cmdline_data* cmd, struct  global_data* gd);
int get_HistZetaMsincos(struct  cmdline_data* cmd,
                         struct  global_data* gd,
                         int m, int type, ErrorMsg errmsg);
int get_HistZetaM_sincos(struct  cmdline_data* cmd,
                         struct  global_data* gd,
                         int m, int type, ErrorMsg errmsg);
//E histograms section

//E PXD functions

#endif	// ! _protodefs_pxd_h
