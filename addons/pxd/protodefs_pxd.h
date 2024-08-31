// Use:
//#include "protodefs_pxd.h"

#ifndef _protodefs_pxd_h
#define _protodefs_pxd_h

//B PXD functions
int get_theta(struct  cmdline_data* cmd,
              real *theta);
int get_sizeHistN(struct  cmdline_data* cmd, int *sizeHistN);
int get_rBins(struct  cmdline_data* cmd,
              struct  global_data* gd);
int get_HistZetaM_sincos(struct  cmdline_data* cmd,
                         struct  global_data* gd,
                         int m, int type, ErrorMsg errmsg);
//E PXD functions

#endif	// ! _protodefs_pxd_h
