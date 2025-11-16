
#ifndef _cballsio_iolib_00_h
#define _cballsio_iolib_00_h

// This tag number must be diferent than others use in I/O like
//#define INGADGET      5
//
#define INCOLUMNS2DTO3D 4
// To read only position columns (ascii)
#define INCOLUMNSPOS    6
#define INMCOLUMNS      7
// check cballsio_cfitsio_00.h for IO-TAG numbers
#define INRADECCOLUMNS  12

//B output tags
// To save only position columns (ascii)
#define OUTCOLUMNSPOS   10
//E

local int inputdata_ascii_pos(struct cmdline_data*, struct  global_data*,
                               string filename, int);
#if NDIM == 3
local int inputdata_ascii_2d_to_3d(struct cmdline_data*, struct  global_data*,
                                    string filename, int);
#endif
local int inputdata_ascii_mcolumns(struct cmdline_data* cmd,
                                   struct  global_data* gd,
                                   string filename, int);
local int inputdata_ascii_ra_dec(struct cmdline_data* cmd, struct  global_data* gd,
                                 string filename, int ifile);
local int outputdata_ascii_pos(struct cmdline_data*, struct  global_data*,
                         bodyptr, INTEGER);

#endif	// ! _cballsio_iolib_00_h
