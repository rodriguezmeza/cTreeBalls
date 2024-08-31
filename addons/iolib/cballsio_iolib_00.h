
#ifndef _cballsio_iolib_00_h
#define _cballsio_iolib_00_h

// This tag number must be diferent than others use in I/O
//#define INGADGET                5
#define INCOLUMNS2DTO3D 4
#define INCOLUMNSPOS    6                           
// To read only
//  position columns (ascii)
#define INMCOLUMNS      7

//B Some global definitions needed for reading formats
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

#endif	// ! _cballsio_iolib_00_h
