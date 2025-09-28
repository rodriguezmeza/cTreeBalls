
#ifndef _cballsio_cfitsio_00_h
#define _cballsio_cfitsio_00_h

#include "fitsio.h"
#include "chealpix.h"

// These tag numbers must be diferent than others used in I/O
#define INCFITSIO           8
#define INFITSRADECTC       80
#define INFITSHEALPIX       81
#define INNUMPYHEALPIX      82

#define OUTCFITSIO          8
#define OUTNUMPYHEALPIX     9

local int inputdata_cfitsio(struct cmdline_data*, struct  global_data*,
                               string filename, int);
local int inputdata_cfitsio_radec_tc(struct cmdline_data*,
                                             struct  global_data*,
                                             string filename, int);
local int inputdata_cfitsio_healpix(struct cmdline_data*,
                                             struct  global_data*,
                                             string filename, int);
local int inputdata_cfitsio_ra_dec(struct cmdline_data* cmd,
                                   struct  global_data* gd,
                                   string filename, int ifile, fitsfile *fptr);
local int inputdata_cfitsio_xyz(struct cmdline_data* cmd,
                                struct  global_data* gd,
                                string filename, int ifile, fitsfile *fptr);
local int inputdata_cfitsio_healpix_map(struct cmdline_data* cmd,
                                struct  global_data* gd,
                                        string filename, int ifile,
                                        fitsfile *fptr);
local int inputdata_cfitsio_healpix_map_mask(struct cmdline_data* cmd,
                                struct  global_data* gd,
                                        string filename, int ifile,
                                             fitsfile *fptr);
local int inputdata_cfitsio_healpix_map_mask_inside(struct cmdline_data* cmd,
                                struct  global_data* gd,
                                        string filename, int ifile,
                                        fitsfile *fptr);
local int outputdata_cfitsio(struct cmdline_data*, struct  global_data*,
                         bodyptr, INTEGER);
void writebintable_kappa(struct cmdline_data* cmd,
                         struct  global_data* gd,
                         bodyptr bodytab, INTEGER nbody,
                         string filename);
void printerror( int status);

//B infileformat: numpy-healpix
local int inputdata_numpy_healpix(struct cmdline_data* cmd,
                                  struct  global_data* gd,
                                  string filename, int ifile);
local int inputdata_numpy_healpix_map(struct cmdline_data* cmd,
                                struct  global_data* gd,
                                        string filename, int ifile);
local int inputdata_numpy_healpix_map_mask(struct cmdline_data* cmd,
                                struct  global_data* gd,
                                        string filename, int ifile);
local int inputdata_numpy_healpix_map_mask_inside(struct cmdline_data* cmd,
                                struct  global_data* gd,
                                        string filename, int ifile);

// outfileformat: numpy-healpix
local int outputdata_numpy_healpix(struct cmdline_data* cmd,
                                   struct  global_data* gd,
                                   bodyptr, INTEGER);
//E

#endif	// ! _cballsio_cfitsio_00_h
