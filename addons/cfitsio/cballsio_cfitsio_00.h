
#ifndef _cballsio_cfitsio_00_h
#define _cballsio_cfitsio_00_h

// This tag number must be diferent than others use in I/O
#define INCFITSIO      8
#define OUTCFITSIO     8

local int inputdata_cfitsio(struct cmdline_data*, struct  global_data*,
                               string filename, int);

local int outputdata_cfitsio(struct cmdline_data*, struct  global_data*,
                         bodyptr, INTEGER);

#endif	// ! _cballsio_cfitsio_00_h
