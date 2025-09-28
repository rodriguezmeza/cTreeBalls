// Use:
//#include "cballsio_cfitsio_09.h"


//B Output file formats:
//#define OUTCOLUMNS       0
//#define OUTCOLUMNSBIN    2
//#define OUTNULL          1
//E

#ifndef _cballsio_cfitsio_09_h
#define _cballsio_cfitsio_09_h

if (strcmp(outfmt_str,"fits") == 0)     *outfmt_int = OUTCFITSIO;
if (strcmp(outfmt_str,"numpy-healpix") == 0) *outfmt_int = OUTNUMPYHEALPIX;

#endif	// ! _cballsio_cfitsio_09_h
