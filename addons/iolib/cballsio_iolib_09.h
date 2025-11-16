// Use:
//#include "cballsio_iolib_09.h"

// included in: addons/addons_include/cballsio/cballsio_include_09.h

//B Output file formats:
//#define OUTCOLUMNS       0
//#define OUTCOLUMNSALL    3
//#define OUTCOLUMNSBIN    2
//#define OUTCOLUMNSBINALL 4
//#define OUTNULL          1
//#define OUTCFITSIO          8
//#define OUTNUMPYHEALPIX     9
//E

#ifndef _cballsio_iolib_09_h
#define _cballsio_iolib_09_h

if (strcmp(outfmt_str,"columns-ascii-pos") == 0)
    *outfmt_int = OUTCOLUMNSPOS;

#endif	// ! _cballsio_iolib_09_h
