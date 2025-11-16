// Use:
//#include "cballsio_iolib_08.h"

// included in: addons/addons_include/source/cballsio/cballsio_indlude_08.h

#ifndef _cballsio_iolib_08_h
#define _cballsio_iolib_08_h

if (strcmp(infmt_str,"columns-ascii-pos") == 0)
    *infmt_int = INCOLUMNSPOS;

if (strcmp(infmt_str,"columns-ascii-2d-to-3d") == 0)
    *infmt_int = INCOLUMNS2DTO3D;

if (strcmp(infmt_str,"multi-columns-ascii") == 0)
    *infmt_int = INMCOLUMNS;

if (strcmp(infmt_str,"ra-dec-ascii") == 0)
    *infmt_int = INRADECCOLUMNS;

#endif	// ! _cballsio_iolib_08_h
