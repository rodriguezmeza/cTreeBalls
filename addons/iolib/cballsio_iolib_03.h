// Use:
//#include "cballsio_iolib_03.h"

#ifndef _cballsio_iolib_03_h
#define _cballsio_iolib_03_h

if (strcmp(infmt_str,"columns-ascii-pos") == 0)     *infmt_int = INCOLUMNSPOS;

if (strcmp(infmt_str,"columns-ascii-2d-to-3d") == 0)
                                                *infmt_int = INCOLUMNS2DTO3D;

if (strcmp(infmt_str,"multi-columns-ascii") == 0) *infmt_int = INMCOLUMNS;

#endif	// ! _cballsio_iolib_03_h
