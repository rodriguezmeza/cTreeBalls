// Use:
//#include "cballsio_cfitsio_08.h"

#ifndef _cballsio_cfitsio_08_h
#define _cballsio_cfitsio_08_h

if (strcmp(infmt_str,"fits") == 0)     *infmt_int = INCFITSIO;
if (strcmp(infmt_str,"fits-radec-tc") == 0)
    *infmt_int = INFITSRADECTC;
if (strcmp(infmt_str,"fits-healpix") == 0)
    *infmt_int = INFITSHEALPIX;

#endif	// ! _cballsio_cfitsio_08_h
