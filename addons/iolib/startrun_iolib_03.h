// Use:
//#include "startrun_iolib_03.h"

#ifndef _startrun_iolib_03_h
#define _startrun_iolib_03_h

if ((int)gd->columns[0] < 1)
    error("CheckParameters: absurd value for columns[0] (%s)\n",
          cmd->columns);
if ((int)gd->columns[1] < 1)
    error("CheckParameters: absurd value for columns[1] (%s)\n",
          cmd->columns);
if ((int)gd->columns[2] < 1)
    error("CheckParameters: absurd value for columns[2] (%s)\n",
          cmd->columns);

#endif	// ! _startrun_iolib_03_h
