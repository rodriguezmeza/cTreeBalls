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
if ((int)gd->columns[3] < 1)
    error("CheckParameters: absurd value for columns[3] (%s)\n",
          cmd->columns);
if ((int)gd->columns[4] < 1)
    error("CheckParameters: absurd value for columns[4] (%s)\n",
          cmd->columns);
if ((int)gd->columns[5] < 1)
    error("CheckParameters: absurd value for columns[5] (%s)\n",
          cmd->columns);
if ((int)gd->columns[6] < 1)
    error("CheckParameters: absurd value for columns[6] (%s)\n",
          cmd->columns);

#endif	// ! _startrun_iolib_03_h
