// Use:
//#include "startrun_balls_omp_03.h"

#ifndef _startrun_balls_omp_03_h
#define _startrun_balls_omp_03_h

if (scanopt(cmd->options, "bodyfound"))
    if (cmd->ntosave < 1 || cmd->ntosave > cmd->nbody)
        error("CheckParameters: absurd value for ntosave\n");
if (cmd->scanLevel < 0)
    error("CheckParameters: absurd value for scanLevel (%d)\n",cmd->scanLevel);
// Root nodes:
if (cmd->scanLevelRoot < 0)
    error("CheckParameters: absurd value for scanLevelRoot (%d)\n",cmd->scanLevelRoot);
if ((int)gd->scanLevelMin[0] > 0)
    error("CheckParameters: absurd value for scanLevelMin[0] (%s)\n",cmd->scanLevelMin);
if (gd->scanLevelMin[1] > 0)
    error("CheckParameters: absurd value for scanLevelMin[1] (%s)\n",cmd->scanLevelMin);
//

#endif	// ! _startrun_balls_omp_03_h
