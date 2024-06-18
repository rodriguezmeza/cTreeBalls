// Use:
//ADDONS:
//#include "startrun_balls_omp_03.h"

#ifndef _startrun_balls_omp_03_h
#define _startrun_balls_omp_03_h

if (scanopt(cmd->options, "bodyfound"))
    if (cmd->ntosave < 1 || cmd->ntosave > cmd->nbody)
        error("CheckParameters: absurd value for ntosave\n");
if (cmd->scanLevel < 0)
    error("CheckParameters: absurd value for scanLevel (%d)\n",cmd->scanLevel);
//    if (cmd.scanLevel > 12)
//        error("CheckParameters: too big value for scanLevel (%d)\n",cmd.scanLevel);
// Root nodes:
if (cmd->scanLevelRoot < 0)
    error("CheckParameters: absurd value for scanLevelRoot (%d)\n",cmd->scanLevelRoot);
//    if (cmd.scanLevelMin > 0)
if ((int)gd->scanLevelMin[0] > 0)
    error("CheckParameters: absurd value for scanLevelMin[0] (%s)\n",cmd->scanLevelMin);
//    if (gd.scanLevelMin[1] < 0 || gd.scanLevelMin[1] > 1)
if (gd->scanLevelMin[1] > 0)
    error("CheckParameters: absurd value for scanLevelMin[1] (%s)\n",cmd->scanLevelMin);
//
//if (gd.rsmooth[0] < 0 || gd.rsmoothFlag==FALSE)
//    error("CheckParameters: absurd value for rsmooth (%s)\n",cmd.rsmooth);


#endif	// ! _startrun_balls_omp_03_h
