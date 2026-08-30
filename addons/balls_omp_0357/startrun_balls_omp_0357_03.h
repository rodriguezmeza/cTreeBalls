// Use:
//#include "startrun_balls_omp_0357_03.h"
//
//
// included in file: addons/addons_include/source/startrun/startrun_include_07.h
//

#ifndef _startrun_balls_omp_0357_03_h
#define _startrun_balls_omp_0357_03_h

/*
//if (scanopt(cmd->options, "bodyfound"))
//    if (cmd->ntosave < 1 || cmd->ntosave > cmd->nbody)
//        error("CheckParameters: absurd value for ntosave\n");
if (cmd->scanLevel < 0)
    error("CheckParameters: absurd value for scanLevel (%d)\n",cmd->scanLevel);
// Root nodes:
if (cmd->scanLevelRoot < 0)
    error("CheckParameters: absurd value for scanLevelRoot (%d)\n",cmd->scanLevelRoot);
if ((int)gd->scanLevelMin[0] > 0)
    error("CheckParameters: absurd value for scanLevelMin[0] (%s)\n",cmd->scanLevelMin);
if (gd->scanLevelMin[1] > 0)
    error("CheckParameters: absurd value for scanLevelMin[1] (%s)\n",cmd->scanLevelMin);
*/
//

if (cmd->scanLevel < 0)
cBALLS_FAIL(cmd, "CheckParameters: absurd value for scanLevel (%d)\n",cmd->scanLevel);
if (cmd->scanLevelRoot < 0)
cBALLS_FAIL(cmd, "CheckParameters: absurd value for scanLevelRoot (%d)\n",cmd->scanLevelRoot);
if (cmd->scanLevelMin >= 1)
cBALLS_FAIL(cmd, "CheckParameters: absurd value for scanLevelMin (%d)\n",
            cmd->scanLevelMin);
//if (gd->scanLevelMin[1] >= 1)
//cBALLS_FAIL(cmd, "CheckParameters: absurd value for scanLevelMin[1] (%d)\n",
//            gd->scanLevelMin[1]);


#endif	// ! _startrun_balls_omp_0357_03_h
