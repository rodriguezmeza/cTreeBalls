// Use:
//#include "startrun_octree_smoothing_07.h"

//
// included in file: addons/addons_include/source/startrun/startrun_include_07.h
//

#ifndef _startrun_balls_omp_07_h
#define _startrun_balls_omp_07_h


if (cmd->scanLevel < 0)
cBALLS_FAIL(cmd, "CheckParameters: absurd value for scanLevel (%d)\n",cmd->scanLevel);
if (cmd->scanLevelRoot < 0)
cBALLS_FAIL(cmd, "CheckParameters: absurd value for scanLevelRoot (%d)\n",cmd->scanLevelRoot);
if (cmd->scanLevelMin >= 1)
cBALLS_FAIL(cmd, "CheckParameters: absurd value for scanLevelMin (%d)\n",
            cmd->scanLevelMin);

#endif	// ! _startrun_balls_omp_07_h
