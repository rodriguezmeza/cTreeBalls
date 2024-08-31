// Use:
//#include "startrun_balls_omp_02.h"

#ifndef _startrun_balls_omp_02_h
#define _startrun_balls_omp_02_h

if (GetParamStat("ntosave") & ARGPARAM)
    cmd->ntosave = GetlParam("ntosave");
if (GetParamStat("scanLevel") & ARGPARAM)
        cmd->scanLevel = GetiParam("scanLevel");
// Root nodes:
    if (GetParamStat("scanLevelRoot") & ARGPARAM)
            cmd->scanLevelRoot = GetiParam("scanLevelRoot");
    if (GetParamStat("scanLevelMin") & ARGPARAM)
            cmd->scanLevelMin = GetParam("scanLevelMin");
//

#endif	// ! _startrun_balls_omp_02_h
