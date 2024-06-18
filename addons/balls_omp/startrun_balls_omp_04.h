// Use:
//ADDONS:
//#include "startrun_balls_omp_04.h"

#ifndef _startrun_balls_omp_04_h
#define _startrun_balls_omp_04_h

LPName(cmd->ntosave,"ntosave");
IPName(cmd->scanLevel,"scanLevel");
// Root nodes:
    IPName(cmd->scanLevelRoot,"scanLevelRoot");
//        IPName(cmd.scanLevelMin,"scanLevelMin");
    SPName(cmd->scanLevelMin,"scanLevelMin",MAXLENGTHOFSTRSCMD);
//    SPName(cmd.rsmooth,"rsmooth",MAXLENGTHOFSTRSCMD);
//

#endif	// ! _startrun_balls_omp_04_h
