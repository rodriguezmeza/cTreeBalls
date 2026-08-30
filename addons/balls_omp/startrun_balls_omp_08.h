// Use:
//#include "startrun_balls_omp_08.h"

//
// included in file: addons/addons_include/source/startrun/startrun_include_08.h
//

#ifndef _startrun_balls_omp_08_h
#define _startrun_balls_omp_08_h

//B
//        fprintf(fdout,FMTIL,"ntosave",cmd->ntosave);    // needed?
fprintf(fdout,FMTI,"scanLevel",cmd->scanLevel);
// Root nodes:
fprintf(fdout,FMTI,"scanLevelRoot",cmd->scanLevelRoot);
//fprintf(fdout,FMTT,"scanLevelMin",cmd->scanLevelMin);
fprintf(fdout,FMTI,"scanLevelMin",cmd->scanLevelMin);
//
//E

#endif	// ! _startrun_balls_omp_08_h
