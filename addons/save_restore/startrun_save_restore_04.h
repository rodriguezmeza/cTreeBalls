// Use:
//NOLSST:
//#include "startrun_save_restore_04.h"

#ifndef _startrun_save_restore_04_h
#define _startrun_save_restore_04_h

SPName(cmd->statefile,"statefile",MAXLENGTHOFSTRSCMD);
SPName(cmd->restorefile,"restorefile",MAXLENGTHOFSTRSCMD);
//B new
gd->statefileFlag = TRUE;
gd->restorefileFlag = TRUE;
//E
#endif	// ! _startrun_save_restore_04_h
