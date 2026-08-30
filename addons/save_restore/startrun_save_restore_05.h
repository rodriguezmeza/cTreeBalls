// Use:
//NOLSST:
//#include "startrun_save_restore_05.h"

#ifndef _startrun_save_restore_05_h
#define _startrun_save_restore_05_h

WRITE_OUTPUT_OR_FAIL(fdout, buf, FMTT,"statefile",cmd->statefile);
WRITE_OUTPUT_OR_FAIL(fdout, buf, FMTT,"restorefile",cmd->restorefile);

#endif	// ! _startrun_save_restore_05_h
