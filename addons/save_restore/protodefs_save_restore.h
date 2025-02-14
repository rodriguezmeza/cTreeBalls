// Use:
//#include "protodefs_save_restore.h"

#ifndef _protodefs_save_restore_h
#define _protodefs_save_restore_h

global void checkstop(struct cmdline_data* cmd,
                      struct  global_data* gd);
global void savestate(struct cmdline_data* cmd,
                      struct  global_data* gd,
                      string);
global int restorestate(struct cmdline_data* cmd,
                        struct  global_data* gd, string);
global int PrintState(struct  cmdline_data *cmd,
                      struct  global_data* gd,
                      char *fname);


#endif	// ! _protodefs_save_restore_h
