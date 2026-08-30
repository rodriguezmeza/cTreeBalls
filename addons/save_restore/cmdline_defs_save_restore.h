// Use:
//#include "cmdline_defs_save_restore.h"

#ifndef _cmdline_defs_save_restore_h
#define _cmdline_defs_save_restore_h

//
// When activate, choose the above explanation instead of the second one
//
//    "stepState=10000",                  ";number of steps to save a state-run info; or to save a state of the run which is not working yet",
//    "stepState=10000",                  ";number of steps to save a state-run info (pivot number completed in the log file)",

//B Work in progress on this section
    "statefile=",                  ";Write run state to a file; use in combination with stepState options above", ":state",
    "restorefile=",                     ";Continue run from state file", ":restore",
//E

#endif	// ! _cmdline_defs_save_restore_h
