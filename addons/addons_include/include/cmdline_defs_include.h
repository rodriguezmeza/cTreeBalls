// Use:
//NOLSST:
//#include "cmdline_defs_include.h"

#ifndef _cmdline_defs_include_h
#define _cmdline_defs_include_h

#ifdef BALLS
#include "cmdline_defs_balls_omp.h"
#endif
#ifdef TREE3PCFDIRECTOMP
#include "cmdline_defs_tree_3pcf_direct_omp.h"
#endif
#ifdef SAVERESTORE
#include "cmdline_defs_save_restore.h"
#endif

#endif	// ! _cmdline_defs_include_h
