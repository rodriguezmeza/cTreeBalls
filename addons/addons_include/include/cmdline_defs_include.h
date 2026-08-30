//=============================================================================
//        1          2          3          4        ^ 5          6          7

// Use:
//#include "cmdline_defs_include.h"

// included in: include/cmdline_defs.h

#ifndef _cmdline_defs_include_h
#define _cmdline_defs_include_h

#ifdef OCTREESMOOTHING
#include "cmdline_defs_octree_smoothing.h"
#endif

#ifndef OCTREESMOOTHING
#ifdef BALLS
#include "cmdline_defs_balls_omp.h"
#endif
#endif

#ifdef BALLS0357
#include "cmdline_defs_balls_omp_0357.h"
#endif

#ifdef IOLIB
#include "cmdline_defs_iolib.h"
#endif

#ifdef LYAFORESTOMP
#include "cmdline_defs_lya_forest_omp.h"
#endif

#ifdef SAVERESTORE
#include "cmdline_defs_save_restore.h"
#endif

/*
 Add your addon item here
 */

#endif	// ! _cmdline_defs_include_h
