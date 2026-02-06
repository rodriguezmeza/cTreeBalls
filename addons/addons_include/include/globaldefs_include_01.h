// Use:
//#include "globaldefs_include_01.h"

#ifndef _globaldefs_include_01_h
#define _globaldefs_include_01_h

//#ifdef BALLS
//#include "globaldefs_balls_omp_01.h"
//#endif

#ifdef OCTREESMOOTHING
#include "cmdline_data_octree_smoothing_01.h"
#endif

#ifdef IOLIB
#include "cmdline_data_iolib.h"
#endif

#ifdef SAVERESTORE
#include "cmdline_data_save_restore.h"
#endif


/*
 Add your addon item here
 */

#endif	// ! _globaldefs_include_01_h
