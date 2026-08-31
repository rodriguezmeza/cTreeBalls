//=============================================================================
//        1          2          3          4        ^ 5          6          7

// Use:
//#include "globaldefs_include_01.h"

// included in: include/cmdline_data.h

#ifndef _globaldefs_include_01_h
#define _globaldefs_include_01_h

//#ifdef BALLS
//#include "globaldefs_balls_omp_01.h"
//#endif

#ifdef OCTREESMOOTHING
#include "cmdline_data_octree_smoothing_01.h"
#endif

#ifdef BALLS0357
#include "globaldefs_balls_omp_0357_01.h"
#endif

#ifndef OCTREESMOOTHING
#ifdef BALLS
#include "cmdline_data_balls_omp.h"
#endif
#endif

#ifdef IOLIB
#include "cmdline_data_iolib.h"
#endif

#if defined(LYAFORESTOMP) || defined(LYAFORESTMPI)
#include "cmdline_data_lya_forest_omp.h"
#endif

#ifdef SAVERESTORE
#include "cmdline_data_save_restore.h"
#endif

#ifdef CLASSLIB
#include "cmdline_data_include_class.h"
#endif

/*
 Add your addon item here
 */

#endif	// ! _globaldefs_include_01_h
