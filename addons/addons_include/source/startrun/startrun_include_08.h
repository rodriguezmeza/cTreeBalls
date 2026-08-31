// Use:
//#include "startrun_include_08.h"

#ifndef _startrun_include_08_h
#define _startrun_include_08_h


//#ifdef BALLS
//#include "startrun_balls_omp_05.h"
//#endif

#ifndef OCTREESMOOTHING
#ifdef BALLS
#include "startrun_balls_omp_08.h"
#endif
#endif

#ifdef OCTREESMOOTHING
#include "startrun_octree_smoothing_08.h"
#endif

#ifdef BALLS0357
#include "startrun_balls_omp_0357_05.h"
#endif

#ifdef IOLIB
#include "startrun_iolib_05.h"
#endif

#if defined(LYAFORESTOMP) || defined(LYAFORESTMPI)
#include "startrun_lya_forest_omp_08.h"
#endif

#ifdef SAVERESTORE
#include "startrun_save_restore_05.h"
#endif


/*
 Add your addon item here
 */

#endif	// ! _startrun_include_08_h
