// Use:
//#include "startrun_include_07.h"

#ifndef _startrun_include_07_h
#define _startrun_include_07_h

#ifndef OCTREESMOOTHING
#ifdef BALLS
#include "startrun_balls_omp_07.h"
#endif
#endif

#ifdef OCTREESMOOTHING
#include "startrun_octree_smoothing_07.h"
#endif

#ifdef BALLS0357
#include "startrun_balls_omp_0357_03.h"
#endif

#ifdef IOLIB
#include "startrun_iolib_03.h"
#endif

#if defined(LYAFORESTOMP) || defined(LYAFORESTMPI)
#include "startrun_lya_forest_omp_07.h"
#endif


/*
 Add your addon item here
 */

#endif	// ! _startrun_include_07_h
