// Use:
//#include "startrun_include_12.h"

#ifndef _startrun_include_12_h
#define _startrun_include_12_h

#ifndef OCTREESMOOTHING
#ifdef BALLS
#include "startrun_balls_omp_12.h"
#endif
#endif

#ifdef OCTREESMOOTHING
#include "startrun_octree_smoothing_12.h"
#endif

#ifdef BALLS0357
#include "startrun_balls_omp_0357_12.h"
#endif

#ifdef IOLIB
#include "startrun_iolib_06.h"
#endif

/*
 Add your addon item here
 */

#endif	// ! _startrun_include_12_h
