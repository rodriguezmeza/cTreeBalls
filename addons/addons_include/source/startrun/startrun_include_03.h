// Use:
//#include "startrun_include_03.h"
//
//  it is included in (the socket):
//      addons/class_lib/input.c
//      in the testParameterFile routine
//

#ifndef _startrun_include_03_h
#define _startrun_include_03_h

//#ifdef BALLS
//#include "startrun_balls_omp_04.h"
//#endif

#ifndef OCTREESMOOTHING
#ifdef BALLS
#include "startrun_balls_omp_03.h"
#endif
#endif

#ifdef OCTREESMOOTHING
#include "startrun_octree_smoothing_03.h"
#endif

#ifdef BALLS0357
#include "startrun_balls_omp_0357_04.h"
#endif

#ifdef IOLIB
#include "startrun_iolib_04.h"
#endif

#ifdef LYAFORESTOMP
#include "startrun_lya_forest_omp_03.h"
#endif

#ifdef SAVERESTORE
#include "startrun_save_restore_04.h"
#endif


/*
 Add your addon item here
 */

#endif	// ! _startrun_include_03_h
