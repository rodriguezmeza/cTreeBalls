// Use:
//#include "globaldefs_include_04.h"

#ifndef _globaldefs_include_04_h
#define _globaldefs_include_04_h

#ifdef GADGETIO
#include "cballsio_gadget_00.h"
#endif

#ifndef OCTREESMOOTHING
#ifdef BALLS
#include "globaldefs_balls_omp_04.h"
#endif
#endif

#ifdef OCTREESMOOTHING
#include "globaldefs_octree_smoothing_04.h"
#endif

#ifdef BALLS0357
#include "globaldefs_octree_smoothing_0357_04.h"
#endif

/*
 Add your addon item here
 */

#endif	// ! _globaldefs_include_04_h
