//=============================================================================
//        1          2          3          4        ^ 5          6          7

// Use:
//#include "globaldefs_include.h"

// included in: include/global_data.h

#ifndef _globaldefs_include_h
#define _globaldefs_include_h

#ifndef OCTREESMOOTHING
#ifdef BALLS
//#include "globaldefs_balls_omp_02.h"
#include "globaldefs_balls_omp_03.h"
#endif
#endif

#ifdef OCTREESMOOTHING
#include "globaldefs_octree_smoothing_03.h"
#endif

#ifdef BALLS0357
#include "globaldefs_balls_omp_0357_03.h"
#endif

#ifdef IOLIB
#include "global_data_iolib.h"
#endif

#ifdef CLASSLIB
#include "global_data_include_class.h"
#endif

#if defined(OCTREE3PCF3DOMP) || defined(OCTREE3PCF3DMPI)
    bool octree3pcf3d_los_ids[MAXITEMS];
#endif


/*
 Add your addon item here
 */

/*
 ############################
 #B Addendum of some not important
 #   modules or that are in
 #   development phase
 # Normally they will be switched OFF
 ############################
 */

#ifdef DIRECTMETHODSIMPLELOOPID
#include "global_data_direct_method_simple_loopId.h"
#endif

/*
 #E Addendum of some not important
 ############################
*/

/*
 ############################
 #B Addendum of some important
 #   modules that will not be part
 #       any longer of the public version
 # Normally they will be switched OFF
 ############################
 */

/*
 #E Addendum of some not important
 #  no longer part of the public version
 ############################
*/



#endif	// ! _globaldefs_include_h
