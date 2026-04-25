// Use:
//#include "globaldefs_include_03.h"

#ifndef _globaldefs_include_03_h
#define _globaldefs_include_03_h

#ifdef BALLS
#include "globaldefs_balls_omp_02.h"
#endif

#ifdef OCTREESMOOTHING
#include "globaldefs_octree_smoothing_03.h"
#endif

#ifdef IOLIB
#include "global_data_iolib.h"
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



#endif	// ! _globaldefs_include_03_h
