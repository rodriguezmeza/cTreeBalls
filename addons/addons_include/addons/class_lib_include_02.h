// Use:
//#include "class_lib_include_02.h"
//
//  it is included in (the socket):
//      addons/class_lib/input.c
//      in the input_default_params routine
//

#ifndef _class_lib_include_02_h
#define _class_lib_include_02_h

#ifdef OCTREESMOOTHING
#include "class_lib_octree_smoothing_02.h"
#endif

#ifdef BALLS0357
#include "input_balls_omp_0357_02.h"
#endif

#ifdef IOLIB
#include "input_iolib_02.h"
#endif

#if defined(LYAFORESTOMP) || defined(LYAFORESTMPI)
#include "input_lya_forest_omp_02.h"
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

#ifdef SAVERESTORE
#include "input_save_restore_02.h"
#endif



/*
 #B Addendum of some not important
 ############################
*/


#endif	// ! _class_lib_include_02_h
