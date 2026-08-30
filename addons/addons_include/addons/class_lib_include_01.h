// Use:
//#include "class_lib_include_01.h"
//
//  it is included in (the socket):
//      addons/class_lib/input.c
//      in the input_read_parameters_general routine
//

#ifndef _class_lib_include_01_h
#define _class_lib_include_01_h

#ifdef OCTREESMOOTHING
#include "class_lib_octree_smoothing_01.h"
#endif

#ifdef BALLS0357
#include "input_balls_omp_0357_01.h"
#endif

#ifdef IOLIB
#include "input_iolib_01.h"
#endif

#ifdef LYAFORESTOMP
#include "input_lya_forest_omp_01.h"
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
#include "input_save_restore_01.h"
#endif


/*
 #B Addendum of some not important
 ############################
*/


#endif	// ! _class_lib_include_01_h
