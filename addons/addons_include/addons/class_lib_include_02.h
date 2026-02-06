// Use:
//#include "class_lib_include_02.h"

#ifndef _class_lib_include_02_h
#define _class_lib_include_02_h

#ifdef OCTREESMOOTHING
#include "class_lib_octree_smoothing_02.h"
#endif

#ifdef IOLIB
#include "input_iolib_02.h"
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
