// Use:
//#include "cballs_include_03.h"

#ifndef _cballs_include_03_h
#define _cballs_include_03_h

#ifdef BALLS
#include "cballs_print_balls_omp.h"
#endif

#ifdef KDTREEOMP
#include "cballs_print_kdtree_omp.h"
#endif

#ifdef OCTREEKKKOMP
#include "cballs_print_octree_kkk_omp.h"
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

#ifdef DIRECTMETHOD
#include "cballs_print_direct_method.h"
#endif

/*
 #B Addendum of some not important
 ############################
*/


#endif	// ! _cballs_include_03_h
