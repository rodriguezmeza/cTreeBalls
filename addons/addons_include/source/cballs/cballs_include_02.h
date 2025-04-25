// Use:
//#include "cballs_include_02.h"
//
//  it is included in (the socket):
//      source/cballs.c
//      in the EvalHist routine
//

#ifndef _cballs_include_02_h
#define _cballs_include_02_h

#ifdef BALLS
#include "cballs_balls_omp.h"
#endif

#ifdef KDTREEOMP
#include "cballs_kdtree_omp.h"
#endif

#ifdef OCTREEKKKOMP
#include "cballs_octree_kkk_omp.h"
#endif

#ifdef OCTREEGGGOMP
#include "cballs_octree_ggg_omp.h"
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
#include "cballs_direct_method.h"
#endif

#ifdef DIRECTMETHODSIMPLE
#include "cballs_direct_method_simple.h"
#endif


/*
 #E Addendum of some not important
 ############################
*/


#endif	// ! _cballs_include_02_h
