// Use:
//#include "protodefs_include.h"
//
//  it is included in (the socket):
//      include/protodefs.h
//      at the end of the file
//

#ifndef _protodefs_include_h
#define _protodefs_include_h

#ifdef BALLS
#include "protodefs_balls_omp.h"
#endif

#ifdef PXD
#include "protodefs_pxd.h"
#endif

#ifdef KDTREEOMP
#include "protodefs_kdtree_omp.h"
#endif

#ifdef OCTREEKKKOMP
#include "protodefs_octree_kkk_omp.h"
#endif

#ifdef OCTREEGGGOMP
#include "protodefs_octree_ggg_omp.h"
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
#include "protodefs_direct_method.h"
#endif

#ifdef DIRECTMETHODSIMPLE
#include "protodefs_direct_method_simple.h"
#endif

#ifdef SAVERESTORE
#include "protodefs_save_restore.h"
#endif

#ifdef KDTREECUTEBOX
#include "protodefs_kdtree_cute_box.h"
#endif

/*
 #E Addendum of some not important
 ############################
*/

#endif	// ! _protodefs_include_h
