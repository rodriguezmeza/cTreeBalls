// Use:
//#include "protodefs_include.h"
//
//  it is included in (the socket):
//      include/protodefs.h
//      at the end of the file
//

#ifndef _protodefs_include_h
#define _protodefs_include_h

#ifdef OCTREEGGGOMP
#include "protodefs_octree_ggg_omp.h"
#endif

#ifdef KDTREEOMP
#include "protodefs_kdtree_omp.h"
#endif

#ifdef KDTREEBOXOMP
#include "protodefs_kdtree_box_omp.h"
#endif

#ifdef NEIGHBORBOXESOMP
#include "protodefs_neighbor_boxes_omp.h"
#endif

#ifdef KDTREECUTEBOX
#include "protodefs_kdtree_cute_box.h"
#endif

#ifdef PXD
#include "protodefs_pxd.h"
#endif

#ifdef DIRECTMETHODSIMPLE
#include "protodefs_direct_method_simple.h"
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

#ifdef SAVERESTORE
#include "protodefs_save_restore.h"
#endif

#ifdef OCTREEGGGOMPTRIANGLES
#include "protodefs_octree_ggg_omp_triangles.h"
#endif

#ifdef OCTREEKKKBALLS4OMP
#include "protodefs_octree_kkk_balls4_omp.h"
#endif

#ifdef OCTREEKKKBALLS4OMPTRIANGLES
#include "protodefs_octree_kkk_balls4_omp_triangles.h"
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

#ifdef OCTREESINCOSOMP
#include "protodefs_octree_sincos_omp.h"
#endif

#ifdef TREEOMPSINCOS
#include "protodefs_tree_omp_sincos.h"
#endif

#ifdef BALLS
#include "protodefs_balls_omp.h"
#endif

#ifdef OCTREEKKKOMP
#include "protodefs_octree_kkk_omp.h"
#endif

#ifdef OCTREEBOXOMP
#include "protodefs_octree_box_omp.h"
#endif

/*
 #E Addendum of some not important
 #  no longer part of the public version
 ############################
*/

#endif	// ! _protodefs_include_h
