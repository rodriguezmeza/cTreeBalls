// Use:
//#include "cballs_include_02.h"
//
//  it is included in (the socket):
//      source/cballs.c
//      in the EvalHist routine
//

#ifndef _cballs_include_02_h
#define _cballs_include_02_h

#ifdef OCTREESINCOSOMP
#include "cballs_octree_sincos_omp.h"
#endif

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

#ifdef KDTREECUTEBOX
#include "cballs_kdtree_cute_box.h"
#endif

#ifdef OCTREEGGGOMPTRIANGLES
#include "cballs_octree_ggg_omp_triangles.h"
#endif

#ifdef OCTREEKKKBALLS4OMP
#include "cballs_octree_kkk_balls4_omp.h"
#endif

#ifdef OCTREEKKKBALLS4OMPTRIANGLES
#include "cballs_octree_kkk_balls4_omp_triangles.h"
#endif

#ifdef KDTREEBOXOMP
#include "cballs_kdtree_box_omp.h"
#endif

#ifdef OCTREEBOXOMP
#include "cballs_octree_box_omp.h"
#endif

#ifdef NEIGHBORBOXESOMP
#include "cballs_neighbor_boxes_omp.h"
#endif


/*
 #E Addendum of some not important
 ############################
*/


#endif	// ! _cballs_include_02_h
