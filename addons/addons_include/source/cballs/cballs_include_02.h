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

#ifdef TREEOMPSINCOS
#include "cballs_tree_omp_sincos.h"
#endif

#ifdef BALLS
#include "cballs_balls_omp.h"
#endif

#ifdef BALLS0357
#include "cballs_balls_omp_0357.h"
#endif

#ifdef KDTREEOMP
#include "cballs_kdtree_omp.h"
#endif

#if defined(BALLTREEOMP) || defined(BALLTREEMPI)
#include "cballs_balltree_omp.h"
#endif

#ifdef BALLTREE2BALLSOMP
#include "cballs_balltree_2balls_omp.h"
#endif

#ifdef BALLTREE2BALLSMPI
#include "cballs_balltree_2balls_mpi.h"
#endif

#ifdef OCTREE2BALLSOMP
#include "cballs_octree_2balls_omp.h"
#endif

#ifdef OCTREE2BALLSMPI
#include "cballs_octree_2balls_mpi.h"
#endif

#ifdef BALLTREE2BALLSOMP3PCF
#include "cballs_balltree_2balls_omp_3pcf.h"
#endif

#ifdef BALLTREE2BALLSMPI3PCF
#include "cballs_balltree_2balls_mpi_3pcf.h"
#endif

#ifdef OCTREEKKKOMP
#include "cballs_octree_kkk_omp.h"
#endif

#ifdef OCTREEGGGOMP
#include "cballs_octree_ggg_omp.h"
#endif

#ifdef OCTREEGGGMPI
#include "cballs_octree_ggg_mpi.h"
#endif

#ifdef OCTREESHEAROMP
#include "cballs_octree_shear_omp.h"
#endif

#ifdef OCTREEGGGCROSSOMP
#include "cballs_octree_ggg_cross_omp.h"
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

#ifdef OCTREEGGG
#include "cballs_octree_ggg.h"
#endif

#ifdef DIRECTMETHODSIMPLE
#include "cballs_direct_method_simple.h"
#endif

#ifdef DIRECTMETHODSIMPLELOOPID
#include "cballs_direct_method_simple_loopId.h"
#endif

#ifdef KDTREECUTEBOX
#include "cballs_kdtree_cute_box.h"
#endif

#ifdef OCTREEGGGOMPTRIANGLES
#include "cballs_octree_ggg_omp_triangles.h"
#endif

#ifdef OCTREEBALLS4OMP
#include "cballs_octree_balls4_omp.h"
#endif
#ifdef OCTREEBALLS4MPI
#include "cballs_octree_balls4_mpi.h"
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

#if defined(OCTREE3PCF3DOMP) || defined(OCTREE3PCF3DMPI)
#include "cballs_octree_3pcf_3d_omp.h"
#endif

#ifdef LYAFORESTOMP
#include "cballs_lya_forest_omp.h"
#endif
#ifdef LYAFORESTMPI
#include "cballs_lya_forest_mpi.h"
#endif


/*
 #E Addendum of some not important
 ############################
*/


#endif	// ! _cballs_include_02_h
