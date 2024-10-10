// Use:
//#include "global_data_include.h"

#ifndef _global_data_include_h
#define _global_data_include_h

#include "global_data_include_NMultipoles.h"

//B Delete all below
#ifdef KDTREEOMP
//#include "global_data_kdtree_omp.h"
#endif

#ifdef OCTREEKKKOMP
//#include "global_data_octree_kkk_omp.h"
#endif

#ifdef OCTREEKKKBALLSOMP
//#include "global_data_octree_kkk_balls_omp.h"
#endif

#ifdef OCTREEKKKBALLS4OMP
//#include "global_data_octree_kkk_balls4_omp.h"
#endif

//E



/*
 Add your addon item here
 */

#endif	// ! _global_data_include_h
