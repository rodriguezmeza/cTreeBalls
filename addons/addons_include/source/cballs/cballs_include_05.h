// Use:
//#include "cballs_include_05.h"

#ifndef _cballs_include_05_h
#define _cballs_include_05_h

#ifdef PXD
#include "cballs_pxd_05.h"
#endif

#include "cballs_include_05_NMultipoles.h"

//B Delete includes below
#ifdef KDTREEOMP
//#include "cballs_kdtree_omp_05.h"
#endif

#ifdef OCTREEKKKOMP
#include "cballs_octree_kkk_omp_05.h"
#endif


#ifdef OCTREEKKKBALLSOMP
//#include "cballs_octree_kkk_balls_omp_05.h"
#endif

#ifdef OCTREEKKKBALLS4OMP
//#include "cballs_octree_kkk_balls4_omp_05.h"
#endif
//E

/*
 Add your addon item here
 */

#endif	// ! _cballs_include_05_h
