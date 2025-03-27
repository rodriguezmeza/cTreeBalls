// Use:
//#include "cballs_include_03.h"

#ifndef _cballs_include_03_h
#define _cballs_include_03_h

#ifdef BALLS
#include "cballs_print_balls_omp.h"
#endif

#ifdef DIRECTMETHOD
#include "cballs_print_direct_method.h"
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

#endif	// ! _cballs_include_03_h
