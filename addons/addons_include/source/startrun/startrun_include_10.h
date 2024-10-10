// Use:
//#include "startrun_include_10.h"

#ifndef _startrun_include_10_h
#define _startrun_include_10_h

#ifdef KDTREEOMP
#include "startrun_kdtree_omp_10.h"
#endif

#ifdef OCTREEKKKOMP
#include "startrun_octree_kkk_omp_10.h"
#endif

#ifdef OCTREEKKKBALLSOMP
#include "startrun_octree_kkk_balls_omp_10.h"
#endif

#ifdef OCTREEKKKBALLS4OMP
#include "startrun_octree_kkk_balls4_omp_10.h"
#endif

#ifdef TC3PCFDIRECTOMP
#include "startrun_tc_3pcf_direct_omp_10.h"
#endif

/*
 Add your addon item here
 */

#endif	// ! _startrun_include_10_h
