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

#ifdef SHEARDIRECTSIMPLE
#include "cballs_print_shear_direct_simple.h"
#endif

#ifdef KDTREEOMP
#include "cballs_print_kdtree_omp.h"
#endif

#ifdef BALLTREEEXP
#include "cballs_print_balltree_exp.h"
#endif

#ifdef OCTREEKKKOMP
#include "cballs_print_octree_kkk_omp.h"
#endif

#ifdef OCTREEKKKBALLSOMP
#include "cballs_print_octree_kkk_balls_omp.h"
#endif

#ifdef OCTREEKKKBALLS4OMP
#include "cballs_print_octree_kkk_balls4_omp.h"
#endif

#ifdef TC3PCFDIRECTOMP
#include "cballs_print_tc_3pcf_direct_omp.h"
#endif

/*
 Add your addon item here
 */

#endif	// ! _cballs_include_03_h
