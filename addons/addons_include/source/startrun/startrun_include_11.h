// Use:
//#include "startrun_include_11.h"

// These header files have defined tag numbers for the searching methods...
// It is recommended to use tag numbers greater than 100

#ifndef _startrun_include_11_h
#define _startrun_include_11_h

#ifdef ADDONSDEVELOP
#include "startrun_02.h"                        // Here several tags: 1 -- 53,
                                                //  some numbers are missing
                                                //  and used below
#endif

#ifdef BALLS
#include "startrun_balls_omp_06.h"              // 46
#endif

#ifdef DIRECTMETHOD
#include "startrun_direct_method.h"             // 19
#endif

#ifdef SHEARDIRECTSIMPLE
#include "startrun_shear_direct_simple.h"       // 100
#endif

#ifdef KDTREEOMP
#include "startrun_kdtree_omp_11.h"             // 59
#endif

#ifdef BALLTREEEXP
#include "startrun_balltree_exp.h"              // 52
#endif

#ifdef TC3PCFDIRECTOMP
#include "startrun_tc_3pcf_direct_omp_11.h"     // 60
#endif

#ifdef OCTREEKKKOMP
#include "startrun_octree_kkk_omp_11.h"             // 61
#endif

#ifdef OCTREEKKKCEXPOMP
#include "startrun_octree_kkk_cexp_omp_11.h"             // 64
#endif

#ifdef OCTREEKKKBALLSOMP
#include "startrun_octree_kkk_balls_omp_11.h"             // 62
#endif

#ifdef OCTREEKKKBALLS4OMP
#include "startrun_octree_kkk_balls4_omp_11.h"             // 63
#endif

/*
 Add your addon item here
 */

#endif	// ! _startrun_include_11_h
