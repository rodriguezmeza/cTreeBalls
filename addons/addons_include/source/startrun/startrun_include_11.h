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

#ifdef KDTREEOMP
#include "startrun_kdtree_omp_11.h"             // 59
#endif

#ifdef OCTREEKKKOMP
#include "startrun_octree_kkk_omp_11.h"             // 61
#endif

#ifdef OCTREEGGGOMP
#include "startrun_octree_ggg_omp_11.h"             // 66
#endif

/*
 Add your addon item here
 */

// it is recommended to use tag numbers above 100...


#endif	// ! _startrun_include_11_h
