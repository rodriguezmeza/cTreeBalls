// Use:
//#include "protodefs_include.h"

#ifndef _protodefs_include_h
#define _protodefs_include_h

#ifdef DIRECTMETHOD
#include "protodefs_direct_method.h"
#endif

#ifdef BALLS
#include "protodefs_balls_omp.h"
#endif

#ifdef PXD
#include "protodefs_pxd.h"
#endif

#ifdef SHEARDIRECTSIMPLE
#include "protodefs_shear_direct_simple.h"
#endif

#ifdef KDTREEOMP
#include "protodefs_kdtree_omp.h"
#endif

#ifdef BALLTREEEXP
#include "protodefs_balltree_exp.h"
#endif

#ifdef OCTREEKKKOMP
#include "protodefs_octree_kkk_omp.h"
#endif

#ifdef OCTREEKKKCEXPOMP
#include "protodefs_octree_kkk_cexp_omp.h"
#endif

#ifdef OCTREEKKKBALLSOMP
#include "protodefs_octree_kkk_balls_omp.h"
#endif

#ifdef OCTREEKKKBALLS4OMP
#include "protodefs_octree_kkk_balls4_omp.h"
#endif

#ifdef TC3PCFDIRECTOMP
#include "protodefs_tc_3pcf_direct_omp.h"
#endif

#ifdef SAVERESTORE
#include "protodefs_save_restore.h"
#endif

/*
 Add your addon item here
 */

#endif	// ! _protodefs_include_h
