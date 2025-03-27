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

#ifdef SAVERESTORE
#include "protodefs_save_restore.h"
#endif

#ifdef KDTREEOMP
#include "protodefs_kdtree_omp.h"
#endif

#ifdef OCTREEKKKOMP
#include "protodefs_octree_kkk_omp.h"
#endif

#ifdef OCTREEGGGOMP
#include "protodefs_octree_ggg_omp.h"
#endif

/*
 Add your addon item here
 */

#endif	// ! _protodefs_include_h
