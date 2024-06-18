// Use:
//NOLSST:
//#include "protodefs_include.h"

#ifndef _protodefs_include_h
#define _protodefs_include_h

#ifdef ADDONSDEVELOP
#include "protodefs_01.h"
#include "protodefs_02.h"
#endif
#ifdef DIRECTMETHOD
#include "protodefs_direct_method.h"
#endif
#ifdef TREEOMP
#include "protodefs_tree_normal_omp.h"
#endif
#ifdef TREE3PCFDIRECTOMP
#include "protodefs_tree_3pcf_direct_omp.h"
#endif
#ifdef BALLS
#include "protodefs_balls_omp.h"
#endif
#ifdef KDTREECUTEBOX
#include "protodefs_kdtree_cute_box.h"
#endif
#ifdef KDTREECUTE
#include "protodefs_kdtree_cute.h"
#endif
#ifdef KDTREE
#include "protodefs_kdtree.h"
#endif

#ifdef TREEOMPSINCOSCROSSCORR
#include "protodefs_tree_omp_sincos_crosscorr.h"
#endif

/*
 Add your addon item here
 */

#endif	// ! _protodefs_include_h
