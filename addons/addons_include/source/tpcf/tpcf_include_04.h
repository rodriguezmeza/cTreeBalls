// Use:
//NOLSST:
//#include "tpcf_include_04.h"

#ifndef _tpcf_include_04_h
#define _tpcf_include_04_h

#ifdef ADDONSDEVELOP
#include "tpcf_03.h"
#endif

#ifdef TREEOMP
#include "tpcf_print_tree_normal_omp.h"
#endif
#ifdef TREE3PCFDIRECTOMP
#include "tpcf_print_tree_3pcf_direct_omp.h"
#endif
#ifdef BALLS
#include "tpcf_print_balls_omp.h"
#endif
#ifdef DIRECTMETHOD
#include "tpcf_print_direct_method.h"
#endif
#ifdef KDTREECUTEBOX
#include "tpcf_print_kdtree_cute_box.h"
#endif
#ifdef KDTREECUTE
#include "tpcf_print_kdtree_cute.h"
#endif
#ifdef KDTREE
#include "tpcf_print_kdtree.h"
#endif

#ifdef TREEOMPSINCOSCROSSCORR
#include "tpcf_print_tree_omp_sincos_crosscorr.h"
#endif

/*
 Add your addon item here
 */

#endif	// ! _tpcf_include_04_h
