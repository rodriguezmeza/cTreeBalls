// Use:
//NOLSST:
//#include "globaldefs_include_01.h"

#ifndef _globaldefs_include_01_h
#define _globaldefs_include_01_h

#ifdef BALLS
#include "globaldefs_balls_omp_01.h"
#endif
//#else // ! ADDONS
//    string scanLevelMin;
#ifdef TREE3PCFDIRECTOMP
#include "globaldefs_tree_3pcf_direct_omp.h"
#endif
#ifdef SAVERESTORE
#include "globaldefs_save_restore.h"
#endif

#endif	// ! _globaldefs_include_01_h
