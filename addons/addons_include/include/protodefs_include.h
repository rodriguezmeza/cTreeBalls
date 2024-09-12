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

/*
 Add your addon item here
 */

#endif	// ! _protodefs_include_h
