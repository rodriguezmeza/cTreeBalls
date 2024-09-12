// Use:
//#include "cballs_include_02.h"

#ifndef _cballs_include_02_h
#define _cballs_include_02_h

#ifdef BALLS
#include "cballs_balls_omp.h"
#endif

#ifdef DIRECTMETHOD
#include "cballs_direct_method.h"
#endif

#ifdef SHEARDIRECTSIMPLE
#include "cballs_shear_direct_simple.h"
#endif

#ifdef KDTREEOMP
#include "cballs_kdtree_omp.h"
#endif

#ifdef BALLTREEEXP
#include "cballs_balltree_exp.h"
#endif


/*
 Add your addon item here
 */

#endif	// ! _cballs_include_02_h
