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

#ifdef IOLIB
#include "protodefs_iolib.h"
#endif


/*
 Add your addon item here
 */

#endif	// ! _protodefs_include_h
