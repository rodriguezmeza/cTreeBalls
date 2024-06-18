// Use:
//NOLSST:
//#include "startrun_patch.h"

#ifndef _startrun_patch_h
#define _startrun_patch_h

//#define DIRECTSIMPLE            12
//#define DIRECTSIMPLESINCOS      19

if (strcmp(method_str,"direct-simple") == 0)            *method_int = 12;

if (strcmp(method_str,"direct-simple-sincos") == 0)     *method_int = 19;

#endif	// ! _startrun_02_h
