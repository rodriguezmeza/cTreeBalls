// Use:
//#include "datastruc_defs_include_01.h"

#ifndef _datastruc_defs_include_01_h
#define _datastruc_defs_include_01_h

#ifdef LYAFORESTOMP
#include "datastruc_defs_lya_forest_omp.h"
#endif

#ifdef OCTREE3PCF3DOMP
    INTEGER octree3pcf3d_los_id;
#define Octree3pcf3dLosId(x) (((bodyptr)(x))->octree3pcf3d_los_id)
#endif

/*
 Add your addon item here
 */

#endif	// ! _datastruc_defs_include_01_h
