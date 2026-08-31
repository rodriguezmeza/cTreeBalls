// Use:
//#include "datastruc_defs_include_01.h"

#ifndef _datastruc_defs_include_01_h
#define _datastruc_defs_include_01_h

#if defined(LYAFORESTOMP) || defined(LYAFORESTMPI)
#include "datastruc_defs_lya_forest_omp.h"
#endif

#if defined(OCTREE3PCF3DOMP) || defined(OCTREE3PCF3DMPI)
    INTEGER octree3pcf3d_los_id;
#define Octree3pcf3dLosId(x) (((bodyptr)(x))->octree3pcf3d_los_id)
#endif

/*
 Add your addon item here
 */

#endif	// ! _datastruc_defs_include_01_h
