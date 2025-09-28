// Use:
//#include "startrun_octree_kkk_omp_10.h"

#ifndef _startrun_octree_kkk_omp_10_h
#define _startrun_octree_kkk_omp_10_h

// check if these allocations are needed!!!...
//  (they may have been allocated already in searchcalc_octree_kkk_omp.c,
//      or not needed any more...)

// 2pcf
gd->histNNSubN2pcf = dvector(1,cmd->sizeHistN);
//B kappa Avg Rmin
gd->histNNSubN2pcftotal = dvector(1,cmd->sizeHistN);
//E
//
gd->histN2pcf = dvector(1,cmd->sizeHistN);

bytes_tot_local += 3*cmd->sizeHistN*sizeof(real);

#endif	// ! _startrun_octree_kkk_omp_10_h
