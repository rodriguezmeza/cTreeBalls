// Use:
//NOLSST:
//#include "startrun_tree_3pcf_direct_omp_02.h"

#ifndef _startrun_tree_3pcf_direct_omp_02_h
#define _startrun_tree_3pcf_direct_omp_02_h

//
//B For 3pcf brute force:
//        gd.deltaTheta = TWOPI/cmd.sizeHistTheta;
        gd->deltaTheta = PI/cmd->sizeHistTheta;
//E

#endif	// ! _startrun_tree_3pcf_direct_omp_02_h
