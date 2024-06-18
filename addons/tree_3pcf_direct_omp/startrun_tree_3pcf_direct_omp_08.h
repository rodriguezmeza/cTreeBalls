// Use:
//NOLSST:
//#include "startrun_tree_3pcf_direct_omp_08.h"

#ifndef _startrun_tree_3pcf_direct_omp_08_h
#define _startrun_tree_3pcf_direct_omp_08_h

gd->histXi3pcf = dmatrix3D(1,cmd->sizeHistN,1,cmd->sizeHistN,1,cmd->sizeHistTheta);
bytes_tot_local +=
            (cmd->sizeHistN*cmd->sizeHistN*cmd->sizeHistTheta)*sizeof(real);

#endif	// ! _startrun_tree_3pcf_direct_omp_08_h
