// Use:
//#include "globaldefs_octree_smoothing_03.h"
//
//
// included in file: addons/addons_include/include/global_data_include.h
//

#ifndef _globaldefs_octree_smoothing_03_h
#define _globaldefs_octree_smoothing_03_h

INTEGER nbodysm;

int scanLevelMin[MAXITEMS];
char nodesfilePath[MAXLENGTHOFFILES];
int nnodescanlev;
int nnodescanlev_root;

char bodiesfilePath[MAXLENGTHOFFILES];

bool flagSmoothCellMin;
bool flagSmooth;
bool flagSetNbNoSel;
//B BUCKET
real rminCell[2];               // used only in treeload's scanLevel
                                //  and in search_balls_omp;
                                //      search_balls_kk_omp;
                                //      search_octree_kk_balls4_omp;
                                //      search_octree_kkk_balls4_omp;
//E


#endif	// ! _globaldefs_octree_smoothing_03_h
