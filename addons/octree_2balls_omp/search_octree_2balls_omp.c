/* TreeCorr traversal over a binary view of the native cTreeBalls octree. */

#include "globaldefs.h"
#include "octree_2balls_tree.h"

#define TREECORR_METHOD_NAME "octree-2balls-omp"
#define TREECORR_NATIVE_BINARY_VIEW 1
#define TREECORR_PUBLISH_NODE_COUNT(gd, catalog, count) ((void)0)
#define fcfc_balltree_build octree_2balls_tree_build
#define fcfc_balltree_frontier octree_2balls_tree_frontier
#define fcfc_balltree_free octree_2balls_tree_free
#define searchcalc_balltree_2balls_omp searchcalc_octree_2balls_omp

#ifdef THREEPCFCONVERGENCE
#define TREECORR_TASK_FRONTIER_ENGINE 1
#define TREECORR_LOG_MULTIPOLE_ENGINE 1
#define TREECORR_BODY_PIVOT_LOG_MULTIPOLE 1
#endif

#include "../balltree_2balls_omp/search_balltree_2balls_omp.c"
