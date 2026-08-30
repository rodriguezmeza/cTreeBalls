/* Distributed LogMultipole traversal over a binary view of the native octree. */

#include "globaldefs.h"
#include "octree_2balls_tree.h"
#include "fcfc_octree_2balls_mpi.h"

#define TREECORR_METHOD_NAME "octree-2balls-mpi"
#define TREECORR_NATIVE_BINARY_VIEW 1
#define TREECORR_DISTRIBUTED_ENGINE 1
#define TREECORR_DISTRIBUTED_TASK_OWNED(task) \
    fcfc_octree_2balls_mpi_task_owned(task)
#define TREECORR_DISTRIBUTED_IS_ROOT() \
    fcfc_octree_2balls_mpi_is_root()
#define TREECORR_DISTRIBUTED_SIZE() \
    fcfc_octree_2balls_mpi_size()
#define TREECORR_DISTRIBUTED_CONSENSUS(cmd, status, operation) \
    fcfc_octree_2balls_mpi_consensus(cmd, status, operation)
#define TREECORR_DISTRIBUTED_REDUCE_REALS(cmd, values, count) \
    fcfc_octree_2balls_mpi_reduce_reals(cmd, values, count)
#define TREECORR_DISTRIBUTED_REDUCE_INTEGERS(cmd, values, count) \
    fcfc_octree_2balls_mpi_reduce_integers(cmd, values, count)
#define TREECORR_PUBLISH_NODE_COUNT(gd, catalog, count) ((void)0)
#define fcfc_balltree_build octree_2balls_tree_build
#define fcfc_balltree_frontier octree_2balls_tree_frontier
#define fcfc_balltree_free octree_2balls_tree_free
#define searchcalc_balltree_2balls_omp searchcalc_octree_2balls_mpi

#ifdef THREEPCFCONVERGENCE
#define TREECORR_TASK_FRONTIER_ENGINE 1
#define TREECORR_LOG_MULTIPOLE_ENGINE 1
#define TREECORR_BODY_PIVOT_LOG_MULTIPOLE 1
#endif

#include "../balltree_2balls_omp/search_balltree_2balls_omp.c"
