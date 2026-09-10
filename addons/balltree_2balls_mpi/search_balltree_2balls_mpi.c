/* Distributed FCFC ball-tree dual-node 2PCF and triple-node 3PCF. */

#include "globaldefs.h"
#include "fcfc_balltree_2balls_mpi.h"

#define DUAL_NODE_METHOD_NAME "balltree-2balls-mpi"
#define BALLTREE_2BALLS_PRIMARY_FEATURES 1
#define BALLTREE_2BALLS_FULL_FUNCTION searchcalc_balltree_2balls_full_mpi
#define BALLTREE_2BALLS_SEARCH_FUNCTION searchcalc_balltree_2balls_mpi
#define BALLTREE_2BALLS_LEGACY_FUNCTION searchcalc_balltree_mpi
#ifdef BALLS4SCANLEV
#define DUAL_NODE_SCAN_LEVEL_FRONTIER 1
#endif
#define DUAL_NODE_DISTRIBUTED_ENGINE 1
#define DUAL_NODE_DISTRIBUTED_TASK_OWNED(task) \
    fcfc_balltree_2balls_mpi_task_owned(task)
#define DUAL_NODE_DISTRIBUTED_IS_ROOT() \
    fcfc_balltree_2balls_mpi_is_root()
#define DUAL_NODE_DISTRIBUTED_SIZE() \
    fcfc_balltree_2balls_mpi_size()
#define DUAL_NODE_DISTRIBUTED_CONSENSUS(cmd, status, operation) \
    fcfc_balltree_2balls_mpi_consensus(cmd, status, operation)
#define DUAL_NODE_DISTRIBUTED_REDUCE_REALS(cmd, values, count) \
    fcfc_balltree_2balls_mpi_reduce_reals(cmd, values, count)
#define DUAL_NODE_DISTRIBUTED_REDUCE_INTEGERS(cmd, values, count) \
    fcfc_balltree_2balls_mpi_reduce_integers(cmd, values, count)
#include "../balltree_2balls_omp/search_balltree_2balls_omp.c"
