/* Distributed dual-node-style scans over the median KD tree. */

#include "globaldefs.h"
#include "fcfc_kdtree_2balls_mpi.h"

#define KDTREE_2BALLS_METHOD_NAME "kdtree-2balls-mpi"
#define KDTREE_2BALLS_FULL_FUNCTION searchcalc_kdtree_2balls_full_mpi
#define KDTREE_2BALLS_SEARCH_FUNCTION searchcalc_kdtree_2balls_mpi
#define DUAL_NODE_DISTRIBUTED_ENGINE 1
#define DUAL_NODE_DISTRIBUTED_DIRECT_TRIPLES 1
#define DUAL_NODE_DISTRIBUTED_TASK_OWNED(task) \
    fcfc_kdtree_2balls_mpi_task_owned(task)
#define DUAL_NODE_DISTRIBUTED_IS_ROOT() \
    fcfc_kdtree_2balls_mpi_is_root()
#define DUAL_NODE_DISTRIBUTED_SIZE() \
    fcfc_kdtree_2balls_mpi_size()
#define DUAL_NODE_DISTRIBUTED_CONSENSUS(cmd, status, operation) \
    fcfc_kdtree_2balls_mpi_consensus((cmd), (status), (operation))
#define DUAL_NODE_DISTRIBUTED_REDUCE_REALS(cmd, values, count) \
    fcfc_kdtree_2balls_mpi_reduce_reals((cmd), (values), (count))
#define DUAL_NODE_DISTRIBUTED_REDUCE_INTEGERS(cmd, values, count) \
    fcfc_kdtree_2balls_mpi_reduce_integers((cmd), (values), (count))

#include "../kdtree_2balls_omp/search_kdtree_2balls_omp.c"
