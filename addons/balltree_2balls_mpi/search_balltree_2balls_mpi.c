/* Distributed FCFC ball-tree dual-node 2PCF and triple-node 3PCF. */

#include "globaldefs.h"
#include "fcfc_balltree_2balls_mpi.h"

#define TREECORR_METHOD_NAME "balltree-2balls-mpi"
#define TREECORR_DISTRIBUTED_ENGINE 1
#define TREECORR_DISTRIBUTED_TASK_OWNED(task) \
    fcfc_balltree_2balls_mpi_task_owned(task)
#define TREECORR_DISTRIBUTED_IS_ROOT() \
    fcfc_balltree_2balls_mpi_is_root()
#define TREECORR_DISTRIBUTED_SIZE() \
    fcfc_balltree_2balls_mpi_size()
#define TREECORR_DISTRIBUTED_CONSENSUS(cmd, status, operation) \
    fcfc_balltree_2balls_mpi_consensus(cmd, status, operation)
#define TREECORR_DISTRIBUTED_REDUCE_REALS(cmd, values, count) \
    fcfc_balltree_2balls_mpi_reduce_reals(cmd, values, count)
#define TREECORR_DISTRIBUTED_REDUCE_INTEGERS(cmd, values, count) \
    fcfc_balltree_2balls_mpi_reduce_integers(cmd, values, count)
#define searchcalc_balltree_2balls_omp searchcalc_balltree_2balls_mpi

#include "../balltree_2balls_omp/search_balltree_2balls_omp.c"
