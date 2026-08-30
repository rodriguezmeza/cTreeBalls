/* Distributed FCFC ball-tree LogMultipole 3PCF traversal. */

#include "globaldefs.h"
#include "fcfc_balltree_2balls_mpi_3pcf.h"

#define TREECORR_TASK_FRONTIER_ENGINE 1
#define TREECORR_LOG_MULTIPOLE_ENGINE 1
#define TREECORR_METHOD_NAME "balltree-2balls-mpi_3pcf"
#define TREECORR_DISTRIBUTED_ENGINE 1
#define TREECORR_DISTRIBUTED_TASK_OWNED(task) \
    fcfc_balltree_2balls_mpi_3pcf_task_owned(task)
#define TREECORR_DISTRIBUTED_IS_ROOT() \
    fcfc_balltree_2balls_mpi_3pcf_is_root()
#define TREECORR_DISTRIBUTED_SIZE() \
    fcfc_balltree_2balls_mpi_3pcf_size()
#define TREECORR_DISTRIBUTED_CONSENSUS(cmd, status, operation) \
    fcfc_balltree_2balls_mpi_3pcf_consensus(cmd, status, operation)
#define TREECORR_DISTRIBUTED_REDUCE_REALS(cmd, values, count) \
    fcfc_balltree_2balls_mpi_3pcf_reduce_reals(cmd, values, count)
#define TREECORR_DISTRIBUTED_REDUCE_INTEGERS(cmd, values, count) \
    fcfc_balltree_2balls_mpi_3pcf_reduce_integers(cmd, values, count)
#define searchcalc_balltree_2balls_omp searchcalc_balltree_2balls_mpi_3pcf

#ifdef TWOPCF
#undef TWOPCF
#endif

#include "../balltree_2balls_omp/search_balltree_2balls_omp.c"
