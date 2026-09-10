/* Distributed LogMultipole traversal over a binary view of the native octree. */

#include "globaldefs.h"
#include "native_octree_pair.h"
#include "octree_2balls_tree.h"
#include "fcfc_octree_2balls_mpi.h"

#define DUAL_NODE_METHOD_NAME "octree-2balls-mpi"
#define DUAL_NODE_NATIVE_BINARY_VIEW 1
#define DUAL_NODE_DISTRIBUTED_ENGINE 1
#define DUAL_NODE_DISTRIBUTED_TASK_OWNED(task) \
    fcfc_octree_2balls_mpi_task_owned(task)
#define DUAL_NODE_DISTRIBUTED_IS_ROOT() \
    fcfc_octree_2balls_mpi_is_root()
#define DUAL_NODE_DISTRIBUTED_SIZE() \
    fcfc_octree_2balls_mpi_size()
#define DUAL_NODE_DISTRIBUTED_CONSENSUS(cmd, status, operation) \
    fcfc_octree_2balls_mpi_consensus(cmd, status, operation)
#define DUAL_NODE_DISTRIBUTED_REDUCE_REALS(cmd, values, count) \
    fcfc_octree_2balls_mpi_reduce_reals(cmd, values, count)
#define DUAL_NODE_DISTRIBUTED_REDUCE_INTEGERS(cmd, values, count) \
    fcfc_octree_2balls_mpi_reduce_integers(cmd, values, count)
#define DUAL_NODE_PUBLISH_NODE_COUNT(gd, catalog, count) ((void)0)
#define fcfc_balltree_build octree_2balls_tree_build
#define fcfc_balltree_frontier octree_2balls_tree_frontier
#define fcfc_balltree_free octree_2balls_tree_free
#define searchcalc_balltree_2balls_omp searchcalc_octree_2balls_full_mpi
#ifdef BALLS4SCANLEV
#define DUAL_NODE_SCAN_LEVEL_FRONTIER 1
#endif

#ifdef THREEPCFCONVERGENCE
#define DUAL_NODE_TASK_FRONTIER_ENGINE 1
#define DUAL_NODE_LOG_MULTIPOLE_ENGINE 1
#define DUAL_NODE_BODY_PIVOT_LOG_MULTIPOLE 1
#endif

#include "../balltree_2balls_omp/search_balltree_2balls_omp.c"

#undef searchcalc_balltree_2balls_omp

static bool octree_2balls_mpi_task_owned(
        struct cmdline_data *cmd, INTEGER task)
{
    (void)cmd;
    return fcfc_octree_2balls_mpi_task_owned(task);
}

static bool octree_2balls_mpi_publish(struct cmdline_data *cmd)
{
    (void)cmd;
    return fcfc_octree_2balls_mpi_is_root();
}

global int searchcalc_octree_2balls_mpi(
        struct cmdline_data *cmd, struct global_data *gd,
        bodyptr *btab, INTEGER *nbody, INTEGER ipmin, INTEGER *ipmax,
        int cat1, int cat2)
{
    if (cballs_opt_legacy_one_ball(cmd)) {
        verb_print(cmd->verbose,
                   "octree-2balls-mpi: dispatching to the distributed "
                   "octree-GGG compatibility kernel\n");
        return searchcalc_octree_ggg_mpi(
            cmd, gd, btab, nbody, ipmin, ipmax, cat1, cat2);
    }

    const cballs_native_pair_parallel parallel = {
        fcfc_octree_2balls_mpi_size(),
        octree_2balls_mpi_task_owned,
        octree_2balls_mpi_publish,
        fcfc_octree_2balls_mpi_consensus,
        fcfc_octree_2balls_mpi_reduce_reals,
        fcfc_octree_2balls_mpi_reduce_integers
    };
    const cballs_native_pair_policy policy = {
        FALSE, FALSE, FALSE, FALSE, &parallel
    };

#ifdef THREEPCFCONVERGENCE
    if (!cballs_opt_only_2pcf(cmd))
        return searchcalc_octree_2balls_full_mpi(
            cmd, gd, btab, nbody, ipmin, ipmax, cat1, cat2);
#endif
    return cballs_native_octree_pair_search(
        cmd, gd, btab, nbody, ipmin, ipmax, cat1, cat2, &policy);
}
