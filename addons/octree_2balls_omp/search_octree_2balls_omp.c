/* dual-node traversal over a binary view of the native cTreeBalls octree. */

#include "globaldefs.h"
#include "native_octree_pair.h"
#include "octree_2balls_tree.h"

#define DUAL_NODE_METHOD_NAME "octree-2balls-omp"
#define DUAL_NODE_NATIVE_BINARY_VIEW 1
#define DUAL_NODE_PUBLISH_NODE_COUNT(gd, catalog, count) ((void)0)
#define fcfc_balltree_build octree_2balls_tree_build
#define fcfc_balltree_frontier octree_2balls_tree_frontier
#define fcfc_balltree_free octree_2balls_tree_free
#define searchcalc_balltree_2balls_omp searchcalc_octree_2balls_full_omp
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

global int searchcalc_octree_2balls_omp(
        struct cmdline_data *cmd, struct global_data *gd,
        bodyptr *btab, INTEGER *nbody, INTEGER ipmin, INTEGER *ipmax,
        int cat1, int cat2)
{
    const cballs_native_pair_policy policy = {
        FALSE, FALSE, FALSE, FALSE, NULL
    };

    if (cballs_opt_legacy_one_ball(cmd)) {
        verb_print(cmd->verbose,
                   "octree-2balls-omp: dispatching to the octree-GGG "
                   "compatibility kernel\n");
        return searchcalc_octree_ggg_omp(
            cmd, gd, btab, nbody, ipmin, ipmax, cat1, cat2);
    }

#ifdef THREEPCFCONVERGENCE
    if (!cballs_opt_only_2pcf(cmd))
        return searchcalc_octree_2balls_full_omp(
            cmd, gd, btab, nbody, ipmin, ipmax, cat1, cat2);
#endif
    return cballs_native_octree_pair_search(
        cmd, gd, btab, nbody, ipmin, ipmax, cat1, cat2, &policy);
}
