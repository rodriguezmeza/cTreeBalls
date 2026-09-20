/* dual-node-style dual-node and LogMultipole scans over a median KD tree. */

#include "globaldefs.h"
#include "kdtree_2balls_tree.h"
#include "protodefs_kdtree_omp.h"

#ifndef KDTREE_2BALLS_METHOD_NAME
#define KDTREE_2BALLS_METHOD_NAME "kdtree-2balls-omp"
#endif
#ifndef KDTREE_2BALLS_FULL_FUNCTION
#define KDTREE_2BALLS_FULL_FUNCTION searchcalc_kdtree_2balls_full_omp
#endif
#ifndef KDTREE_2BALLS_SEARCH_FUNCTION
#define KDTREE_2BALLS_SEARCH_FUNCTION searchcalc_kdtree_2balls_omp
#endif

static inline int kdtree_2balls_prepare_pivots(
        struct cmdline_data *cmd, struct global_data *gd,
        bodyptr *body_table, INTEGER *body_count,
        INTEGER pivot_minimum, INTEGER *pivot_maximum,
        int pivot_catalog, int neighbor_catalog)
{
#ifdef SMOOTHPIVOT
    return prepare_smooth_pivots(
        cmd, gd, body_table, body_count, pivot_minimum, pivot_maximum,
        pivot_catalog, neighbor_catalog);
#else
    (void)cmd;
    (void)gd;
    (void)body_table;
    (void)body_count;
    (void)pivot_minimum;
    (void)pivot_maximum;
    (void)pivot_catalog;
    (void)neighbor_catalog;
    return SUCCESS;
#endif
}

#define DUAL_NODE_METHOD_NAME KDTREE_2BALLS_METHOD_NAME
#define DUAL_NODE_PUBLISH_NODE_COUNT(gd, catalog, count) ((void)0)
#define DUAL_NODE_PREPARE_PIVOTS kdtree_2balls_prepare_pivots
#define DUAL_NODE_BUILD_PIVOT_TREE kdtree_2balls_tree_build_pivot
#define DUAL_NODE_BUILD_NEIGHBOR_TREE kdtree_2balls_tree_build_neighbor
#define DUAL_NODE_SEPARATE_PIVOT_TREE(cmd) cballs_opt_smooth_pivot(cmd)
#ifdef SMOOTHPIVOT
#define DUAL_NODE_PIVOT_IS_ACTIVE(context, pivot) \
    ((!cballs_opt_read_mask((context)->cmd) \
      || Mask(pivot) != MASK_NODE_MASKED) \
     && (!cballs_opt_smooth_pivot((context)->cmd) || Update(pivot) != FALSE))
#define DUAL_NODE_PIVOT_FIELD(context, pivot) \
    (cballs_opt_smooth_pivot((context)->cmd) \
     ? KappaRmin(pivot) : dual_node_body_weighted_field((context), (pivot)))
#define DUAL_NODE_PIVOT_NORMALIZATION(context, pivot) \
    (cballs_opt_smooth_pivot((context)->cmd) \
     ? (cballs_opt_weights_norm((context)->cmd) \
        ? WeightRmin(pivot) : (real)MAX(NbRmin(pivot), 1)) \
     : dual_node_body_normalization((context), (pivot)))
#else
#define DUAL_NODE_PIVOT_IS_ACTIVE(context, pivot) \
    (!cballs_opt_read_mask((context)->cmd) || Mask(pivot) != MASK_NODE_MASKED)
#endif
#define DUAL_NODE_CONTEXT_WEIGHTED(cmd) \
    (cballs_opt_weights_norm(cmd) || cballs_opt_smooth_pivot(cmd))
#ifdef BALLS4SCANLEV
#define DUAL_NODE_SCAN_LEVEL_FRONTIER 1
#endif
#define fcfc_balltree_frontier kdtree_2balls_tree_frontier
#define fcfc_balltree_free kdtree_2balls_tree_free
#define searchcalc_balltree_2balls_omp KDTREE_2BALLS_FULL_FUNCTION
#define DUAL_NODE_ADAPTIVE_PAIR_LEAVES 1
#define DUAL_NODE_PAIR_BATCH_SIZE 256
#define DUAL_NODE_USE_NATURAL_LOG_BINS 1
#if defined(__APPLE__)
#define DUAL_NODE_USE_ACCELERATE_VFORCE 1
#endif

#ifdef THREEPCFCONVERGENCE
#define DUAL_NODE_TASK_FRONTIER_ENGINE 1
#define DUAL_NODE_LOG_MULTIPOLE_ENGINE 1
#define DUAL_NODE_BODY_PIVOT_LOG_MULTIPOLE 1
#define DUAL_NODE_PERSISTENT_PARTIAL_FRONTIER 1
#ifndef DUAL_NODE_DISTRIBUTED_ENGINE
#define DUAL_NODE_PIVOT_PROGRESS 1
#endif
#endif

#include "../balltree_2balls_omp/search_balltree_2balls_omp.c"

#undef searchcalc_balltree_2balls_omp

global int KDTREE_2BALLS_SEARCH_FUNCTION(
        struct cmdline_data *cmd, struct global_data *gd,
        bodyptr *body_table, INTEGER *body_count,
        INTEGER pivot_minimum, INTEGER *pivot_maximum,
        int pivot_catalog, int neighbor_catalog)
{
#ifdef DUAL_NODE_DISTRIBUTED_ENGINE
    if (cballs_opt_legacy_one_ball(cmd)) {
        if (cballs_opt_no_two_balls(cmd)
            || scanopt(cmd->options, "dual-node-bin-slop")
            || scanopt(cmd->options, "dual-node-direct-triples")) {
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     "%s: legacy-one-ball cannot be combined with "
                     "no-two-balls, dual-node-bin-slop, or "
                     "dual-node-direct-triples",
                     cmd->searchMethod);
            return FAILURE;
        }
        verb_print(cmd->verbose,
                   "%s: dispatching to the distributed kdtree legacy kernel\n",
                   cmd->searchMethod);
        return searchcalc_kdtree_omp(
            cmd, gd, body_table, body_count, pivot_minimum, pivot_maximum,
            pivot_catalog, neighbor_catalog);
    }
#else
    if (cballs_opt_legacy_one_ball(cmd)) {
        if (cballs_opt_no_two_balls(cmd)
            || scanopt(cmd->options, "dual-node-bin-slop")
            || scanopt(cmd->options, "dual-node-direct-triples")) {
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     "%s: legacy-one-ball cannot be combined with "
                     "no-two-balls, dual-node-bin-slop, or "
                     "dual-node-direct-triples",
                     cmd->searchMethod);
            return FAILURE;
        }
        verb_print(cmd->verbose,
                   "%s: dispatching to the kdtree-omp legacy kernel\n",
                   cmd->searchMethod);
        return searchcalc_kdtree_omp(
            cmd, gd, body_table, body_count, pivot_minimum, pivot_maximum,
            pivot_catalog, neighbor_catalog);
    }
#endif
#ifdef SMOOTHPIVOT
    if (cballs_opt_smooth_pivot(cmd)
        && scanopt(cmd->options, "dual-node-direct-triples")) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: dual-node-direct-triples does not support smooth-pivot; "
                 "use no-smooth-pivot for that validation oracle",
                 cmd->searchMethod);
        return FAILURE;
    }
#endif
    return KDTREE_2BALLS_FULL_FUNCTION(
        cmd, gd, body_table, body_count, pivot_minimum, pivot_maximum,
        pivot_catalog, neighbor_catalog);
}
