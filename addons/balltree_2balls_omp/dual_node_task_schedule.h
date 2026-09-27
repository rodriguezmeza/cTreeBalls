/* Work-estimated deterministic frontier and scheduling policy. */
#ifdef DUAL_NODE_TASK_FRONTIER_ENGINE
enum {
    DUAL_NODE_TASK_AUTO3 = 0,
    DUAL_NODE_TASK_TWO_ONE = 1,
    DUAL_NODE_TASK_THREE = 2,
    DUAL_NODE_TASK_CROSS12 = 3
};

typedef struct {
    int kind;
    fcfc_balltreeptr trees[3];
    INTEGER nodes[3];
    unsigned orientation_mask;
    long double estimated_work;
} dual_node_triple_task;

static unsigned dual_node_popcount(unsigned value)
{
    unsigned count = 0;

    while (value != 0U) {
        count += value & 1U;
        value >>= 1;
    }
    return count;
}

static long double dual_node_task_work(const dual_node_triple_task *task)
{
    const long double n0 = (long double)dual_node_node_count(
        &task->trees[0]->nodes[task->nodes[0]]);

    switch (task->kind) {
        case DUAL_NODE_TASK_AUTO3:
            return n0 >= 3.0L ? n0 * (n0 - 1.0L) * (n0 - 2.0L) : 0.0L;
        case DUAL_NODE_TASK_TWO_ONE: {
            const long double n1 = (long double)dual_node_node_count(
                &task->trees[0]->nodes[task->nodes[1]]);
            return n0 >= 2.0L
                ? 3.0L * n0 * (n0 - 1.0L) * n1 : 0.0L;
        }
        case DUAL_NODE_TASK_THREE: {
            const long double n1 = (long double)dual_node_node_count(
                &task->trees[1]->nodes[task->nodes[1]]);
            const long double n2 = (long double)dual_node_node_count(
                &task->trees[2]->nodes[task->nodes[2]]);
            return n0 * n1 * n2
                * (long double)dual_node_popcount(task->orientation_mask);
        }
        case DUAL_NODE_TASK_CROSS12: {
            const long double n1 = (long double)dual_node_node_count(
                &task->trees[1]->nodes[task->nodes[1]]);
            return n1 >= 2.0L ? n0 * n1 * (n1 - 1.0L) : 0.0L;
        }
        default:
            return 0.0L;
    }
}

static dual_node_triple_task dual_node_make_task(
        int kind, fcfc_balltreeptr tree0, INTEGER node0,
        fcfc_balltreeptr tree1, INTEGER node1,
        fcfc_balltreeptr tree2, INTEGER node2,
        unsigned orientation_mask)
{
    dual_node_triple_task task;

    task.kind = kind;
    task.trees[0] = tree0;
    task.trees[1] = tree1;
    task.trees[2] = tree2;
    task.nodes[0] = node0;
    task.nodes[1] = node1;
    task.nodes[2] = node2;
    task.orientation_mask = orientation_mask;
    task.estimated_work = 0.0L;
    task.estimated_work = dual_node_task_work(&task);
    return task;
}

static int dual_node_split_task(const dual_node_triple_task *task,
                               dual_node_triple_task children[4])
{
    const fcfc_ballnode *node0 =
        &task->trees[0]->nodes[task->nodes[0]];

    switch (task->kind) {
        case DUAL_NODE_TASK_AUTO3:
            if (dual_node_node_is_leaf(node0)) return 0;
            children[0] = dual_node_make_task(
                DUAL_NODE_TASK_AUTO3, task->trees[0], node0->left,
                NULL, 0, NULL, 0, 0U);
            children[1] = dual_node_make_task(
                DUAL_NODE_TASK_AUTO3, task->trees[0], node0->right,
                NULL, 0, NULL, 0, 0U);
            children[2] = dual_node_make_task(
                DUAL_NODE_TASK_TWO_ONE, task->trees[0], node0->left,
                task->trees[0], node0->right, NULL, 0, 0U);
            children[3] = dual_node_make_task(
                DUAL_NODE_TASK_TWO_ONE, task->trees[0], node0->right,
                task->trees[0], node0->left, NULL, 0, 0U);
            return 4;

        case DUAL_NODE_TASK_TWO_ONE: {
            const fcfc_ballnode *one =
                &task->trees[0]->nodes[task->nodes[1]];

            if (dual_node_node_count(node0) < 2
                || (dual_node_node_is_leaf(node0)
                    && dual_node_node_is_leaf(one)))
                return 0;
            if (!dual_node_node_is_leaf(one)
                && (dual_node_node_is_leaf(node0)
                    || one->radius > node0->radius)) {
                children[0] = dual_node_make_task(
                    DUAL_NODE_TASK_TWO_ONE,
                    task->trees[0], task->nodes[0],
                    task->trees[0], one->left, NULL, 0, 0U);
                children[1] = dual_node_make_task(
                    DUAL_NODE_TASK_TWO_ONE,
                    task->trees[0], task->nodes[0],
                    task->trees[0], one->right, NULL, 0, 0U);
                return 2;
            }
            children[0] = dual_node_make_task(
                DUAL_NODE_TASK_TWO_ONE, task->trees[0], node0->left,
                task->trees[0], task->nodes[1], NULL, 0, 0U);
            children[1] = dual_node_make_task(
                DUAL_NODE_TASK_TWO_ONE, task->trees[0], node0->right,
                task->trees[0], task->nodes[1], NULL, 0, 0U);
            children[2] = dual_node_make_task(
                DUAL_NODE_TASK_THREE, task->trees[0], node0->left,
                task->trees[0], node0->right,
                task->trees[0], task->nodes[1], 0x3fU);
            return 3;
        }

        case DUAL_NODE_TASK_CROSS12: {
            const fcfc_ballnode *pair =
                &task->trees[1]->nodes[task->nodes[1]];

            if (dual_node_node_count(pair) < 2
                || (dual_node_node_is_leaf(node0)
                    && dual_node_node_is_leaf(pair)))
                return 0;
            if (!dual_node_node_is_leaf(node0)
                && (dual_node_node_is_leaf(pair)
                    || node0->radius > pair->radius)) {
                children[0] = dual_node_make_task(
                    DUAL_NODE_TASK_CROSS12, task->trees[0], node0->left,
                    task->trees[1], task->nodes[1], NULL, 0, 0U);
                children[1] = dual_node_make_task(
                    DUAL_NODE_TASK_CROSS12, task->trees[0], node0->right,
                    task->trees[1], task->nodes[1], NULL, 0, 0U);
                return 2;
            }
            children[0] = dual_node_make_task(
                DUAL_NODE_TASK_CROSS12, task->trees[0], task->nodes[0],
                task->trees[1], pair->left, NULL, 0, 0U);
            children[1] = dual_node_make_task(
                DUAL_NODE_TASK_CROSS12, task->trees[0], task->nodes[0],
                task->trees[1], pair->right, NULL, 0, 0U);
            children[2] = dual_node_make_task(
                DUAL_NODE_TASK_THREE, task->trees[0], task->nodes[0],
                task->trees[1], pair->left,
                task->trees[1], pair->right, 0x03U);
            return 3;
        }

        case DUAL_NODE_TASK_THREE: {
            int split = -1;
            real largest = -1.0;

            for (int slot = 0; slot < 3; slot++) {
                const fcfc_ballnode *ball_node =
                    &task->trees[slot]->nodes[task->nodes[slot]];
                if (!dual_node_node_is_leaf(ball_node)
                    && (real)ball_node->radius > largest) {
                    split = slot;
                    largest = (real)ball_node->radius;
                }
            }
            if (split < 0) return 0;
            for (int child = 0; child < 2; child++) {
                INTEGER nodes[3] = {
                    task->nodes[0], task->nodes[1], task->nodes[2]
                };
                const fcfc_ballnode *ball_node =
                    &task->trees[split]->nodes[task->nodes[split]];
                nodes[split] = child == 0
                    ? ball_node->left : ball_node->right;
                children[child] = dual_node_make_task(
                    DUAL_NODE_TASK_THREE,
                    task->trees[0], nodes[0],
                    task->trees[1], nodes[1],
                    task->trees[2], nodes[2],
                    task->orientation_mask);
            }
            return 2;
        }
        default:
            return 0;
    }
}

static INTEGER dual_node_task_target(struct cmdline_data *cmd,
                                    size_t stride, int orders)
{
    size_t values;
    size_t bytes;
    INTEGER target = DUAL_NODE_TRIPLE_TASK_TARGET;

    if (!dual_node_triple_values(stride, orders, dual_node_window_orders(cmd),
                                &values))
        return 1;
    if (values > SIZE_MAX / sizeof(real)) return 1;
    bytes = values * sizeof(real);
    if (bytes > 0 && (size_t)target > DUAL_NODE_TRIPLE_TASK_MEMORY / bytes)
        target = (INTEGER)MAX(
            (size_t)1, DUAL_NODE_TRIPLE_TASK_MEMORY / bytes);
    return target;
}

static int dual_node_build_task_frontier(
        struct cmdline_data *cmd, fcfc_balltreeptr tree1,
        fcfc_balltreeptr tree2, bool auto_correlation, INTEGER target,
        dual_node_triple_task **result, INTEGER *result_count)
{
    dual_node_triple_task *tasks;
    INTEGER count = 1;

    *result = NULL;
    *result_count = 0;
    if (target <= 0 || (size_t)target > SIZE_MAX / sizeof(*tasks)) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: invalid task-frontier size", DUAL_NODE_METHOD_NAME);
        return FAILURE;
    }
    tasks = calloc((size_t)target, sizeof(*tasks));
    if (tasks == NULL) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: task-frontier allocation failed", DUAL_NODE_METHOD_NAME);
        return FAILURE;
    }
    tasks[0] = auto_correlation
        ? dual_node_make_task(DUAL_NODE_TASK_AUTO3, tree1, 0,
                             NULL, 0, NULL, 0, 0U)
        : dual_node_make_task(DUAL_NODE_TASK_CROSS12, tree1, 0,
                             tree2, 0, NULL, 0, 0U);

    while (count < target) {
        INTEGER best = -1;
        int best_children = 0;
        long double largest_work = -1.0L;

        for (INTEGER i = 0; i < count; i++) {
            dual_node_triple_task children[4];
            const int nchildren = dual_node_split_task(&tasks[i], children);

            if (nchildren <= 0 || count + nchildren - 1 > target) continue;
            if (tasks[i].estimated_work > largest_work) {
                best = i;
                best_children = nchildren;
                largest_work = tasks[i].estimated_work;
            }
        }
        if (best < 0) break;
        {
            dual_node_triple_task children[4];
            const int nchildren = dual_node_split_task(&tasks[best], children);

            if (nchildren != best_children) {
                free(tasks);
                snprintf(cmd->error_message, _ERRORMSGSIZE_,
                         "%s: inconsistent task split", DUAL_NODE_METHOD_NAME);
                return FAILURE;
            }
            tasks[best] = children[0];
            for (int child = 1; child < nchildren; child++)
                tasks[count++] = children[child];
        }
    }
    *result = tasks;
    *result_count = count;
    return SUCCESS;
}

static void dual_node_run_task(const dual_node_search_context *context,
                              const dual_node_triple_task *task,
                              dual_node_triple_histogram *hist)
{
    switch (task->kind) {
        case DUAL_NODE_TASK_AUTO3:
            dual_node_process3_auto(
                context, task->trees[0], task->nodes[0], hist);
            break;
        case DUAL_NODE_TASK_TWO_ONE:
            dual_node_process21_auto(
                context, task->trees[0], task->nodes[0],
                task->nodes[1], hist);
            break;
        case DUAL_NODE_TASK_THREE:
            dual_node_process111(
                context, task->trees, task->nodes,
                task->orientation_mask, hist);
            break;
        case DUAL_NODE_TASK_CROSS12:
            dual_node_process12_cross(
                context, task->trees[0], task->nodes[0],
                task->trees[1], task->nodes[1], hist);
            break;
    }
}

#ifdef DUAL_NODE_LOG_MULTIPOLE_ENGINE
static int dual_node_search_direct_triples(
        struct cmdline_data *, struct global_data *,
        bodyptr *, INTEGER *, INTEGER, INTEGER *, int, int);
#ifdef THREEPCFCONVERGENCE

#include "dual_node_multipole.h"

#include "dual_node_search_policy.h"
#endif /* THREEPCFCONVERGENCE */

global int searchcalc_balltree_2balls_omp(
        struct cmdline_data *cmd, struct global_data *gd,
        bodyptr *btab, INTEGER *nbody, INTEGER ipmin, INTEGER *ipmax,
        int cat1, int cat2)
{
#if defined(DUAL_NODE_DISTRIBUTED_ENGINE) \
    && !defined(DUAL_NODE_DISTRIBUTED_DIRECT_TRIPLES)
    if (scanopt(cmd->options, "dual-node-direct-triples")) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s does not distribute dual-node-direct-triples; "
                 "use %s without that validation option",
                 DUAL_NODE_METHOD_NAME, DUAL_NODE_METHOD_NAME);
        return FAILURE;
    }
#else
    if (scanopt(cmd->options, "dual-node-direct-triples"))
        return dual_node_search_direct_triples(
            cmd, gd, btab, nbody, ipmin, ipmax, cat1, cat2);
#endif
#ifdef THREEPCFCONVERGENCE
    return dual_node_search_log_multipole(
        cmd, gd, btab, nbody, ipmin, ipmax, cat1, cat2);
#else
    return dual_node_search_direct_triples(
        cmd, gd, btab, nbody, ipmin, ipmax, cat1, cat2);
#endif
}
#endif /* DUAL_NODE_LOG_MULTIPOLE_ENGINE */
#endif /* DUAL_NODE_TASK_FRONTIER_ENGINE */
