/* Binary traversal view of the native cTreeBalls octree. */

#include <limits.h>
#include <stdint.h>
#include <stdlib.h>

#include "globaldefs.h"
#include "octree_2balls_tree.h"

#define OCTREE_2BALLS_ROOT 0

typedef struct {
    nodeptr node;
    INTEGER count;
} octree_2balls_child;

static real octree_2balls_field(bodyptr p)
{
#ifdef KappaAvgON
    return KappaAvg(p);
#else
    return Kappa(p);
#endif
}

static real octree_2balls_distance(const cballs_storage_real *a,
                                   const cballs_storage_real *b)
{
    real distance2 = 0.0;
    int k;

    DO_COORD(k) {
        const real difference = (real)a[k] - (real)b[k];
        distance2 += difference * difference;
    }
    return rsqrt(distance2);
}

static INTEGER octree_2balls_native_count(struct cmdline_data *cmd,
                                          nodeptr source)
{
    if (source == NULL) return 0;
    if (cballs_opt_read_mask(cmd) && Mask(source) == MASK_NODE_MASKED)
        return 0;
    if (Type(source) == BODY || Type(source) == BODY3)
        return !cballs_opt_read_mask(cmd) || Mask(source) == MASK_NODE_VALID;
    if (Type(source) == CELL && Nb(source) > 0) return Nb(source);
    return 0;
}

static int octree_2balls_reserve_node(struct cmdline_data *cmd,
                                      fcfc_balltreeptr tree,
                                      INTEGER *result)
{
    if (tree->nnode >= tree->capacity) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "octree-2balls-omp: binary octree capacity exceeded");
        return FAILURE;
    }
    *result = tree->nnode++;
    return SUCCESS;
}

static void octree_2balls_finish_parent(struct cmdline_data *cmd,
                                        fcfc_balltreeptr tree,
                                        INTEGER parent_index,
                                        INTEGER left_index,
                                        INTEGER right_index)
{
    fcfc_ballnode *parent = &tree->nodes[parent_index];
    const fcfc_ballnode *left = &tree->nodes[left_index];
    const fcfc_ballnode *right = &tree->nodes[right_index];
    const INTEGER left_count = left->last - left->first + 1;
    const INTEGER right_count = right->last - right->first + 1;
    const real count = (real)left_count + (real)right_count;
    const real mass = left->weight + right->weight;
    real radius = 0.0;
    real aggregate_radius = 0.0;
    int k;

    parent->first = left->first;
    parent->last = right->last;
    parent->left = left_index;
    parent->right = right_index;
    DO_COORD(k) {
        parent->center[k] = (cballs_storage_real)
            (((real)left->center[k] * (real)left_count
              + (real)right->center[k] * (real)right_count) / count);
        parent->cmpos[k] = mass > 0.0
            ? (cballs_storage_real)
                (((real)left->cmpos[k] * left->weight
                  + (real)right->cmpos[k] * right->weight) / mass)
            : parent->center[k];
    }

    radius = MAX(octree_2balls_distance(parent->center, left->center)
                     + (real)left->radius,
                 octree_2balls_distance(parent->center, right->center)
                     + (real)right->radius);
    aggregate_radius = MAX(
        octree_2balls_distance(parent->cmpos, left->center)
            + (real)left->radius,
        octree_2balls_distance(parent->cmpos, right->center)
            + (real)right->radius);
    parent->radius = cballs_store_search_bound(radius);
    parent->aggregate_radius = cmd->theta > 0.0
        ? cballs_store_search_bound(aggregate_radius / cmd->theta)
        : cballs_store_upper_bound(MAX_REAL_NUMBER);
    parent->weight = mass;
    parent->kappa_sum = left->kappa_sum + right->kappa_sum;
    parent->kappa_sq_sum = left->kappa_sq_sum + right->kappa_sq_sum;
    parent->field_weight_sum =
        left->field_weight_sum + right->field_weight_sum;
    parent->field_weight_sq_sum =
        left->field_weight_sq_sum + right->field_weight_sq_sum;
    parent->weighted_kappa_sum =
        left->weighted_kappa_sum + right->weighted_kappa_sum;
    parent->weighted_kappa_sq_sum =
        left->weighted_kappa_sq_sum + right->weighted_kappa_sq_sum;
    parent->kappa = parent->kappa_sum / count;
}

static int octree_2balls_build_native(struct cmdline_data *,
                                      fcfc_balltreeptr, nodeptr, int,
                                      INTEGER *);

static int octree_2balls_build_group(struct cmdline_data *cmd,
                                     fcfc_balltreeptr tree,
                                     const octree_2balls_child *children,
                                     int first, int last, int depth,
                                     INTEGER *result)
{
    INTEGER parent_index;
    INTEGER left_index;
    INTEGER right_index;
    INTEGER total = 0;
    INTEGER prefix = 0;
    INTEGER best_difference;
    int split = first;
    int i;

    if (first == last)
        return octree_2balls_build_native(
            cmd, tree, children[first].node, depth, result);
    if (octree_2balls_reserve_node(cmd, tree, &parent_index) == FAILURE)
        return FAILURE;
    if (depth > tree->max_depth) tree->max_depth = depth;

    for (i = first; i <= last; i++) total += children[i].count;
    best_difference = total;
    for (i = first; i < last; i++) {
        INTEGER difference;
        prefix += children[i].count;
        difference = prefix > total - prefix
            ? prefix - (total - prefix) : (total - prefix) - prefix;
        if (difference < best_difference) {
            best_difference = difference;
            split = i;
        }
    }

    if (octree_2balls_build_group(cmd, tree, children, first, split,
                                  depth + 1, &left_index) == FAILURE
        || octree_2balls_build_group(cmd, tree, children, split + 1, last,
                                     depth + 1, &right_index) == FAILURE)
        return FAILURE;
    octree_2balls_finish_parent(
        cmd, tree, parent_index, left_index, right_index);
    *result = parent_index;
    return SUCCESS;
}

static int octree_2balls_build_leaf(struct cmdline_data *cmd,
                                    fcfc_balltreeptr tree, bodyptr body,
                                    int depth, INTEGER *result)
{
    fcfc_ballnode *target;
    INTEGER index;
    const real field = octree_2balls_field(body);
    const real field_weight = cballs_opt_weights_norm(cmd) ? Weight(body) : 0.0;
    int k;

    if (cballs_opt_read_mask(cmd) && Mask(body) != MASK_NODE_VALID) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "octree-2balls: masked body reached the binary search tree");
        return FAILURE;
    }
    if (tree->npoint >= tree->capacity / 2) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "octree-2balls-omp: point capacity exceeded");
        return FAILURE;
    }
    if (octree_2balls_reserve_node(cmd, tree, &index) == FAILURE)
        return FAILURE;
    target = &tree->nodes[index];
    target->first = tree->npoint;
    target->last = tree->npoint;
    target->left = -1;
    target->right = -1;
    tree->bptr[tree->npoint++] = body;
    if (depth > tree->max_depth) tree->max_depth = depth;
    DO_COORD(k) {
        target->center[k] = Pos(body)[k];
        target->cmpos[k] = Pos(body)[k];
    }
    target->radius = 0.0;
    target->aggregate_radius = 0.0;
    target->weight = Mass(body);
    target->kappa = field;
    target->kappa_sum = field;
    target->kappa_sq_sum = field * field;
    target->field_weight_sum = field_weight;
    target->field_weight_sq_sum = field_weight * field_weight;
    target->weighted_kappa_sum = field_weight * field;
    target->weighted_kappa_sq_sum =
        (field_weight * field) * (field_weight * field);
    *result = index;
    return SUCCESS;
}

static int octree_2balls_build_native(struct cmdline_data *cmd,
                                      fcfc_balltreeptr tree, nodeptr source,
                                      int depth, INTEGER *result)
{
    octree_2balls_child children[NSUB];
    int child_count = 0;
    int i;

    if (Type(source) == BODY || Type(source) == BODY3)
        return octree_2balls_build_leaf(
            cmd, tree, (bodyptr)source, depth, result);
    if (Type(source) != CELL) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "octree-2balls-omp: unsupported native node type %d",
                 (int)Type(source));
        return FAILURE;
    }

    for (i = 0; i < NSUB; i++) {
        const INTEGER count = octree_2balls_native_count(cmd, Subp(source)[i]);
        if (count <= 0) continue;
        children[child_count].node = Subp(source)[i];
        children[child_count].count = count;
        child_count++;
    }
    if (child_count == 0) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "octree-2balls-omp: native octree contains an empty live cell");
        return FAILURE;
    }
    return octree_2balls_build_group(
        cmd, tree, children, 0, child_count - 1, depth, result);
}

int octree_2balls_tree_build(struct cmdline_data *cmd,
                             struct global_data *gd, bodyptr btab,
                             INTEGER nbody, int leaf_capacity,
                             fcfc_balltreeptr *result)
{
    fcfc_balltreeptr tree = NULL;
    INTEGER root = -1;
    INTEGER live_count;
    int catalog = -1;
    int i;

    (void)leaf_capacity;
    if (result == NULL || btab == NULL || nbody <= 0) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "octree-2balls-omp: invalid tree dimensions");
        return FAILURE;
    }
    *result = NULL;
    for (i = 0; i < gd->ninfiles; i++) {
        if (bodytable[i] == btab) {
            catalog = i;
            break;
        }
    }
    if (catalog < 0 || roottable[catalog] == NULL) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "octree-2balls-omp: native octree is not available for this catalog");
        return FAILURE;
    }
    /* The traversal and all its moments contain valid bodies only. Native
     * mixed cells are opened, never copied as aggregate search nodes. */
    live_count = octree_2balls_native_count(cmd, (nodeptr)roottable[catalog]);
    if (live_count <= 0 || live_count > nbody) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "octree-2balls: catalog %d has no unmasked bodies "
                 "or an invalid native body count", catalog + 1);
        return FAILURE;
    }
    if ((uintmax_t)live_count >
#ifdef LONGINT
        (uintmax_t)LONG_MAX / 2
#else
        (uintmax_t)INT_MAX / 2
#endif
        || (uintmax_t)live_count > (uintmax_t)SIZE_MAX / (2 * sizeof(fcfc_ballnode))
        || (uintmax_t)live_count > (uintmax_t)SIZE_MAX / sizeof(bodyptr)
#ifdef SINGLEP
        || (uintmax_t)live_count > (uintmax_t)SIZE_MAX / sizeof(fcfc_ballpoint)
#endif
        ) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "octree-2balls-omp: tree allocation size overflow");
        return FAILURE;
    }

    tree = calloc(1, sizeof(*tree));
    if (tree == NULL) goto allocation_failure;
    tree->capacity = 2 * live_count;
    tree->bptr = malloc((size_t)live_count * sizeof(*tree->bptr));
    tree->nodes = calloc((size_t)tree->capacity, sizeof(*tree->nodes));
    if (tree->bptr == NULL || tree->nodes == NULL) goto allocation_failure;
    if (octree_2balls_build_native(cmd, tree, (nodeptr)roottable[catalog],
                                   0, &root) == FAILURE)
        goto failure;
    if (root != OCTREE_2BALLS_ROOT || tree->npoint != live_count) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "octree-2balls-omp: invalid binary octree root");
        goto failure;
    }

#ifdef SINGLEP
    tree->packed_points = malloc(
        (size_t)tree->npoint * sizeof(*tree->packed_points));
    if (tree->packed_points == NULL) goto allocation_failure;
    for (INTEGER point = 0; point < tree->npoint; point++) {
        SETV(tree->packed_points[point].pos, Pos(tree->bptr[point]));
        tree->packed_points[point].kappa = Kappa(tree->bptr[point]);
    }
#endif

    gd->bytes_tot += sizeof(*tree)
        + (size_t)live_count * sizeof(*tree->bptr)
        + (size_t)tree->capacity * sizeof(*tree->nodes)
#ifdef SINGLEP
        + (size_t)tree->npoint * sizeof(*tree->packed_points)
#endif
        ;
    if (cballs_opt_read_mask(cmd))
        verb_print(cmd->verbose,
                   "octree-2balls: catalog %d mask keeps %" INTEGER_FMT
                   " of %" INTEGER_FMT " bodies; binary capacity=%" INTEGER_FMT
                   " nodes\n", catalog + 1, live_count, nbody, tree->capacity);
    *result = tree;
    return SUCCESS;

allocation_failure:
    snprintf(cmd->error_message, _ERRORMSGSIZE_,
             "octree-2balls-omp: memory allocation failed");
failure:
    octree_2balls_tree_free(tree);
    return FAILURE;
}

int octree_2balls_tree_frontier(struct cmdline_data *cmd,
                                fcfc_balltreeptr tree,
                                INTEGER minimum_count, INTEGER **result,
                                INTEGER *result_count)
{
    INTEGER target;
    INTEGER count = 1;
    INTEGER *frontier;

    if (tree == NULL || tree->nodes == NULL || tree->nnode < 1
        || minimum_count < 1 || result == NULL || result_count == NULL) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "octree-2balls-omp: invalid task-frontier arguments");
        return FAILURE;
    }
    *result = NULL;
    *result_count = 0;
    target = MIN(minimum_count, tree->npoint);
    if ((size_t)target > SIZE_MAX / sizeof(*frontier)) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "octree-2balls-omp: task-frontier size overflow");
        return FAILURE;
    }
    frontier = malloc((size_t)target * sizeof(*frontier));
    if (frontier == NULL) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "octree-2balls-omp: task-frontier allocation failed");
        return FAILURE;
    }
    frontier[0] = OCTREE_2BALLS_ROOT;

    while (count < target) {
        INTEGER selected = -1;
        INTEGER selected_points = -1;
        real selected_radius = -1.0;
        INTEGER i;

        for (i = 0; i < count; i++) {
            const fcfc_ballnode *node = &tree->nodes[frontier[i]];
            const INTEGER points = node->last - node->first + 1;
            if (node->left < 0) continue;
            if (points > selected_points
                || (points == selected_points
                    && (real)node->radius > selected_radius)) {
                selected = i;
                selected_points = points;
                selected_radius = (real)node->radius;
            }
        }
        if (selected < 0) break;
        {
            const fcfc_ballnode *node = &tree->nodes[frontier[selected]];
            frontier[selected] = node->left;
            frontier[count++] = node->right;
        }
    }

    *result = frontier;
    *result_count = count;
    return SUCCESS;
}

void octree_2balls_tree_free(fcfc_balltreeptr tree)
{
    if (tree == NULL) return;
#ifdef SINGLEP
    free(tree->packed_points);
#endif
    free(tree->nodes);
    free(tree->bptr);
    free(tree);
}
