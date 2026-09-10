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

static real octree_2balls_distance_squared(const cballs_storage_real *a,
                                           const cballs_storage_real *b)
{
    real distance2 = 0.0;
    int k;

    DO_COORD(k) {
        const real difference = (real)a[k] - (real)b[k];
        distance2 += difference * difference;
    }
    return distance2;
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
    const real mass = left->weight + right->weight;
    real radius_squared = 0.0;
    INTEGER point;
    int k;

    parent->first = left->first;
    parent->last = right->last;
    parent->left = left_index;
    parent->right = right_index;
    DO_COORD(k) {
        parent->cmpos[k] = mass > 0.0
            ? (cballs_storage_real)
                (((real)left->cmpos[k] * left->weight
                  + (real)right->cmpos[k] * right->weight) / mass)
            : (cballs_storage_real)
                (((real)left->cmpos[k] * (real)left_count
                  + (real)right->cmpos[k] * (real)right_count)
                 / ((real)left_count + (real)right_count));
        parent->center[k] = parent->cmpos[k];
    }

    /* Match dual-node's cell geometry: aggregate at the centroid and measure
     * the exact maximum point displacement from the stored centroid. */
    for (point = parent->first; point <= parent->last; point++)
        radius_squared = MAX(radius_squared,
            octree_2balls_distance_squared(
                parent->center, Pos(tree->bptr[point])));
    parent->radius = cballs_store_search_bound(rsqrt(radius_squared));
    parent->aggregate_radius = cmd->theta > 0.0
        ? cballs_store_search_bound(rsqrt(radius_squared) / cmd->theta)
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
    parent->kappa = parent->kappa_sum
        / ((real)left_count + (real)right_count);
}

static int octree_2balls_build_native(struct cmdline_data *,
                                      fcfc_balltreeptr, nodeptr, int, int,
                                      INTEGER *);

static void octree_2balls_sort_children(octree_2balls_child *children,
                                        int child_count)
{
    real minimum[NDIM];
    real maximum[NDIM];
    real largest_span = -1.0;
    int axis = 0;
    int i;
    int k;

    DO_COORD(k) {
        minimum[k] = maximum[k] = (real)Pos(children[0].node)[k];
        for (i = 1; i < child_count; i++) {
            minimum[k] = MIN(minimum[k], (real)Pos(children[i].node)[k]);
            maximum[k] = MAX(maximum[k], (real)Pos(children[i].node)[k]);
        }
        if (maximum[k] - minimum[k] > largest_span) {
            largest_span = maximum[k] - minimum[k];
            axis = k;
        }
    }

    for (i = 1; i < child_count; i++) {
        const octree_2balls_child value = children[i];
        int j = i;

        while (j > 0
               && (real)Pos(children[j - 1].node)[axis]
                    > (real)Pos(value.node)[axis]) {
            children[j] = children[j - 1];
            j--;
        }
        children[j] = value;
    }
}

static int octree_2balls_build_group(struct cmdline_data *cmd,
                                     fcfc_balltreeptr tree,
                                     const octree_2balls_child *children,
                                     int first, int last, int depth,
                                     int leaf_capacity,
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
            cmd, tree, children[first].node, depth, leaf_capacity, result);
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
                                  depth + 1, leaf_capacity,
                                  &left_index) == FAILURE
        || octree_2balls_build_group(cmd, tree, children, split + 1, last,
                                     depth + 1, leaf_capacity,
                                     &right_index) == FAILURE)
        return FAILURE;
    octree_2balls_finish_parent(
        cmd, tree, parent_index, left_index, right_index);
    *result = parent_index;
    return SUCCESS;
}

static int octree_2balls_collect_bodies(struct cmdline_data *cmd,
                                        fcfc_balltreeptr tree,
                                        nodeptr source)
{
    int i;

    if (source == NULL) return SUCCESS;
    if (cballs_opt_read_mask(cmd) && Mask(source) == MASK_NODE_MASKED)
        return SUCCESS;
    if (Type(source) == BODY || Type(source) == BODY3) {
        if (cballs_opt_read_mask(cmd) && Mask(source) != MASK_NODE_VALID) {
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     "octree-2balls: non-valid body reached the binary tree");
            return FAILURE;
        }
        if (tree->npoint >= tree->capacity / 2) {
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     "octree-2balls-omp: point capacity exceeded");
            return FAILURE;
        }
        tree->bptr[tree->npoint++] = (bodyptr)source;
        return SUCCESS;
    }
    if (Type(source) != CELL) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "octree-2balls-omp: unsupported native node type %d",
                 (int)Type(source));
        return FAILURE;
    }
    for (i = 0; i < NSUB; i++)
        if (octree_2balls_collect_bodies(
                cmd, tree, Subp(source)[i]) == FAILURE)
            return FAILURE;
    return SUCCESS;
}

static int octree_2balls_build_leaf(struct cmdline_data *cmd,
                                    fcfc_balltreeptr tree, nodeptr source,
                                    INTEGER expected_count, int depth,
                                    INTEGER *result)
{
    compute_vector cmpos_sum;
    compute_vector geometric_center;
    fcfc_ballnode *target;
    INTEGER index;
    INTEGER i;
    real farthest2;
    int k;

    if (octree_2balls_reserve_node(cmd, tree, &index) == FAILURE)
        return FAILURE;
    target = &tree->nodes[index];
    target->first = tree->npoint;
    if (octree_2balls_collect_bodies(cmd, tree, source) == FAILURE)
        return FAILURE;
    target->last = tree->npoint - 1;
    if (target->last - target->first + 1 != expected_count) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "octree-2balls-omp: native/body count mismatch");
        return FAILURE;
    }
    target->left = -1;
    target->right = -1;
    if (depth > tree->max_depth) tree->max_depth = depth;

    CLRV(cmpos_sum);
    CLRV(geometric_center);
    target->weight = 0.0;
    target->kappa_sum = 0.0;
    target->kappa_sq_sum = 0.0;
    target->field_weight_sum = 0.0;
    target->field_weight_sq_sum = 0.0;
    target->weighted_kappa_sum = 0.0;
    target->weighted_kappa_sq_sum = 0.0;
    for (i = target->first; i <= target->last; i++) {
        bodyptr body = tree->bptr[i];
        const real mass = Mass(body);
        const real field = octree_2balls_field(body);
        const real field_weight = Weight(body);

        DO_COORD(k) {
            cmpos_sum[k] += mass * (real)Pos(body)[k];
            geometric_center[k] += (real)Pos(body)[k];
        }
        target->weight += mass;
        target->kappa_sum += field;
        target->kappa_sq_sum += field * field;
        target->field_weight_sum += field_weight;
        target->field_weight_sq_sum += field_weight * field_weight;
        target->weighted_kappa_sum += field_weight * field;
        target->weighted_kappa_sq_sum +=
            (field_weight * field) * (field_weight * field);
    }
    DO_COORD(k)
        target->cmpos[k] = (cballs_storage_real)
            (target->weight > 0.0
             ? cmpos_sum[k] / target->weight
             : geometric_center[k] / (real)expected_count);
    target->kappa = target->kappa_sum / (real)expected_count;

    SETV(target->center, target->cmpos);
    farthest2 = 0.0;
    for (i = target->first; i <= target->last; i++)
        farthest2 = MAX(farthest2, octree_2balls_distance_squared(
            target->center, Pos(tree->bptr[i])));
    target->radius = cballs_store_search_bound(rsqrt(farthest2));
    target->aggregate_radius = cmd->theta > 0.0
        ? cballs_store_search_bound(rsqrt(farthest2) / cmd->theta)
        : cballs_store_upper_bound(MAX_REAL_NUMBER);
    *result = index;
    return SUCCESS;
}

static int octree_2balls_build_native(struct cmdline_data *cmd,
                                      fcfc_balltreeptr tree, nodeptr source,
                                      int depth, int leaf_capacity,
                                      INTEGER *result)
{
    octree_2balls_child children[NSUB];
    const INTEGER source_count = octree_2balls_native_count(cmd, source);
    int child_count = 0;
    int i;

    if (source_count <= 0) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "octree-2balls-omp: empty node reached the binary tree");
        return FAILURE;
    }
    if (Type(source) == BODY || Type(source) == BODY3
        || source_count <= leaf_capacity)
        return octree_2balls_build_leaf(
            cmd, tree, source, source_count, depth, result);
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
    octree_2balls_sort_children(children, child_count);
    return octree_2balls_build_group(
        cmd, tree, children, 0, child_count - 1, depth,
        leaf_capacity, result);
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

    if (result == NULL || btab == NULL || nbody <= 0) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "octree-2balls-omp: invalid tree dimensions");
        return FAILURE;
    }
    if (leaf_capacity < 1) leaf_capacity = 1;
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
        || (uintmax_t)live_count > (uintmax_t)SIZE_MAX / sizeof(fcfc_ballpoint)
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
                                   0, leaf_capacity, &root) == FAILURE)
        goto failure;
    if (root != OCTREE_2BALLS_ROOT || tree->npoint != live_count) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "octree-2balls-omp: invalid binary octree root");
        goto failure;
    }

    tree->packed_points = malloc(
        (size_t)tree->npoint * sizeof(*tree->packed_points));
    if (tree->packed_points == NULL) goto allocation_failure;
    for (INTEGER point = 0; point < tree->npoint; point++) {
        const real field = octree_2balls_field(tree->bptr[point]);
        const real weight = Weight(tree->bptr[point]);
        SETV(tree->packed_points[point].pos, Pos(tree->bptr[point]));
        tree->packed_points[point].kappa = field;
        tree->packed_points[point].weight = weight;
        tree->packed_points[point].weighted_kappa = weight * field;
        tree->packed_points[point].source = tree->bptr[point];
    }

    gd->bytes_tot += sizeof(*tree)
        + (size_t)live_count * sizeof(*tree->bptr)
        + (size_t)tree->capacity * sizeof(*tree->nodes)
        + (size_t)tree->npoint * sizeof(*tree->packed_points)
        ;
    if (cballs_opt_read_mask(cmd))
        verb_print(cmd->verbose,
                   "octree-2balls: catalog %d mask keeps %" INTEGER_FMT
                   " of %" INTEGER_FMT " bodies; binary capacity=%" INTEGER_FMT
                   " nodes\n", catalog + 1, live_count, nbody, tree->capacity);
    verb_print(cmd->verbose >= 2,
               "octree-2balls: compact binary view has %" INTEGER_FMT
               " nodes for %" INTEGER_FMT " bodies (leaf capacity %d)\n",
               tree->nnode, tree->npoint, leaf_capacity);
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
    free(tree->packed_points);
    free(tree->nodes);
    free(tree->bptr);
    free(tree);
}
