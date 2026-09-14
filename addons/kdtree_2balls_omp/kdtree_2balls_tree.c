/* Median-split KD tree used by the dual-node-style two-ball traversal. */

#include <limits.h>
#include <stdint.h>
#include <stdlib.h>

#ifdef OPENMPCODE
#include <omp.h>
#endif

#include "globaldefs.h"
#include "kdtree_2balls_tree.h"

#define KDTREE_2BALLS_ROOT 0
#define KDTREE_2BALLS_PARALLEL_BUILD_MIN ((INTEGER)262144)
#define KDTREE_2BALLS_PARALLEL_SUBTREE_MIN ((INTEGER)65536)
#define KDTREE_2BALLS_SHAPE_CACHE_SIZE 256

typedef enum {
    KDTREE_2BALLS_PIVOT,
    KDTREE_2BALLS_NEIGHBOR
} kdtree_2balls_role;

static real kdtree_2balls_raw_field(bodyptr body)
{
#ifdef KappaAvgON
    return KappaAvg(body);
#else
    return Kappa(body);
#endif
}

static bool kdtree_2balls_valid_body(
        const struct cmdline_data *cmd, bodyptr body,
        kdtree_2balls_role role)
{
    if (cballs_opt_read_mask(cmd) && Mask(body) == MASK_NODE_MASKED)
        return FALSE;
#ifdef SMOOTHPIVOT
    if (role == KDTREE_2BALLS_PIVOT
        && cballs_opt_smooth_pivot(cmd) && Update(body) == FALSE)
        return FALSE;
#else
    (void)role;
#endif
    return TRUE;
}

static void kdtree_2balls_point_values(
        const struct cmdline_data *cmd, bodyptr body,
        kdtree_2balls_role role, real *field, real *normalization)
{
#ifdef SMOOTHPIVOT
    if (role == KDTREE_2BALLS_PIVOT && cballs_opt_smooth_pivot(cmd)) {
        if (cballs_opt_weights_norm(cmd)) {
            *normalization = WeightRmin(body);
            *field = KappaRmin(body);
        } else {
            *normalization = (real)MAX(NbRmin(body), 1);
            *field = KappaRmin(body);
        }
        return;
    }
#else
    (void)role;
#endif
    *normalization = cballs_opt_weights_norm(cmd) ? Weight(body) : 1.0;
    *field = *normalization * kdtree_2balls_raw_field(body);
}

static int kdtree_2balls_before(bodyptr left, bodyptr right, int axis)
{
    const real a = (real)Pos(left)[axis];
    const real b = (real)Pos(right)[axis];

    if (a < b) return TRUE;
    if (a > b) return FALSE;
    return (uintptr_t)left < (uintptr_t)right;
}

static void kdtree_2balls_swap(bodyptr *left, bodyptr *right)
{
    bodyptr temporary = *left;
    *left = *right;
    *right = temporary;
}

static void kdtree_2balls_select(
        bodyptr *points, INTEGER lo, INTEGER hi, INTEGER target, int axis)
{
    while (lo < hi) {
        INTEGER i = lo;
        INTEGER j = hi;
        bodyptr pivot = points[lo + (hi - lo) / 2];

        while (i <= j) {
            while (i <= hi && kdtree_2balls_before(points[i], pivot, axis))
                i++;
            while (j >= lo && kdtree_2balls_before(pivot, points[j], axis))
                j--;
            if (i <= j) {
                kdtree_2balls_swap(&points[i], &points[j]);
                i++;
                j--;
            }
        }
        if (target <= j)
            hi = j;
        else if (target >= i)
            lo = i;
        else
            return;
    }
}

static real kdtree_2balls_distance_squared(
        const cballs_storage_real *left, const cballs_storage_real *right)
{
    real distance2 = 0.0;
    int axis;

    DO_COORD(axis)
        distance2 += rsqr((real)left[axis] - (real)right[axis]);
    return distance2;
}

static void kdtree_2balls_aggregate(
        const struct cmdline_data *cmd, fcfc_balltreeptr tree,
        fcfc_ballnode *node, kdtree_2balls_role role)
{
    compute_vector weighted_position;
    compute_vector geometric_position;
    real position_weight = 0.0;
    real radius2 = 0.0;
    INTEGER point;
    int axis;

    CLRV(weighted_position);
    CLRV(geometric_position);
    node->weight = 0.0;
    node->kappa_sum = 0.0;
    node->kappa_sq_sum = 0.0;
    node->field_weight_sum = 0.0;
    node->field_weight_sq_sum = 0.0;
    node->weighted_kappa_sum = 0.0;
    node->weighted_kappa_sq_sum = 0.0;

    for (point = node->first; point <= node->last; point++) {
        bodyptr body = tree->bptr[point];
        real field;
        real normalization;
        const real mass = Mass(body);
        const real raw_field = kdtree_2balls_raw_field(body);

        kdtree_2balls_point_values(cmd, body, role, &field, &normalization);
        DO_COORD(axis) {
            weighted_position[axis] += mass * (real)Pos(body)[axis];
            geometric_position[axis] += (real)Pos(body)[axis];
        }
        position_weight += mass;
        node->weight += mass;
        node->kappa_sum += raw_field;
        node->kappa_sq_sum += raw_field * raw_field;
        node->field_weight_sum += normalization;
        node->field_weight_sq_sum += normalization * normalization;
        node->weighted_kappa_sum += field;
        node->weighted_kappa_sq_sum += field * field;
    }

    DO_COORD(axis)
        node->cmpos[axis] = (cballs_storage_real)
            (position_weight > 0.0
             ? weighted_position[axis] / position_weight
             : geometric_position[axis]
               / (real)(node->last - node->first + 1));
    SETV(node->center, node->cmpos);
    node->kappa = node->kappa_sum
        / (real)(node->last - node->first + 1);

    for (point = node->first; point <= node->last; point++)
        radius2 = MAX(radius2, kdtree_2balls_distance_squared(
            node->center, Pos(tree->bptr[point])));
    node->radius = cballs_store_search_bound(rsqrt(radius2));
    node->aggregate_radius = cmd->theta > 0.0
        ? cballs_store_search_bound(rsqrt(radius2) / cmd->theta)
        : cballs_store_upper_bound(MAX_REAL_NUMBER);
}

static int kdtree_2balls_build_node_serial(
        struct cmdline_data *cmd, fcfc_balltreeptr tree,
        INTEGER first, INTEGER last, int leaf_capacity, int depth,
        kdtree_2balls_role role, INTEGER *result)
{
    cballs_storage_real minimum[NDIM];
    cballs_storage_real maximum[NDIM];
    fcfc_ballnode *node;
    INTEGER index;
    INTEGER point;
    int split_axis = 0;
    int axis;

    if (tree->nnode >= tree->capacity) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "kdtree-2balls: node capacity exceeded");
        return FAILURE;
    }
    index = tree->nnode++;
    node = &tree->nodes[index];
    node->first = first;
    node->last = last;
    node->left = -1;
    node->right = -1;
    if (depth > tree->max_depth) tree->max_depth = depth;

    DO_COORD(axis)
        minimum[axis] = maximum[axis] = Pos(tree->bptr[first])[axis];
    for (point = first + 1; point <= last; point++)
        DO_COORD(axis) {
            minimum[axis] = MIN(minimum[axis], Pos(tree->bptr[point])[axis]);
            maximum[axis] = MAX(maximum[axis], Pos(tree->bptr[point])[axis]);
        }
    DO_COORD(axis)
        if ((real)maximum[axis] - (real)minimum[axis]
            > (real)maximum[split_axis] - (real)minimum[split_axis])
            split_axis = axis;

    kdtree_2balls_aggregate(cmd, tree, node, role);
    if (last - first + 1 > leaf_capacity) {
        const INTEGER middle = first + (last - first + 1) / 2;

        kdtree_2balls_select(tree->bptr, first, last, middle, split_axis);
        if (kdtree_2balls_build_node_serial(
                cmd, tree, first, middle - 1, leaf_capacity, depth + 1,
                role, &node->left) == FAILURE
            || kdtree_2balls_build_node_serial(
                cmd, tree, middle, last, leaf_capacity, depth + 1,
                role, &node->right) == FAILURE)
            return FAILURE;
    }
    *result = index;
    return SUCCESS;
}

typedef struct {
    INTEGER points;
    INTEGER nodes;
    int depth;
} kdtree_2balls_shape;

typedef struct {
    kdtree_2balls_shape entries[KDTREE_2BALLS_SHAPE_CACHE_SIZE];
    INTEGER *nodes_by_points;
    int count;
    int leaf_capacity;
} kdtree_2balls_shape_cache;

static int kdtree_2balls_shape_find(
        const kdtree_2balls_shape_cache *cache, INTEGER points)
{
    for (int i = 0; i < cache->count; i++)
        if (cache->entries[i].points == points) return i;
    return -1;
}

static int kdtree_2balls_shape_add(
        kdtree_2balls_shape_cache *cache, INTEGER points)
{
    int found = kdtree_2balls_shape_find(cache, points);
    INTEGER nodes = 1;
    int depth = 0;

    if (found >= 0) return found;
    if (points > cache->leaf_capacity) {
        const INTEGER left_points = points / 2;
        const INTEGER right_points = points - left_points;
        const int left = kdtree_2balls_shape_add(cache, left_points);
        const int right = kdtree_2balls_shape_add(cache, right_points);

        if (left < 0 || right < 0) return -1;
        nodes += cache->entries[left].nodes + cache->entries[right].nodes;
        depth = 1 + MAX(cache->entries[left].depth,
                        cache->entries[right].depth);
    }
    if (cache->count >= KDTREE_2BALLS_SHAPE_CACHE_SIZE) return -1;
    found = cache->count++;
    cache->entries[found].points = points;
    cache->entries[found].nodes = nodes;
    cache->entries[found].depth = depth;
    cache->nodes_by_points[points] = nodes;
    return found;
}

static INTEGER kdtree_2balls_shape_nodes(
        const kdtree_2balls_shape_cache *cache, INTEGER points)
{
    return cache->nodes_by_points[points];
}

static int kdtree_2balls_build_node(
        struct cmdline_data *cmd, fcfc_balltreeptr tree,
        INTEGER first, INTEGER last, int leaf_capacity, int depth,
        kdtree_2balls_role role, INTEGER index,
        const kdtree_2balls_shape_cache *shape, bool parallel_build)
{
    cballs_storage_real minimum[NDIM];
    cballs_storage_real maximum[NDIM];
    fcfc_ballnode *node;
    INTEGER point;
    int split_axis = 0;
    int axis;

    if (index < 0 || index >= tree->capacity) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "kdtree-2balls: node capacity exceeded");
        return FAILURE;
    }
    node = &tree->nodes[index];
    node->first = first;
    node->last = last;
    node->left = -1;
    node->right = -1;

    DO_COORD(axis)
        minimum[axis] = maximum[axis] = Pos(tree->bptr[first])[axis];
    for (point = first + 1; point <= last; point++)
        DO_COORD(axis) {
            minimum[axis] = MIN(minimum[axis], Pos(tree->bptr[point])[axis]);
            maximum[axis] = MAX(maximum[axis], Pos(tree->bptr[point])[axis]);
        }
    DO_COORD(axis)
        if ((real)maximum[axis] - (real)minimum[axis]
            > (real)maximum[split_axis] - (real)minimum[split_axis])
            split_axis = axis;

    kdtree_2balls_aggregate(cmd, tree, node, role);
    if (last - first + 1 > leaf_capacity) {
        const INTEGER points = last - first + 1;
        const INTEGER middle = first + points / 2;
        const INTEGER left_points = middle - first;
        const INTEGER left_index = index + 1;
        const INTEGER right_index = left_index
            + kdtree_2balls_shape_nodes(shape, left_points);
        int left_status = SUCCESS;
        int right_status = SUCCESS;

        kdtree_2balls_select(tree->bptr, first, last, middle, split_axis);
        node->left = left_index;
        node->right = right_index;
#ifdef OPENMPCODE
        if (parallel_build && points >= KDTREE_2BALLS_PARALLEL_SUBTREE_MIN) {
#pragma omp task shared(left_status)
            left_status = kdtree_2balls_build_node(
                cmd, tree, first, middle - 1, leaf_capacity, depth + 1,
                role, left_index, shape, parallel_build);
#pragma omp task shared(right_status)
            right_status = kdtree_2balls_build_node(
                cmd, tree, middle, last, leaf_capacity, depth + 1,
                role, right_index, shape, parallel_build);
#pragma omp taskwait
        } else
#endif
        {
            left_status = kdtree_2balls_build_node(
                cmd, tree, first, middle - 1, leaf_capacity, depth + 1,
                role, left_index, shape, parallel_build);
            right_status = kdtree_2balls_build_node(
                cmd, tree, middle, last, leaf_capacity, depth + 1,
                role, right_index, shape, parallel_build);
        }
        if (left_status == FAILURE || right_status == FAILURE)
            return FAILURE;
    }
    return SUCCESS;
}

void kdtree_2balls_tree_free(fcfc_balltreeptr tree)
{
    if (tree == NULL) return;
    free(tree->packed_points);
    free(tree->nodes);
    free(tree->bptr);
    free(tree);
}

static int kdtree_2balls_tree_build(
        struct cmdline_data *cmd, struct global_data *gd,
        bodyptr body_table, INTEGER body_count, int leaf_capacity,
        kdtree_2balls_role role, fcfc_balltreeptr *result)
{
    fcfc_balltreeptr tree = NULL;
    kdtree_2balls_shape_cache shape = {0};
    INTEGER valid_count = 0;
    INTEGER root = -1;
    INTEGER source;
    int root_shape = -1;
    int build_status = FAILURE;
    bool parallel_build = FALSE;

    if (result == NULL || body_table == NULL
        || body_count <= 0 || leaf_capacity <= 0) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "kdtree-2balls: invalid tree dimensions");
        return FAILURE;
    }
    *result = NULL;
    for (source = 0; source < body_count; source++)
        if (kdtree_2balls_valid_body(cmd, nthBody(body_table, source), role))
            valid_count++;
    if (valid_count <= 0) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "kdtree-2balls: mask/smoothing selected no bodies");
        return FAILURE;
    }
    if ((uintmax_t)valid_count >
#ifdef LONGINT
        (uintmax_t)LONG_MAX / 2 ||
#else
        (uintmax_t)INT_MAX / 2 ||
#endif
        (uintmax_t)valid_count > (uintmax_t)SIZE_MAX / (2 * sizeof(fcfc_ballnode))
        || (uintmax_t)valid_count > (uintmax_t)SIZE_MAX / sizeof(bodyptr)
        || (uintmax_t)valid_count > (uintmax_t)SIZE_MAX / sizeof(fcfc_ballpoint)) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "kdtree-2balls: tree size overflow");
        return FAILURE;
    }

#ifdef OPENMPCODE
    parallel_build = valid_count >= KDTREE_2BALLS_PARALLEL_BUILD_MIN
                  && omp_get_max_threads() > 1;
#endif
    if (parallel_build) {
        shape.leaf_capacity = leaf_capacity;
        shape.nodes_by_points = calloc(
            (size_t)valid_count + 1, sizeof(*shape.nodes_by_points));
        if (shape.nodes_by_points == NULL) {
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     "kdtree-2balls: tree-shape allocation failed");
            return FAILURE;
        }
        root_shape = kdtree_2balls_shape_add(&shape, valid_count);
        if (root_shape < 0 || shape.entries[root_shape].nodes <= 0) {
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     "kdtree-2balls: tree-shape construction failed");
            free(shape.nodes_by_points);
            return FAILURE;
        }
    }

    tree = calloc(1, sizeof(*tree));
    if (tree == NULL) goto allocation_failure;
    tree->npoint = valid_count;
    tree->nnode = parallel_build ? shape.entries[root_shape].nodes : 0;
    tree->capacity = parallel_build ? tree->nnode : 2 * valid_count;
    tree->max_depth = parallel_build ? shape.entries[root_shape].depth : 0;
    tree->bptr = malloc((size_t)valid_count * sizeof(*tree->bptr));
    tree->nodes = calloc((size_t)tree->capacity, sizeof(*tree->nodes));
    tree->packed_points = malloc(
        (size_t)valid_count * sizeof(*tree->packed_points));
    if (tree->bptr == NULL || tree->nodes == NULL
        || tree->packed_points == NULL)
        goto allocation_failure;

    valid_count = 0;
    for (source = 0; source < body_count; source++) {
        bodyptr body = nthBody(body_table, source);
        if (kdtree_2balls_valid_body(cmd, body, role))
            tree->bptr[valid_count++] = body;
    }
#ifdef OPENMPCODE
    if (parallel_build) {
#pragma omp parallel shared(build_status)
        {
#pragma omp single
            build_status = kdtree_2balls_build_node(
                cmd, tree, 0, valid_count - 1, leaf_capacity, 0,
                role, KDTREE_2BALLS_ROOT, &shape, parallel_build);
        }
    } else
#endif
    {
        build_status = kdtree_2balls_build_node_serial(
            cmd, tree, 0, valid_count - 1, leaf_capacity, 0,
            role, &root);
        if (root != KDTREE_2BALLS_ROOT) build_status = FAILURE;
    }
    if (build_status == FAILURE) {
        if (cmd->error_message[0] == '\0')
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     "kdtree-2balls: node construction failed");
        kdtree_2balls_tree_free(tree);
        free(shape.nodes_by_points);
        return FAILURE;
    }
    free(shape.nodes_by_points);
    shape.nodes_by_points = NULL;

    for (source = 0; source < valid_count; source++) {
        real field;
        real normalization;
        bodyptr body = tree->bptr[source];

        kdtree_2balls_point_values(cmd, body, role, &field, &normalization);
        SETV(tree->packed_points[source].pos, Pos(body));
        tree->packed_points[source].kappa = normalization != 0.0
            ? field / normalization : 0.0;
        tree->packed_points[source].weight = normalization;
        tree->packed_points[source].weighted_kappa = field;
        tree->packed_points[source].source = body;
    }

    gd->bytes_tot += sizeof(*tree)
        + (size_t)valid_count * sizeof(*tree->bptr)
        + (size_t)tree->capacity * sizeof(*tree->nodes)
        + (size_t)valid_count * sizeof(*tree->packed_points);
    *result = tree;
    return SUCCESS;

allocation_failure:
    snprintf(cmd->error_message, _ERRORMSGSIZE_,
             "kdtree-2balls: memory allocation failed");
    kdtree_2balls_tree_free(tree);
    free(shape.nodes_by_points);
    return FAILURE;
}

int kdtree_2balls_tree_build_pivot(
        struct cmdline_data *cmd, struct global_data *gd,
        bodyptr body_table, INTEGER body_count, int leaf_capacity,
        fcfc_balltreeptr *result)
{
    return kdtree_2balls_tree_build(
        cmd, gd, body_table, body_count, leaf_capacity,
        KDTREE_2BALLS_PIVOT, result);
}

int kdtree_2balls_tree_build_neighbor(
        struct cmdline_data *cmd, struct global_data *gd,
        bodyptr body_table, INTEGER body_count, int leaf_capacity,
        fcfc_balltreeptr *result)
{
    return kdtree_2balls_tree_build(
        cmd, gd, body_table, body_count, leaf_capacity,
        KDTREE_2BALLS_NEIGHBOR, result);
}

static int kdtree_2balls_minimum_leaf_depth(
        const fcfc_balltreeptr tree, INTEGER node_index, int depth)
{
    const fcfc_ballnode *node = &tree->nodes[node_index];

    if (node->left < 0) return depth;
    const int left = kdtree_2balls_minimum_leaf_depth(
        tree, node->left, depth + 1);
    const int right = kdtree_2balls_minimum_leaf_depth(
        tree, node->right, depth + 1);
    return MIN(left, right);
}

static void kdtree_2balls_gather_frontier(
        const fcfc_balltreeptr tree, INTEGER node_index,
        int depth, int requested_depth, INTEGER *frontier, INTEGER *count)
{
    const fcfc_ballnode *node = &tree->nodes[node_index];

    if (depth == requested_depth || node->left < 0) {
        frontier[(*count)++] = node_index;
        return;
    }
    kdtree_2balls_gather_frontier(
        tree, node->left, depth + 1, requested_depth, frontier, count);
    kdtree_2balls_gather_frontier(
        tree, node->right, depth + 1, requested_depth, frontier, count);
}

int kdtree_2balls_tree_frontier(
        struct cmdline_data *cmd, fcfc_balltreeptr tree,
        INTEGER minimum_count, INTEGER **result, INTEGER *result_count)
{
    INTEGER capacity = 1;
    INTEGER count = 0;
    INTEGER *frontier;
    int depth = 0;
    int maximum_depth;

    if (tree == NULL || tree->nodes == NULL || tree->nnode < 1
        || minimum_count < 1 || result == NULL || result_count == NULL) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "kdtree-2balls: invalid frontier arguments");
        return FAILURE;
    }
    *result = NULL;
    *result_count = 0;
    maximum_depth = kdtree_2balls_minimum_leaf_depth(
        tree, KDTREE_2BALLS_ROOT, 0);
    while (capacity < minimum_count && depth < maximum_depth) {
        if (capacity >
#ifdef LONGINT
            LONG_MAX / 2
#else
            INT_MAX / 2
#endif
            ) {
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     "kdtree-2balls: frontier size overflow");
            return FAILURE;
        }
        capacity *= 2;
        depth++;
    }
    if ((size_t)capacity > SIZE_MAX / sizeof(*frontier)) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "kdtree-2balls: frontier allocation overflow");
        return FAILURE;
    }
    frontier = malloc((size_t)capacity * sizeof(*frontier));
    if (frontier == NULL) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "kdtree-2balls: frontier allocation failed");
        return FAILURE;
    }
    kdtree_2balls_gather_frontier(
        tree, KDTREE_2BALLS_ROOT, 0, depth, frontier, &count);
    if (count <= 0) {
        free(frontier);
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "kdtree-2balls: empty task frontier");
        return FAILURE;
    }
    *result = frontier;
    *result_count = count;
    return SUCCESS;
}
