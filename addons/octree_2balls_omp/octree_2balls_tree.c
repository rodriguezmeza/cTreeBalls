/* Binary traversal view of the native cTreeBalls octree. */

#include <limits.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "globaldefs.h"
#include "octree_2balls_tree.h"

#define OCTREE_2BALLS_ROOT 0
#define OCTREE_2BALLS_PARALLEL_BUILD_CUTOFF ((INTEGER)32768)
#define OCTREE_2BALLS_CACHE_SLOTS 2
#define OCTREE_2BALLS_HASH_CHUNKS 64
#define OCTREE_2BALLS_PARALLEL_HASH_CUTOFF ((INTEGER)262144)

typedef struct {
    nodeptr node;
    INTEGER count;
} octree_2balls_child;

typedef struct {
    uint64_t fingerprint;
    uint64_t stamp;
    INTEGER nbody;
    int leaf_capacity;
    bool read_mask;
    int users;
    fcfc_balltreeptr tree;
} octree_2balls_cache_entry;

static octree_2balls_cache_entry
    octree_2balls_cache[OCTREE_2BALLS_CACHE_SLOTS];
static uint64_t octree_2balls_cache_stamp;
static bool octree_2balls_cache_registered;

static real octree_2balls_field(bodyptr p)
{
#ifdef KappaAvgON
    return KappaAvg(p);
#else
    return Kappa(p);
#endif
}

static uint64_t octree_2balls_hash_word(uint64_t hash, uint64_t value)
{
    value ^= value >> 30;
    value *= UINT64_C(0xbf58476d1ce4e5b9);
    value ^= value >> 27;
    value *= UINT64_C(0x94d049bb133111eb);
    value ^= value >> 31;
    return (hash ^ value) * UINT64_C(1099511628211);
}

static uint64_t octree_2balls_hash_real(uint64_t hash, real value)
{
    uint64_t bits = 0;

    memcpy(&bits, &value, sizeof(value));
    return octree_2balls_hash_word(hash, bits);
}

static uint64_t octree_2balls_catalog_range_fingerprint(
        bodyptr btab, INTEGER first, INTEGER last, bool read_mask,
        uint64_t seed)
{
    uint64_t hash = seed;

    for (INTEGER i = first; i < last; i++) {
        bodyptr body = btab + i;
        const real mass = Mass(body);
        const real field = octree_2balls_field(body);
        const real weight = Weight(body);

        for (int k = 0; k < NDIM; k++) {
            const real position = (real)Pos(body)[k];

            hash = octree_2balls_hash_real(hash, position);
        }
        hash = octree_2balls_hash_real(hash, mass);
        hash = octree_2balls_hash_real(hash, field);
        hash = octree_2balls_hash_real(hash, weight);
        if (read_mask) {
            const int mask = Mask(body);

            hash = octree_2balls_hash_word(hash, (uint64_t)(unsigned int)mask);
        }
    }
    return hash;
}

static uint64_t octree_2balls_catalog_fingerprint(
        struct cmdline_data *cmd, bodyptr btab, INTEGER nbody,
        int leaf_capacity)
{
    uint64_t partial[OCTREE_2BALLS_HASH_CHUNKS];
    uint64_t hash = UINT64_C(1469598103934665603);
    const bool read_mask = cballs_opt_read_mask(cmd);
    const int chunks = (int)MIN(
        (INTEGER)OCTREE_2BALLS_HASH_CHUNKS, nbody);
    const INTEGER base = nbody / chunks;
    const INTEGER remainder = nbody % chunks;

    hash = octree_2balls_hash_word(hash, (uint64_t)nbody);
    hash = octree_2balls_hash_word(hash, (uint64_t)leaf_capacity);
    hash = octree_2balls_hash_word(hash, (uint64_t)read_mask);
#ifdef OPENMPCODE
#pragma omp parallel for schedule(static) \
    if(nbody >= OCTREE_2BALLS_PARALLEL_HASH_CUTOFF && cmd->numthreads > 1)
#endif
    for (int chunk = 0; chunk < chunks; chunk++) {
        const INTEGER first = (INTEGER)chunk * base
            + MIN((INTEGER)chunk, remainder);
        const INTEGER count = base + ((INTEGER)chunk < remainder);
        const uint64_t seed = octree_2balls_hash_word(
            UINT64_C(1469598103934665603), (uint64_t)chunk);

        partial[chunk] = octree_2balls_catalog_range_fingerprint(
            btab, first, first + count, read_mask, seed);
    }
    for (int chunk = 0; chunk < chunks; chunk++)
        hash = octree_2balls_hash_word(hash, partial[chunk]);
    return hash;
}

static void octree_2balls_cache_clear(void)
{
    for (int i = 0; i < OCTREE_2BALLS_CACHE_SLOTS; i++) {
        octree_2balls_tree_free(octree_2balls_cache[i].tree);
        octree_2balls_cache[i].tree = NULL;
        octree_2balls_cache[i].users = 0;
    }
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

static void octree_2balls_link_parent(fcfc_balltreeptr tree,
                                      INTEGER parent_index,
                                      INTEGER left_index,
                                      INTEGER right_index)
{
    fcfc_ballnode *parent = &tree->nodes[parent_index];
    const fcfc_ballnode *left = &tree->nodes[left_index];
    const fcfc_ballnode *right = &tree->nodes[right_index];

    parent->first = left->first;
    parent->last = right->last;
    parent->left = left_index;
    parent->right = right_index;
}

static void octree_2balls_finish_parent_geometry(struct cmdline_data *cmd,
                                                 fcfc_balltreeptr tree,
                                                 INTEGER parent_index)
{
    fcfc_ballnode *parent = &tree->nodes[parent_index];
    const INTEGER left_index = parent->left;
    const INTEGER right_index = parent->right;
    const fcfc_ballnode *left = &tree->nodes[left_index];
    const fcfc_ballnode *right = &tree->nodes[right_index];
    const INTEGER left_count = left->last - left->first + 1;
    const INTEGER right_count = right->last - right->first + 1;
    const real mass = left->weight + right->weight;
    real radius_squared = 0.0;
    INTEGER point;
    int k;

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
                parent->center, tree->packed_points[point].pos));
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

static void octree_2balls_finish_leaf_geometry(struct cmdline_data *cmd,
                                               fcfc_balltreeptr tree,
                                               INTEGER index)
{
    compute_vector cmpos_sum;
    compute_vector geometric_center;
    fcfc_ballnode *target = &tree->nodes[index];
    real farthest2 = 0.0;
    int k;

    CLRV(cmpos_sum);
    CLRV(geometric_center);
    target->weight = 0.0;
    target->kappa_sum = 0.0;
    target->kappa_sq_sum = 0.0;
    target->field_weight_sum = 0.0;
    target->field_weight_sq_sum = 0.0;
    target->weighted_kappa_sum = 0.0;
    target->weighted_kappa_sq_sum = 0.0;
    for (INTEGER point = target->first; point <= target->last; point++) {
        bodyptr body = tree->bptr[point];
        fcfc_ballpoint *packed = &tree->packed_points[point];
        const real mass = Mass(body);
        const real field = octree_2balls_field(body);
        const real field_weight = Weight(body);

        SETV(packed->pos, Pos(body));
        packed->kappa = field;
        packed->weight = field_weight;
        packed->weighted_kappa = field_weight * field;
        packed->source = body;
        DO_COORD(k) {
            cmpos_sum[k] += mass * (real)packed->pos[k];
            geometric_center[k] += (real)packed->pos[k];
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
             : geometric_center[k]
               / (real)(target->last - target->first + 1));
    target->kappa = target->kappa_sum
                  / (real)(target->last - target->first + 1);
    SETV(target->center, target->cmpos);
    for (INTEGER point = target->first; point <= target->last; point++)
        farthest2 = MAX(farthest2, octree_2balls_distance_squared(
            target->center, tree->packed_points[point].pos));
    target->radius = cballs_store_search_bound(rsqrt(farthest2));
    target->aggregate_radius = cmd->theta > 0.0
        ? cballs_store_search_bound(rsqrt(farthest2) / cmd->theta)
        : cballs_store_upper_bound(MAX_REAL_NUMBER);
}

static void octree_2balls_finish_subtree(struct cmdline_data *cmd,
                                         fcfc_balltreeptr tree,
                                         INTEGER index)
{
    fcfc_ballnode *node = &tree->nodes[index];

    if (node->left < 0) {
        octree_2balls_finish_leaf_geometry(cmd, tree, index);
        return;
    }
#ifdef OPENMPCODE
    if (node->last - node->first + 1 >= OCTREE_2BALLS_PARALLEL_BUILD_CUTOFF) {
        const INTEGER left = node->left;
        const INTEGER right = node->right;

#pragma omp taskgroup
        {
#pragma omp task firstprivate(left) shared(cmd, tree)
            octree_2balls_finish_subtree(cmd, tree, left);
#pragma omp task firstprivate(right) shared(cmd, tree)
            octree_2balls_finish_subtree(cmd, tree, right);
        }
    } else
#endif
    {
        octree_2balls_finish_subtree(cmd, tree, node->left);
        octree_2balls_finish_subtree(cmd, tree, node->right);
    }
    octree_2balls_finish_parent_geometry(cmd, tree, index);
}

static int octree_2balls_build_native(struct cmdline_data *,
                                      fcfc_balltreeptr, nodeptr, int, int,
                                      bool, INTEGER *);

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
                                     int leaf_capacity, bool finish_inline,
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
            cmd, tree, children[first].node, depth, leaf_capacity,
            finish_inline, result);
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
                                  depth + 1, leaf_capacity, finish_inline,
                                  &left_index) == FAILURE
        || octree_2balls_build_group(cmd, tree, children, split + 1, last,
                                     depth + 1, leaf_capacity, finish_inline,
                                     &right_index) == FAILURE)
        return FAILURE;
    octree_2balls_link_parent(tree, parent_index, left_index, right_index);
    if (finish_inline)
        octree_2balls_finish_parent_geometry(cmd, tree, parent_index);
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
                                    bool finish_inline, INTEGER *result)
{
    fcfc_ballnode *target;
    INTEGER index;

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
    if (finish_inline)
        octree_2balls_finish_leaf_geometry(cmd, tree, index);
    *result = index;
    return SUCCESS;
}

static int octree_2balls_build_native(struct cmdline_data *cmd,
                                      fcfc_balltreeptr tree, nodeptr source,
                                      int depth, int leaf_capacity,
                                      bool finish_inline, INTEGER *result)
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
            cmd, tree, source, source_count, depth, finish_inline, result);
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
        leaf_capacity, finish_inline, result);
}

int octree_2balls_tree_build(struct cmdline_data *cmd,
                             struct global_data *gd, bodyptr btab,
                             INTEGER nbody, int leaf_capacity,
                             fcfc_balltreeptr *result)
{
    fcfc_balltreeptr tree = NULL;
    INTEGER root = -1;
    INTEGER live_count;
    bool finish_inline;
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
    tree->packed_points = malloc(
        (size_t)live_count * sizeof(*tree->packed_points));
    if (tree->bptr == NULL || tree->nodes == NULL
        || tree->packed_points == NULL)
        goto allocation_failure;
#ifdef OPENMPCODE
    finish_inline = live_count < OCTREE_2BALLS_PARALLEL_BUILD_CUTOFF
                 || cmd->numthreads <= 1;
#else
    finish_inline = TRUE;
#endif
    if (octree_2balls_build_native(cmd, tree, (nodeptr)roottable[catalog],
                                   0, leaf_capacity, finish_inline,
                                   &root) == FAILURE)
        goto failure;
    if (root != OCTREE_2BALLS_ROOT || tree->npoint != live_count) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "octree-2balls-omp: invalid binary octree root");
        goto failure;
    }

    if (!finish_inline) {
#ifdef OPENMPCODE
#pragma omp parallel
        {
#pragma omp single nowait
            octree_2balls_finish_subtree(cmd, tree, root);
        }
#else
        octree_2balls_finish_subtree(cmd, tree, root);
#endif
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

int octree_2balls_tree_build_cached(struct cmdline_data *cmd,
                                    struct global_data *gd, bodyptr btab,
                                    INTEGER nbody, int leaf_capacity,
                                    fcfc_balltreeptr *result, bool *cache_hit)
{
    uint64_t fingerprint;
    bool read_mask;
    fcfc_balltreeptr built = NULL;
    fcfc_balltreeptr found = NULL;
    int install_slot = -1;
    bool reused = FALSE;

    if (cmd == NULL || btab == NULL || nbody < 1 || leaf_capacity < 1
        || result == NULL || cache_hit == NULL) {
        if (cmd == NULL) return FAILURE;
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "octree-2balls-omp: invalid compact-tree cache result");
        return FAILURE;
    }
    *result = NULL;
    *cache_hit = FALSE;
    read_mask = cballs_opt_read_mask(cmd);
    fingerprint = octree_2balls_catalog_fingerprint(
        cmd, btab, nbody, leaf_capacity);

#ifdef OPENMPCODE
#pragma omp critical(octree_2balls_tree_cache)
#endif
    {
        for (int i = 0; i < OCTREE_2BALLS_CACHE_SLOTS; i++) {
            octree_2balls_cache_entry *entry = &octree_2balls_cache[i];

            if (entry->tree != NULL
                && entry->fingerprint == fingerprint
                && entry->nbody == nbody
                && entry->leaf_capacity == leaf_capacity
                && entry->read_mask == read_mask) {
                entry->users++;
                entry->stamp = ++octree_2balls_cache_stamp;
                found = entry->tree;
                reused = TRUE;
                break;
            }
        }
    }
    if (found != NULL) {
        *result = found;
        *cache_hit = TRUE;
        return SUCCESS;
    }

    if (octree_2balls_tree_build(
            cmd, gd, btab, nbody, leaf_capacity, &built) == FAILURE)
        return FAILURE;

#ifdef OPENMPCODE
#pragma omp critical(octree_2balls_tree_cache)
#endif
    {
        uint64_t oldest_stamp = UINT64_MAX;

        /* A concurrent caller may have installed the same catalog while this
         * thread built its private candidate. */
        for (int i = 0; i < OCTREE_2BALLS_CACHE_SLOTS; i++) {
            octree_2balls_cache_entry *entry = &octree_2balls_cache[i];

            if (entry->tree != NULL
                && entry->fingerprint == fingerprint
                && entry->nbody == nbody
                && entry->leaf_capacity == leaf_capacity
                && entry->read_mask == read_mask) {
                entry->users++;
                entry->stamp = ++octree_2balls_cache_stamp;
                found = entry->tree;
                reused = TRUE;
                break;
            }
            if (entry->users == 0
                && (entry->tree == NULL || entry->stamp < oldest_stamp)) {
                install_slot = i;
                oldest_stamp = entry->tree == NULL ? 0 : entry->stamp;
            }
        }
        if (found == NULL && install_slot >= 0) {
            octree_2balls_cache_entry *entry =
                &octree_2balls_cache[install_slot];

            octree_2balls_tree_free(entry->tree);
            free(built->bptr);
            built->bptr = NULL;
            for (INTEGER i = 0; i < built->npoint; i++)
                built->packed_points[i].source = NULL;
            entry->fingerprint = fingerprint;
            entry->stamp = ++octree_2balls_cache_stamp;
            entry->nbody = nbody;
            entry->leaf_capacity = leaf_capacity;
            entry->read_mask = read_mask;
            entry->users = 1;
            entry->tree = built;
            found = built;
            built = NULL;
            if (!octree_2balls_cache_registered) {
                atexit(octree_2balls_cache_clear);
                octree_2balls_cache_registered = TRUE;
            }
        }
    }

    if (built != NULL && found != NULL) octree_2balls_tree_free(built);
    if (found != NULL) {
        *result = found;
        *cache_hit = reused;
    } else {
        *result = built;
    }
    return SUCCESS;
}

bool octree_2balls_tree_cache_contains(struct cmdline_data *cmd,
                                       bodyptr btab, INTEGER nbody,
                                       int leaf_capacity)
{
    uint64_t fingerprint;
    bool found = FALSE;
    bool read_mask;

    if (cmd == NULL || btab == NULL || nbody < 1 || leaf_capacity < 1)
        return FALSE;
    read_mask = cballs_opt_read_mask(cmd);
    fingerprint = octree_2balls_catalog_fingerprint(
        cmd, btab, nbody, leaf_capacity);
#ifdef OPENMPCODE
#pragma omp critical(octree_2balls_tree_cache)
#endif
    {
        for (int i = 0; i < OCTREE_2BALLS_CACHE_SLOTS; i++) {
            const octree_2balls_cache_entry *entry = &octree_2balls_cache[i];

            if (entry->tree != NULL
                && entry->fingerprint == fingerprint
                && entry->nbody == nbody
                && entry->leaf_capacity == leaf_capacity
                && entry->read_mask == read_mask) {
                found = TRUE;
                break;
            }
        }
    }
    return found;
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

void octree_2balls_tree_release(fcfc_balltreeptr tree)
{
    bool cached = FALSE;

    if (tree == NULL) return;
#ifdef OPENMPCODE
#pragma omp critical(octree_2balls_tree_cache)
#endif
    {
        for (int i = 0; i < OCTREE_2BALLS_CACHE_SLOTS; i++) {
            octree_2balls_cache_entry *entry = &octree_2balls_cache[i];

            if (entry->tree != tree) continue;
            if (entry->users > 0) entry->users--;
            cached = TRUE;
            break;
        }
    }
    if (!cached) octree_2balls_tree_free(tree);
}
