/* Octree forest discovery followed by exact per-forest radial-tree queries.
 * This is a neighbor backend for the 3D estimators, not a radial-only CF.
 */
#include "globaldefs.h"
#include "lya_forest_defs.h"
#include "lya_forest_los_tree.h"
#include <float.h>

#define LYA_LOS_LEAF_SIZE 8
#define LYA_LOS_MIXED SIZE_MAX

typedef struct {
    nodeptr node;
    size_t escape;
    size_t forest;
} lya_los_octant;

typedef struct {
    size_t first, end, escape;
} lya_los_radial_node;

typedef struct {
    size_t root;
    long double direction[NDIM], origin[NDIM];
    long double deviation, scale;
} lya_los_forest;

struct lya_los_index {
    bodyptr table;
    bodyptr *points;
    size_t *body_forest;
    lya_los_forest *forests;
    lya_los_octant *octants;
    lya_los_radial_node *radial;
    size_t forest_count, octant_count, radial_count;
};

static int lya_los_compare(const void *a, const void *b)
{
    bodyptr p = *(bodyptr const *)a, q = *(bodyptr const *)b;
    if (LyaForestId(p) != LyaForestId(q))
        return LyaForestId(p) < LyaForestId(q) ? -1 : 1;
    if (LyaDistance(p) != LyaDistance(q))
        return LyaDistance(p) < LyaDistance(q) ? -1 : 1;
    return Id(p) < Id(q) ? -1 : Id(p) > Id(q);
}

static size_t lya_los_radial_size(size_t count)
{
    if (count <= LYA_LOS_LEAF_SIZE) return 1;
    return 1 + lya_los_radial_size(count / 2)
             + lya_los_radial_size(count - count / 2);
}

static void lya_los_build_radial(lya_los_index *index, size_t first, size_t end)
{
    size_t id = index->radial_count++;
    lya_los_radial_node *node = index->radial + id;
    node->first = first;
    node->end = end;
    if (end - first > LYA_LOS_LEAF_SIZE) {
        size_t middle = first + (end - first) / 2;
        lya_los_build_radial(index, first, middle);
        lya_los_build_radial(index, middle, end);
    }
    node->escape = index->radial_count;
}

static void lya_los_pack_octree(lya_los_index *index, nodeptr q)
{
    size_t id = index->octant_count++;
    lya_los_octant *node = index->octants + id;
    node->node = q;
    if (Type(q) == CELL) {
        nodeptr child;
        int first = TRUE;
        node->forest = LYA_LOS_MIXED;
        for (child = More(q); child != Next(q); child = Next(child)) {
            size_t child_id = index->octant_count;
            lya_los_pack_octree(index, child);
            if (first) node->forest = index->octants[child_id].forest;
            else if (node->forest != index->octants[child_id].forest)
                node->forest = LYA_LOS_MIXED;
            first = FALSE;
        }
    } else {
        node->forest = index->body_forest[(bodyptr)q - index->table];
    }
    node->escape = index->octant_count;
}

void lya_los_free(lya_los_index *index)
{
    if (index == NULL) return;
    free(index->points);
    free(index->body_forest);
    free(index->forests);
    free(index->octants);
    free(index->radial);
    free(index);
}

size_t lya_los_forest_count(const lya_los_index *index)
{
    return index->forest_count;
}

int lya_los_build(lya_los_index **result, bodyptr table, INTEGER count,
                  nodeptr root, ErrorMsg error_message)
{
    lya_los_index *index = NULL;
    size_t i, first, n, nodes = 0, radial_nodes = 0;
    nodeptr q;
    *result = NULL;
    if (count <= 0 || (uintmax_t)count > SIZE_MAX || root == NULL) {
        snprintf(error_message, _ERRORMSGSIZE_, "invalid LOS-tree catalog");
        return FAILURE;
    }
    n = (size_t)count;
#define LYA_LOS_ALLOC(pointer, length) \
    if (cballs_calloc_checked((void **)&(pointer), (length), sizeof(*(pointer)), \
        "Ly-alpha LOS tree", error_message, _ERRORMSGSIZE_) == FAILURE) goto fail
    LYA_LOS_ALLOC(index, 1);
    index->table = table;
    LYA_LOS_ALLOC(index->points, n);
    LYA_LOS_ALLOC(index->body_forest, n);
    for (i = 0; i < n; i++) index->points[i] = table + i;
    qsort(index->points, n, sizeof(*index->points), lya_los_compare);
    for (first = 0; first < n; first = i) {
        for (i = first + 1; i < n
             && LyaForestId(index->points[i]) == LyaForestId(index->points[first]); i++);
        index->forest_count++;
        /* With leaves of at least four pixels, the node count is <= n. */
        radial_nodes += lya_los_radial_size(i - first);
    }
    LYA_LOS_ALLOC(index->forests, index->forest_count);
    LYA_LOS_ALLOC(index->radial, radial_nodes);
    size_t forest = 0;
    for (first = 0; first < n; first = i, forest++) {
        lya_los_forest *f = index->forests + forest;
        bodyptr reference = index->points[first];
        long double norm = 0.0L;
        int axis;
        for (axis = 0; axis < NDIM; axis++)
            norm += (long double)LyaLOS(reference)[axis] * LyaLOS(reference)[axis];
        norm = sqrtl(norm);
        f->scale = 1.0L;
        for (axis = 0; axis < NDIM; axis++) {
            f->direction[axis] = LyaLOS(reference)[axis] / norm;
            /* Use actual stored positions: Pos() may have been recentered or
             * rounded to float. The residual bound below covers both. */
            f->origin[axis] = (long double)Pos(reference)[axis]
                           - LyaDistance(reference) * f->direction[axis];
        }
        for (i = first; i < n
             && LyaForestId(index->points[i]) == LyaForestId(reference); i++) {
            bodyptr p = index->points[i];
            long double residual2 = 0.0L;
            index->body_forest[p - table] = forest;
            f->scale = fmaxl(f->scale, fabsl(LyaDistance(p)));
            for (axis = 0; axis < NDIM; axis++) {
                long double residual = (long double)Pos(p)[axis] - f->origin[axis]
                                     - LyaDistance(p) * f->direction[axis];
                residual2 += residual * residual;
                f->scale = fmaxl(f->scale, fabsl(Pos(p)[axis]));
            }
            f->deviation = fmaxl(f->deviation, sqrtl(residual2));
        }
        f->root = index->radial_count;
        lya_los_build_radial(index, first, i);
    }
    for (q = root; q != Next(root); q = Type(q) == CELL ? More(q) : Next(q)) {
        if (nodes == SIZE_MAX) {
            snprintf(error_message, _ERRORMSGSIZE_, "LOS octree size overflow");
            goto fail;
        }
        nodes++;
    }
    LYA_LOS_ALLOC(index->octants, nodes);
    lya_los_pack_octree(index, root);
#undef LYA_LOS_ALLOC
    *result = index;
    return SUCCESS;
fail:
    lya_los_free(index);
    return FAILURE;
}

int lya_los_workspace_init(const lya_los_index *index,
                            lya_los_workspace *workspace,
                            ErrorMsg error_message)
{
    memset(workspace, 0, sizeof(*workspace));
    if (cballs_calloc_checked((void **)&workspace->seen, index->forest_count,
            sizeof(*workspace->seen), "LOS forest stamps", error_message,
            _ERRORMSGSIZE_) == FAILURE
        || cballs_calloc_checked((void **)&workspace->forests, index->forest_count,
            sizeof(*workspace->forests), "LOS forest frontier", error_message,
            _ERRORMSGSIZE_) == FAILURE) {
        lya_los_workspace_free(workspace);
        return FAILURE;
    }
    return SUCCESS;
}

void lya_los_workspace_free(lya_los_workspace *workspace)
{
    free(workspace->seen);
    free(workspace->forests);
    memset(workspace, 0, sizeof(*workspace));
}

static void lya_los_interval(const lya_los_forest *forest, bodyptr pivot,
                              REAL cutoff, long double *lower, long double *upper)
{
    long double d[NDIM], a = 0.0L, b2 = 0.0L;
    long double scale = forest->scale + fabsl(cutoff), padding, radius, half;
    int axis;
    for (axis = 0; axis < NDIM; axis++) {
        d[axis] = (long double)Pos(pivot)[axis] - forest->origin[axis];
        a += d[axis] * forest->direction[axis];
        scale += fabsl(d[axis]) + fabsl(Pos(pivot)[axis]);
    }
    /* Every actual pixel is within deviation of origin + chi * direction.
     * Inflating the sphere therefore remains conservative even for a forest
     * with inconsistent sightlines. Exact 3D checks are still made per pixel. */
    padding = 128.0L * DBL_EPSILON * (scale + forest->deviation);
    radius = (long double)cutoff + forest->deviation + padding;
    for (axis = 0; axis < NDIM; axis++) {
        long double perpendicular = d[axis] - a * forest->direction[axis];
        b2 += perpendicular * perpendicular;
    }
    /* A witness was already found inside the sphere. In the exceptional case
     * of inconsistent numeric bounds, retain the whole forest. */
    if (!isfinite(radius) || !isfinite(b2) || b2 > radius * radius) {
        *lower = -LDBL_MAX;
        *upper = LDBL_MAX;
        return;
    }
    half = sqrtl(fmaxl(0.0L, radius * radius - b2));
    *lower = a - half - padding;
    *upper = a + half + padding;
}

int lya_los_query(const lya_los_index *index, lya_los_workspace *workspace,
                  bodyptr pivot, REAL cutoff, lya_los_visit visit,
                  void *context, ErrorMsg error_message)
{
    size_t node = 0, found = 0, f;
    const size_t own = index->body_forest[pivot - index->table];
    if (++workspace->epoch == 0) {
        memset(workspace->seen, 0, index->forest_count * sizeof(*workspace->seen));
        workspace->epoch = 1;
    }
    workspace->seen[own] = workspace->epoch;
    while (node < index->octant_count) {
        const lya_los_octant *entry = index->octants + node;
        nodeptr q = entry->node;
        REAL d2, distance;
        compute_vector displacement;
        workspace->octree_nodes++;
        if (entry->forest != LYA_LOS_MIXED
            && workspace->seen[entry->forest] == workspace->epoch) {
            workspace->forest_skips++;
            node = entry->escape;
            continue;
        }
        DOTPSUBV(d2, displacement, Pos(pivot), Pos(q));
        distance = rsqrt(d2);
        if (Type(q) == CELL) {
            if (distance >= cutoff + Size(q) * rsqrt((REAL)NDIM))
                node = entry->escape;
            else node++;
        } else {
            if (distance < cutoff && Update(q) != FALSE && Mask(q) == MASK_NODE_VALID) {
                workspace->seen[entry->forest] = workspace->epoch;
                workspace->forests[found++] = entry->forest;
            }
            node++;
        }
    }
    workspace->forest_hits += found;
    for (f = 0; f < found; f++) {
        const lya_los_forest *forest = index->forests + workspace->forests[f];
        long double lower, upper;
        size_t end = index->radial[forest->root].escape;
        lya_los_interval(forest, pivot, cutoff, &lower, &upper);
        node = forest->root;
        while (node < end) {
            const lya_los_radial_node *entry = index->radial + node;
            workspace->radial_nodes++;
            if (LyaDistance(index->points[entry->end - 1]) < lower
                || LyaDistance(index->points[entry->first]) > upper) {
                node = entry->escape;
                continue;
            }
            if (entry->end - entry->first <= LYA_LOS_LEAF_SIZE) {
                size_t i;
                for (i = entry->first; i < entry->end; i++) {
                    bodyptr q = index->points[i];
                    REAL d2, distance;
                    compute_vector displacement;
                    if (LyaDistance(q) < lower || LyaDistance(q) > upper
                        || Update(q) == FALSE || Mask(q) != MASK_NODE_VALID) continue;
                    workspace->pixel_tests++;
                    DOTPSUBV(d2, displacement, Pos(pivot), Pos(q));
                    distance = rsqrt(d2);
                    if (distance < cutoff
                        && visit(q, distance, context, error_message) == FAILURE)
                        return FAILURE;
                }
            }
            node++;
        }
    }
    return SUCCESS;
}
