/* Median KD tree with full-sky spin-2 moments. */

#include <float.h>
#include <limits.h>
#include <stdint.h>
#include <stdlib.h>

#include "globaldefs.h"
#include "kdtree_shear_sphere_tree.h"

typedef struct {
    real re;
    real im;
} kd_shear_complex;

static kd_shear_complex kd_shear_make(real re, real im)
{
    kd_shear_complex value = {re, im};
    return value;
}

static kd_shear_complex kd_shear_add(kd_shear_complex a, kd_shear_complex b)
{
    return kd_shear_make(a.re + b.re, a.im + b.im);
}

static kd_shear_complex kd_shear_mul(kd_shear_complex a, kd_shear_complex b)
{
    return kd_shear_make(a.re*b.re - a.im*b.im,
                         a.re*b.im + a.im*b.re);
}

static kd_shear_complex kd_shear_scale(kd_shear_complex value, real scale)
{
    return kd_shear_make(value.re*scale, value.im*scale);
}

static real kd_shear_dot3(const real *a, const real *b)
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

static bool kd_shear_unit3(const real *position, compute_vector unit)
{
    const real norm2 = kd_shear_dot3(position, position);
    real inverse_norm;

    if (!(norm2 > 0.0) || !isfinite(norm2)) return FALSE;
    inverse_norm = 1.0/rsqrt(norm2);
    unit[0] = position[0]*inverse_norm;
    unit[1] = position[1]*inverse_norm;
    unit[2] = position[2]*inverse_norm;
    return isfinite(unit[0]) && isfinite(unit[1]) && isfinite(unit[2]);
}

static bool kd_shear_basis(const real *unit, compute_vector east,
                           compute_vector north)
{
    const real equatorial_norm = rsqrt(unit[0]*unit[0] + unit[1]*unit[1]);

    if (equatorial_norm > 64.0*DBL_EPSILON) {
        east[0] = -unit[1]/equatorial_norm;
        east[1] = unit[0]/equatorial_norm;
        east[2] = 0.0;
    } else {
        east[0] = 1.0;
        east[1] = 0.0;
        east[2] = 0.0;
    }
    north[0] = unit[1]*east[2] - unit[2]*east[1];
    north[1] = unit[2]*east[0] - unit[0]*east[2];
    north[2] = unit[0]*east[1] - unit[1]*east[0];
    return isfinite(north[0]) && isfinite(north[1]) && isfinite(north[2]);
}

static bool kd_shear_rotation(const real *target_unit,
                              const real *target_east,
                              const real *target_north,
                              const real *source_position,
                              kd_shear_complex *rotation)
{
    compute_vector source_unit;
    compute_vector source_east;
    compute_vector source_north;
    compute_vector transported_east;
    compute_vector transported_north;
    real denominator;
    real c;
    real s;
    real norm;
    int axis;

    if (!kd_shear_unit3(source_position, source_unit)
        || !kd_shear_basis(source_unit, source_east, source_north))
        return FALSE;
    (void)source_north;
    denominator = 1.0 + kd_shear_dot3(target_unit, source_unit);
    if (!(denominator > 64.0*DBL_EPSILON) || !isfinite(denominator))
        return FALSE;
    for (axis = 0; axis < 3; axis++) {
        transported_east[axis] = target_east[axis]
            - kd_shear_dot3(target_east, source_unit)/denominator
              *(target_unit[axis] + source_unit[axis]);
        transported_north[axis] = target_north[axis]
            - kd_shear_dot3(target_north, source_unit)/denominator
              *(target_unit[axis] + source_unit[axis]);
    }
    c = kd_shear_dot3(source_east, transported_east);
    s = kd_shear_dot3(source_east, transported_north);
    norm = rsqrt(c*c + s*s);
    if (!(norm > 64.0*DBL_EPSILON) || !isfinite(norm)) return FALSE;
    c /= norm;
    s /= norm;
    *rotation = kd_shear_make(c*c - s*s, 2.0*c*s);
    return isfinite(rotation->re) && isfinite(rotation->im);
}

static bool kd_shear_valid(const struct cmdline_data *cmd, bodyptr body,
                           bool pivot_role)
{
    if (cballs_opt_read_mask(cmd) && Mask(body) == MASK_NODE_MASKED)
        return FALSE;
#ifdef SMOOTHPIVOT
    if (pivot_role && cballs_opt_smooth_pivot(cmd) && !Update(body))
        return FALSE;
#else
    (void)pivot_role;
#endif
    return TRUE;
}

static real kd_shear_body_weight(const struct cmdline_data *cmd, bodyptr body,
                                 bool pivot_role)
{
#ifdef SMOOTHPIVOT
    if (pivot_role && cballs_opt_smooth_pivot(cmd))
        return WeightRmin(body);
#else
    (void)cmd;
    (void)pivot_role;
#endif
    return Weight(body);
}

static kd_shear_complex kd_shear_body_weighted_gamma(
        const struct cmdline_data *cmd, bodyptr body, bool pivot_role)
{
#ifdef SMOOTHPIVOT
    if (pivot_role && cballs_opt_smooth_pivot(cmd))
        return kd_shear_make(Gamma1Rmin(body), Gamma2Rmin(body));
#else
    (void)cmd;
    (void)pivot_role;
#endif
    return kd_shear_scale(kd_shear_make(Gamma1(body), Gamma2(body)),
                          Weight(body));
}

static int kd_shear_before(bodyptr left, bodyptr right, int axis)
{
    const real a = (real)Pos(left)[axis];
    const real b = (real)Pos(right)[axis];

    if (a < b) return TRUE;
    if (a > b) return FALSE;
    return (uintptr_t)left < (uintptr_t)right;
}

static void kd_shear_swap(bodyptr *left, bodyptr *right)
{
    bodyptr temporary = *left;
    *left = *right;
    *right = temporary;
}

static void kd_shear_select(bodyptr *points, INTEGER lo, INTEGER hi,
                            INTEGER target, int axis)
{
    while (lo < hi) {
        INTEGER i = lo;
        INTEGER j = hi;
        bodyptr pivot = points[lo + (hi - lo)/2];

        while (i <= j) {
            while (i <= hi && kd_shear_before(points[i], pivot, axis)) i++;
            while (j >= lo && kd_shear_before(pivot, points[j], axis)) j--;
            if (i <= j) {
                kd_shear_swap(&points[i], &points[j]);
                i++;
                j--;
            }
        }
        if (target <= j) hi = j;
        else if (target >= i) lo = i;
        else return;
    }
}

static real kd_shear_distance2(const real *left, const real *right)
{
    real value = 0.0;
    int axis;
    DO_COORD(axis) value += rsqr(left[axis] - right[axis]);
    return value;
}

static int kd_shear_aggregate(struct cmdline_data *cmd,
                              fcfc_balltreeptr tree, fcfc_ballnode *node,
                              bool pivot_role)
{
    compute_vector center_sum;
    compute_vector center_unit;
    compute_vector center_east;
    compute_vector center_north;
    real center_weight = 0.0;
    real radius2 = 0.0;
    INTEGER point;
    int axis;

    CLRV(center_sum);
    for (point = node->first; point <= node->last; point++) {
        bodyptr body = tree->bptr[point];
        real weight = kd_shear_body_weight(cmd, body, pivot_role);
        real position_weight = weight > 0.0 ? weight : 1.0;

        DO_COORD(axis) center_sum[axis] += position_weight*Pos(body)[axis];
        center_weight += position_weight;
    }
    if (!(center_weight > 0.0)) return FAILURE;
    DO_COORD(axis) center_sum[axis] /= center_weight;
    if (!kd_shear_unit3(center_sum, center_unit))
        if (!kd_shear_unit3(Pos(tree->bptr[node->first]), center_unit))
            return FAILURE;
    if (!kd_shear_basis(center_unit, center_east, center_north)) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "kdtree-shear-sphere-2balls-omp: node tangent frame is undefined");
        return FAILURE;
    }
    DO_COORD(axis) {
        node->center[axis] = (cballs_storage_real)center_unit[axis];
        node->cmpos[axis] = (cballs_storage_real)center_unit[axis];
    }

    node->weight = 0.0;
    node->shear_gamma_re = 0.0;
    node->shear_gamma_im = 0.0;
    node->shear_gamma2_re = 0.0;
    node->shear_gamma2_im = 0.0;
    node->shear_gamma_abs2 = 0.0;
    node->shear_weight2 = 0.0;
    node->shear_transport_error = 0.0;
    for (point = node->first; point <= node->last; point++) {
        bodyptr body = tree->bptr[point];
        const real weight = kd_shear_body_weight(cmd, body, pivot_role);
        kd_shear_complex weighted_gamma =
            kd_shear_body_weighted_gamma(cmd, body, pivot_role);
        kd_shear_complex rotation;
        kd_shear_complex transported;

        if (!(weight >= 0.0) || !isfinite(weight)
            || !kd_shear_rotation(center_unit, center_east, center_north,
                                  Pos(body), &rotation)) {
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     "kdtree-shear-sphere-2balls-omp: invalid spin-2 node member");
            return FAILURE;
        }
        transported = kd_shear_mul(weighted_gamma, rotation);
        node->weight += weight;
        node->shear_gamma_re += transported.re;
        node->shear_gamma_im += transported.im;
        transported = kd_shear_mul(transported, transported);
        node->shear_gamma2_re += transported.re;
        node->shear_gamma2_im += transported.im;
        node->shear_gamma_abs2 += weighted_gamma.re*weighted_gamma.re
                                + weighted_gamma.im*weighted_gamma.im;
        node->shear_weight2 += weight*weight;
        radius2 = MAX(radius2, kd_shear_distance2(center_unit, Pos(body)));
    }
    node->radius = cballs_store_search_bound(rsqrt(radius2));
    node->aggregate_radius = node->radius;
    return SUCCESS;
}

static int kd_shear_build_node(struct cmdline_data *cmd,
                               fcfc_balltreeptr tree, INTEGER first,
                               INTEGER last, int leaf_capacity, int depth,
                               bool pivot_role, INTEGER *result)
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
                 "kdtree-shear-sphere-2balls-omp: node capacity exceeded");
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

    if (kd_shear_aggregate(cmd, tree, node, pivot_role) == FAILURE)
        return FAILURE;
    if (last - first + 1 > leaf_capacity) {
        const INTEGER middle = first + (last - first + 1)/2;

        kd_shear_select(tree->bptr, first, last, middle, split_axis);
        if (kd_shear_build_node(cmd, tree, first, middle - 1,
                                leaf_capacity, depth + 1, pivot_role,
                                &node->left) == FAILURE
            || kd_shear_build_node(cmd, tree, middle, last,
                                   leaf_capacity, depth + 1, pivot_role,
                                   &node->right) == FAILURE)
            return FAILURE;
    }
    *result = index;
    return SUCCESS;
}

void kdtree_shear_sphere_free(fcfc_balltreeptr tree)
{
    if (tree == NULL) return;
    free(tree->packed_points);
    free(tree->nodes);
    free(tree->bptr);
    free(tree);
}

int kdtree_shear_sphere_build(
        struct cmdline_data *cmd, struct global_data *gd,
        bodyptr body_table, INTEGER body_count, int leaf_capacity,
        bool pivot_role, fcfc_balltreeptr *result)
{
    fcfc_balltreeptr tree = NULL;
    INTEGER valid_count = 0;
    INTEGER source;
    INTEGER root = -1;
#ifdef LONGINT
    const uintmax_t integer_max = (uintmax_t)LONG_MAX;
#else
    const uintmax_t integer_max = (uintmax_t)INT_MAX;
#endif

    if (result == NULL || body_table == NULL || body_count <= 0
        || leaf_capacity <= 0) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "kdtree-shear-sphere-2balls-omp: invalid tree dimensions");
        return FAILURE;
    }
    *result = NULL;
    for (source = 0; source < body_count; source++)
        if (kd_shear_valid(cmd, nthBody(body_table, source), pivot_role))
            valid_count++;
    if (valid_count <= 0) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "kdtree-shear-sphere-2balls-omp: mask/smoothing selected no bodies");
        return FAILURE;
    }
    if ((uintmax_t)valid_count > integer_max/2
        || (uintmax_t)valid_count > (uintmax_t)SIZE_MAX/(2*sizeof(fcfc_ballnode))
        || (uintmax_t)valid_count > (uintmax_t)SIZE_MAX/sizeof(bodyptr)) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "kdtree-shear-sphere-2balls-omp: tree size overflow");
        return FAILURE;
    }

    tree = calloc(1, sizeof(*tree));
    if (tree == NULL) goto allocation_failure;
    tree->npoint = valid_count;
    tree->capacity = 2*valid_count;
    tree->bptr = malloc((size_t)valid_count*sizeof(*tree->bptr));
    tree->nodes = calloc((size_t)tree->capacity, sizeof(*tree->nodes));
    if (tree->bptr == NULL || tree->nodes == NULL)
        goto allocation_failure;

    valid_count = 0;
    for (source = 0; source < body_count; source++) {
        bodyptr body = nthBody(body_table, source);
        if (kd_shear_valid(cmd, body, pivot_role))
            tree->bptr[valid_count++] = body;
    }
    if (kd_shear_build_node(cmd, tree, 0, valid_count - 1, leaf_capacity,
                            0, pivot_role, &root) == FAILURE || root != 0) {
        kdtree_shear_sphere_free(tree);
        return FAILURE;
    }
    gd->bytes_tot += sizeof(*tree)
        + (size_t)valid_count*sizeof(*tree->bptr)
        + (size_t)tree->capacity*sizeof(*tree->nodes);
    *result = tree;
    return SUCCESS;

allocation_failure:
    snprintf(cmd->error_message, _ERRORMSGSIZE_,
             "kdtree-shear-sphere-2balls-omp: memory allocation failed");
    kdtree_shear_sphere_free(tree);
    return FAILURE;
}

int kdtree_shear_sphere_frontier(
        struct cmdline_data *cmd, fcfc_balltreeptr tree, INTEGER target,
        INTEGER **result, INTEGER *result_count)
{
    INTEGER *frontier;
    INTEGER count = 1;

    *result = NULL;
    *result_count = 0;
    if (tree == NULL || tree->nnode <= 0 || target <= 0
        || (size_t)target > SIZE_MAX/sizeof(*frontier)) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "kdtree-shear-sphere-2balls-omp: invalid frontier request");
        return FAILURE;
    }
    frontier = malloc((size_t)target*sizeof(*frontier));
    if (frontier == NULL) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "kdtree-shear-sphere-2balls-omp: frontier allocation failed");
        return FAILURE;
    }
    frontier[0] = 0;
    while (count < target) {
        INTEGER selected = -1;
        INTEGER largest = -1;
        INTEGER index;

        for (index = 0; index < count; index++) {
            fcfc_ballnode *node = &tree->nodes[frontier[index]];
            INTEGER size = node->last - node->first + 1;

            if (node->left >= 0 && node->right >= 0 && size > largest) {
                selected = index;
                largest = size;
            }
        }
        if (selected < 0) break;
        {
            fcfc_ballnode *node = &tree->nodes[frontier[selected]];
            frontier[selected] = node->left;
            frontier[count++] = node->right;
        }
    }
    *result = frontier;
    *result_count = count;
    return SUCCESS;
}
