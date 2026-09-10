/*
 * FCFC-style ball tree for cTreeBalls.
 *
 * The PCA split and bounding-sphere construction are adapted from FCFC:
 * https://github.com/cheng-zhao/FCFC
 * Copyright (c) 2020--2022 Cheng Zhao.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include <float.h>
#include <limits.h>
#include <stdint.h>
#include <stdlib.h>

#include "globaldefs.h"
#include "fcfc_balltree.h"

#define FCFC_BALLTREE_ROOT 0

static real fcfc_balltree_field(bodyptr point)
{
#ifdef KappaAvgON
    return KappaAvg(point);
#else
    return Kappa(point);
#endif
}

static int fcfc_balltree_pack_points(const struct global_data *gd)
{
#ifdef SINGLEP
    (void)gd;
    return TRUE;
#else
#ifdef BALLTREE2BALLSOMP
    if (gd->searchMethod_int == BALLTREE2BALLSMETHOD) return TRUE;
#endif
#ifdef BALLTREE2BALLSMPI
    if (gd->searchMethod_int == BALLTREE2BALLSMPIMETHOD) return TRUE;
#endif
    return FALSE;
#endif
}

static real projection(bodyptr p, const real axis[NDIM])
{
    real value = 0.0;
    int k;

    DO_COORD(k)
        value += Pos(p)[k] * axis[k];
    return value;
}

static int point_before(bodyptr a, bodyptr b, const real axis[NDIM])
{
    const real pa = projection(a, axis);
    const real pb = projection(b, axis);

    if (pa < pb) return TRUE;
    if (pa > pb) return FALSE;
    return (uintptr_t)a < (uintptr_t)b;
}

static void swap_body(bodyptr *a, bodyptr *b)
{
    bodyptr tmp = *a;
    *a = *b;
    *b = tmp;
}

/* Partition around the requested median without allocating projection arrays. */
static void select_median(bodyptr *points, INTEGER lo, INTEGER hi,
                          INTEGER target, const real axis[NDIM])
{
    while (lo < hi) {
        INTEGER i = lo;
        INTEGER j = hi;
        bodyptr pivot = points[lo + (hi - lo) / 2];

        while (i <= j) {
            while (point_before(points[i], pivot, axis)) i++;
            while (point_before(pivot, points[j], axis)) j--;
            if (i <= j) {
                swap_body(&points[i], &points[j]);
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

/* Jacobi diagonalisation of the small symmetric covariance matrix. */
static void principal_axes(bodyptr *points, INTEGER lo, INTEGER hi,
                           real axes[NDIM][NDIM])
{
    double mean[NDIM] = {0};
    double cov[NDIM][NDIM] = {{0}};
    double eigvec[NDIM][NDIM] = {{0}};
    double eigval[NDIM];
    const double count = (double)(hi - lo + 1);
    INTEGER i;
    int j, k, iteration;

    for (i = lo; i <= hi; i++)
        DO_COORD(j)
            mean[j] += Pos(points[i])[j];
    DO_COORD(j)
        mean[j] /= count;

    for (i = lo; i <= hi; i++) {
        DO_COORD(j) {
            const double dj = Pos(points[i])[j] - mean[j];
            DO_COORD(k)
                cov[j][k] += dj * (Pos(points[i])[k] - mean[k]);
        }
    }
    DO_COORD(j)
        eigvec[j][j] = 1.0;

    for (iteration = 0; iteration < 32; iteration++) {
        int p = 0;
        int q = NDIM > 1 ? 1 : 0;
        double largest = 0.0;

        DO_COORD(j) {
            for (k = j + 1; k < NDIM; k++) {
                const double value = fabs(cov[j][k]);
                if (value > largest) {
                    largest = value;
                    p = j;
                    q = k;
                }
            }
        }
        if (largest <= DBL_EPSILON) break;

        const double angle = 0.5 * atan2(2.0 * cov[p][q],
                                         cov[q][q] - cov[p][p]);
        const double c = cos(angle);
        const double s = sin(angle);
        const double app = cov[p][p];
        const double aqq = cov[q][q];
        const double apq = cov[p][q];

        cov[p][p] = c*c*app - 2.0*s*c*apq + s*s*aqq;
        cov[q][q] = s*s*app + 2.0*s*c*apq + c*c*aqq;
        cov[p][q] = cov[q][p] = 0.0;
        DO_COORD(j) {
            if (j != p && j != q) {
                const double ajp = cov[j][p];
                const double ajq = cov[j][q];
                cov[j][p] = cov[p][j] = c*ajp - s*ajq;
                cov[j][q] = cov[q][j] = s*ajp + c*ajq;
            }
            const double vjp = eigvec[j][p];
            const double vjq = eigvec[j][q];
            eigvec[j][p] = c*vjp - s*vjq;
            eigvec[j][q] = s*vjp + c*vjq;
        }
    }

    DO_COORD(j)
        eigval[j] = cov[j][j];
    DO_COORD(j) {
        int best = j;
        for (k = j + 1; k < NDIM; k++)
            if (eigval[k] > eigval[best]) best = k;
        if (best != j) {
            const double tmpval = eigval[j];
            eigval[j] = eigval[best];
            eigval[best] = tmpval;
            DO_COORD(k) {
                const double tmpvec = eigvec[k][j];
                eigvec[k][j] = eigvec[k][best];
                eigvec[k][best] = tmpvec;
            }
        }
        DO_COORD(k)
            axes[j][k] = (real)eigvec[k][j];
    }
}

static real storage_distance_squared(const cballs_storage_real a[NDIM],
                                     const cballs_storage_real b[NDIM])
{
    real value = 0.0;
    int k;

    DO_COORD(k)
        value += rsqr((real)a[k] - (real)b[k]);
    return value;
}

static real mixed_distance_squared(const compute_vector a,
                                   const cballs_storage_real b[NDIM])
{
    real value = 0.0;
    int k;

    DO_COORD(k)
        value += rsqr(a[k] - (real)b[k]);
    return value;
}

/* FCFC seeds the sphere from extremes along the first two PCA directions. */
static void enclosing_sphere(bodyptr *points, INTEGER lo, INTEGER hi,
                             real axes[NDIM][NDIM], compute_vector center,
                             real *radius)
{
    INTEGER extreme[4];
    int nextreme = 0;
    int ndirection = NDIM > 1 ? 2 : 1;
    int direction;

    for (direction = 0; direction < ndirection; direction++) {
        INTEGER imin = lo;
        INTEGER imax = lo;
        real pmin = projection(points[lo], axes[direction]);
        real pmax = pmin;
        INTEGER i;

        for (i = lo + 1; i <= hi; i++) {
            const real value = projection(points[i], axes[direction]);
            if (value < pmin) { pmin = value; imin = i; }
            if (value > pmax) { pmax = value; imax = i; }
        }
        extreme[nextreme++] = imin;
        extreme[nextreme++] = imax;
    }

    INTEGER ia = extreme[0];
    INTEGER ib = extreme[0];
    real farthest = -1.0;
    int a, b, k;
    for (a = 0; a < nextreme; a++) {
        for (b = a + 1; b < nextreme; b++) {
            const real d2 = storage_distance_squared(
                Pos(points[extreme[a]]), Pos(points[extreme[b]]));
            if (d2 > farthest) {
                farthest = d2;
                ia = extreme[a];
                ib = extreme[b];
            }
        }
    }

    DO_COORD(k)
        center[k] = 0.5 * ((real)Pos(points[ia])[k]
                         + (real)Pos(points[ib])[k]);
    *radius = farthest > 0.0 ? 0.5 * rsqrt(farthest) : 0.0;

    /* Ritter growth makes the PCA seed an enclosing sphere. */
    INTEGER i;
    for (i = lo; i <= hi; i++) {
        const real d2 = mixed_distance_squared(center, Pos(points[i]));
        if (d2 > rsqr(*radius)) {
            const real distance = rsqrt(d2);
            const real grown = 0.5 * (*radius + distance);
            const real shift = distance > 0.0
                ? (distance - grown) / distance : 0.0;
            DO_COORD(k)
                center[k] += ((real)Pos(points[i])[k] - center[k]) * shift;
            *radius = grown;
        }
    }

}

static bool fcfc_scalar_body_valid(const struct cmdline_data *cmd,
                                   bodyptr body, bool filter_active,
                                   bool pivot_role)
{
    if (!filter_active) return TRUE;
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

static void fcfc_scalar_body_values(const struct cmdline_data *cmd,
                                    bodyptr body, bool pivot_role,
                                    real *field, real *normalization,
                                    real *raw_field)
{
#ifdef SMOOTHPIVOT
    if (pivot_role && cballs_opt_smooth_pivot(cmd)) {
        *normalization = cballs_opt_weights_norm(cmd)
            ? WeightRmin(body) : (real)MAX(NbRmin(body), 1);
        *field = KappaRmin(body);
        *raw_field = *normalization != 0.0
            ? *field / *normalization : 0.0;
        return;
    }
#else
    (void)pivot_role;
#endif
    *raw_field = fcfc_balltree_field(body);
    *normalization = cballs_opt_weights_norm(cmd) ? Weight(body) : 1.0;
    *field = *normalization * *raw_field;
}

static void aggregate_node(struct cmdline_data *cmd, fcfc_ballnode *node,
                           bodyptr *points, bool pivot_role)
{
    compute_vector cmpos_sum;
    compute_vector geometric_center;
    real kappa_sum = 0.0;
    INTEGER i;
    int k;

    CLRV(cmpos_sum);
    CLRV(geometric_center);
    node->weight = 0.0;
    node->kappa_sq_sum = 0.0;
    node->field_weight_sum = 0.0;
    node->field_weight_sq_sum = 0.0;
    node->weighted_kappa_sum = 0.0;
    node->weighted_kappa_sq_sum = 0.0;
    for (i = node->first; i <= node->last; i++) {
        const real mass = Mass(points[i]);
        real field_weight;
        real weighted_field;
        real field;

        fcfc_scalar_body_values(cmd, points[i], pivot_role,
                                &weighted_field, &field_weight, &field);
        DO_COORD(k) {
            cmpos_sum[k] += mass * (real)Pos(points[i])[k];
            geometric_center[k] += (real)Pos(points[i])[k];
        }
        node->weight += mass;
        kappa_sum += field;
        node->kappa_sq_sum += field * field;
        node->field_weight_sum += field_weight;
        node->field_weight_sq_sum += field_weight * field_weight;
        node->weighted_kappa_sum += weighted_field;
        node->weighted_kappa_sq_sum += weighted_field * weighted_field;
    }
    if (node->weight > 0.0) {
        DO_COORD(k)
            node->cmpos[k] = (cballs_storage_real)
                (cmpos_sum[k] / node->weight);
    } else {
        const real count = (real)(node->last - node->first + 1);
        DO_COORD(k)
            node->cmpos[k] = (cballs_storage_real)
                (geometric_center[k] / count);
    }
    node->kappa_sum = kappa_sum;
    node->kappa = kappa_sum / (real)(node->last - node->first + 1);

    real max2 = 0.0;
    for (i = node->first; i <= node->last; i++) {
        const real d2 = storage_distance_squared(node->cmpos,
                                                  Pos(points[i]));
        if (d2 > max2) max2 = d2;
    }
    if (cmd->theta == 0.0) {
        node->aggregate_radius = cballs_store_upper_bound(MAX_REAL_NUMBER);
    } else {
        const real radius = rsqrt(max2) / cmd->theta;
        node->aggregate_radius = cballs_store_search_bound(radius);
    }
}

static int build_node(struct cmdline_data *cmd, fcfc_balltreeptr tree,
                      INTEGER lo, INTEGER hi, int nleaf, int depth,
                      bool pivot_role, INTEGER *result)
{
    if (tree->nnode >= tree->capacity) return FAILURE;

    const INTEGER index = tree->nnode++;
    fcfc_ballnode *node = &tree->nodes[index];
    real axes[NDIM][NDIM];
    compute_vector center;
    real radius;
    INTEGER i;
    int k;

    node->first = lo;
    node->last = hi;
    node->left = node->right = -1;
    if (depth > tree->max_depth) tree->max_depth = depth;

    principal_axes(tree->bptr, lo, hi, axes);
    enclosing_sphere(tree->bptr, lo, hi, axes, center, &radius);
    DO_COORD(k)
        node->center[k] = (cballs_storage_real)center[k];

    /* Recompute after the storage cast so the stored sphere still encloses. */
    radius = 0.0;
    for (i = lo; i <= hi; i++) {
        const real d2 = storage_distance_squared(node->center,
                                                  Pos(tree->bptr[i]));
        if (d2 > radius) radius = d2;
    }
    radius = rsqrt(radius);
    node->radius = cballs_store_search_bound(radius);
    aggregate_node(cmd, node, tree->bptr, pivot_role);

    if (hi - lo + 1 > nleaf) {
        const INTEGER median = lo + (hi - lo + 1) / 2;
        select_median(tree->bptr, lo, hi, median, axes[0]);
        if (build_node(cmd, tree, lo, median - 1, nleaf, depth + 1,
                       pivot_role,
                       &node->left) == FAILURE ||
            build_node(cmd, tree, median, hi, nleaf, depth + 1,
                       pivot_role,
                       &node->right) == FAILURE)
            return FAILURE;
    }

    *result = index;
    return SUCCESS;
}

static int fcfc_balltree_build_internal(
        struct cmdline_data *cmd, struct global_data *gd,
        bodyptr btab, INTEGER nbody, int nleaf,
        bool filter_active, bool pivot_role, fcfc_balltreeptr *result)
{
    fcfc_balltreeptr tree = NULL;
    INTEGER valid_count = 0;
    INTEGER root = -1;
    INTEGER i;

    if (result == NULL || btab == NULL || nbody <= 0 || nleaf <= 0) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "fcfc_balltree_build: invalid tree dimensions");
        return FAILURE;
    }
    *result = NULL;
    for (i = 0; i < nbody; i++)
        if (fcfc_scalar_body_valid(cmd, nthBody(btab, i),
                                   filter_active, pivot_role))
            valid_count++;
    if (valid_count <= 0) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "fcfc_balltree_build: mask/smoothing selected no bodies");
        return FAILURE;
    }
    if ((uintmax_t)valid_count >
#ifdef LONGINT
        (uintmax_t)LONG_MAX / 2 ||
#else
        (uintmax_t)INT_MAX / 2 ||
#endif
        (uintmax_t)valid_count
            > (uintmax_t)SIZE_MAX / (2 * sizeof(fcfc_ballnode)) ||
        (uintmax_t)valid_count > (uintmax_t)SIZE_MAX / sizeof(bodyptr) ||
        (uintmax_t)valid_count > (uintmax_t)SIZE_MAX / sizeof(fcfc_ballpoint)
        ) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "fcfc_balltree_build: tree size overflows size_t");
        return FAILURE;
    }

    tree = calloc(1, sizeof(*tree));
    if (tree == NULL) goto allocation_failure;
    tree->capacity = 2 * valid_count;
    tree->npoint = valid_count;
    tree->bptr = malloc((size_t)valid_count * sizeof(*tree->bptr));
    tree->nodes = calloc((size_t)tree->capacity, sizeof(*tree->nodes));
    if (tree->bptr == NULL || tree->nodes == NULL) goto allocation_failure;

    valid_count = 0;
    for (i = 0; i < nbody; i++) {
        bodyptr body = nthBody(btab, i);

        if (fcfc_scalar_body_valid(cmd, body, filter_active, pivot_role))
            tree->bptr[valid_count++] = body;
    }

    if (build_node(cmd, tree, 0, valid_count - 1, nleaf, 0,
                   pivot_role, &root) == FAILURE ||
        root != FCFC_BALLTREE_ROOT) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "fcfc_balltree_build: node construction failed");
        fcfc_balltree_free(tree);
        return FAILURE;
    }

    if (filter_active || fcfc_balltree_pack_points(gd)) {
        tree->packed_points = malloc(
            (size_t)valid_count * sizeof(*tree->packed_points));
        if (tree->packed_points == NULL) goto allocation_failure;
        for (i = 0; i < valid_count; i++) {
            real field;
            real weight;
            real raw_field;

            fcfc_scalar_body_values(cmd, tree->bptr[i], pivot_role,
                                    &field, &weight, &raw_field);
            SETV(tree->packed_points[i].pos, Pos(tree->bptr[i]));
            tree->packed_points[i].kappa = raw_field;
            tree->packed_points[i].weight = weight;
            tree->packed_points[i].weighted_kappa = field;
            tree->packed_points[i].source = tree->bptr[i];
        }
    }

    gd->bytes_tot += sizeof(*tree) +
        (size_t)valid_count * sizeof(*tree->bptr) +
        (size_t)tree->capacity * sizeof(*tree->nodes) +
        (tree->packed_points == NULL ? 0 :
         (size_t)valid_count * sizeof(*tree->packed_points));
    *result = tree;
    return SUCCESS;

allocation_failure:
    snprintf(cmd->error_message, _ERRORMSGSIZE_,
             "fcfc_balltree_build: memory allocation failed");
    fcfc_balltree_free(tree);
    return FAILURE;
}

int fcfc_balltree_build(struct cmdline_data *cmd, struct global_data *gd,
                        bodyptr btab, INTEGER nbody, int nleaf,
                        fcfc_balltreeptr *result)
{
    return fcfc_balltree_build_internal(
        cmd, gd, btab, nbody, nleaf, FALSE, FALSE, result);
}

int fcfc_balltree_build_scalar_role(
        struct cmdline_data *cmd, struct global_data *gd,
        bodyptr btab, INTEGER nbody, int nleaf,
        bool pivot_role, fcfc_balltreeptr *result)
{
    return fcfc_balltree_build_internal(
        cmd, gd, btab, nbody, nleaf, TRUE, pivot_role, result);
}

#ifdef BALLTREESHEARSPHERE2BALLSOMP

typedef struct {
    real re;
    real im;
} fcfc_shear_complex;

static fcfc_shear_complex fcfc_shear_make(real re, real im)
{
    fcfc_shear_complex value = {re, im};
    return value;
}

static fcfc_shear_complex fcfc_shear_mul(fcfc_shear_complex a,
                                         fcfc_shear_complex b)
{
    return fcfc_shear_make(a.re*b.re - a.im*b.im,
                           a.re*b.im + a.im*b.re);
}

static fcfc_shear_complex fcfc_shear_scale(fcfc_shear_complex value,
                                           real scale)
{
    return fcfc_shear_make(value.re*scale, value.im*scale);
}

static real fcfc_shear_dot3(const real *a, const real *b)
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

static bool fcfc_shear_unit3(const real *position, compute_vector unit)
{
    const real norm2 = fcfc_shear_dot3(position, position);
    real inverse_norm;

    if (!(norm2 > 0.0) || !isfinite(norm2)) return FALSE;
    inverse_norm = 1.0/rsqrt(norm2);
    unit[0] = position[0]*inverse_norm;
    unit[1] = position[1]*inverse_norm;
    unit[2] = position[2]*inverse_norm;
    return isfinite(unit[0]) && isfinite(unit[1]) && isfinite(unit[2]);
}

static bool fcfc_shear_basis(const real *unit, compute_vector east,
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

static bool fcfc_shear_rotation(const real *target_unit,
                                const real *target_east,
                                const real *target_north,
                                const real *source_position,
                                fcfc_shear_complex *rotation)
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

    if (!fcfc_shear_unit3(source_position, source_unit)
        || !fcfc_shear_basis(source_unit, source_east, source_north))
        return FALSE;
    (void)source_north;
    denominator = 1.0 + fcfc_shear_dot3(target_unit, source_unit);
    if (!(denominator > 64.0*DBL_EPSILON) || !isfinite(denominator))
        return FALSE;
    for (axis = 0; axis < 3; axis++) {
        transported_east[axis] = target_east[axis]
            - fcfc_shear_dot3(target_east, source_unit)/denominator
              *(target_unit[axis] + source_unit[axis]);
        transported_north[axis] = target_north[axis]
            - fcfc_shear_dot3(target_north, source_unit)/denominator
              *(target_unit[axis] + source_unit[axis]);
    }
    c = fcfc_shear_dot3(source_east, transported_east);
    s = fcfc_shear_dot3(source_east, transported_north);
    norm = rsqrt(c*c + s*s);
    if (!(norm > 64.0*DBL_EPSILON) || !isfinite(norm)) return FALSE;
    c /= norm;
    s /= norm;
    *rotation = fcfc_shear_make(c*c - s*s, 2.0*c*s);
    return isfinite(rotation->re) && isfinite(rotation->im);
}

static bool fcfc_shear_valid(const struct cmdline_data *cmd, bodyptr body,
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

static real fcfc_shear_body_weight(const struct cmdline_data *cmd,
                                   bodyptr body, bool pivot_role)
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

static fcfc_shear_complex fcfc_shear_body_gamma(
        const struct cmdline_data *cmd, bodyptr body, bool pivot_role)
{
#ifdef SMOOTHPIVOT
    if (pivot_role && cballs_opt_smooth_pivot(cmd))
        return fcfc_shear_make(Gamma1Rmin(body), Gamma2Rmin(body));
#else
    (void)cmd;
    (void)pivot_role;
#endif
    return fcfc_shear_scale(
        fcfc_shear_make(Gamma1(body), Gamma2(body)), Weight(body));
}

static int fcfc_shear_aggregate(struct cmdline_data *cmd,
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
        const real weight = fcfc_shear_body_weight(cmd, body, pivot_role);
        const real position_weight = weight > 0.0 ? weight : 1.0;

        DO_COORD(axis)
            center_sum[axis] += position_weight*(real)Pos(body)[axis];
        center_weight += position_weight;
    }
    if (center_weight > 0.0)
        DO_COORD(axis)
            center_sum[axis] /= center_weight;
    if (!fcfc_shear_unit3(center_sum, center_unit)
        && !fcfc_shear_unit3(Pos(tree->bptr[node->first]), center_unit)) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "balltree-shear-sphere-2balls-omp: invalid node center");
        return FAILURE;
    }
    if (!fcfc_shear_basis(center_unit, center_east, center_north)) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "balltree-shear-sphere-2balls-omp: node tangent frame is undefined");
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
        const real weight = fcfc_shear_body_weight(cmd, body, pivot_role);
        fcfc_shear_complex weighted_gamma =
            fcfc_shear_body_gamma(cmd, body, pivot_role);
        fcfc_shear_complex rotation;
        fcfc_shear_complex transported;
        real distance2 = 0.0;

        if (!(weight >= 0.0) || !isfinite(weight)
            || !fcfc_shear_rotation(center_unit, center_east, center_north,
                                    Pos(body), &rotation)) {
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     "balltree-shear-sphere-2balls-omp: invalid spin-2 node member");
            return FAILURE;
        }
        transported = fcfc_shear_mul(weighted_gamma, rotation);
        node->weight += weight;
        node->shear_gamma_re += transported.re;
        node->shear_gamma_im += transported.im;
        transported = fcfc_shear_mul(transported, transported);
        node->shear_gamma2_re += transported.re;
        node->shear_gamma2_im += transported.im;
        node->shear_gamma_abs2 += weighted_gamma.re*weighted_gamma.re
                                + weighted_gamma.im*weighted_gamma.im;
        node->shear_weight2 += weight*weight;
        DO_COORD(axis)
            distance2 += rsqr(center_unit[axis] - (real)Pos(body)[axis]);
        radius2 = MAX(radius2, distance2);
    }
    node->radius = cballs_store_search_bound(rsqrt(radius2));
    node->aggregate_radius = node->radius;
    return SUCCESS;
}

static int fcfc_shear_build_node(struct cmdline_data *cmd,
                                 fcfc_balltreeptr tree, INTEGER first,
                                 INTEGER last, int leaf_capacity, int depth,
                                 bool pivot_role, INTEGER *result)
{
    real axes[NDIM][NDIM];
    fcfc_ballnode *node;
    INTEGER index;

    if (tree->nnode >= tree->capacity) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "balltree-shear-sphere-2balls-omp: node capacity exceeded");
        return FAILURE;
    }
    index = tree->nnode++;
    node = &tree->nodes[index];
    node->first = first;
    node->last = last;
    node->left = -1;
    node->right = -1;
    if (depth > tree->max_depth) tree->max_depth = depth;

    principal_axes(tree->bptr, first, last, axes);
    if (fcfc_shear_aggregate(cmd, tree, node, pivot_role) == FAILURE)
        return FAILURE;
    if (last - first + 1 > leaf_capacity) {
        const INTEGER middle = first + (last - first + 1)/2;

        select_median(tree->bptr, first, last, middle, axes[0]);
        if (fcfc_shear_build_node(cmd, tree, first, middle - 1,
                                  leaf_capacity, depth + 1, pivot_role,
                                  &node->left) == FAILURE
            || fcfc_shear_build_node(cmd, tree, middle, last,
                                     leaf_capacity, depth + 1, pivot_role,
                                     &node->right) == FAILURE)
            return FAILURE;
    }
    *result = index;
    return SUCCESS;
}

int fcfc_balltree_build_shear_sphere(
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
                 "balltree-shear-sphere-2balls-omp: invalid tree dimensions");
        return FAILURE;
    }
    *result = NULL;
    for (source = 0; source < body_count; source++)
        if (fcfc_shear_valid(cmd, nthBody(body_table, source), pivot_role))
            valid_count++;
    if (valid_count <= 0) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "balltree-shear-sphere-2balls-omp: mask/smoothing selected no bodies");
        return FAILURE;
    }
    if ((uintmax_t)valid_count > integer_max/2
        || (uintmax_t)valid_count
             > (uintmax_t)SIZE_MAX/(2*sizeof(fcfc_ballnode))
        || (uintmax_t)valid_count > (uintmax_t)SIZE_MAX/sizeof(bodyptr)) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "balltree-shear-sphere-2balls-omp: tree size overflow");
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
        if (fcfc_shear_valid(cmd, body, pivot_role))
            tree->bptr[valid_count++] = body;
    }
    if (fcfc_shear_build_node(cmd, tree, 0, valid_count - 1,
                              leaf_capacity, 0, pivot_role, &root) == FAILURE
        || root != FCFC_BALLTREE_ROOT) {
        fcfc_balltree_free(tree);
        return FAILURE;
    }
    gd->bytes_tot += sizeof(*tree)
        + (size_t)valid_count*sizeof(*tree->bptr)
        + (size_t)tree->capacity*sizeof(*tree->nodes);
    *result = tree;
    return SUCCESS;

allocation_failure:
    snprintf(cmd->error_message, _ERRORMSGSIZE_,
             "balltree-shear-sphere-2balls-omp: memory allocation failed");
    fcfc_balltree_free(tree);
    return FAILURE;
}

#endif /* BALLTREESHEARSPHERE2BALLSOMP */

static int minimum_leaf_depth(const fcfc_balltreeptr tree, INTEGER inode,
                              int depth)
{
    const fcfc_ballnode *node = &tree->nodes[inode];

    if (node->left < 0) return depth;
    const int left = minimum_leaf_depth(tree, node->left, depth + 1);
    const int right = minimum_leaf_depth(tree, node->right, depth + 1);
    return left < right ? left : right;
}

static void gather_frontier(const fcfc_balltreeptr tree, INTEGER inode,
                            int depth, int requested_depth,
                            INTEGER *frontier, INTEGER *count)
{
    const fcfc_ballnode *node = &tree->nodes[inode];

    if (depth == requested_depth || node->left < 0) {
        frontier[(*count)++] = inode;
        return;
    }
    gather_frontier(tree, node->left, depth + 1, requested_depth,
                    frontier, count);
    gather_frontier(tree, node->right, depth + 1, requested_depth,
                    frontier, count);
}

/*
 * Return a complete breadth-first frontier.  Limiting the requested depth to
 * the shallowest leaf ensures that every point belongs to exactly one task.
 */
int fcfc_balltree_frontier(struct cmdline_data *cmd, fcfc_balltreeptr tree,
                           INTEGER minimum_count, INTEGER **result,
                           INTEGER *result_count)
{
    INTEGER capacity = 1;
    INTEGER count = 0;
    int depth = 0;
    int max_frontier_depth;
    INTEGER *frontier;

    if (tree == NULL || tree->nodes == NULL || tree->nnode < 1
        || minimum_count < 1 || result == NULL || result_count == NULL) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "fcfc_balltree_frontier: invalid arguments");
        return FAILURE;
    }
    *result = NULL;
    *result_count = 0;

    max_frontier_depth = minimum_leaf_depth(tree, FCFC_BALLTREE_ROOT, 0);
    while (capacity < minimum_count && depth < max_frontier_depth) {
        if (capacity >
#ifdef LONGINT
            LONG_MAX / 2
#else
            INT_MAX / 2
#endif
            ) {
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     "fcfc_balltree_frontier: task count overflow");
            return FAILURE;
        }
        capacity *= 2;
        depth++;
    }
    if ((size_t)capacity > SIZE_MAX / sizeof(*frontier)) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "fcfc_balltree_frontier: allocation size overflow");
        return FAILURE;
    }
    frontier = malloc((size_t)capacity * sizeof(*frontier));
    if (frontier == NULL) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "fcfc_balltree_frontier: memory allocation failed");
        return FAILURE;
    }

    gather_frontier(tree, FCFC_BALLTREE_ROOT, 0, depth, frontier, &count);
    if (count < 1) {
        free(frontier);
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "fcfc_balltree_frontier: empty task frontier");
        return FAILURE;
    }
    *result = frontier;
    *result_count = count;
    return SUCCESS;
}

void fcfc_balltree_free(fcfc_balltreeptr tree)
{
    if (tree == NULL) return;
    free(tree->packed_points);
    free(tree->nodes);
    free(tree->bptr);
    free(tree);
}
