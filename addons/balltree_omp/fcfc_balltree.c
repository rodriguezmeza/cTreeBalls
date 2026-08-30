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

static void aggregate_node(struct cmdline_data *cmd, fcfc_ballnode *node,
                           bodyptr *points)
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
        const real field_weight = cballs_opt_weights_norm(cmd)
            ? Weight(points[i]) : 0.0;
#ifdef KappaAvgON
        const real field = KappaAvg(points[i]);
#else
        const real field = Kappa(points[i]);
#endif
        DO_COORD(k) {
            cmpos_sum[k] += mass * (real)Pos(points[i])[k];
            geometric_center[k] += (real)Pos(points[i])[k];
        }
        node->weight += mass;
        kappa_sum += field;
        node->kappa_sq_sum += field * field;
        node->field_weight_sum += field_weight;
        node->field_weight_sq_sum += field_weight * field_weight;
        node->weighted_kappa_sum += field_weight * field;
        node->weighted_kappa_sq_sum +=
            (field_weight * field) * (field_weight * field);
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
                      INTEGER *result)
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
    aggregate_node(cmd, node, tree->bptr);

    if (hi - lo + 1 > nleaf) {
        const INTEGER median = lo + (hi - lo + 1) / 2;
        select_median(tree->bptr, lo, hi, median, axes[0]);
        if (build_node(cmd, tree, lo, median - 1, nleaf, depth + 1,
                       &node->left) == FAILURE ||
            build_node(cmd, tree, median, hi, nleaf, depth + 1,
                       &node->right) == FAILURE)
            return FAILURE;
    }

    *result = index;
    return SUCCESS;
}

int fcfc_balltree_build(struct cmdline_data *cmd, struct global_data *gd,
                        bodyptr btab, INTEGER nbody, int nleaf,
                        fcfc_balltreeptr *result)
{
    fcfc_balltreeptr tree = NULL;
    INTEGER root = -1;
    INTEGER i;

    if (result == NULL || btab == NULL || nbody <= 0 || nleaf <= 0) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "fcfc_balltree_build: invalid tree dimensions");
        return FAILURE;
    }
    *result = NULL;
    if ((uintmax_t)nbody >
#ifdef LONGINT
        (uintmax_t)LONG_MAX / 2 ||
#else
        (uintmax_t)INT_MAX / 2 ||
#endif
        (uintmax_t)nbody > (uintmax_t)SIZE_MAX / (2 * sizeof(fcfc_ballnode)) ||
        (uintmax_t)nbody > (uintmax_t)SIZE_MAX / sizeof(bodyptr)
#ifdef SINGLEP
        || (uintmax_t)nbody >
            (uintmax_t)SIZE_MAX / sizeof(fcfc_ballpoint)
#endif
        ) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "fcfc_balltree_build: tree size overflows size_t");
        return FAILURE;
    }

    tree = calloc(1, sizeof(*tree));
    if (tree == NULL) goto allocation_failure;
    tree->capacity = 2 * nbody;
    tree->npoint = nbody;
    tree->bptr = malloc((size_t)nbody * sizeof(*tree->bptr));
    tree->nodes = calloc((size_t)tree->capacity, sizeof(*tree->nodes));
    if (tree->bptr == NULL || tree->nodes == NULL) goto allocation_failure;

    for (i = 0; i < nbody; i++)
        tree->bptr[i] = nthBody(btab, i);

    if (build_node(cmd, tree, 0, nbody - 1, nleaf, 0, &root) == FAILURE ||
        root != FCFC_BALLTREE_ROOT) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "fcfc_balltree_build: node construction failed");
        fcfc_balltree_free(tree);
        return FAILURE;
    }

#ifdef SINGLEP
    tree->packed_points = malloc(
        (size_t)nbody * sizeof(*tree->packed_points));
    if (tree->packed_points == NULL) goto allocation_failure;
    for (i = 0; i < nbody; i++) {
        SETV(tree->packed_points[i].pos, Pos(tree->bptr[i]));
        tree->packed_points[i].kappa = Kappa(tree->bptr[i]);
    }
#endif

    gd->bytes_tot += sizeof(*tree) +
        (size_t)nbody * sizeof(*tree->bptr) +
        (size_t)tree->capacity * sizeof(*tree->nodes)
#ifdef SINGLEP
        + (size_t)nbody * sizeof(*tree->packed_points)
#endif
        ;
    *result = tree;
    return SUCCESS;

allocation_failure:
    snprintf(cmd->error_message, _ERRORMSGSIZE_,
             "fcfc_balltree_build: memory allocation failed");
    fcfc_balltree_free(tree);
    return FAILURE;
}

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
#ifdef SINGLEP
    free(tree->packed_points);
#endif
    free(tree->nodes);
    free(tree->bptr);
    free(tree);
}
