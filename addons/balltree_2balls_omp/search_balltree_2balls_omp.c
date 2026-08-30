/*
 * TreeCorr-style dual-ball 2PCF and triple-node scalar 3PCF search.
 *
 * The node-pair recursion and split heuristic are adapted from TreeCorr:
 * Copyright (c) 2003-2024 Mike Jarvis, used under its BSD-style license.
 * The ball-tree builder is the FCFC-derived implementation shared with
 * balltree-omp; its source file carries the full FCFC MIT notice.
 */

#include <float.h>
#include <limits.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "globaldefs.h"
#include "fcfc_balltree.h"

#ifndef TREECORR_METHOD_NAME
#define TREECORR_METHOD_NAME "balltree-2balls-omp"
#endif
#ifndef TREECORR_PUBLISH_NODE_COUNT
#define TREECORR_PUBLISH_NODE_COUNT(gd, catalog, count) \
    ((gd)->ncellTable[(catalog)] = (count))
#endif

#define TREECORR_SPLIT_FACTOR ((real)0.585)
#define TREECORR_PAIR_FRONTIER_TARGET ((INTEGER)256)
#define TREECORR_TRIPLE_FRONTIER_TARGET ((INTEGER)32)
#define TREECORR_LARGE_DIRECT_3PCF ((INTEGER)1000000)
#define TREECORR_TRIPLE_TASK_TARGET ((INTEGER)256)
#define TREECORR_TRIPLE_TASK_MEMORY ((size_t)128 * 1024 * 1024)

#ifdef TWOPCF
typedef struct {
    real *pair_count;
    real *weight_product;
    real *field_product;
    INTEGER body_pairs;
    INTEGER cell_pairs;
} treecorr_histogram;
#endif

typedef struct {
    struct cmdline_data *cmd;
    struct global_data *gd;
    bool use_two_balls;
    bool use_three_cells;
    bool weighted;
    real angular_tolerance;
    real max_angular_ratio;
} treecorr_search_context;

#ifdef TREECORR_DISTRIBUTED_ENGINE
static inline bool treecorr_distributed_task_owned(INTEGER task)
{
    return TREECORR_DISTRIBUTED_TASK_OWNED(task);
}

static inline bool treecorr_distributed_publish(void)
{
    return TREECORR_DISTRIBUTED_IS_ROOT();
}

static inline int treecorr_distributed_consensus(
        struct cmdline_data *cmd, int status, const char *operation)
{
    return TREECORR_DISTRIBUTED_CONSENSUS(cmd, status, operation);
}

static inline int treecorr_distributed_reduce_reals(
        struct cmdline_data *cmd, real *values, size_t count)
{
    return TREECORR_DISTRIBUTED_REDUCE_REALS(cmd, values, count);
}

static inline int treecorr_distributed_reduce_integers(
        struct cmdline_data *cmd, INTEGER *values, size_t count)
{
    return TREECORR_DISTRIBUTED_REDUCE_INTEGERS(cmd, values, count);
}
#else
static inline bool treecorr_distributed_task_owned(INTEGER task)
{
    (void)task;
    return TRUE;
}

static inline bool treecorr_distributed_publish(void)
{
    return TRUE;
}

static inline int treecorr_distributed_consensus(
        struct cmdline_data *cmd, int status, const char *operation)
{
    (void)cmd;
    (void)operation;
    return status;
}

static inline int treecorr_distributed_reduce_reals(
        struct cmdline_data *cmd, real *values, size_t count)
{
    (void)cmd;
    (void)values;
    (void)count;
    return SUCCESS;
}

static inline int treecorr_distributed_reduce_integers(
        struct cmdline_data *cmd, INTEGER *values, size_t count)
{
    (void)cmd;
    (void)values;
    (void)count;
    return SUCCESS;
}
#endif

#ifdef TWOPCF
static int treecorr_run_pair_tasks(
        const treecorr_search_context *, fcfc_balltreeptr,
        fcfc_balltreeptr, bool, const INTEGER *, INTEGER,
        const INTEGER *, INTEGER);
#endif

#ifdef THREEPCFCONVERGENCE
static inline bool treecorr_normalize_3pcf(const struct cmdline_data *cmd)
{
    return !scanopt(cmd->options, "no-normalize-HistZeta");
}
#endif

static inline real treecorr_body_field(bodyptr p)
{
#ifdef KappaAvgON
    return KappaAvg(p);
#else
    return Kappa(p);
#endif
}

static inline INTEGER treecorr_node_count(const fcfc_ballnode *ball_node)
{
    return ball_node->last - ball_node->first + 1;
}

static inline bool treecorr_node_is_leaf(const fcfc_ballnode *ball_node)
{
    return ball_node->left < 0;
}

static inline real treecorr_center_distance(struct cmdline_data *cmd,
                                            struct global_data *gd,
                                            const fcfc_ballnode *node1,
                                            const fcfc_ballnode *node2)
{
    compute_vector dr;
    real distance2;

    DOTPSUBV(distance2, dr, node1->center, node2->center);
    if (cmd->usePeriodic) {
        VWrapAll(dr);
        DOTVP(distance2, dr, dr);
    }
    return rsqrt(distance2);
}

static inline real treecorr_body_distance(struct cmdline_data *cmd,
                                          struct global_data *gd,
                                          bodyptr p, bodyptr q)
{
    compute_vector dr;
    real distance2;

    DOTPSUBV(distance2, dr, Pos(p), Pos(q));
    if (cmd->usePeriodic) {
        VWrapAll(dr);
        DOTVP(distance2, dr, dr);
    }
    return rsqrt(distance2);
}

static inline int treecorr_bin_index(const treecorr_search_context *context,
                                     real distance)
{
    const struct cmdline_data *cmd = context->cmd;
    const struct global_data *gd = context->gd;
    int n;

    if (!(distance > cmd->rminHist && distance < cmd->rangeN))
        return -1;
    if (cmd->useLogHist) {
        if (cmd->rminHist == 0.0) {
            n = (int)(cmd->logHistBinsPD
                * (rlog10(distance) - rlog10(cmd->rangeN))
                + cmd->sizeHistN) + 1;
        } else {
            n = (int)(rlog10(distance / cmd->rminHist)
                * gd->i_deltaR) + 1;
        }
    } else {
        n = (int)((distance - cmd->rminHist) * gd->i_deltaR) + 1;
    }
    return n >= 1 && n <= cmd->sizeHistN ? n : -1;
}

/* TreeCorr's b = bin_slop * bin_size, expressed in distance units. */
static inline real treecorr_bin_slop_width(
        const treecorr_search_context *context, real distance)
{
    if (context->cmd->useLogHist)
        return context->cmd->theta * rlog(10.0)
             * context->gd->deltaR * distance;
    return context->cmd->theta * context->gd->deltaR;
}

#ifdef TWOPCF
static inline void treecorr_accumulate_body_pair(
        const treecorr_search_context *context, treecorr_histogram *hist,
        bodyptr p, bodyptr q)
{
    const int n = treecorr_bin_index(
        context, treecorr_body_distance(context->cmd, context->gd, p, q));
    real denominator;
    real numerator;

    if (n < 0) return;
    if (context->weighted) {
        const real wp = Weight(p);
        const real wq = Weight(q);
        denominator = wp * wq;
        numerator = (wp * treecorr_body_field(p))
                  * (wq * treecorr_body_field(q));
    } else {
        denominator = 1.0;
        numerator = treecorr_body_field(p) * treecorr_body_field(q);
    }
    hist->pair_count[n] += 1.0;
    hist->weight_product[n] += denominator;
    hist->field_product[n] += numerator;
    hist->body_pairs++;
}

static inline bool treecorr_pair_outside_range(
        const treecorr_search_context *context,
        const fcfc_ballnode *node1, const fcfc_ballnode *node2,
        real distance)
{
    const real size = (real)node1->radius + (real)node2->radius;

    if (distance + size <= context->cmd->rminHist) return TRUE;
    if (distance - size >= context->cmd->rangeN) return TRUE;
    return FALSE;
}

/* Return the common bin, or -1 when the pair must be split. */
static int treecorr_two_ball_bin(const treecorr_search_context *context,
                                 const fcfc_ballnode *node1,
                                 const fcfc_ballnode *node2, real distance)
{
    const struct cmdline_data *cmd = context->cmd;
    const real size = (real)node1->radius + (real)node2->radius;
    const real lower = distance - size;
    const real upper = distance + size;
    int center_bin;

    if (!context->use_two_balls || !(distance > 0.0)
        || !(cmd->theta > 0.0))
        return -1;

    if (size > treecorr_bin_slop_width(context, distance)) return -1;
    if (!(lower > cmd->rminHist && upper < cmd->rangeN)) return -1;
    center_bin = treecorr_bin_index(context, distance);
    if (center_bin < 0
        || treecorr_bin_index(context, lower) != center_bin
        || treecorr_bin_index(context, upper) != center_bin)
        return -1;
    return center_bin;
}

static inline void treecorr_accumulate_node_pair(
        const treecorr_search_context *context, treecorr_histogram *hist,
        const fcfc_ballnode *node1, const fcfc_ballnode *node2, int n)
{
    const real count = (real)treecorr_node_count(node1)
                     * (real)treecorr_node_count(node2);
    real denominator;
    real numerator;

    if (context->weighted) {
        denominator = node1->field_weight_sum * node2->field_weight_sum;
        numerator = node1->weighted_kappa_sum
                  * node2->weighted_kappa_sum;
    } else {
        denominator = count;
        numerator = node1->kappa_sum * node2->kappa_sum;
    }
    hist->pair_count[n] += count;
    hist->weight_product[n] += denominator;
    hist->field_product[n] += numerator;
    hist->cell_pairs++;
}

static void treecorr_process_pair(const treecorr_search_context *context,
                                  const fcfc_balltreeptr tree1, INTEGER i1,
                                  const fcfc_balltreeptr tree2, INTEGER i2,
                                  treecorr_histogram *hist)
{
    const fcfc_ballnode *node1 = &tree1->nodes[i1];
    const fcfc_ballnode *node2 = &tree2->nodes[i2];
    const bool leaf1 = treecorr_node_is_leaf(node1);
    const bool leaf2 = treecorr_node_is_leaf(node2);
    const real distance = treecorr_center_distance(
        context->cmd, context->gd, node1, node2);
    int n;

    if (treecorr_pair_outside_range(context, node1, node2, distance)) return;

    n = treecorr_two_ball_bin(context, node1, node2, distance);
    if (n >= 0) {
        treecorr_accumulate_node_pair(context, hist, node1, node2, n);
        return;
    }

    if (leaf1 && leaf2) {
        INTEGER p1;
        INTEGER p2;
        for (p1 = node1->first; p1 <= node1->last; p1++)
            for (p2 = node2->first; p2 <= node2->last; p2++)
                treecorr_accumulate_body_pair(
                    context, hist, tree1->bptr[p1], tree2->bptr[p2]);
        return;
    }

    if (leaf1) {
        treecorr_process_pair(context, tree1, i1, tree2, node2->left, hist);
        treecorr_process_pair(context, tree1, i1, tree2, node2->right, hist);
        return;
    }
    if (leaf2) {
        treecorr_process_pair(context, tree1, node1->left, tree2, i2, hist);
        treecorr_process_pair(context, tree1, node1->right, tree2, i2, hist);
        return;
    }

    {
        bool split1 = FALSE;
        bool split2 = FALSE;
        const real size1 = (real)node1->radius;
        const real size2 = (real)node2->radius;
        real effective_width;

        if (context->cmd->useLogHist)
            effective_width = context->cmd->theta * rlog(10.0)
                * context->gd->deltaR * distance;
        else
            effective_width = context->cmd->theta * context->gd->deltaR;

        if (size2 > size1) {
            split2 = TRUE;
            if (!(size2 > 2.0 * size1))
                split1 = size1 > TREECORR_SPLIT_FACTOR * effective_width;
        } else {
            split1 = TRUE;
            if (!(size1 > 2.0 * size2))
                split2 = size2 > TREECORR_SPLIT_FACTOR * effective_width;
        }

        if (split1 && split2) {
            treecorr_process_pair(context, tree1, node1->left,
                                  tree2, node2->left, hist);
            treecorr_process_pair(context, tree1, node1->left,
                                  tree2, node2->right, hist);
            treecorr_process_pair(context, tree1, node1->right,
                                  tree2, node2->left, hist);
            treecorr_process_pair(context, tree1, node1->right,
                                  tree2, node2->right, hist);
        } else if (split1) {
            treecorr_process_pair(context, tree1, node1->left,
                                  tree2, i2, hist);
            treecorr_process_pair(context, tree1, node1->right,
                                  tree2, i2, hist);
        } else {
            treecorr_process_pair(context, tree1, i1,
                                  tree2, node2->left, hist);
            treecorr_process_pair(context, tree1, i1,
                                  tree2, node2->right, hist);
        }
    }
}

static void treecorr_process_auto(const treecorr_search_context *context,
                                  const fcfc_balltreeptr tree, INTEGER inode,
                                  treecorr_histogram *hist)
{
    const fcfc_ballnode *node = &tree->nodes[inode];

    if (2.0 * (real)node->radius <= context->cmd->rminHist) return;
    if (treecorr_node_is_leaf(node)) {
        INTEGER i;
        INTEGER j;
        for (i = node->first; i <= node->last; i++)
            for (j = i + 1; j <= node->last; j++)
                treecorr_accumulate_body_pair(
                    context, hist, tree->bptr[i], tree->bptr[j]);
        return;
    }

    treecorr_process_auto(context, tree, node->left, hist);
    treecorr_process_auto(context, tree, node->right, hist);
    treecorr_process_pair(context, tree, node->left,
                          tree, node->right, hist);
}
#endif /* TWOPCF */

#ifdef THREEPCFCONVERGENCE

enum {
    TREECORR_ZETA_COS = 0,
    TREECORR_ZETA_SIN = 1,
    TREECORR_ZETA_SINCOS = 2,
    TREECORR_ZETA_COSSIN = 3,
    TREECORR_ZETA_COMPONENTS = 4
};

enum {
    TREECORR_TRIPLE_OUTSIDE = 0,
    TREECORR_TRIPLE_SPLIT = 1,
    TREECORR_TRIPLE_ACCEPT = 2
};

typedef struct {
    real *components;
    real *normalization;
    size_t stride;
    size_t plane;
    size_t order_plane;
    int orders;
    INTEGER body_triples;
    INTEGER cell_triples;
} treecorr_triple_histogram;

static const unsigned char treecorr_permutations[6][3] = {
    {0, 1, 2}, {0, 2, 1},
    {1, 0, 2}, {1, 2, 0},
    {2, 0, 1}, {2, 1, 0}
};

static inline real *treecorr_zeta_component(
        treecorr_triple_histogram *hist, int component, int order)
{
    return hist->components
        + ((size_t)component * (size_t)hist->orders + (size_t)order)
        * hist->plane;
}

static inline real treecorr_position_distance(
        const treecorr_search_context *context,
        const cballs_storage_real *p, const cballs_storage_real *q,
        compute_vector dr)
{
    struct global_data *gd = context->gd;
    real distance2;

    DOTPSUBV(distance2, dr, p, q);
    if (context->cmd->usePeriodic) {
        VWrapAll(dr);
        DOTVP(distance2, dr, dr);
    }
    return distance2 > 0.0 ? rsqrt(distance2) : 0.0;
}

static inline bool treecorr_leg_outside(
        const treecorr_search_context *context, real distance, real size)
{
    return distance + size <= context->cmd->rminHist
        || distance - size >= context->cmd->rangeN;
}

static int treecorr_sloppy_bin(const treecorr_search_context *context,
                               real distance, real size)
{
    if (size > treecorr_bin_slop_width(context, distance)) return -1;
    return treecorr_bin_index(context, distance);
}

static bool treecorr_polar_coordinates_from_displacement(
        const treecorr_search_context *context,
        const cballs_storage_real *pivot,
        const cballs_storage_real *neighbor,
        const compute_vector dr, real distance,
        real *cosphi, real *sinphi)
{
    (void)context;
    (void)neighbor;
#if NDIM == 2
    *cosphi = -dr[0] / distance;
    *sinphi = -dr[1] / distance;
#else
#ifdef POLARAXIS
    {
        compute_vector q0 = {0.0, 0.0, 1.0};
        compute_vector dr0;
        compute_vector v21;
        compute_vector v31;
        compute_vector cross;
        real chord2;
        real chord;
        real a;
        real b;
        real denominator;
        real orientation;

        DOTPSUBV(chord2, dr0, pivot, q0);
        chord = chord2 > 0.0 ? rsqrt(chord2) : 0.0;
        b = 2.0 * rasin(MIN(1.0, 0.5 * chord));
#ifdef NOLIMBER
        a = 2.0 * rasin(MIN(1.0, 0.5 * distance));
#else
        a = distance;
#endif
        denominator = a * rsin(b);
        if (rabs(denominator) <= 16.0 * DBL_EPSILON) return FALSE;
        *cosphi = ((real)neighbor[2]
            - (1.0 - 0.5 * rsqr(a)) * rcos(b)) / denominator;
        *cosphi = MAX(-1.0, MIN(1.0, *cosphi));
        *sinphi = rsqrt(MAX(0.0, 1.0 - rsqr(*cosphi)));
        SUBV(v21, q0, pivot);
        SUBV(v31, neighbor, pivot);
        CROSSVP(cross, v21, v31);
        DOTVP(orientation, cross, pivot);
        if (!(orientation < 0.0)) *sinphi = -*sinphi;
    }
#else
    {
        compute_vector q0;
        compute_vector dr0;
        compute_vector cross;
        real reference2;
        real projection;
        real orientation;

        dRotation3D(pivot, ROTANGLE, ROTANGLE, ROTANGLE, q0);
        DOTPSUBV(reference2, dr0, pivot, q0);
        if (!(reference2 > 0.0)) return FALSE;
        DOTVP(projection, dr, dr0);
        *cosphi = projection / (distance * rsqrt(reference2));
        *cosphi = MAX(-1.0, MIN(1.0, *cosphi));
        *sinphi = rsqrt(MAX(0.0, 1.0 - rsqr(*cosphi)));
        CROSSVP(cross, dr0, pivot);
        DOTVP(orientation, dr, cross);
        if (orientation < 0.0) *sinphi = -*sinphi;
    }
#endif /* POLARAXIS */
#endif /* NDIM */
    return TRUE;
}

static bool treecorr_polar_coordinates(
        const treecorr_search_context *context,
        const cballs_storage_real *pivot,
        const cballs_storage_real *neighbor,
        real *distance, real *cosphi, real *sinphi)
{
    compute_vector dr;

    *distance = treecorr_position_distance(context, pivot, neighbor, dr);
    return *distance > 0.0
        && treecorr_polar_coordinates_from_displacement(
            context, pivot, neighbor, dr, *distance, cosphi, sinphi);
}

static void treecorr_initialize_angular_tolerance(
        treecorr_search_context *context)
{
    const real max_order = (real)MAX(context->cmd->mChebyshev, 0);

    /* TreeCorr's LogMultipole default is pi / (2 * max_n + 1). */
    context->angular_tolerance =
        context->cmd->theta * PI / (2.0 * max_order + 1.0);
    context->max_angular_ratio = context->angular_tolerance >= 0.5 * PI
        ? 1.0 : rsin(MAX(0.0, context->angular_tolerance));
}

static int treecorr_node_orientation_status(
        const treecorr_search_context *context,
        const fcfc_ballnode *pivot, const fcfc_ballnode *q,
        const fcfc_ballnode *r, int *nq, int *nr,
        real *cosq, real *sinq, real *cosr, real *sinr)
{
    compute_vector unused;
    real dq = treecorr_position_distance(context, pivot->center, q->center,
                                         unused);
    real dr = treecorr_position_distance(context, pivot->center, r->center,
                                         unused);
    const real sq = (real)pivot->radius + (real)q->radius;
    const real sr = (real)pivot->radius + (real)r->radius;
    real ignored;

    if (treecorr_leg_outside(context, dq, sq)
        || treecorr_leg_outside(context, dr, sr))
        return TREECORR_TRIPLE_OUTSIDE;
    if (!context->use_three_cells || !(dq > sq) || !(dr > sr))
        return TREECORR_TRIPLE_SPLIT;

    *nq = treecorr_sloppy_bin(context, dq, sq);
    *nr = treecorr_sloppy_bin(context, dr, sr);
    if (*nq < 0 || *nr < 0) return TREECORR_TRIPLE_SPLIT;

    if (sq / dq > context->max_angular_ratio
        || sr / dr > context->max_angular_ratio)
        return TREECORR_TRIPLE_SPLIT;
    if (!treecorr_polar_coordinates(context, pivot->center, q->center,
                                    &ignored, cosq, sinq)
        || !treecorr_polar_coordinates(context, pivot->center, r->center,
                                       &ignored, cosr, sinr))
        return TREECORR_TRIPLE_SPLIT;
    return TREECORR_TRIPLE_ACCEPT;
}

static void treecorr_accumulate_modes(treecorr_triple_histogram *hist,
                                      int nq, int nr, real numerator,
                                      real denominator,
                                      real cosq1, real sinq1,
                                      real cosr1, real sinr1)
{
    const size_t index = (size_t)nq * hist->stride + (size_t)nr;
    real cosq = 1.0;
    real sinq = 0.0;
    real cosr = 1.0;
    real sinr = 0.0;

    hist->normalization[index] += denominator;
    for (int order = 0; order < hist->orders; order++) {
        treecorr_zeta_component(hist, TREECORR_ZETA_COS, order)[index]
            += numerator * cosq * cosr;
        treecorr_zeta_component(hist, TREECORR_ZETA_SIN, order)[index]
            += numerator * sinq * sinr;
        treecorr_zeta_component(hist, TREECORR_ZETA_SINCOS, order)[index]
            += numerator * sinq * cosr;
        treecorr_zeta_component(hist, TREECORR_ZETA_COSSIN, order)[index]
            += numerator * cosq * sinr;
        {
            const real next_cosq = cosq * cosq1 - sinq * sinq1;
            const real next_sinq = sinq * cosq1 + cosq * sinq1;
            const real next_cosr = cosr * cosr1 - sinr * sinr1;
            const real next_sinr = sinr * cosr1 + cosr * sinr1;
            cosq = next_cosq;
            sinq = next_sinq;
            cosr = next_cosr;
            sinr = next_sinr;
        }
    }
}

static inline real treecorr_node_field_sum(
        const treecorr_search_context *context, const fcfc_ballnode *ball_node)
{
    return context->weighted
        ? ball_node->weighted_kappa_sum : ball_node->kappa_sum;
}

static inline real treecorr_node_normalization_sum(
        const treecorr_search_context *context, const fcfc_ballnode *ball_node)
{
    return context->weighted ? ball_node->field_weight_sum
                             : (real)treecorr_node_count(ball_node);
}

static void treecorr_accumulate_node_orientation(
        const treecorr_search_context *context,
        treecorr_triple_histogram *hist,
        const fcfc_ballnode *pivot, const fcfc_ballnode *q,
        const fcfc_ballnode *r, int nq, int nr,
        real cosq, real sinq, real cosr, real sinr)
{
    const real numerator = treecorr_node_field_sum(context, pivot)
                         * treecorr_node_field_sum(context, q)
                         * treecorr_node_field_sum(context, r);
    const real denominator = treecorr_node_normalization_sum(context, pivot)
                           * treecorr_node_normalization_sum(context, q)
                           * treecorr_node_normalization_sum(context, r);

    treecorr_accumulate_modes(hist, nq, nr, numerator, denominator,
                              cosq, sinq, cosr, sinr);
    hist->cell_triples++;
}

static void treecorr_accumulate_body_orientation(
        const treecorr_search_context *context,
        treecorr_triple_histogram *hist, bodyptr pivot, bodyptr q, bodyptr r)
{
    real dq;
    real dr;
    real cosq;
    real sinq;
    real cosr;
    real sinr;
    real fp;
    real fq;
    real fr;
    real denominator;
    int nq;
    int nr;

    if (!treecorr_polar_coordinates(context, Pos(pivot), Pos(q),
                                    &dq, &cosq, &sinq)
        || !treecorr_polar_coordinates(context, Pos(pivot), Pos(r),
                                       &dr, &cosr, &sinr))
        return;
    nq = treecorr_bin_index(context, dq);
    nr = treecorr_bin_index(context, dr);
    if (nq < 0 || nr < 0) return;

    fp = treecorr_body_field(pivot);
    fq = treecorr_body_field(q);
    fr = treecorr_body_field(r);
    if (context->weighted) {
        fp *= Weight(pivot);
        fq *= Weight(q);
        fr *= Weight(r);
        denominator = Weight(pivot) * Weight(q) * Weight(r);
    } else {
        denominator = 1.0;
    }
    treecorr_accumulate_modes(hist, nq, nr, fp * fq * fr, denominator,
                              cosq, sinq, cosr, sinr);
    hist->body_triples++;
}

static void treecorr_accumulate_body_mask(
        const treecorr_search_context *context,
        treecorr_triple_histogram *hist, bodyptr bodies[3], unsigned mask)
{
    for (int permutation = 0; permutation < 6; permutation++) {
        if (mask & (1U << permutation)) {
            const unsigned char *p = treecorr_permutations[permutation];
            treecorr_accumulate_body_orientation(
                context, hist, bodies[p[0]], bodies[p[1]], bodies[p[2]]);
        }
    }
}

#ifdef TREECORR_TASK_FRONTIER_ENGINE
static unsigned treecorr_request_leg_split(
        const fcfc_ballnode *pivot, int pivot_slot,
        const fcfc_ballnode *neighbor, int neighbor_slot)
{
    const bool pivot_leaf = treecorr_node_is_leaf(pivot);
    const bool neighbor_leaf = treecorr_node_is_leaf(neighbor);

    if (pivot_leaf && neighbor_leaf) return 0U;
    if (pivot_leaf) return 1U << neighbor_slot;
    if (neighbor_leaf) return 1U << pivot_slot;
    return (real)pivot->radius >= (real)neighbor->radius
        ? 1U << pivot_slot : 1U << neighbor_slot;
}

/* Select the cells responsible for a failed two-ball leg criterion. */
static unsigned treecorr_orientation_split_mask(
        const treecorr_search_context *context,
        const fcfc_ballnode *cells[3], const unsigned char permutation[3])
{
    const int pivot_slot = permutation[0];
    const fcfc_ballnode *pivot = cells[pivot_slot];
    unsigned split_mask = 0U;

    for (int leg = 1; leg <= 2; leg++) {
        const int neighbor_slot = permutation[leg];
        const fcfc_ballnode *neighbor = cells[neighbor_slot];
        compute_vector unused;
        const real distance = treecorr_position_distance(
            context, pivot->center, neighbor->center, unused);
        const real size = (real)pivot->radius + (real)neighbor->radius;
        bool leg_needs_split = !context->use_three_cells
            || !(distance > size)
            || treecorr_sloppy_bin(context, distance, size) < 0;

        if (!leg_needs_split) {
            leg_needs_split = size / distance
                > context->max_angular_ratio;
        }
        if (leg_needs_split)
            split_mask |= treecorr_request_leg_split(
                pivot, pivot_slot, neighbor, neighbor_slot);
    }
    return split_mask;
}

static void treecorr_process111_children(
        const treecorr_search_context *context,
        const fcfc_balltreeptr trees[3], const INTEGER nodes[3],
        unsigned orientation_mask, unsigned split_mask, int slot,
        treecorr_triple_histogram *hist);
#endif

static void treecorr_process111(const treecorr_search_context *context,
                                const fcfc_balltreeptr trees[3],
                                const INTEGER nodes[3], unsigned mask,
                                treecorr_triple_histogram *hist)
{
    const fcfc_ballnode *cells[3] = {
        &trees[0]->nodes[nodes[0]],
        &trees[1]->nodes[nodes[1]],
        &trees[2]->nodes[nodes[2]]
    };
    unsigned remaining = mask;
#ifdef TREECORR_TASK_FRONTIER_ENGINE
    unsigned requested_splits = 0U;
#endif

    for (int permutation = 0; permutation < 6; permutation++) {
        int nq;
        int nr;
        real cosq;
        real sinq;
        real cosr;
        real sinr;
        int orientation_status;
        const unsigned bit = 1U << permutation;
        const unsigned char *p;

        if (!(remaining & bit)) continue;
        p = treecorr_permutations[permutation];
        orientation_status = treecorr_node_orientation_status(
            context, cells[p[0]], cells[p[1]], cells[p[2]],
            &nq, &nr, &cosq, &sinq, &cosr, &sinr);
        if (orientation_status == TREECORR_TRIPLE_OUTSIDE) {
            remaining &= ~bit;
        } else if (orientation_status == TREECORR_TRIPLE_ACCEPT) {
            treecorr_accumulate_node_orientation(
                context, hist, cells[p[0]], cells[p[1]], cells[p[2]],
                nq, nr, cosq, sinq, cosr, sinr);
            remaining &= ~bit;
#ifdef TREECORR_TASK_FRONTIER_ENGINE
        } else {
            requested_splits |= treecorr_orientation_split_mask(
                context, cells, p);
#endif
        }
    }
    if (remaining == 0) return;

    if (treecorr_node_is_leaf(cells[0])
        && treecorr_node_is_leaf(cells[1])
        && treecorr_node_is_leaf(cells[2])) {
        for (INTEGER i = cells[0]->first; i <= cells[0]->last; i++)
            for (INTEGER j = cells[1]->first; j <= cells[1]->last; j++)
                for (INTEGER k = cells[2]->first; k <= cells[2]->last; k++) {
                    bodyptr bodies[3] = {
                        trees[0]->bptr[i], trees[1]->bptr[j], trees[2]->bptr[k]
                    };
                    treecorr_accumulate_body_mask(context, hist, bodies, remaining);
                }
        return;
    }

    {
#ifdef TREECORR_TASK_FRONTIER_ENGINE
        unsigned split_mask = 0U;
        real largest_requested = -1.0;

        for (int slot = 0; slot < 3; slot++) {
            if ((requested_splits & (1U << slot))
                && !treecorr_node_is_leaf(cells[slot]))
                largest_requested = MAX(
                    largest_requested, (real)cells[slot]->radius);
        }
        if (largest_requested >= 0.0) {
            for (int slot = 0; slot < 3; slot++) {
                if ((requested_splits & (1U << slot))
                    && !treecorr_node_is_leaf(cells[slot])
                    && (real)cells[slot]->radius
                        >= TREECORR_SPLIT_FACTOR * largest_requested)
                    split_mask |= 1U << slot;
            }
        }
        if (split_mask != 0U) {
            treecorr_process111_children(
                context, trees, nodes, remaining, split_mask, 0, hist);
            return;
        }
#endif
        int split = -1;
        real largest = -1.0;
        for (int slot = 0; slot < 3; slot++) {
            if (!treecorr_node_is_leaf(cells[slot])
                && (real)cells[slot]->radius > largest) {
                split = slot;
                largest = (real)cells[slot]->radius;
            }
        }
        if (split >= 0) {
            INTEGER child_nodes[3] = {nodes[0], nodes[1], nodes[2]};
            child_nodes[split] = cells[split]->left;
            treecorr_process111(context, trees, child_nodes, remaining, hist);
            child_nodes[split] = cells[split]->right;
            treecorr_process111(context, trees, child_nodes, remaining, hist);
        }
    }
}

#ifdef TREECORR_TASK_FRONTIER_ENGINE
static void treecorr_process111_children(
        const treecorr_search_context *context,
        const fcfc_balltreeptr trees[3], const INTEGER nodes[3],
        unsigned orientation_mask, unsigned split_mask, int slot,
        treecorr_triple_histogram *hist)
{
    if (slot == 3) {
        treecorr_process111(
            context, trees, nodes, orientation_mask, hist);
        return;
    }
    if (split_mask & (1U << slot)) {
        INTEGER child_nodes[3] = {nodes[0], nodes[1], nodes[2]};
        const fcfc_ballnode *ball_node =
            &trees[slot]->nodes[nodes[slot]];

        child_nodes[slot] = ball_node->left;
        treecorr_process111_children(
            context, trees, child_nodes, orientation_mask,
            split_mask, slot + 1, hist);
        child_nodes[slot] = ball_node->right;
        treecorr_process111_children(
            context, trees, child_nodes, orientation_mask,
            split_mask, slot + 1, hist);
    } else {
        treecorr_process111_children(
            context, trees, nodes, orientation_mask,
            split_mask, slot + 1, hist);
    }
}
#endif

static void treecorr_process21_auto(
        const treecorr_search_context *context, const fcfc_balltreeptr tree,
        INTEGER two_index, INTEGER one_index,
        treecorr_triple_histogram *hist)
{
    const fcfc_ballnode *two = &tree->nodes[two_index];
    const fcfc_ballnode *one = &tree->nodes[one_index];

    if (treecorr_node_count(two) < 2) return;
    {
        const real distance = treecorr_center_distance(
            context->cmd, context->gd, two, one);
        const real size = (real)two->radius + (real)one->radius;
        if (treecorr_leg_outside(context, distance, size)) return;
    }
    if (treecorr_node_is_leaf(two) && treecorr_node_is_leaf(one)) {
        for (INTEGER i = two->first; i < two->last; i++)
            for (INTEGER j = i + 1; j <= two->last; j++)
                for (INTEGER k = one->first; k <= one->last; k++) {
                    bodyptr bodies[3] = {
                        tree->bptr[i], tree->bptr[j], tree->bptr[k]
                    };
                    treecorr_accumulate_body_mask(context, hist, bodies, 0x3fU);
                }
        return;
    }

    if (!treecorr_node_is_leaf(one)
        && (treecorr_node_is_leaf(two) || one->radius > two->radius)) {
        treecorr_process21_auto(context, tree, two_index, one->left, hist);
        treecorr_process21_auto(context, tree, two_index, one->right, hist);
    } else {
        const fcfc_balltreeptr trees[3] = {tree, tree, tree};
        const INTEGER nodes[3] = {two->left, two->right, one_index};
        treecorr_process21_auto(context, tree, two->left, one_index, hist);
        treecorr_process21_auto(context, tree, two->right, one_index, hist);
        treecorr_process111(context, trees, nodes, 0x3fU, hist);
    }
}

static void treecorr_process3_auto(const treecorr_search_context *context,
                                   const fcfc_balltreeptr tree, INTEGER inode,
                                   treecorr_triple_histogram *hist)
{
    const fcfc_ballnode *ball_node = &tree->nodes[inode];

    if (treecorr_node_count(ball_node) < 3) return;
    if (2.0 * (real)ball_node->radius <= context->cmd->rminHist) return;
    if (treecorr_node_is_leaf(ball_node)) {
        for (INTEGER i = ball_node->first; i < ball_node->last - 1; i++)
            for (INTEGER j = i + 1; j < ball_node->last; j++)
                for (INTEGER k = j + 1; k <= ball_node->last; k++) {
                    bodyptr bodies[3] = {
                        tree->bptr[i], tree->bptr[j], tree->bptr[k]
                    };
                    treecorr_accumulate_body_mask(context, hist, bodies, 0x3fU);
                }
        return;
    }
    treecorr_process3_auto(context, tree, ball_node->left, hist);
    treecorr_process3_auto(context, tree, ball_node->right, hist);
    treecorr_process21_auto(
        context, tree, ball_node->left, ball_node->right, hist);
    treecorr_process21_auto(
        context, tree, ball_node->right, ball_node->left, hist);
}

static void treecorr_process12_cross(
        const treecorr_search_context *context,
        const fcfc_balltreeptr pivot_tree, INTEGER pivot_index,
        const fcfc_balltreeptr pair_tree, INTEGER pair_index,
        treecorr_triple_histogram *hist)
{
    const fcfc_ballnode *pivot = &pivot_tree->nodes[pivot_index];
    const fcfc_ballnode *pair = &pair_tree->nodes[pair_index];

    if (treecorr_node_count(pair) < 2) return;
    {
        const real distance = treecorr_center_distance(
            context->cmd, context->gd, pivot, pair);
        const real size = (real)pivot->radius + (real)pair->radius;
        if (treecorr_leg_outside(context, distance, size)) return;
    }
    if (treecorr_node_is_leaf(pivot) && treecorr_node_is_leaf(pair)) {
        for (INTEGER i = pivot->first; i <= pivot->last; i++)
            for (INTEGER j = pair->first; j < pair->last; j++)
                for (INTEGER k = j + 1; k <= pair->last; k++) {
                    bodyptr bodies[3] = {
                        pivot_tree->bptr[i], pair_tree->bptr[j], pair_tree->bptr[k]
                    };
                    treecorr_accumulate_body_mask(context, hist, bodies, 0x03U);
                }
        return;
    }

    if (!treecorr_node_is_leaf(pivot)
        && (treecorr_node_is_leaf(pair) || pivot->radius > pair->radius)) {
        treecorr_process12_cross(context, pivot_tree, pivot->left,
                                 pair_tree, pair_index, hist);
        treecorr_process12_cross(context, pivot_tree, pivot->right,
                                 pair_tree, pair_index, hist);
    } else {
        const fcfc_balltreeptr trees[3] = {
            pivot_tree, pair_tree, pair_tree
        };
        const INTEGER nodes[3] = {pivot_index, pair->left, pair->right};
        treecorr_process12_cross(context, pivot_tree, pivot_index,
                                 pair_tree, pair->left, hist);
        treecorr_process12_cross(context, pivot_tree, pivot_index,
                                 pair_tree, pair->right, hist);
        treecorr_process111(context, trees, nodes, 0x03U, hist);
    }
}

static int treecorr_allocate_triple_histograms(
        struct cmdline_data *cmd, INTEGER task_count, size_t stride,
        int orders, real **histograms, size_t *values_per_task,
        INTEGER **body_triples, INTEGER **cell_triples)
{
    size_t tasks;
    size_t plane;
    size_t component_values;
    size_t values;

    *histograms = NULL;
    *values_per_task = 0;
    *body_triples = NULL;
    *cell_triples = NULL;
    if (task_count <= 0 || orders <= 0
        || stride > SIZE_MAX / stride)
        goto invalid_size;
    tasks = (size_t)task_count;
    plane = stride * stride;
    if ((size_t)orders > SIZE_MAX / TREECORR_ZETA_COMPONENTS
        || (size_t)orders * TREECORR_ZETA_COMPONENTS > SIZE_MAX / plane)
        goto invalid_size;
    component_values = (size_t)orders * TREECORR_ZETA_COMPONENTS * plane;
    if (component_values > SIZE_MAX - plane)
        goto invalid_size;
    *values_per_task = component_values + plane;
    if (tasks > SIZE_MAX / *values_per_task
        || (values = tasks * *values_per_task) > SIZE_MAX / sizeof(**histograms))
        goto invalid_size;
    *histograms = calloc(values, sizeof(**histograms));
    *body_triples = calloc(tasks, sizeof(**body_triples));
    *cell_triples = calloc(tasks, sizeof(**cell_triples));
    if (*histograms == NULL || *body_triples == NULL || *cell_triples == NULL) {
        free(*cell_triples);
        free(*body_triples);
        free(*histograms);
        *cell_triples = NULL;
        *body_triples = NULL;
        *histograms = NULL;
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: triple histogram allocation failed",
                 TREECORR_METHOD_NAME);
        return FAILURE;
    }
    return SUCCESS;

invalid_size:
    snprintf(cmd->error_message, _ERRORMSGSIZE_,
             "%s: triple histogram size overflow", TREECORR_METHOD_NAME);
    return FAILURE;
}

#ifdef TREECORR_TASK_FRONTIER_ENGINE
enum {
    TREECORR_TASK_AUTO3 = 0,
    TREECORR_TASK_TWO_ONE = 1,
    TREECORR_TASK_THREE = 2,
    TREECORR_TASK_CROSS12 = 3
};

typedef struct {
    int kind;
    fcfc_balltreeptr trees[3];
    INTEGER nodes[3];
    unsigned orientation_mask;
    long double estimated_work;
} treecorr_triple_task;

static unsigned treecorr_popcount(unsigned value)
{
    unsigned count = 0;

    while (value != 0U) {
        count += value & 1U;
        value >>= 1;
    }
    return count;
}

static long double treecorr_task_work(const treecorr_triple_task *task)
{
    const long double n0 = (long double)treecorr_node_count(
        &task->trees[0]->nodes[task->nodes[0]]);

    switch (task->kind) {
        case TREECORR_TASK_AUTO3:
            return n0 >= 3.0L ? n0 * (n0 - 1.0L) * (n0 - 2.0L) : 0.0L;
        case TREECORR_TASK_TWO_ONE: {
            const long double n1 = (long double)treecorr_node_count(
                &task->trees[0]->nodes[task->nodes[1]]);
            return n0 >= 2.0L
                ? 3.0L * n0 * (n0 - 1.0L) * n1 : 0.0L;
        }
        case TREECORR_TASK_THREE: {
            const long double n1 = (long double)treecorr_node_count(
                &task->trees[1]->nodes[task->nodes[1]]);
            const long double n2 = (long double)treecorr_node_count(
                &task->trees[2]->nodes[task->nodes[2]]);
            return n0 * n1 * n2
                * (long double)treecorr_popcount(task->orientation_mask);
        }
        case TREECORR_TASK_CROSS12: {
            const long double n1 = (long double)treecorr_node_count(
                &task->trees[1]->nodes[task->nodes[1]]);
            return n1 >= 2.0L ? n0 * n1 * (n1 - 1.0L) : 0.0L;
        }
        default:
            return 0.0L;
    }
}

static treecorr_triple_task treecorr_make_task(
        int kind, fcfc_balltreeptr tree0, INTEGER node0,
        fcfc_balltreeptr tree1, INTEGER node1,
        fcfc_balltreeptr tree2, INTEGER node2,
        unsigned orientation_mask)
{
    treecorr_triple_task task;

    task.kind = kind;
    task.trees[0] = tree0;
    task.trees[1] = tree1;
    task.trees[2] = tree2;
    task.nodes[0] = node0;
    task.nodes[1] = node1;
    task.nodes[2] = node2;
    task.orientation_mask = orientation_mask;
    task.estimated_work = 0.0L;
    task.estimated_work = treecorr_task_work(&task);
    return task;
}

static int treecorr_split_task(const treecorr_triple_task *task,
                               treecorr_triple_task children[4])
{
    const fcfc_ballnode *node0 =
        &task->trees[0]->nodes[task->nodes[0]];

    switch (task->kind) {
        case TREECORR_TASK_AUTO3:
            if (treecorr_node_is_leaf(node0)) return 0;
            children[0] = treecorr_make_task(
                TREECORR_TASK_AUTO3, task->trees[0], node0->left,
                NULL, 0, NULL, 0, 0U);
            children[1] = treecorr_make_task(
                TREECORR_TASK_AUTO3, task->trees[0], node0->right,
                NULL, 0, NULL, 0, 0U);
            children[2] = treecorr_make_task(
                TREECORR_TASK_TWO_ONE, task->trees[0], node0->left,
                task->trees[0], node0->right, NULL, 0, 0U);
            children[3] = treecorr_make_task(
                TREECORR_TASK_TWO_ONE, task->trees[0], node0->right,
                task->trees[0], node0->left, NULL, 0, 0U);
            return 4;

        case TREECORR_TASK_TWO_ONE: {
            const fcfc_ballnode *one =
                &task->trees[0]->nodes[task->nodes[1]];

            if (treecorr_node_count(node0) < 2
                || (treecorr_node_is_leaf(node0)
                    && treecorr_node_is_leaf(one)))
                return 0;
            if (!treecorr_node_is_leaf(one)
                && (treecorr_node_is_leaf(node0)
                    || one->radius > node0->radius)) {
                children[0] = treecorr_make_task(
                    TREECORR_TASK_TWO_ONE,
                    task->trees[0], task->nodes[0],
                    task->trees[0], one->left, NULL, 0, 0U);
                children[1] = treecorr_make_task(
                    TREECORR_TASK_TWO_ONE,
                    task->trees[0], task->nodes[0],
                    task->trees[0], one->right, NULL, 0, 0U);
                return 2;
            }
            children[0] = treecorr_make_task(
                TREECORR_TASK_TWO_ONE, task->trees[0], node0->left,
                task->trees[0], task->nodes[1], NULL, 0, 0U);
            children[1] = treecorr_make_task(
                TREECORR_TASK_TWO_ONE, task->trees[0], node0->right,
                task->trees[0], task->nodes[1], NULL, 0, 0U);
            children[2] = treecorr_make_task(
                TREECORR_TASK_THREE, task->trees[0], node0->left,
                task->trees[0], node0->right,
                task->trees[0], task->nodes[1], 0x3fU);
            return 3;
        }

        case TREECORR_TASK_CROSS12: {
            const fcfc_ballnode *pair =
                &task->trees[1]->nodes[task->nodes[1]];

            if (treecorr_node_count(pair) < 2
                || (treecorr_node_is_leaf(node0)
                    && treecorr_node_is_leaf(pair)))
                return 0;
            if (!treecorr_node_is_leaf(node0)
                && (treecorr_node_is_leaf(pair)
                    || node0->radius > pair->radius)) {
                children[0] = treecorr_make_task(
                    TREECORR_TASK_CROSS12, task->trees[0], node0->left,
                    task->trees[1], task->nodes[1], NULL, 0, 0U);
                children[1] = treecorr_make_task(
                    TREECORR_TASK_CROSS12, task->trees[0], node0->right,
                    task->trees[1], task->nodes[1], NULL, 0, 0U);
                return 2;
            }
            children[0] = treecorr_make_task(
                TREECORR_TASK_CROSS12, task->trees[0], task->nodes[0],
                task->trees[1], pair->left, NULL, 0, 0U);
            children[1] = treecorr_make_task(
                TREECORR_TASK_CROSS12, task->trees[0], task->nodes[0],
                task->trees[1], pair->right, NULL, 0, 0U);
            children[2] = treecorr_make_task(
                TREECORR_TASK_THREE, task->trees[0], task->nodes[0],
                task->trees[1], pair->left,
                task->trees[1], pair->right, 0x03U);
            return 3;
        }

        case TREECORR_TASK_THREE: {
            int split = -1;
            real largest = -1.0;

            for (int slot = 0; slot < 3; slot++) {
                const fcfc_ballnode *ball_node =
                    &task->trees[slot]->nodes[task->nodes[slot]];
                if (!treecorr_node_is_leaf(ball_node)
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
                children[child] = treecorr_make_task(
                    TREECORR_TASK_THREE,
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

static INTEGER treecorr_task_target(size_t stride, int orders)
{
    const size_t plane = stride <= SIZE_MAX / stride
        ? stride * stride : SIZE_MAX;
    size_t values;
    size_t bytes;
    INTEGER target = TREECORR_TRIPLE_TASK_TARGET;

    if (orders <= 0 || plane == SIZE_MAX
        || (size_t)orders > SIZE_MAX / TREECORR_ZETA_COMPONENTS
        || (size_t)orders * TREECORR_ZETA_COMPONENTS > SIZE_MAX / plane)
        return 1;
    values = (size_t)orders * TREECORR_ZETA_COMPONENTS * plane;
    if (values > SIZE_MAX - plane) return 1;
    values += plane;
    if (values > SIZE_MAX / sizeof(real)) return 1;
    bytes = values * sizeof(real);
    if (bytes > 0 && (size_t)target > TREECORR_TRIPLE_TASK_MEMORY / bytes)
        target = (INTEGER)MAX(
            (size_t)1, TREECORR_TRIPLE_TASK_MEMORY / bytes);
    return target;
}

static int treecorr_build_task_frontier(
        struct cmdline_data *cmd, fcfc_balltreeptr tree1,
        fcfc_balltreeptr tree2, bool auto_correlation, INTEGER target,
        treecorr_triple_task **result, INTEGER *result_count)
{
    treecorr_triple_task *tasks;
    INTEGER count = 1;

    *result = NULL;
    *result_count = 0;
    if (target <= 0 || (size_t)target > SIZE_MAX / sizeof(*tasks)) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: invalid task-frontier size", TREECORR_METHOD_NAME);
        return FAILURE;
    }
    tasks = calloc((size_t)target, sizeof(*tasks));
    if (tasks == NULL) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: task-frontier allocation failed", TREECORR_METHOD_NAME);
        return FAILURE;
    }
    tasks[0] = auto_correlation
        ? treecorr_make_task(TREECORR_TASK_AUTO3, tree1, 0,
                             NULL, 0, NULL, 0, 0U)
        : treecorr_make_task(TREECORR_TASK_CROSS12, tree1, 0,
                             tree2, 0, NULL, 0, 0U);

    while (count < target) {
        INTEGER best = -1;
        int best_children = 0;
        long double largest_work = -1.0L;

        for (INTEGER i = 0; i < count; i++) {
            treecorr_triple_task children[4];
            const int nchildren = treecorr_split_task(&tasks[i], children);

            if (nchildren <= 0 || count + nchildren - 1 > target) continue;
            if (tasks[i].estimated_work > largest_work) {
                best = i;
                best_children = nchildren;
                largest_work = tasks[i].estimated_work;
            }
        }
        if (best < 0) break;
        {
            treecorr_triple_task children[4];
            const int nchildren = treecorr_split_task(&tasks[best], children);

            if (nchildren != best_children) {
                free(tasks);
                snprintf(cmd->error_message, _ERRORMSGSIZE_,
                         "%s: inconsistent task split", TREECORR_METHOD_NAME);
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

static void treecorr_run_task(const treecorr_search_context *context,
                              const treecorr_triple_task *task,
                              treecorr_triple_histogram *hist)
{
    switch (task->kind) {
        case TREECORR_TASK_AUTO3:
            treecorr_process3_auto(
                context, task->trees[0], task->nodes[0], hist);
            break;
        case TREECORR_TASK_TWO_ONE:
            treecorr_process21_auto(
                context, task->trees[0], task->nodes[0],
                task->nodes[1], hist);
            break;
        case TREECORR_TASK_THREE:
            treecorr_process111(
                context, task->trees, task->nodes,
                task->orientation_mask, hist);
            break;
        case TREECORR_TASK_CROSS12:
            treecorr_process12_cross(
                context, task->trees[0], task->nodes[0],
                task->trees[1], task->nodes[1], hist);
            break;
    }
}

#ifdef TREECORR_LOG_MULTIPOLE_ENGINE
static int treecorr_search_direct_triples(
        struct cmdline_data *, struct global_data *,
        bodyptr *, INTEGER *, INTEGER, INTEGER *, int, int);
#ifdef THREEPCFCONVERGENCE

typedef struct {
    real *field_cos;
    real *field_sin;
    real *self_coscos;
    real *self_sinsin;
    real *self_sincos;
    real *normalization;
    real *normalization_sq;
    size_t stride;
    size_t values;
    int orders;
    INTEGER body_visits;
    INTEGER accepted_nodes;
    INTEGER pair_tests;
    INTEGER pivot_restarts;
    INTEGER pivot_finishes;
} treecorr_multipole_scratch;

typedef struct {
    const cballs_storage_real *position;
    real radius;
    real field_sum;
    real normalization_sum;
    INTEGER first;
    INTEGER last;
    bool can_split;
} treecorr_multipole_pivot;

static void treecorr_initialize_multipole_scratch(
        treecorr_multipole_scratch *, real *, size_t, int, size_t);

static inline real treecorr_node_field_sq_sum(
        const treecorr_search_context *context,
        const fcfc_ballnode *ball_node)
{
    return context->weighted ? ball_node->weighted_kappa_sq_sum
                             : ball_node->kappa_sq_sum;
}

static inline real treecorr_node_normalization_sq_sum(
        const treecorr_search_context *context,
        const fcfc_ballnode *ball_node)
{
    return context->weighted ? ball_node->field_weight_sq_sum
                             : (real)treecorr_node_count(ball_node);
}

static inline real treecorr_body_normalization(
        const treecorr_search_context *context, bodyptr p)
{
    return context->weighted ? Weight(p) : 1.0;
}

static inline real treecorr_body_weighted_field(
        const treecorr_search_context *context, bodyptr p)
{
    const real field = treecorr_body_field(p);

    return context->weighted ? Weight(p) * field : field;
}

static void treecorr_multipole_clear(treecorr_multipole_scratch *scratch)
{
    memset(scratch->field_cos, 0, scratch->values * sizeof(real));
}

static void treecorr_multipole_add(
        treecorr_multipole_scratch *scratch, int radial_bin,
        real field_sum, real field_sq_sum,
        real normalization_sum, real normalization_sq_sum,
        real cosphi1, real sinphi1)
{
    real cosphi = cosphi1;
    real sinphi = sinphi1;
    size_t index = (size_t)radial_bin;

    scratch->normalization[radial_bin] += normalization_sum;
    scratch->normalization_sq[radial_bin] += normalization_sq_sum;
    scratch->field_cos[index] += field_sum;
    scratch->self_coscos[index] += field_sq_sum;
    for (int order = 1; order < scratch->orders; order++) {
        const real field_cos = field_sum * cosphi;
        const real field_sin = field_sum * sinphi;

        index += scratch->stride;
        scratch->field_cos[index] += field_cos;
        scratch->field_sin[index] += field_sin;
        scratch->self_coscos[index] += field_sq_sum * cosphi * cosphi;
        scratch->self_sinsin[index] += field_sq_sum * sinphi * sinphi;
        scratch->self_sincos[index] += field_sq_sum * sinphi * cosphi;
        {
            const real next_cos = cosphi * cosphi1 - sinphi * sinphi1;
            const real next_sin = sinphi * cosphi1 + cosphi * sinphi1;
            cosphi = next_cos;
            sinphi = next_sin;
        }
    }
}

static int treecorr_multipole_pair_status_limited(
        const treecorr_search_context *context,
        const treecorr_multipole_pivot *pivot,
        const cballs_storage_real *neighbor_position, real neighbor_radius,
        real upper_limit, int *radial_bin,
        real *cosphi, real *sinphi, real *upper_extent)
{
    compute_vector dr;
    real distance;
    const real size = pivot->radius + neighbor_radius;
    const real lower_limit = context->cmd->rminHist;

    distance = treecorr_position_distance(
        context, pivot->position, neighbor_position, dr);
    *upper_extent = MIN(upper_limit, distance + size);
    if (!(distance > 0.0))
        return TREECORR_TRIPLE_SPLIT;
    if (distance + size <= lower_limit
        || distance - size >= upper_limit)
        return TREECORR_TRIPLE_OUTSIDE;
    if (!context->use_two_balls || !(distance > size)
        || size > treecorr_bin_slop_width(context, distance)
        || size / distance > context->max_angular_ratio)
        return TREECORR_TRIPLE_SPLIT;
    *radial_bin = treecorr_bin_index(context, distance);
    if (*radial_bin < 0 || !(distance < upper_limit))
        return TREECORR_TRIPLE_OUTSIDE;
    if (!treecorr_polar_coordinates_from_displacement(
            context, pivot->position, neighbor_position,
            dr, distance, cosphi, sinphi))
        return TREECORR_TRIPLE_SPLIT;
    return TREECORR_TRIPLE_ACCEPT;
}

static int treecorr_multipole_pair_status(
        const treecorr_search_context *context,
        const treecorr_multipole_pivot *pivot,
        const cballs_storage_real *neighbor_position, real neighbor_radius,
        int *radial_bin, real *cosphi, real *sinphi)
{
    real upper_extent;

    return treecorr_multipole_pair_status_limited(
        context, pivot, neighbor_position, neighbor_radius,
        context->cmd->rangeN, radial_bin, cosphi, sinphi, &upper_extent);
}

static int treecorr_multipole_add_body(
        const treecorr_search_context *context,
        const treecorr_multipole_pivot *pivot,
        treecorr_multipole_scratch *scratch, bodyptr neighbor)
{
    real distance;
    real cosphi;
    real sinphi;
    int radial_bin;

    if (pivot->radius > 0.0) {
        scratch->pair_tests++;
        const int pair_status = treecorr_multipole_pair_status(
            context, pivot, Pos(neighbor), 0.0,
            &radial_bin, &cosphi, &sinphi);
        if (pair_status == TREECORR_TRIPLE_OUTSIDE) return SUCCESS;
        if (pair_status != TREECORR_TRIPLE_ACCEPT) return FAILURE;
    } else {
        if (!treecorr_polar_coordinates(
                context, pivot->position, Pos(neighbor),
                &distance, &cosphi, &sinphi))
            return SUCCESS;
        radial_bin = treecorr_bin_index(context, distance);
        if (radial_bin < 0) return SUCCESS;
    }
    {
        const real field = treecorr_body_weighted_field(context, neighbor);
        const real normalization = treecorr_body_normalization(
            context, neighbor);

        treecorr_multipole_add(
            scratch, radial_bin, field, field * field,
            normalization, normalization * normalization,
            cosphi, sinphi);
    }
    scratch->body_visits++;
    return SUCCESS;
}

static bool treecorr_ranges_overlap(INTEGER first1, INTEGER last1,
                                    INTEGER first2, INTEGER last2)
{
    return first1 <= last2 && first2 <= last1;
}

/*
 * Scan one neighbor tree for a pivot ball.  FAILURE means that the neighbor
 * side can no longer resolve the requested radial/angular accuracy, so the
 * caller must discard this scratch buffer and split the pivot ball.
 */
static int treecorr_multipole_scan_neighbors(
        const treecorr_search_context *context,
        const treecorr_multipole_pivot *pivot,
        const fcfc_balltreeptr neighbor_tree, INTEGER neighbor_index,
        bool same_tree, treecorr_multipole_scratch *scratch)
{
    const fcfc_ballnode *neighbor = &neighbor_tree->nodes[neighbor_index];
    int radial_bin;
    real cosphi;
    real sinphi;
    int pair_status;

    if (same_tree && treecorr_ranges_overlap(
            pivot->first, pivot->last, neighbor->first, neighbor->last)) {
        if (neighbor->first >= pivot->first
            && neighbor->last <= pivot->last)
            return SUCCESS;
        if (treecorr_node_is_leaf(neighbor)) {
            for (INTEGER i = neighbor->first; i <= neighbor->last; i++) {
                if (i >= pivot->first && i <= pivot->last) continue;
                if (treecorr_multipole_add_body(
                        context, pivot, scratch,
                        neighbor_tree->bptr[i]) == FAILURE)
                    return FAILURE;
            }
            return SUCCESS;
        }
        if (treecorr_multipole_scan_neighbors(
                context, pivot, neighbor_tree, neighbor->left,
                same_tree, scratch) == FAILURE)
            return FAILURE;
        return treecorr_multipole_scan_neighbors(
            context, pivot, neighbor_tree, neighbor->right,
            same_tree, scratch);
    }

    scratch->pair_tests++;
    pair_status = treecorr_multipole_pair_status(
        context, pivot, neighbor->center, (real)neighbor->radius,
        &radial_bin, &cosphi, &sinphi);
    if (pair_status == TREECORR_TRIPLE_OUTSIDE) return SUCCESS;
    if (pair_status == TREECORR_TRIPLE_ACCEPT) {
        treecorr_multipole_add(
            scratch, radial_bin,
            treecorr_node_field_sum(context, neighbor),
            treecorr_node_field_sq_sum(context, neighbor),
            treecorr_node_normalization_sum(context, neighbor),
            treecorr_node_normalization_sq_sum(context, neighbor),
            cosphi, sinphi);
        scratch->accepted_nodes++;
        return SUCCESS;
    }

    if (!treecorr_node_is_leaf(neighbor)) {
        if (pivot->can_split
            && pivot->radius > (real)neighbor->radius)
            return FAILURE;
        if (treecorr_multipole_scan_neighbors(
                context, pivot, neighbor_tree, neighbor->left,
                same_tree, scratch) == FAILURE)
            return FAILURE;
        return treecorr_multipole_scan_neighbors(
            context, pivot, neighbor_tree, neighbor->right,
            same_tree, scratch);
    }

    for (INTEGER i = neighbor->first; i <= neighbor->last; i++) {
        if (same_tree && i >= pivot->first && i <= pivot->last) continue;
        if (treecorr_multipole_add_body(
                context, pivot, scratch, neighbor_tree->bptr[i]) == FAILURE)
            return FAILURE;
    }
    return SUCCESS;
}

static void treecorr_multipole_finish_pivot_range(
        treecorr_triple_histogram *hist,
        const treecorr_multipole_pivot *pivot,
        const treecorr_multipole_scratch *scratch,
        int first_radial_bin, int radial_bin_limit, int radial_bins)
{
    for (int n1 = first_radial_bin; n1 < radial_bin_limit; n1++) {
        for (int n2 = n1; n2 <= radial_bins; n2++) {
            const size_t hist_index = (size_t)n1 * hist->stride
                                    + (size_t)n2;
            const size_t transpose_index = (size_t)n2 * hist->stride
                                         + (size_t)n1;
            real denominator = scratch->normalization[n1]
                             * scratch->normalization[n2];

            if (n1 == n2) denominator -= scratch->normalization_sq[n1];
            hist->normalization[hist_index] +=
                pivot->normalization_sum * denominator;
            if (n1 != n2)
                hist->normalization[transpose_index] +=
                    pivot->normalization_sum * denominator;
            {
                real coscos = scratch->field_cos[n1]
                            * scratch->field_cos[n2];

                if (n1 == n2)
                    coscos -= scratch->self_coscos[n1];
                treecorr_zeta_component(
                    hist, TREECORR_ZETA_COS, 0)[hist_index]
                    += pivot->field_sum * coscos;
                if (n1 != n2)
                    treecorr_zeta_component(
                        hist, TREECORR_ZETA_COS, 0)[transpose_index]
                        += pivot->field_sum * coscos;
            }
            for (int order = 1; order < scratch->orders; order++) {
                const size_t index1 = (size_t)order * scratch->stride
                                    + (size_t)n1;
                const size_t index2 = (size_t)order * scratch->stride
                                    + (size_t)n2;
                real coscos = scratch->field_cos[index1]
                            * scratch->field_cos[index2];
                real sinsin = scratch->field_sin[index1]
                            * scratch->field_sin[index2];
                real sincos = scratch->field_sin[index1]
                            * scratch->field_cos[index2];
                real cossin = scratch->field_cos[index1]
                            * scratch->field_sin[index2];

                if (n1 == n2) {
                    coscos -= scratch->self_coscos[index1];
                    sinsin -= scratch->self_sinsin[index1];
                    sincos -= scratch->self_sincos[index1];
                    cossin -= scratch->self_sincos[index1];
                }
                treecorr_zeta_component(
                    hist, TREECORR_ZETA_COS, order)[hist_index]
                    += pivot->field_sum * coscos;
                treecorr_zeta_component(
                    hist, TREECORR_ZETA_SIN, order)[hist_index]
                    += pivot->field_sum * sinsin;
                treecorr_zeta_component(
                    hist, TREECORR_ZETA_SINCOS, order)[hist_index]
                    += pivot->field_sum * sincos;
                treecorr_zeta_component(
                    hist, TREECORR_ZETA_COSSIN, order)[hist_index]
                    += pivot->field_sum * cossin;
                if (n1 != n2) {
                    treecorr_zeta_component(
                        hist, TREECORR_ZETA_COS, order)[transpose_index]
                        += pivot->field_sum * coscos;
                    treecorr_zeta_component(
                        hist, TREECORR_ZETA_SIN, order)[transpose_index]
                        += pivot->field_sum * sinsin;
                    treecorr_zeta_component(
                        hist, TREECORR_ZETA_SINCOS, order)[transpose_index]
                        += pivot->field_sum * cossin;
                    treecorr_zeta_component(
                        hist, TREECORR_ZETA_COSSIN, order)[transpose_index]
                        += pivot->field_sum * sincos;
                }
            }
        }
    }
}

static void treecorr_multipole_finish_pivot(
        treecorr_triple_histogram *hist,
        const treecorr_multipole_pivot *pivot,
        treecorr_multipole_scratch *scratch, int radial_bins)
{
    treecorr_multipole_finish_pivot_range(
        hist, pivot, scratch, 1, radial_bins + 1, radial_bins);
    scratch->pivot_finishes++;
}

static void treecorr_multipole_process_body(
        const treecorr_search_context *context,
        const fcfc_balltreeptr pivot_tree, INTEGER pivot_index,
        const fcfc_balltreeptr neighbor_tree, bool same_tree,
        treecorr_multipole_scratch *scratch,
        treecorr_triple_histogram *hist)
{
    bodyptr pivot_body = pivot_tree->bptr[pivot_index];
    treecorr_multipole_pivot pivot;

    pivot.position = Pos(pivot_body);
    pivot.radius = 0.0;
    pivot.field_sum = treecorr_body_weighted_field(context, pivot_body);
    pivot.normalization_sum = treecorr_body_normalization(
        context, pivot_body);
    pivot.first = pivot_index;
    pivot.last = pivot_index;
    pivot.can_split = FALSE;
    treecorr_multipole_clear(scratch);
    if (treecorr_multipole_scan_neighbors(
            context, &pivot, neighbor_tree, 0,
            same_tree, scratch) == SUCCESS)
        treecorr_multipole_finish_pivot(
            hist, &pivot, scratch, context->cmd->sizeHistN);
}

#ifdef TREECORR_BODY_PIVOT_LOG_MULTIPOLE
static void treecorr_multipole_process_body_pivots(
        const treecorr_search_context *context,
        const fcfc_balltreeptr pivot_tree, INTEGER pivot_node_index,
        const fcfc_balltreeptr neighbor_tree, bool same_tree,
        treecorr_multipole_scratch *scratch,
        treecorr_triple_histogram *hist)
{
    const fcfc_ballnode *pivot_node =
        &pivot_tree->nodes[pivot_node_index];

    if (!treecorr_node_is_leaf(pivot_node)) {
        treecorr_multipole_process_body_pivots(
            context, pivot_tree, pivot_node->left,
            neighbor_tree, same_tree, scratch, hist);
        treecorr_multipole_process_body_pivots(
            context, pivot_tree, pivot_node->right,
            neighbor_tree, same_tree, scratch, hist);
        return;
    }
    for (INTEGER i = pivot_node->first; i <= pivot_node->last; i++)
        treecorr_multipole_process_body(
            context, pivot_tree, i, neighbor_tree,
            same_tree, scratch, hist);
}
#endif

static void treecorr_multipole_process_pivots(
        const treecorr_search_context *context,
        const fcfc_balltreeptr pivot_tree, INTEGER pivot_node_index,
        const fcfc_balltreeptr neighbor_tree, bool same_tree,
        treecorr_multipole_scratch *scratch,
        treecorr_triple_histogram *hist)
{
    const fcfc_ballnode *pivot_node =
        &pivot_tree->nodes[pivot_node_index];

    if (!context->use_two_balls
        || (same_tree
            && 2.0 * (real)pivot_node->radius > context->cmd->rminHist)) {
        if (!treecorr_node_is_leaf(pivot_node)) {
            treecorr_multipole_process_pivots(
                context, pivot_tree, pivot_node->left,
                neighbor_tree, same_tree, scratch, hist);
            treecorr_multipole_process_pivots(
                context, pivot_tree, pivot_node->right,
                neighbor_tree, same_tree, scratch, hist);
        } else {
            for (INTEGER i = pivot_node->first; i <= pivot_node->last; i++)
                treecorr_multipole_process_body(
                    context, pivot_tree, i, neighbor_tree,
                    same_tree, scratch, hist);
        }
        return;
    }

    {
        treecorr_multipole_pivot pivot;

        pivot.position = pivot_node->center;
        pivot.radius = (real)pivot_node->radius;
        pivot.field_sum = treecorr_node_field_sum(context, pivot_node);
        pivot.normalization_sum =
            treecorr_node_normalization_sum(context, pivot_node);
        pivot.first = pivot_node->first;
        pivot.last = pivot_node->last;
        pivot.can_split = !treecorr_node_is_leaf(pivot_node);
        treecorr_multipole_clear(scratch);
        if (treecorr_multipole_scan_neighbors(
                context, &pivot, neighbor_tree, 0,
                same_tree, scratch) == SUCCESS) {
            treecorr_multipole_finish_pivot(
                hist, &pivot, scratch, context->cmd->sizeHistN);
            return;
        }
        scratch->pivot_restarts++;
    }

    if (!treecorr_node_is_leaf(pivot_node)) {
        treecorr_multipole_process_pivots(
            context, pivot_tree, pivot_node->left,
            neighbor_tree, same_tree, scratch, hist);
        treecorr_multipole_process_pivots(
            context, pivot_tree, pivot_node->right,
            neighbor_tree, same_tree, scratch, hist);
    } else {
        for (INTEGER i = pivot_node->first; i <= pivot_node->last; i++)
            treecorr_multipole_process_body(
                context, pivot_tree, i, neighbor_tree,
                same_tree, scratch, hist);
    }
}

static real treecorr_radial_upper_edge(
        const treecorr_search_context *context, int radial_bin)
{
    const struct cmdline_data *cmd = context->cmd;

    if (radial_bin >= cmd->sizeHistN) return cmd->rangeN;
    if (radial_bin <= 0) return cmd->rminHist;
    if (!cmd->useLogHist)
        return cmd->rminHist + (real)radial_bin * context->gd->deltaR;
    if (cmd->rminHist > 0.0)
        return cmd->rminHist
             * rpow(10.0, (real)radial_bin * context->gd->deltaR);
    return rpow(10.0,
        ((real)radial_bin - (real)cmd->sizeHistN) / cmd->logHistBinsPD
        + rlog10(cmd->rangeN));
}

static int treecorr_unresolved_radial_bin(
        const treecorr_search_context *context,
        real upper_extent, int radial_bin_limit)
{
    int radial_bin;

    if (!(upper_extent > context->cmd->rminHist)) return 0;
    if (upper_extent >= context->cmd->rangeN)
        return MIN(context->cmd->sizeHistN, radial_bin_limit - 1);
    radial_bin = treecorr_bin_index(context, upper_extent);
    if (radial_bin < 1) radial_bin = radial_bin_limit - 1;
    return MIN(radial_bin, radial_bin_limit - 1);
}

static void treecorr_multipole_mark_unresolved(
        const treecorr_search_context *context,
        real upper_extent, real *max_unresolved_extent)
{
    if (upper_extent > context->cmd->rminHist)
        *max_unresolved_extent = MAX(*max_unresolved_extent, upper_extent);
}

static void treecorr_multipole_clear_radial_through(
        treecorr_multipole_scratch *scratch, int radial_bin)
{
    const size_t count = (size_t)radial_bin + 1;
    real *order_arrays[5] = {
        scratch->field_cos, scratch->field_sin,
        scratch->self_coscos, scratch->self_sinsin,
        scratch->self_sincos
    };

    for (size_t array = 0; array < 5; array++)
        for (int order = 0; order < scratch->orders; order++)
            memset(order_arrays[array] + (size_t)order * scratch->stride,
                   0, count * sizeof(real));
    memset(scratch->normalization, 0, count * sizeof(real));
    memset(scratch->normalization_sq, 0, count * sizeof(real));
}

static void treecorr_multipole_copy_values(
        treecorr_multipole_scratch *destination,
        const treecorr_multipole_scratch *source)
{
    memcpy(destination->field_cos, source->field_cos,
           source->values * sizeof(real));
}

static void treecorr_multipole_scan_body_partial(
        const treecorr_search_context *context,
        const treecorr_multipole_pivot *pivot,
        treecorr_multipole_scratch *scratch,
        treecorr_multipole_scratch *statistics,
        bodyptr neighbor, real active_upper,
        real *max_unresolved_extent)
{
    real distance;
    real cosphi;
    real sinphi;
    real upper_extent;
    int radial_bin;

    if (pivot->radius > 0.0) {
        const int pair_status = treecorr_multipole_pair_status_limited(
            context, pivot, Pos(neighbor), 0.0, active_upper,
            &radial_bin, &cosphi, &sinphi, &upper_extent);

        statistics->pair_tests++;
        if (pair_status == TREECORR_TRIPLE_OUTSIDE) return;
        if (pair_status != TREECORR_TRIPLE_ACCEPT) {
            treecorr_multipole_mark_unresolved(
                context, upper_extent, max_unresolved_extent);
            return;
        }
    } else {
        if (!treecorr_polar_coordinates(
                context, pivot->position, Pos(neighbor),
                &distance, &cosphi, &sinphi)
            || !(distance < active_upper))
            return;
        radial_bin = treecorr_bin_index(context, distance);
        if (radial_bin < 0) return;
    }
    {
        const real field = treecorr_body_weighted_field(context, neighbor);
        const real normalization = treecorr_body_normalization(
            context, neighbor);

        treecorr_multipole_add(
            scratch, radial_bin, field, field * field,
            normalization, normalization * normalization,
            cosphi, sinphi);
    }
    statistics->body_visits++;
}

static void treecorr_multipole_scan_neighbors_partial(
        const treecorr_search_context *context,
        const treecorr_multipole_pivot *pivot,
        const fcfc_balltreeptr neighbor_tree, INTEGER neighbor_index,
        bool same_tree, treecorr_multipole_scratch *scratch,
        treecorr_multipole_scratch *statistics,
        real active_upper, real *max_unresolved_extent)
{
    const fcfc_ballnode *neighbor = &neighbor_tree->nodes[neighbor_index];
    int radial_bin;
    real cosphi;
    real sinphi;
    real upper_extent;
    int pair_status;

    if (same_tree && treecorr_ranges_overlap(
            pivot->first, pivot->last, neighbor->first, neighbor->last)) {
        if (neighbor->first >= pivot->first
            && neighbor->last <= pivot->last) {
            treecorr_multipole_mark_unresolved(
                context, MIN(active_upper, 2.0 * pivot->radius),
                max_unresolved_extent);
            return;
        }
        if (treecorr_node_is_leaf(neighbor)) {
            for (INTEGER i = neighbor->first; i <= neighbor->last; i++) {
                if (i >= pivot->first && i <= pivot->last) continue;
                treecorr_multipole_scan_body_partial(
                    context, pivot, scratch, statistics,
                    neighbor_tree->bptr[i], active_upper,
                    max_unresolved_extent);
            }
            return;
        }
        treecorr_multipole_scan_neighbors_partial(
            context, pivot, neighbor_tree, neighbor->left,
            same_tree, scratch, statistics,
            active_upper, max_unresolved_extent);
        treecorr_multipole_scan_neighbors_partial(
            context, pivot, neighbor_tree, neighbor->right,
            same_tree, scratch, statistics,
            active_upper, max_unresolved_extent);
        return;
    }

    statistics->pair_tests++;
    pair_status = treecorr_multipole_pair_status_limited(
        context, pivot, neighbor->center, (real)neighbor->radius,
        active_upper, &radial_bin, &cosphi, &sinphi, &upper_extent);
    if (pair_status == TREECORR_TRIPLE_OUTSIDE) return;
    if (pair_status == TREECORR_TRIPLE_ACCEPT) {
        treecorr_multipole_add(
            scratch, radial_bin,
            treecorr_node_field_sum(context, neighbor),
            treecorr_node_field_sq_sum(context, neighbor),
            treecorr_node_normalization_sum(context, neighbor),
            treecorr_node_normalization_sq_sum(context, neighbor),
            cosphi, sinphi);
        statistics->accepted_nodes++;
        return;
    }

    if (!treecorr_node_is_leaf(neighbor)
        && (!pivot->can_split
            || (real)neighbor->radius > pivot->radius)) {
        treecorr_multipole_scan_neighbors_partial(
            context, pivot, neighbor_tree, neighbor->left,
            same_tree, scratch, statistics,
            active_upper, max_unresolved_extent);
        treecorr_multipole_scan_neighbors_partial(
            context, pivot, neighbor_tree, neighbor->right,
            same_tree, scratch, statistics,
            active_upper, max_unresolved_extent);
        return;
    }
    if (pivot->can_split) {
        treecorr_multipole_mark_unresolved(
            context, upper_extent, max_unresolved_extent);
        return;
    }

    for (INTEGER i = neighbor->first; i <= neighbor->last; i++) {
        if (same_tree && i >= pivot->first && i <= pivot->last) continue;
        treecorr_multipole_scan_body_partial(
            context, pivot, scratch, statistics,
            neighbor_tree->bptr[i], active_upper,
            max_unresolved_extent);
    }
}

static int treecorr_balltree_depth(
        const fcfc_balltreeptr tree, INTEGER node_index)
{
    const fcfc_ballnode *ball_node = &tree->nodes[node_index];

    if (treecorr_node_is_leaf(ball_node)) return 1;
    return 1 + MAX(treecorr_balltree_depth(tree, ball_node->left),
                   treecorr_balltree_depth(tree, ball_node->right));
}

/*
 * Complete radial bins from large to small.  Before splitting a pivot, retain
 * its completed high-radius moments and clear only the unresolved low bins;
 * each child then inherits the completed work, as in TreeCorr's LogMultipole
 * traversal.
 */
static void treecorr_multipole_process_pivots_partial(
        const treecorr_search_context *context,
        const fcfc_balltreeptr pivot_tree, INTEGER pivot_node_index,
        const fcfc_balltreeptr neighbor_tree, bool same_tree,
        treecorr_multipole_scratch *scratch,
        treecorr_multipole_scratch *statistics,
        real *scratch_levels, size_t scratch_values_per_level,
        int depth,
        real active_upper, int radial_bin_limit,
        treecorr_triple_histogram *hist)
{
    const fcfc_ballnode *pivot_node = &pivot_tree->nodes[pivot_node_index];
    treecorr_multipole_pivot pivot;
    real max_unresolved_extent = 0.0;
    int unresolved_bin;

    pivot.position = pivot_node->center;
    pivot.radius = (real)pivot_node->radius;
    pivot.field_sum = treecorr_node_field_sum(context, pivot_node);
    pivot.normalization_sum =
        treecorr_node_normalization_sum(context, pivot_node);
    pivot.first = pivot_node->first;
    pivot.last = pivot_node->last;
    pivot.can_split = !treecorr_node_is_leaf(pivot_node);

    treecorr_multipole_scan_neighbors_partial(
        context, &pivot, neighbor_tree, 0, same_tree,
        scratch, statistics, active_upper, &max_unresolved_extent);
    unresolved_bin = treecorr_unresolved_radial_bin(
        context, max_unresolved_extent, radial_bin_limit);

    if (unresolved_bin + 1 < radial_bin_limit) {
        treecorr_multipole_finish_pivot_range(
            hist, &pivot, scratch, unresolved_bin + 1,
            radial_bin_limit, context->cmd->sizeHistN);
        statistics->pivot_finishes++;
    }
    if (unresolved_bin == 0) return;

    statistics->pivot_restarts++;
    treecorr_multipole_clear_radial_through(scratch, unresolved_bin);
    active_upper = treecorr_radial_upper_edge(context, unresolved_bin);
    radial_bin_limit = unresolved_bin + 1;

    if (!treecorr_node_is_leaf(pivot_node)) {
        const INTEGER children[2] = {pivot_node->left, pivot_node->right};

        for (int child_index = 0; child_index < 2; child_index++) {
            treecorr_multipole_scratch child_scratch;
            real *child_base = scratch_levels
                + (size_t)(depth + 1) * scratch_values_per_level;

            treecorr_initialize_multipole_scratch(
                &child_scratch, child_base, scratch->stride,
                scratch->orders, scratch_values_per_level);
            treecorr_multipole_copy_values(&child_scratch, scratch);
            treecorr_multipole_process_pivots_partial(
                context, pivot_tree, children[child_index],
                neighbor_tree, same_tree, &child_scratch, statistics,
                scratch_levels, scratch_values_per_level,
                depth + 1,
                active_upper, radial_bin_limit, hist);
        }
    } else {
        for (INTEGER i = pivot_node->first; i <= pivot_node->last; i++) {
            treecorr_multipole_scratch body_scratch;
            treecorr_multipole_pivot body_pivot;
            bodyptr pivot_body = pivot_tree->bptr[i];
            real body_unresolved_extent = 0.0;
            real *body_base = scratch_levels
                + (size_t)(depth + 1) * scratch_values_per_level;

            treecorr_initialize_multipole_scratch(
                &body_scratch, body_base, scratch->stride,
                scratch->orders, scratch_values_per_level);
            treecorr_multipole_copy_values(&body_scratch, scratch);
            body_pivot.position = Pos(pivot_body);
            body_pivot.radius = 0.0;
            body_pivot.field_sum =
                treecorr_body_weighted_field(context, pivot_body);
            body_pivot.normalization_sum =
                treecorr_body_normalization(context, pivot_body);
            body_pivot.first = i;
            body_pivot.last = i;
            body_pivot.can_split = FALSE;
            treecorr_multipole_scan_neighbors_partial(
                context, &body_pivot, neighbor_tree, 0, same_tree,
                &body_scratch, statistics, active_upper,
                &body_unresolved_extent);
            treecorr_multipole_finish_pivot_range(
                hist, &body_pivot, &body_scratch,
                1, radial_bin_limit, context->cmd->sizeHistN);
            statistics->pivot_finishes++;
        }
    }
}

static int treecorr_allocate_multipole_scratch(
        struct cmdline_data *cmd, INTEGER task_count,
        size_t stride, int orders, int levels, real **storage,
        size_t *values_per_level, size_t *values_per_task)
{
    const size_t order_arrays = 5;
    size_t order_values;
    size_t values;
    size_t tasks;

    *storage = NULL;
    *values_per_level = 0;
    *values_per_task = 0;
    if (task_count <= 0 || orders <= 0 || levels <= 0
        || (size_t)orders > SIZE_MAX / order_arrays
        || (size_t)orders * order_arrays > SIZE_MAX / stride)
        goto invalid_size;
    order_values = order_arrays * (size_t)orders * stride;
    if (stride > (SIZE_MAX - order_values) / 2) goto invalid_size;
    *values_per_level = order_values + 2 * stride;
    if ((size_t)levels > SIZE_MAX / *values_per_level)
        goto invalid_size;
    *values_per_task = (size_t)levels * *values_per_level;
    tasks = (size_t)task_count;
    if (tasks > SIZE_MAX / *values_per_task
        || (values = tasks * *values_per_task) > SIZE_MAX / sizeof(**storage))
        goto invalid_size;
    *storage = calloc(values, sizeof(**storage));
    if (*storage == NULL) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: multipole scratch allocation failed", TREECORR_METHOD_NAME);
        return FAILURE;
    }
    return SUCCESS;

invalid_size:
    snprintf(cmd->error_message, _ERRORMSGSIZE_,
             "%s: multipole scratch size overflow", TREECORR_METHOD_NAME);
    return FAILURE;
}

static void treecorr_initialize_multipole_scratch(
        treecorr_multipole_scratch *scratch, real *base,
        size_t stride, int orders, size_t values)
{
    const size_t order_values = (size_t)orders * stride;

    scratch->field_cos = base;
    scratch->field_sin = scratch->field_cos + order_values;
    scratch->self_coscos = scratch->field_sin + order_values;
    scratch->self_sinsin = scratch->self_coscos + order_values;
    scratch->self_sincos = scratch->self_sinsin + order_values;
    scratch->normalization = scratch->self_sincos + order_values;
    scratch->normalization_sq = scratch->normalization + stride;
    scratch->stride = stride;
    scratch->values = values;
    scratch->orders = orders;
    scratch->body_visits = 0;
    scratch->accepted_nodes = 0;
    scratch->pair_tests = 0;
    scratch->pivot_restarts = 0;
    scratch->pivot_finishes = 0;
}

static int treecorr_search_log_multipole(
        struct cmdline_data *cmd, struct global_data *gd,
        bodyptr *btab, INTEGER *nbody, INTEGER ipmin, INTEGER *ipmax,
        int cat1, int cat2)
{
    const bool auto_correlation = cat1 == cat2;
    const bool only_2pcf = scanopt(cmd->options, "only-2pcf");
    const bool only_3pcf = scanopt(cmd->options, "only-3pcf");
#ifdef TWOPCF
    const bool run_2pcf = !only_3pcf;
#else
    const bool run_2pcf = FALSE;
#endif
    const bool run_3pcf = !only_2pcf;
    const size_t stride = (size_t)cmd->sizeHistN + 1;
    const int orders = cmd->mChebyshev + 1;
    const INTEGER target_tasks = run_3pcf
        ? treecorr_task_target(stride, orders)
        : TREECORR_PAIR_FRONTIER_TARGET;
    const int leaf_capacity = scanopt(cmd->options, "treecorr-singleton-leaves")
        ? 1 : cmd->nsmooth;
    treecorr_search_context context;
    fcfc_balltreeptr tree1 = NULL;
    fcfc_balltreeptr tree2 = NULL;
    INTEGER *frontier = NULL;
    INTEGER *pair_frontier2 = NULL;
    INTEGER task_count = 0;
    INTEGER pair_frontier_count2 = 0;
    real *task_histograms = NULL;
    real *task_scratch = NULL;
    INTEGER *task_body_counts = NULL;
    INTEGER *task_cell_counts = NULL;
    size_t hist_values_per_task = 0;
    size_t scratch_values_per_level = 0;
    size_t scratch_values_per_task = 0;
    int scratch_level_count = 0;
    INTEGER body_total = 0;
    INTEGER cell_total = 0;
    INTEGER pair_test_total = 0;
    INTEGER pivot_restart_total = 0;
    INTEGER pivot_finish_total = 0;
    INTEGER distributed_statistics[3] = {0, 0, 0};
    int operation_status;
    int reduction_status = SUCCESS;
    int status = FAILURE;
    const double cpustart = CPUTIME;

    if (only_2pcf && only_3pcf) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: only-2pcf and only-3pcf are mutually exclusive",
                 TREECORR_METHOD_NAME);
        return FAILURE;
    }
    if (only_2pcf && !run_2pcf) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: only-2pcf requires TWOPCFON=1",
                 TREECORR_METHOD_NAME);
        return FAILURE;
    }
    if (ipmin != 1 || ipmax[cat1] != nbody[cat1]) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s requires the complete pivot catalog",
                 TREECORR_METHOD_NAME);
        return FAILURE;
    }
    if (cmd->nsmooth <= 0) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s requires nsmooth > 0", TREECORR_METHOD_NAME);
        return FAILURE;
    }

    verb_print(cmd->verbose, "Search: Running %s", cmd->searchMethod);
#ifdef TWOPCF
    if (run_2pcf) verb_print(cmd->verbose, " with dual-node 2PCF");
#endif
    if (run_3pcf)
        verb_print(cmd->verbose, " with LogMultipole 3PCF");
    verb_print(cmd->verbose, "\n");
#ifdef TREECORR_DISTRIBUTED_ENGINE
    verb_print(cmd->verbose,
               "%s: %d ranks with deterministic cyclic frontier ownership\n",
               TREECORR_METHOD_NAME, TREECORR_DISTRIBUTED_SIZE());
#endif
    if (cballs_opt_no_two_balls(cmd))
        verb_print(cmd->verbose,
                   "no-two-balls: exact body pivot-neighbor accumulation\n");
    else
        verb_print(cmd->verbose,
                   "two-ball radial/angular node acceptance enabled; theta=%g\n",
                   cmd->theta);
#ifdef TREECORR_NATIVE_BINARY_VIEW
    verb_print(cmd->verbose,
               "native-octree binary view with exact body pivots\n");
#else
    verb_print(cmd->verbose,
               "ball-tree leaf capacity = %d\n", leaf_capacity);
#endif
    if (run_3pcf)
        verb_print(cmd->verbose, "3PCF multipoles: %s\n",
                   treecorr_normalize_3pcf(cmd)
                   ? (cballs_opt_weights_norm(cmd)
                      ? "normalized by distinct-triplet weight sum"
                      : "normalized by distinct-triplet count")
                   : "raw distinct-triplet sums");

#ifdef OPENMPCODE
    ThreadCount(cmd, gd, nbody[cat1], cat1);
#endif

    operation_status = search_init_gd_hist_sincos(cmd, gd);
    if (treecorr_distributed_consensus(
            cmd, operation_status,
            "two-ball histogram initialization") == FAILURE)
        goto cleanup;
    operation_status = fcfc_balltree_build(
        cmd, gd, btab[cat1], nbody[cat1], leaf_capacity, &tree1);
    if (treecorr_distributed_consensus(
            cmd, operation_status,
            "two-ball first tree construction") == FAILURE)
        goto cleanup;
    TREECORR_PUBLISH_NODE_COUNT(gd, cat1, tree1->nnode);
    if (auto_correlation) {
        tree2 = tree1;
    } else {
        operation_status = fcfc_balltree_build(
            cmd, gd, btab[cat2], nbody[cat2], leaf_capacity, &tree2);
        if (treecorr_distributed_consensus(
                cmd, operation_status,
                "two-ball second tree construction") == FAILURE)
            goto cleanup;
        TREECORR_PUBLISH_NODE_COUNT(gd, cat2, tree2->nnode);
    }

    context.cmd = cmd;
    context.gd = gd;
    context.use_two_balls = !cballs_opt_no_two_balls(cmd);
    context.use_three_cells = context.use_two_balls;
    context.weighted = cballs_opt_weights_norm(cmd);
    treecorr_initialize_angular_tolerance(&context);

    operation_status = fcfc_balltree_frontier(
        cmd, tree1, target_tasks, &frontier, &task_count);
    if (treecorr_distributed_consensus(
            cmd, operation_status,
            "two-ball pivot-frontier construction") == FAILURE)
        goto cleanup;
    if (auto_correlation) {
        pair_frontier2 = frontier;
        pair_frontier_count2 = task_count;
#ifdef TWOPCF
    } else if (run_2pcf) {
        operation_status = fcfc_balltree_frontier(
            cmd, tree2, target_tasks,
            &pair_frontier2, &pair_frontier_count2);
        if (treecorr_distributed_consensus(
                cmd, operation_status,
                "two-ball neighbor-frontier construction") == FAILURE)
            goto cleanup;
#endif
    }

#ifdef TWOPCF
    if (run_2pcf
        && treecorr_run_pair_tasks(
            &context, tree1, tree2, auto_correlation,
            frontier, task_count,
            pair_frontier2, pair_frontier_count2) == FAILURE)
        goto cleanup;
#endif

    if (run_3pcf) {
#ifdef TREECORR_BODY_PIVOT_LOG_MULTIPOLE
    scratch_level_count = 1;
#else
    scratch_level_count = treecorr_balltree_depth(tree1, 0) + 1;
#endif
    operation_status = treecorr_allocate_triple_histograms(
        cmd, task_count, stride, orders,
        &task_histograms, &hist_values_per_task,
        &task_body_counts, &task_cell_counts);
    if (treecorr_distributed_consensus(
            cmd, operation_status,
            "two-ball multipole histogram allocation") == FAILURE)
        goto cleanup;
    operation_status = treecorr_allocate_multipole_scratch(
        cmd, task_count, stride, orders, scratch_level_count,
        &task_scratch, &scratch_values_per_level,
        &scratch_values_per_task);
    if (treecorr_distributed_consensus(
            cmd, operation_status,
            "two-ball multipole scratch allocation") == FAILURE)
        goto cleanup;

#pragma omp parallel for schedule(dynamic,1) \
    reduction(+:pair_test_total,pivot_restart_total,pivot_finish_total)
    for (INTEGER itask = 0; itask < task_count; itask++) {
        real *hist_base = task_histograms
                        + (size_t)itask * hist_values_per_task;
        real *scratch_base = task_scratch
                           + (size_t)itask * scratch_values_per_task;
        treecorr_triple_histogram hist;
        treecorr_multipole_scratch scratch;

        if (!treecorr_distributed_task_owned(itask)) continue;

        hist.components = hist_base;
        hist.stride = stride;
        hist.plane = stride * stride;
        hist.order_plane = (size_t)orders * hist.plane;
        hist.orders = orders;
        hist.normalization = hist_base
            + TREECORR_ZETA_COMPONENTS * hist.order_plane;
        hist.body_triples = 0;
        hist.cell_triples = 0;
        treecorr_initialize_multipole_scratch(
            &scratch, scratch_base, stride, orders,
            scratch_values_per_level);
#ifdef TREECORR_BODY_PIVOT_LOG_MULTIPOLE
        treecorr_multipole_process_body_pivots(
            &context, tree1, frontier[itask], tree2,
            auto_correlation, &scratch, &hist);
#else
        if (context.use_two_balls) {
            treecorr_multipole_clear(&scratch);
            treecorr_multipole_process_pivots_partial(
                &context, tree1, frontier[itask], tree2,
                auto_correlation, &scratch, &scratch,
                scratch_base, scratch_values_per_level,
                0, cmd->rangeN,
                cmd->sizeHistN + 1, &hist);
        } else {
            treecorr_multipole_process_pivots(
                &context, tree1, frontier[itask], tree2,
                auto_correlation, &scratch, &hist);
        }
#endif
        pair_test_total += scratch.pair_tests;
        pivot_restart_total += scratch.pivot_restarts;
        pivot_finish_total += scratch.pivot_finishes;
        task_body_counts[itask] = scratch.body_visits;
        task_cell_counts[itask] = scratch.accepted_nodes;
    }

    distributed_statistics[0] = pair_test_total;
    distributed_statistics[1] = pivot_restart_total;
    distributed_statistics[2] = pivot_finish_total;
    if (treecorr_distributed_reduce_reals(
            cmd, task_histograms,
            (size_t)task_count * hist_values_per_task) == FAILURE)
        reduction_status = FAILURE;
    if (treecorr_distributed_reduce_integers(
            cmd, task_body_counts, (size_t)task_count) == FAILURE)
        reduction_status = FAILURE;
    if (treecorr_distributed_reduce_integers(
            cmd, task_cell_counts, (size_t)task_count) == FAILURE)
        reduction_status = FAILURE;
    if (treecorr_distributed_reduce_integers(
            cmd, distributed_statistics, 3) == FAILURE)
        reduction_status = FAILURE;
    if (treecorr_distributed_consensus(
            cmd, reduction_status,
            "two-ball multipole reduction") == FAILURE)
        goto cleanup;

    if (treecorr_distributed_publish()) {
    pair_test_total = distributed_statistics[0];
    pivot_restart_total = distributed_statistics[1];
    pivot_finish_total = distributed_statistics[2];
    for (INTEGER itask = 0; itask < task_count; itask++) {
        const real *base = task_histograms
                         + (size_t)itask * hist_values_per_task;
        const size_t plane = stride * stride;

        for (int order = 0; order < orders; order++) {
            const real *zcos = base
                + ((size_t)TREECORR_ZETA_COS * (size_t)orders
                   + (size_t)order) * plane;
            const real *zsin = base
                + ((size_t)TREECORR_ZETA_SIN * (size_t)orders
                   + (size_t)order) * plane;
            const real *zsincos = base
                + ((size_t)TREECORR_ZETA_SINCOS * (size_t)orders
                   + (size_t)order) * plane;
            const real *zcossin = base
                + ((size_t)TREECORR_ZETA_COSSIN * (size_t)orders
                   + (size_t)order) * plane;

            for (int n1 = 1; n1 <= cmd->sizeHistN; n1++)
                for (int n2 = 1; n2 <= cmd->sizeHistN; n2++) {
                    const size_t index = (size_t)n1 * stride + (size_t)n2;
                    gd->histZetaMcos[order + 1][n1][n2] += zcos[index];
                    gd->histZetaMsin[order + 1][n1][n2] += zsin[index];
                    gd->histZetaMsincos[order + 1][n1][n2] += zsincos[index];
                    gd->histZetaMcossin[order + 1][n1][n2] += zcossin[index];
                }
        }
        body_total += task_body_counts[itask];
        cell_total += task_cell_counts[itask];
    }

    if (treecorr_normalize_3pcf(cmd)) {
        for (int n1 = 1; n1 <= cmd->sizeHistN; n1++)
            for (int n2 = 1; n2 <= cmd->sizeHistN; n2++) {
            const size_t index = (size_t)n1 * stride + (size_t)n2;
            real denominator = 0.0;

            for (INTEGER itask = 0; itask < task_count; itask++) {
                const real *base = task_histograms
                    + (size_t)itask * hist_values_per_task;
                denominator += base[
                    TREECORR_ZETA_COMPONENTS
                    * (size_t)orders * stride * stride + index];
            }
            for (int order = 1; order <= orders; order++) {
                gd->histZetaMcos[order][n1][n2] = cballs_normalize_or_zero(
                    gd->histZetaMcos[order][n1][n2], denominator);
                gd->histZetaMsin[order][n1][n2] = cballs_normalize_or_zero(
                    gd->histZetaMsin[order][n1][n2], denominator);
                gd->histZetaMsincos[order][n1][n2] = cballs_normalize_or_zero(
                    gd->histZetaMsincos[order][n1][n2], denominator);
                gd->histZetaMcossin[order][n1][n2] = cballs_normalize_or_zero(
                    gd->histZetaMcossin[order][n1][n2], denominator);
            }
        }
    }
    }
    }

    gd->cpusearch = CPUTIME - cpustart;
#ifdef TWOPCF
    if (run_2pcf)
        verb_print(cmd->verbose,
                   "%s: nbbcalc = %" INTEGER_FMT
                   ", nbccalc = %" INTEGER_FMT
                   ", frontier tasks = %" INTEGER_FMT "\n",
                   cmd->searchMethod, gd->nbbcalc, gd->nbccalc, task_count);
#endif
    if (run_3pcf) {
        verb_print(cmd->verbose,
                   "%s: neighbor body visits = %" INTEGER_FMT
                   ", accepted cell orientations = %" INTEGER_FMT
                   ", pivot tasks = %" INTEGER_FMT "\n",
                   cmd->searchMethod, body_total, cell_total, task_count);
        verb_print(cmd->verbose,
                   "%s: pair tests = %" INTEGER_FMT
                   ", restarted pivot scans = %" INTEGER_FMT
                   ", finished pivots = %" INTEGER_FMT "\n",
                   cmd->searchMethod, pair_test_total,
                   pivot_restart_total, pivot_finish_total);
    }
    verb_print(cmd->verbose, "Going out: CPU time = %lf\n", gd->cpusearch);
    status = SUCCESS;

cleanup:
    free(task_scratch);
    free(task_cell_counts);
    free(task_body_counts);
    free(task_histograms);
    if (!auto_correlation) free(pair_frontier2);
    free(frontier);
    if (tree2 != tree1) fcfc_balltree_free(tree2);
    fcfc_balltree_free(tree1);
    return status;
}
#endif /* THREEPCFCONVERGENCE */

global int searchcalc_balltree_2balls_omp(
        struct cmdline_data *cmd, struct global_data *gd,
        bodyptr *btab, INTEGER *nbody, INTEGER ipmin, INTEGER *ipmax,
        int cat1, int cat2)
{
#ifdef TREECORR_DISTRIBUTED_ENGINE
    if (scanopt(cmd->options, "treecorr-direct-triples")) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s does not distribute treecorr-direct-triples; "
                 "use %s without that validation option",
                 TREECORR_METHOD_NAME, TREECORR_METHOD_NAME);
        return FAILURE;
    }
#else
    if (scanopt(cmd->options, "treecorr-direct-triples"))
        return treecorr_search_direct_triples(
            cmd, gd, btab, nbody, ipmin, ipmax, cat1, cat2);
#endif
#ifdef THREEPCFCONVERGENCE
    return treecorr_search_log_multipole(
        cmd, gd, btab, nbody, ipmin, ipmax, cat1, cat2);
#else
    return treecorr_search_direct_triples(
        cmd, gd, btab, nbody, ipmin, ipmax, cat1, cat2);
#endif
}
#endif /* TREECORR_LOG_MULTIPOLE_ENGINE */
#endif /* TREECORR_TASK_FRONTIER_ENGINE */

#endif /* THREEPCFCONVERGENCE */

#ifdef TWOPCF
static int treecorr_allocate_task_histograms(
        struct cmdline_data *cmd, INTEGER task_count, size_t stride,
        real **histograms, INTEGER **body_pairs, INTEGER **cell_pairs)
{
    size_t tasks;
    size_t values;

    *histograms = NULL;
    *body_pairs = NULL;
    *cell_pairs = NULL;
    if (task_count <= 0) goto invalid_size;
    tasks = (size_t)task_count;
    if (tasks > SIZE_MAX / (3 * stride)
        || (values = tasks * 3 * stride) > SIZE_MAX / sizeof(**histograms)
        || tasks > SIZE_MAX / sizeof(**body_pairs))
        goto invalid_size;

    *histograms = calloc(values, sizeof(**histograms));
    *body_pairs = calloc(tasks, sizeof(**body_pairs));
    *cell_pairs = calloc(tasks, sizeof(**cell_pairs));
    if (*histograms == NULL || *body_pairs == NULL || *cell_pairs == NULL) {
        free(*cell_pairs);
        free(*body_pairs);
        free(*histograms);
        *histograms = NULL;
        *body_pairs = NULL;
        *cell_pairs = NULL;
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: task histogram allocation failed", TREECORR_METHOD_NAME);
        return FAILURE;
    }
    return SUCCESS;

invalid_size:
    snprintf(cmd->error_message, _ERRORMSGSIZE_,
             "%s: task histogram size overflow", TREECORR_METHOD_NAME);
    return FAILURE;
}

static int treecorr_run_pair_tasks(
        const treecorr_search_context *context,
        fcfc_balltreeptr tree1, fcfc_balltreeptr tree2,
        bool auto_correlation, const INTEGER *frontier1,
        INTEGER frontier_count1, const INTEGER *frontier2,
        INTEGER frontier_count2)
{
    struct cmdline_data *cmd = context->cmd;
    struct global_data *gd = context->gd;
    const size_t stride = (size_t)cmd->sizeHistN + 1;
    real *task_histograms = NULL;
    INTEGER *task_body_counts = NULL;
    INTEGER *task_cell_counts = NULL;
    int allocation_status;
    int reduction_status = SUCCESS;
    int status = FAILURE;

    allocation_status = treecorr_allocate_task_histograms(
        cmd, frontier_count1, stride, &task_histograms,
        &task_body_counts, &task_cell_counts);
    if (treecorr_distributed_consensus(
            cmd, allocation_status,
            "two-ball pair-task allocation") == FAILURE)
        goto cleanup;

#pragma omp parallel for schedule(dynamic,1)
    for (INTEGER itask = 0; itask < frontier_count1; itask++) {
        real *base = task_histograms + (size_t)itask * 3 * stride;
        treecorr_histogram hist;

        if (!treecorr_distributed_task_owned(itask)) continue;

        hist.pair_count = base;
        hist.weight_product = base + stride;
        hist.field_product = base + 2 * stride;
        hist.body_pairs = 0;
        hist.cell_pairs = 0;

        if (auto_correlation) {
            treecorr_process_auto(
                context, tree1, frontier1[itask], &hist);
            for (INTEGER jtask = itask + 1;
                 jtask < frontier_count2; jtask++)
                treecorr_process_pair(
                    context, tree1, frontier1[itask],
                    tree2, frontier2[jtask], &hist);
        } else {
            for (INTEGER jtask = 0; jtask < frontier_count2; jtask++)
                treecorr_process_pair(
                    context, tree1, frontier1[itask],
                    tree2, frontier2[jtask], &hist);
        }
        task_body_counts[itask] = hist.body_pairs;
        task_cell_counts[itask] = hist.cell_pairs;
    }

    if (treecorr_distributed_reduce_reals(
            cmd, task_histograms,
            (size_t)frontier_count1 * 3 * stride) == FAILURE)
        reduction_status = FAILURE;
    if (treecorr_distributed_reduce_integers(
            cmd, task_body_counts, (size_t)frontier_count1) == FAILURE)
        reduction_status = FAILURE;
    if (treecorr_distributed_reduce_integers(
            cmd, task_cell_counts, (size_t)frontier_count1) == FAILURE)
        reduction_status = FAILURE;
    if (treecorr_distributed_consensus(
            cmd, reduction_status,
            "two-ball pair-task reduction") == FAILURE)
        goto cleanup;
    if (!treecorr_distributed_publish()) {
        status = SUCCESS;
        goto cleanup;
    }

    for (INTEGER itask = 0; itask < frontier_count1; itask++) {
        const real *base = task_histograms + (size_t)itask * 3 * stride;
        const real *pair_count = base;
        const real *field_product = base + 2 * stride;

        for (int n = 1; n <= cmd->sizeHistN; n++) {
            gd->histNN[n] += pair_count[n];
            gd->histNNSubXi2pcf[n] += pair_count[n];
            gd->histXi2pcf[n] += field_product[n];
        }
        gd->nbbcalc += task_body_counts[itask];
        gd->nbccalc += task_cell_counts[itask];
    }

    for (int n = 1; n <= cmd->sizeHistN; n++) {
        real denominator = 0.0;

        for (INTEGER itask = 0; itask < frontier_count1; itask++) {
            const real *base = task_histograms
                + (size_t)itask * 3 * stride;
            denominator += base[stride + (size_t)n];
        }
        gd->histXi2pcf[n] = cballs_normalize_or_zero(
            gd->histXi2pcf[n], denominator);
    }

    if (cballs_opt_compute_histn(cmd) && cballs_opt_and_cf(cmd)) {
        for (int n = 1; n <= cmd->sizeHistN; n++) gd->histNN[n] *= 2.0;
#ifdef LONGINT
        if (tree1->npoint > INT_MAX) {
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     "%s: and-CF body count exceeds int",
                     TREECORR_METHOD_NAME);
            goto cleanup;
        }
#endif
        if (search_compute_HistN(cmd, gd, (int)tree1->npoint) == FAILURE)
            goto cleanup;
    }
    status = SUCCESS;

cleanup:
    status = treecorr_distributed_consensus(
        cmd, status, "two-ball pair publication");
    free(task_cell_counts);
    free(task_body_counts);
    free(task_histograms);
    return status;
}
#endif /* TWOPCF */

#ifndef TREECORR_TASK_FRONTIER_ENGINE
global int searchcalc_balltree_2balls_omp(
        struct cmdline_data *cmd, struct global_data *gd,
        bodyptr *btab, INTEGER *nbody, INTEGER ipmin, INTEGER *ipmax,
        int cat1, int cat2)
{
    const bool auto_correlation = cat1 == cat2;
    const bool only_2pcf = scanopt(cmd->options, "only-2pcf");
    const bool only_3pcf = scanopt(cmd->options, "only-3pcf");
#ifdef TWOPCF
    const bool run_2pcf = !only_3pcf;
#else
    const bool run_2pcf = FALSE;
#endif
#ifdef THREEPCFCONVERGENCE
    const bool run_3pcf = !only_2pcf;
#else
    const bool run_3pcf = FALSE;
#endif
    const size_t stride = (size_t)cmd->sizeHistN + 1;
    treecorr_search_context context;
    fcfc_balltreeptr tree1 = NULL;
    fcfc_balltreeptr tree2 = NULL;
    INTEGER *frontier1 = NULL;
    INTEGER *frontier2 = NULL;
    INTEGER frontier_count1 = 0;
    INTEGER frontier_count2 = 0;
#ifdef THREEPCFCONVERGENCE
    real *triple_task_histograms = NULL;
    INTEGER *triple_task_body_counts = NULL;
    INTEGER *triple_task_cell_counts = NULL;
    size_t triple_values_per_task = 0;
    const int triple_orders = cmd->mChebyshev + 1;
    INTEGER triple_body_total = 0;
    INTEGER triple_cell_total = 0;
    int triple_reduction_status = SUCCESS;
#endif
    const INTEGER target_tasks = run_3pcf
        ? TREECORR_TRIPLE_FRONTIER_TARGET : TREECORR_PAIR_FRONTIER_TARGET;
    int operation_status;
    int status = FAILURE;
    const double cpustart = CPUTIME;

#if !defined(TWOPCF) && !defined(THREEPCFCONVERGENCE)
    snprintf(cmd->error_message, _ERRORMSGSIZE_,
             "%s was built with TWOPCFON=0 and TPCFON=0",
             TREECORR_METHOD_NAME);
    return FAILURE;
#endif

    if (only_2pcf && only_3pcf) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: only-2pcf and only-3pcf are mutually exclusive",
                 TREECORR_METHOD_NAME);
        return FAILURE;
    }
    if (only_2pcf && !run_2pcf) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: only-2pcf requires TWOPCFON=1", TREECORR_METHOD_NAME);
        return FAILURE;
    }
    if (only_3pcf && !run_3pcf) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: only-3pcf requires TPCFON=1", TREECORR_METHOD_NAME);
        return FAILURE;
    }

    if (ipmin != 1 || ipmax[cat1] != nbody[cat1]) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s requires the complete pivot catalog",
                 TREECORR_METHOD_NAME);
        return FAILURE;
    }
    if (cmd->nsmooth <= 0) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s requires nsmooth > 0", TREECORR_METHOD_NAME);
        return FAILURE;
    }

    verb_print(cmd->verbose, "Search: Running %s", cmd->searchMethod);
#ifdef TWOPCF
    if (run_2pcf) verb_print(cmd->verbose, " with dual-node 2PCF");
#endif
#ifdef THREEPCFCONVERGENCE
    if (run_3pcf) verb_print(cmd->verbose, " with triple-node 3PCF");
#endif
    verb_print(cmd->verbose, "\n");
#ifdef TREECORR_DISTRIBUTED_ENGINE
    verb_print(cmd->verbose,
               "%s: %d ranks with deterministic cyclic frontier ownership\n",
               TREECORR_METHOD_NAME, TREECORR_DISTRIBUTED_SIZE());
#endif
#ifdef THREEPCFCONVERGENCE
    if (run_3pcf && nbody[cat1] >= TREECORR_LARGE_DIRECT_3PCF)
        verb_print(cmd->verbose,
                   "warning: direct triple-node 3PCF on %" INTEGER_FMT
                   " pivots can be extremely expensive; use only-2pcf or "
                   "reduce the catalog when 3PCF is not required\n",
                   nbody[cat1]);
#endif
    if (cballs_opt_no_two_balls(cmd))
        verb_print(cmd->verbose,
                   "no-two-balls: forcing exact body-pair and body-triplet accumulation\n");
    else
        verb_print(cmd->verbose,
                   "two/three-cell aggregation enabled; theta=%g\n", cmd->theta);
    if (cballs_opt_weights_norm(cmd))
        verb_print(cmd->verbose,
                   "using TreeCorr-style pair/triplet-weight normalization\n");

#ifdef OPENMPCODE
    ThreadCount(cmd, gd, nbody[cat1], cat1);
#endif

    operation_status = search_init_gd_hist_sincos(cmd, gd);
    if (treecorr_distributed_consensus(
            cmd, operation_status,
            "two-ball histogram initialization") == FAILURE)
        goto cleanup;
    operation_status = fcfc_balltree_build(
        cmd, gd, btab[cat1], nbody[cat1], cmd->nsmooth, &tree1);
    if (treecorr_distributed_consensus(
            cmd, operation_status,
            "two-ball first tree construction") == FAILURE)
        goto cleanup;
    TREECORR_PUBLISH_NODE_COUNT(gd, cat1, tree1->nnode);

    if (auto_correlation) {
        tree2 = tree1;
    } else {
        operation_status = fcfc_balltree_build(
            cmd, gd, btab[cat2], nbody[cat2], cmd->nsmooth, &tree2);
        if (treecorr_distributed_consensus(
                cmd, operation_status,
                "two-ball second tree construction") == FAILURE)
            goto cleanup;
        TREECORR_PUBLISH_NODE_COUNT(gd, cat2, tree2->nnode);
    }

    operation_status = fcfc_balltree_frontier(
        cmd, tree1, target_tasks, &frontier1, &frontier_count1);
    if (treecorr_distributed_consensus(
            cmd, operation_status,
            "two-ball pivot-frontier construction") == FAILURE)
        goto cleanup;
    if (auto_correlation) {
        frontier2 = frontier1;
        frontier_count2 = frontier_count1;
    } else {
        operation_status = fcfc_balltree_frontier(
            cmd, tree2, target_tasks, &frontier2, &frontier_count2);
        if (treecorr_distributed_consensus(
                cmd, operation_status,
                "two-ball neighbor-frontier construction") == FAILURE)
            goto cleanup;
    }

    context.cmd = cmd;
    context.gd = gd;
    context.use_two_balls = !cballs_opt_no_two_balls(cmd);
    context.use_three_cells = !cballs_opt_no_two_balls(cmd);
    context.weighted = cballs_opt_weights_norm(cmd);
#ifdef THREEPCFCONVERGENCE
    treecorr_initialize_angular_tolerance(&context);
#endif

#ifdef TWOPCF
    if (run_2pcf) {
        if (treecorr_run_pair_tasks(
                &context, tree1, tree2, auto_correlation,
                frontier1, frontier_count1,
                frontier2, frontier_count2) == FAILURE)
            goto cleanup;
    }
#endif /* TWOPCF */

#ifdef THREEPCFCONVERGENCE
    if (run_3pcf) {
    operation_status = treecorr_allocate_triple_histograms(
        cmd, frontier_count1, stride, triple_orders,
        &triple_task_histograms, &triple_values_per_task,
        &triple_task_body_counts, &triple_task_cell_counts);
    if (treecorr_distributed_consensus(
            cmd, operation_status,
            "two-ball triple histogram allocation") == FAILURE)
        goto cleanup;

#pragma omp parallel for schedule(dynamic,1)
    for (INTEGER itask = 0; itask < frontier_count1; itask++) {
        real *base = triple_task_histograms
                   + (size_t)itask * triple_values_per_task;
        treecorr_triple_histogram hist;

        hist.components = base;
        hist.stride = stride;
        hist.plane = stride * stride;
        hist.order_plane = (size_t)triple_orders * hist.plane;
        hist.orders = triple_orders;
        hist.normalization = base
            + TREECORR_ZETA_COMPONENTS * hist.order_plane;
        hist.body_triples = 0;
        hist.cell_triples = 0;

        if (!treecorr_distributed_task_owned(itask)) continue;

        if (auto_correlation) {
            const fcfc_balltreeptr trees[3] = {tree1, tree1, tree1};

            treecorr_process3_auto(&context, tree1, frontier1[itask], &hist);
            for (INTEGER jtask = itask + 1;
                 jtask < frontier_count1; jtask++) {
                treecorr_process21_auto(&context, tree1,
                                        frontier1[itask], frontier1[jtask],
                                        &hist);
                treecorr_process21_auto(&context, tree1,
                                        frontier1[jtask], frontier1[itask],
                                        &hist);
                for (INTEGER ktask = jtask + 1;
                     ktask < frontier_count1; ktask++) {
                    const INTEGER nodes[3] = {
                        frontier1[itask], frontier1[jtask], frontier1[ktask]
                    };
                    treecorr_process111(&context, trees, nodes, 0x3fU, &hist);
                }
            }
        } else {
            const fcfc_balltreeptr trees[3] = {tree1, tree2, tree2};

            for (INTEGER jtask = 0; jtask < frontier_count2; jtask++) {
                treecorr_process12_cross(
                    &context, tree1, frontier1[itask],
                    tree2, frontier2[jtask], &hist);
                for (INTEGER ktask = jtask + 1;
                     ktask < frontier_count2; ktask++) {
                    const INTEGER nodes[3] = {
                        frontier1[itask], frontier2[jtask], frontier2[ktask]
                    };
                    treecorr_process111(&context, trees, nodes, 0x03U, &hist);
                }
            }
        }
        triple_task_body_counts[itask] = hist.body_triples;
        triple_task_cell_counts[itask] = hist.cell_triples;
    }

    if (treecorr_distributed_reduce_reals(
            cmd, triple_task_histograms,
            (size_t)frontier_count1 * triple_values_per_task) == FAILURE)
        triple_reduction_status = FAILURE;
    if (treecorr_distributed_reduce_integers(
            cmd, triple_task_body_counts,
            (size_t)frontier_count1) == FAILURE)
        triple_reduction_status = FAILURE;
    if (treecorr_distributed_reduce_integers(
            cmd, triple_task_cell_counts,
            (size_t)frontier_count1) == FAILURE)
        triple_reduction_status = FAILURE;
    if (treecorr_distributed_consensus(
            cmd, triple_reduction_status,
            "two-ball triple histogram reduction") == FAILURE)
        goto cleanup;

    if (treecorr_distributed_publish()) {
    for (INTEGER itask = 0; itask < frontier_count1; itask++) {
        const real *base = triple_task_histograms
                         + (size_t)itask * triple_values_per_task;
        const size_t plane = stride * stride;

        for (int m = 0; m < triple_orders; m++) {
            const real *zcos = base
                + ((size_t)TREECORR_ZETA_COS * (size_t)triple_orders
                   + (size_t)m)
                * plane;
            const real *zsin = base
                + ((size_t)TREECORR_ZETA_SIN * (size_t)triple_orders
                   + (size_t)m)
                * plane;
            const real *zsincos = base
                + ((size_t)TREECORR_ZETA_SINCOS * (size_t)triple_orders
                   + (size_t)m)
                * plane;
            const real *zcossin = base
                + ((size_t)TREECORR_ZETA_COSSIN * (size_t)triple_orders
                   + (size_t)m)
                * plane;
            for (int n1 = 1; n1 <= cmd->sizeHistN; n1++)
                for (int n2 = 1; n2 <= cmd->sizeHistN; n2++) {
                    const size_t index = (size_t)n1 * stride + (size_t)n2;
                    gd->histZetaMcos[m + 1][n1][n2] += zcos[index];
                    gd->histZetaMsin[m + 1][n1][n2] += zsin[index];
                    gd->histZetaMsincos[m + 1][n1][n2] += zsincos[index];
                    gd->histZetaMcossin[m + 1][n1][n2] += zcossin[index];
                }
        }
        triple_body_total += triple_task_body_counts[itask];
        triple_cell_total += triple_task_cell_counts[itask];
    }

    if (treecorr_normalize_3pcf(cmd)) {
        for (int n1 = 1; n1 <= cmd->sizeHistN; n1++)
            for (int n2 = 1; n2 <= cmd->sizeHistN; n2++) {
            const size_t index = (size_t)n1 * stride + (size_t)n2;
            real denominator = 0.0;
            for (INTEGER itask = 0; itask < frontier_count1; itask++) {
                const real *base = triple_task_histograms
                    + (size_t)itask * triple_values_per_task;
                denominator += base[
                    TREECORR_ZETA_COMPONENTS
                    * (size_t)triple_orders * stride * stride + index];
            }
            for (int m = 1; m <= triple_orders; m++) {
                gd->histZetaMcos[m][n1][n2] = cballs_normalize_or_zero(
                    gd->histZetaMcos[m][n1][n2], denominator);
                gd->histZetaMsin[m][n1][n2] = cballs_normalize_or_zero(
                    gd->histZetaMsin[m][n1][n2], denominator);
                gd->histZetaMsincos[m][n1][n2] = cballs_normalize_or_zero(
                    gd->histZetaMsincos[m][n1][n2], denominator);
                gd->histZetaMcossin[m][n1][n2] = cballs_normalize_or_zero(
                    gd->histZetaMcossin[m][n1][n2], denominator);
            }
        }
    }
    }
    }
#endif /* THREEPCFCONVERGENCE */

    gd->cpusearch = CPUTIME - cpustart;
#ifdef TWOPCF
    if (run_2pcf)
        verb_print(cmd->verbose,
                   "%s: nbbcalc = %" INTEGER_FMT ", nbccalc = %" INTEGER_FMT
                   ", frontier tasks = %" INTEGER_FMT "\n",
                   cmd->searchMethod, gd->nbbcalc, gd->nbccalc,
                   frontier_count1);
#endif
#ifdef THREEPCFCONVERGENCE
    if (run_3pcf)
        verb_print(cmd->verbose,
                   "%s: body-triplet evaluations = %" INTEGER_FMT
                   ", accepted cell orientations = %" INTEGER_FMT "\n",
                   cmd->searchMethod, triple_body_total, triple_cell_total);
#endif
    verb_print(cmd->verbose, "Going out: CPU time = %lf\n", gd->cpusearch);
    status = SUCCESS;

cleanup:
#ifdef THREEPCFCONVERGENCE
    free(triple_task_cell_counts);
    free(triple_task_body_counts);
    free(triple_task_histograms);
#endif
    if (!auto_correlation) free(frontier2);
    free(frontier1);
    if (tree2 != tree1) fcfc_balltree_free(tree2);
    fcfc_balltree_free(tree1);
    return status;
}
#endif /* !TREECORR_TASK_FRONTIER_ENGINE */

#ifdef TREECORR_TASK_FRONTIER_ENGINE
#ifdef TREECORR_LOG_MULTIPOLE_ENGINE
static int treecorr_search_direct_triples(
#else
global int searchcalc_balltree_2balls_omp(
#endif
        struct cmdline_data *cmd, struct global_data *gd,
        bodyptr *btab, INTEGER *nbody, INTEGER ipmin, INTEGER *ipmax,
        int cat1, int cat2)
{
#ifdef THREEPCFCONVERGENCE
    const bool auto_correlation = cat1 == cat2;
    const bool only_2pcf = scanopt(cmd->options, "only-2pcf");
    const bool only_3pcf = scanopt(cmd->options, "only-3pcf");
#ifdef TWOPCF
    const bool run_2pcf = !only_3pcf;
#else
    const bool run_2pcf = FALSE;
#endif
    const bool run_3pcf = !only_2pcf;
    const size_t stride = (size_t)cmd->sizeHistN + 1;
    const int triple_orders = cmd->mChebyshev + 1;
    const INTEGER target_tasks = treecorr_task_target(stride, triple_orders);
    const int leaf_capacity = scanopt(cmd->options, "treecorr-bucket-leaves")
        ? cmd->nsmooth : 1;
    treecorr_search_context context;
    fcfc_balltreeptr tree1 = NULL;
    fcfc_balltreeptr tree2 = NULL;
    treecorr_triple_task *tasks = NULL;
    INTEGER *pair_frontier1 = NULL;
    INTEGER *pair_frontier2 = NULL;
    INTEGER task_count = 0;
    INTEGER pair_frontier_count1 = 0;
    INTEGER pair_frontier_count2 = 0;
    real *task_histograms = NULL;
    INTEGER *task_body_counts = NULL;
    INTEGER *task_cell_counts = NULL;
    size_t values_per_task = 0;
    INTEGER body_total = 0;
    INTEGER cell_total = 0;
    int status = FAILURE;
    const double cpustart = CPUTIME;

    if (only_2pcf && only_3pcf) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: only-2pcf and only-3pcf are mutually exclusive",
                 TREECORR_METHOD_NAME);
        return FAILURE;
    }
    if (only_2pcf && !run_2pcf) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: only-2pcf requires TWOPCFON=1",
                 TREECORR_METHOD_NAME);
        return FAILURE;
    }
    if (ipmin != 1 || ipmax[cat1] != nbody[cat1]) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s requires the complete pivot catalog",
                 TREECORR_METHOD_NAME);
        return FAILURE;
    }
    if (cmd->nsmooth <= 0) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s requires nsmooth > 0", TREECORR_METHOD_NAME);
        return FAILURE;
    }

    verb_print(cmd->verbose, "Search: Running %s", cmd->searchMethod);
#ifdef TWOPCF
    if (run_2pcf) verb_print(cmd->verbose, " with dual-node 2PCF");
#endif
    if (run_3pcf)
        verb_print(cmd->verbose, " with direct triple-node 3PCF");
    verb_print(cmd->verbose, "\n");
    if (cballs_opt_no_two_balls(cmd))
        verb_print(cmd->verbose,
                   "no-two-balls: forcing exact body-triplet accumulation\n");
    else
        verb_print(cmd->verbose,
                   "two-ball radial/angular acceptance enabled; theta=%g\n",
                   cmd->theta);
#ifdef TREECORR_NATIVE_BINARY_VIEW
    verb_print(cmd->verbose,
               "native-octree binary view with singleton body leaves\n");
#else
    verb_print(cmd->verbose,
               "ball-tree leaf capacity = %d%s\n", leaf_capacity,
               leaf_capacity == 1 ? " (TreeCorr singleton policy)"
                                  : " (nsmooth bucket policy)");
#endif
    if (cballs_opt_weights_norm(cmd))
        verb_print(cmd->verbose,
                   "using TreeCorr-style triplet-weight normalization\n");

#ifdef OPENMPCODE
    ThreadCount(cmd, gd, nbody[cat1], cat1);
#endif

    if (search_init_gd_hist_sincos(cmd, gd) == FAILURE) goto cleanup;
    if (fcfc_balltree_build(cmd, gd, btab[cat1], nbody[cat1],
                            leaf_capacity, &tree1) == FAILURE)
        goto cleanup;
    TREECORR_PUBLISH_NODE_COUNT(gd, cat1, tree1->nnode);

    if (auto_correlation) {
        tree2 = tree1;
    } else {
        if (fcfc_balltree_build(cmd, gd, btab[cat2], nbody[cat2],
                                leaf_capacity, &tree2) == FAILURE)
            goto cleanup;
        TREECORR_PUBLISH_NODE_COUNT(gd, cat2, tree2->nnode);
    }

    context.cmd = cmd;
    context.gd = gd;
    context.use_two_balls = !cballs_opt_no_two_balls(cmd);
    context.use_three_cells = !cballs_opt_no_two_balls(cmd);
    context.weighted = cballs_opt_weights_norm(cmd);
    treecorr_initialize_angular_tolerance(&context);

 #ifdef TWOPCF
    if (run_2pcf) {
        if (fcfc_balltree_frontier(
                cmd, tree1, TREECORR_PAIR_FRONTIER_TARGET,
                &pair_frontier1, &pair_frontier_count1) == FAILURE)
            goto cleanup;
        if (auto_correlation) {
            pair_frontier2 = pair_frontier1;
            pair_frontier_count2 = pair_frontier_count1;
        } else if (fcfc_balltree_frontier(
                cmd, tree2, TREECORR_PAIR_FRONTIER_TARGET,
                &pair_frontier2, &pair_frontier_count2) == FAILURE) {
            goto cleanup;
        }
        if (treecorr_run_pair_tasks(
                &context, tree1, tree2, auto_correlation,
                pair_frontier1, pair_frontier_count1,
                pair_frontier2, pair_frontier_count2) == FAILURE)
            goto cleanup;
    }
#endif

    if (run_3pcf) {
    if (treecorr_build_task_frontier(
            cmd, tree1, tree2, auto_correlation, target_tasks,
            &tasks, &task_count) == FAILURE)
        goto cleanup;
    if (treecorr_allocate_triple_histograms(
            cmd, task_count, stride, triple_orders,
            &task_histograms, &values_per_task,
            &task_body_counts, &task_cell_counts) == FAILURE)
        goto cleanup;

#pragma omp parallel for schedule(dynamic,1)
    for (INTEGER itask = 0; itask < task_count; itask++) {
        real *base = task_histograms
                   + (size_t)itask * values_per_task;
        treecorr_triple_histogram hist;

        hist.components = base;
        hist.stride = stride;
        hist.plane = stride * stride;
        hist.order_plane = (size_t)triple_orders * hist.plane;
        hist.orders = triple_orders;
        hist.normalization = base
            + TREECORR_ZETA_COMPONENTS * hist.order_plane;
        hist.body_triples = 0;
        hist.cell_triples = 0;
        treecorr_run_task(&context, &tasks[itask], &hist);
        task_body_counts[itask] = hist.body_triples;
        task_cell_counts[itask] = hist.cell_triples;
    }

    for (INTEGER itask = 0; itask < task_count; itask++) {
        const real *base = task_histograms
                         + (size_t)itask * values_per_task;
        const size_t plane = stride * stride;

        for (int m = 0; m < triple_orders; m++) {
            const real *zcos = base
                + ((size_t)TREECORR_ZETA_COS * (size_t)triple_orders
                   + (size_t)m)
                * plane;
            const real *zsin = base
                + ((size_t)TREECORR_ZETA_SIN * (size_t)triple_orders
                   + (size_t)m)
                * plane;
            const real *zsincos = base
                + ((size_t)TREECORR_ZETA_SINCOS * (size_t)triple_orders
                   + (size_t)m)
                * plane;
            const real *zcossin = base
                + ((size_t)TREECORR_ZETA_COSSIN * (size_t)triple_orders
                   + (size_t)m)
                * plane;

            for (int n1 = 1; n1 <= cmd->sizeHistN; n1++)
                for (int n2 = 1; n2 <= cmd->sizeHistN; n2++) {
                    const size_t index = (size_t)n1 * stride + (size_t)n2;
                    gd->histZetaMcos[m + 1][n1][n2] += zcos[index];
                    gd->histZetaMsin[m + 1][n1][n2] += zsin[index];
                    gd->histZetaMsincos[m + 1][n1][n2] += zsincos[index];
                    gd->histZetaMcossin[m + 1][n1][n2] += zcossin[index];
                }
        }
        body_total += task_body_counts[itask];
        cell_total += task_cell_counts[itask];
    }

    if (treecorr_normalize_3pcf(cmd)) {
        for (int n1 = 1; n1 <= cmd->sizeHistN; n1++)
            for (int n2 = 1; n2 <= cmd->sizeHistN; n2++) {
            const size_t index = (size_t)n1 * stride + (size_t)n2;
            real denominator = 0.0;

            for (INTEGER itask = 0; itask < task_count; itask++) {
                const real *base = task_histograms
                    + (size_t)itask * values_per_task;
                denominator += base[
                    TREECORR_ZETA_COMPONENTS
                    * (size_t)triple_orders * stride * stride + index];
            }
            for (int m = 1; m <= triple_orders; m++) {
                gd->histZetaMcos[m][n1][n2] = cballs_normalize_or_zero(
                    gd->histZetaMcos[m][n1][n2], denominator);
                gd->histZetaMsin[m][n1][n2] = cballs_normalize_or_zero(
                    gd->histZetaMsin[m][n1][n2], denominator);
                gd->histZetaMsincos[m][n1][n2] = cballs_normalize_or_zero(
                    gd->histZetaMsincos[m][n1][n2], denominator);
                gd->histZetaMcossin[m][n1][n2] = cballs_normalize_or_zero(
                    gd->histZetaMcossin[m][n1][n2], denominator);
            }
        }
    }
    }

    gd->cpusearch = CPUTIME - cpustart;
#ifdef TWOPCF
    if (run_2pcf)
        verb_print(cmd->verbose,
                   "%s: nbbcalc = %" INTEGER_FMT
                   ", nbccalc = %" INTEGER_FMT
                   ", frontier tasks = %" INTEGER_FMT "\n",
                   cmd->searchMethod, gd->nbbcalc, gd->nbccalc,
                   pair_frontier_count1);
#endif
    if (run_3pcf)
        verb_print(cmd->verbose,
                   "%s: body-triplet evaluations = %" INTEGER_FMT
                   ", accepted cell orientations = %" INTEGER_FMT
                   ", task frontier = %" INTEGER_FMT "\n",
                   cmd->searchMethod, body_total, cell_total, task_count);
    verb_print(cmd->verbose, "Going out: CPU time = %lf\n", gd->cpusearch);
    status = SUCCESS;

cleanup:
    free(task_cell_counts);
    free(task_body_counts);
    free(task_histograms);
    free(tasks);
    if (!auto_correlation) free(pair_frontier2);
    free(pair_frontier1);
    if (tree2 != tree1) fcfc_balltree_free(tree2);
    fcfc_balltree_free(tree1);
    return status;
#else
    (void)gd;
    (void)btab;
    (void)nbody;
    (void)ipmin;
    (void)ipmax;
    (void)cat1;
    (void)cat2;
    snprintf(cmd->error_message, _ERRORMSGSIZE_,
             "%s requires TPCFON=1", TREECORR_METHOD_NAME);
    return FAILURE;
#endif /* THREEPCFCONVERGENCE */
}
#endif /* TREECORR_TASK_FRONTIER_ENGINE */
