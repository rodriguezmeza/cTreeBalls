/*
 * dual-node-style dual-ball 2PCF and body-pivot LogMultipole scalar 3PCF.
 * Direct triples remain available only as an explicit validation option.
 *
 * The node-pair recursion and split heuristic are adapted from dual-node:
 * Copyright (c) 2003-2024 Mike Jarvis, used under its BSD-style license.
 * The ball-tree builder is the FCFC-derived implementation shared with
 * balltree-omp; its source file carries the full FCFC MIT notice.
 */

#include <float.h>
#include <limits.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#ifdef OPENMPCODE
#include <omp.h>
#endif

#include "globaldefs.h"
#include "fcfc_balltree.h"

#if defined(__APPLE__) && !defined(DUAL_NODE_METHOD_NAME)
#define DUAL_NODE_USE_ACCELERATE_VFORCE 1
#endif
#ifdef DUAL_NODE_USE_ACCELERATE_VFORCE
extern void vvlog(double *, const double *, const int *);
extern void vvlogf(float *, const float *, const int *);
#endif
#ifdef DUAL_NODE_USE_SLEEF
#include "dual_node_sleef_log.h"
#endif

#ifndef DUAL_NODE_METHOD_NAME
#define DUAL_NODE_METHOD_NAME "balltree-2balls-omp"
#define BALLTREE_2BALLS_PRIMARY_FEATURES 1
#define BALLTREE_2BALLS_FULL_FUNCTION searchcalc_balltree_2balls_full_omp
#define BALLTREE_2BALLS_SEARCH_FUNCTION searchcalc_balltree_2balls_omp
#define BALLTREE_2BALLS_LEGACY_FUNCTION searchcalc_balltree_omp
#endif

#ifdef BALLTREE_2BALLS_LEGACY_FUNCTION
#include "protodefs_balltree_omp.h"
#endif

#ifdef BALLTREE_2BALLS_SEARCH_FUNCTION
#define searchcalc_balltree_2balls_omp BALLTREE_2BALLS_FULL_FUNCTION
#endif

#ifdef BALLTREE_2BALLS_PRIMARY_FEATURES
static inline int balltree_2balls_prepare_pivots(
        struct cmdline_data *cmd, struct global_data *gd,
        bodyptr *body_table, INTEGER *body_count,
        INTEGER pivot_minimum, INTEGER *pivot_maximum,
        int pivot_catalog, int neighbor_catalog)
{
#ifdef SMOOTHPIVOT
    return prepare_smooth_pivots(
        cmd, gd, body_table, body_count, pivot_minimum, pivot_maximum,
        pivot_catalog, neighbor_catalog);
#else
    (void)cmd;
    (void)gd;
    (void)body_table;
    (void)body_count;
    (void)pivot_minimum;
    (void)pivot_maximum;
    (void)pivot_catalog;
    (void)neighbor_catalog;
    return SUCCESS;
#endif
}

#define DUAL_NODE_PREPARE_PIVOTS balltree_2balls_prepare_pivots
#define DUAL_NODE_BUILD_PIVOT_TREE(cmd, gd, btab, nbody, leaf, result) \
    fcfc_balltree_build_scalar_role(                                  \
        (cmd), (gd), (btab), (nbody), (leaf), TRUE, (result))
#define DUAL_NODE_BUILD_NEIGHBOR_TREE(cmd, gd, btab, nbody, leaf, result) \
    fcfc_balltree_build_scalar_role(                                     \
        (cmd), (gd), (btab), (nbody), (leaf), FALSE, (result))
#define DUAL_NODE_SEPARATE_PIVOT_TREE(cmd) cballs_opt_smooth_pivot(cmd)
#ifdef SMOOTHPIVOT
#define DUAL_NODE_PIVOT_IS_ACTIVE(context, pivot) \
    ((!cballs_opt_read_mask((context)->cmd) \
      || Mask(pivot) != MASK_NODE_MASKED) \
     && (!cballs_opt_smooth_pivot((context)->cmd) || Update(pivot) != FALSE))
#define DUAL_NODE_PIVOT_FIELD(context, pivot) \
    (cballs_opt_smooth_pivot((context)->cmd) \
     ? KappaRmin(pivot) : dual_node_body_weighted_field((context), (pivot)))
#define DUAL_NODE_PIVOT_NORMALIZATION(context, pivot) \
    (cballs_opt_smooth_pivot((context)->cmd) \
     ? (cballs_opt_weights_norm((context)->cmd) \
        ? WeightRmin(pivot) : (real)MAX(NbRmin(pivot), 1)) \
     : dual_node_body_normalization((context), (pivot)))
#else
#define DUAL_NODE_PIVOT_IS_ACTIVE(context, pivot) \
    (!cballs_opt_read_mask((context)->cmd) || Mask(pivot) != MASK_NODE_MASKED)
#endif
#define DUAL_NODE_CONTEXT_WEIGHTED(cmd) \
    (cballs_opt_weights_norm(cmd) || cballs_opt_smooth_pivot(cmd))
#define DUAL_NODE_ADAPTIVE_PAIR_LEAVES 1
#define DUAL_NODE_PAIR_BATCH_SIZE 256
#define DUAL_NODE_USE_NATURAL_LOG_BINS 1
#if !defined(DUAL_NODE_DISTRIBUTED_ENGINE)
#define DUAL_NODE_SCALAR_TREE_CACHE 1
#define DUAL_NODE_RELEASE_TREE(tree) fcfc_balltree_release(tree)
#endif
#ifdef THREEPCFCONVERGENCE
#define DUAL_NODE_TASK_FRONTIER_ENGINE 1
#define DUAL_NODE_LOG_MULTIPOLE_ENGINE 1
#define DUAL_NODE_BODY_PIVOT_LOG_MULTIPOLE 1
#define DUAL_NODE_PERSISTENT_NEIGHBOR_FRONTIER 1
#ifndef DUAL_NODE_DISTRIBUTED_ENGINE
#define DUAL_NODE_PIVOT_PROGRESS 1
#endif
#endif
#endif

#ifdef BALLS4SCANLEV
#define DUAL_NODE_SCAN_LEVEL_FRONTIER 1
#endif
#ifndef DUAL_NODE_PUBLISH_NODE_COUNT
#define DUAL_NODE_PUBLISH_NODE_COUNT(gd, catalog, count) \
    ((gd)->ncellTable[(catalog)] = (count))
#endif
#ifndef DUAL_NODE_PREPARE_PIVOTS
#define DUAL_NODE_PREPARE_PIVOTS(cmd, gd, btab, nbody, ipmin, ipmax, cat1, cat2) \
    SUCCESS
#endif
#ifndef DUAL_NODE_BUILD_PIVOT_TREE
#define DUAL_NODE_BUILD_PIVOT_TREE(cmd, gd, btab, nbody, leaf, result) \
    fcfc_balltree_build((cmd), (gd), (btab), (nbody), (leaf), (result))
#endif
#ifndef DUAL_NODE_BUILD_NEIGHBOR_TREE
#define DUAL_NODE_BUILD_NEIGHBOR_TREE(cmd, gd, btab, nbody, leaf, result) \
    fcfc_balltree_build((cmd), (gd), (btab), (nbody), (leaf), (result))
#endif
#ifndef DUAL_NODE_SEPARATE_PIVOT_TREE
#define DUAL_NODE_SEPARATE_PIVOT_TREE(cmd) FALSE
#endif
#ifndef DUAL_NODE_PIVOT_IS_ACTIVE
#define DUAL_NODE_PIVOT_IS_ACTIVE(context, pivot) TRUE
#endif
#ifndef DUAL_NODE_PIVOT_FIELD
#define DUAL_NODE_PIVOT_FIELD(context, pivot) \
    dual_node_body_weighted_field((context), (pivot))
#endif
#ifndef DUAL_NODE_PIVOT_NORMALIZATION
#define DUAL_NODE_PIVOT_NORMALIZATION(context, pivot) \
    dual_node_body_normalization((context), (pivot))
#endif
#ifndef DUAL_NODE_CONTEXT_WEIGHTED
#define DUAL_NODE_CONTEXT_WEIGHTED(cmd) cballs_opt_weights_norm(cmd)
#endif
#ifndef DUAL_NODE_RELEASE_TREE
#define DUAL_NODE_RELEASE_TREE(tree) fcfc_balltree_free(tree)
#endif

#define DUAL_NODE_SPLIT_FACTOR ((real)0.585)
#define DUAL_NODE_PAIR_FRONTIER_TARGET ((INTEGER)256)
#define DUAL_NODE_TRIPLE_FRONTIER_TARGET ((INTEGER)32)
#define DUAL_NODE_LARGE_DIRECT_3PCF ((INTEGER)1000000)
#define DUAL_NODE_TRIPLE_TASK_TARGET ((INTEGER)256)
#define DUAL_NODE_TRIPLE_TASK_MEMORY ((size_t)128 * 1024 * 1024)

typedef struct {
    double build;
    double frontier;
    double pair_traversal;
    double pivot_transport_thread;
    double scratch_clear_thread;
    double multipole_products_thread;
    double reduction;
} dual_node_phase_timers;

#ifdef DUAL_NODE_PIVOT_PROGRESS
typedef struct {
    bool enabled;
    INTEGER total;
    INTEGER completed;
    INTEGER reported;
    INTEGER interval;
    INTEGER batch;
    double started;
} dual_node_pivot_progress;
#endif

static inline double dual_node_timer_now(void)
{
#ifdef OPENMPCODE
    return omp_get_wtime();
#else
    return CPUTIME;
#endif
}

#ifdef TWOPCF
typedef struct {
    real *pair_count;
    real *weight_product;
    real *field_product;
    INTEGER body_pairs;
    INTEGER cell_pairs;
} dual_node_histogram;
#endif

typedef struct {
    struct cmdline_data *cmd;
    struct global_data *gd;
    bool reuse_enabled;
    real reuse_bin_theta, reuse_phase_ratio, reuse_max_relative_width;
    bool use_two_balls;
    bool use_bin_theta;
    bool use_three_cells;
    bool weighted;
    real angular_tolerance;
    real max_angular_ratio;
    real minimum2;
    real maximum2;
    real natural_log_scale;
    real logarithmic_bin_size;
    real logarithmic_minimum2;
    real logarithmic_maximum;
    real bin_size;
    real bin_theta;
    real theta2;
    real bin_theta2;
    real half_bin_plus_slop2;
    bool profile;
    dual_node_phase_timers *timers;
#ifdef DUAL_NODE_PIVOT_PROGRESS
    dual_node_pivot_progress *progress;
#endif
} dual_node_search_context;

#ifdef DUAL_NODE_PIVOT_PROGRESS
static void dual_node_print_pivot_progress(
        const dual_node_search_context *context)
{
    const dual_node_pivot_progress *progress = context->progress;

    verb_print_min_info(
        context->cmd->verbose, context->cmd->verbose_log, context->gd->outlog,
        "%s: 3PCF progress: completed pivots %" INTEGER_FMT
        " / %" INTEGER_FMT " (%.1f%%); elapsed %.2f s\n",
        DUAL_NODE_METHOD_NAME, progress->completed, progress->total,
        100.0 * (double)progress->completed / (double)progress->total,
        dual_node_timer_now() - progress->started);
}

static void dual_node_publish_pivot_progress(
        const dual_node_search_context *context, INTEGER count)
{
    dual_node_pivot_progress *progress = context->progress;

    if (progress == NULL || !progress->enabled || count == 0) return;
    /* Publish small local batches, serializing counters and flushed output
     * together so dynamic OpenMP scheduling cannot reorder progress lines. */
#pragma omp critical(cballs_two_balls_progress)
    {
        progress->completed += count;
        if (progress->completed - progress->reported >= progress->interval
            || progress->completed == progress->total) {
            dual_node_print_pivot_progress(context);
            progress->reported = progress->completed;
        }
    }
}
#endif

static void dual_node_initialize_radial_context(
        dual_node_search_context *context, struct cmdline_data *cmd,
        struct global_data *gd)
{
    const real log10_natural = rlog(10.0);

    context->cmd = cmd;
    context->gd = gd;
    context->minimum2 = rsqr(cmd->rminHist);
    context->maximum2 = rsqr(cmd->rangeN);
    context->natural_log_scale = gd->i_deltaR / log10_natural;
    context->logarithmic_bin_size = log10_natural * gd->deltaR;
    context->logarithmic_minimum2 = cmd->rminHist > 0.0
        ? rlog(context->minimum2) : 0.0;
    context->logarithmic_maximum = rlog(cmd->rangeN);
    context->bin_size = cmd->useLogHist
        ? context->logarithmic_bin_size : gd->deltaR;
    context->bin_theta = cmd->theta * context->bin_size;
    context->theta2 = rsqr(cmd->theta);
    context->bin_theta2 = rsqr(context->bin_theta);
    context->half_bin_plus_slop2 =
        0.25 * rsqr(context->bin_size + context->bin_theta);
}

#ifdef DUAL_NODE_DISTRIBUTED_ENGINE
static inline bool dual_node_distributed_task_owned(INTEGER task)
{
    return DUAL_NODE_DISTRIBUTED_TASK_OWNED(task);
}

static inline bool dual_node_distributed_publish(void)
{
    return DUAL_NODE_DISTRIBUTED_IS_ROOT();
}

static inline int dual_node_distributed_consensus(
        struct cmdline_data *cmd, int status, const char *operation)
{
    return DUAL_NODE_DISTRIBUTED_CONSENSUS(cmd, status, operation);
}

static inline int dual_node_distributed_reduce_reals(
        struct cmdline_data *cmd, real *values, size_t count)
{
    return DUAL_NODE_DISTRIBUTED_REDUCE_REALS(cmd, values, count);
}

static inline int dual_node_distributed_reduce_integers(
        struct cmdline_data *cmd, INTEGER *values, size_t count)
{
    return DUAL_NODE_DISTRIBUTED_REDUCE_INTEGERS(cmd, values, count);
}
#else
static inline bool dual_node_distributed_task_owned(INTEGER task)
{
    (void)task;
    return TRUE;
}

static inline bool dual_node_distributed_publish(void)
{
    return TRUE;
}

static inline int dual_node_distributed_consensus(
        struct cmdline_data *cmd, int status, const char *operation)
{
    (void)cmd;
    (void)operation;
    return status;
}

static inline int dual_node_distributed_reduce_reals(
        struct cmdline_data *cmd, real *values, size_t count)
{
    (void)cmd;
    (void)values;
    (void)count;
    return SUCCESS;
}

static inline int dual_node_distributed_reduce_integers(
        struct cmdline_data *cmd, INTEGER *values, size_t count)
{
    (void)cmd;
    (void)values;
    (void)count;
    return SUCCESS;
}
#endif

static INTEGER dual_node_pair_frontier_target(
        const struct cmdline_data *cmd)
{
    INTEGER workers = 1;
    INTEGER target;

#ifdef OPENMPCODE
    if (cmd->numthreads > 0) workers = (INTEGER)cmd->numthreads;
#else
    (void)cmd;
#endif
#ifdef DUAL_NODE_DISTRIBUTED_ENGINE
    if (DUAL_NODE_DISTRIBUTED_SIZE() > 1) {
        const INTEGER ranks = (INTEGER)DUAL_NODE_DISTRIBUTED_SIZE();
        if (workers <= DUAL_NODE_PAIR_FRONTIER_TARGET / ranks)
            workers *= ranks;
        else
            workers = DUAL_NODE_PAIR_FRONTIER_TARGET;
    }
#endif
#ifdef DUAL_NODE_SCAN_LEVEL_FRONTIER
    target = workers <= DUAL_NODE_PAIR_FRONTIER_TARGET / 4
        ? 4 * workers : DUAL_NODE_PAIR_FRONTIER_TARGET;
    target = MAX((INTEGER)64, target);
#else
    target = workers <= DUAL_NODE_PAIR_FRONTIER_TARGET / 4
        ? 4 * workers : DUAL_NODE_PAIR_FRONTIER_TARGET;
    target = MAX((INTEGER)16, target);
#endif
    return MIN(DUAL_NODE_PAIR_FRONTIER_TARGET, target);
}

static INTEGER dual_node_triple_frontier_target(
        const struct cmdline_data *cmd)
{
#ifdef DUAL_NODE_SCAN_LEVEL_FRONTIER
    return dual_node_pair_frontier_target(cmd);
#else
    (void)cmd;
    return DUAL_NODE_TRIPLE_FRONTIER_TARGET;
#endif
}

#ifdef TWOPCF
static int dual_node_run_pair_tasks(
        const dual_node_search_context *, fcfc_balltreeptr,
        fcfc_balltreeptr, bool, const INTEGER *, INTEGER,
        const INTEGER *, INTEGER);
#endif

#ifdef THREEPCFCONVERGENCE
static inline int dual_node_window_orders(const struct cmdline_data *cmd)
{
    if (!cballs_opt_edge_corrections(cmd)) return 0;
    if (cmd->mChebyshev < 0 || cmd->mChebyshev > (INT_MAX - 1) / 2)
        return -1;
    return 2 * cmd->mChebyshev + 1;
}

static inline bool dual_node_normalize_3pcf(const struct cmdline_data *cmd)
{
    return !scanopt(cmd->options, "no-normalize-HistZeta");
}
#endif

static inline real dual_node_body_field(bodyptr p)
{
#ifdef KappaAvgON
    return KappaAvg(p);
#else
    return Kappa(p);
#endif
}

static inline INTEGER dual_node_node_count(const fcfc_ballnode *ball_node)
{
    return ball_node->last - ball_node->first + 1;
}

static inline bool dual_node_node_is_leaf(const fcfc_ballnode *ball_node)
{
    return ball_node->left < 0;
}

static inline real dual_node_center_distance_squared(
        struct cmdline_data *cmd, struct global_data *gd,
        const fcfc_ballnode *node1, const fcfc_ballnode *node2)
{
    compute_vector dr;
    real distance2;

    DOTPSUBV(distance2, dr, node1->center, node2->center);
    if (cmd->usePeriodic) {
        VWrapAll(dr);
        DOTVP(distance2, dr, dr);
    }
    return distance2;
}

static inline real dual_node_center_distance(struct cmdline_data *cmd,
                                            struct global_data *gd,
                                            const fcfc_ballnode *node1,
                                            const fcfc_ballnode *node2)
{
    return rsqrt(dual_node_center_distance_squared(cmd, gd, node1, node2));
}

static inline real dual_node_point_distance_squared(struct cmdline_data *cmd,
                                                   struct global_data *gd,
                                                   const fcfc_ballpoint *p,
                                                   const fcfc_ballpoint *q)
{
    compute_vector dr;
    real distance2;

    DOTPSUBV(distance2, dr, p->pos, q->pos);
    if (cmd->usePeriodic) {
        VWrapAll(dr);
        DOTVP(distance2, dr, dr);
    }
    return distance2;
}

#include "dual_node_radial_bins.h"

#ifdef TWOPCF
static inline void dual_node_accumulate_body_pair(
        const dual_node_search_context *context, dual_node_histogram *hist,
        const fcfc_ballpoint *p, const fcfc_ballpoint *q)
{
    if (p->source != NULL && p->source == q->source) return;
    const int n = dual_node_bin_index_squared(
        context, dual_node_point_distance_squared(
            context->cmd, context->gd, p, q));
    real denominator;
    real numerator;

    if (n < 0) return;
    if (context->weighted) {
        denominator = p->weight * q->weight;
        numerator = p->weighted_kappa * q->weighted_kappa;
    } else {
        denominator = 1.0;
        numerator = p->kappa * q->kappa;
    }
    hist->pair_count[n] += 1.0;
    hist->weight_product[n] += denominator;
    hist->field_product[n] += numerator;
    hist->body_pairs++;
}

#ifdef DUAL_NODE_PAIR_BATCH_SIZE
static void dual_node_vector_log(real *output, const real *input, int count)
{
#ifdef DUAL_NODE_USE_ACCELERATE_VFORCE
#if defined(DOUBLEPREC)
    vvlog(output, input, &count);
#else
    vvlogf(output, input, &count);
#endif
#elif defined(DUAL_NODE_USE_SLEEF)
#if defined(DOUBLEPREC)
    dual_node_sleef_log_double(output, input, (size_t)count);
#else
    dual_node_sleef_log_float(output, input, (size_t)count);
#endif
#else
    for (int i = 0; i < count; i++)
        output[i] = rlog(input[i]);
#endif
}

static void dual_node_accumulate_body_pair_batch_flush(
        const dual_node_search_context *context, dual_node_histogram *hist,
        real *distance2, real *denominator, real *numerator, int count)
{
    int bins[DUAL_NODE_PAIR_BATCH_SIZE];
    const struct cmdline_data *cmd = context->cmd;

    if (count <= 0) return;
    if (cmd->useLogHist && cmd->rminHist > 0.0) {
        real logarithms[DUAL_NODE_PAIR_BATCH_SIZE];
        const real scale = 0.5 * context->natural_log_scale;
#if defined(DUAL_NODE_USE_ACCELERATE_VFORCE) || defined(DUAL_NODE_USE_SLEEF)
        if (count >= 8) {
            dual_node_vector_log(logarithms, distance2, count);
        } else
#endif
        {
            for (int i = 0; i < count; i++)
                logarithms[i] = rlog(distance2[i]);
        }
        for (int i = 0; i < count; i++) {
            const int n = (int)((logarithms[i]
                - context->logarithmic_minimum2) * scale) + 1;
            bins[i] = n >= 1 && n <= cmd->sizeHistN ? n : -1;
        }
    } else {
        for (int i = 0; i < count; i++)
            bins[i] = dual_node_bin_index_squared(context, distance2[i]);
    }

    for (int i = 0; i < count; i++) {
        const int n = bins[i];

        if (n < 0) continue;
        hist->pair_count[n] += 1.0;
        hist->weight_product[n] += denominator[i];
        hist->field_product[n] += numerator[i];
        hist->body_pairs++;
    }
}

static void dual_node_accumulate_body_pair_batch(
        const dual_node_search_context *context, dual_node_histogram *hist,
        const fcfc_balltreeptr tree1, INTEGER first1, INTEGER last1,
        const fcfc_balltreeptr tree2, INTEGER first2, INTEGER last2,
        bool triangular)
{
    real distance2[DUAL_NODE_PAIR_BATCH_SIZE];
    real denominator[DUAL_NODE_PAIR_BATCH_SIZE];
    real numerator[DUAL_NODE_PAIR_BATCH_SIZE];
    int count = 0;

    for (INTEGER i = first1; i <= last1; i++) {
        const INTEGER begin = triangular ? i + 1 : first2;
        const fcfc_ballpoint *p = &tree1->packed_points[i];

        for (INTEGER j = begin; j <= last2; j++) {
            const fcfc_ballpoint *q = &tree2->packed_points[j];
            real d2;

            if (p->source != NULL && p->source == q->source) continue;
            d2 = dual_node_point_distance_squared(
                context->cmd, context->gd, p, q);
            if (!(d2 > context->minimum2 && d2 < context->maximum2))
                continue;
            distance2[count] = d2;
            if (context->weighted) {
                denominator[count] = p->weight * q->weight;
                numerator[count] = p->weighted_kappa * q->weighted_kappa;
            } else {
                denominator[count] = 1.0;
                numerator[count] = p->kappa * q->kappa;
            }
            count++;
            if (count == DUAL_NODE_PAIR_BATCH_SIZE) {
                dual_node_accumulate_body_pair_batch_flush(
                    context, hist, distance2, denominator, numerator, count);
                count = 0;
            }
        }
    }
    dual_node_accumulate_body_pair_batch_flush(
        context, hist, distance2, denominator, numerator, count);
}
#endif

#include "dual_node_pair_acceptance.h"

static inline void dual_node_accumulate_node_pair(
        const dual_node_search_context *context, dual_node_histogram *hist,
        const fcfc_ballnode *node1, const fcfc_ballnode *node2, int n)
{
    const real count = (real)dual_node_node_count(node1)
                     * (real)dual_node_node_count(node2);
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

static void dual_node_process_pair(const dual_node_search_context *context,
                                  const fcfc_balltreeptr tree1, INTEGER i1,
                                  const fcfc_balltreeptr tree2, INTEGER i2,
                                  dual_node_histogram *hist)
{
    const fcfc_ballnode *node1 = &tree1->nodes[i1];
    const fcfc_ballnode *node2 = &tree2->nodes[i2];
    const bool leaf1 = dual_node_node_is_leaf(node1);
    const bool leaf2 = dual_node_node_is_leaf(node2);
    const real distance2 = dual_node_center_distance_squared(
        context->cmd, context->gd, node1, node2);
    int n;

    if (dual_node_pair_outside_range(context, node1, node2, distance2)) return;

    n = dual_node_two_ball_bin(context, node1, node2, distance2);
    if (n == -2) return;
    if (n >= 0) {
        dual_node_accumulate_node_pair(context, hist, node1, node2, n);
        return;
    }

    if (leaf1 && leaf2) {
#ifdef DUAL_NODE_PAIR_BATCH_SIZE
        if (dual_node_node_count(node1) * dual_node_node_count(node2) > 64) {
            dual_node_accumulate_body_pair_batch(
                context, hist, tree1, node1->first, node1->last,
                tree2, node2->first, node2->last, FALSE);
            return;
        }
#endif
        INTEGER p1;
        INTEGER p2;
        for (p1 = node1->first; p1 <= node1->last; p1++)
            for (p2 = node2->first; p2 <= node2->last; p2++)
                dual_node_accumulate_body_pair(
                    context, hist, &tree1->packed_points[p1],
                    &tree2->packed_points[p2]);
        return;
    }

    if (leaf1) {
        dual_node_process_pair(context, tree1, i1, tree2, node2->left, hist);
        dual_node_process_pair(context, tree1, i1, tree2, node2->right, hist);
        return;
    }
    if (leaf2) {
        dual_node_process_pair(context, tree1, node1->left, tree2, i2, hist);
        dual_node_process_pair(context, tree1, node1->right, tree2, i2, hist);
        return;
    }

    {
        bool split1 = FALSE;
        bool split2 = FALSE;
        const real size1 = (real)node1->radius;
        const real size2 = (real)node2->radius;
        real effective_width2;

        if (context->cmd->useLogHist)
            effective_width2 = rsqr(context->cmd->theta * rlog(10.0)
                * context->gd->deltaR) * distance2;
        else
            effective_width2 = rsqr(
                context->cmd->theta * context->gd->deltaR);

        if (size2 > size1) {
            split2 = TRUE;
            if (!(size2 > 2.0 * size1))
                split1 = size1 * size1
                    > rsqr(DUAL_NODE_SPLIT_FACTOR) * effective_width2;
        } else {
            split1 = TRUE;
            if (!(size1 > 2.0 * size2))
                split2 = size2 * size2
                    > rsqr(DUAL_NODE_SPLIT_FACTOR) * effective_width2;
        }

        if (split1 && split2) {
            dual_node_process_pair(context, tree1, node1->left,
                                  tree2, node2->left, hist);
            dual_node_process_pair(context, tree1, node1->left,
                                  tree2, node2->right, hist);
            dual_node_process_pair(context, tree1, node1->right,
                                  tree2, node2->left, hist);
            dual_node_process_pair(context, tree1, node1->right,
                                  tree2, node2->right, hist);
        } else if (split1) {
            dual_node_process_pair(context, tree1, node1->left,
                                  tree2, i2, hist);
            dual_node_process_pair(context, tree1, node1->right,
                                  tree2, i2, hist);
        } else {
            dual_node_process_pair(context, tree1, i1,
                                  tree2, node2->left, hist);
            dual_node_process_pair(context, tree1, i1,
                                  tree2, node2->right, hist);
        }
    }
}

static void dual_node_process_auto(const dual_node_search_context *context,
                                  const fcfc_balltreeptr tree, INTEGER inode,
                                  dual_node_histogram *hist)
{
    const fcfc_ballnode *node = &tree->nodes[inode];

    if (2.0 * (real)node->radius <= context->cmd->rminHist) return;
    if (dual_node_node_is_leaf(node)) {
#ifdef DUAL_NODE_PAIR_BATCH_SIZE
        if (dual_node_node_count(node) > 12) {
            dual_node_accumulate_body_pair_batch(
                context, hist, tree, node->first, node->last,
                tree, node->first, node->last, TRUE);
            return;
        }
#endif
        INTEGER i;
        INTEGER j;
        for (i = node->first; i <= node->last; i++)
            for (j = i + 1; j <= node->last; j++)
                dual_node_accumulate_body_pair(
                    context, hist, &tree->packed_points[i],
                    &tree->packed_points[j]);
        return;
    }

    dual_node_process_auto(context, tree, node->left, hist);
    dual_node_process_auto(context, tree, node->right, hist);
    dual_node_process_pair(context, tree, node->left,
                          tree, node->right, hist);
}
#endif /* TWOPCF */

#ifdef THREEPCFCONVERGENCE

enum {
    DUAL_NODE_ZETA_COS = 0,
    DUAL_NODE_ZETA_SIN = 1,
    DUAL_NODE_ZETA_SINCOS = 2,
    DUAL_NODE_ZETA_COSSIN = 3,
    DUAL_NODE_ZETA_COMPONENTS = 4
};

enum {
    DUAL_NODE_TRIPLE_OUTSIDE = 0,
    DUAL_NODE_TRIPLE_SPLIT = 1,
    DUAL_NODE_TRIPLE_ACCEPT = 2
};

typedef struct {
    real *components;
    real *normalization;
    real *window_re;
    real *window_im;
    size_t stride;
    size_t plane;
    size_t order_plane;
    int orders;
    int window_orders;
    INTEGER body_triples;
    INTEGER cell_triples;
} dual_node_triple_histogram;

#include "dual_node_edge_correction.h"

static void dual_node_initialize_triple_histogram(
        dual_node_triple_histogram *hist, real *base, size_t stride,
        int orders, int window_orders)
{
    hist->components = base;
    hist->stride = stride;
    hist->plane = stride * stride;
    hist->order_plane = (size_t)orders * hist->plane;
    hist->orders = orders;
    hist->normalization = base + DUAL_NODE_ZETA_COMPONENTS * hist->order_plane;
    hist->window_orders = window_orders;
    hist->window_re = window_orders ? hist->normalization + hist->plane : NULL;
    hist->window_im = window_orders
        ? hist->window_re + (size_t)window_orders * hist->plane : NULL;
    hist->body_triples = 0;
    hist->cell_triples = 0;
}

static const unsigned char dual_node_permutations[6][3] = {
    {0, 1, 2}, {0, 2, 1},
    {1, 0, 2}, {1, 2, 0},
    {2, 0, 1}, {2, 1, 0}
};

static inline real *dual_node_zeta_component(
        dual_node_triple_histogram *hist, int component, int order)
{
    return hist->components
        + ((size_t)component * (size_t)hist->orders + (size_t)order)
        * hist->plane;
}

static inline real dual_node_position_distance(
        const dual_node_search_context *context,
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

static inline bool dual_node_leg_outside(
        const dual_node_search_context *context, real distance, real size)
{
    return distance + size <= context->cmd->rminHist
        || distance - size >= context->cmd->rangeN;
}

static int dual_node_sloppy_bin(const dual_node_search_context *context,
                               real distance, real size)
{
    if (size > dual_node_bin_theta_width(context, distance)) return -1;
    return dual_node_bin_index(context, distance);
}

static bool dual_node_polar_coordinates_from_displacement(
        const dual_node_search_context *context,
        const cballs_storage_real *pivot,
        const cballs_storage_real *neighbor,
        const compute_vector dr, real distance,
        real *cosphi, real *sinphi)
{
    (void)context;
    (void)neighbor;
    (void)distance;
    return cballs_angular_phase(pivot, dr, cosphi, sinphi);
}

static bool dual_node_polar_coordinates(
        const dual_node_search_context *context,
        const cballs_storage_real *pivot,
        const cballs_storage_real *neighbor,
        real *distance, real *cosphi, real *sinphi)
{
    compute_vector dr;

    *distance = dual_node_position_distance(context, pivot, neighbor, dr);
    return *distance > 0.0
        && dual_node_polar_coordinates_from_displacement(
            context, pivot, neighbor, dr, *distance, cosphi, sinphi);
}

static void dual_node_initialize_angular_tolerance(
        dual_node_search_context *context)
{
    const real max_order = cballs_opt_edge_corrections(context->cmd)
        ? (real)(dual_node_window_orders(context->cmd) - 1)
        : (real)MAX(context->cmd->mChebyshev, 0);

    /* dual-node's LogMultipole default is pi / (2 * max_n + 1). */
    context->angular_tolerance =
        context->cmd->theta * PI / (2.0 * max_order + 1.0);
    context->max_angular_ratio = context->angular_tolerance >= 0.5 * PI
        ? 1.0 : rsin(MAX(0.0, context->angular_tolerance));
}

static int dual_node_node_orientation_status(
        const dual_node_search_context *context,
        const fcfc_ballnode *pivot, const fcfc_ballnode *q,
        const fcfc_ballnode *r, int *nq, int *nr,
        real *cosq, real *sinq, real *cosr, real *sinr)
{
    compute_vector unused;
    real dq = dual_node_position_distance(context, pivot->center, q->center,
                                         unused);
    real dr = dual_node_position_distance(context, pivot->center, r->center,
                                         unused);
    const real sq = (real)pivot->radius + (real)q->radius;
    const real sr = (real)pivot->radius + (real)r->radius;
    real ignored;

    if (dual_node_leg_outside(context, dq, sq)
        || dual_node_leg_outside(context, dr, sr))
        return DUAL_NODE_TRIPLE_OUTSIDE;
    if (!context->use_three_cells || !(dq > sq) || !(dr > sr))
        return DUAL_NODE_TRIPLE_SPLIT;

    *nq = dual_node_sloppy_bin(context, dq, sq);
    *nr = dual_node_sloppy_bin(context, dr, sr);
    if (*nq < 0 || *nr < 0) return DUAL_NODE_TRIPLE_SPLIT;

    compute_vector legq, legr;
    SUBV(legq, pivot->center, q->center);
    SUBV(legr, pivot->center, r->center);
    if (cballs_angular_extent(pivot->center, legq, pivot->radius, q->radius)
            > context->max_angular_ratio
        || cballs_angular_extent(pivot->center, legr, pivot->radius, r->radius)
            > context->max_angular_ratio)
        return DUAL_NODE_TRIPLE_SPLIT;
    if (!dual_node_polar_coordinates(context, pivot->center, q->center,
                                    &ignored, cosq, sinq)
        || !dual_node_polar_coordinates(context, pivot->center, r->center,
                                       &ignored, cosr, sinr))
        return DUAL_NODE_TRIPLE_SPLIT;
    return DUAL_NODE_TRIPLE_ACCEPT;
}

static void dual_node_accumulate_modes(dual_node_triple_histogram *hist,
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
    for (int order = 0; order < MAX(hist->orders, hist->window_orders); order++) {
        if (order < hist->orders) {
            dual_node_zeta_component(hist, DUAL_NODE_ZETA_COS, order)[index]
                += numerator * cosq * cosr;
            dual_node_zeta_component(hist, DUAL_NODE_ZETA_SIN, order)[index]
                += numerator * sinq * sinr;
            dual_node_zeta_component(hist, DUAL_NODE_ZETA_SINCOS, order)[index]
                += numerator * sinq * cosr;
            dual_node_zeta_component(hist, DUAL_NODE_ZETA_COSSIN, order)[index]
                += numerator * cosq * sinr;
        }
        if (order < hist->window_orders) {
            const size_t wi = (size_t)order * hist->plane + index;
            hist->window_re[wi] += denominator * (cosq * cosr + sinq * sinr);
            hist->window_im[wi] += denominator * (sinq * cosr - cosq * sinr);
        }
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

static inline real dual_node_node_field_sum(
        const dual_node_search_context *context, const fcfc_ballnode *ball_node)
{
    return context->weighted
        ? ball_node->weighted_kappa_sum : ball_node->kappa_sum;
}

static inline real dual_node_node_normalization_sum(
        const dual_node_search_context *context, const fcfc_ballnode *ball_node)
{
    return context->weighted ? ball_node->field_weight_sum
                             : (real)dual_node_node_count(ball_node);
}

static void dual_node_accumulate_node_orientation(
        const dual_node_search_context *context,
        dual_node_triple_histogram *hist,
        const fcfc_ballnode *pivot, const fcfc_ballnode *q,
        const fcfc_ballnode *r, int nq, int nr,
        real cosq, real sinq, real cosr, real sinr)
{
    const real numerator = dual_node_node_field_sum(context, pivot)
                         * dual_node_node_field_sum(context, q)
                         * dual_node_node_field_sum(context, r);
    const real denominator = dual_node_node_normalization_sum(context, pivot)
                           * dual_node_node_normalization_sum(context, q)
                           * dual_node_node_normalization_sum(context, r);

    dual_node_accumulate_modes(hist, nq, nr, numerator, denominator,
                              cosq, sinq, cosr, sinr);
    hist->cell_triples++;
}

static void dual_node_accumulate_body_orientation(
        const dual_node_search_context *context,
        dual_node_triple_histogram *hist, bodyptr pivot, bodyptr q, bodyptr r)
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

    if (!dual_node_polar_coordinates(context, Pos(pivot), Pos(q),
                                    &dq, &cosq, &sinq)
        || !dual_node_polar_coordinates(context, Pos(pivot), Pos(r),
                                       &dr, &cosr, &sinr))
        return;
    nq = dual_node_bin_index(context, dq);
    nr = dual_node_bin_index(context, dr);
    if (nq < 0 || nr < 0) return;

    fp = dual_node_body_field(pivot);
    fq = dual_node_body_field(q);
    fr = dual_node_body_field(r);
    if (context->weighted) {
        fp *= Weight(pivot);
        fq *= Weight(q);
        fr *= Weight(r);
        denominator = Weight(pivot) * Weight(q) * Weight(r);
    } else {
        denominator = 1.0;
    }
    dual_node_accumulate_modes(hist, nq, nr, fp * fq * fr, denominator,
                              cosq, sinq, cosr, sinr);
    hist->body_triples++;
}

static void dual_node_accumulate_body_mask(
        const dual_node_search_context *context,
        dual_node_triple_histogram *hist, bodyptr bodies[3], unsigned mask)
{
    for (int permutation = 0; permutation < 6; permutation++) {
        if (mask & (1U << permutation)) {
            const unsigned char *p = dual_node_permutations[permutation];
            dual_node_accumulate_body_orientation(
                context, hist, bodies[p[0]], bodies[p[1]], bodies[p[2]]);
        }
    }
}

#ifdef DUAL_NODE_TASK_FRONTIER_ENGINE
static unsigned dual_node_request_leg_split(
        const fcfc_ballnode *pivot, int pivot_slot,
        const fcfc_ballnode *neighbor, int neighbor_slot)
{
    const bool pivot_leaf = dual_node_node_is_leaf(pivot);
    const bool neighbor_leaf = dual_node_node_is_leaf(neighbor);

    if (pivot_leaf && neighbor_leaf) return 0U;
    if (pivot_leaf) return 1U << neighbor_slot;
    if (neighbor_leaf) return 1U << pivot_slot;
    return (real)pivot->radius >= (real)neighbor->radius
        ? 1U << pivot_slot : 1U << neighbor_slot;
}

/* Select the cells responsible for a failed two-ball leg criterion. */
static unsigned dual_node_orientation_split_mask(
        const dual_node_search_context *context,
        const fcfc_ballnode *cells[3], const unsigned char permutation[3])
{
    const int pivot_slot = permutation[0];
    const fcfc_ballnode *pivot = cells[pivot_slot];
    unsigned split_mask = 0U;

    for (int leg = 1; leg <= 2; leg++) {
        const int neighbor_slot = permutation[leg];
        const fcfc_ballnode *neighbor = cells[neighbor_slot];
        compute_vector unused;
        const real distance = dual_node_position_distance(
            context, pivot->center, neighbor->center, unused);
        const real size = (real)pivot->radius + (real)neighbor->radius;
        bool leg_needs_split = !context->use_three_cells
            || !(distance > size)
            || dual_node_sloppy_bin(context, distance, size) < 0;

        if (!leg_needs_split) {
            leg_needs_split = size / distance
                > context->max_angular_ratio;
        }
        if (leg_needs_split)
            split_mask |= dual_node_request_leg_split(
                pivot, pivot_slot, neighbor, neighbor_slot);
    }
    return split_mask;
}

static void dual_node_process111_children(
        const dual_node_search_context *context,
        const fcfc_balltreeptr trees[3], const INTEGER nodes[3],
        unsigned orientation_mask, unsigned split_mask, int slot,
        dual_node_triple_histogram *hist);
#endif

static void dual_node_process111(const dual_node_search_context *context,
                                const fcfc_balltreeptr trees[3],
                                const INTEGER nodes[3], unsigned mask,
                                dual_node_triple_histogram *hist)
{
    const fcfc_ballnode *cells[3] = {
        &trees[0]->nodes[nodes[0]],
        &trees[1]->nodes[nodes[1]],
        &trees[2]->nodes[nodes[2]]
    };
    unsigned remaining = mask;
#ifdef DUAL_NODE_TASK_FRONTIER_ENGINE
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
        p = dual_node_permutations[permutation];
        orientation_status = dual_node_node_orientation_status(
            context, cells[p[0]], cells[p[1]], cells[p[2]],
            &nq, &nr, &cosq, &sinq, &cosr, &sinr);
        if (orientation_status == DUAL_NODE_TRIPLE_OUTSIDE) {
            remaining &= ~bit;
        } else if (orientation_status == DUAL_NODE_TRIPLE_ACCEPT) {
            dual_node_accumulate_node_orientation(
                context, hist, cells[p[0]], cells[p[1]], cells[p[2]],
                nq, nr, cosq, sinq, cosr, sinr);
            remaining &= ~bit;
#ifdef DUAL_NODE_TASK_FRONTIER_ENGINE
        } else {
            requested_splits |= dual_node_orientation_split_mask(
                context, cells, p);
#endif
        }
    }
    if (remaining == 0) return;

    if (dual_node_node_is_leaf(cells[0])
        && dual_node_node_is_leaf(cells[1])
        && dual_node_node_is_leaf(cells[2])) {
        for (INTEGER i = cells[0]->first; i <= cells[0]->last; i++)
            for (INTEGER j = cells[1]->first; j <= cells[1]->last; j++)
                for (INTEGER k = cells[2]->first; k <= cells[2]->last; k++) {
                    bodyptr bodies[3] = {
                        trees[0]->bptr[i], trees[1]->bptr[j], trees[2]->bptr[k]
                    };
                    dual_node_accumulate_body_mask(context, hist, bodies, remaining);
                }
        return;
    }

    {
#ifdef DUAL_NODE_TASK_FRONTIER_ENGINE
        unsigned split_mask = 0U;
        real largest_requested = -1.0;

        for (int slot = 0; slot < 3; slot++) {
            if ((requested_splits & (1U << slot))
                && !dual_node_node_is_leaf(cells[slot]))
                largest_requested = MAX(
                    largest_requested, (real)cells[slot]->radius);
        }
        if (largest_requested >= 0.0) {
            for (int slot = 0; slot < 3; slot++) {
                if ((requested_splits & (1U << slot))
                    && !dual_node_node_is_leaf(cells[slot])
                    && (real)cells[slot]->radius
                        >= DUAL_NODE_SPLIT_FACTOR * largest_requested)
                    split_mask |= 1U << slot;
            }
        }
        if (split_mask != 0U) {
            dual_node_process111_children(
                context, trees, nodes, remaining, split_mask, 0, hist);
            return;
        }
#endif
        int split = -1;
        real largest = -1.0;
        for (int slot = 0; slot < 3; slot++) {
            if (!dual_node_node_is_leaf(cells[slot])
                && (real)cells[slot]->radius > largest) {
                split = slot;
                largest = (real)cells[slot]->radius;
            }
        }
        if (split >= 0) {
            INTEGER child_nodes[3] = {nodes[0], nodes[1], nodes[2]};
            child_nodes[split] = cells[split]->left;
            dual_node_process111(context, trees, child_nodes, remaining, hist);
            child_nodes[split] = cells[split]->right;
            dual_node_process111(context, trees, child_nodes, remaining, hist);
        }
    }
}

#ifdef DUAL_NODE_TASK_FRONTIER_ENGINE
static void dual_node_process111_children(
        const dual_node_search_context *context,
        const fcfc_balltreeptr trees[3], const INTEGER nodes[3],
        unsigned orientation_mask, unsigned split_mask, int slot,
        dual_node_triple_histogram *hist)
{
    if (slot == 3) {
        dual_node_process111(
            context, trees, nodes, orientation_mask, hist);
        return;
    }
    if (split_mask & (1U << slot)) {
        INTEGER child_nodes[3] = {nodes[0], nodes[1], nodes[2]};
        const fcfc_ballnode *ball_node =
            &trees[slot]->nodes[nodes[slot]];

        child_nodes[slot] = ball_node->left;
        dual_node_process111_children(
            context, trees, child_nodes, orientation_mask,
            split_mask, slot + 1, hist);
        child_nodes[slot] = ball_node->right;
        dual_node_process111_children(
            context, trees, child_nodes, orientation_mask,
            split_mask, slot + 1, hist);
    } else {
        dual_node_process111_children(
            context, trees, nodes, orientation_mask,
            split_mask, slot + 1, hist);
    }
}
#endif

static void dual_node_process21_auto(
        const dual_node_search_context *context, const fcfc_balltreeptr tree,
        INTEGER two_index, INTEGER one_index,
        dual_node_triple_histogram *hist)
{
    const fcfc_ballnode *two = &tree->nodes[two_index];
    const fcfc_ballnode *one = &tree->nodes[one_index];

    if (dual_node_node_count(two) < 2) return;
    {
        const real distance = dual_node_center_distance(
            context->cmd, context->gd, two, one);
        const real size = (real)two->radius + (real)one->radius;
        if (dual_node_leg_outside(context, distance, size)) return;
    }
    if (dual_node_node_is_leaf(two) && dual_node_node_is_leaf(one)) {
        for (INTEGER i = two->first; i < two->last; i++)
            for (INTEGER j = i + 1; j <= two->last; j++)
                for (INTEGER k = one->first; k <= one->last; k++) {
                    bodyptr bodies[3] = {
                        tree->bptr[i], tree->bptr[j], tree->bptr[k]
                    };
                    dual_node_accumulate_body_mask(context, hist, bodies, 0x3fU);
                }
        return;
    }

    if (!dual_node_node_is_leaf(one)
        && (dual_node_node_is_leaf(two) || one->radius > two->radius)) {
        dual_node_process21_auto(context, tree, two_index, one->left, hist);
        dual_node_process21_auto(context, tree, two_index, one->right, hist);
    } else {
        const fcfc_balltreeptr trees[3] = {tree, tree, tree};
        const INTEGER nodes[3] = {two->left, two->right, one_index};
        dual_node_process21_auto(context, tree, two->left, one_index, hist);
        dual_node_process21_auto(context, tree, two->right, one_index, hist);
        dual_node_process111(context, trees, nodes, 0x3fU, hist);
    }
}

static void dual_node_process3_auto(const dual_node_search_context *context,
                                   const fcfc_balltreeptr tree, INTEGER inode,
                                   dual_node_triple_histogram *hist)
{
    const fcfc_ballnode *ball_node = &tree->nodes[inode];

    if (dual_node_node_count(ball_node) < 3) return;
    if (2.0 * (real)ball_node->radius <= context->cmd->rminHist) return;
    if (dual_node_node_is_leaf(ball_node)) {
        for (INTEGER i = ball_node->first; i < ball_node->last - 1; i++)
            for (INTEGER j = i + 1; j < ball_node->last; j++)
                for (INTEGER k = j + 1; k <= ball_node->last; k++) {
                    bodyptr bodies[3] = {
                        tree->bptr[i], tree->bptr[j], tree->bptr[k]
                    };
                    dual_node_accumulate_body_mask(context, hist, bodies, 0x3fU);
                }
        return;
    }
    dual_node_process3_auto(context, tree, ball_node->left, hist);
    dual_node_process3_auto(context, tree, ball_node->right, hist);
    dual_node_process21_auto(
        context, tree, ball_node->left, ball_node->right, hist);
    dual_node_process21_auto(
        context, tree, ball_node->right, ball_node->left, hist);
}

static void dual_node_process12_cross(
        const dual_node_search_context *context,
        const fcfc_balltreeptr pivot_tree, INTEGER pivot_index,
        const fcfc_balltreeptr pair_tree, INTEGER pair_index,
        dual_node_triple_histogram *hist)
{
    const fcfc_ballnode *pivot = &pivot_tree->nodes[pivot_index];
    const fcfc_ballnode *pair = &pair_tree->nodes[pair_index];

    if (dual_node_node_count(pair) < 2) return;
    {
        const real distance = dual_node_center_distance(
            context->cmd, context->gd, pivot, pair);
        const real size = (real)pivot->radius + (real)pair->radius;
        if (dual_node_leg_outside(context, distance, size)) return;
    }
    if (dual_node_node_is_leaf(pivot) && dual_node_node_is_leaf(pair)) {
        for (INTEGER i = pivot->first; i <= pivot->last; i++)
            for (INTEGER j = pair->first; j < pair->last; j++)
                for (INTEGER k = j + 1; k <= pair->last; k++) {
                    bodyptr bodies[3] = {
                        pivot_tree->bptr[i], pair_tree->bptr[j], pair_tree->bptr[k]
                    };
                    dual_node_accumulate_body_mask(context, hist, bodies, 0x03U);
                }
        return;
    }

    if (!dual_node_node_is_leaf(pivot)
        && (dual_node_node_is_leaf(pair) || pivot->radius > pair->radius)) {
        dual_node_process12_cross(context, pivot_tree, pivot->left,
                                 pair_tree, pair_index, hist);
        dual_node_process12_cross(context, pivot_tree, pivot->right,
                                 pair_tree, pair_index, hist);
    } else {
        const fcfc_balltreeptr trees[3] = {
            pivot_tree, pair_tree, pair_tree
        };
        const INTEGER nodes[3] = {pivot_index, pair->left, pair->right};
        dual_node_process12_cross(context, pivot_tree, pivot_index,
                                 pair_tree, pair->left, hist);
        dual_node_process12_cross(context, pivot_tree, pivot_index,
                                 pair_tree, pair->right, hist);
        dual_node_process111(context, trees, nodes, 0x03U, hist);
    }
}

static int dual_node_allocate_triple_histograms(
        struct cmdline_data *cmd, INTEGER task_count, size_t stride,
        int orders, real **histograms, size_t *values_per_task,
        INTEGER **body_triples, INTEGER **cell_triples)
{
    size_t tasks;
    size_t values;

    *histograms = NULL;
    *values_per_task = 0;
    *body_triples = NULL;
    *cell_triples = NULL;
    if (task_count <= 0
        || !dual_node_triple_values(stride, orders, dual_node_window_orders(cmd),
                                   values_per_task))
        goto invalid_size;
    tasks = (size_t)task_count;
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
                 DUAL_NODE_METHOD_NAME);
        return FAILURE;
    }
    return SUCCESS;

invalid_size:
    snprintf(cmd->error_message, _ERRORMSGSIZE_,
             "%s: triple histogram size overflow", DUAL_NODE_METHOD_NAME);
    return FAILURE;
}

#include "dual_node_task_schedule.h"

#endif /* THREEPCFCONVERGENCE */

#ifdef TWOPCF
static int dual_node_allocate_task_histograms(
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
                 "%s: task histogram allocation failed", DUAL_NODE_METHOD_NAME);
        return FAILURE;
    }
    return SUCCESS;

invalid_size:
    snprintf(cmd->error_message, _ERRORMSGSIZE_,
             "%s: task histogram size overflow", DUAL_NODE_METHOD_NAME);
    return FAILURE;
}

static int dual_node_run_pair_tasks(
        const dual_node_search_context *context,
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
    double phase_started = 0.0;

    if (tree1 == NULL || tree2 == NULL
        || tree1->packed_points == NULL || tree2->packed_points == NULL) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: compact point storage is unavailable",
                 DUAL_NODE_METHOD_NAME);
        goto cleanup;
    }

    allocation_status = dual_node_allocate_task_histograms(
        cmd, frontier_count1, stride, &task_histograms,
        &task_body_counts, &task_cell_counts);
    if (dual_node_distributed_consensus(
            cmd, allocation_status,
            "two-ball pair-task allocation") == FAILURE)
        goto cleanup;

    if (context->profile) phase_started = dual_node_timer_now();
#pragma omp parallel for schedule(dynamic,1)
    for (INTEGER itask = 0; itask < frontier_count1; itask++) {
        real *base = task_histograms + (size_t)itask * 3 * stride;
        dual_node_histogram hist;

        if (!dual_node_distributed_task_owned(itask)) continue;

        hist.pair_count = base;
        hist.weight_product = base + stride;
        hist.field_product = base + 2 * stride;
        hist.body_pairs = 0;
        hist.cell_pairs = 0;

        if (auto_correlation) {
            dual_node_process_auto(
                context, tree1, frontier1[itask], &hist);
            for (INTEGER jtask = itask + 1;
                 jtask < frontier_count2; jtask++)
                dual_node_process_pair(
                    context, tree1, frontier1[itask],
                    tree2, frontier2[jtask], &hist);
        } else {
            for (INTEGER jtask = 0; jtask < frontier_count2; jtask++)
                dual_node_process_pair(
                    context, tree1, frontier1[itask],
                    tree2, frontier2[jtask], &hist);
        }
        task_body_counts[itask] = hist.body_pairs;
        task_cell_counts[itask] = hist.cell_pairs;
    }
    if (context->profile)
        context->timers->pair_traversal +=
            dual_node_timer_now() - phase_started;

    if (context->profile) phase_started = dual_node_timer_now();
    if (dual_node_distributed_reduce_reals(
            cmd, task_histograms,
            (size_t)frontier_count1 * 3 * stride) == FAILURE)
        reduction_status = FAILURE;
    if (dual_node_distributed_reduce_integers(
            cmd, task_body_counts, (size_t)frontier_count1) == FAILURE)
        reduction_status = FAILURE;
    if (dual_node_distributed_reduce_integers(
            cmd, task_cell_counts, (size_t)frontier_count1) == FAILURE)
        reduction_status = FAILURE;
    if (dual_node_distributed_consensus(
            cmd, reduction_status,
            "two-ball pair-task reduction") == FAILURE)
        goto cleanup;
    if (!dual_node_distributed_publish()) {
        if (context->profile)
            context->timers->reduction +=
                dual_node_timer_now() - phase_started;
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

    gd->histogram_products |= CBALLS_PRODUCT_NN | CBALLS_PRODUCT_XI;
    if (cballs_opt_compute_histn(cmd) && cballs_opt_and_cf(cmd)) {
        for (int n = 1; n <= cmd->sizeHistN; n++) gd->histNN[n] *= 2.0;
        if (search_compute_HistN(cmd, gd, tree1->npoint) == FAILURE)
            goto cleanup;
    }
    if (context->profile)
        context->timers->reduction += dual_node_timer_now() - phase_started;
    status = SUCCESS;

cleanup:
    status = dual_node_distributed_consensus(
        cmd, status, "two-ball pair publication");
    free(task_cell_counts);
    free(task_body_counts);
    free(task_histograms);
    return status;
}
#endif /* TWOPCF */

#ifndef DUAL_NODE_TASK_FRONTIER_ENGINE
global int searchcalc_balltree_2balls_omp(
        struct cmdline_data *cmd, struct global_data *gd,
        bodyptr *btab, INTEGER *nbody, INTEGER ipmin, INTEGER *ipmax,
        int cat1, int cat2)
{
    const bool catalog_auto_correlation = cat1 == cat2;
    const bool auto_correlation = catalog_auto_correlation
        && !DUAL_NODE_SEPARATE_PIVOT_TREE(cmd);
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
    dual_node_search_context context = {0};
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
        ? dual_node_triple_frontier_target(cmd)
        : dual_node_pair_frontier_target(cmd);
    int operation_status;
    int status = FAILURE;
    const double cpustart = CPUTIME;

    gd->cpu_edge_correction = 0.0;
    gd->wall_edge_correction = 0.0;

#if !defined(TWOPCF) && !defined(THREEPCFCONVERGENCE)
    snprintf(cmd->error_message, _ERRORMSGSIZE_,
             "%s was built with TWOPCFON=0 and TPCFON=0",
             DUAL_NODE_METHOD_NAME);
    return FAILURE;
#endif

    if (only_2pcf && only_3pcf) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: only-2pcf and only-3pcf are mutually exclusive",
                 DUAL_NODE_METHOD_NAME);
        return FAILURE;
    }
    if (only_2pcf && !run_2pcf) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: only-2pcf requires TWOPCFON=1", DUAL_NODE_METHOD_NAME);
        return FAILURE;
    }
    if (only_3pcf && !run_3pcf) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: only-3pcf requires TPCFON=1", DUAL_NODE_METHOD_NAME);
        return FAILURE;
    }

    if (ipmin != 1 || ipmax[cat1] != nbody[cat1]) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s requires the complete pivot catalog",
                 DUAL_NODE_METHOD_NAME);
        return FAILURE;
    }
    if (cmd->nsmooth <= 0) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s requires nsmooth > 0", DUAL_NODE_METHOD_NAME);
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
#ifdef DUAL_NODE_DISTRIBUTED_ENGINE
    verb_print(cmd->verbose,
               "%s: %d ranks with deterministic cyclic frontier ownership\n",
               DUAL_NODE_METHOD_NAME, DUAL_NODE_DISTRIBUTED_SIZE());
#endif
#ifdef DUAL_NODE_SCAN_LEVEL_FRONTIER
    verb_print(cmd->verbose,
               "balanced scan-level task frontier enabled by BALLS4SCANLEV\n");
#endif
#ifdef THREEPCFCONVERGENCE
    if (run_3pcf && nbody[cat1] >= DUAL_NODE_LARGE_DIRECT_3PCF)
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
                   "using dual-node-style pair/triplet-weight normalization\n");

#ifdef OPENMPCODE
    ThreadCount(cmd, gd, nbody[cat1], cat1);
#endif

    operation_status = DUAL_NODE_PREPARE_PIVOTS(
        cmd, gd, btab, nbody, ipmin, ipmax, cat1, cat2);
    if (dual_node_distributed_consensus(
            cmd, operation_status,
            "two-ball pivot preparation") == FAILURE)
        goto cleanup;
    operation_status = search_init_gd_hist_sincos(cmd, gd);
    if (dual_node_distributed_consensus(
            cmd, operation_status,
            "two-ball histogram initialization") == FAILURE)
        goto cleanup;
    operation_status = DUAL_NODE_BUILD_PIVOT_TREE(
        cmd, gd, btab[cat1], nbody[cat1], cmd->nsmooth, &tree1);
    if (dual_node_distributed_consensus(
            cmd, operation_status,
            "two-ball first tree construction") == FAILURE)
        goto cleanup;
    DUAL_NODE_PUBLISH_NODE_COUNT(gd, cat1, tree1->nnode);

    if (auto_correlation) {
        tree2 = tree1;
    } else {
        operation_status = DUAL_NODE_BUILD_NEIGHBOR_TREE(
            cmd, gd, btab[cat2], nbody[cat2], cmd->nsmooth, &tree2);
        if (dual_node_distributed_consensus(
                cmd, operation_status,
                "two-ball second tree construction") == FAILURE)
            goto cleanup;
        DUAL_NODE_PUBLISH_NODE_COUNT(gd, cat2, tree2->nnode);
    }

    operation_status = fcfc_balltree_frontier(
        cmd, tree1, target_tasks, &frontier1, &frontier_count1);
    if (dual_node_distributed_consensus(
            cmd, operation_status,
            "two-ball pivot-frontier construction") == FAILURE)
        goto cleanup;
    if (auto_correlation) {
        frontier2 = frontier1;
        frontier_count2 = frontier_count1;
    } else {
        operation_status = fcfc_balltree_frontier(
            cmd, tree2, target_tasks, &frontier2, &frontier_count2);
        if (dual_node_distributed_consensus(
                cmd, operation_status,
                "two-ball neighbor-frontier construction") == FAILURE)
            goto cleanup;
    }

    dual_node_initialize_radial_context(&context, cmd, gd);
    context.use_two_balls = !cballs_opt_no_two_balls(cmd);
    context.use_bin_theta = scanopt(cmd->options, "dual-node-bin-theta");
    context.use_three_cells = !cballs_opt_no_two_balls(cmd);
    context.weighted = DUAL_NODE_CONTEXT_WEIGHTED(cmd);
#ifdef THREEPCFCONVERGENCE
    dual_node_initialize_angular_tolerance(&context);
#endif

#ifdef TWOPCF
    if (run_2pcf) {
        if (dual_node_run_pair_tasks(
                &context, tree1, tree2, auto_correlation,
                frontier1, frontier_count1,
                frontier2, frontier_count2) == FAILURE)
            goto cleanup;
    }
#endif /* TWOPCF */

#ifdef THREEPCFCONVERGENCE
    if (run_3pcf) {
    operation_status = dual_node_allocate_triple_histograms(
        cmd, frontier_count1, stride, triple_orders,
        &triple_task_histograms, &triple_values_per_task,
        &triple_task_body_counts, &triple_task_cell_counts);
    if (dual_node_distributed_consensus(
            cmd, operation_status,
            "two-ball triple histogram allocation") == FAILURE)
        goto cleanup;

#pragma omp parallel for schedule(dynamic,1)
    for (INTEGER itask = 0; itask < frontier_count1; itask++) {
        real *base = triple_task_histograms
                   + (size_t)itask * triple_values_per_task;
        dual_node_triple_histogram hist;

        dual_node_initialize_triple_histogram(
            &hist, base, stride, triple_orders, dual_node_window_orders(cmd));

        if (!dual_node_distributed_task_owned(itask)) continue;

        if (auto_correlation) {
            const fcfc_balltreeptr trees[3] = {tree1, tree1, tree1};

            dual_node_process3_auto(&context, tree1, frontier1[itask], &hist);
            for (INTEGER jtask = itask + 1;
                 jtask < frontier_count1; jtask++) {
                dual_node_process21_auto(&context, tree1,
                                        frontier1[itask], frontier1[jtask],
                                        &hist);
                dual_node_process21_auto(&context, tree1,
                                        frontier1[jtask], frontier1[itask],
                                        &hist);
                for (INTEGER ktask = jtask + 1;
                     ktask < frontier_count1; ktask++) {
                    const INTEGER nodes[3] = {
                        frontier1[itask], frontier1[jtask], frontier1[ktask]
                    };
                    dual_node_process111(&context, trees, nodes, 0x3fU, &hist);
                }
            }
        } else {
            const fcfc_balltreeptr trees[3] = {tree1, tree2, tree2};

            for (INTEGER jtask = 0; jtask < frontier_count2; jtask++) {
                dual_node_process12_cross(
                    &context, tree1, frontier1[itask],
                    tree2, frontier2[jtask], &hist);
                for (INTEGER ktask = jtask + 1;
                     ktask < frontier_count2; ktask++) {
                    const INTEGER nodes[3] = {
                        frontier1[itask], frontier2[jtask], frontier2[ktask]
                    };
                    dual_node_process111(&context, trees, nodes, 0x03U, &hist);
                }
            }
        }
        triple_task_body_counts[itask] = hist.body_triples;
        triple_task_cell_counts[itask] = hist.cell_triples;
    }

    if (dual_node_distributed_reduce_reals(
            cmd, triple_task_histograms,
            (size_t)frontier_count1 * triple_values_per_task) == FAILURE)
        triple_reduction_status = FAILURE;
    if (dual_node_distributed_reduce_integers(
            cmd, triple_task_body_counts,
            (size_t)frontier_count1) == FAILURE)
        triple_reduction_status = FAILURE;
    if (dual_node_distributed_reduce_integers(
            cmd, triple_task_cell_counts,
            (size_t)frontier_count1) == FAILURE)
        triple_reduction_status = FAILURE;
    if (dual_node_distributed_consensus(
            cmd, triple_reduction_status,
            "two-ball triple histogram reduction") == FAILURE)
        goto cleanup;

    if (dual_node_distributed_publish()) {
    for (INTEGER itask = 0; itask < frontier_count1; itask++) {
        const real *base = triple_task_histograms
                         + (size_t)itask * triple_values_per_task;
        const size_t plane = stride * stride;

        for (int m = 0; m < triple_orders; m++) {
            const real *zcos = base
                + ((size_t)DUAL_NODE_ZETA_COS * (size_t)triple_orders
                   + (size_t)m)
                * plane;
            const real *zsin = base
                + ((size_t)DUAL_NODE_ZETA_SIN * (size_t)triple_orders
                   + (size_t)m)
                * plane;
            const real *zsincos = base
                + ((size_t)DUAL_NODE_ZETA_SINCOS * (size_t)triple_orders
                   + (size_t)m)
                * plane;
            const real *zcossin = base
                + ((size_t)DUAL_NODE_ZETA_COSSIN * (size_t)triple_orders
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

    }
    operation_status = dual_node_distributed_publish()
        ? dual_node_publish_edge(cmd, gd, triple_task_histograms, frontier_count1,
                                triple_values_per_task, stride, triple_orders)
        : SUCCESS;
    if (dual_node_distributed_consensus(
            cmd, operation_status, "two-ball edge correction") == FAILURE)
        goto cleanup;
    if (dual_node_distributed_publish()) {
    if (dual_node_normalize_3pcf(cmd)) {
        for (int n1 = 1; n1 <= cmd->sizeHistN; n1++)
            for (int n2 = 1; n2 <= cmd->sizeHistN; n2++) {
            const size_t index = (size_t)n1 * stride + (size_t)n2;
            real denominator = 0.0;
            for (INTEGER itask = 0; itask < frontier_count1; itask++) {
                const real *base = triple_task_histograms
                    + (size_t)itask * triple_values_per_task;
                denominator += base[
                    DUAL_NODE_ZETA_COMPONENTS
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
    if (tree2 != tree1) DUAL_NODE_RELEASE_TREE(tree2);
    DUAL_NODE_RELEASE_TREE(tree1);
    return status;
}
#endif /* !DUAL_NODE_TASK_FRONTIER_ENGINE */

#ifdef DUAL_NODE_TASK_FRONTIER_ENGINE
#ifdef DUAL_NODE_LOG_MULTIPOLE_ENGINE
static int dual_node_search_direct_triples(
#else
global int searchcalc_balltree_2balls_omp(
#endif
        struct cmdline_data *cmd, struct global_data *gd,
        bodyptr *btab, INTEGER *nbody, INTEGER ipmin, INTEGER *ipmax,
        int cat1, int cat2)
{
#ifdef THREEPCFCONVERGENCE
    const bool catalog_auto_correlation = cat1 == cat2;
    const bool auto_correlation = catalog_auto_correlation
        && !DUAL_NODE_SEPARATE_PIVOT_TREE(cmd);
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
    const INTEGER target_tasks = dual_node_task_target(cmd, stride, triple_orders);
    const int leaf_capacity = scanopt(cmd->options, "dual-node-bucket-leaves")
        ? cmd->nsmooth : 1;
    dual_node_search_context context = {0};
    fcfc_balltreeptr tree1 = NULL;
    fcfc_balltreeptr tree2 = NULL;
    dual_node_triple_task *tasks = NULL;
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
    int operation_status;
    int reduction_status = SUCCESS;
    int status = FAILURE;
    const double cpustart = CPUTIME;

    gd->cpu_edge_correction = 0.0;
    gd->wall_edge_correction = 0.0;

    if (only_2pcf && only_3pcf) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: only-2pcf and only-3pcf are mutually exclusive",
                 DUAL_NODE_METHOD_NAME);
        return FAILURE;
    }
    if (only_2pcf && !run_2pcf) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: only-2pcf requires TWOPCFON=1",
                 DUAL_NODE_METHOD_NAME);
        return FAILURE;
    }
    if (ipmin != 1 || ipmax[cat1] != nbody[cat1]) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s requires the complete pivot catalog",
                 DUAL_NODE_METHOD_NAME);
        return FAILURE;
    }
    if (cmd->nsmooth <= 0) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s requires nsmooth > 0", DUAL_NODE_METHOD_NAME);
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
#ifdef DUAL_NODE_NATIVE_BINARY_VIEW
    verb_print(cmd->verbose,
               "native-octree binary view with leaf capacity %d\n",
               leaf_capacity);
#else
    verb_print(cmd->verbose,
               "ball-tree leaf capacity = %d%s\n", leaf_capacity,
               leaf_capacity == 1 ? " (dual-node singleton policy)"
                                  : " (nsmooth bucket policy)");
#endif
    if (cballs_opt_weights_norm(cmd))
        verb_print(cmd->verbose,
                   "using dual-node-style triplet-weight normalization\n");

#ifdef OPENMPCODE
    ThreadCount(cmd, gd, nbody[cat1], cat1);
#endif

    operation_status = DUAL_NODE_PREPARE_PIVOTS(
        cmd, gd, btab, nbody, ipmin, ipmax, cat1, cat2);
    if (dual_node_distributed_consensus(
            cmd, operation_status,
            "direct-triple pivot preparation") == FAILURE)
        goto cleanup;
    operation_status = search_init_gd_hist_sincos(cmd, gd);
    if (dual_node_distributed_consensus(
            cmd, operation_status,
            "direct-triple histogram initialization") == FAILURE)
        goto cleanup;
    operation_status = DUAL_NODE_BUILD_PIVOT_TREE(
        cmd, gd, btab[cat1], nbody[cat1], leaf_capacity, &tree1);
    if (dual_node_distributed_consensus(
            cmd, operation_status,
            "direct-triple pivot-tree construction") == FAILURE)
        goto cleanup;
    DUAL_NODE_PUBLISH_NODE_COUNT(gd, cat1, tree1->nnode);

    if (auto_correlation) {
        tree2 = tree1;
    } else {
        operation_status = DUAL_NODE_BUILD_NEIGHBOR_TREE(
            cmd, gd, btab[cat2], nbody[cat2], leaf_capacity, &tree2);
        if (dual_node_distributed_consensus(
                cmd, operation_status,
                "direct-triple neighbor-tree construction") == FAILURE)
            goto cleanup;
        DUAL_NODE_PUBLISH_NODE_COUNT(gd, cat2, tree2->nnode);
    }

    dual_node_initialize_radial_context(&context, cmd, gd);
    context.use_two_balls = !cballs_opt_no_two_balls(cmd);
    context.use_bin_theta = scanopt(cmd->options, "dual-node-bin-theta");
    context.use_three_cells = !cballs_opt_no_two_balls(cmd);
    context.weighted = DUAL_NODE_CONTEXT_WEIGHTED(cmd);
    dual_node_initialize_angular_tolerance(&context);

 #ifdef TWOPCF
    if (run_2pcf) {
        operation_status = fcfc_balltree_frontier(
            cmd, tree1, dual_node_pair_frontier_target(cmd),
            &pair_frontier1, &pair_frontier_count1);
        if (dual_node_distributed_consensus(
                cmd, operation_status,
                "direct-triple pair pivot-frontier construction") == FAILURE)
            goto cleanup;
        if (auto_correlation) {
            pair_frontier2 = pair_frontier1;
            pair_frontier_count2 = pair_frontier_count1;
        } else {
            operation_status = fcfc_balltree_frontier(
                cmd, tree2, dual_node_pair_frontier_target(cmd),
                &pair_frontier2, &pair_frontier_count2);
            if (dual_node_distributed_consensus(
                    cmd, operation_status,
                    "direct-triple pair neighbor-frontier construction")
                == FAILURE)
                goto cleanup;
        }
        if (dual_node_run_pair_tasks(
                &context, tree1, tree2, auto_correlation,
                pair_frontier1, pair_frontier_count1,
                pair_frontier2, pair_frontier_count2) == FAILURE)
            goto cleanup;
    }
#endif

    if (run_3pcf) {
    operation_status = dual_node_build_task_frontier(
        cmd, tree1, tree2, auto_correlation, target_tasks,
        &tasks, &task_count);
    if (dual_node_distributed_consensus(
            cmd, operation_status,
            "direct-triple task-frontier construction") == FAILURE)
        goto cleanup;
    operation_status = dual_node_allocate_triple_histograms(
        cmd, task_count, stride, triple_orders,
        &task_histograms, &values_per_task,
        &task_body_counts, &task_cell_counts);
    if (dual_node_distributed_consensus(
            cmd, operation_status,
            "direct-triple task-histogram allocation") == FAILURE)
        goto cleanup;

#pragma omp parallel for schedule(dynamic,1)
    for (INTEGER itask = 0; itask < task_count; itask++) {
        real *base = task_histograms
                   + (size_t)itask * values_per_task;
        dual_node_triple_histogram hist;

        if (!dual_node_distributed_task_owned(itask)) continue;

        dual_node_initialize_triple_histogram(
            &hist, base, stride, triple_orders, dual_node_window_orders(cmd));
        dual_node_run_task(&context, &tasks[itask], &hist);
        task_body_counts[itask] = hist.body_triples;
        task_cell_counts[itask] = hist.cell_triples;
    }

    if (dual_node_distributed_reduce_reals(
            cmd, task_histograms,
            (size_t)task_count * values_per_task) == FAILURE)
        reduction_status = FAILURE;
    if (dual_node_distributed_reduce_integers(
            cmd, task_body_counts, (size_t)task_count) == FAILURE)
        reduction_status = FAILURE;
    if (dual_node_distributed_reduce_integers(
            cmd, task_cell_counts, (size_t)task_count) == FAILURE)
        reduction_status = FAILURE;
    if (dual_node_distributed_consensus(
            cmd, reduction_status,
            "direct-triple task reduction") == FAILURE)
        goto cleanup;

    if (dual_node_distributed_publish()) {
    for (INTEGER itask = 0; itask < task_count; itask++) {
        const real *base = task_histograms
                         + (size_t)itask * values_per_task;
        const size_t plane = stride * stride;

        for (int m = 0; m < triple_orders; m++) {
            const real *zcos = base
                + ((size_t)DUAL_NODE_ZETA_COS * (size_t)triple_orders
                   + (size_t)m)
                * plane;
            const real *zsin = base
                + ((size_t)DUAL_NODE_ZETA_SIN * (size_t)triple_orders
                   + (size_t)m)
                * plane;
            const real *zsincos = base
                + ((size_t)DUAL_NODE_ZETA_SINCOS * (size_t)triple_orders
                   + (size_t)m)
                * plane;
            const real *zcossin = base
                + ((size_t)DUAL_NODE_ZETA_COSSIN * (size_t)triple_orders
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
    }

    operation_status = dual_node_distributed_publish()
        ? dual_node_publish_edge(cmd, gd, task_histograms, task_count,
                                values_per_task, stride, triple_orders)
        : SUCCESS;
    if (dual_node_distributed_consensus(
            cmd, operation_status,
            "direct-triple edge correction") == FAILURE)
        goto cleanup;
    if (dual_node_distributed_publish() && dual_node_normalize_3pcf(cmd)) {
        for (int n1 = 1; n1 <= cmd->sizeHistN; n1++)
            for (int n2 = 1; n2 <= cmd->sizeHistN; n2++) {
            const size_t index = (size_t)n1 * stride + (size_t)n2;
            real denominator = 0.0;

            for (INTEGER itask = 0; itask < task_count; itask++) {
                const real *base = task_histograms
                    + (size_t)itask * values_per_task;
                denominator += base[
                    DUAL_NODE_ZETA_COMPONENTS
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
    if (tree2 != tree1) DUAL_NODE_RELEASE_TREE(tree2);
    DUAL_NODE_RELEASE_TREE(tree1);
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
             "%s requires TPCFON=1", DUAL_NODE_METHOD_NAME);
    return FAILURE;
#endif /* THREEPCFCONVERGENCE */
}
#endif /* DUAL_NODE_TASK_FRONTIER_ENGINE */

#ifdef BALLTREE_2BALLS_SEARCH_FUNCTION
#undef searchcalc_balltree_2balls_omp

global int BALLTREE_2BALLS_SEARCH_FUNCTION(
        struct cmdline_data *cmd, struct global_data *gd,
        bodyptr *btab, INTEGER *nbody, INTEGER ipmin, INTEGER *ipmax,
        int cat1, int cat2)
{
    /* The legacy ball kernel has no scalar angular-window solver. */
    if (cballs_opt_legacy_one_ball(cmd) && cballs_opt_edge_corrections(cmd)) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: legacy-one-ball does not support edge-corrections; "
                 "remove legacy-one-ball to use the two-ball window solver",
                 cmd->searchMethod);
        return FAILURE;
    }
#ifdef DUAL_NODE_DISTRIBUTED_ENGINE
    if (cballs_opt_legacy_one_ball(cmd)) {
#ifdef BALLTREE_2BALLS_LEGACY_FUNCTION
        if (cballs_opt_no_two_balls(cmd)
            || scanopt(cmd->options, "dual-node-bin-theta")
            || scanopt(cmd->options, "dual-node-direct-triples")) {
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     "%s: legacy-one-ball cannot be combined with "
                     "no-two-balls, dual-node-bin-theta, or "
                     "dual-node-direct-triples",
                     cmd->searchMethod);
            return FAILURE;
        }
        verb_print(cmd->verbose,
                   "%s: dispatching to the distributed balltree legacy "
                   "kernel\n",
                   cmd->searchMethod);
        return BALLTREE_2BALLS_LEGACY_FUNCTION(
            cmd, gd, btab, nbody, ipmin, ipmax, cat1, cat2);
#else
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s does not support legacy-one-ball; use "
                 "search=balltree-mpi for the distributed legacy kernel",
                 cmd->searchMethod);
        return FAILURE;
#endif
    }
#else
    if (cballs_opt_legacy_one_ball(cmd)) {
        if (cballs_opt_no_two_balls(cmd)
            || scanopt(cmd->options, "dual-node-bin-theta")
            || scanopt(cmd->options, "dual-node-direct-triples")) {
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     "%s: legacy-one-ball cannot be combined with "
                     "no-two-balls, dual-node-bin-theta, or "
                     "dual-node-direct-triples",
                     cmd->searchMethod);
            return FAILURE;
        }
        verb_print(cmd->verbose,
                   "%s: dispatching to the balltree-omp legacy kernel\n",
                   cmd->searchMethod);
        return BALLTREE_2BALLS_LEGACY_FUNCTION(
            cmd, gd, btab, nbody, ipmin, ipmax, cat1, cat2);
    }
#endif
#ifdef SMOOTHPIVOT
    if (cballs_opt_smooth_pivot(cmd)
        && scanopt(cmd->options, "dual-node-direct-triples")) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: dual-node-direct-triples does not support smooth-pivot; "
                 "use no-smooth-pivot for that validation oracle",
                 cmd->searchMethod);
        return FAILURE;
    }
#endif
    return BALLTREE_2BALLS_FULL_FUNCTION(
        cmd, gd, btab, nbody, ipmin, ipmax, cat1, cat2);
}
#endif
