/*
 * Shared deterministic 2PCF traversal over a compact binary view of the
 * production cTreeBalls octree.
 *
 * The dual-node split and bin-slop criteria are adapted from dual-node:
 * Copyright (c) 2003-2024 Mike Jarvis, under its BSD-style license.
 */

#include <limits.h>
#include <stdint.h>
#include <stdlib.h>

#ifdef OPENMPCODE
#include <omp.h>
#endif

#include "globaldefs.h"
#include "native_octree_pair.h"
#include "octree_2balls_tree.h"

#define NATIVE_PAIR_SPLIT_FACTOR ((real)0.585)
#define NATIVE_PAIR_MAX_FRONTIER ((INTEGER)256)
#define NATIVE_PAIR_BATCH_SIZE 256

#if defined(__APPLE__)
extern void vvlog(double *, const double *, const int *);
extern void vvlogf(float *, const float *, const int *);
#elif defined(CBALLS_HAVE_SLEEF)
#include <sleef.h>
#endif

#if defined(__APPLE__) || defined(CBALLS_HAVE_SLEEF) \
    || defined(CBALLS_COMPILER_VECTOR_LOG)
#define NATIVE_PAIR_VECTOR_LOG 1
#endif

#ifdef TWOPCF
typedef struct {
    real *pair_count;
    real *weight_product;
    real *field_product;
    INTEGER body_pairs;
    INTEGER cell_pairs;
#ifdef NATIVE_PAIR_VECTOR_LOG
    real batch_distance2[NATIVE_PAIR_BATCH_SIZE];
    real batch_denominator[NATIVE_PAIR_BATCH_SIZE];
    real batch_numerator[NATIVE_PAIR_BATCH_SIZE];
    int batch_count;
#endif
} native_pair_histogram;

typedef struct {
    double build;
    double frontier;
    double traversal;
    double reduction;
} native_pair_phase_timers;

typedef struct {
    struct cmdline_data *cmd;
    struct global_data *gd;
    bool use_two_balls;
    bool use_bin_slop;
    bool weighted_signal;
    bool weighted_normalization;
    real pair_scale;
    real minimum2;
    real maximum2;
    real theta2;
    real bin_size;
    real inverse_bin_size;
    real bin_slop;
    real bin_slop2;
    real half_bin_plus_slop2;
    real log_acceptance_limit2;
    bool profile;
    native_pair_phase_timers *timers;
} native_pair_context;

static inline double native_pair_timer_now(void)
{
#ifdef OPENMPCODE
    return omp_get_wtime();
#else
    return CPUTIME;
#endif
}

static const char *native_pair_method(const struct cmdline_data *cmd)
{
    return cmd->searchMethod != NULL ? cmd->searchMethod
                                     : "native-octree-pair";
}

static bool native_pair_task_owned(
        const cballs_native_pair_policy *policy,
        struct cmdline_data *cmd, INTEGER task)
{
    const cballs_native_pair_parallel *parallel = policy->parallel;

    return parallel == NULL || parallel->task_owned == NULL
        || parallel->task_owned(cmd, task);
}

static bool native_pair_publish(
        const cballs_native_pair_policy *policy,
        struct cmdline_data *cmd)
{
    const cballs_native_pair_parallel *parallel = policy->parallel;

    return parallel == NULL || parallel->publish == NULL
        || parallel->publish(cmd);
}

static int native_pair_consensus(
        const cballs_native_pair_policy *policy,
        struct cmdline_data *cmd, int status, const char *operation)
{
    const cballs_native_pair_parallel *parallel = policy->parallel;

    return parallel != NULL && parallel->consensus != NULL
        ? parallel->consensus(cmd, status, operation) : status;
}

static int native_pair_reduce_reals(
        const cballs_native_pair_policy *policy,
        struct cmdline_data *cmd, real *values, size_t count)
{
    const cballs_native_pair_parallel *parallel = policy->parallel;

    return parallel != NULL && parallel->reduce_reals != NULL
        ? parallel->reduce_reals(cmd, values, count) : SUCCESS;
}

static int native_pair_reduce_integers(
        const cballs_native_pair_policy *policy,
        struct cmdline_data *cmd, INTEGER *values, size_t count)
{
    const cballs_native_pair_parallel *parallel = policy->parallel;

    return parallel != NULL && parallel->reduce_integers != NULL
        ? parallel->reduce_integers(cmd, values, count) : SUCCESS;
}

static INTEGER native_pair_frontier_target(
        const struct cmdline_data *cmd,
        const cballs_native_pair_policy *policy)
{
    INTEGER workers = 1;
    INTEGER ranks = 1;
    INTEGER target;

#ifdef OPENMPCODE
    if (cmd->numthreads > 0) workers = (INTEGER)cmd->numthreads;
#else
    (void)cmd;
#endif
    if (policy->parallel != NULL && policy->parallel->rank_count > 1)
        ranks = (INTEGER)policy->parallel->rank_count;
    if (workers <= NATIVE_PAIR_MAX_FRONTIER / ranks)
        workers *= ranks;
    else
        workers = NATIVE_PAIR_MAX_FRONTIER;
    target = workers <= NATIVE_PAIR_MAX_FRONTIER / 4
        ? 4 * workers : NATIVE_PAIR_MAX_FRONTIER;
#ifdef BALLS4SCANLEV
    target = MAX((INTEGER)64, target);
#else
    target = MAX((INTEGER)16, target);
#endif
    return MIN(NATIVE_PAIR_MAX_FRONTIER, target);
}

int cballs_native_pair_leaf_capacity(
        const struct cmdline_data *cmd, const struct global_data *gd,
        INTEGER point_count)
{
    const int base = MAX(1, cmd->nsmooth);
    const int adaptive_base = MIN(64, base);
    const int minimum = MAX(1, (3 * adaptive_base + 3) / 4);
    const int maximum = MIN(
        64, adaptive_base + MAX(1, adaptive_base / 4));
    real bin_width;
    real target_radius;
    real surface_density;
    real expected_occupancy;
    int capacity;

    if (scanopt(cmd->options, "dual-node-singleton-leaves")) return 1;
    if (cballs_opt_no_two_balls(cmd) || !(cmd->theta > 0.0)
        || point_count <= 0)
        return base;

    /* The unit-sphere pair workload is dominated by the widest, outermost
     * annulus.  Estimate the number of bodies in a cell whose radius is half
     * the accepted theta-scaled bin width.  The 75--125 percent calibration
     * band around nsmooth prevents body-pair traffic from exploding when a
     * dense catalog predicts an overly coarse leaf. */
    if (cmd->useLogHist) {
        if (!(gd->deltaR > 0.0) || !(cmd->rangeN > 0.0)) return base;
        bin_width = cmd->rangeN * (rpow(10.0, gd->deltaR) - 1.0);
    } else {
        bin_width = gd->deltaR;
    }
    target_radius = 0.5 * cmd->theta * bin_width;
    surface_density = (real)point_count / (4.0 * PI);
    expected_occupancy = surface_density * PI * target_radius * target_radius;
    if (!isfinite(expected_occupancy) || !(expected_occupancy > 0.0))
        return base;
    capacity = (int)rfloor(expected_occupancy + 0.5);
    return MAX(minimum, MIN(maximum, capacity));
}

static inline INTEGER native_pair_node_count(const fcfc_ballnode *node)
{
    return node->last - node->first + 1;
}

static inline bool native_pair_node_is_leaf(const fcfc_ballnode *node)
{
    return node->left < 0;
}

static inline real native_pair_center_distance_squared(
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

static inline real native_pair_point_distance_squared(
        struct cmdline_data *cmd, struct global_data *gd,
        const fcfc_ballpoint *point1, const fcfc_ballpoint *point2)
{
    compute_vector dr;
    real distance2;

    DOTPSUBV(distance2, dr, point1->pos, point2->pos);
    if (cmd->usePeriodic) {
        VWrapAll(dr);
        DOTVP(distance2, dr, dr);
    }
    return distance2;
}

static inline int native_pair_bin_index(
        const native_pair_context *context, real distance)
{
    const struct cmdline_data *cmd = context->cmd;
    const struct global_data *gd = context->gd;
    int bin;

    if (!(distance > cmd->rminHist && distance < cmd->rangeN)) return -1;
    if (cmd->useLogHist) {
        if (cmd->rminHist == 0.0) {
            bin = (int)(cmd->logHistBinsPD
                * (rlog10(distance) - rlog10(cmd->rangeN))
                + cmd->sizeHistN) + 1;
        } else {
            bin = (int)(rlog10(distance / cmd->rminHist)
                * gd->i_deltaR) + 1;
        }
    } else {
        bin = (int)((distance - cmd->rminHist) * gd->i_deltaR) + 1;
    }
    return bin >= 1 && bin <= cmd->sizeHistN ? bin : -1;
}

static inline int native_pair_bin_index_squared(
        const native_pair_context *context, real distance2)
{
    const struct cmdline_data *cmd = context->cmd;
    const struct global_data *gd = context->gd;
    int bin;

    if (!(distance2 > context->minimum2
          && distance2 < context->maximum2)) return -1;
    if (!cmd->useLogHist)
        return native_pair_bin_index(context, rsqrt(distance2));
    if (cmd->rminHist == 0.0) {
        bin = (int)(cmd->logHistBinsPD
            * (0.5 * rlog10(distance2) - rlog10(cmd->rangeN))
            + cmd->sizeHistN) + 1;
    } else {
        bin = (int)(0.5 * rlog10(distance2 / context->minimum2)
            * gd->i_deltaR) + 1;
    }
    return bin >= 1 && bin <= cmd->sizeHistN ? bin : -1;
}

static int native_pair_initialize_acceptance(native_pair_context *context)
{
    struct cmdline_data *cmd = context->cmd;
    const struct global_data *gd = context->gd;

    context->minimum2 = rsqr(cmd->rminHist);
    context->maximum2 = rsqr(cmd->rangeN);
    context->theta2 = rsqr(cmd->theta);
    context->bin_size = cmd->useLogHist
        ? rlog(10.0) * gd->deltaR : gd->deltaR;
    context->bin_slop = cmd->theta * context->bin_size;
    context->bin_slop2 = rsqr(context->bin_slop);
    context->half_bin_plus_slop2 =
        0.25 * rsqr(context->bin_size + context->bin_slop);
    context->log_acceptance_limit2 = MIN(
        context->theta2,
        MAX(context->bin_slop2, context->half_bin_plus_slop2));
    if (!(context->bin_size > 0.0)) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: radial bin width must be positive",
                 native_pair_method(cmd));
        return FAILURE;
    }
    context->inverse_bin_size = 1.0 / context->bin_size;
    return SUCCESS;
}

#ifdef NATIVE_PAIR_VECTOR_LOG
static void native_pair_accumulate_body_batch_flush(
        const native_pair_context *, native_pair_histogram *,
        real *, real *, real *, int);
#endif

static inline void native_pair_accumulate_body(
        const native_pair_context *context, native_pair_histogram *hist,
        const fcfc_ballpoint *point1, const fcfc_ballpoint *point2)
{
#ifdef NATIVE_PAIR_VECTOR_LOG
    const real distance2 = native_pair_point_distance_squared(
        context->cmd, context->gd, point1, point2);
    const int index = hist->batch_count;

    if (!(distance2 > context->minimum2
          && distance2 < context->maximum2)) return;
    hist->batch_distance2[index] = distance2;
    hist->batch_denominator[index] = context->weighted_normalization
        ? point1->weight * point2->weight : 1.0;
    hist->batch_numerator[index] = context->weighted_signal
        ? point1->weighted_kappa * point2->weighted_kappa
        : point1->kappa * point2->kappa;
    hist->batch_count++;
    if (hist->batch_count == NATIVE_PAIR_BATCH_SIZE) {
        native_pair_accumulate_body_batch_flush(
            context, hist, hist->batch_distance2, hist->batch_denominator,
            hist->batch_numerator, hist->batch_count);
        hist->batch_count = 0;
    }
#else
    const int bin = native_pair_bin_index_squared(
        context, native_pair_point_distance_squared(
            context->cmd, context->gd, point1, point2));
    real denominator;
    real numerator;

    if (bin < 0) return;
    denominator = context->weighted_normalization
        ? point1->weight * point2->weight : 1.0;
    numerator = context->weighted_signal
        ? point1->weighted_kappa * point2->weighted_kappa
        : point1->kappa * point2->kappa;
    hist->pair_count[bin] += context->pair_scale;
    hist->weight_product[bin] += context->pair_scale * denominator;
    hist->field_product[bin] += context->pair_scale * numerator;
    hist->body_pairs++;
#endif
}

#ifdef NATIVE_PAIR_VECTOR_LOG
static void native_pair_vector_log(real *output, const real *input, int count)
{
#if defined(__APPLE__)
#if defined(DOUBLEPREC)
    vvlog(output, input, &count);
#else
    vvlogf(output, input, &count);
#endif
#elif defined(CBALLS_HAVE_SLEEF)
#pragma omp simd
    for (int i = 0; i < count; i++) {
#if defined(DOUBLEPREC)
        output[i] = Sleef_log_u10(input[i]);
#else
        output[i] = Sleef_logf_u10(input[i]);
#endif
    }
#else
#pragma omp simd
    for (int i = 0; i < count; i++) output[i] = rlog(input[i]);
#endif
}

static void native_pair_accumulate_body_batch_flush(
        const native_pair_context *context, native_pair_histogram *hist,
        real *distance2, real *denominator, real *numerator, int count)
{
    int bins[NATIVE_PAIR_BATCH_SIZE];
    const struct cmdline_data *cmd = context->cmd;

    if (count <= 0) return;
    if (cmd->useLogHist && cmd->rminHist > 0.0) {
        real logarithms[NATIVE_PAIR_BATCH_SIZE];
        const real minimum2 = cmd->rminHist * cmd->rminHist;
        const real scale = 0.5 * context->gd->i_deltaR / rlog(10.0);

        for (int i = 0; i < count; i++) distance2[i] /= minimum2;
        if (count >= 8)
            native_pair_vector_log(logarithms, distance2, count);
        else
        {
            for (int i = 0; i < count; i++)
                logarithms[i] = rlog(distance2[i]);
        }
        for (int i = 0; i < count; i++) {
            const int bin = (int)(logarithms[i] * scale) + 1;

            bins[i] = bin >= 1 && bin <= cmd->sizeHistN ? bin : -1;
        }
    } else {
        for (int i = 0; i < count; i++)
            bins[i] = native_pair_bin_index_squared(context, distance2[i]);
    }

    for (int i = 0; i < count; i++) {
        const int bin = bins[i];

        if (bin < 0) continue;
        hist->pair_count[bin] += context->pair_scale;
        hist->weight_product[bin] += context->pair_scale * denominator[i];
        hist->field_product[bin] += context->pair_scale * numerator[i];
        hist->body_pairs++;
    }
}
#endif

static inline bool native_pair_outside_range(
        const native_pair_context *context,
        const fcfc_ballnode *node1, const fcfc_ballnode *node2,
        real distance2)
{
    const real size = (real)node1->radius + (real)node2->radius;
    const real minimum = context->cmd->rminHist;
    const real maximum = context->cmd->rangeN;

    if (size <= minimum && distance2 <= rsqr(minimum - size)) return TRUE;
    if (distance2 >= rsqr(maximum + size)) return TRUE;
    return FALSE;
}

/* Return a common bin, -1 when nodes must split, and -2 when an accepted
 * approximate pair has its center outside the histogram domain. */
static int native_pair_two_ball_bin(
        const native_pair_context *context,
        const fcfc_ballnode *node1, const fcfc_ballnode *node2,
        real distance2)
{
    const struct cmdline_data *cmd = context->cmd;
    const real size = (real)node1->radius + (real)node2->radius;
    const real size2 = size * size;
    bool have_coordinate = FALSE;
    real coordinate = 0.0;
    int center_bin;

    if (!context->use_two_balls || !(distance2 > 0.0)
        || !(cmd->theta > 0.0))
        return -1;

    if (context->use_bin_slop) {
        real fraction;

        if (cmd->useLogHist) {
            if (!(cmd->rminHist > 0.0)) return -1;
            if (size2 > context->log_acceptance_limit2 * distance2)
                return -1;
            if (size2 > context->bin_slop2 * distance2) {
                const real relative_size2 = size2 / distance2;

                if (size2 > context->half_bin_plus_slop2 * distance2)
                    return -1;
                coordinate = 0.5 * rlog(distance2 / context->minimum2)
                           * context->inverse_bin_size;
                have_coordinate = TRUE;
                fraction = coordinate - rfloor(coordinate);
                if (fraction > 0.5) fraction = 1.0 - fraction;
                if (size2 > rsqr(fraction * context->bin_size
                                 + context->bin_slop) * distance2)
                    return -1;
                fraction = coordinate - rfloor(coordinate);
                if (size2
                    > rsqr(fraction * context->bin_size
                           + context->bin_slop - relative_size2)
                      * distance2)
                    return -1;
            }
        } else {
            if (size2 > context->theta2 * distance2) return -1;
            if (size2 > context->bin_slop2) {
                if (size2 > context->half_bin_plus_slop2) return -1;
                coordinate = (rsqrt(distance2) - cmd->rminHist)
                           / context->bin_size;
                fraction = coordinate - rfloor(coordinate);
                if (fraction > 0.5) fraction = 1.0 - fraction;
                if (size > fraction * context->bin_size
                         + context->bin_slop) return -1;
            }
        }
        if (!(distance2 > context->minimum2
              && distance2 < context->maximum2))
            return -2;
        if (!have_coordinate)
            coordinate = 0.5 * rlog(distance2 / context->minimum2)
                       * context->inverse_bin_size;
        center_bin = (int)coordinate + 1;
        if (center_bin < 1 || center_bin > cmd->sizeHistN) return -2;
        return center_bin < 0 ? -2 : center_bin;
    }

    {
        const real distance = rsqrt(distance2);
        const real lower = distance - size;
        const real upper = distance + size;

        center_bin = native_pair_bin_index(context, distance);
        if (center_bin < 0
            || native_pair_bin_index(context, lower) != center_bin
            || native_pair_bin_index(context, upper) != center_bin)
            return -1;
        return center_bin;
    }
}

static inline void native_pair_accumulate_nodes(
        const native_pair_context *context, native_pair_histogram *hist,
        const fcfc_ballnode *node1, const fcfc_ballnode *node2, int bin)
{
    const real count = (real)native_pair_node_count(node1)
                     * (real)native_pair_node_count(node2);
    const real denominator = context->weighted_normalization
        ? node1->field_weight_sum * node2->field_weight_sum : count;
    const real numerator = context->weighted_signal
        ? node1->weighted_kappa_sum * node2->weighted_kappa_sum
        : node1->kappa_sum * node2->kappa_sum;

    hist->pair_count[bin] += context->pair_scale * count;
    hist->weight_product[bin] += context->pair_scale * denominator;
    hist->field_product[bin] += context->pair_scale * numerator;
    hist->cell_pairs++;
}

static void native_pair_process_pair(
        const native_pair_context *context,
        const fcfc_balltreeptr tree1, INTEGER index1,
        const fcfc_balltreeptr tree2, INTEGER index2,
        native_pair_histogram *hist)
{
    const fcfc_ballnode *node1 = &tree1->nodes[index1];
    const fcfc_ballnode *node2 = &tree2->nodes[index2];
    const bool leaf1 = native_pair_node_is_leaf(node1);
    const bool leaf2 = native_pair_node_is_leaf(node2);
    const real distance2 = native_pair_center_distance_squared(
        context->cmd, context->gd, node1, node2);
    int bin;

    if (native_pair_outside_range(context, node1, node2, distance2)) return;
    bin = native_pair_two_ball_bin(context, node1, node2, distance2);
    if (bin == -2) return;
    if (bin >= 0) {
        native_pair_accumulate_nodes(context, hist, node1, node2, bin);
        return;
    }

    if (leaf1 && leaf2) {
        INTEGER point1;
        INTEGER point2;

        for (point1 = node1->first; point1 <= node1->last; point1++)
            for (point2 = node2->first; point2 <= node2->last; point2++)
                native_pair_accumulate_body(
                    context, hist, &tree1->packed_points[point1],
                    &tree2->packed_points[point2]);
        return;
    }
    if (leaf1) {
        native_pair_process_pair(
            context, tree1, index1, tree2, node2->left, hist);
        native_pair_process_pair(
            context, tree1, index1, tree2, node2->right, hist);
        return;
    }
    if (leaf2) {
        native_pair_process_pair(
            context, tree1, node1->left, tree2, index2, hist);
        native_pair_process_pair(
            context, tree1, node1->right, tree2, index2, hist);
        return;
    }

    {
        bool split1 = FALSE;
        bool split2 = FALSE;
        const real size1 = (real)node1->radius;
        const real size2 = (real)node2->radius;
        real effective_width2;

        if (context->cmd->useLogHist)
            effective_width2 = context->bin_slop2 * distance2;
        else
            effective_width2 = context->bin_slop2;

        if (size2 > size1) {
            split2 = TRUE;
            if (!(size2 > 2.0 * size1))
                split1 = size1 * size1
                    > rsqr(NATIVE_PAIR_SPLIT_FACTOR) * effective_width2;
        } else {
            split1 = TRUE;
            if (!(size1 > 2.0 * size2))
                split2 = size2 * size2
                    > rsqr(NATIVE_PAIR_SPLIT_FACTOR) * effective_width2;
        }

        if (split1 && split2) {
            native_pair_process_pair(context, tree1, node1->left,
                                     tree2, node2->left, hist);
            native_pair_process_pair(context, tree1, node1->left,
                                     tree2, node2->right, hist);
            native_pair_process_pair(context, tree1, node1->right,
                                     tree2, node2->left, hist);
            native_pair_process_pair(context, tree1, node1->right,
                                     tree2, node2->right, hist);
        } else if (split1) {
            native_pair_process_pair(context, tree1, node1->left,
                                     tree2, index2, hist);
            native_pair_process_pair(context, tree1, node1->right,
                                     tree2, index2, hist);
        } else {
            native_pair_process_pair(context, tree1, index1,
                                     tree2, node2->left, hist);
            native_pair_process_pair(context, tree1, index1,
                                     tree2, node2->right, hist);
        }
    }
}

static void native_pair_process_auto(
        const native_pair_context *context,
        const fcfc_balltreeptr tree, INTEGER index,
        native_pair_histogram *hist)
{
    const fcfc_ballnode *node = &tree->nodes[index];

    if (2.0 * (real)node->radius <= context->cmd->rminHist) return;
    if (native_pair_node_is_leaf(node)) {
        INTEGER point1;
        INTEGER point2;

        for (point1 = node->first; point1 <= node->last; point1++)
            for (point2 = point1 + 1; point2 <= node->last; point2++)
                native_pair_accumulate_body(
                    context, hist, &tree->packed_points[point1],
                    &tree->packed_points[point2]);
        return;
    }
    native_pair_process_auto(context, tree, node->left, hist);
    native_pair_process_auto(context, tree, node->right, hist);
    native_pair_process_pair(
        context, tree, node->left, tree, node->right, hist);
}

static int native_pair_allocate_tasks(
        struct cmdline_data *cmd, INTEGER task_count, size_t stride,
        real **histograms, INTEGER **body_pairs, INTEGER **cell_pairs)
{
    size_t tasks;
    size_t values;

    *histograms = NULL;
    *body_pairs = NULL;
    *cell_pairs = NULL;
    if (task_count <= 0 || stride == 0 || stride > SIZE_MAX / 3)
        goto invalid_size;
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
                 "%s: pair-task allocation failed", native_pair_method(cmd));
        return FAILURE;
    }
    return SUCCESS;

invalid_size:
    snprintf(cmd->error_message, _ERRORMSGSIZE_,
             "%s: pair-task size overflow", native_pair_method(cmd));
    return FAILURE;
}

static int native_pair_compute_balls4_cf(
        struct cmdline_data *cmd, struct global_data *gd, INTEGER nbody)
{
    real volume = 1.0;
    const real body_count = (real)nbody;
    int coordinate;

    DO_COORD(coordinate) volume *= gd->Box[coordinate];
    for (int bin = 1; bin <= cmd->sizeHistN; bin++) {
        real r0;
        real r1;

        if (gd->histNN[bin] == 0.0 || body_count == 0.0 || volume == 0.0) {
            gd->histCF[bin] = 0.0;
            continue;
        }
        if (cmd->useLogHist) {
            if (cmd->rminHist == 0.0) {
                r0 = rpow(10.0, ((real)(bin-cmd->sizeHistN))
                    /cmd->logHistBinsPD + rlog10(cmd->rangeN));
                r1 = rpow(10.0, ((real)(bin+1-cmd->sizeHistN))
                    /cmd->logHistBinsPD + rlog10(cmd->rangeN));
            } else {
                r0 = rpow(10.0, rlog10(cmd->rminHist)
                    + (real)bin*gd->deltaR);
                r1 = rpow(10.0, rlog10(cmd->rminHist)
                    + (real)(bin+1)*gd->deltaR);
            }
        } else {
            r0 = cmd->rminHist + (real)(bin-1)*gd->deltaR;
            r1 = cmd->rminHist + (real)bin*gd->deltaR;
        }
#if NDIM == 3
        if (cballs_opt_cute_box(cmd)) {
            const real shell_volume =
                4.0*PI*(r1*r1*r1-r0*r0*r0)/3.0;
            const real expected =
                body_count*body_count*shell_volume/volume;
            gd->histCF[bin] = cballs_normalize_or_zero(
                gd->histNN[bin], expected) - 1.0;
        } else {
            const real norm = volume
                /(2.0*PI*rpow(gd->deltaR,3.0)*body_count*body_count);
            gd->histCF[bin] =
                gd->histNN[bin]*norm/rsqr((int)bin-0.5) - 1.0;
        }
#else
        {
            const real norm = volume
                /(PI*rpow(gd->deltaR,2.0)*body_count*body_count);
            gd->histCF[bin] =
                gd->histNN[bin]*norm/((real)bin-0.5) - 1.0;
        }
#endif
    }
    return SUCCESS;
}

static int native_pair_run_tasks(
        const native_pair_context *context,
        const cballs_native_pair_policy *policy,
        fcfc_balltreeptr tree1, fcfc_balltreeptr tree2,
        bool auto_correlation, const INTEGER *frontier1,
        INTEGER frontier_count1, const INTEGER *frontier2,
        INTEGER frontier_count2)
{
    struct cmdline_data *cmd = context->cmd;
    struct global_data *gd = context->gd;
    const size_t stride = (size_t)cmd->sizeHistN + 1;
    real *thread_histograms = NULL;
    INTEGER *thread_body_counts = NULL;
    INTEGER *thread_cell_counts = NULL;
    int thread_count = 1;
    int allocation_status;
    int reduction_status = SUCCESS;
    int status = FAILURE;
    double phase_started = 0.0;

    if (tree1 == NULL || tree2 == NULL
        || tree1->packed_points == NULL || tree2->packed_points == NULL) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: compact point storage is unavailable",
                 native_pair_method(cmd));
        goto cleanup;
    }
#ifdef OPENMPCODE
    thread_count = MAX(1, omp_get_max_threads());
#endif
    allocation_status = native_pair_allocate_tasks(
        cmd, (INTEGER)thread_count, stride, &thread_histograms,
        &thread_body_counts, &thread_cell_counts);
    if (native_pair_consensus(
            policy, cmd, allocation_status,
            "native-octree pair-task allocation") == FAILURE)
        goto cleanup;

    if (context->profile) phase_started = native_pair_timer_now();
#ifdef OPENMPCODE
#pragma omp parallel
#endif
    {
#ifdef OPENMPCODE
        const int thread = omp_get_thread_num();
#else
        const int thread = 0;
#endif
        real *base = thread_histograms + (size_t)thread * 3 * stride;
        native_pair_histogram hist;

        hist.pair_count = base;
        hist.weight_product = base + stride;
        hist.field_product = base + 2 * stride;
        hist.body_pairs = 0;
        hist.cell_pairs = 0;
#ifdef NATIVE_PAIR_VECTOR_LOG
        hist.batch_count = 0;
#endif

#ifdef OPENMPCODE
#pragma omp for schedule(dynamic,1)
#endif
        for (INTEGER task = 0; task < frontier_count1; task++) {
            if (!native_pair_task_owned(policy, cmd, task)) continue;
            if (auto_correlation) {
                native_pair_process_auto(
                    context, tree1, frontier1[task], &hist);
                for (INTEGER other = task + 1;
                     other < frontier_count2; other++)
                    native_pair_process_pair(
                        context, tree1, frontier1[task],
                        tree2, frontier2[other], &hist);
            } else {
                for (INTEGER other = 0; other < frontier_count2; other++)
                    native_pair_process_pair(
                        context, tree1, frontier1[task],
                        tree2, frontier2[other], &hist);
            }
        }
#ifdef NATIVE_PAIR_VECTOR_LOG
        native_pair_accumulate_body_batch_flush(
            context, &hist, hist.batch_distance2, hist.batch_denominator,
            hist.batch_numerator, hist.batch_count);
#endif
        thread_body_counts[thread] = hist.body_pairs;
        thread_cell_counts[thread] = hist.cell_pairs;
    }
    if (context->profile)
        context->timers->traversal += native_pair_timer_now() - phase_started;

    if (context->profile) phase_started = native_pair_timer_now();

#ifdef OPENMPCODE
#pragma omp parallel
    {
        for (int merge_span = 1; merge_span < thread_count;
             merge_span *= 2) {
            const int blocks =
                (thread_count + 2 * merge_span - 1) / (2 * merge_span);

#pragma omp for schedule(static)
            for (int block = 0; block < blocks; block++) {
                const int destination = block * 2 * merge_span;
                const int source = destination + merge_span;

                if (source >= thread_count) continue;
                for (size_t value = 0; value < 3 * stride; value++)
                    thread_histograms[(size_t)destination * 3 * stride + value]
                        += thread_histograms[(size_t)source * 3 * stride + value];
                thread_body_counts[destination] += thread_body_counts[source];
                thread_cell_counts[destination] += thread_cell_counts[source];
            }
        }
    }
#endif
    if (native_pair_reduce_reals(
            policy, cmd, thread_histograms, 3 * stride) == FAILURE)
        reduction_status = FAILURE;
    if (native_pair_reduce_integers(
            policy, cmd, thread_body_counts, 1) == FAILURE)
        reduction_status = FAILURE;
    if (native_pair_reduce_integers(
            policy, cmd, thread_cell_counts, 1) == FAILURE)
        reduction_status = FAILURE;
    if (native_pair_consensus(
            policy, cmd, reduction_status,
            "native-octree pair-task reduction") == FAILURE)
        goto cleanup;
    if (!native_pair_publish(policy, cmd)) {
        status = SUCCESS;
        goto cleanup;
    }

    for (int bin = 1; bin <= cmd->sizeHistN; bin++) {
        gd->histNN[bin] += thread_histograms[bin];
        gd->histNNSubXi2pcf[bin] += thread_histograms[bin];
        gd->histXi2pcf[bin] += thread_histograms[2 * stride + (size_t)bin];
    }
    gd->nbbcalc += thread_body_counts[0];
    gd->nbccalc += thread_cell_counts[0];

    for (int bin = 1; bin <= cmd->sizeHistN; bin++) {
        gd->histXi2pcf[bin] = cballs_normalize_or_zero(
            gd->histXi2pcf[bin], thread_histograms[stride + (size_t)bin]);
    }
    if (context->profile)
        context->timers->reduction += native_pair_timer_now() - phase_started;

    if (cballs_opt_compute_histn(cmd) && cballs_opt_and_cf(cmd)) {
        if (policy->balls4_density_normalization) {
            if (native_pair_compute_balls4_cf(
                    cmd, gd, tree1->npoint) == FAILURE)
                goto cleanup;
            status = SUCCESS;
            goto cleanup;
        }
        for (int bin = 1; bin <= cmd->sizeHistN; bin++) gd->histNN[bin] *= 2.0;
#ifdef LONGINT
        if (tree1->npoint > INT_MAX) {
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     "%s: and-CF body count exceeds int",
                     native_pair_method(cmd));
            goto cleanup;
        }
#endif
        if (search_compute_HistN(cmd, gd, (int)tree1->npoint) == FAILURE)
            goto cleanup;
    }
    status = SUCCESS;

cleanup:
    status = native_pair_consensus(
        policy, cmd, status, "native-octree pair publication");
    free(thread_cell_counts);
    free(thread_body_counts);
    free(thread_histograms);
    return status;
}
#endif /* TWOPCF */

int cballs_native_octree_pair_search(
        struct cmdline_data *cmd, struct global_data *gd,
        bodyptr *body_tables, INTEGER *body_counts,
        INTEGER ipmin, INTEGER *ipmax, int cat1, int cat2,
        const cballs_native_pair_policy *policy)
{
#ifdef TWOPCF
    const bool auto_correlation = cat1 == cat2;
    int leaf_capacity1;
    int leaf_capacity2;
    INTEGER target;
    native_pair_context context;
    fcfc_balltreeptr tree1 = NULL;
    fcfc_balltreeptr tree2 = NULL;
    INTEGER *frontier1 = NULL;
    INTEGER *frontier2 = NULL;
    INTEGER frontier_count1 = 0;
    INTEGER frontier_count2 = 0;
    bool tree_cache_hit1 = FALSE;
    bool tree_cache_hit2 = FALSE;
    int operation_status;
    int status = FAILURE;
    const double cpustart = CPUTIME;
    double phase_started = 0.0;
    native_pair_phase_timers phase_timers = {0};

    if (policy == NULL) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: missing native-pair policy", native_pair_method(cmd));
        return FAILURE;
    }
    target = native_pair_frontier_target(cmd, policy);
    if (scanopt(cmd->options, "only-3pcf")) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: compact pair search cannot run only-3pcf",
                 native_pair_method(cmd));
        return FAILURE;
    }
    if (ipmin != 1 || ipmax[cat1] != body_counts[cat1]) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s requires the complete pivot catalog",
                 native_pair_method(cmd));
        return FAILURE;
    }
    if (cmd->nsmooth <= 0) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s requires nsmooth > 0", native_pair_method(cmd));
        return FAILURE;
    }
    context.cmd = cmd;
    context.gd = gd;
    context.profile = scanopt(cmd->options, "dual-node-profile");
    context.timers = &phase_timers;
    leaf_capacity1 = cballs_native_pair_leaf_capacity(
        cmd, gd, body_counts[cat1]);
    leaf_capacity2 = auto_correlation ? leaf_capacity1
        : cballs_native_pair_leaf_capacity(cmd, gd, body_counts[cat2]);

    verb_print(cmd->verbose,
               "Search: Running %s with shared compact native-octree 2PCF\n",
               native_pair_method(cmd));
#ifdef BALLS4SCANLEV
    verb_print(cmd->verbose,
               "%s: BALLS4SCANLEV balanced native frontier enabled\n",
               native_pair_method(cmd));
#endif
    if (policy->parallel != NULL && policy->parallel->rank_count > 1)
        verb_print(cmd->verbose,
                   "%s: %d ranks with deterministic cyclic frontier ownership\n",
                   native_pair_method(cmd), policy->parallel->rank_count);
    if (cballs_opt_no_two_balls(cmd)
        || (policy->no_one_ball_is_exact && cballs_opt_no_one_ball(cmd)))
        verb_print(cmd->verbose,
                   "exact body-pair accumulation enabled\n");
    else if (scanopt(cmd->options, "dual-node-bin-slop"))
        verb_print(cmd->verbose,
                   "dual-node-compatible controlled 2PCF bin slop enabled\n");
    else
        verb_print(cmd->verbose,
                   "conservative same-bin dual-node acceptance enabled\n");
    verb_print(cmd->verbose >= 2 || context.profile,
               "%s: adaptive 2PCF leaf capacities = %d, %d "
               "(configured nsmooth=%d)\n",
               native_pair_method(cmd), leaf_capacity1, leaf_capacity2,
               cmd->nsmooth);

#ifdef OPENMPCODE
    ThreadCount(cmd, gd, body_counts[cat1], cat1);
#endif
    operation_status = search_init_gd_hist_sincos(cmd, gd);
    if (native_pair_consensus(
            policy, cmd, operation_status,
            "native-octree histogram initialization") == FAILURE)
        goto cleanup;
    operation_status = native_pair_initialize_acceptance(&context);
    if (native_pair_consensus(
            policy, cmd, operation_status,
            "native-octree pair-acceptance initialization") == FAILURE)
        goto cleanup;
    if (context.profile) phase_started = native_pair_timer_now();
    if (scanopt(cmd->options, "no-native-tree-cache"))
        operation_status = octree_2balls_tree_build(
            cmd, gd, body_tables[cat1], body_counts[cat1], leaf_capacity1,
            &tree1);
    else
        operation_status = octree_2balls_tree_build_cached(
            cmd, gd, body_tables[cat1], body_counts[cat1], leaf_capacity1,
            &tree1, &tree_cache_hit1);
    if (native_pair_consensus(
            policy, cmd, operation_status,
            "native-octree first compact-tree construction") == FAILURE)
        goto cleanup;
    if (auto_correlation) {
        tree2 = tree1;
    } else {
        if (scanopt(cmd->options, "no-native-tree-cache"))
            operation_status = octree_2balls_tree_build(
                cmd, gd, body_tables[cat2], body_counts[cat2],
                leaf_capacity2, &tree2);
        else
            operation_status = octree_2balls_tree_build_cached(
                cmd, gd, body_tables[cat2], body_counts[cat2],
                leaf_capacity2, &tree2, &tree_cache_hit2);
        if (native_pair_consensus(
                policy, cmd, operation_status,
                "native-octree second compact-tree construction") == FAILURE)
            goto cleanup;
    }
    if (context.profile)
        phase_timers.build += native_pair_timer_now() - phase_started;
    verb_print(context.profile,
               "%s: compact-tree cache = %s, %s\n",
               native_pair_method(cmd),
               scanopt(cmd->options, "no-native-tree-cache") ? "disabled" :
                   (tree_cache_hit1 ? "hit" : "miss"),
               auto_correlation ? "shared" :
                   (scanopt(cmd->options, "no-native-tree-cache") ? "disabled" :
                    (tree_cache_hit2 ? "hit" : "miss")));

    if (context.profile) phase_started = native_pair_timer_now();
    operation_status = octree_2balls_tree_frontier(
        cmd, tree1, target, &frontier1, &frontier_count1);
    if (native_pair_consensus(
            policy, cmd, operation_status,
            "native-octree pivot-frontier construction") == FAILURE)
        goto cleanup;
    if (auto_correlation) {
        frontier2 = frontier1;
        frontier_count2 = frontier_count1;
    } else {
        operation_status = octree_2balls_tree_frontier(
            cmd, tree2, target, &frontier2, &frontier_count2);
        if (native_pair_consensus(
                policy, cmd, operation_status,
                "native-octree neighbor-frontier construction") == FAILURE)
            goto cleanup;
    }
    if (context.profile)
        phase_timers.frontier += native_pair_timer_now() - phase_started;

    context.use_two_balls = !cballs_opt_no_two_balls(cmd)
        && !(policy->no_one_ball_is_exact && cballs_opt_no_one_ball(cmd));
    context.use_bin_slop = scanopt(cmd->options, "dual-node-bin-slop");
    context.weighted_normalization = cballs_opt_weights_norm(cmd);
    context.weighted_signal = policy->weighted_signal
        || context.weighted_normalization;
    context.pair_scale = auto_correlation && policy->honor_asymmetric
        && cballs_opt_asymmetric(cmd) ? 2.0 : 1.0;

    if (native_pair_run_tasks(
            &context, policy, tree1, tree2, auto_correlation,
            frontier1, frontier_count1,
            frontier2, frontier_count2) == FAILURE)
        goto cleanup;

    gd->cpusearch = CPUTIME - cpustart;
    if (native_pair_publish(policy, cmd))
        verb_print(cmd->verbose,
                   "%s: nbbcalc = %" INTEGER_FMT
                   ", nbccalc = %" INTEGER_FMT
                   ", frontier tasks = %" INTEGER_FMT "\n",
                   native_pair_method(cmd), gd->nbbcalc, gd->nbccalc,
                   frontier_count1);
    if (context.profile && native_pair_publish(policy, cmd))
        verb_print(TRUE,
                   "%s phase-timers: build_wall=%.9g, frontier_wall=%.9g, "
                   "pair_traversal_wall=%.9g, reduction_wall=%.9g\n",
                   native_pair_method(cmd), phase_timers.build,
                   phase_timers.frontier, phase_timers.traversal,
                   phase_timers.reduction);
    verb_print(cmd->verbose, "Going out: CPU time = %lf\n", gd->cpusearch);
    status = SUCCESS;

cleanup:
    if (!auto_correlation) free(frontier2);
    free(frontier1);
    if (!auto_correlation) octree_2balls_tree_release(tree2);
    octree_2balls_tree_release(tree1);
    return native_pair_consensus(
        policy, cmd, status, "native-octree pair cleanup");
#else
    (void)gd;
    (void)body_tables;
    (void)body_counts;
    (void)ipmin;
    (void)ipmax;
    (void)cat1;
    (void)cat2;
    (void)policy;
    snprintf(cmd->error_message, _ERRORMSGSIZE_,
             "%s requires TWOPCFON=1",
             cmd->searchMethod != NULL ? cmd->searchMethod
                                       : "native-octree-pair");
    return FAILURE;
#endif
}
