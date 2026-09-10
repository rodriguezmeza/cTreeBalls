/* Deterministic masked LogMultipole/window estimator for the KD tree. */
#include <float.h>
#include <stdint.h>
#include "globaldefs.h"
#include "kdtree.h"
#include "kdtree_parallel.h"

#if defined(TPCF) && NDIM == 3
#define DUAL_NODE_METHOD_NAME (cmd->searchMethod)
#define DUAL_NODE_ZETA_COMPONENTS 4
static int dual_node_window_orders(const struct cmdline_data *cmd)
{
    return 2 * cmd->mChebyshev + 1;
}
#include "../balltree_2balls_omp/dual_node_edge_correction.h"

#define KD_EDGE_BLOCK 32
#define KD_EDGE_BATCH 8

typedef struct {
    struct cmdline_data *cmd;
    struct global_data *gd;
    ballxptr tree;
    size_t stride, plane, values, triple_values, work_values;
    int orders, window_orders;
    bool weighted, aggregate, run_2pcf;
} kd_edge_context;

static int kd_edge_bin(const kd_edge_context *ctx, real distance)
{
    if (!(distance > ctx->cmd->rminHist && distance < ctx->cmd->rangeN))
        return 0;
    const real value = ctx->cmd->rminHist == 0
        ? ctx->cmd->logHistBinsPD * log10(distance / ctx->cmd->rangeN)
            + ctx->cmd->sizeHistN
        : (ctx->cmd->useLogHist
            ? log10(distance / ctx->cmd->rminHist) * ctx->gd->i_deltaR
            : (distance - ctx->cmd->rminHist) * ctx->gd->i_deltaR);
    if (!isfinite(value) || value < 0 || value >= ctx->cmd->sizeHistN)
        return 0;
    return (int)floor(value) + 1;
}

static bool kd_edge_accept_node(const kd_edge_context *ctx, bodyptr pivot,
                                const ballnode *node, real distance,
                                const compute_vector dr, int bin)
{
    const real radius = node->bnd.geometric_radius;
    if (!ctx->aggregate || bin == 0 || !(distance > radius)
        || kdtree_node_contains_body(ctx->tree, node, pivot)
        || kd_edge_bin(ctx, distance - radius) != bin
        || kd_edge_bin(ctx, distance + radius) != bin)
        return FALSE;
    const real extent = cballs_angular_extent(Pos(pivot), dr, 0.0, radius);
    return extent < 1.0
        && (ctx->window_orders - 1) * asin(extent)
            < 0.05 * MAX(ctx->cmd->theta, DBL_EPSILON);
}

static void kd_edge_add(const kd_edge_context *ctx, bodyptr pivot,
                        real count, real pair_weight, real pair_field,
                        real weight, real field, real weight2, real field2,
                        int bin, const compute_vector dr, bool cell,
                        real *output, real *work, INTEGER *encounters)
{
    const size_t stride = ctx->stride;
    const size_t moment = (size_t)ctx->window_orders * stride;
    real *sr = work;
    real *si = sr + moment;
    real *wr = si + moment;
    real *wi = wr + moment;
    real *scc = wi + moment;
    real *sss = scc + (size_t)ctx->orders * stride;
    real *ssc = sss + (size_t)ctx->orders * stride;
    real *wself = ssc + (size_t)ctx->orders * stride;
    real *pairs = output + ctx->triple_values;

    if (ctx->run_2pcf) {
        pairs[bin] += count;
        pairs[stride + bin] += Weight(pivot) * pair_weight;
#ifdef SMOOTHPIVOT
        const real pivot_pair_field = cballs_opt_smooth_pivot(ctx->cmd)
            ? KappaRmin(pivot) : Weight(pivot)*Kappa(pivot);
#else
        const real pivot_pair_field = Weight(pivot)*Kappa(pivot);
#endif
        pairs[2*stride + bin] += pivot_pair_field * pair_field;
    }
    encounters[cell ? 1 : 0]++;

    real cosine, sine;
    if (!cballs_angular_phase(Pos(pivot), dr, &cosine, &sine)) return;
    real cm = 1.0, sm = 0.0;
    for (int m = 0; m < ctx->window_orders; m++) {
        const size_t k = (size_t)m * stride + bin;
        sr[k] += field*cm;
        si[k] += field*sm;
        wr[k] += weight*cm;
        wi[k] += weight*sm;
        if (m < ctx->orders) {
            scc[k] += field2*cm*cm;
            sss[k] += field2*sm*sm;
            ssc[k] += field2*sm*cm;
        }
        const real next = cm*cosine - sm*sine;
        sm = sm*cosine + cm*sine;
        cm = next;
    }
    wself[bin] += weight2;
}

static int kd_edge_pivot(const kd_edge_context *ctx, bodyptr pivot,
                         real *output, real *work, INTEGER *encounters)
{
    struct global_data *gd = ctx->gd;
    const size_t stride = ctx->stride;
    const size_t moment = (size_t)ctx->window_orders * stride;
    real *sr = work, *si = sr + moment, *wr = si + moment, *wi = wr + moment;
    real *scc = wi + moment;
    real *sss = scc + (size_t)ctx->orders * stride;
    real *ssc = sss + (size_t)ctx->orders * stride;
    real *wself = ssc + (size_t)ctx->orders * stride;
    memset(work, 0, ctx->work_values*sizeof(*work));

    INTEGER cp = KDROOT;
    do {
        const ballnode *node = &ctx->tree->ntab[cp];
        if (node->valid_count == 0) {
            SetNext(cp);
            continue;
        }
        compute_vector dr;
        real distance2;
        DOTPSUBV(distance2, dr, Pos(pivot), node->cmpos);
        if (ctx->cmd->usePeriodic) {
            VWrapAll(dr);
            DOTVP(distance2, dr, dr);
        }
        const real distance = sqrt(distance2);
        const real radius = node->bnd.geometric_radius;
        if (distance - radius >= ctx->cmd->rangeN
            || distance + radius <= ctx->cmd->rminHist) {
            SetNext(cp);
            continue;
        }
        const int bin = kd_edge_bin(ctx, distance);
        if (cp < ctx->tree->nsplit) {
            if (!kd_edge_accept_node(ctx, pivot, node, distance, dr, bin)) {
                cp = Lower(cp);
                continue;
            }
            const real weight = ctx->weighted ? node->weight_sum
                                              : (real)node->valid_count;
            const real field = ctx->weighted ? node->weighted_kappa_sum
                                             : node->kappa_sum;
            const real weight2 = ctx->weighted ? node->weight_sq_sum
                                               : (real)node->valid_count;
            const real field2 = ctx->weighted ? node->weighted_kappa_sq_sum
                                              : node->kappa_sq_sum;
            kd_edge_add(ctx, pivot, (real)node->valid_count, node->weight_sum,
                        node->weighted_kappa_sum, weight, field, weight2,
                        field2, bin, dr, TRUE, output, work, encounters);
            SetNext(cp);
            continue;
        }

        for (INTEGER i = node->first; i <= node->last; i++) {
            bodyptr q = ctx->tree->bptr[i];
            if (q == pivot || (cballs_opt_read_mask(ctx->cmd)
                && Mask(q) == MASK_NODE_MASKED)) continue;
            DOTPSUBV(distance2, dr, Pos(pivot), Pos(q));
            if (ctx->cmd->usePeriodic) {
                VWrapAll(dr);
                DOTVP(distance2, dr, dr);
            }
            const real d = sqrt(distance2);
            const int b = kd_edge_bin(ctx, d);
            if (!b) continue;
            const real w = ctx->weighted ? Weight(q) : 1.0;
            const real f = w*Kappa(q);
            kd_edge_add(ctx, pivot, 1.0, Weight(q), Weight(q)*Kappa(q),
                        w, f, w*w, f*f, b, dr, FALSE,
                        output, work, encounters);
        }
        SetNext(cp);
    } while (cp != KDROOT);

#ifdef SMOOTHPIVOT
    const real pivot_field = cballs_opt_smooth_pivot(ctx->cmd)
        ? KappaRmin(pivot)/MAX((real)NbRmin(pivot), 1.0)
        : (ctx->weighted ? Weight(pivot)*Kappa(pivot) : Kappa(pivot));
#else
    const real pivot_field = ctx->weighted
        ? Weight(pivot)*Kappa(pivot) : Kappa(pivot);
#endif
    const real pivot_weight = ctx->weighted ? Weight(pivot) : 1.0;
    const size_t component = (size_t)ctx->orders * ctx->plane;
    real *window_re = output + 4*component + ctx->plane;
    real *window_im = window_re + (size_t)ctx->window_orders*ctx->plane;
    for (int i = 1; i <= ctx->cmd->sizeHistN; i++)
        for (int j = 1; j <= ctx->cmd->sizeHistN; j++) {
            const size_t ij = (size_t)i*stride + j;
            for (int m = 0; m < ctx->window_orders; m++) {
                const size_t a = (size_t)m*stride+i;
                const size_t b = (size_t)m*stride+j;
                const size_t k = (size_t)m*ctx->plane+ij;
                const real diagonal = i == j ? wself[i] : 0.0;
                window_re[k] += pivot_weight*(wr[a]*wr[b]+wi[a]*wi[b]-diagonal);
                window_im[k] += pivot_weight*(wi[a]*wr[b]-wr[a]*wi[b]);
                if (m < ctx->orders) {
                    output[k] += pivot_field*(sr[a]*sr[b]-(i == j ? scc[a] : 0));
                    output[component+k] += pivot_field*(si[a]*si[b]-(i == j ? sss[a] : 0));
                    output[2*component+k] += pivot_field*(si[a]*sr[b]-(i == j ? ssc[a] : 0));
                    output[3*component+k] += pivot_field*(sr[a]*si[b]-(i == j ? ssc[a] : 0));
                }
            }
        }
    return SUCCESS;
}

static int kd_edge_output(struct cmdline_data *cmd, struct global_data *gd)
{
    const char *names[4] = {"cos", "sin", "sincos", "cossin"};
    real ***matrices[4] = {gd->histZetaMcos, gd->histZetaMsin,
                           gd->histZetaMsincos, gd->histZetaMcossin};
    for (int c = 0; c < 4; c++)
        for (int m = 0; m <= cmd->mChebyshev; m++)
            if (dual_node_write_edge_matrix(cmd, gd, names[c], m, NULL,
                    matrices[c][m+1], (size_t)cmd->sizeHistN+1) == FAILURE)
                return FAILURE;
    return SUCCESS;
}

global int kdtree_edge_search(struct cmdline_data *cmd, struct global_data *gd,
                              bodyptr *btable, INTEGER *nbody, INTEGER ipmin,
                              INTEGER *ipmax, int cat1, int cat2, ballxptr tree)
{
    kd_edge_context ctx = {0};
    bodyptr *pivots = NULL;
    real *total = NULL, *batch = NULL, *work = NULL;
    size_t pivot_count = 0;
    INTEGER counters[2] = {0,0};
    INTEGER batch_counts[KD_EDGE_BATCH][2] = {{0}};
    int status = FAILURE, failed = FALSE;
    ctx.cmd = cmd; ctx.gd = gd; ctx.tree = tree;
    ctx.orders = cmd->mChebyshev + 1;
    ctx.window_orders = dual_node_window_orders(cmd);
    ctx.stride = (size_t)cmd->sizeHistN + 1;
    ctx.plane = ctx.stride*ctx.stride;
    ctx.weighted = cballs_opt_weights_norm(cmd);
    ctx.aggregate = !cballs_opt_no_one_ball(cmd);
#ifdef TWOPCF
    ctx.run_2pcf = !cballs_opt_only_3pcf(cmd);
#endif
    if (ipmin < 1 || ipmax[cat1] < ipmin || ipmax[cat1] > nbody[cat1]
        || !dual_node_triple_values(ctx.stride, ctx.orders,
                                    ctx.window_orders, &ctx.triple_values))
        goto setup_fail;
    if (ctx.run_2pcf && ctx.triple_values > SIZE_MAX - 3*ctx.stride)
        goto setup_fail;
    ctx.values = ctx.triple_values + (ctx.run_2pcf ? 3*ctx.stride : 0);
    if ((size_t)ctx.window_orders > (SIZE_MAX/ctx.stride-1)/7)
        goto setup_fail;
    ctx.work_values = (4*(size_t)ctx.window_orders
                      +3*(size_t)ctx.orders+1)*ctx.stride;
    if ((size_t)nbody[cat1] > SIZE_MAX/sizeof(*pivots)
        || ctx.values > SIZE_MAX/KD_EDGE_BATCH/sizeof(*batch)
        || ctx.work_values > SIZE_MAX/KD_EDGE_BATCH/sizeof(*work))
        goto setup_fail;
    pivots = malloc((size_t)nbody[cat1]*sizeof(*pivots));
    total = calloc(ctx.values, sizeof(*total));
    batch = calloc(KD_EDGE_BATCH*ctx.values, sizeof(*batch));
    work = calloc(KD_EDGE_BATCH*ctx.work_values, sizeof(*work));
    if (!pivots || !total || !batch || !work) goto setup_fail;
    for (bodyptr p = btable[cat1]+ipmin-1; p < btable[cat1]+ipmax[cat1]; p++) {
        if (cballs_opt_read_mask(cmd) && Mask(p) == MASK_NODE_MASKED) continue;
#ifdef SMOOTHPIVOT
        if (cballs_opt_smooth_pivot(cmd) && Update(p) == FALSE) continue;
#endif
        pivots[pivot_count++] = p;
    }
    if (pivot_count == 0) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: mask/smoothing selected no pivots", cmd->searchMethod);
        goto setup_done;
    }
    status = SUCCESS;
    goto setup_done;
setup_fail:
    snprintf(cmd->error_message, _ERRORMSGSIZE_,
             "%s: invalid bounds, size overflow, or workspace allocation failure",
             cmd->searchMethod);
setup_done:
    status = kdtree_consensus(cmd, status, "KDTREE edge setup");
    if (status == FAILURE) goto cleanup;

    search_init_gd_hist_sincos(cmd, gd);
    const size_t tasks = (pivot_count+KD_EDGE_BLOCK-1)/KD_EDGE_BLOCK;
    for (size_t first = 0; first < tasks; first += KD_EDGE_BATCH) {
        const int slots = (int)MIN((size_t)KD_EDGE_BATCH, tasks-first);
        memset(batch, 0, KD_EDGE_BATCH*ctx.values*sizeof(*batch));
        memset(batch_counts, 0, sizeof(batch_counts));
        failed = FALSE;
#pragma omp parallel for schedule(dynamic) reduction(|:failed)
        for (int slot = 0; slot < slots; slot++) {
            const size_t task = first+(size_t)slot;
            if (!kdtree_task_owned(cmd, (INTEGER)task)) continue;
            const size_t last = MIN(pivot_count, (task+1)*KD_EDGE_BLOCK);
            for (size_t p = task*KD_EDGE_BLOCK; p < last; p++)
                if (kd_edge_pivot(&ctx, pivots[p],
                        batch+(size_t)slot*ctx.values,
                        work+(size_t)slot*ctx.work_values,
                        batch_counts[slot]) == FAILURE) failed = TRUE;
        }
        status = kdtree_consensus(cmd, failed ? FAILURE : SUCCESS,
                                   "KDTREE edge pivots");
        if (status == FAILURE) goto cleanup;
        if (kdtree_reduce(cmd, batch, (size_t)slots*ctx.values) == FAILURE
            || kdtree_reduce_counts(cmd, &batch_counts[0][0],
                                    (size_t)slots*2) == FAILURE) {
            status = FAILURE;
            goto cleanup;
        }
        if (kdtree_publish(cmd))
            for (int slot = 0; slot < slots; slot++) {
                for (size_t i = 0; i < ctx.values; i++)
                    total[i] += batch[(size_t)slot*ctx.values+i];
                counters[0] += batch_counts[slot][0];
                counters[1] += batch_counts[slot][1];
            }
    }
    if (kdtree_publish(cmd)) {
        real ***matrices[4] = {gd->histZetaMcos, gd->histZetaMsin,
                               gd->histZetaMsincos, gd->histZetaMcossin};
        for (int c = 0; c < 4; c++)
            for (int m = 0; m < ctx.orders; m++)
                for (int i = 1; i <= cmd->sizeHistN; i++)
                    for (int j = 1; j <= cmd->sizeHistN; j++)
                        matrices[c][m+1][i][j] = total[
                            ((size_t)c*ctx.orders+m)*ctx.plane
                            +(size_t)i*ctx.stride+j];
#ifdef TWOPCF
        if (ctx.run_2pcf) {
            const real *pairs = total+ctx.triple_values;
            const real symmetry = cballs_opt_asymmetric(cmd) ? 1.0 : 0.5;
            for (int n = 1; n <= cmd->sizeHistN; n++) {
                gd->histNN[n] = symmetry*pairs[n];
                gd->histNNSubXi2pcf[n] = symmetry*pairs[n];
                gd->histXi2pcf[n] = cballs_normalize_or_zero(
                    pairs[2*ctx.stride+n],
                    ctx.weighted ? pairs[ctx.stride+n] : pairs[n]);
            }
        }
#endif
        gd->nbbcalc = counters[0];
        gd->nbccalc = counters[1];
        gd->ncccalc = 0;
        if (!cballs_opt_no_out_hist(cmd)) status = kd_edge_output(cmd, gd);
        if (status == SUCCESS)
            status = dual_node_publish_edge(cmd, gd, total, 1, ctx.values,
                                            ctx.stride, ctx.orders);
    }
    status = kdtree_consensus(cmd, status, "KDTREE edge publication");
cleanup:
    free(work); free(batch); free(total); free(pivots);
    return status;
}
#else
global int kdtree_edge_search(struct cmdline_data *cmd, struct global_data *gd,
                              bodyptr *btable, INTEGER *nbody, INTEGER ipmin,
                              INTEGER *ipmax, int cat1, int cat2, ballxptr tree)
{
    (void)gd; (void)btable; (void)nbody; (void)ipmin; (void)ipmax;
    (void)cat1; (void)cat2; (void)tree;
    cBALLS_FAIL(cmd, "kdtree edge corrections require NDIM=3 and TPCFON=1\n");
}
#endif
