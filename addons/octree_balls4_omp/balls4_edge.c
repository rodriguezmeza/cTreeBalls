/* Masked, distinct-neighbor angular multipoles over the native B4 partition.
 * Body pivots avoid mixing reference axes when correcting the survey window.
 * Fixed pivot blocks are reduced in order, independent of threads/MPI ranks.
 */
#include <float.h>
#include <stdint.h>
#include "globaldefs.h"
#include "balls4_parallel.h"

#if defined(THREEPCFCONVERGENCE) && NDIM == 3
#define TREECORR_METHOD_NAME (cmd->searchMethod)
#define TREECORR_ZETA_COMPONENTS 4
static int treecorr_window_orders(const struct cmdline_data *cmd)
{
    return 2 * cmd->mChebyshev + 1;
}
#include "../balltree_2balls_omp/treecorr_edge_correction.h"

#define B4_EDGE_BLOCK 32
#define B4_EDGE_BATCH 8

typedef struct {
    nodeptr source;
    size_t end;
    real radius;
    real count, weight, field, weight2, field2;
    real pair_weight, pair_field;
} b4_edge_node;

typedef struct {
    struct cmdline_data *cmd;
    struct global_data *gd;
    b4_edge_node *nodes;
    size_t count, capacity, stride, plane, values, triple_values, work_values;
    int orders, window_orders;
    bool weighted, aggregate;
} b4_edge_context;

static size_t b4_edge_copy(b4_edge_context *ctx, nodeptr source)
{
    const size_t index = ctx->count;
    if (index >= ctx->capacity) return SIZE_MAX;
    b4_edge_node *node = &ctx->nodes[ctx->count++];
    node->source = source;
    if (Type(source) == CELL) {
        for (nodeptr child = More(source); child != Next(source); child = Next(child)) {
            const size_t ci = b4_edge_copy(ctx, child);
            if (ci == SIZE_MAX) return SIZE_MAX;
            const b4_edge_node *q = &ctx->nodes[ci];
            node->count += q->count;
            node->weight += q->weight;
            node->field += q->field;
            node->weight2 += q->weight2;
            node->field2 += q->field2;
            node->pair_weight += q->pair_weight;
            node->pair_field += q->pair_field;
            /* Native Radius is theta-scaled, not a geometric enclosing radius. */
            if (q->count > 0) {
                compute_vector dr;
                real distance2;
                DOTPSUBV(distance2, dr, Pos(source), Pos(child));
                const real epsilon = sizeof(real) == sizeof(float) ? FLT_EPSILON : DBL_EPSILON;
                const real enclosing = (sqrt(distance2) + q->radius) * (1.0 + 8.0*epsilon);
                node->radius = MAX(node->radius, enclosing);
            }
        }
    } else if (!cballs_opt_read_mask(ctx->cmd) || Mask(source) != MASK_NODE_MASKED) {
        node->count = 1;
        node->pair_weight = Weight(source);
        node->pair_field = Weight(source) * Kappa(source);
        node->weight = ctx->weighted ? Weight(source) : 1.0;
        node->field = node->weight * Kappa(source);
        node->weight2 = node->weight * node->weight;
        node->field2 = node->field * node->field;
    }
    node->end = ctx->count;
    return index;
}

static int b4_edge_pivots(struct cmdline_data *cmd, nodeptr node,
                          bodyptr first, bodyptr last, bodyptr *pivots,
                          size_t capacity, size_t *count)
{
    if (cballs_opt_read_mask(cmd) && Mask(node) == MASK_NODE_MASKED) return SUCCESS;
    if (Type(node) == CELL) {
        for (nodeptr q = More(node); q != Next(node); q = Next(q))
            if (b4_edge_pivots(cmd, q, first, last, pivots, capacity, count) == FAILURE)
                return FAILURE;
    } else if ((bodyptr)node >= first && (bodyptr)node < last) {
        if (*count >= capacity) return FAILURE;
        pivots[(*count)++] = (bodyptr)node;
    }
    return SUCCESS;
}

static int b4_edge_bin(const b4_edge_context *ctx, real distance)
{
    if (!(distance > ctx->cmd->rminHist && distance < ctx->cmd->rangeN)) return 0;
    const real bin = ctx->cmd->rminHist == 0
        ? ctx->cmd->logHistBinsPD * log10(distance / ctx->cmd->rangeN)
            + ctx->cmd->sizeHistN
        : log10(distance / ctx->cmd->rminHist) * ctx->gd->i_deltaR;
    if (!isfinite(bin) || bin < 0 || bin >= ctx->cmd->sizeHistN) return 0;
    return (int)floor(bin) + 1;
}

static int b4_edge_pivot(const b4_edge_context *ctx, bodyptr pivot,
                         real *output, real *work, INTEGER *encounters)
{
    const int orders = ctx->orders, worders = ctx->window_orders;
    const size_t stride = ctx->stride, moment = (size_t)worders * stride;
    real *sr = work, *si = sr + moment, *wr = si + moment, *wi = wr + moment;
    real *scc = wi + moment, *sss = scc + (size_t)orders * stride;
    real *ssc = sss + (size_t)orders * stride;
    real *wself = ssc + (size_t)orders * stride;
    real *pairs = output + ctx->triple_values;
    memset(work, 0, ctx->work_values * sizeof(*work));
    const real pivot_weight = ctx->weighted ? Weight(pivot) : 1.0;
    const real pivot_field = pivot_weight * Kappa(pivot);

    for (size_t qi = 0; qi < ctx->count;) {
        const b4_edge_node *q = &ctx->nodes[qi];
        nodeptr source = q->source;
        const bool cell = Type(source) == CELL;
        compute_vector dr;
        real distance2;
        if (q->count == 0 || source == (nodeptr)pivot) {
            qi = q->end;
            continue;
        }
        DOTPSUBV(distance2, dr, Pos(pivot), Pos(source));
        const real distance = sqrt(distance2);
        const real radius = q->radius;
        if (distance - radius >= ctx->cmd->rangeN
            || distance + radius <= ctx->cmd->rminHist) {
            qi = q->end;
            continue;
        }
        const int bin = b4_edge_bin(ctx, distance);
        if (cell) {
            const bool pure = !cballs_opt_read_mask(ctx->cmd)
                || mask_node_can_approximate(Mask(source));
            const real angular_extent = cballs_angular_extent(Pos(pivot), dr, 0.0, radius);
            const bool accept = ctx->aggregate && pure && bin != 0
                && radius < distance
                && b4_edge_bin(ctx, distance - radius) == bin
                && b4_edge_bin(ctx, distance + radius) == bin
                && angular_extent < 1.0
                && (worders - 1) * asin(angular_extent) < 0.05 * ctx->cmd->theta;
            if (!accept) { qi++; continue; }
        }
        qi = q->end;
        if (!bin || !(distance > 0)) continue;
        pairs[bin] += q->count;
        pairs[stride + bin] += Weight(pivot) * q->pair_weight;
        pairs[2*stride + bin] += Weight(pivot)*Kappa(pivot) * q->pair_field;
        encounters[cell ? 1 : 0]++;
        real cosine, sine;
        if (!cballs_angular_phase(Pos(pivot), dr, &cosine, &sine)) continue;
        real cm = 1.0, sm = 0.0;
        for (int m = 0; m < worders; m++) {
            const size_t k = (size_t)m * stride + bin;
            sr[k] += q->field * cm; si[k] += q->field * sm;
            wr[k] += q->weight * cm; wi[k] += q->weight * sm;
            if (m < orders) {
                scc[k] += q->field2 * cm*cm;
                sss[k] += q->field2 * sm*sm;
                ssc[k] += q->field2 * sm*cm;
            }
            const real next = cm*cosine - sm*sine;
            sm = sm*cosine + cm*sine; cm = next;
        }
        wself[bin] += q->weight2;
    }

    const size_t plane = ctx->plane, component = (size_t)orders * plane;
    real *window_re = output + 4*component + plane;
    real *window_im = window_re + (size_t)worders*plane;
    for (int i = 1; i <= ctx->cmd->sizeHistN; i++) {
        for (int j = 1; j <= ctx->cmd->sizeHistN; j++) {
            const size_t bin = (size_t)i * stride + j;
            for (int m = 0; m < worders; m++) {
                const size_t a = (size_t)m * stride + i, b = (size_t)m * stride + j;
                const size_t k = (size_t)m * plane + bin;
                const real diagonal = i == j ? wself[i] : 0.0;
                window_re[k] += pivot_weight * (wr[a]*wr[b] + wi[a]*wi[b] - diagonal);
                window_im[k] += pivot_weight * (wi[a]*wr[b] - wr[a]*wi[b]);
                if (m < orders) {
                    output[k] += pivot_field * (sr[a]*sr[b] - (i == j ? scc[a] : 0));
                    output[component+k] += pivot_field * (si[a]*si[b] - (i == j ? sss[a] : 0));
                    output[2*component+k] += pivot_field * (si[a]*sr[b] - (i == j ? ssc[a] : 0));
                    output[3*component+k] += pivot_field * (sr[a]*si[b] - (i == j ? ssc[a] : 0));
                }
            }
        }
    }
    return SUCCESS;
}

static int b4_edge_output(struct cmdline_data *cmd, struct global_data *gd)
{
    string routineName = "BALLS4 edge output";
    stream output = NULL;
    const char *names[4] = {"cos", "sin", "sincos", "cossin"};
    real ***matrices[4] = {gd->histZetaMcos, gd->histZetaMsin,
                           gd->histZetaMsincos, gd->histZetaMcossin};
    for (int c = 0; c < 4; c++)
        for (int m = 0; m <= cmd->mChebyshev; m++)
            if (treecorr_write_edge_matrix(cmd, gd, names[c], m, NULL,
                    matrices[c][m+1], (size_t)cmd->sizeHistN+1) == FAILURE) return FAILURE;
#ifdef TWOPCF
    char xi_path[MAXLENGTHOFFILES + 80];
    if (format_checked(xi_path, sizeof(xi_path), "BALLS4 Xi path", "%s%s%s",
            gd->fpfnamehistXi2pcfFileName, cmd->suffixOutFiles, EXTFILES) != 0)
        return FAILURE;
    char *paths[3] = {gd->fpfnamehistNNFileName, gd->fpfnamehistCFFileName,
                            xi_path};
    real *vectors[3] = {gd->histNN, gd->histCF, gd->histXi2pcf};
    for (int v = 0; v < 3; v++) {
        if ((v == 0 && !cballs_opt_compute_histn(cmd))
            || (v == 1 && !cballs_opt_and_cf(cmd))) continue;
        OPEN_OUTPUT_OR_FAIL(output, paths[v], "w!");
        for (int n = 1; n <= cmd->sizeHistN; n++) {
            real radius = cmd->rminHist == 0
                ? cmd->rangeN * pow(10.0, (n-0.5-cmd->sizeHistN)/cmd->logHistBinsPD)
                : cmd->rminHist * pow(10.0, (n-0.5)*gd->deltaR);
            if (cballs_opt_rbin_arcmin(cmd)) radius *= RADTOARCMIN;
            else if (cballs_opt_rbin_degree(cmd)) radius *= RADTOARCMIN/60.0;
            WRITE_OUTPUT_OR_FAIL(output, paths[v], "%.17g %.17g\n", radius, vectors[v][n]);
        }
        CLOSE_OUTPUT_OR_FAIL(output, paths[v]);
    }
#else
    (void)output; (void)routineName;
#endif
    return SUCCESS;
}

global int balls4_edge_search(struct cmdline_data *cmd, struct global_data *gd,
                              bodyptr *btable, INTEGER *nbody, INTEGER ipmin,
                              INTEGER *ipmax, int cat1, int cat2)
{
    b4_edge_context ctx = {0};
    bodyptr *pivots = NULL;
    real *total = NULL, *batch = NULL, *work = NULL;
    size_t pivot_count = 0;
    INTEGER counters[2] = {0,0};
    INTEGER batch_counts[B4_EDGE_BATCH][2] = {{0}};
    int status = FAILURE, failed = 0;
    const double started = CPUTIME;
    ctx.cmd = cmd; ctx.gd = gd;
    ctx.orders = cmd->mChebyshev + 1;
    ctx.window_orders = treecorr_window_orders(cmd);
    ctx.stride = (size_t)cmd->sizeHistN + 1;
    ctx.weighted = cballs_opt_weights_norm(cmd);
    ctx.aggregate = !cballs_opt_no_one_ball(cmd) && !cballs_opt_no_two_balls(cmd);
    if (ipmin < 1 || ipmax[cat1] < ipmin || ipmax[cat1] > nbody[cat1]
        || !treecorr_triple_values(ctx.stride, ctx.orders, ctx.window_orders,
                                    &ctx.triple_values)) goto allocation_fail;
    ctx.plane = ctx.stride * ctx.stride;
    if (ctx.triple_values > SIZE_MAX - 3*ctx.stride) goto allocation_fail;
    ctx.values = ctx.triple_values + 3*ctx.stride;
    if ((size_t)nbody[cat2] > SIZE_MAX - (size_t)gd->ncellTable[cat2]) goto allocation_fail;
    ctx.capacity = (size_t)nbody[cat2] + (size_t)gd->ncellTable[cat2];
    if (ctx.capacity > SIZE_MAX / sizeof(*ctx.nodes)
        || (size_t)nbody[cat1] > SIZE_MAX / sizeof(*pivots)
        || (size_t)ctx.window_orders > (SIZE_MAX / ctx.stride - 1) / 7)
        goto allocation_fail;
    ctx.work_values = (4*(size_t)ctx.window_orders + 3*(size_t)ctx.orders + 1)*ctx.stride;
    if (ctx.values > SIZE_MAX / B4_EDGE_BATCH / sizeof(real)
        || ctx.work_values > SIZE_MAX / B4_EDGE_BATCH / sizeof(real)) goto allocation_fail;
    ctx.nodes = calloc(ctx.capacity, sizeof(*ctx.nodes));
    pivots = malloc((size_t)nbody[cat1] * sizeof(*pivots));
    total = calloc(ctx.values, sizeof(*total));
    batch = calloc(B4_EDGE_BATCH*ctx.values, sizeof(*batch));
    work = calloc(B4_EDGE_BATCH*ctx.work_values, sizeof(*work));
    if (!ctx.nodes || !pivots || !total || !batch || !work) goto allocation_fail;
    if (b4_edge_copy(&ctx, (nodeptr)roottable[cat2]) == SIZE_MAX) goto allocation_fail;
    for (INTEGER i = 0; i < gd->nnodescanlevTableB4[cat1]; i++)
        if (b4_edge_pivots(cmd, nodetablescanlevB4[cat1][i],
                btable[cat1]+ipmin-1, btable[cat1]+ipmax[cat1],
                pivots, (size_t)nbody[cat1], &pivot_count) == FAILURE)
            goto allocation_fail;
    status = SUCCESS;
    goto allocated;
allocation_fail:
    snprintf(cmd->error_message, _ERRORMSGSIZE_,
             "%s: invalid bounds, size overflow, or edge workspace allocation failure",
             cmd->searchMethod);
allocated:
    status = balls4_consensus(cmd, status, "BALLS4 edge setup");
    if (status == FAILURE) goto cleanup;
    const size_t tasks = pivot_count / B4_EDGE_BLOCK + (pivot_count % B4_EDGE_BLOCK != 0);
    for (size_t first = 0; first < tasks; first += B4_EDGE_BATCH) {
        const int slots = (int)MIN((size_t)B4_EDGE_BATCH, tasks-first);
        memset(batch, 0, B4_EDGE_BATCH*ctx.values*sizeof(*batch));
        memset(batch_counts, 0, sizeof(batch_counts));
        failed = 0;
#pragma omp parallel for schedule(dynamic) reduction(|:failed)
        for (int slot = 0; slot < slots; slot++) {
            const size_t task = first + (size_t)slot;
            if (!balls4_task_owned(cmd, (INTEGER)task)) continue;
            const size_t last = MIN(pivot_count, (task+1)*B4_EDGE_BLOCK);
            for (size_t p = task*B4_EDGE_BLOCK; p < last; p++)
                if (b4_edge_pivot(&ctx, pivots[p], batch+(size_t)slot*ctx.values,
                        work+(size_t)slot*ctx.work_values, batch_counts[slot]) == FAILURE)
                    failed = 1;
        }
        if (failed) snprintf(cmd->error_message, _ERRORMSGSIZE_,
                              "%s: undefined angular reference at a pivot", cmd->searchMethod);
        status = balls4_consensus(cmd, failed ? FAILURE : SUCCESS, "BALLS4 edge pivots");
        if (status == FAILURE) goto cleanup;
        if (balls4_reduce(cmd, batch, (size_t)slots*ctx.values) == FAILURE
            || balls4_reduce_counts(cmd, &batch_counts[0][0], (size_t)slots*2) == FAILURE) {
            status = FAILURE; goto cleanup;
        }
        if (balls4_publish(cmd))
            for (int slot = 0; slot < slots; slot++) {
                for (size_t i = 0; i < ctx.values; i++)
                    total[i] += batch[(size_t)slot*ctx.values+i];
                counters[0] += batch_counts[slot][0]; counters[1] += batch_counts[slot][1];
            }
    }
    if (balls4_publish(cmd)) {
        real ***matrices[4] = {gd->histZetaMcos, gd->histZetaMsin,
                               gd->histZetaMsincos, gd->histZetaMcossin};
        for (int c = 0; c < 4; c++)
            for (int m = 0; m < ctx.orders; m++)
                for (int i = 1; i <= cmd->sizeHistN; i++)
                    for (int j = 1; j <= cmd->sizeHistN; j++)
                        matrices[c][m+1][i][j] =
                            total[((size_t)c*ctx.orders+m)*ctx.plane+(size_t)i*ctx.stride+j];
#ifdef TWOPCF
        const real *pairs = total + ctx.triple_values;
        const real symmetry = cballs_opt_asymmetric(cmd) ? 1.0 : 0.5;
        const real count_scale = cballs_opt_compute_histn(cmd) ? symmetry : 1.0;
        const real valid = (real)pivot_count;
        real volume = 1;
        for (int k = 0; k < NDIM; k++) volume *= gd->Box[k];
        for (int n = 1; n <= cmd->sizeHistN; n++) {
            gd->histNN[n] = count_scale*pairs[n];
            gd->histNNSubXi2pcf[n] = symmetry*pairs[n];
            gd->histXi2pcf[n] = cballs_normalize_or_zero(pairs[2*ctx.stride+n],
                ctx.weighted ? pairs[ctx.stride+n] : pairs[n]);
            gd->histCF[n] = 0;
            if (valid > 0 && gd->histNN[n] > 0 && cballs_opt_and_cf(cmd))
                gd->histCF[n] = gd->histNN[n]*volume /
                    (2.0*PI*gd->deltaR*gd->deltaR*gd->deltaR*valid*valid*(n-0.5)*(n-0.5))-1;
        }
#endif
        gd->nbbcalc = counters[0]; gd->nbccalc = counters[1]; gd->ncccalc = 0;
        if (!cballs_opt_no_out_hist(cmd)) status = b4_edge_output(cmd, gd);
        if (status == SUCCESS && cballs_opt_edge_corrections(cmd))
            status = treecorr_publish_edge(cmd, gd, total, 1, ctx.values, ctx.stride, ctx.orders);
    }
    status = balls4_consensus(cmd, status, "BALLS4 edge publication");
cleanup:
    free(work); free(batch); free(total); free(pivots); free(ctx.nodes);
    gd->flagPrint = FALSE;
    gd->cpusearch = CPUTIME - started;
    return status;
}
#else
global int balls4_edge_search(struct cmdline_data *cmd, struct global_data *gd,
                              bodyptr *btable, INTEGER *nbody, INTEGER ipmin,
                              INTEGER *ipmax, int cat1, int cat2)
{
    (void)gd; (void)btable; (void)nbody; (void)ipmin; (void)ipmax; (void)cat1; (void)cat2;
    cBALLS_FAIL(cmd, "BALLS4 edge corrections require 3D and TPCFON=1\n");
}
#endif
