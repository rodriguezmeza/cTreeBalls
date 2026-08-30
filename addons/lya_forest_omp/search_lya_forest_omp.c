/* Exact weighted Lyman-alpha forest 2PCF and anisotropic 3PCF.
 *
 * This implements equations (2.4) and (2.8) of the Lya2pcf paper.  Octree
 * cells are used only to reject regions outside the largest requested scale;
 * all accepted contributions are evaluated from individual forest pixels.
 */

#include "globaldefs.h"
#include "lya_forest_defs.h"

#include <errno.h>
#include <float.h>

#define LYA_PI 3.141592653589793238462643383279502884

typedef struct {
    REAL *num2;
    REAL *den2;
    REAL *num3;
    REAL *den3;
    size_t *touched2;
    size_t touched2_count;
    size_t *touched3;
    size_t touched3_count;
    bodyptr *neighbors;
    size_t neighbor_count;
    size_t neighbor_capacity;
    INTEGER accepted_visits;
    INTEGER pair_count;
    INTEGER ordered_triplet_count;
} lya_worker_hist;

local int lya_size_mul(size_t left, size_t right, size_t *result)
{
    if (left != 0 && right > SIZE_MAX / left) return FAILURE;
    *result = left * right;
    return SUCCESS;
}

local REAL lya_clamp(REAL value, REAL lower, REAL upper)
{
    return value < lower ? lower : (value > upper ? upper : value);
}

local int lya_bin_positive(REAL value, REAL maximum, int bins)
{
    int bin;
    if (!(value >= 0.0) || value >= maximum) return -1;
    bin = (int)(value / maximum * (REAL)bins);
    return bin >= 0 && bin < bins ? bin : -1;
}

local int lya_bin_theta(REAL theta, int bins)
{
    int bin;
    theta = lya_clamp(theta, 0.0, (REAL)LYA_PI);
    bin = (int)(theta / (REAL)LYA_PI * (REAL)bins);
    return bin == bins ? bins - 1 : bin;
}

local int lya_bin_mu(REAL mu, int bins)
{
    int bin;
    mu = lya_clamp(mu, -1.0, 1.0);
    bin = (int)((mu + 1.0) * 0.5 * (REAL)bins);
    return bin == bins ? bins - 1 : bin;
}

local size_t lya_index3(int b1, int b2, int t1, int t2, int mu,
                        int radial_bins, int theta_bins, int mu_bins)
{
    size_t index = (size_t)b1;
    index = index * (size_t)radial_bins + (size_t)b2;
    index = index * (size_t)theta_bins + (size_t)t1;
    index = index * (size_t)theta_bins + (size_t)t2;
    index = index * (size_t)mu_bins + (size_t)mu;
    return index;
}

local int lya_worker_init(struct cmdline_data *cmd, lya_worker_hist *worker,
                          size_t bins2, size_t bins3, int compute_2pcf,
                          int compute_3pcf, ErrorMsg error_message)
{
    memset(worker, 0, sizeof(*worker));

    if (compute_2pcf) {
        if (cballs_calloc_checked((void **)&worker->num2, bins2,
                                  sizeof(*worker->num2), "Ly-alpha 2PCF numerator",
                                  error_message, _ERRORMSGSIZE_) == FAILURE
            || cballs_calloc_checked((void **)&worker->den2, bins2,
                                     sizeof(*worker->den2), "Ly-alpha 2PCF denominator",
                                     error_message, _ERRORMSGSIZE_) == FAILURE
            || cballs_calloc_checked((void **)&worker->touched2, bins2,
                                     sizeof(*worker->touched2), "Ly-alpha 2PCF touched bins",
                                     error_message, _ERRORMSGSIZE_) == FAILURE)
            goto fail;
    }
    if (compute_3pcf) {
        if (cballs_calloc_checked((void **)&worker->num3, bins3,
                                  sizeof(*worker->num3), "Ly-alpha 3PCF numerator",
                                  error_message, _ERRORMSGSIZE_) == FAILURE
            || cballs_calloc_checked((void **)&worker->den3, bins3,
                                     sizeof(*worker->den3), "Ly-alpha 3PCF denominator",
                                     error_message, _ERRORMSGSIZE_) == FAILURE
            || cballs_calloc_checked((void **)&worker->touched3, bins3,
                                     sizeof(*worker->touched3), "Ly-alpha 3PCF touched bins",
                                     error_message, _ERRORMSGSIZE_) == FAILURE)
            goto fail;
    }
    return SUCCESS;

fail:
    free(worker->num2); worker->num2 = NULL;
    free(worker->den2); worker->den2 = NULL;
    free(worker->num3); worker->num3 = NULL;
    free(worker->den3); worker->den3 = NULL;
    free(worker->touched2); worker->touched2 = NULL;
    free(worker->touched3); worker->touched3 = NULL;
    return FAILURE;
}

local void lya_worker_free(lya_worker_hist *worker)
{
    free(worker->num2);
    free(worker->den2);
    free(worker->num3);
    free(worker->den3);
    free(worker->touched2);
    free(worker->touched3);
    free(worker->neighbors);
    memset(worker, 0, sizeof(*worker));
}

local int lya_append_neighbor(lya_worker_hist *worker, bodyptr neighbor,
                              ErrorMsg error_message)
{
    if (worker->neighbor_count == worker->neighbor_capacity) {
        size_t new_capacity = worker->neighbor_capacity == 0
                            ? 128 : worker->neighbor_capacity * 2;
        bodyptr *resized;
        if (new_capacity < worker->neighbor_capacity
            || new_capacity > SIZE_MAX / sizeof(*resized)) {
            snprintf(error_message, _ERRORMSGSIZE_,
                     "Ly-alpha neighbor-list size overflow");
            return FAILURE;
        }
        resized = realloc(worker->neighbors,
                          new_capacity * sizeof(*resized));
        if (resized == NULL) {
            snprintf(error_message, _ERRORMSGSIZE_,
                     "not enough memory growing Ly-alpha neighbor list");
            return FAILURE;
        }
        worker->neighbors = resized;
        worker->neighbor_capacity = new_capacity;
    }
    worker->neighbors[worker->neighbor_count++] = neighbor;
    return SUCCESS;
}

local void lya_accumulate_2pcf(struct cmdline_data *cmd, bodyptr p, bodyptr q,
                               lya_worker_hist *worker)
{
    REAL cos_sight;
    REAL cos_half;
    REAL sin_half;
    REAL rp;
    REAL rt;
    REAL numerator;
    REAL denominator;
    int bp;
    int bt;
    size_t index;

    if (Id(q) <= Id(p) || LyaForestId(q) == LyaForestId(p)) return;

    DOTVP(cos_sight, LyaLOS(p), LyaLOS(q));
    cos_sight = lya_clamp(cos_sight, -1.0, 1.0);
    cos_half = rsqrt(MAX(0.0, 0.5 * (1.0 + cos_sight)));
    sin_half = rsqrt(MAX(0.0, 0.5 * (1.0 - cos_sight)));
    rp = rabs(LyaDistance(p) - LyaDistance(q)) * cos_half;
    rt = (LyaDistance(p) + LyaDistance(q)) * sin_half;
    bp = lya_bin_positive(rp, cmd->lya2RpMax, cmd->lya2RpBins);
    bt = lya_bin_positive(rt, cmd->lya2RtMax, cmd->lya2RtBins);
    if (bp < 0 || bt < 0) return;

    index = (size_t)bp * (size_t)cmd->lya2RtBins + (size_t)bt;
    denominator = Weight(p) * Weight(q);
    numerator = (Weight(p) * Kappa(p)) * (Weight(q) * Kappa(q));
    worker->pair_count++;
    if (denominator <= 0.0) return;
    if (worker->den2[index] == 0.0)
        worker->touched2[worker->touched2_count++] = index;
    worker->den2[index] += denominator;
    worker->num2[index] += numerator;
}

local int lya_walktree(struct cmdline_data *cmd, bodyptr p, nodeptr q,
                       REAL cutoff, int compute_2pcf, int compute_3pcf,
                       lya_worker_hist *worker, ErrorMsg error_message)
{
    REAL distance2;
    REAL distance;
    compute_vector displacement;
    nodeptr child;

    if ((nodeptr)p == q) return SUCCESS;
    DOTPSUBV(distance2, displacement, Pos(p), Pos(q));
    distance = rsqrt(distance2);

    if (Type(q) == CELL) {
        if (distance >= cutoff + Radius(q)) return SUCCESS;
        for (child = More(q); child != Next(q); child = Next(child)) {
            if (lya_walktree(cmd, p, child, cutoff, compute_2pcf,
                             compute_3pcf, worker, error_message) == FAILURE)
                return FAILURE;
        }
        return SUCCESS;
    }

    if (distance >= cutoff || Update(q) == FALSE
        || Mask(q) != MASK_NODE_VALID)
        return SUCCESS;

    worker->accepted_visits++;
    if (compute_2pcf)
        lya_accumulate_2pcf(cmd, p, (bodyptr)q, worker);
    if (compute_3pcf && distance < cmd->lya3RMax
        && LyaForestId(p) != LyaForestId((bodyptr)q))
        return lya_append_neighbor(worker, (bodyptr)q, error_message);
    return SUCCESS;
}

local void lya_accumulate_3pcf(struct cmdline_data *cmd, bodyptr p,
                               lya_worker_hist *worker)
{
    size_t iq;
    size_t ir;
    for (iq = 0; iq < worker->neighbor_count; iq++) {
        bodyptr q = worker->neighbors[iq];
        compute_vector dq;
        REAL rq2;
        REAL rq;
        REAL cos_theta_q;
        int bq;
        int tq;
        DOTPSUBV(rq2, dq, Pos(q), Pos(p));
        rq = rsqrt(rq2);
        bq = lya_bin_positive(rq, cmd->lya3RMax, cmd->lya3RBins);
        if (bq < 0 || rq <= 0.0) continue;
        DOTVP(cos_theta_q, dq, LyaLOS(p));
        tq = lya_bin_theta(racos(lya_clamp(cos_theta_q / rq, -1.0, 1.0)),
                           cmd->lya3ThetaBins);

        for (ir = iq + 1; ir < worker->neighbor_count; ir++) {
            bodyptr r = worker->neighbors[ir];
            compute_vector dr;
            REAL rr2;
            REAL rr;
            REAL cos_theta_r;
            REAL dot_qr;
            REAL mu;
            REAL numerator;
            REAL denominator;
            int br;
            int tr;
            int bm;
            size_t forward;
            size_t reverse;

            if (LyaForestId(q) == LyaForestId(r)) continue;
            DOTPSUBV(rr2, dr, Pos(r), Pos(p));
            rr = rsqrt(rr2);
            br = lya_bin_positive(rr, cmd->lya3RMax, cmd->lya3RBins);
            if (br < 0 || rr <= 0.0) continue;
            DOTVP(cos_theta_r, dr, LyaLOS(p));
            tr = lya_bin_theta(racos(lya_clamp(cos_theta_r / rr, -1.0, 1.0)),
                               cmd->lya3ThetaBins);
            DOTVP(dot_qr, dq, dr);
            mu = lya_clamp(dot_qr / (rq * rr), -1.0, 1.0);
            bm = lya_bin_mu(mu, cmd->lya3MuBins);

            denominator = Weight(p) * Weight(q) * Weight(r);
            numerator = (Weight(p) * Kappa(p))
                      * (Weight(q) * Kappa(q))
                      * (Weight(r) * Kappa(r));
            worker->ordered_triplet_count += 2;
            if (denominator <= 0.0) continue;
            forward = lya_index3(bq, br, tq, tr, bm,
                                 cmd->lya3RBins, cmd->lya3ThetaBins,
                                 cmd->lya3MuBins);
            reverse = lya_index3(br, bq, tr, tq, bm,
                                 cmd->lya3RBins, cmd->lya3ThetaBins,
                                 cmd->lya3MuBins);
            if (worker->den3[forward] == 0.0)
                worker->touched3[worker->touched3_count++] = forward;
            worker->num3[forward] += numerator;
            worker->den3[forward] += denominator;
            if (worker->den3[reverse] == 0.0)
                worker->touched3[worker->touched3_count++] = reverse;
            worker->num3[reverse] += numerator;
            worker->den3[reverse] += denominator;
        }
    }
}

local void lya_commit_worker(lya_worker_hist *worker,
                             REAL *num2, REAL *den2,
                             REAL *num3, REAL *den3,
                             int compute_2pcf, int compute_3pcf)
{
    size_t i;
    if (compute_2pcf) {
        for (i = 0; i < worker->touched2_count; i++) {
            size_t index = worker->touched2[i];
            num2[index] += worker->num2[index];
            den2[index] += worker->den2[index];
            worker->num2[index] = 0.0;
            worker->den2[index] = 0.0;
        }
        worker->touched2_count = 0;
    }
    if (compute_3pcf) {
        for (i = 0; i < worker->touched3_count; i++) {
            size_t index = worker->touched3[i];
            num3[index] += worker->num3[index];
            den3[index] += worker->den3[index];
            worker->num3[index] = 0.0;
            worker->den3[index] = 0.0;
        }
        worker->touched3_count = 0;
    }
}

local int lya_write_2pcf(struct cmdline_data *cmd, struct global_data *gd,
                         const REAL *num, const REAL *den,
                         INTEGER pair_count)
{
    char path[MAXLENGTHOFFILES];
    FILE *stream = NULL;
    int write_failed;
    int bp;
    int bt;

    if (format_checked(path, sizeof(path), "Ly-alpha 2PCF path", "%s_lya%s",
                       gd->fpfnamehistXi2pcfFileName, EXTFILES) != 0) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "Ly-alpha 2PCF output path is too long");
        return FAILURE;
    }
    stream = fopen(path, "w");
    if (stream == NULL) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "cannot open Ly-alpha 2PCF output '%s': %s",
                 path, strerror(errno));
        return FAILURE;
    }
    fprintf(stream, "# Weighted Lyman-alpha forest 2PCF (paper equation 2.4)\n");
    fprintf(stream, "# distinct-forest unordered pairs: %" INTEGER_FMT "\n",
            pair_count);
    fprintf(stream, "# columns: bp bt rp_center rt_center xi numerator denominator\n");
    for (bp = 0; bp < cmd->lya2RpBins; bp++) {
        for (bt = 0; bt < cmd->lya2RtBins; bt++) {
            size_t index = (size_t)bp * (size_t)cmd->lya2RtBins + (size_t)bt;
            REAL rp = ((REAL)bp + 0.5) * cmd->lya2RpMax
                    / (REAL)cmd->lya2RpBins;
            REAL rt = ((REAL)bt + 0.5) * cmd->lya2RtMax
                    / (REAL)cmd->lya2RtBins;
            fprintf(stream, "%d %d %.17g %.17g %.17g %.17g %.17g\n",
                    bp, bt, rp, rt,
                    cballs_normalize_or_zero(num[index], den[index]),
                    num[index], den[index]);
        }
    }
    write_failed = ferror(stream);
    if (fclose(stream) != 0) write_failed = TRUE;
    if (write_failed) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "failed writing Ly-alpha 2PCF output '%s'", path);
        return FAILURE;
    }
    return SUCCESS;
}

local int lya_write_3pcf(struct cmdline_data *cmd, struct global_data *gd,
                         const REAL *num, const REAL *den,
                         INTEGER ordered_triplet_count)
{
    char path[MAXLENGTHOFFILES];
    FILE *stream = NULL;
    int write_failed;
    int b1, b2, t1, t2, bm;
    int output_empty = cballs_opt_lya_output_empty_bins(cmd);

    if (format_checked(path, sizeof(path), "Ly-alpha 3PCF path", "%s_lya5d%s",
                       gd->fpfnamehistZetaMFileName, EXTFILES) != 0) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "Ly-alpha 3PCF output path is too long");
        return FAILURE;
    }
    stream = fopen(path, "w");
    if (stream == NULL) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "cannot open Ly-alpha 3PCF output '%s': %s",
                 path, strerror(errno));
        return FAILURE;
    }
    fprintf(stream, "# Weighted anisotropic Lyman-alpha forest 3PCF (paper equation 2.8)\n");
    fprintf(stream, "# distinct-forest ordered triplets: %" INTEGER_FMT "\n",
            ordered_triplet_count);
    fprintf(stream, "# zero-denominator policy: zeta=0; empty bins are %s\n",
            output_empty ? "included" : "omitted");
    fprintf(stream, "# columns: b1 b2 t1 t2 bmu r1 r2 theta1 theta2 mu zeta numerator denominator\n");
    for (b1 = 0; b1 < cmd->lya3RBins; b1++)
        for (b2 = 0; b2 < cmd->lya3RBins; b2++)
            for (t1 = 0; t1 < cmd->lya3ThetaBins; t1++)
                for (t2 = 0; t2 < cmd->lya3ThetaBins; t2++)
                    for (bm = 0; bm < cmd->lya3MuBins; bm++) {
                        size_t index = lya_index3(
                            b1, b2, t1, t2, bm, cmd->lya3RBins,
                            cmd->lya3ThetaBins, cmd->lya3MuBins);
                        REAL r1, r2, theta1, theta2, mu;
                        if (!output_empty && den[index] == 0.0) continue;
                        r1 = ((REAL)b1 + 0.5) * cmd->lya3RMax
                           / (REAL)cmd->lya3RBins;
                        r2 = ((REAL)b2 + 0.5) * cmd->lya3RMax
                           / (REAL)cmd->lya3RBins;
                        theta1 = ((REAL)t1 + 0.5) * (REAL)LYA_PI
                               / (REAL)cmd->lya3ThetaBins;
                        theta2 = ((REAL)t2 + 0.5) * (REAL)LYA_PI
                               / (REAL)cmd->lya3ThetaBins;
                        mu = -1.0 + ((REAL)bm + 0.5) * 2.0
                           / (REAL)cmd->lya3MuBins;
                        fprintf(stream,
                            "%d %d %d %d %d %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g\n",
                            b1, b2, t1, t2, bm, r1, r2, theta1, theta2,
                            mu, cballs_normalize_or_zero(num[index], den[index]),
                            num[index], den[index]);
                    }
    write_failed = ferror(stream);
    if (fclose(stream) != 0) write_failed = TRUE;
    if (write_failed) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "failed writing Ly-alpha 3PCF output '%s'", path);
        return FAILURE;
    }
    return SUCCESS;
}

global int searchcalc_lya_forest_omp(struct cmdline_data *cmd,
                                     struct global_data *gd,
                                     bodyptr *btable, INTEGER *nbody,
                                     INTEGER ipmin, INTEGER *ipmax,
                                     int cat, int compute_2pcf,
                                     int compute_3pcf)
{
    double cpustart = CPUTIME;
    size_t bins2 = 0;
    size_t bins3 = 0;
    REAL *num2 = NULL;
    REAL *den2 = NULL;
    REAL *num3 = NULL;
    REAL *den3 = NULL;
    REAL cutoff2 = 0.0;
    REAL cutoff;
    INTEGER accepted_visits = 0;
    INTEGER pair_count = 0;
    INTEGER ordered_triplet_count = 0;
    int allocation_failed = FALSE;
    ErrorMsg worker_error = "";
    int status = FAILURE;

#if NDIM != 3
    snprintf(cmd->error_message, _ERRORMSGSIZE_,
             "Ly-alpha forest correlations require NDIM=3");
    return FAILURE;
#endif
#ifndef OPENMPCODE
    snprintf(cmd->error_message, _ERRORMSGSIZE_,
             "Ly-alpha forest addon requires OPENMPMACHINE=1");
    return FAILURE;
#endif

    if (compute_2pcf
        && lya_size_mul((size_t)cmd->lya2RpBins,
                        (size_t)cmd->lya2RtBins, &bins2) == FAILURE)
        goto size_error;
    if (compute_3pcf) {
        bins3 = (size_t)cmd->lya3RBins;
        if (lya_size_mul(bins3, (size_t)cmd->lya3RBins, &bins3) == FAILURE
            || lya_size_mul(bins3, (size_t)cmd->lya3ThetaBins, &bins3) == FAILURE
            || lya_size_mul(bins3, (size_t)cmd->lya3ThetaBins, &bins3) == FAILURE
            || lya_size_mul(bins3, (size_t)cmd->lya3MuBins, &bins3) == FAILURE)
            goto size_error;
    }

    if (compute_2pcf) {
        if (cballs_calloc_checked((void **)&num2, bins2, sizeof(*num2),
                                  "global Ly-alpha 2PCF numerator",
                                  cmd->error_message,
                                  sizeof(cmd->error_message)) == FAILURE
            || cballs_calloc_checked((void **)&den2, bins2, sizeof(*den2),
                                     "global Ly-alpha 2PCF denominator",
                                     cmd->error_message,
                                     sizeof(cmd->error_message)) == FAILURE)
            goto cleanup;
        cutoff2 = hypot(cmd->lya2RpMax, cmd->lya2RtMax);
    }
    if (compute_3pcf) {
        if (cballs_calloc_checked((void **)&num3, bins3, sizeof(*num3),
                                  "global Ly-alpha 3PCF numerator",
                                  cmd->error_message,
                                  sizeof(cmd->error_message)) == FAILURE
            || cballs_calloc_checked((void **)&den3, bins3, sizeof(*den3),
                                     "global Ly-alpha 3PCF denominator",
                                     cmd->error_message,
                                     sizeof(cmd->error_message)) == FAILURE)
            goto cleanup;
    }
    cutoff = MAX(cutoff2, compute_3pcf ? cmd->lya3RMax : 0.0);

    ThreadCount(cmd, gd, nbody[cat], cat);
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
        "\n%s: exact Ly-alpha estimator; 2PCF=%d 3PCF=%d cutoff=%g\n",
        cmd->searchMethod, compute_2pcf, compute_3pcf, cutoff);

#pragma omp parallel shared(allocation_failed,worker_error,num2,den2,num3,den3,accepted_visits,pair_count,ordered_triplet_count)
    {
        lya_worker_hist worker;
        ErrorMsg local_error = "";
        int worker_ready = lya_worker_init(cmd, &worker, bins2, bins3,
                                           compute_2pcf, compute_3pcf,
                                           local_error) == SUCCESS;
        int worker_failed = !worker_ready;
        bodyptr p;
        if (!worker_ready) {
#pragma omp critical(lya_failure)
            {
                if (!allocation_failed)
                    snprintf(worker_error, sizeof(worker_error), "%s", local_error);
                allocation_failed = TRUE;
            }
        }

#pragma omp barrier
#pragma omp for schedule(static,1) ordered
        for (p = btable[cat] + ipmin - 1;
             p < btable[cat] + ipmax[cat]; p++) {
            worker.accepted_visits = 0;
            worker.pair_count = 0;
            worker.ordered_triplet_count = 0;
            worker.neighbor_count = 0;
            if (!worker_failed && Update(p) != FALSE
                && Mask(p) == MASK_NODE_VALID) {
                if (lya_walktree(cmd, p, (nodeptr)roottable[cat], cutoff,
                                 compute_2pcf, compute_3pcf, &worker,
                                 local_error) == FAILURE) {
                    worker_failed = TRUE;
#pragma omp critical(lya_failure)
                    {
                        if (!allocation_failed)
                            snprintf(worker_error, sizeof(worker_error), "%s",
                                     local_error);
                        allocation_failed = TRUE;
                    }
                } else if (compute_3pcf) {
                    lya_accumulate_3pcf(cmd, p, &worker);
                }
            }

#pragma omp ordered
            {
                if (worker_ready && !allocation_failed && !worker_failed) {
                    lya_commit_worker(&worker, num2, den2, num3, den3,
                                      compute_2pcf, compute_3pcf);
                    accepted_visits += worker.accepted_visits;
                    pair_count += worker.pair_count;
                    ordered_triplet_count += worker.ordered_triplet_count;
                }
            }
        }
        if (worker_ready) lya_worker_free(&worker);
    }

    if (allocation_failed) {
        if (cmd->error_message[0] == '\0')
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     "Ly-alpha worker failed: %s",
                     worker_error[0] != '\0' ? worker_error : "unknown error");
        goto cleanup;
    }

    gd->nbbcalc += accepted_visits;
    if (!cballs_opt_no_out_hist(cmd)) {
        if (compute_2pcf
            && lya_write_2pcf(cmd, gd, num2, den2, pair_count) == FAILURE)
            goto cleanup;
        if (compute_3pcf
            && lya_write_3pcf(cmd, gd, num3, den3,
                              ordered_triplet_count) == FAILURE)
            goto cleanup;
    }

    gd->cpusearch = CPUTIME - cpustart;
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
        "%s: accepted=%" INTEGER_FMT " pairs=%" INTEGER_FMT
        " ordered_triplets=%" INTEGER_FMT " CPU=%g\n",
        cmd->searchMethod, accepted_visits, pair_count,
        ordered_triplet_count, gd->cpusearch);
    status = SUCCESS;
    goto cleanup;

size_error:
    snprintf(cmd->error_message, _ERRORMSGSIZE_,
             "Ly-alpha histogram dimensions overflow size_t");

cleanup:
    free(num2);
    free(den2);
    free(num3);
    free(den3);
    return status;
}

#undef LYA_PI
