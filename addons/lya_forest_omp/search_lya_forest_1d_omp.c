/* Radial-only weighted Lyman-alpha forest 2PCF and 3PCF.
 *
 * Catalogs retain their observer-centered 3D coordinates, but this search
 * deliberately ignores angular and transverse separation.  A sorted radial
 * index provides exact one-dimensional range rejection.  Contributions are
 * always evaluated from individual pixels belonging to distinct forests.
 */

#include "globaldefs.h"
#include "lya_forest_defs.h"

#include <errno.h>
#include <limits.h>

#ifndef LYA1D_OMP_PIVOT_BLOCK_SIZE
#define LYA1D_OMP_PIVOT_BLOCK_SIZE 256
#endif
#if LYA1D_OMP_PIVOT_BLOCK_SIZE < 1
#error "LYA1D_OMP_PIVOT_BLOCK_SIZE must be positive"
#endif

typedef struct {
    REAL *num2;
    REAL *den2;
    REAL *num3;
    REAL *den3;
    size_t *touched2;
    size_t touched2_count;
    size_t *touched3;
    size_t touched3_count;
    INTEGER accepted_visits;
    INTEGER pair_count;
    INTEGER ordered_triplet_count;
} lya1d_worker_hist;

local int lya1d_size_mul(size_t left, size_t right, size_t *result)
{
    if (left != 0 && right > SIZE_MAX / left) return FAILURE;
    *result = left * right;
    return SUCCESS;
}

local int lya1d_radial_order(const void *left_ptr, const void *right_ptr)
{
    bodyptr left = *(bodyptr const *)left_ptr;
    bodyptr right = *(bodyptr const *)right_ptr;

    if (LyaDistance(left) < LyaDistance(right)) return -1;
    if (LyaDistance(left) > LyaDistance(right)) return 1;
    if (Id(left) < Id(right)) return -1;
    if (Id(left) > Id(right)) return 1;
    return 0;
}

local int lya1d_positive_bin(REAL value, REAL maximum, int bins)
{
    int bin;

    if (!isfinite(value) || !(value >= 0.0) || value >= maximum) return -1;
    bin = (int)(value / maximum * (REAL)bins);
    return bin >= 0 && bin < bins ? bin : -1;
}

local int lya1d_signed_bin(REAL value, REAL maximum, int bins_per_side)
{
    int bin;
    int total_bins = 2 * bins_per_side;

    if (!isfinite(value) || !(value > -maximum && value < maximum)) return -1;
    bin = (int)((value + maximum) / maximum * (REAL)bins_per_side);
    return bin >= 0 && bin < total_bins ? bin : -1;
}

local int lya1d_worker_init(lya1d_worker_hist *worker,
                            size_t bins2, size_t bins3,
                            int compute_2pcf, int compute_3pcf,
                            ErrorMsg error_message)
{
    memset(worker, 0, sizeof(*worker));

    if (compute_2pcf) {
        if (cballs_calloc_checked((void **)&worker->num2, bins2,
                                  sizeof(*worker->num2),
                                  "radial Ly-alpha 2PCF numerator",
                                  error_message, _ERRORMSGSIZE_) == FAILURE
            || cballs_calloc_checked((void **)&worker->den2, bins2,
                                     sizeof(*worker->den2),
                                     "radial Ly-alpha 2PCF denominator",
                                     error_message, _ERRORMSGSIZE_) == FAILURE
            || cballs_calloc_checked((void **)&worker->touched2, bins2,
                                     sizeof(*worker->touched2),
                                     "radial Ly-alpha 2PCF touched bins",
                                     error_message, _ERRORMSGSIZE_) == FAILURE)
            goto fail;
    }
    if (compute_3pcf) {
        if (cballs_calloc_checked((void **)&worker->num3, bins3,
                                  sizeof(*worker->num3),
                                  "radial Ly-alpha 3PCF numerator",
                                  error_message, _ERRORMSGSIZE_) == FAILURE
            || cballs_calloc_checked((void **)&worker->den3, bins3,
                                     sizeof(*worker->den3),
                                     "radial Ly-alpha 3PCF denominator",
                                     error_message, _ERRORMSGSIZE_) == FAILURE
            || cballs_calloc_checked((void **)&worker->touched3, bins3,
                                     sizeof(*worker->touched3),
                                     "radial Ly-alpha 3PCF touched bins",
                                     error_message, _ERRORMSGSIZE_) == FAILURE)
            goto fail;
    }
    return SUCCESS;

fail:
    free(worker->num2);
    free(worker->den2);
    free(worker->num3);
    free(worker->den3);
    free(worker->touched2);
    free(worker->touched3);
    memset(worker, 0, sizeof(*worker));
    return FAILURE;
}

local void lya1d_worker_free(lya1d_worker_hist *worker)
{
    free(worker->num2);
    free(worker->den2);
    free(worker->num3);
    free(worker->den3);
    free(worker->touched2);
    free(worker->touched3);
    memset(worker, 0, sizeof(*worker));
}

local void lya1d_add_2pcf(struct cmdline_data *cmd, bodyptr p, bodyptr q,
                          lya1d_worker_hist *worker)
{
    REAL separation;
    REAL denominator;
    REAL numerator;
    int bin;

    if (Id(q) <= Id(p) || LyaForestId(q) == LyaForestId(p)) return;
    separation = rabs(LyaDistance(q) - LyaDistance(p));
    bin = lya1d_positive_bin(separation, cmd->lya2RpMax,
                             cmd->lya2RpBins);
    if (bin < 0) return;

    worker->pair_count++;
    denominator = Weight(p) * Weight(q);
    if (denominator <= 0.0) return;
    numerator = (Weight(p) * Kappa(p)) * (Weight(q) * Kappa(q));
    if (worker->den2[bin] == 0.0)
        worker->touched2[worker->touched2_count++] = (size_t)bin;
    worker->num2[bin] += numerator;
    worker->den2[bin] += denominator;
}

local void lya1d_add_3pcf(struct cmdline_data *cmd,
                          bodyptr p, bodyptr q, bodyptr r,
                          lya1d_worker_hist *worker)
{
    REAL lag_q;
    REAL lag_r;
    REAL denominator;
    REAL numerator;
    int bins_per_side = cmd->lya3RBins;
    int total_bins = 2 * bins_per_side;
    int bq;
    int br;
    size_t forward;
    size_t reverse;

    if (LyaForestId(p) == LyaForestId(q)
        || LyaForestId(p) == LyaForestId(r)
        || LyaForestId(q) == LyaForestId(r))
        return;

    lag_q = LyaDistance(q) - LyaDistance(p);
    lag_r = LyaDistance(r) - LyaDistance(p);
    bq = lya1d_signed_bin(lag_q, cmd->lya3RMax, bins_per_side);
    br = lya1d_signed_bin(lag_r, cmd->lya3RMax, bins_per_side);
    if (bq < 0 || br < 0) return;

    worker->ordered_triplet_count += 2;
    denominator = Weight(p) * Weight(q) * Weight(r);
    if (denominator <= 0.0) return;
    numerator = (Weight(p) * Kappa(p))
              * (Weight(q) * Kappa(q))
              * (Weight(r) * Kappa(r));
    forward = (size_t)bq * (size_t)total_bins + (size_t)br;
    reverse = (size_t)br * (size_t)total_bins + (size_t)bq;

    if (worker->den3[forward] == 0.0)
        worker->touched3[worker->touched3_count++] = forward;
    worker->num3[forward] += numerator;
    worker->den3[forward] += denominator;
    if (worker->den3[reverse] == 0.0)
        worker->touched3[worker->touched3_count++] = reverse;
    worker->num3[reverse] += numerator;
    worker->den3[reverse] += denominator;
}

local void lya1d_accumulate_pivot(struct cmdline_data *cmd,
                                  bodyptr *radial_order,
                                  size_t pivot_index,
                                  size_t first, size_t end,
                                  int compute_2pcf, int compute_3pcf,
                                  lya1d_worker_hist *worker)
{
    bodyptr p = radial_order[pivot_index];
    size_t iq;

    for (iq = first; iq < end; iq++) {
        bodyptr q = radial_order[iq];
        REAL lag_q;
        size_t ir;

        if (iq == pivot_index || Update(q) == FALSE
            || Mask(q) != MASK_NODE_VALID)
            continue;
        worker->accepted_visits++;
        if (compute_2pcf) lya1d_add_2pcf(cmd, p, q, worker);
        if (!compute_3pcf || LyaForestId(q) == LyaForestId(p)) continue;

        lag_q = LyaDistance(q) - LyaDistance(p);
        if (!(rabs(lag_q) < cmd->lya3RMax)) continue;
        for (ir = iq + 1; ir < end; ir++) {
            bodyptr r = radial_order[ir];
            REAL lag_r;

            if (ir == pivot_index || Update(r) == FALSE
                || Mask(r) != MASK_NODE_VALID)
                continue;
            lag_r = LyaDistance(r) - LyaDistance(p);
            if (!(rabs(lag_r) < cmd->lya3RMax)) continue;
            lya1d_add_3pcf(cmd, p, q, r, worker);
        }
    }
}

local void lya1d_commit_worker(lya1d_worker_hist *worker,
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

local int lya1d_write_2pcf(struct cmdline_data *cmd,
                            struct global_data *gd,
                            const REAL *num, const REAL *den,
                            INTEGER pair_count)
{
    char path[MAXLENGTHOFFILES];
    FILE *stream;
    int bin;
    int write_failed;

    if (format_checked(path, sizeof(path), "radial Ly-alpha 2PCF path",
                       "%s_lya1d%s", gd->fpfnamehistXi2pcfFileName,
                       EXTFILES) != 0) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "radial Ly-alpha 2PCF output path is too long");
        return FAILURE;
    }
    stream = fopen(path, "w");
    if (stream == NULL) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "cannot open radial Ly-alpha 2PCF output '%s': %s",
                 path, strerror(errno));
        return FAILURE;
    }

    fprintf(stream, "# Radial-only weighted Lyman-alpha forest 2PCF\n");
    fprintf(stream, "# distinct-forest unordered pairs: %" INTEGER_FMT "\n",
            pair_count);
    fprintf(stream, "# transverse separation is ignored\n");
    fprintf(stream, "# columns: bin radial_separation xi numerator denominator\n");
    for (bin = 0; bin < cmd->lya2RpBins; bin++) {
        REAL separation = ((REAL)bin + 0.5) * cmd->lya2RpMax
                        / (REAL)cmd->lya2RpBins;
        fprintf(stream, "%d %.17g %.17g %.17g %.17g\n", bin, separation,
                cballs_normalize_or_zero(num[bin], den[bin]),
                num[bin], den[bin]);
    }
    write_failed = ferror(stream);
    if (fclose(stream) != 0) write_failed = TRUE;
    if (write_failed) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "failed writing radial Ly-alpha 2PCF output '%s'", path);
        return FAILURE;
    }
    return SUCCESS;
}

local int lya1d_write_3pcf(struct cmdline_data *cmd,
                            struct global_data *gd,
                            const REAL *num, const REAL *den,
                            INTEGER ordered_triplet_count)
{
    char path[MAXLENGTHOFFILES];
    FILE *stream;
    int bins_per_side = cmd->lya3RBins;
    int total_bins = 2 * bins_per_side;
    int b1;
    int b2;
    int output_empty = cballs_opt_lya_output_empty_bins(cmd);
    int write_failed;

    if (format_checked(path, sizeof(path), "radial Ly-alpha 3PCF path",
                       "%s_lya1d%s", gd->fpfnamehistZetaMFileName,
                       EXTFILES) != 0) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "radial Ly-alpha 3PCF output path is too long");
        return FAILURE;
    }
    stream = fopen(path, "w");
    if (stream == NULL) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "cannot open radial Ly-alpha 3PCF output '%s': %s",
                 path, strerror(errno));
        return FAILURE;
    }

    fprintf(stream, "# Radial-only weighted Lyman-alpha forest 3PCF\n");
    fprintf(stream, "# distinct-forest ordered triplets: %" INTEGER_FMT "\n",
            ordered_triplet_count);
    fprintf(stream, "# transverse separation is ignored; lags are signed about the pivot\n");
    fprintf(stream, "# zero-denominator policy: zeta=0; empty bins are %s\n",
            output_empty ? "included" : "omitted");
    fprintf(stream, "# columns: b1 b2 lag1 lag2 zeta numerator denominator\n");
    for (b1 = 0; b1 < total_bins; b1++) {
        for (b2 = 0; b2 < total_bins; b2++) {
            size_t index = (size_t)b1 * (size_t)total_bins + (size_t)b2;
            REAL lag1;
            REAL lag2;

            if (!output_empty && den[index] == 0.0) continue;
            lag1 = -cmd->lya3RMax
                 + ((REAL)b1 + 0.5) * cmd->lya3RMax
                 / (REAL)bins_per_side;
            lag2 = -cmd->lya3RMax
                 + ((REAL)b2 + 0.5) * cmd->lya3RMax
                 / (REAL)bins_per_side;
            fprintf(stream, "%d %d %.17g %.17g %.17g %.17g %.17g\n",
                    b1, b2, lag1, lag2,
                    cballs_normalize_or_zero(num[index], den[index]),
                    num[index], den[index]);
        }
    }
    write_failed = ferror(stream);
    if (fclose(stream) != 0) write_failed = TRUE;
    if (write_failed) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "failed writing radial Ly-alpha 3PCF output '%s'", path);
        return FAILURE;
    }
    return SUCCESS;
}

global int searchcalc_lya_forest_1d_omp(struct cmdline_data *cmd,
                                        struct global_data *gd,
                                        bodyptr table, INTEGER nbody,
                                        int compute_2pcf,
                                        int compute_3pcf)
{
    double cpustart = CPUTIME;
    size_t count;
    size_t bins2 = 0;
    size_t bins3 = 0;
    size_t total_radial_bins = 0;
    bodyptr *radial_order = NULL;
    size_t *first = NULL;
    size_t *end = NULL;
    REAL *num2 = NULL;
    REAL *den2 = NULL;
    REAL *num3 = NULL;
    REAL *den3 = NULL;
    REAL cutoff;
    INTEGER accepted_visits = 0;
    INTEGER pair_count = 0;
    INTEGER ordered_triplet_count = 0;
    INTEGER block_count;
    int allocation_failed = FALSE;
    ErrorMsg worker_error = "";
    int status = FAILURE;
    size_t i;
    size_t left_cursor = 0;
    size_t right_cursor = 0;

#if NDIM != 3
    snprintf(cmd->error_message, _ERRORMSGSIZE_,
             "radial Ly-alpha correlations require 3D catalog storage");
    return FAILURE;
#endif
#ifndef OPENMPCODE
    snprintf(cmd->error_message, _ERRORMSGSIZE_,
             "radial Ly-alpha addon requires OPENMPMACHINE=1");
    return FAILURE;
#endif

    if (nbody < 1) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "radial Ly-alpha search received an empty catalog");
        return FAILURE;
    }
    if (compute_3pcf && cmd->lya3RBins > INT_MAX / 2) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "radial Ly-alpha 3PCF bin count is too large");
        return FAILURE;
    }
    block_count = 1 + (nbody - 1) / (INTEGER)LYA1D_OMP_PIVOT_BLOCK_SIZE;

    count = (size_t)nbody;
    if (compute_2pcf) bins2 = (size_t)cmd->lya2RpBins;
    if (compute_3pcf) {
        total_radial_bins = 2 * (size_t)cmd->lya3RBins;
        if (lya1d_size_mul(total_radial_bins, total_radial_bins,
                           &bins3) == FAILURE)
            goto size_error;
    }

    if (cballs_calloc_checked((void **)&radial_order, count,
                              sizeof(*radial_order),
                              "radial Ly-alpha sorted index",
                              cmd->error_message,
                              sizeof(cmd->error_message)) == FAILURE
        || cballs_calloc_checked((void **)&first, count, sizeof(*first),
                                 "radial Ly-alpha lower limits",
                                 cmd->error_message,
                                 sizeof(cmd->error_message)) == FAILURE
        || cballs_calloc_checked((void **)&end, count, sizeof(*end),
                                 "radial Ly-alpha upper limits",
                                 cmd->error_message,
                                 sizeof(cmd->error_message)) == FAILURE)
        goto cleanup;
    if (compute_2pcf
        && (cballs_calloc_checked((void **)&num2, bins2, sizeof(*num2),
                                  "global radial Ly-alpha 2PCF numerator",
                                  cmd->error_message,
                                  sizeof(cmd->error_message)) == FAILURE
            || cballs_calloc_checked((void **)&den2, bins2, sizeof(*den2),
                                     "global radial Ly-alpha 2PCF denominator",
                                     cmd->error_message,
                                     sizeof(cmd->error_message)) == FAILURE))
        goto cleanup;
    if (compute_3pcf
        && (cballs_calloc_checked((void **)&num3, bins3, sizeof(*num3),
                                  "global radial Ly-alpha 3PCF numerator",
                                  cmd->error_message,
                                  sizeof(cmd->error_message)) == FAILURE
            || cballs_calloc_checked((void **)&den3, bins3, sizeof(*den3),
                                     "global radial Ly-alpha 3PCF denominator",
                                     cmd->error_message,
                                     sizeof(cmd->error_message)) == FAILURE))
        goto cleanup;

    for (i = 0; i < count; i++) radial_order[i] = table + i;
    qsort(radial_order, count, sizeof(*radial_order), lya1d_radial_order);
    cutoff = MAX(compute_2pcf ? cmd->lya2RpMax : 0.0,
                 compute_3pcf ? cmd->lya3RMax : 0.0);
    for (i = 0; i < count; i++) {
        REAL pivot_distance = LyaDistance(radial_order[i]);

        while (left_cursor < i
               && pivot_distance - LyaDistance(radial_order[left_cursor])
                  >= cutoff)
            left_cursor++;
        if (right_cursor < i + 1) right_cursor = i + 1;
        while (right_cursor < count
               && LyaDistance(radial_order[right_cursor]) - pivot_distance
                  < cutoff)
            right_cursor++;
        first[i] = left_cursor;
        end[i] = right_cursor;
    }

    ThreadCount(cmd, gd, nbody, 0);
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
        "\n%s: exact radial-only Ly-alpha estimator; 2PCF=%d 3PCF=%d "
        "cutoff=%g pivot_block=%d blocks=%" INTEGER_FMT "\n",
        cmd->searchMethod, compute_2pcf, compute_3pcf, cutoff,
        LYA1D_OMP_PIVOT_BLOCK_SIZE, block_count);

#pragma omp parallel shared(allocation_failed,worker_error,num2,den2,num3,den3,accepted_visits,pair_count,ordered_triplet_count)
    {
        lya1d_worker_hist worker;
        ErrorMsg local_error = "";
        int worker_ready = lya1d_worker_init(
            &worker, bins2, bins3, compute_2pcf, compute_3pcf,
            local_error) == SUCCESS;
        int worker_failed = !worker_ready;
        INTEGER block_index;

        if (!worker_ready) {
#pragma omp critical(lya1d_failure)
            {
                if (!allocation_failed)
                    snprintf(worker_error, sizeof(worker_error), "%s",
                             local_error);
                allocation_failed = TRUE;
            }
        }

#pragma omp barrier
#pragma omp for schedule(static,1) ordered
        for (block_index = 0; block_index < block_count; block_index++) {
            size_t block_first = (size_t)block_index
                               * (size_t)LYA1D_OMP_PIVOT_BLOCK_SIZE;
            size_t block_end = MIN(
                block_first + (size_t)LYA1D_OMP_PIVOT_BLOCK_SIZE, count);
            size_t pivot;

            worker.accepted_visits = 0;
            worker.pair_count = 0;
            worker.ordered_triplet_count = 0;
            if (!worker_failed && !allocation_failed) {
                for (pivot = block_first; pivot < block_end; pivot++) {
                    bodyptr p = radial_order[pivot];

                    if (Update(p) != FALSE && Mask(p) == MASK_NODE_VALID)
                        lya1d_accumulate_pivot(
                            cmd, radial_order, pivot,
                            first[pivot], end[pivot],
                            compute_2pcf, compute_3pcf, &worker);
                }
            }

#pragma omp ordered
            {
                if (worker_ready && !allocation_failed && !worker_failed) {
                    lya1d_commit_worker(&worker, num2, den2, num3, den3,
                                        compute_2pcf, compute_3pcf);
                    accepted_visits += worker.accepted_visits;
                    pair_count += worker.pair_count;
                    ordered_triplet_count += worker.ordered_triplet_count;
                }
            }
        }
        if (worker_ready) lya1d_worker_free(&worker);
    }

    if (allocation_failed) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "radial Ly-alpha worker failed: %.2000s",
                 worker_error[0] != '\0' ? worker_error : "unknown error");
        goto cleanup;
    }

    gd->nbbcalc += accepted_visits;
    if (!cballs_opt_no_out_hist(cmd)) {
        if (compute_2pcf
            && lya1d_write_2pcf(cmd, gd, num2, den2,
                                pair_count) == FAILURE)
            goto cleanup;
        if (compute_3pcf
            && lya1d_write_3pcf(cmd, gd, num3, den3,
                                ordered_triplet_count) == FAILURE)
            goto cleanup;
    }

    gd->cpusearch = CPUTIME - cpustart;
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
        "%s: radial candidates=%" INTEGER_FMT " pairs=%" INTEGER_FMT
        " ordered_triplets=%" INTEGER_FMT " CPU=%g\n",
        cmd->searchMethod, accepted_visits, pair_count,
        ordered_triplet_count, gd->cpusearch);
    status = SUCCESS;
    goto cleanup;

size_error:
    snprintf(cmd->error_message, _ERRORMSGSIZE_,
             "radial Ly-alpha histogram dimensions overflow size_t");

cleanup:
    free(radial_order);
    free(first);
    free(end);
    free(num2);
    free(den2);
    free(num3);
    free(den3);
    return status;
}
