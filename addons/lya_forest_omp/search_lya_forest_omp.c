/* Exact weighted Lyman-alpha forest 2PCF and anisotropic 3PCF.
 *
 * This implements equations (2.4) and (2.8) of the Lya2pcf paper.  Octree
 * cells reject regions outside the largest requested scale in the reference
 * walker. Optional persistent 2PCF cell pairs certify rp/rt bins and reuse
 * observer-angle bounds, with exact radial-window and pixel fallbacks.
 * 3PCF contributions use certified same-bin forest segment sums with exact
 * pixel fallback. The opt-in anisotropic multipole method computes exact
 * Legendre moments and explicitly approximate mu-bin reconstruction.
 */

#include "globaldefs.h"
#include "lya_forest_defs.h"
#include "lya_forest_parallel.h"
#include "lya_forest_los_tree.h"

#include <errno.h>
#include <float.h>

#define LYA_PI 3.141592653589793238462643383279502884

/* Pivot-dependent quantities are invariant across the neighbor-pair loop.
 * Keep its inputs contiguous rather than repeatedly dereferencing full bodies.
 * Retain displacement/radius (not unit vectors) to preserve the original mu
 * arithmetic at bin boundaries. */
typedef struct {
    compute_vector displacement;
    REAL radius, weight, weighted_delta;
    INTEGER forest_id;
    size_t first_index, second_index, ordinal, leg_bin;
} lya_neighbor;

#include "lya_triplet_types.h"

typedef struct {
    REAL *moments;
    size_t *moment_bins;
    unsigned char *moment_seen;
    int aggregation_safe, moment_mode;
    REAL polar_edges[65];
    int polar_lookup;
    lya_segment *segments;
    size_t *segment_roots, segment_count, segment_capacity;
    unsigned long long aggregated_pairs, direct_pairs, segment_accepts, approximate_pairs;
    REAL *num2;
    REAL *den2;
    REAL *num3;
    REAL *den3;
    size_t *touched2;
    size_t touched2_count;
    size_t *touched3;
    size_t touched3_count;
    lya_neighbor *neighbors;
    size_t neighbor_count;
    size_t neighbor_capacity;
    size_t histogram_plan_bytes;
    size_t scratch_workers;
    /* Zero multiplicity means the unmodified pixel-pivot path. */
    INTEGER pivot_multiplicity;
    REAL pivot_weight, pivot_field;
    bodyptr *pivot_members;
    INTEGER pivot_idlo, pivot_idhi;
    INTEGER accepted_visits;
    INTEGER pair_count;
    INTEGER ordered_triplet_count;
} lya_worker_hist;

local REAL lya_pivot_weight(const lya_worker_hist *w,bodyptr p)
{return w->pivot_multiplicity?w->pivot_weight:Weight(p);}

local REAL lya_pivot_field(const lya_worker_hist *w,bodyptr p)
{return w->pivot_multiplicity?w->pivot_field:Weight(p)*Kappa(p);}

local INTEGER lya_pivot_multiplicity(const lya_worker_hist *w)
{return w->pivot_multiplicity?w->pivot_multiplicity:1;}

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
                          int compute_3pcf, size_t histogram_plan_bytes,
                          ErrorMsg error_message)
{
    memset(worker, 0, sizeof(*worker));
    worker->moment_mode = lya_forest_is_multipole_method(cmd->searchMethod);
    worker->histogram_plan_bytes = histogram_plan_bytes;
    worker->scratch_workers = (size_t)MAX(1, cmd->numthreads);
#ifndef __FAST_MATH__
    worker->polar_lookup=sizeof(REAL)==sizeof(double) && sizeof(real)==sizeof(double)
        && cmd->lya3ThetaBins<=64;
#endif
    if(compute_3pcf && worker->polar_lookup)
        for(int j=1;j<cmd->lya3ThetaBins;j++)
            worker->polar_edges[j]=cos((REAL)LYA_PI*j/cmd->lya3ThetaBins);

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
    free(worker->moments);free(worker->moment_bins);free(worker->moment_seen);
    free(worker->segments);
    free(worker->segment_roots);
    free(worker->neighbors);
    memset(worker, 0, sizeof(*worker));
}

#include "lya_neighbor_geometry.h"

local int lya_append_neighbor(struct cmdline_data *cmd, bodyptr pivot,
                              lya_worker_hist *worker, bodyptr neighbor,
                              REAL distance, ErrorMsg error_message)
{
    lya_neighbor entry;
    REAL cos_theta;
    int radial_bin, theta_bin;
    SUBV(entry.displacement,Pos(neighbor),Pos(pivot));
    entry.radius = distance;
    radial_bin = lya_bin_positive(entry.radius, cmd->lya3RMax, cmd->lya3RBins);
    if (radial_bin < 0 || entry.radius <= 0.0) return SUCCESS;
    DOTVP(cos_theta, entry.displacement, LyaLOS(pivot));
    cos_theta=lya_clamp(cos_theta/entry.radius,-1.,1.);
    theta_bin = worker->polar_lookup
        ? lya_polar_bin(cos_theta,cmd->lya3ThetaBins,worker->polar_edges)
        : lya_bin_theta(racos(cos_theta),cmd->lya3ThetaBins);
    entry.leg_bin=(size_t)radial_bin*cmd->lya3ThetaBins+theta_bin;
    entry.ordinal=worker->neighbor_count;
    entry.forest_id = LyaForestId(neighbor);
    entry.weight = Weight(neighbor);
    entry.weighted_delta = Weight(neighbor) * Kappa(neighbor);
    /* Hard-bin offsets are only used by exact kernels. The multipole grid
     * was checked with L+1 instead of MuBins, so do not form unused MuBins
     * products here. leg_bin supplies the common sorting/grouping key. */
    if (worker->moment_mode) entry.first_index=entry.second_index=0;
    else {
    /* bins3 was checked before worker allocation, so these partial indices
     * and their eventual sum fit size_t. */
    entry.first_index = (((size_t)radial_bin * (size_t)cmd->lya3RBins
                         * (size_t)cmd->lya3ThetaBins + (size_t)theta_bin)
                        * (size_t)cmd->lya3ThetaBins) * (size_t)cmd->lya3MuBins;
    entry.second_index = ((size_t)radial_bin * (size_t)cmd->lya3ThetaBins
                          * (size_t)cmd->lya3ThetaBins + (size_t)theta_bin)
                         * (size_t)cmd->lya3MuBins;
    }
    if (worker->neighbor_count == worker->neighbor_capacity) {
        size_t new_capacity = worker->neighbor_capacity == 0
                            ? 128 : worker->neighbor_capacity * 2;
        size_t bytes, all_scratch, all_bytes, old_segments;
        lya_neighbor *resized;
        if (new_capacity < worker->neighbor_capacity
            || !cballs_size_mul(new_capacity, sizeof(*resized), &bytes)
            || !cballs_size_mul(worker->segment_capacity,2*sizeof(lya_segment)+sizeof(size_t),&old_segments)
            || !cballs_size_mul(bytes,2,&all_scratch)
            || !cballs_size_add(all_scratch,old_segments,&all_scratch)
            || !cballs_size_mul(all_scratch, worker->scratch_workers, &all_scratch)
            || !cballs_size_add(worker->histogram_plan_bytes, all_scratch, &all_bytes)) {
            snprintf(error_message, _ERRORMSGSIZE_,
                     "Ly-alpha neighbor-list size overflow");
            return FAILURE;
        }
        if (cballs_memory_preflight(all_bytes, "Ly-alpha histograms and neighbor scratch",
                                    error_message, _ERRORMSGSIZE_) == FAILURE)
            return FAILURE;
        resized = realloc(worker->neighbors, bytes);
        if (resized == NULL) {
            snprintf(error_message, _ERRORMSGSIZE_,
                     "not enough memory growing Ly-alpha neighbor list");
            return FAILURE;
        }
        worker->neighbors = resized;
        worker->neighbor_capacity = new_capacity;
    }
    worker->neighbors[worker->neighbor_count++] = entry;
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

    if (LyaForestId(q) == LyaForestId(p)) return;
    REAL pivot_weight=lya_pivot_weight(worker,p),pivot_field=lya_pivot_field(worker,p);
    INTEGER multiplicity=lya_pivot_multiplicity(worker);
    if(worker->pivot_multiplicity) {
        if(Id(q)<=worker->pivot_idlo) return;
        if(Id(q)<=worker->pivot_idhi) {
            /* Preserve original row-ID ownership even when a smoothed group
             * straddles q in a shuffled catalog. Do not compare only the
             * representative's ID or multiply by an ineligible member. */
            long double weight=0,field=0; multiplicity=0;
            for(INTEGER j=0;j<worker->pivot_multiplicity;j++) {
                bodyptr member=worker->pivot_members[j];
                if(Id(member)<Id(q)) {
                    weight+=Weight(member);field+=(REAL)(Weight(member)*Kappa(member));multiplicity++;
                }
            }
            if(!multiplicity) return;
            pivot_weight=(REAL)weight;pivot_field=(REAL)field;
        }
    } else if(Id(q)<=Id(p)) return;

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
    denominator = pivot_weight * Weight(q);
    numerator = pivot_field * (Weight(q) * Kappa(q));
    worker->pair_count+=multiplicity;
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
    /* Pair ownership and forest exclusion do not depend on geometry. Do not
     * spend a distance calculation on the reverse orientation. Smoothed
     * groups retain their original per-member ownership check. */
    if(Type(q)!=CELL && !compute_3pcf && !worker->pivot_multiplicity
        && (Id(q)<=Id(p) || LyaForestId((bodyptr)q)==LyaForestId(p))) return SUCCESS;
    DOTPSUBV(distance2, displacement, Pos(p), Pos(q));
    distance = rsqrt(distance2);

    if (Type(q) == CELL) {
        /* Radius is theta-scaled. The cell diagonal bounds both geometric
         * and center-of-mass centers without changing the exact cutoff. */
        const REAL enclosing = Size(q) * rsqrt((REAL)NDIM);
        if (distance >= cutoff + enclosing) return SUCCESS;
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
        return lya_append_neighbor(cmd, p, worker, (bodyptr)q, distance, error_message);
    return SUCCESS;
}

local void lya_accumulate_3pcf(struct cmdline_data *cmd, bodyptr p,
                               lya_worker_hist *worker)
{
    size_t iq;
    size_t ir;
    const REAL pivot_weight = lya_pivot_weight(worker,p);
    const REAL pivot_weighted_delta = lya_pivot_field(worker,p);
    for (iq = 0; iq < worker->neighbor_count; iq++) {
        const lya_neighbor *q = &worker->neighbors[iq];
        const REAL pq_weight = pivot_weight * q->weight;
        const REAL pq_weighted_delta = pivot_weighted_delta * q->weighted_delta;

        for (ir = iq + 1; ir < worker->neighbor_count; ir++) {
            const lya_neighbor *r = &worker->neighbors[ir];
            REAL dot_qr;
            REAL mu;
            REAL numerator;
            REAL denominator;
            int bm;
            size_t forward;
            size_t reverse;

            if (q->forest_id == r->forest_id) continue;
            DOTVP(dot_qr, q->displacement, r->displacement);
            mu = lya_clamp(dot_qr / (q->radius * r->radius), -1.0, 1.0);
            bm = lya_bin_mu(mu, cmd->lya3MuBins);

            denominator = pq_weight * r->weight;
            numerator = pq_weighted_delta * r->weighted_delta;
            worker->ordered_triplet_count += 2;
            if (denominator <= 0.0) continue;
            forward = q->first_index + r->second_index + (size_t)bm;
            reverse = r->first_index + q->second_index + (size_t)bm;
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

#include "lya_triplet_exact.h"
#include "lya_triplet_multipole.h"

typedef struct {
    struct cmdline_data *cmd;
    bodyptr pivot;
    lya_worker_hist *worker;
    int compute_2pcf, compute_3pcf;
} lya_los_accumulator;

local int lya_los_accumulate(bodyptr q, REAL distance, void *context,
                              ErrorMsg error_message)
{
    lya_los_accumulator *acc = context;
    acc->worker->accepted_visits++;
    if (acc->compute_2pcf)
        lya_accumulate_2pcf(acc->cmd, acc->pivot, q, acc->worker);
    if (acc->compute_3pcf && distance < acc->cmd->lya3RMax)
        return lya_append_neighbor(acc->cmd, acc->pivot, acc->worker, q, distance, error_message);
    return SUCCESS;
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
    fprintf(stream,"# pivot_frontier level=%d radius=%.17g max_pixels=%d; pivot_geometry_approximate=%d; neighbors=original_pixels\n",
        cmd->lyaScanLevel,(double)cmd->lyaPivotRadius,cmd->lyaPivotMax,cmd->lyaPivotRadius>0 && cmd->lya2Kernel==0);
    fprintf(stream, "# pair_geometry kernel=%d rp_slop=%.17g rt_slop=%.17g approximate=%d\n",
        cmd->lya2Kernel,(double)cmd->lya2RpSlop,(double)cmd->lya2RtSlop,
        cmd->lya2RpSlop>0 || cmd->lya2RtSlop>0);
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
    fprintf(stream,"# pivot_frontier level=%d radius=%.17g max_pixels=%d; pivot_geometry_approximate=%d; neighbors=original_pixels\n",
        cmd->lyaScanLevel,(double)cmd->lyaPivotRadius,cmd->lyaPivotMax,cmd->lyaPivotRadius>0);
    fprintf(stream, "# geometry_slop mu=%.17g radial=%.17g polar=%.17g; kernel=%d; approximate=%d\n",
        (double)cmd->lya3MuSlop,(double)cmd->lya3RadialSlop,(double)cmd->lya3PolarSlop,cmd->lya3Kernel,
        cmd->lya3MuSlop>0 || cmd->lya3RadialSlop>0 || cmd->lya3PolarSlop>0);
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

// Persistent-node kernels share the estimator and publication helpers above.
#include "lya_forest_cells.h"
#include "lya_pair_cells.h"
#include "lya_triplet_cells.h"
#include "lya_pivot_frontier.h"

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
    const int multipole = lya_forest_is_multipole_method(cmd->searchMethod);
    const int persistent = compute_3pcf && !multipole && cmd->lya3Kernel>=3;
    const int pair_cells = compute_2pcf && cmd->lya2Kernel==1;
    const int legacy2 = compute_2pcf && !pair_cells;
    const int use_los_tree = !persistent && (compute_3pcf || legacy2) && lya_forest_is_los_tree_method(cmd->searchMethod);
    lya_los_index *los_index = NULL;
    lya_pivot_frontier frontier={0};
    const int use_frontier=cmd->lyaScanLevel || cmd->lyaPivotRadius>0;
    lya_los_workspace los_totals = {0};
    double los_build_cpu = 0.0;
    unsigned long long aggregated_pairs=0, direct_pairs=0, segment_accepts=0, approximate_pairs=0;

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
            || lya_size_mul(bins3, (size_t)(multipole ? cmd->lya3LMax+1 : cmd->lya3MuBins), &bins3) == FAILURE)
            goto size_error;
    }

    size_t grid_cells, grid_bytes, worker_bytes, all_bytes, base_bytes;
    if (cballs_resource_base(cmd,gd,&base_bytes)==FAILURE) goto setup_done;
    if (!cballs_size_add(bins2,bins3,&grid_cells)
        || !cballs_size_mul(grid_cells,2*sizeof(real),&grid_bytes)
        || !cballs_size_mul(grid_cells,2*sizeof(real)+sizeof(size_t),&worker_bytes)
        || !cballs_size_mul(worker_bytes,(size_t)MAX(1,cmd->numthreads),&worker_bytes)
        || !cballs_size_add(grid_bytes,worker_bytes,&all_bytes)
        || !cballs_size_add(all_bytes,base_bytes,&all_bytes)) goto size_error;
    if (cballs_memory_preflight(all_bytes,"Ly-alpha global and worker histogram plan",
                               cmd->error_message,sizeof(cmd->error_message)) == FAILURE)
        goto setup_done;

    if (compute_2pcf) {
        if (cballs_calloc_checked((void **)&num2, bins2, sizeof(*num2),
                                  "global Ly-alpha 2PCF numerator",
                                  cmd->error_message,
                                  sizeof(cmd->error_message)) == FAILURE
            || cballs_calloc_checked((void **)&den2, bins2, sizeof(*den2),
                                     "global Ly-alpha 2PCF denominator",
                                     cmd->error_message,
                                     sizeof(cmd->error_message)) == FAILURE)
            goto setup_done;
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
            goto setup_done;
    }
    cutoff = MAX(cutoff2, compute_3pcf ? cmd->lya3RMax : 0.0);

    if (use_los_tree && !pair_cells) {
        double start = CPUTIME;
        if (lya_los_build(&los_index, btable[cat], nbody[cat],
                          (nodeptr)roottable[cat], cmd->error_message) == FAILURE)
            goto setup_done;
        los_build_cpu = CPUTIME - start;
    }

    status = SUCCESS;
setup_done:
    status = lya_parallel_consensus(cmd, status, "Ly-alpha 3D setup");
    if (status == FAILURE) goto cleanup;

    ThreadCount(cmd, gd, nbody[cat], cat);
    const INTEGER first_task = (INTEGER)lya_parallel_first(cmd);
    const INTEGER task_stride = (INTEGER)lya_parallel_stride(cmd);
    /* Fixed grouping, independent of thread count. Small runs need enough
     * blocks to distribute dense pivots; large runs amortize ordered commits.
     * The explicit override supports workload calibration without rebuilding.
     * Original MPI rank ownership continues to use individual pivots. */
    const INTEGER pivot_count = ipmax[cat] - (ipmin - 1);
    const INTEGER chosen_block = cmd->lya3PivotBlock ? cmd->lya3PivotBlock
                                : (pivot_count <= 16384 ? 8 : 64);
    const INTEGER pivot_block = use_los_tree || !lya_forest_is_mpi_method(cmd->searchMethod)
                              ? (compute_3pcf ? chosen_block : 64) : 1;
    INTEGER blocks = pivot_count / pivot_block
                         + (pivot_count % pivot_block != 0);
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
        "\n%s: %s; 2PCF=%d 3PCF=%d cutoff=%g pivot_block=%ld\n",
        cmd->searchMethod, multipole ? "anisotropic moments; approximate mu reconstruction" : (cmd->lyaPivotRadius>0 ? "approximate forest-local pivot geometry; original neighbors" : ((cmd->lya3MuSlop>0 || cmd->lya3RadialSlop>0 || cmd->lya3PolarSlop>0 || cmd->lya2RpSlop>0 || cmd->lya2RtSlop>0) ? "approximate geometry; exact forest exclusions/cutoffs" : "exact Ly-alpha estimator")), compute_2pcf, compute_3pcf, cutoff, (long)pivot_block);

    if (pair_cells && !persistent) {
        lya_cells pair_tree;
        double start=CPUTIME;
        allocation_failed=lya_cells_build(&pair_tree,cmd,btable[cat],nbody[cat],
            ipmin-1,ipmax[cat],all_bytes,cmd->error_message)==FAILURE;
        if(!allocation_failed) {
            verb_print_normal_info(cmd->verbose,cmd->verbose_log,gd->outlog,
                "Ly-alpha 2PCF forest build: nodes=%zu build_CPU=%g shared_with_3pcf=0\n",
                pair_tree.nodes_count,CPUTIME-start);
            allocation_failed=lya_pairs_run(&pair_tree,cmd,gd,btable[cat],ipmin-1,ipmax[cat],
                bins2,pair_tree.plan,num2,den2,&pair_count,&accepted_visits)==FAILURE;
            lya_cells_free(&pair_tree);
        }
        if(allocation_failed || !compute_3pcf) goto workers_done;
        /* The 3PCF discovery now uses only its own sphere. */
        cutoff=cmd->lya3RMax;
        /* The pair hierarchy is already freed. Avoid overlapping it with the
         * independent legacy LOS index in combined kernels 0/1/2. */
        if(use_los_tree) {
            double start=CPUTIME;
            allocation_failed=lya_los_build(&los_index,btable[cat],nbody[cat],
                (nodeptr)roottable[cat],cmd->error_message)==FAILURE;
            los_build_cpu=CPUTIME-start;
            if(allocation_failed) goto workers_done;
        }
    }

    if (persistent) {
        allocation_failed=lya_cells_run(cmd,gd,btable[cat],nbody[cat],ipmin-1,ipmax[cat],
            (nodeptr)roottable[cat],compute_2pcf,bins2,bins3,all_bytes,num2,den2,num3,den3,
            &accepted_visits,&pair_count,&ordered_triplet_count,
            &aggregated_pairs,&direct_pairs,&segment_accepts,&approximate_pairs)==FAILURE;
        goto workers_done;
    }

    if(use_frontier) {
        double start=CPUTIME;
        if(lya_pivot_frontier_build(&frontier,cmd,btable[cat],nbody[cat],ipmin-1,ipmax[cat],
            (nodeptr)roottable[cat],(size_t)pivot_block,compute_3pcf,all_bytes,cmd->error_message)==FAILURE) {
            allocation_failed=TRUE;goto workers_done;
        }
        all_bytes=frontier.plan;blocks=(INTEGER)frontier.blocks;
        verb_print_normal_info(cmd->verbose,cmd->verbose_log,gd->outlog,
            "Ly-alpha pivot frontier: level=%d active=%zu representatives=%zu tasks=%zu radius=%g max_actual_radius=%g build_CPU=%g\n",
            cmd->lyaScanLevel,frontier.active,frontier.count,frontier.blocks,
            (double)cmd->lyaPivotRadius,(double)frontier.max_radius,CPUTIME-start);
    }

#pragma omp parallel shared(allocation_failed,worker_error,num2,den2,num3,den3,accepted_visits,pair_count,ordered_triplet_count)
    {
        lya_worker_hist worker;
        lya_los_workspace los_workspace = {0};
        ErrorMsg local_error = "";
        int worker_ready = lya_worker_init(cmd, &worker, bins2, bins3,
                                           legacy2, compute_3pcf,
                                           all_bytes, local_error) == SUCCESS;
        if (worker_ready && use_los_tree
            && lya_los_workspace_init(los_index, &los_workspace,
                                      local_error) == FAILURE) {
            lya_worker_free(&worker);
            worker_ready = FALSE;
        }
        int worker_failed = !worker_ready;
        INTEGER block;
        if (!worker_ready) {
#pragma omp critical(lya_failure)
            {
                if (!allocation_failed)
                    snprintf(worker_error, sizeof(worker_error), "%s", local_error);
                allocation_failed = TRUE;
            }
        }

#pragma omp barrier
#pragma omp for schedule(dynamic,1) ordered
        for (block = first_task; block < blocks; block += task_stride) {
            INTEGER first = use_frontier?(INTEGER)frontier.offsets[block]:ipmin-1+block*pivot_block;
            INTEGER end = use_frontier?(INTEGER)frontier.offsets[block+1]:first+MIN(pivot_block,ipmax[cat]-first);
            INTEGER pivot;
            worker.accepted_visits = 0;
            worker.pair_count = 0;
            worker.ordered_triplet_count = 0;
            for (pivot = first; pivot < end; pivot++) {
                const lya_pivot_group *group=frontier.groups?&frontier.groups[frontier.order[pivot]]:NULL;
                bodyptr p=group?group->pivot:btable[cat]+(use_frontier?(INTEGER)frontier.order[pivot]:pivot);
                worker.pivot_multiplicity=(group && cmd->lyaPivotRadius>0)?group->count:0;
                if(worker.pivot_multiplicity) {
                    worker.pivot_weight=group->weight;worker.pivot_field=group->field;
                    worker.pivot_members=group->members;worker.pivot_idlo=group->idlo;worker.pivot_idhi=group->idhi;
                }
                worker.neighbor_count = 0;
                if (!worker_failed && Update(p) != FALSE
                    && Mask(p) == MASK_NODE_VALID) {
                    lya_los_accumulator acc = {cmd, p, &worker,
                                               legacy2, compute_3pcf};
                    int walk_status = use_los_tree
                        ? lya_los_query(los_index, &los_workspace, p, cutoff,
                                         !compute_3pcf && !worker.pivot_multiplicity?Id(p):0,
                                         lya_los_accumulate, &acc, local_error)
                        : lya_walktree(cmd, p, (nodeptr)roottable[cat], cutoff,
                                        legacy2, compute_3pcf, &worker, local_error);
                    if (walk_status == FAILURE) {
                        worker_failed = TRUE;
#pragma omp critical(lya_failure)
                        {
                            if (!allocation_failed)
                                snprintf(worker_error, sizeof(worker_error), "%s",
                                         local_error);
                            allocation_failed = TRUE;
                        }
                    } else if (compute_3pcf) {
                        int triplet_status=SUCCESS;
                        INTEGER before_triplets=worker.ordered_triplet_count;
                        if (multipole) triplet_status=lya_accumulate_multipoles(cmd,p,&worker,local_error);
                        else if (cmd->lya3Kernel==1) lya_accumulate_3pcf(cmd,p,&worker);
                        else triplet_status=lya_accumulate_segments(cmd,p,&worker,local_error);
                        /* Scale represented counts once per pivot, outside the
                         * hot neighbor-pair loop. The frontier preflight bounds
                         * the full multiplicity-weighted count before execution. */
                        if(worker.pivot_multiplicity>1)
                            worker.ordered_triplet_count +=
                                (worker.ordered_triplet_count-before_triplets)*(worker.pivot_multiplicity-1);
                        if (triplet_status==FAILURE) {
                            worker_failed=TRUE;
#pragma omp critical(lya_failure)
                            {
                                if (!allocation_failed) snprintf(worker_error,sizeof(worker_error),"%s",local_error);
                                allocation_failed=TRUE;
                            }
                        }
                    }
                }
            }

#pragma omp ordered
            {
                if (worker_ready && !worker_failed) {
                    lya_commit_worker(&worker, num2, den2, num3, den3,
                                      legacy2, compute_3pcf);
                    accepted_visits += worker.accepted_visits;
                    pair_count += worker.pair_count;
                    ordered_triplet_count += worker.ordered_triplet_count;
                }
            }
        }
#pragma omp critical(lya_segment_counters)
        { aggregated_pairs+=worker.aggregated_pairs; direct_pairs+=worker.direct_pairs; segment_accepts+=worker.segment_accepts; approximate_pairs+=worker.approximate_pairs; }
        if (worker_ready) lya_worker_free(&worker);
        if (use_los_tree) {
#pragma omp critical(lya_los_counters)
            {
                los_totals.octree_nodes += los_workspace.octree_nodes;
                los_totals.forest_skips += los_workspace.forest_skips;
                los_totals.forest_hits += los_workspace.forest_hits;
                los_totals.radial_nodes += los_workspace.radial_nodes;
                los_totals.pixel_tests += los_workspace.pixel_tests;
            }
            lya_los_workspace_free(&los_workspace);
        }
    }

workers_done:
    if (allocation_failed) {
        if (cmd->error_message[0] == '\0')
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     "Ly-alpha worker failed: %s",
                     worker_error[0] != '\0' ? worker_error : "unknown error");
    }

    status = lya_parallel_consensus(cmd, allocation_failed ? FAILURE : SUCCESS,
                                    "Ly-alpha 3D workers");
    if (status == FAILURE) goto cleanup;
    INTEGER counters[3] = {accepted_visits, pair_count, ordered_triplet_count};
    if (lya_parallel_reduce_reals(cmd, num2, bins2) == FAILURE
        || lya_parallel_reduce_reals(cmd, den2, bins2) == FAILURE
        || lya_parallel_reduce_reals(cmd, num3, bins3) == FAILURE
        || lya_parallel_reduce_reals(cmd, den3, bins3) == FAILURE
        || lya_parallel_reduce_integers(cmd, counters, 3) == FAILURE) {
        status = FAILURE;
        goto cleanup;
    }
    if (!lya_parallel_publish(cmd)) {
        status = SUCCESS;
        goto publication;
    }
    accepted_visits = counters[0];
    pair_count = counters[1];
    ordered_triplet_count = counters[2];
    status = FAILURE;

    gd->nbbcalc += accepted_visits;
    if (!cballs_opt_no_out_hist(cmd)) {
        if (compute_2pcf
            && lya_write_2pcf(cmd, gd, num2, den2, pair_count) == FAILURE)
            goto publication;
        if (compute_3pcf
            && (multipole ? lya_write_multipoles(cmd,gd,num3,den3,ordered_triplet_count)
                : lya_write_3pcf(cmd, gd, num3, den3,ordered_triplet_count)) == FAILURE)
            goto publication;
    }

    if (compute_3pcf && !multipole && cmd->lya3Kernel!=1) verb_print_normal_info(cmd->verbose,cmd->verbose_log,gd->outlog,
        "Ly-alpha 3PCF kernel=%d aggregated_pairs=%llu direct_pairs=%llu segment_accepts=%llu approximate_pairs=%llu (rank-local)\n",
        cmd->lya3Kernel,aggregated_pairs,direct_pairs,segment_accepts,approximate_pairs);
    gd->cpusearch = CPUTIME - cpustart;
    if (use_los_tree)
        verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
            "LOS-tree: forests=%zu build_CPU=%g octree_nodes=%llu "
            "forest_skips=%llu forest_hits=%llu radial_nodes=%llu pixel_tests=%llu\n",
            lya_los_forest_count(los_index), los_build_cpu,
            los_totals.octree_nodes, los_totals.forest_skips,
            los_totals.forest_hits, los_totals.radial_nodes, los_totals.pixel_tests);
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
        "%s: accepted=%" INTEGER_FMT " pairs=%" INTEGER_FMT
        " ordered_triplets=%" INTEGER_FMT " CPU=%g\n",
        cmd->searchMethod, accepted_visits, pair_count,
        ordered_triplet_count, gd->cpusearch);
    status = SUCCESS;
    goto publication;

size_error:
    snprintf(cmd->error_message, _ERRORMSGSIZE_,
             "Ly-alpha histogram dimensions overflow size_t");

    goto setup_done;

publication:
    status = lya_parallel_consensus(cmd, status, "Ly-alpha 3D output");
cleanup:
    if(status==SUCCESS && lya_parallel_publish(cmd)) {
        size_t pair_shape[2]={(size_t)cmd->lya2RpBins,(size_t)cmd->lya2RtBins};
        size_t triple_shape[5]={(size_t)cmd->lya3RBins,(size_t)cmd->lya3RBins,(size_t)cmd->lya3ThetaBins,(size_t)cmd->lya3ThetaBins,(size_t)(multipole?cmd->lya3LMax+1:cmd->lya3MuBins)};
        cballs_result_adopt("pair_numerator",(void **)&num2,0,2,pair_shape);
        cballs_result_adopt("pair_denominator",(void **)&den2,0,2,pair_shape);
        cballs_result_adopt(multipole?"moments_numerator":"triple_numerator",(void **)&num3,0,5,triple_shape);
        cballs_result_adopt(multipole?"moments_denominator":"triple_denominator",(void **)&den3,0,5,triple_shape);
    }
    lya_pivot_frontier_free(&frontier);
    lya_los_free(los_index);
    gd->cpusearch = CPUTIME - cpustart;
    free(num2);
    free(den2);
    free(num3);
    free(den3);
    return status;
}

#undef LYA_PI
