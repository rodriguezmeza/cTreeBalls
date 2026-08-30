/* =============================================================================
 MODULE: search_octree_3pcf_3d_omp.c      [cTreeBalls]
 Purpose: 3D scalar-field 3PCF multipoles using the cTreeBalls octree traversal.

 Estimator:
   For every primary i, construct shell-binned spherical harmonic coefficients

       a_lm^i(b) = sum_{j in shell b} w_j delta_j Y_lm^*(rhat_ij)

   and accumulate

       N_l(b1,b2) += w_i delta_i (2l+1) sum_{j != k} w_j delta_j w_k delta_k
                                      P_l(rhat_ij . rhat_ik)

   implemented as

       4*pi * w_i delta_i * [ sum_m a_lm(b1) a_lm^*(b2)
                              - delta_b1b2 (2l+1)/(4*pi) sum_j (w_j delta_j)^2 ]

   The denominator is the analogous weighted pair count without delta fields.

   The same neighbor traversal can also accumulate the scalar 2PCF

       xi(b) = sum_(i,j in b) w_i delta_i w_j delta_j
             / sum_(i,j in b) w_i w_j.

   With no explicit mode option the addon computes the 3PCF only, preserving
   its original behavior. Use compute-2pcf-3d and compute-3pcf-3d together,
   or use only-2pcf-3d/only-3pcf-3d.

 Notes:
   - This is an exact-in-angle leaf traversal by default. The octree is used for
     pruning outside Rcut; cells are not collapsed into approximate angular nodes.
   - Harmonics use the standard orthonormal Condon--Shortley convention and the
     same m>=0 triangular layout used by ENCORE: n = ell*(ell+1)/2 + m.
   - The code supports ASCII/FITS catalog readers already present in cBalls;
     Kappa(p) is interpreted as delta, Weight(p) as w.
   - exclude-same-los (aliases exclude-los and exclude-pivot-los) removes
     pivot-neighbor pairs with equal LOS IDs. It therefore removes i=j_LOS and
     i=k_LOS terms, but does not subtract j=k_LOS terms from the harmonic
     product.
 =============================================================================*/

#include "globaldefs.h"

#ifndef OCTREE3PCF3D_MAX_LMAX
#define OCTREE3PCF3D_MAX_LMAX 32
#endif

#define CB3D_PI      3.14159265358979323846264338327950288419716939937510
#define CB3D_FOURPI  12.566370614359172953850573533118011536788677597500

#define CB3D_LM_INDEX(ell,m) ((ell)*((ell)+1)/2 + (m))
#define CB3D_ACC_INDEX(ell,b1,b2,nb) ((((ell)*(nb) + (b1))*(nb)) + (b2))
#define CB3D_ALM_INDEX(b,nlm,n) (((b)*(nlm)) + (n))

typedef struct {
    int nbins;
    int lmax;
    int nlm;
    REAL *alm_re;
    REAL *alm_im;
    REAL *ylm_norm;
    REAL *ylm_a;
    REAL *ylm_b;
    REAL *ylm_diag;
    REAL *shell_w;
    REAL *shell_w2;
    REAL *shell_f2;
    REAL *num;
    REAL *den;
    REAL *xi_num;
    REAL *xi_den;
    int *active_bins;
    unsigned char *bin_active;
    int nactive;
    int compute_xi;
    int compute_zeta;
    int exclude_same_los;
    int read_mask;
    int use_log_hist;
    int use_periodic;
    REAL rmin;
    REAL rmin2;
    REAL rcut;
    REAL rcut2;
    REAL i_delta_r;
    REAL log_bins_per_decade;
    REAL log_range_n;
    REAL cell_radius_scale;
    REAL pivot_weight;
    REAL pivot_field;
    INTEGER pivot_los;
    INTEGER nbbcalcthread;
    INTEGER npivots_thread;
} cb3d_hist, *cb3d_histptr;

local int cb3d_init_hist(struct cmdline_data* cmd, struct global_data* gd,
                         cb3d_histptr h,
                         int compute_xi, int compute_zeta,
                         int exclude_same_los, ErrorMsg allocation_error);
local void cb3d_free_hist(cb3d_histptr h);
local void cb3d_clear_pivot(cb3d_histptr h);
local void cb3d_activate_bin(cb3d_histptr h, int b);
local int cb3d_bin_index(cb3d_histptr h, REAL r);
local int cb3d_reject_cell(struct global_data* gd, bodyptr p, nodeptr q,
                           cb3d_histptr h);
local void cb3d_walktree(struct global_data* gd,
                         bodyptr p, nodeptr q, cb3d_histptr h);
local void cb3d_add_body(struct global_data* gd,
                         bodyptr p, bodyptr q, cb3d_histptr h);
local void cb3d_accumulate_pivot(cb3d_histptr h);
local int cb3d_get_modes(struct cmdline_data* cmd,
                         int *compute_xi, int *compute_zeta);
local int cb3d_exclude_same_los_option(struct cmdline_data* cmd);
local void cb3d_init_ylm_norm(int lmax, REAL *norm);
local void cb3d_init_ylm_recurrence(cb3d_histptr h);
local void cb3d_accumulate_ylm_cartesian(cb3d_histptr h, int b, REAL field,
                                         double xhat, double yhat,
                                         double zhat);
local int cb3d_print_zeta(struct cmdline_data* cmd, struct global_data* gd,
                          REAL *num, REAL *den, INTEGER npivots, INTEGER nbb);
local int cb3d_print_xi(struct cmdline_data* cmd, struct global_data* gd,
                        REAL *num, REAL *den, INTEGER npivots, INTEGER nbb);
local REAL cb3d_bin_center(struct cmdline_data* cmd, struct global_data* gd,
                           int bzero);
local void cb3d_add_arrays(REAL *dst, const REAL *src, size_t n);
local void cb3d_zero_array(REAL *dst, size_t n);

/* ------------------------------------------------------------------------- */

global int searchcalc_octree_3pcf_3d_omp(struct cmdline_data* cmd,
                                         struct global_data* gd,
                                         bodyptr *btable, INTEGER *nbody,
                                         INTEGER ipmin, INTEGER *ipmax,
                                         int cat1, int cat2)
{
    string routineName = "searchcalc_octree_3pcf_3d_omp";
    double cpustart = CPUTIME;
    int nb = cmd->sizeHistN;
    int lmax = cmd->mChebyshev;
    int nlm;
    size_t nacc;
    REAL *num_global = NULL;
    REAL *den_global = NULL;
    REAL *xi_num_global = NULL;
    REAL *xi_den_global = NULL;
    int compute_xi = FALSE;
    int compute_zeta = FALSE;
    int exclude_same_los = FALSE;
    INTEGER nbb_global = 0;
    INTEGER npivots_global = 0;
    int allocation_failed = FALSE;
    int status = FAILURE;
    ErrorMsg allocation_error = "";

#if NDIM != 3
    snprintf(cmd->error_message, _ERRORMSGSIZE_,
             "%s: this estimator requires NDIM=3", routineName);
    return FAILURE;
#endif

    if (cb3d_get_modes(cmd, &compute_xi, &compute_zeta) == FAILURE)
        return FAILURE;
    exclude_same_los = cb3d_exclude_same_los_option(cmd);

    if (nb <= 0) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: sizeHistN must be positive", routineName);
        return FAILURE;
    }
    if (lmax < 0 || lmax > OCTREE3PCF3D_MAX_LMAX) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: invalid mChebyshev/lmax=%d; supported range is 0..%d",
                 routineName, lmax, OCTREE3PCF3D_MAX_LMAX);
        return FAILURE;
    }
    nlm = (lmax + 1)*(lmax + 2)/2;
    nacc = (size_t)(lmax + 1) * (size_t)nb * (size_t)nb;

    if (cmd->useLogHist && cmd->rminHist <= 0.0) {
        verb_print_warning(cmd->verbose,
            "%s: log histograms with rminHist<=0 use cBalls legacy bin centers; consider rminHist>0 for 3D.",
            routineName);
    }

    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
        "\n%s: 3D scalar estimator, nbins=%d, lmax=%d, nlm=%d, 2pcf=%d, 3pcf=%d, exclude_same_los=%d\n",
        routineName, nb, lmax, nlm, compute_xi, compute_zeta,
        exclude_same_los);
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
        "%s: input fields: Pos=(x,y,z), delta=Kappa, weight=Weight\n",
        routineName);
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
        "%s: physical cutoff and exact-leaf pruning use rangeN=%g"
        " (global Rcut=%g)\n\n",
        routineName, cmd->rangeN, gd->Rcut);

#ifdef OPENMPCODE
    ThreadCount(cmd, gd, nbody[cat1], cat1);
#else
#error `OPENMPMACHINE` is not defined. Switch it on in Makefile_settings
#endif

    if (compute_zeta) {
        if (cballs_calloc_checked((void **)&num_global, nacc,
                                  sizeof(*num_global), "3D 3PCF numerator",
                                  cmd->error_message,
                                  sizeof(cmd->error_message)) == FAILURE)
            goto cleanup;
        if (cballs_calloc_checked((void **)&den_global, nacc,
                                  sizeof(*den_global), "3D 3PCF denominator",
                                  cmd->error_message,
                                  sizeof(cmd->error_message)) == FAILURE)
            goto cleanup;
    }
    if (compute_xi) {
        if (cballs_calloc_checked((void **)&xi_num_global, (size_t)nb,
                                  sizeof(*xi_num_global), "3D 2PCF numerator",
                                  cmd->error_message,
                                  sizeof(cmd->error_message)) == FAILURE)
            goto cleanup;
        if (cballs_calloc_checked((void **)&xi_den_global, (size_t)nb,
                                  sizeof(*xi_den_global), "3D 2PCF denominator",
                                  cmd->error_message,
                                  sizeof(cmd->error_message)) == FAILURE)
            goto cleanup;
    }

    verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                        "\nRunning 3D scalar 2PCF/3PCF... completed pivot node:\n");

#pragma omp parallel default(none) \
    shared(cmd,gd,btable,nbody,roottable,ipmin,ipmax,cat1,cat2, \
           num_global,den_global,xi_num_global,xi_den_global,nacc,nb, \
           compute_xi,compute_zeta,exclude_same_los, \
           nbb_global,npivots_global,allocation_failed,allocation_error)
    {
        bodyptr p;
        INTEGER ip;
        cb3d_hist hist;
        ErrorMsg thread_error = "";
        int hist_ready =
            cb3d_init_hist(cmd, gd, &hist, compute_xi, compute_zeta,
                           exclude_same_los, thread_error) == SUCCESS;

        if (!hist_ready) {
#pragma omp critical(cb3d_allocation_failure)
            {
                if (!allocation_failed)
                    snprintf(allocation_error, sizeof(allocation_error),
                             "%s", thread_error);
                allocation_failed = TRUE;
            }
        }

#pragma omp barrier

#pragma omp for schedule(dynamic)
        for (p = btable[cat1] + ipmin - 1; p < btable[cat1] + ipmax[cat1]; p++) {
            if (allocation_failed) continue;
            if (Update(p) == FALSE) continue;
            if (hist.read_mask && Mask(p) == FALSE) continue;

            cb3d_clear_pivot(&hist);
            hist.pivot_weight = Weight(p);
            hist.pivot_field = Weight(p) * Kappa(p);
            if (hist.exclude_same_los)
                hist.pivot_los = Octree3pcf3dLosId(p);
            cb3d_walktree(gd, p, ((nodeptr) roottable[cat2]), &hist);
            cb3d_accumulate_pivot(&hist);
            hist.npivots_thread++;

            ip = p - btable[cat1] + 1;
            if (cmd->stepState > 0 && ip % cmd->stepState == 0)
                verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog, "%ld ", ip);
        }

        if (hist_ready) {
            if (!allocation_failed) {
#pragma omp critical(cb3d_histogram_reduction)
                {
                    if (compute_zeta) {
                        cb3d_add_arrays(num_global, hist.num, nacc);
                        cb3d_add_arrays(den_global, hist.den, nacc);
                    }
                    if (compute_xi) {
                        cb3d_add_arrays(xi_num_global, hist.xi_num,
                                        (size_t)nb);
                        cb3d_add_arrays(xi_den_global, hist.xi_den,
                                        (size_t)nb);
                    }
                    nbb_global += hist.nbbcalcthread;
                    npivots_global += hist.npivots_thread;
                }
            }
            cb3d_free_hist(&hist);
        }
    }

    if (allocation_failed) {
        snprintf(cmd->error_message, sizeof(cmd->error_message),
                 "%s: %s", routineName,
                 allocation_error[0] != '\0'
                     ? allocation_error : "worker allocation failed");
        goto cleanup;
    }

    verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog, "\n");
    if (compute_zeta
        && cb3d_print_zeta(cmd, gd, num_global, den_global,
                           npivots_global, nbb_global) == FAILURE)
        goto cleanup;
    if (compute_xi
        && cb3d_print_xi(cmd, gd, xi_num_global, xi_den_global,
                         npivots_global, nbb_global) == FAILURE)
        goto cleanup;

    gd->cpusearch = CPUTIME - cpustart;
    verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                        "\nGoing out: CPU time = %lf %s\n",
                        CPUTIME-cpustart, PRNUNITOFTIMEUSED);

    status = SUCCESS;

cleanup:
    free(num_global);
    free(den_global);
    free(xi_num_global);
    free(xi_den_global);
    return status;
}

/* ------------------------------------------------------------------------- */

local int cb3d_get_modes(struct cmdline_data* cmd,
                         int *compute_xi, int *compute_zeta)
{
    int has_xi = scanopt(cmd->options, "compute-2pcf-3d");
    int has_zeta = scanopt(cmd->options, "compute-3pcf-3d");
    int only_xi = scanopt(cmd->options, "only-2pcf-3d");
    int only_zeta = scanopt(cmd->options, "only-3pcf-3d");

    if (only_xi && only_zeta) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "octree-3pcf-3d-omp: only-2pcf-3d and only-3pcf-3d are incompatible");
        return FAILURE;
    }
    if (only_xi) {
        *compute_xi = TRUE;
        *compute_zeta = FALSE;
    } else if (only_zeta) {
        *compute_xi = FALSE;
        *compute_zeta = TRUE;
    } else if (has_xi || has_zeta) {
        *compute_xi = has_xi;
        *compute_zeta = has_zeta;
    } else {
        *compute_xi = FALSE;
        *compute_zeta = TRUE;
    }
    return SUCCESS;
}

local int cb3d_exclude_same_los_option(struct cmdline_data* cmd)
{
    return scanopt(cmd->options, "exclude-same-los")
        || scanopt(cmd->options, "exclude-los")
        || scanopt(cmd->options, "exclude-pivot-los");
}

local int cb3d_init_hist(struct cmdline_data* cmd, struct global_data* gd,
                         cb3d_histptr h,
                         int compute_xi, int compute_zeta,
                         int exclude_same_los, ErrorMsg allocation_error)
{
    int nb = cmd->sizeHistN;
    int lmax = cmd->mChebyshev;
    int nlm = (lmax + 1)*(lmax + 2)/2;
    size_t nalm = (size_t)nb * (size_t)nlm;
    size_t nacc = (size_t)(lmax + 1) * (size_t)nb * (size_t)nb;

    memset(h, 0, sizeof(*h));
    h->nbins = nb;
    h->lmax = lmax;
    h->nlm = nlm;
    h->compute_xi = compute_xi;
    h->compute_zeta = compute_zeta;
    h->exclude_same_los = exclude_same_los;
    h->read_mask = cballs_opt_read_mask(cmd);
    h->use_log_hist = cmd->useLogHist;
    h->use_periodic = cmd->usePeriodic;
    h->rmin = cmd->rminHist;
    h->rmin2 = cmd->rminHist * cmd->rminHist;
    /* rangeN is the estimator's physical radial domain.  gd->Rcut is a
     * legacy tree-search tuning value and may be rangeN/theta; using it here
     * would silently discard valid bodies before radial binning. */
    h->rcut = cmd->rangeN;
    h->rcut2 = cmd->rangeN * cmd->rangeN;
    h->i_delta_r = gd->i_deltaR;
    h->log_bins_per_decade = cmd->logHistBinsPD;
    h->log_range_n = cmd->rangeN > 0.0 ? rlog10(cmd->rangeN) : 0.0;
    h->cell_radius_scale = cmd->theta > 0.0 ? cmd->theta : 1.0;
#define CB3D_ALLOC(field,count)                                           \
    do {                                                                  \
        if (cballs_malloc_checked((void **)&h->field, (count),            \
                                  sizeof(*h->field), #field,              \
                                  allocation_error, _ERRORMSGSIZE_)       \
            == FAILURE)                                                   \
            goto fail;                                                    \
    } while (0)

    if (compute_zeta) {
        CB3D_ALLOC(alm_re, nalm);
        CB3D_ALLOC(alm_im, nalm);
        CB3D_ALLOC(ylm_norm, (size_t)nlm);
        CB3D_ALLOC(ylm_a, (size_t)nlm);
        CB3D_ALLOC(ylm_b, (size_t)nlm);
        CB3D_ALLOC(ylm_diag, (size_t)lmax + 1);
        CB3D_ALLOC(shell_w, (size_t)nb);
        CB3D_ALLOC(shell_w2, (size_t)nb);
        CB3D_ALLOC(shell_f2, (size_t)nb);
        CB3D_ALLOC(num, nacc);
        CB3D_ALLOC(den, nacc);
        CB3D_ALLOC(active_bins, (size_t)nb);
        CB3D_ALLOC(bin_active, (size_t)nb);
        memset(h->bin_active, 0, (size_t)nb * sizeof(*h->bin_active));
        cb3d_init_ylm_norm(lmax, h->ylm_norm);
        cb3d_init_ylm_recurrence(h);
        cb3d_zero_array(h->num, nacc);
        cb3d_zero_array(h->den, nacc);
    }
    if (compute_xi) {
        CB3D_ALLOC(xi_num, (size_t)nb);
        CB3D_ALLOC(xi_den, (size_t)nb);
        cb3d_zero_array(h->xi_num, (size_t)nb);
        cb3d_zero_array(h->xi_den, (size_t)nb);
    }
    h->nbbcalcthread = 0;
    h->npivots_thread = 0;
    cb3d_clear_pivot(h);
    return SUCCESS;

fail:
    cb3d_free_hist(h);
    return FAILURE;

#undef CB3D_ALLOC
}

local void cb3d_free_hist(cb3d_histptr h)
{
    free(h->alm_re);
    free(h->alm_im);
    free(h->ylm_norm);
    free(h->ylm_a);
    free(h->ylm_b);
    free(h->ylm_diag);
    free(h->shell_w);
    free(h->shell_w2);
    free(h->shell_f2);
    free(h->num);
    free(h->den);
    free(h->xi_num);
    free(h->xi_den);
    free(h->active_bins);
    free(h->bin_active);
}

local void cb3d_clear_pivot(cb3d_histptr h)
{
    int i;
    if (!h->compute_zeta) return;
    for (i = 0; i < h->nactive; i++)
        h->bin_active[h->active_bins[i]] = 0;
    h->nactive = 0;
}

local void cb3d_activate_bin(cb3d_histptr h, int b)
{
    size_t offset;

    if (h->bin_active[b]) return;
    h->bin_active[b] = 1;
    h->active_bins[h->nactive++] = b;
    offset = (size_t)b * (size_t)h->nlm;
    memset(h->alm_re + offset, 0, (size_t)h->nlm * sizeof(*h->alm_re));
    memset(h->alm_im + offset, 0, (size_t)h->nlm * sizeof(*h->alm_im));
    h->shell_w[b] = 0.0;
    h->shell_w2[b] = 0.0;
    h->shell_f2[b] = 0.0;
}

local int cb3d_bin_index(cb3d_histptr h, REAL r)
{
    int n;
#ifndef NORMALHISTSCALE
    if (h->use_log_hist) {
        if (h->rmin == 0)
            n = (int)(h->log_bins_per_decade
                      * (rlog10(r) - h->log_range_n)
                      + h->nbins) + 1;
        else
            n = (int)(rlog10(r/h->rmin) * h->i_delta_r) + 1;
    } else {
        n = (int)((r - h->rmin) * h->i_delta_r) + 1;
    }
#else
    n = (int)((r - h->rmin) * h->i_delta_r) + 1;
#endif
    if (n < 1 || n > h->nbins) return -1;
    return n - 1;
}

/* Radius(q) is divided by theta when cell properties are finalized.  Restore
 * the geometric bound before rejecting a whole cell, or theta > 1 can discard
 * bodies that are still inside Rcut. */
local int cb3d_reject_cell(struct global_data* gd, bodyptr p, nodeptr q,
                           cb3d_histptr h)
{
    REAL distance_squared;
    REAL cutoff = h->rcut + Radius(q) * h->cell_radius_scale;
    compute_vector dr;

    DOTPSUBV(distance_squared, dr, Pos(p), Pos(q));
    if (h->use_periodic) {
        VWrapAll(dr);
        DOTVP(distance_squared, dr, dr);
    }
    return distance_squared >= cutoff * cutoff;
}

local void cb3d_walktree(struct global_data* gd,
                         bodyptr p, nodeptr q, cb3d_histptr h)
{
    nodeptr child;

    if (((nodeptr)p) == q) return;
    if (Type(q) == CELL) {
        if ((h->read_mask && Mask(q) == MASK_NODE_MASKED)
            || cb3d_reject_cell(gd, p, q, h))
            return;
        for (child = More(q); child != Next(q); child = Next(child))
            cb3d_walktree(gd, p, child, h);
    } else {
        cb3d_add_body(gd, p, (bodyptr)q, h);
    }
}

local void cb3d_add_body(struct global_data* gd,
                         bodyptr p, bodyptr q, cb3d_histptr h)
{
    REAL dr2;
    REAL dr1;
    compute_vector dr;
    int b;
    REAL f, w;

    if (p == q) return;
    if (h->read_mask && Mask(q) == FALSE) return;
    if (Update(q) == FALSE) return;
    if (h->exclude_same_los
        && h->pivot_los == Octree3pcf3dLosId(q)) return;

    DOTPSUBV(dr2, dr, Pos(p), Pos(q));
    if (h->use_periodic) {
        VWrapAll(dr);
        DOTVP(dr2, dr, dr);
    }
    if (dr2 <= h->rmin2 || dr2 >= h->rcut2) return;
    dr1 = rsqrt(dr2);
    b = cb3d_bin_index(h, dr1);
    if (b < 0) return;

    w = Weight(q);
    f = w * Kappa(q);

    if (h->compute_xi) {
        h->xi_num[b] += h->pivot_field * f;
        h->xi_den[b] += h->pivot_weight * w;
    }

    if (h->compute_zeta) {
        double invr = 1.0 / (double)dr1;
        double xhat = (double)dr[0] * invr;
        double yhat = (double)dr[1] * invr;
        double zhat = (double)dr[2] * invr;

        cb3d_activate_bin(h, b);
        cb3d_accumulate_ylm_cartesian(h, b, f, xhat, yhat, zhat);
        h->shell_w[b]  += w;
        h->shell_w2[b] += w*w;
        h->shell_f2[b] += f*f;
    }

    h->nbbcalcthread++;
}

local void cb3d_accumulate_pivot(cb3d_histptr h)
{
    int ell, m, i1_active, i2_active;
    REAL wp_delta;
    REAL wp_count;

    if (!h->compute_zeta) return;
    wp_delta = h->pivot_field;
    wp_count = h->pivot_weight;

    for (ell = 0; ell <= h->lmax; ell++) {
        for (i1_active = 0; i1_active < h->nactive; i1_active++) {
            int b1 = h->active_bins[i1_active];
            for (i2_active = i1_active;
                 i2_active < h->nactive; i2_active++) {
                int b2 = h->active_bins[i2_active];
                REAL power = 0.0;
                REAL d;
                REAL diag_f = 0.0;
                REAL diag_w = 0.0;
                size_t idx = CB3D_ACC_INDEX(ell, b1, b2, h->nbins);
                REAL num_increment;
                REAL den_increment;

                for (m = 0; m <= ell; m++) {
                    int n = CB3D_LM_INDEX(ell, m);
                    size_t i1 = CB3D_ALM_INDEX(b1, h->nlm, n);
                    size_t i2 = CB3D_ALM_INDEX(b2, h->nlm, n);
                    REAL reprod = h->alm_re[i1]*h->alm_re[i2]
                                + h->alm_im[i1]*h->alm_im[i2];
                    power += (m == 0) ? reprod : 2.0*reprod;
                }

                d = h->shell_w[b1] * h->shell_w[b2];
                if (b1 == b2) {
                    diag_f = ((REAL)(2*ell + 1) / (REAL)CB3D_FOURPI) * h->shell_f2[b1];
                    diag_w = h->shell_w2[b1];
                }

                num_increment = wp_delta
                              * ((REAL)CB3D_FOURPI * (power - diag_f));
                den_increment = wp_count * (d - diag_w);
                h->num[idx] += num_increment;
                h->den[idx] += den_increment;
                if (b1 != b2) {
                    idx = CB3D_ACC_INDEX(ell, b2, b1, h->nbins);
                    h->num[idx] += num_increment;
                    h->den[idx] += den_increment;
                }
            }
        }
    }
}

local void cb3d_init_ylm_norm(int lmax, REAL *norm)
{
    int ell, m;
    for (ell = 0; ell <= lmax; ell++) {
        for (m = 0; m <= ell; m++) {
            int n = CB3D_LM_INDEX(ell, m);
            double logfac = lgamma((double)(ell - m) + 1.0)
                          - lgamma((double)(ell + m) + 1.0);
            norm[n] = (REAL)sqrt(
                ((2.0*(double)ell + 1.0)/(4.0*CB3D_PI)) * exp(logfac));
        }
    }
}

local void cb3d_init_ylm_recurrence(cb3d_histptr h)
{
    int ell, m;

    cb3d_zero_array(h->ylm_a, (size_t)h->nlm);
    cb3d_zero_array(h->ylm_b, (size_t)h->nlm);
    cb3d_zero_array(h->ylm_diag, (size_t)h->lmax + 1);
    for (m = 0; m < h->lmax; m++) {
        int diagonal = CB3D_LM_INDEX(m, m);
        int upward = CB3D_LM_INDEX(m + 1, m);
        int next_diagonal = CB3D_LM_INDEX(m + 1, m + 1);
        double scale = (double)(2*m + 1);

        h->ylm_a[upward] = (REAL)(scale
            * (double)h->ylm_norm[upward]
            / (double)h->ylm_norm[diagonal]);
        h->ylm_diag[m] = (REAL)(-scale
            * (double)h->ylm_norm[next_diagonal]
            / (double)h->ylm_norm[diagonal]);

        for (ell = m + 2; ell <= h->lmax; ell++) {
            int n = CB3D_LM_INDEX(ell, m);
            int previous = CB3D_LM_INDEX(ell - 1, m);
            int previous2 = CB3D_LM_INDEX(ell - 2, m);
            double denominator = (double)(ell - m);

            h->ylm_a[n] = (REAL)(
                ((2.0*(double)ell - 1.0) / denominator)
                * (double)h->ylm_norm[n]
                / (double)h->ylm_norm[previous]);
            h->ylm_b[n] = (REAL)(
                ((double)(ell + m - 1) / denominator)
                * (double)h->ylm_norm[n]
                / (double)h->ylm_norm[previous2]);
        }
    }
}

/* Cartesian normalized spherical harmonics for m>=0. Accumulate directly into
 * the selected radial shell, as octree-GGG does with its recurrence-backed
 * thread histograms. This removes two maximum-lmax temporary arrays and the
 * second pass over all modes for every accepted neighbor. */
local void cb3d_accumulate_ylm_cartesian(cb3d_histptr h, int radial_bin,
                                         REAL field, double xhat,
                                         double yhat, double zhat)
{
    int m, ell;
    size_t offset = (size_t)radial_bin * (size_t)h->nlm;
    double ymm_re = (double)h->ylm_norm[0];
    double ymm_im = 0.0;
    double weighted_field = (double)field;

    for (m = 0; m <= h->lmax; m++) {
        int n = CB3D_LM_INDEX(m, m);
        size_t index = offset + (size_t)n;
        h->alm_re[index] += (REAL)(weighted_field * ymm_re);
        h->alm_im[index] -= (REAL)(weighted_field * ymm_im);

        if (m < h->lmax) {
            double y_lm2_re = ymm_re;
            double y_lm2_im = ymm_im;
            double y_lm1_re;
            double y_lm1_im;
            double next_ymm_re;
            double next_ymm_im;

            n = CB3D_LM_INDEX(m + 1, m);
            y_lm1_re = (double)h->ylm_a[n] * zhat * ymm_re;
            y_lm1_im = (double)h->ylm_a[n] * zhat * ymm_im;
            index = offset + (size_t)n;
            h->alm_re[index] += (REAL)(weighted_field * y_lm1_re);
            h->alm_im[index] -= (REAL)(weighted_field * y_lm1_im);

            for (ell = m + 2; ell <= h->lmax; ell++) {
                double ylm_re;
                double ylm_im;

                n = CB3D_LM_INDEX(ell, m);
                ylm_re = (double)h->ylm_a[n] * zhat * y_lm1_re
                       - (double)h->ylm_b[n] * y_lm2_re;
                ylm_im = (double)h->ylm_a[n] * zhat * y_lm1_im
                       - (double)h->ylm_b[n] * y_lm2_im;
                index = offset + (size_t)n;
                h->alm_re[index] += (REAL)(weighted_field * ylm_re);
                h->alm_im[index] -= (REAL)(weighted_field * ylm_im);
                y_lm2_re = y_lm1_re;
                y_lm2_im = y_lm1_im;
                y_lm1_re = ylm_re;
                y_lm1_im = ylm_im;
            }

            next_ymm_re = (double)h->ylm_diag[m]
                        * (ymm_re*xhat - ymm_im*yhat);
            next_ymm_im = (double)h->ylm_diag[m]
                        * (ymm_re*yhat + ymm_im*xhat);
            ymm_re = next_ymm_re;
            ymm_im = next_ymm_im;
        }
    }
}

local REAL cb3d_bin_center(struct cmdline_data* cmd, struct global_data* gd,
                           int bzero)
{
    int n = bzero + 1;
    REAL rbin;
#ifndef NORMALHISTSCALE
    if (cmd->useLogHist) {
        REAL rbinlog;
        if (cmd->rminHist == 0)
            rbinlog = ((REAL)(n - cmd->sizeHistN))/cmd->logHistBinsPD
                    + rlog10(cmd->rangeN);
        else
            rbinlog = rlog10(cmd->rminHist) + ((REAL)n - 0.5)*gd->deltaR;
        rbin = rpow(10.0, rbinlog);
    } else {
        rbin = cmd->rminHist + ((REAL)n - 0.5)/gd->i_deltaR;
    }
#else
    rbin = cmd->rminHist + ((REAL)n - 0.5)/gd->i_deltaR;
#endif
    return rbin;
}

local int cb3d_print_zeta(struct cmdline_data* cmd, struct global_data* gd,
                          REAL *num, REAL *den, INTEGER npivots, INTEGER nbb)
{
    char namebuf[MAXLENGTHOFFILES];
    stream outstr = NULL;
    int write_failed;
    int ell, b1, b2;
    int nb = cmd->sizeHistN;

    if (format_checked(namebuf, sizeof(namebuf), "3D 3PCF output path",
                       "%s_3d%s", gd->fpfnamehistZetaMFileName,
                       EXTFILES) != 0) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "3D 3PCF output path is too long");
        return FAILURE;
    }
    outstr = fopen(namebuf, "w");
    if (outstr == NULL) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "cannot open 3D 3PCF output '%s': %s",
                 namebuf, strerror(errno));
        return FAILURE;
    }

    fprintf(outstr, "# cBalls 3D scalar 3PCF multipoles\n");
    fprintf(outstr, "# estimator: zeta_l = numerator/denominator; no random/edge correction\n");
    fprintf(outstr, "# npivots %" INTEGER_FMT
                    "  accepted_neighbor_visits %" INTEGER_FMT "\n",
            npivots, nbb);
    fprintf(outstr, "# columns: ell bin1 bin2 r1_center r2_center zeta numerator denominator\n");

    for (ell = 0; ell <= cmd->mChebyshev; ell++) {
        for (b1 = 0; b1 < nb; b1++) {
            for (b2 = 0; b2 < nb; b2++) {
                size_t idx = CB3D_ACC_INDEX(ell, b1, b2, nb);
                REAL zeta = cballs_normalize_or_zero(num[idx], den[idx]);
                fprintf(outstr, "%d %d %d %.12e %.12e %.12e %.12e %.12e\n",
                        ell, b1+1, b2+1,
                        cb3d_bin_center(cmd, gd, b1),
                        cb3d_bin_center(cmd, gd, b2),
                        zeta, num[idx], den[idx]);
            }
        }
    }
    write_failed = ferror(outstr);
    if (fclose(outstr) != 0) write_failed = TRUE;
    if (write_failed) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "failed writing 3D 3PCF output '%s'", namebuf);
        return FAILURE;
    }

    verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                        "\nPrinted 3D scalar 3PCF multipoles to %s\n", namebuf);
    return SUCCESS;
}

local int cb3d_print_xi(struct cmdline_data* cmd, struct global_data* gd,
                        REAL *num, REAL *den, INTEGER npivots, INTEGER nbb)
{
    char namebuf[MAXLENGTHOFFILES];
    stream outstr = NULL;
    int write_failed;
    int b;

    if (format_checked(namebuf, sizeof(namebuf), "3D 2PCF output path",
                       "%s_3d%s", gd->fpfnamehistXi2pcfFileName,
                       EXTFILES) != 0) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "3D 2PCF output path is too long");
        return FAILURE;
    }
    outstr = fopen(namebuf, "w");
    if (outstr == NULL) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "cannot open 3D 2PCF output '%s': %s",
                 namebuf, strerror(errno));
        return FAILURE;
    }

    fprintf(outstr, "# cBalls 3D scalar 2PCF\n");
    fprintf(outstr, "# estimator: xi = numerator/denominator; no random/edge correction\n");
    fprintf(outstr, "# directed pivot-neighbor visits give the same auto-catalog ratio as unordered pairs\n");
    fprintf(outstr, "# npivots %" INTEGER_FMT
                    "  accepted_neighbor_visits %" INTEGER_FMT "\n",
            npivots, nbb);
    fprintf(outstr, "# columns: bin r_center xi numerator denominator\n");
    for (b = 0; b < cmd->sizeHistN; b++) {
        fprintf(outstr, "%d %.12e %.12e %.12e %.12e\n",
                b + 1, cb3d_bin_center(cmd, gd, b),
                cballs_normalize_or_zero(num[b], den[b]), num[b], den[b]);
    }

    write_failed = ferror(outstr);
    if (fclose(outstr) != 0) write_failed = TRUE;
    if (write_failed) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "failed writing 3D 2PCF output '%s'", namebuf);
        return FAILURE;
    }

    verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                        "\nPrinted 3D scalar 2PCF to %s\n", namebuf);
    return SUCCESS;
}

local void cb3d_add_arrays(REAL *dst, const REAL *src, size_t n)
{
    size_t i;
    for (i = 0; i < n; i++) dst[i] += src[i];
}

local void cb3d_zero_array(REAL *dst, size_t n)
{
    size_t i;
    for (i = 0; i < n; i++) dst[i] = 0.0;
}

#undef CB3D_LM_INDEX
#undef CB3D_ACC_INDEX
#undef CB3D_ALM_INDEX
