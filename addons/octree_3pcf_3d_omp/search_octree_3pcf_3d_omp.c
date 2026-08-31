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

#include <limits.h>

#ifndef CB3D_OMP_PIVOT_BLOCK_SIZE
#define CB3D_OMP_PIVOT_BLOCK_SIZE 64
#endif
#if CB3D_OMP_PIVOT_BLOCK_SIZE < 1
#error CB3D_OMP_PIVOT_BLOCK_SIZE must be positive
#endif

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

typedef struct {
    REAL *num;
    REAL *den;
    REAL *xi_num;
    REAL *xi_den;
    INTEGER nbb;
    INTEGER npivots;
} cb3d_result, *cb3d_resultptr;

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
local int cb3d_init_result(struct cmdline_data* cmd, cb3d_resultptr result,
                           int compute_xi, int compute_zeta);
local void cb3d_free_result(cb3d_resultptr result);
local int cb3d_measure(struct cmdline_data* cmd, struct global_data* gd,
                       bodyptr *btable, INTEGER *nbody,
                       INTEGER ipmin, INTEGER ipmax, int cat,
                       int compute_xi, int compute_zeta,
                       int exclude_same_los, cb3d_resultptr result);
local int cb3d_validate(struct cmdline_data* cmd,
                        int *compute_xi, int *compute_zeta,
                        int survey_mode);
local int cb3d_print_survey_zeta(struct cmdline_data* cmd,
                                 struct global_data* gd,
                                 const cb3d_resultptr numerator,
                                 const cb3d_resultptr randoms,
                                 REAL random_scale);
local int cb3d_print_survey_xi(struct cmdline_data* cmd,
                               struct global_data* gd,
                               const cb3d_resultptr numerator,
                               const cb3d_resultptr randoms,
                               REAL random_scale);
local double cb3d_wigner_3j_zero_squared(int l1, int l2, int l3);
local int cb3d_solve_real_system(double *matrix, double *rhs, int n,
                                  double *condition_proxy);

/* ------------------------------------------------------------------------- */

global int cb3d_prepare_common_frame(struct cmdline_data* cmd,
                                     struct global_data* gd,
                                     bodyptr *btable, INTEGER *nbody)
{
    compute_vector minimum;
    compute_vector maximum;
    compute_vector center;
    bodyptr p;
    int ifile;
    int axis;
    int initialized = FALSE;

    if (cmd == NULL)
        return FAILURE;
    if (gd == NULL || btable == NULL || nbody == NULL || gd->ninfiles != 2) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "octree-3pcf-3d-omp: survey common-frame preparation "
                 "requires two catalogs");
        return FAILURE;
    }
    for (ifile = 0; ifile < gd->ninfiles; ifile++) {
        if (btable[ifile] == NULL || nbody[ifile] <= 0) {
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     "octree-3pcf-3d-omp: empty survey catalog %d", ifile);
            return FAILURE;
        }
        DO_BODY(p, btable[ifile], btable[ifile] + nbody[ifile]) {
            DO_COORD(axis) {
                REAL coordinate = Pos(p)[axis];
                if (!isfinite((double)coordinate)) {
                    snprintf(cmd->error_message, _ERRORMSGSIZE_,
                             "octree-3pcf-3d-omp: non-finite coordinate in "
                             "survey catalog %d", ifile);
                    return FAILURE;
                }
                if (!initialized) {
                    minimum[axis] = coordinate;
                    maximum[axis] = coordinate;
                } else {
                    if (coordinate < minimum[axis])
                        minimum[axis] = coordinate;
                    if (coordinate > maximum[axis])
                        maximum[axis] = coordinate;
                }
            }
            initialized = TRUE;
        }
    }
    DO_COORD(axis)
        center[axis] = 0.5*(minimum[axis] + maximum[axis]);
    for (ifile = 0; ifile < gd->ninfiles; ifile++)
        DO_BODY(p, btable[ifile], btable[ifile] + nbody[ifile])
            DO_COORD(axis)
                Pos(p)[axis] -= center[axis];

    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
        "octree-3pcf-3d-omp: centered data and random catalogs in one "
        "common Cartesian frame\n");
    return SUCCESS;
}

global int cb3d_prepare_survey_catalogs(struct cmdline_data* cmd,
                                        struct global_data* gd,
                                        bodyptr *btable, INTEGER *nbody,
                                        int data_cat, int random_cat,
                                        int *dmr_cat, int *normalized_random_cat,
                                        REAL *random_scale)
{
    INTEGER i;
    INTEGER nd;
    INTEGER nr;
    size_t total;
    size_t first_bytes;
    size_t second_bytes;
    size_t integer_max;
    int first_scratch;
    double data_weight = 0.0;
    double random_weight = 0.0;
    double alpha;

    if (data_cat < 0 || data_cat >= gd->ninfiles
        || random_cat < 0 || random_cat >= gd->ninfiles
        || data_cat == random_cat) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "octree-3pcf-3d-omp: invalid data/random catalog indices");
        return FAILURE;
    }
    if (cballs_opt_remove_mean(cmd)) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "octree-3pcf-3d-omp: remove-mean is incompatible with "
                 "survey-estimator-3d");
        return FAILURE;
    }
    if (cb3d_exclude_same_los_option(cmd)) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "octree-3pcf-3d-omp: LOS exclusion is not defined for the "
                 "survey point-process estimator");
        return FAILURE;
    }

    nd = nbody[data_cat];
    nr = nbody[random_cat];
    if (nd < 2 || nr < 2) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "octree-3pcf-3d-omp: survey mode requires at least two data "
                 "and two random points");
        return FAILURE;
    }
    if ((size_t)nd > SIZE_MAX - (size_t)nr) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "octree-3pcf-3d-omp: combined catalog size overflow");
        return FAILURE;
    }
    total = (size_t)nd + (size_t)nr;
#ifdef LONGINT
    integer_max = (size_t)LONG_MAX;
#else
    integer_max = (size_t)INT_MAX;
#endif
    if (total > integer_max || (size_t)nr > integer_max
        || total > SIZE_MAX/sizeof(body)
        || (size_t)nr > SIZE_MAX/sizeof(body)) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "octree-3pcf-3d-omp: survey workspace size is not "
                 "representable");
        return FAILURE;
    }
    first_bytes = total*sizeof(body);
    second_bytes = (size_t)nr*sizeof(body);
    if (first_bytes > integer_max
        || second_bytes > integer_max - first_bytes
        || gd->bytes_tot < 0
        || (size_t)gd->bytes_tot > integer_max - first_bytes - second_bytes) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "octree-3pcf-3d-omp: survey workspace byte accounting "
                 "would overflow");
        return FAILURE;
    }
    first_scratch = gd->ninfiles;
    if (first_scratch < 0 || first_scratch + 1 >= MAXITEMS) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "octree-3pcf-3d-omp: no catalog slots available for survey workspaces");
        return FAILURE;
    }

    for (i = 0; i < nd; i++) {
        double weight = (double)Weight(btable[data_cat] + i);
        if (!isfinite(weight) || weight < 0.0) {
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     "octree-3pcf-3d-omp: data weight %" INTEGER_FMT
                     " is negative or non-finite", i + 1);
            return FAILURE;
        }
        data_weight += weight;
    }
    for (i = 0; i < nr; i++) {
        double weight = (double)Weight(btable[random_cat] + i);
        if (!isfinite(weight) || weight < 0.0) {
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     "octree-3pcf-3d-omp: random weight %" INTEGER_FMT
                     " is negative or non-finite", i + 1);
            return FAILURE;
        }
        random_weight += weight;
    }
    if (!isfinite(data_weight) || !isfinite(random_weight)
        || data_weight <= 0.0 || random_weight <= 0.0) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "octree-3pcf-3d-omp: data and random weight sums must be "
                 "finite and positive");
        return FAILURE;
    }
    alpha = data_weight/random_weight;

    if (cballs_calloc_checked((void **)&btable[first_scratch], total,
                              sizeof(body), "3D survey D-R catalog",
                              cmd->error_message, _ERRORMSGSIZE_) == FAILURE)
        return FAILURE;
    nbody[first_scratch] = (INTEGER)total;
    gd->rSizeTable[first_scratch] = 1.0;
    gd->ninfiles = first_scratch + 1;
    gd->bodytable_allocated = TRUE;
    gd->bytes_tot += (INTEGER)first_bytes;

    if (cballs_calloc_checked((void **)&btable[first_scratch + 1], (size_t)nr,
                              sizeof(body), "3D survey random catalog",
                              cmd->error_message, _ERRORMSGSIZE_) == FAILURE)
        return FAILURE;
    nbody[first_scratch + 1] = nr;
    gd->rSizeTable[first_scratch + 1] = 1.0;
    gd->ninfiles = first_scratch + 2;
    gd->bytes_tot += (INTEGER)second_bytes;

    for (i = 0; i < nd; i++) {
        bodyptr source = btable[data_cat] + i;
        bodyptr target = btable[first_scratch] + i;
        int k;
        Type(target) = BODY;
        Update(target) = TRUE;
        Mask(target) = MASK_NODE_VALID;
        Mass(target) = 1.0;
        Kappa(target) = 1.0;
        Weight(target) = Weight(source);
        Id(target) = i + 1;
        Octree3pcf3dLosId(target) = i + 1;
        DO_COORD(k) Pos(target)[k] = Pos(source)[k];
    }
    for (i = 0; i < nr; i++) {
        bodyptr source = btable[random_cat] + i;
        bodyptr target_dmr = btable[first_scratch] + nd + i;
        bodyptr target_random = btable[first_scratch + 1] + i;
        REAL normalized_weight = (REAL)(alpha*(double)Weight(source));
        int k;

        Type(target_dmr) = Type(target_random) = BODY;
        Update(target_dmr) = Update(target_random) = TRUE;
        Mask(target_dmr) = Mask(target_random) = MASK_NODE_VALID;
        Mass(target_dmr) = Mass(target_random) = 1.0;
        Kappa(target_dmr) = Kappa(target_random) = 1.0;
        Weight(target_dmr) = -normalized_weight;
        Weight(target_random) = normalized_weight;
        Id(target_dmr) = nd + i + 1;
        Id(target_random) = i + 1;
        Octree3pcf3dLosId(target_dmr) = nd + i + 1;
        Octree3pcf3dLosId(target_random) = i + 1;
        DO_COORD(k) {
            Pos(target_dmr)[k] = Pos(source)[k];
            Pos(target_random)[k] = Pos(source)[k];
        }
    }

    *dmr_cat = first_scratch;
    *normalized_random_cat = first_scratch + 1;
    *random_scale = (REAL)alpha;
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
        "octree-3pcf-3d-omp survey weights: sum(D)=%.12e, "
        "sum(R_raw)=%.12e, alpha=%.12e, sum(D-R)=%.12e\n",
        data_weight, random_weight, alpha,
        data_weight-alpha*random_weight);
    return SUCCESS;
}

global int searchcalc_octree_3pcf_3d_omp(struct cmdline_data* cmd,
                                         struct global_data* gd,
                                         bodyptr *btable, INTEGER *nbody,
                                         INTEGER ipmin, INTEGER *ipmax,
                                         int cat1, int cat2)
{
    double cpustart = CPUTIME;
    cb3d_result result = {0};
    int compute_xi = FALSE;
    int compute_zeta = FALSE;
    int status = FAILURE;

    if (cat1 != cat2) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: the 3D scalar multipole estimator requires an auto-catalog",
                 cmd->searchMethod);
        goto setup_done;
    }
    if (cb3d_validate(cmd, &compute_xi, &compute_zeta, FALSE) == FAILURE)
        goto setup_done;
    if (cb3d_init_result(cmd, &result, compute_xi, compute_zeta) == FAILURE)
        goto setup_done;
    status = SUCCESS;

setup_done:
    status = cb3d_parallel_consensus(cmd, status, "3D scalar histogram setup");
    if (status == FAILURE) goto cleanup;
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
        "\n%s: exact 3D scalar estimator; nbins=%d lmax=%d 2pcf=%d 3pcf=%d; "
        "physical cutoff=%g (global Rcut=%g)\n",
        cmd->searchMethod, cmd->sizeHistN, cmd->mChebyshev,
        compute_xi, compute_zeta, cmd->rangeN, gd->Rcut);
    status = cb3d_measure(cmd, gd, btable, nbody, ipmin, ipmax[cat1], cat1,
                          compute_xi, compute_zeta,
                          cb3d_exclude_same_los_option(cmd), &result);
    if (status == FAILURE) goto cleanup;
    if (cb3d_parallel_publish(cmd) && !cballs_opt_no_out_hist(cmd)) {
        if (compute_zeta)
            status = cb3d_print_zeta(cmd, gd, result.num, result.den,
                                     result.npivots, result.nbb);
        if (status == SUCCESS && compute_xi)
            status = cb3d_print_xi(cmd, gd, result.xi_num, result.xi_den,
                                   result.npivots, result.nbb);
    }
    status = cb3d_parallel_consensus(cmd, status, "3D scalar output");

cleanup:
    gd->cpusearch = CPUTIME - cpustart;
    cb3d_free_result(&result);
    return status;
}

global int searchcalc_octree_3pcf_3d_omp_survey(
    struct cmdline_data* cmd, struct global_data* gd,
    bodyptr *btable, INTEGER *nbody,
    int dmr_cat, int normalized_random_cat, REAL random_scale)
{
    double cpustart = CPUTIME;
    cb3d_result numerator = {0};
    cb3d_result randoms = {0};
    int compute_xi = FALSE;
    int compute_zeta = FALSE;
    int status = FAILURE;

    if (cb3d_validate(cmd, &compute_xi, &compute_zeta, TRUE) == FAILURE)
        goto setup_done;
    if (cb3d_init_result(cmd, &numerator, compute_xi, compute_zeta) == FAILURE)
        goto setup_done;
    if (cb3d_init_result(cmd, &randoms, compute_xi, compute_zeta) == FAILURE)
        goto setup_done;
    status = SUCCESS;

setup_done:
    status = cb3d_parallel_consensus(cmd, status, "3D survey histogram setup");
    if (status == FAILURE) goto cleanup;
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
        "\n%s: ENCORE-style survey estimator; measuring D-R and random "
        "multipoles through ell=%d\n", cmd->searchMethod, cmd->mChebyshev);
    status = cb3d_measure(cmd, gd, btable, nbody, 1, nbody[dmr_cat], dmr_cat,
                          compute_xi, compute_zeta, FALSE, &numerator);
    if (status == FAILURE) goto cleanup;
    status = cb3d_measure(cmd, gd, btable, nbody, 1,
                          nbody[normalized_random_cat], normalized_random_cat,
                          compute_xi, compute_zeta, FALSE, &randoms);
    if (status == FAILURE) goto cleanup;

    /* Solve the window system only after both globally summed measurements. */
    if (cb3d_parallel_publish(cmd) && !cballs_opt_no_out_hist(cmd)) {
        if (compute_zeta)
            status = cb3d_print_survey_zeta(cmd, gd, &numerator, &randoms,
                                            random_scale);
        if (status == SUCCESS && compute_xi)
            status = cb3d_print_survey_xi(cmd, gd, &numerator, &randoms,
                                          random_scale);
    }
    status = cb3d_parallel_consensus(cmd, status, "3D survey edge correction/output");

cleanup:
    gd->cpusearch = CPUTIME - cpustart;
    cb3d_free_result(&numerator);
    cb3d_free_result(&randoms);
    return status;
}

/* ------------------------------------------------------------------------- */

local int cb3d_validate(struct cmdline_data* cmd,
                        int *compute_xi, int *compute_zeta,
                        int survey_mode)
{
    const char *routine_name = "octree-3pcf-3d-omp";
    size_t nb;
    size_t lcount;

#if NDIM != 3
    snprintf(cmd->error_message, _ERRORMSGSIZE_,
             "%s: this estimator requires NDIM=3", routine_name);
    return FAILURE;
#endif

    if (cb3d_get_modes(cmd, compute_xi, compute_zeta) == FAILURE)
        return FAILURE;
    if (cmd->sizeHistN <= 0) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: sizeHistN must be positive", routine_name);
        return FAILURE;
    }
    if (cmd->mChebyshev < 0
        || cmd->mChebyshev > OCTREE3PCF3D_MAX_LMAX) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: invalid mChebyshev/lmax=%d; supported range is 0..%d",
                 routine_name, cmd->mChebyshev, OCTREE3PCF3D_MAX_LMAX);
        return FAILURE;
    }
    if (survey_mode && *compute_zeta && cmd->mChebyshev < 1) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: survey 3PCF requires mChebyshev>=1 so one extra "
                 "window-coupling multipole can be measured", routine_name);
        return FAILURE;
    }
    if (!isfinite((double)cmd->rangeN) || cmd->rangeN <= 0.0
        || !isfinite((double)cmd->rminHist) || cmd->rminHist < 0.0
        || cmd->rminHist >= cmd->rangeN) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: require finite 0<=rminHist<rangeN", routine_name);
        return FAILURE;
    }

    nb = (size_t)cmd->sizeHistN;
    lcount = (size_t)cmd->mChebyshev + 1;
    if (nb > SIZE_MAX/nb || lcount > SIZE_MAX/(nb*nb)) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: histogram dimensions overflow size_t", routine_name);
        return FAILURE;
    }
    if (cmd->useLogHist && cmd->rminHist <= 0.0) {
        verb_print_warning(cmd->verbose,
            "%s: log histograms with rminHist<=0 use cBalls legacy bin "
            "centers; consider rminHist>0 for 3D.", routine_name);
    }
    return SUCCESS;
}

local int cb3d_init_result(struct cmdline_data* cmd, cb3d_resultptr result,
                           int compute_xi, int compute_zeta)
{
    size_t nb = (size_t)cmd->sizeHistN;
    size_t nacc = ((size_t)cmd->mChebyshev + 1)*nb*nb;

    memset(result, 0, sizeof(*result));
    if (compute_zeta) {
        if (cballs_calloc_checked((void **)&result->num, nacc,
                                  sizeof(*result->num), "3D 3PCF numerator",
                                  cmd->error_message,
                                  _ERRORMSGSIZE_) == FAILURE)
            goto fail;
        if (cballs_calloc_checked((void **)&result->den, nacc,
                                  sizeof(*result->den), "3D 3PCF denominator",
                                  cmd->error_message,
                                  _ERRORMSGSIZE_) == FAILURE)
            goto fail;
    }
    if (compute_xi) {
        if (cballs_calloc_checked((void **)&result->xi_num, nb,
                                  sizeof(*result->xi_num), "3D 2PCF numerator",
                                  cmd->error_message,
                                  _ERRORMSGSIZE_) == FAILURE)
            goto fail;
        if (cballs_calloc_checked((void **)&result->xi_den, nb,
                                  sizeof(*result->xi_den), "3D 2PCF denominator",
                                  cmd->error_message,
                                  _ERRORMSGSIZE_) == FAILURE)
            goto fail;
    }
    return SUCCESS;

fail:
    cb3d_free_result(result);
    return FAILURE;
}

local void cb3d_free_result(cb3d_resultptr result)
{
    if (result == NULL) return;
    free(result->num);
    free(result->den);
    free(result->xi_num);
    free(result->xi_den);
    memset(result, 0, sizeof(*result));
}

local int cb3d_measure(struct cmdline_data* cmd, struct global_data* gd,
                       bodyptr *btable, INTEGER *nbody,
                       INTEGER ipmin, INTEGER ipmax, int cat,
                       int compute_xi, int compute_zeta,
                       int exclude_same_los, cb3d_resultptr result)
{
    const size_t nb = (size_t)cmd->sizeHistN;
    const size_t nacc = ((size_t)cmd->mChebyshev + 1)*nb*nb;
    const INTEGER first_block = cb3d_parallel_first(cmd);
    const INTEGER block_stride = cb3d_parallel_stride(cmd);
    int allocation_failed = FALSE;
    int status = SUCCESS;
    ErrorMsg allocation_error = "";

    if (cat < 0 || cat >= gd->ninfiles || btable[cat] == NULL
        || nbody[cat] <= 0 || roottable[cat] == NULL) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "cb3d_measure: catalog %d is not ready for tree traversal", cat);
        status = FAILURE;
    } else if (ipmin < 1 || ipmax < ipmin || ipmax > nbody[cat]) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "cb3d_measure: invalid pivot interval");
        status = FAILURE;
    }
    status = cb3d_parallel_consensus(cmd, status, "3D traversal validation");
    if (status == FAILURE) return FAILURE;

#ifdef OPENMPCODE
    ThreadCount(cmd, gd, nbody[cat], cat);
#else
#error OPENMPMACHINE is not defined. Switch it on in Makefile_settings
#endif
    const INTEGER pivot_count = ipmax - ipmin + 1;
    const INTEGER block_count = 1 + (pivot_count - 1)/CB3D_OMP_PIVOT_BLOCK_SIZE;

#pragma omp parallel default(none) \
    shared(cmd,gd,btable,roottable,ipmin,ipmax,cat,result,nacc,nb, \
           compute_xi,compute_zeta,exclude_same_los,allocation_failed, \
           allocation_error,first_block,block_stride,block_count)
    {
        cb3d_hist hist;
        ErrorMsg thread_error = "";
        int hist_ready = cb3d_init_hist(cmd, gd, &hist, compute_xi, compute_zeta,
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
#pragma omp for schedule(dynamic,1) ordered
        for (INTEGER block = first_block; block < block_count; block += block_stride) {
            if (!allocation_failed) {
                const INTEGER first = ipmin - 1 + block*CB3D_OMP_PIVOT_BLOCK_SIZE;
                const INTEGER count = MIN((INTEGER)CB3D_OMP_PIVOT_BLOCK_SIZE,
                                           ipmax - first);
                if (compute_zeta) {
                    cb3d_zero_array(hist.num, nacc);
                    cb3d_zero_array(hist.den, nacc);
                }
                if (compute_xi) {
                    cb3d_zero_array(hist.xi_num, nb);
                    cb3d_zero_array(hist.xi_den, nb);
                }
                hist.nbbcalcthread = 0;
                hist.npivots_thread = 0;
                for (INTEGER offset = 0; offset < count; offset++) {
                    bodyptr p = btable[cat] + first + offset;
                    if (Update(p) == FALSE) continue;
                    if (hist.read_mask && Mask(p) == FALSE) continue;
                    cb3d_clear_pivot(&hist);
                    hist.pivot_weight = Weight(p);
                    hist.pivot_field = Weight(p)*Kappa(p);
                    if (hist.exclude_same_los)
                        hist.pivot_los = Octree3pcf3dLosId(p);
                    cb3d_walktree(gd, p, (nodeptr)roottable[cat], &hist);
                    cb3d_accumulate_pivot(&hist);
                    hist.npivots_thread++;
                }
            }
#pragma omp ordered
            {
                if (!allocation_failed) {
                    if (compute_zeta) {
                        cb3d_add_arrays(result->num, hist.num, nacc);
                        cb3d_add_arrays(result->den, hist.den, nacc);
                    }
                    if (compute_xi) {
                        cb3d_add_arrays(result->xi_num, hist.xi_num, nb);
                        cb3d_add_arrays(result->xi_den, hist.xi_den, nb);
                    }
                    result->nbb += hist.nbbcalcthread;
                    result->npivots += hist.npivots_thread;
                }
            }
        }
        if (hist_ready) cb3d_free_hist(&hist);
    }
    if (allocation_failed)
        snprintf(cmd->error_message, _ERRORMSGSIZE_, "cb3d_measure: %s",
                 allocation_error[0] ? allocation_error : "worker allocation failed");
    status = cb3d_parallel_consensus(cmd, allocation_failed ? FAILURE : SUCCESS,
                                     "3D traversal workers");
    if (status == FAILURE) return FAILURE;
    INTEGER counters[2] = {result->nbb, result->npivots};
    if ((compute_zeta
         && (cb3d_parallel_reduce_reals(cmd, result->num, nacc) == FAILURE
             || cb3d_parallel_reduce_reals(cmd, result->den, nacc) == FAILURE))
        || (compute_xi
            && (cb3d_parallel_reduce_reals(cmd, result->xi_num, nb) == FAILURE
                || cb3d_parallel_reduce_reals(cmd, result->xi_den, nb) == FAILURE))
        || cb3d_parallel_reduce_integers(cmd, counters, 2) == FAILURE)
        return FAILURE;
    result->nbb = counters[0];
    result->npivots = counters[1];
    return SUCCESS;
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

local double cb3d_wigner_3j_zero_squared(int l1, int l2, int l3)
{
    int sum = l1 + l2 + l3;
    int g;
    double log_value;

    if (l1 < 0 || l2 < 0 || l3 < 0 || (sum & 1)
        || l1 + l2 < l3 || l1 + l3 < l2 || l2 + l3 < l1)
        return 0.0;
    g = sum/2;
    log_value =
          2.0*lgamma((double)g + 1.0)
        - 2.0*lgamma((double)(g-l1) + 1.0)
        - 2.0*lgamma((double)(g-l2) + 1.0)
        - 2.0*lgamma((double)(g-l3) + 1.0)
        + lgamma((double)(sum-2*l1) + 1.0)
        + lgamma((double)(sum-2*l2) + 1.0)
        + lgamma((double)(sum-2*l3) + 1.0)
        - lgamma((double)sum + 2.0);
    return exp(log_value);
}

local int cb3d_solve_real_system(double *matrix, double *rhs, int n,
                                  double *condition_proxy)
{
    int i, j, k;
    double matrix_scale = 0.0;
    double largest_pivot = 0.0;
    double smallest_pivot = DBL_MAX;
    double tolerance;

    if (matrix == NULL || rhs == NULL || n <= 0)
        return FAILURE;
    for (i = 0; i < n*n; i++) {
        double value = fabs(matrix[i]);
        if (!isfinite(value)) return FAILURE;
        if (value > matrix_scale) matrix_scale = value;
    }
    for (i = 0; i < n; i++)
        if (!isfinite(rhs[i])) return FAILURE;
    tolerance = 128.0*DBL_EPSILON*(double)n
              * (matrix_scale > 1.0 ? matrix_scale : 1.0);

    for (k = 0; k < n; k++) {
        int pivot_row = k;
        double pivot_size = fabs(matrix[k*n+k]);

        for (i = k + 1; i < n; i++) {
            double candidate = fabs(matrix[i*n+k]);
            if (candidate > pivot_size) {
                pivot_size = candidate;
                pivot_row = i;
            }
        }
        if (!isfinite(pivot_size) || pivot_size <= tolerance)
            return FAILURE;
        if (pivot_row != k) {
            for (j = k; j < n; j++) {
                double temporary = matrix[k*n+j];
                matrix[k*n+j] = matrix[pivot_row*n+j];
                matrix[pivot_row*n+j] = temporary;
            }
            {
                double temporary = rhs[k];
                rhs[k] = rhs[pivot_row];
                rhs[pivot_row] = temporary;
            }
        }

        pivot_size = fabs(matrix[k*n+k]);
        if (pivot_size > largest_pivot) largest_pivot = pivot_size;
        if (pivot_size < smallest_pivot) smallest_pivot = pivot_size;
        for (i = k + 1; i < n; i++) {
            double factor = matrix[i*n+k]/matrix[k*n+k];
            matrix[i*n+k] = 0.0;
            for (j = k + 1; j < n; j++)
                matrix[i*n+j] -= factor*matrix[k*n+j];
            rhs[i] -= factor*rhs[k];
        }
    }

    for (i = n - 1; i >= 0; i--) {
        double value = rhs[i];
        for (j = i + 1; j < n; j++)
            value -= matrix[i*n+j]*rhs[j];
        if (fabs(matrix[i*n+i]) <= tolerance)
            return FAILURE;
        rhs[i] = value/matrix[i*n+i];
        if (!isfinite(rhs[i])) return FAILURE;
    }
    if (condition_proxy != NULL)
        *condition_proxy = smallest_pivot > 0.0
            ? largest_pivot/smallest_pivot : HUGE_VAL;
    return SUCCESS;
}

local int cb3d_print_survey_zeta(struct cmdline_data* cmd,
                                 struct global_data* gd,
                                 const cb3d_resultptr numerator,
                                 const cb3d_resultptr randoms,
                                 REAL random_scale)
{
    char namebuf[MAXLENGTHOFFILES];
    stream outstr = NULL;
    double *matrix = NULL;
    double *rhs = NULL;
    double *n_encore = NULL;
    double *r_encore = NULL;
    double *window = NULL;
    int nb = cmd->sizeHistN;
    int measured_lmax = cmd->mChebyshev;
    int nmode = measured_lmax + 1;
    int output_lmax = scanopt(cmd->options, "survey-keep-top-multipole")
                    ? measured_lmax : measured_lmax - 1;
    int b1, b2, ell, ell_window, ell_signal;
    int invalid_pairs = 0;
    int write_failed = FALSE;
    int status = FAILURE;

    if (format_checked(namebuf, sizeof(namebuf), "3D survey 3PCF output path",
                       "%s_3d_survey%s", gd->fpfnamehistZetaMFileName,
                       EXTFILES) != 0) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "3D survey 3PCF output path is too long");
        return FAILURE;
    }
    if (cballs_malloc_checked((void **)&matrix, (size_t)nmode*(size_t)nmode,
                              sizeof(*matrix), "3D survey coupling matrix",
                              cmd->error_message, _ERRORMSGSIZE_) == FAILURE)
        goto cleanup;
    if (cballs_malloc_checked((void **)&rhs, (size_t)nmode, sizeof(*rhs),
                              "3D survey coupling right-hand side",
                              cmd->error_message, _ERRORMSGSIZE_) == FAILURE)
        goto cleanup;
    if (cballs_malloc_checked((void **)&n_encore, (size_t)nmode,
                              sizeof(*n_encore), "3D survey D-R multipoles",
                              cmd->error_message, _ERRORMSGSIZE_) == FAILURE)
        goto cleanup;
    if (cballs_malloc_checked((void **)&r_encore, (size_t)nmode,
                              sizeof(*r_encore), "3D survey random multipoles",
                              cmd->error_message, _ERRORMSGSIZE_) == FAILURE)
        goto cleanup;
    if (cballs_malloc_checked((void **)&window, (size_t)nmode,
                              sizeof(*window), "3D survey window multipoles",
                              cmd->error_message, _ERRORMSGSIZE_) == FAILURE)
        goto cleanup;

    outstr = fopen(namebuf, "w");
    if (outstr == NULL) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "cannot open 3D survey 3PCF output '%s': %s",
                 namebuf, strerror(errno));
        goto cleanup;
    }
    fprintf(outstr, "# cBalls ENCORE-style survey 3PCF edge correction\n");
    fprintf(outstr, "# estimator: N=(D-alpha*R)^3; window=R^3; "
                    "M zeta_ENCORE=N_ENCORE/R0_ENCORE\n");
    fprintf(outstr, "# alpha %.17e  measured_lmax %d  output_lmax %d\n",
            (double)random_scale, measured_lmax, output_lmax);
    fprintf(outstr, "# only distinct radial shells bin1<bin2 are emitted\n");
    fprintf(outstr, "# standard coefficient multiplies P_ell(mu); "
                    "ENCORE coefficient multiplies "
                    "(-1)^ell*sqrt(2ell+1)*P_ell(mu)/(4*pi)\n");
    fprintf(outstr, "# D-R pivots %" INTEGER_FMT
                    " visits %" INTEGER_FMT
                    "  R pivots %" INTEGER_FMT
                    " visits %" INTEGER_FMT "\n",
            numerator->npivots, numerator->nbb,
            randoms->npivots, randoms->nbb);
    fprintf(outstr, "# zero policy: corrected multipoles are 0 and valid=0 "
                    "when R0 is zero/non-finite or M is singular\n");
    fprintf(outstr, "# columns: ell bin1 bin2 r1_center r2_center "
                    "zeta_legendre zeta_encore N_legendre R_legendre "
                    "N_encore R_encore R0_encore pivot_condition valid\n");

    for (b1 = 0; b1 < nb; b1++) {
        for (b2 = b1 + 1; b2 < nb; b2++) {
            double r0_encore;
            double condition_proxy = HUGE_VAL;
            int valid = TRUE;

            for (ell = 0; ell <= measured_lmax; ell++) {
                size_t index = CB3D_ACC_INDEX(ell, b1, b2, nb);
                double sign = (ell & 1) ? -1.0 : 1.0;
                double conversion =
                    CB3D_FOURPI*sign*sqrt(2.0*(double)ell + 1.0);
                n_encore[ell] = (double)numerator->num[index]/conversion;
                r_encore[ell] = (double)randoms->num[index]/conversion;
            }
            r0_encore = r_encore[0];
            if (!isfinite(r0_encore) || r0_encore == 0.0) {
                valid = FALSE;
            } else {
                for (ell = 0; ell <= measured_lmax; ell++)
                    window[ell] = r_encore[ell]/r0_encore;
                for (ell_signal = 0; ell_signal <= measured_lmax;
                     ell_signal++) {
                    rhs[ell_signal] = n_encore[ell_signal]/r0_encore;
                    for (ell = 0; ell <= measured_lmax; ell++) {
                        double coupling = 0.0;
                        for (ell_window = 0;
                             ell_window <= measured_lmax; ell_window++) {
                            double wigner2 = cb3d_wigner_3j_zero_squared(
                                ell_signal, ell_window, ell);
                            if (wigner2 == 0.0) continue;
                            coupling += window[ell_window]
                                * wigner2
                                * sqrt((2.0*(double)ell_signal + 1.0)
                                     *(2.0*(double)ell_window + 1.0)
                                     *(2.0*(double)ell + 1.0))
                                / CB3D_FOURPI;
                        }
                        matrix[ell_signal*nmode + ell] = coupling;
                    }
                }
                if (cb3d_solve_real_system(matrix, rhs, nmode,
                                           &condition_proxy) == FAILURE)
                    valid = FALSE;
            }

            if (!valid) {
                invalid_pairs++;
                for (ell = 0; ell <= measured_lmax; ell++) rhs[ell] = 0.0;
            }
            for (ell = 0; ell <= output_lmax; ell++) {
                size_t index = CB3D_ACC_INDEX(ell, b1, b2, nb);
                double sign = (ell & 1) ? -1.0 : 1.0;
                double encore_to_legendre =
                    sign*sqrt(2.0*(double)ell + 1.0)/CB3D_FOURPI;
                fprintf(outstr,
                        "%d %d %d %.12e %.12e %.12e %.12e "
                        "%.12e %.12e %.12e %.12e %.12e %.12e %d\n",
                        ell, b1 + 1, b2 + 1,
                        cb3d_bin_center(cmd, gd, b1),
                        cb3d_bin_center(cmd, gd, b2),
                        rhs[ell]*encore_to_legendre, rhs[ell],
                        (double)numerator->num[index],
                        (double)randoms->num[index],
                        n_encore[ell], r_encore[ell], r0_encore,
                        condition_proxy, valid);
            }
        }
    }

    write_failed = ferror(outstr);
    if (fclose(outstr) != 0) write_failed = TRUE;
    outstr = NULL;
    if (write_failed) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "failed writing 3D survey 3PCF output '%s'", namebuf);
        goto cleanup;
    }
    if (invalid_pairs > 0)
        verb_print_warning(cmd->verbose,
            "octree-3pcf-3d-omp: %d radial-bin pairs had an empty or "
            "singular random-window correction and were written as invalid",
            invalid_pairs);
    verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                        "\nPrinted edge-corrected 3D survey 3PCF to %s\n",
                        namebuf);
    status = SUCCESS;

cleanup:
    if (outstr != NULL) fclose(outstr);
    free(matrix);
    free(rhs);
    free(n_encore);
    free(r_encore);
    free(window);
    return status;
}

local int cb3d_print_survey_xi(struct cmdline_data* cmd,
                               struct global_data* gd,
                               const cb3d_resultptr numerator,
                               const cb3d_resultptr randoms,
                               REAL random_scale)
{
    char namebuf[MAXLENGTHOFFILES];
    stream outstr = NULL;
    int b;
    int write_failed;

    if (format_checked(namebuf, sizeof(namebuf), "3D survey 2PCF output path",
                       "%s_3d_survey%s", gd->fpfnamehistXi2pcfFileName,
                       EXTFILES) != 0) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "3D survey 2PCF output path is too long");
        return FAILURE;
    }
    outstr = fopen(namebuf, "w");
    if (outstr == NULL) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "cannot open 3D survey 2PCF output '%s': %s",
                 namebuf, strerror(errno));
        return FAILURE;
    }
    fprintf(outstr, "# cBalls ENCORE/Landy-Szalay survey 2PCF\n");
    fprintf(outstr, "# estimator: xi=(D-alpha*R)^2/RR\n");
    fprintf(outstr, "# alpha %.17e\n", (double)random_scale);
    fprintf(outstr, "# D-R pivots %" INTEGER_FMT
                    " visits %" INTEGER_FMT
                    "  R pivots %" INTEGER_FMT
                    " visits %" INTEGER_FMT "\n",
            numerator->npivots, numerator->nbb,
            randoms->npivots, randoms->nbb);
    fprintf(outstr, "# zero policy: xi=0 and valid=0 when RR is zero or "
                    "non-finite\n");
    fprintf(outstr, "# columns: bin r_center xi N_DminusR_squared RR valid\n");
    for (b = 0; b < cmd->sizeHistN; b++) {
        double nvalue = (double)numerator->xi_num[b];
        double rvalue = (double)randoms->xi_num[b];
        int valid = isfinite(nvalue) && isfinite(rvalue)
                 && rvalue != 0.0;
        double xi = valid ? nvalue/rvalue : 0.0;
        fprintf(outstr, "%d %.12e %.12e %.12e %.12e %d\n",
                b + 1, cb3d_bin_center(cmd, gd, b),
                xi, nvalue, rvalue, valid);
    }

    write_failed = ferror(outstr);
    if (fclose(outstr) != 0) write_failed = TRUE;
    if (write_failed) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "failed writing 3D survey 2PCF output '%s'", namebuf);
        return FAILURE;
    }
    verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                        "\nPrinted edge-corrected 3D survey 2PCF to %s\n",
                        namebuf);
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
