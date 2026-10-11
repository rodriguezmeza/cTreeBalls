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
     pivot-neighbor pairs with equal LOS IDs.
   - exclude-all-same-los additionally subtracts neighbor-neighbor products
     with equal LOS IDs. It is the Ly-alpha contract in which all three
     members of an accepted triplet belong to distinct forests.
 =============================================================================*/

#include "globaldefs.h"
#include "distinct_products.h"
#include "legendre_window_solver.h"
#include "physical_survey_config.h"

#include <limits.h>
#ifdef BALLS4SCANLEV
#include "octree_scan_frontier.h"
#endif

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
    bodyptr body;
    int bin;
} cb3d_neighbor;

typedef struct {
    int nbins;
    int lmax;
    int nlm;
    REAL *alm_re;
    REAL *alm_im;
    REAL *alm_distinct;
    REAL *ylm_norm;
    REAL *ylm_a;
    REAL *ylm_b;
    REAL *ylm_diag;
    REAL *shell_w;
    REAL *shell_weight_distinct;
    REAL *distinct_los_power;
    REAL *distinct_los_weight;
    REAL *group_alm_re;
    REAL *group_alm_im;
    REAL *group_shell_w;
    REAL *num;
    REAL *den;
    REAL *xi_num;
    REAL *xi_den;
    int *active_bins;
    int *group_active_bins;
    unsigned char *bin_active;
    unsigned char *group_bin_active;
    cb3d_neighbor *neighbors;
    size_t neighbor_count;
    size_t neighbor_capacity;
    int nactive;
    int group_nactive;
    int compute_xi;
    int compute_zeta;
    int exclude_same_los;
    int exclude_all_same_los;
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
    int lmax;
    REAL *num;
    REAL *den;
    REAL *xi_num;
    REAL *xi_den;
    INTEGER nbb;
    INTEGER npivots;
} cb3d_result, *cb3d_resultptr;

typedef struct {
    REAL *zeta, *encore, *valid, *status, *rcond, *residual;
    REAL *xi, *xi_valid, *xi_status;
} cb3d_correction;
local int cb3d_correct_survey(struct cmdline_data *, const cb3d_resultptr,
    const cb3d_resultptr, cb3d_correction *);
local void cb3d_free_correction(cb3d_correction *);

local int cb3d_init_hist(struct cmdline_data* cmd, struct global_data* gd,
                         cb3d_histptr h,
                         int compute_xi, int compute_zeta,
                         int exclude_same_los, int exclude_all_same_los, int lmax,
                         ErrorMsg allocation_error);
local void cb3d_free_hist(cb3d_histptr h);
local void cb3d_clear_pivot(cb3d_histptr h);
local void cb3d_activate_bin(cb3d_histptr h, int b);
local int cb3d_bin_index(cb3d_histptr h, REAL r);
local int cb3d_reject_cell(struct global_data* gd, bodyptr p, nodeptr q,
                           cb3d_histptr h);
local int cb3d_walktree(struct global_data* gd,
                        bodyptr p, nodeptr q, cb3d_histptr h,
                        ErrorMsg allocation_error);
local int cb3d_add_body(struct global_data* gd,
                        bodyptr p, bodyptr q, cb3d_histptr h,
                        ErrorMsg allocation_error);
local int cb3d_prepare_same_los_correction(struct global_data* gd,
                                           bodyptr p, cb3d_histptr h,
                                           ErrorMsg allocation_error);
local void cb3d_accumulate_pivot(cb3d_histptr h);
local int cb3d_get_modes(struct cmdline_data* cmd,
                         int *compute_xi, int *compute_zeta);
local int cb3d_exclude_same_los_option(struct cmdline_data* cmd);
local int cb3d_exclude_all_same_los_option(struct cmdline_data* cmd);
local void cb3d_init_ylm_norm(int lmax, REAL *norm);
local void cb3d_init_ylm_recurrence(cb3d_histptr h);
local void cb3d_accumulate_ylm_cartesian(cb3d_histptr h, int b, REAL field,
                                         double xhat, double yhat,
                                         double zhat);
local void cb3d_accumulate_ylm_to(cb3d_histptr h, REAL *alm_re, REAL *alm_im,
                                  int b, REAL field, double xhat,
                                  double yhat, double zhat);
local int cb3d_print_zeta(struct cmdline_data* cmd, struct global_data* gd,
                          REAL *num, REAL *den, INTEGER npivots, INTEGER nbb);
local int cb3d_print_xi(struct cmdline_data* cmd, struct global_data* gd,
                        REAL *num, REAL *den, INTEGER npivots, INTEGER nbb);
local REAL cb3d_bin_center(struct cmdline_data* cmd, struct global_data* gd,
                           int bzero);
local void cb3d_add_arrays(REAL *dst, const REAL *src, size_t n);
local void cb3d_zero_array(REAL *dst, size_t n);
local int cb3d_init_result(struct cmdline_data* cmd, cb3d_resultptr result,
                           int compute_xi, int compute_zeta, int lmax);
local void cb3d_free_result(cb3d_resultptr result);
local int cb3d_measure(struct cmdline_data* cmd, struct global_data* gd,
                       bodyptr *btable, INTEGER *nbody,
                       INTEGER ipmin, INTEGER ipmax, int cat,
                       int compute_xi, int compute_zeta,
                       int exclude_same_los, int exclude_all_same_los,
                       cb3d_resultptr result);
local int cb3d_validate(struct cmdline_data* cmd,
                        int *compute_xi, int *compute_zeta,
                        int survey_mode);
local int cb3d_print_survey_zeta(struct cmdline_data* cmd,
                                 struct global_data* gd,
                                 const cb3d_resultptr numerator,
                                 const cb3d_resultptr randoms,
                                 REAL random_scale, const cb3d_correction *corrected);
local int cb3d_print_survey_xi(struct cmdline_data* cmd,
                               struct global_data* gd,
                               const cb3d_resultptr numerator,
                               const cb3d_resultptr randoms,
                               REAL random_scale, const cb3d_correction *corrected);
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
    if (cb3d_init_result(cmd, &result, compute_xi, compute_zeta, cmd->mChebyshev) == FAILURE)
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
                          cb3d_exclude_same_los_option(cmd),
                          cb3d_exclude_all_same_los_option(cmd), &result);
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
    if(status==SUCCESS && cb3d_parallel_publish(cmd)) {
        size_t pair_shape[1]={(size_t)cmd->sizeHistN};
        size_t triple_shape[3]={(size_t)cmd->mChebyshev+1,(size_t)cmd->sizeHistN,(size_t)cmd->sizeHistN};
        cballs_result_adopt("pair_numerator",(void **)&result.xi_num,0,1,pair_shape);
        cballs_result_adopt("pair_denominator",(void **)&result.xi_den,0,1,pair_shape);
        cballs_result_adopt("triple_numerator",(void **)&result.num,0,3,triple_shape);
        cballs_result_adopt("triple_denominator",(void **)&result.den,0,3,triple_shape);
    }
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
    cb3d_correction corrected = {0};
    int compute_xi = FALSE;
    int compute_zeta = FALSE;
    int status = FAILURE;

    if (cb3d_validate(cmd, &compute_xi, &compute_zeta, TRUE) == FAILURE)
        goto setup_done;
    if (cb3d_init_result(cmd, &numerator, compute_xi, compute_zeta, cb_survey_solve_lmax(cmd)) == FAILURE)
        goto setup_done;
    if (cb3d_init_result(cmd, &randoms, compute_xi, compute_zeta, cb_survey_window_lmax(cmd)) == FAILURE)
        goto setup_done;
    status = SUCCESS;

setup_done:
    status = cb3d_parallel_consensus(cmd, status, "3D survey histogram setup");
    if (status == FAILURE) goto cleanup;
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
        "\n%s: ENCORE-style survey estimator; measuring D-R and random "
        "signal solve through ell=%d, window through ell=%d\n", cmd->searchMethod, numerator.lmax, randoms.lmax);
    status = cb3d_measure(cmd, gd, btable, nbody, 1, nbody[dmr_cat], dmr_cat,
                          compute_xi, compute_zeta, FALSE, FALSE, &numerator);
    if (status == FAILURE) goto cleanup;
    status = cb3d_measure(cmd, gd, btable, nbody, 1,
                          nbody[normalized_random_cat], normalized_random_cat,
                          compute_xi, compute_zeta, FALSE, FALSE, &randoms);
    if (status == FAILURE) goto cleanup;

    /* Correction is computation: retain it even when serialization is disabled. */
    if (cb3d_parallel_publish(cmd))
        status=cb3d_correct_survey(cmd,&numerator,&randoms,&corrected);
    status=cb3d_parallel_consensus(cmd,status,"3D survey correction");
    if(status==FAILURE) goto cleanup;
    /* Serialize only already computed products. */
    if (cb3d_parallel_publish(cmd) && !cballs_opt_no_out_hist(cmd)) {
        if (compute_zeta)
            status = cb3d_print_survey_zeta(cmd, gd, &numerator, &randoms,
                                            random_scale, &corrected);
        if (status == SUCCESS && compute_xi)
            status = cb3d_print_survey_xi(cmd, gd, &numerator, &randoms,
                                          random_scale, &corrected);
    }
    status = cb3d_parallel_consensus(cmd, status, "3D survey edge correction/output");

cleanup:
    gd->cpusearch = CPUTIME - cpustart;
    if(status==SUCCESS && cb3d_parallel_publish(cmd)) {
        size_t pair_shape[1]={(size_t)cmd->sizeHistN};
        size_t triple_shape[3]={(size_t)numerator.lmax+1,(size_t)cmd->sizeHistN,(size_t)cmd->sizeHistN};
        size_t random_shape[3]={(size_t)randoms.lmax+1,(size_t)cmd->sizeHistN,(size_t)cmd->sizeHistN};
        size_t corrected_shape[3]={(size_t)cb_survey_signal_lmax(cmd)+1,(size_t)cmd->sizeHistN,(size_t)cmd->sizeHistN};
        size_t shell_shape[2]={(size_t)cmd->sizeHistN,(size_t)cmd->sizeHistN};
        cballs_result_adopt("data_pair_numerator",(void **)&numerator.xi_num,0,1,pair_shape);
        cballs_result_adopt("data_pair_denominator",(void **)&numerator.xi_den,0,1,pair_shape);
        cballs_result_adopt("data_triple_numerator",(void **)&numerator.num,0,3,triple_shape);
        cballs_result_adopt("data_triple_denominator",(void **)&numerator.den,0,3,triple_shape);
        cballs_result_adopt("random_pair_numerator",(void **)&randoms.xi_num,0,1,pair_shape);
        cballs_result_adopt("random_pair_denominator",(void **)&randoms.xi_den,0,1,pair_shape);
        cballs_result_adopt("random_triple_numerator",(void **)&randoms.num,0,3,random_shape);
        cballs_result_adopt("random_triple_denominator",(void **)&randoms.den,0,3,random_shape);
#define CB3D_ADOPT(name,field,rank,shape) cballs_result_adopt(name,(void **)&corrected.field,0,rank,shape)
        CB3D_ADOPT("zeta_corrected_legendre",zeta,3,corrected_shape);
        CB3D_ADOPT("zeta_corrected_encore",encore,3,corrected_shape);
        CB3D_ADOPT("zeta_valid",valid,2,shell_shape);
        CB3D_ADOPT("zeta_status",status,2,shell_shape);
        CB3D_ADOPT("zeta_rcond",rcond,2,shell_shape);
        CB3D_ADOPT("zeta_relative_residual",residual,2,shell_shape);
        CB3D_ADOPT("xi_corrected",xi,1,pair_shape);
        CB3D_ADOPT("xi_valid",xi_valid,1,pair_shape);
        CB3D_ADOPT("xi_status",xi_status,1,pair_shape);
#undef CB3D_ADOPT
    }
    cb3d_free_correction(&corrected);
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
    if (survey_mode) {
        /* Bound each input before addition/multiplication in the order helpers. */
        if(cmd->surveySignalLmax < -1 || cmd->surveySignalLmax > CB_SURVEY_MAX_SOLVE
            || cmd->surveySolvePadding < 0 || cmd->surveySolvePadding > CB_SURVEY_MAX_SOLVE
            || cmd->surveyWindowLmax < -1 || cmd->surveyWindowLmax > CB_SURVEY_MAX_WINDOW
            || cb_survey_signal_lmax(cmd)<0 || cb_survey_solve_lmax(cmd)>CB_SURVEY_MAX_SOLVE
            || cb_survey_window_lmax(cmd)<2*cb_survey_solve_lmax(cmd)
            || !isfinite(cmd->surveyRcondMin) || cmd->surveyRcondMin<=0. || cmd->surveyRcondMin>1.) {
            snprintf(cmd->error_message,_ERRORMSGSIZE_,
                "physical survey requires signal>=0, padding>=0, signal+padding<=32, "
                "2*(signal+padding)<=window<=64 and 0<surveyRcondMin<=1; -1 selects automatic orders");
            return FAILURE;
        }
    }
    if (!isfinite((double)cmd->rangeN) || cmd->rangeN <= 0.0
        || !isfinite((double)cmd->rminHist) || cmd->rminHist < 0.0
        || cmd->rminHist >= cmd->rangeN) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: require finite 0<=rminHist<rangeN", routine_name);
        return FAILURE;
    }

    nb = (size_t)cmd->sizeHistN;
    lcount = (size_t)(survey_mode ? cb_survey_window_lmax(cmd) : cmd->mChebyshev) + 1;
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
                           int compute_xi, int compute_zeta, int lmax)
{
    size_t nb = (size_t)cmd->sizeHistN;
    size_t nacc = ((size_t)lmax + 1)*nb*nb;

    memset(result, 0, sizeof(*result));
    result->lmax=lmax;
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
                       int exclude_same_los, int exclude_all_same_los,
                       cb3d_resultptr result)
{
    const size_t nb = (size_t)cmd->sizeHistN;
    const size_t nacc = ((size_t)result->lmax + 1)*nb*nb;
    const INTEGER first_block = cb3d_parallel_first(cmd);
    const INTEGER block_stride = cb3d_parallel_stride(cmd);
    int allocation_failed = FALSE;
    int status = SUCCESS;
    ErrorMsg allocation_error = "";
#ifdef BALLS4SCANLEV
    bodyptr *pivot_order = NULL;
    INTEGER ordered_pivot_count = 0;
#endif

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
#ifdef BALLS4SCANLEV
    status = cballs_octree_scan_body_order(
        cmd, gd, btable[cat] + ipmin - 1, btable[cat] + ipmax,
        cat, &pivot_order, &ordered_pivot_count);
    status = cb3d_parallel_consensus(
        cmd, status, "3D scan-level body frontier");
    if (status == FAILURE) {
        free(pivot_order);
        return FAILURE;
    }
    verb_print(cmd->verbose,
               "octree-3pcf-3d: BALLS4SCANLEV spatially ordered %"
               INTEGER_FMT " exact body pivots\n",
               ordered_pivot_count);
    const INTEGER pivot_count = ordered_pivot_count;
#else
    const INTEGER pivot_count = ipmax - ipmin + 1;
#endif
    const INTEGER block_count = 1 + (pivot_count - 1)/CB3D_OMP_PIVOT_BLOCK_SIZE;

#ifdef BALLS4SCANLEV
#define CB3D_SCAN_SHARED ,pivot_order
#else
#define CB3D_SCAN_SHARED
#endif
#pragma omp parallel default(none) \
    shared(cmd,gd,btable,roottable,ipmin,ipmax,cat,result,nacc,nb, \
           compute_xi,compute_zeta,exclude_same_los,exclude_all_same_los, \
           allocation_failed, \
           allocation_error,first_block,block_stride,block_count,pivot_count \
           CB3D_SCAN_SHARED)
    {
        cb3d_hist hist;
        ErrorMsg thread_error = "";
        int hist_ready = cb3d_init_hist(cmd, gd, &hist, compute_xi, compute_zeta,
                                        exclude_same_los, exclude_all_same_los, result->lmax,
                                        thread_error) == SUCCESS;
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
#ifdef BALLS4SCANLEV
                const INTEGER first = block*CB3D_OMP_PIVOT_BLOCK_SIZE;
#else
                const INTEGER first = ipmin - 1 + block*CB3D_OMP_PIVOT_BLOCK_SIZE;
#endif
                const INTEGER count = MIN((INTEGER)CB3D_OMP_PIVOT_BLOCK_SIZE,
                                           pivot_count - block*CB3D_OMP_PIVOT_BLOCK_SIZE);
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
#ifdef BALLS4SCANLEV
                    bodyptr p = pivot_order[first + offset];
#else
                    bodyptr p = btable[cat] + first + offset;
#endif
                    if (Update(p) == FALSE) continue;
                    if (hist.read_mask && Mask(p) == FALSE) continue;
                    cb3d_clear_pivot(&hist);
                    hist.pivot_weight = Weight(p);
                    hist.pivot_field = cballs_product(Weight(p),Kappa(p));
                    if (hist.exclude_same_los)
                        hist.pivot_los = Octree3pcf3dLosId(p);
                    if (cb3d_walktree(gd, p, (nodeptr)roottable[cat], &hist,
                                      thread_error) == FAILURE
                        || cb3d_prepare_same_los_correction(
                               gd, p, &hist, thread_error) == FAILURE) {
#pragma omp critical(cb3d_allocation_failure)
                        {
                            if (!allocation_failed)
                                snprintf(allocation_error,
                                         sizeof(allocation_error), "%s",
                                         thread_error);
                            allocation_failed = TRUE;
                        }
                        break;
                    }
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
#undef CB3D_SCAN_SHARED
#ifdef BALLS4SCANLEV
    free(pivot_order);
#endif
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
    status=SUCCESS;
    if (cb3d_parallel_publish(cmd)) {
        if (compute_zeta) for(size_t k=0;k<nacc;k++)
            if (!isfinite(result->num[k]) || !isfinite(result->den[k])
                || !isfinite(cballs_checked_ratio(result->num[k],result->den[k]))) status=FAILURE;
        if (compute_xi) for(size_t k=0;k<nb;k++)
            if (!isfinite(result->xi_num[k]) || !isfinite(result->xi_den[k])
                || !isfinite(cballs_checked_ratio(result->xi_num[k],result->xi_den[k]))) status=FAILURE;
        if (status==FAILURE) snprintf(cmd->error_message,_ERRORMSGSIZE_,
            "numerical range: raw physical numerator/denominator or normalized result is not representable (overflow or underflow)");
    }
    status=cb3d_parallel_consensus(cmd,status,"physical raw arithmetic");
    if (status==FAILURE) return FAILURE;
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
        || scanopt(cmd->options, "exclude-pivot-los")
        || cb3d_exclude_all_same_los_option(cmd);
}

local int cb3d_exclude_all_same_los_option(struct cmdline_data* cmd)
{
    return scanopt(cmd->options, "exclude-all-same-los")
        || scanopt(cmd->options, "lya-distinct-forests");
}

local int cb3d_init_hist(struct cmdline_data* cmd, struct global_data* gd,
                         cb3d_histptr h,
                         int compute_xi, int compute_zeta,
                         int exclude_same_los, int exclude_all_same_los, int lmax,
                         ErrorMsg allocation_error)
{
    int nb = cmd->sizeHistN;
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
    h->exclude_all_same_los = exclude_all_same_los;
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
        CB3D_ALLOC(alm_distinct, nalm);
        CB3D_ALLOC(alm_im, nalm);
        CB3D_ALLOC(ylm_norm, (size_t)nlm);
        CB3D_ALLOC(ylm_a, (size_t)nlm);
        CB3D_ALLOC(ylm_b, (size_t)nlm);
        CB3D_ALLOC(ylm_diag, (size_t)lmax + 1);
        CB3D_ALLOC(shell_w, (size_t)nb);
        CB3D_ALLOC(shell_weight_distinct, (size_t)nb);
        if (exclude_all_same_los) {
            CB3D_ALLOC(distinct_los_power, nacc);
            CB3D_ALLOC(distinct_los_weight, (size_t)nb*(size_t)nb);
            CB3D_ALLOC(group_alm_re, nalm);
            CB3D_ALLOC(group_alm_im, nalm);
            CB3D_ALLOC(group_shell_w, (size_t)nb);
            CB3D_ALLOC(group_active_bins, (size_t)nb);
            CB3D_ALLOC(group_bin_active, (size_t)nb);
            memset(h->group_bin_active, 0,
                   (size_t)nb*sizeof(*h->group_bin_active));
        }
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
    free(h->alm_distinct);
    free(h->alm_im);
    free(h->ylm_norm);
    free(h->ylm_a);
    free(h->ylm_b);
    free(h->ylm_diag);
    free(h->shell_w);
    free(h->shell_weight_distinct);
    free(h->distinct_los_power);
    free(h->distinct_los_weight);
    free(h->group_alm_re);
    free(h->group_alm_im);
    free(h->group_shell_w);
    free(h->num);
    free(h->den);
    free(h->xi_num);
    free(h->xi_den);
    free(h->active_bins);
    free(h->group_active_bins);
    free(h->bin_active);
    free(h->group_bin_active);
    free(h->neighbors);
    memset(h, 0, sizeof(*h));
}

local void cb3d_clear_pivot(cb3d_histptr h)
{
    int i;
    if (!h->compute_zeta) return;
    for (i = 0; i < h->nactive; i++)
        h->bin_active[h->active_bins[i]] = 0;
    h->nactive = 0;
    h->neighbor_count = 0;
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
    memset(h->alm_distinct + offset, 0, (size_t)h->nlm * sizeof(*h->alm_distinct));
    h->shell_weight_distinct[b] = 0.0;
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

local int cb3d_walktree(struct global_data* gd,
                        bodyptr p, nodeptr q, cb3d_histptr h,
                        ErrorMsg allocation_error)
{
    nodeptr child;

    if (((nodeptr)p) == q) return SUCCESS;
    if (Type(q) == CELL) {
        if ((h->read_mask && Mask(q) == MASK_NODE_MASKED)
            || cb3d_reject_cell(gd, p, q, h))
            return SUCCESS;
        for (child = More(q); child != Next(q); child = Next(child)) {
            if (cb3d_walktree(gd, p, child, h, allocation_error) == FAILURE)
                return FAILURE;
        }
    } else {
        return cb3d_add_body(gd, p, (bodyptr)q, h, allocation_error);
    }
    return SUCCESS;
}

local int cb3d_add_body(struct global_data* gd,
                        bodyptr p, bodyptr q, cb3d_histptr h,
                        ErrorMsg allocation_error)
{
    REAL dr2;
    REAL dr1;
    compute_vector dr;
    int b;
    REAL f, w;

    if (p == q) return SUCCESS;
    if (h->read_mask && Mask(q) == FALSE) return SUCCESS;
    if (Update(q) == FALSE) return SUCCESS;
    if (h->exclude_same_los
        && h->pivot_los == Octree3pcf3dLosId(q)) return SUCCESS;

    DOTPSUBV(dr2, dr, Pos(p), Pos(q));
    if (h->use_periodic) {
        VWrapAll(dr);
        DOTVP(dr2, dr, dr);
    }
    if (dr2 <= h->rmin2 || dr2 >= h->rcut2) return SUCCESS;
    dr1 = rsqrt(dr2);
    b = cb3d_bin_index(h, dr1);
    if (b < 0) return SUCCESS;

    w = Weight(q);
    f = cballs_product(w,Kappa(q));

    if (h->compute_xi) {
        h->xi_num[b] += cballs_product(h->pivot_field,f);
        h->xi_den[b] += cballs_product(h->pivot_weight,w);
    }

    if (h->compute_zeta) {
        double invr = 1.0 / (double)dr1;
        double xhat = (double)dr[0] * invr;
        double yhat = (double)dr[1] * invr;
        double zhat = (double)dr[2] * invr;

        cb3d_activate_bin(h, b);
        cb3d_accumulate_ylm_cartesian(h, b, f, xhat, yhat, zhat);
        h->shell_weight_distinct[b] = cballs_distinct_merge(h->shell_weight_distinct[b],0.0,h->shell_w[b],w);
        h->shell_w[b] += w;
        if (h->exclude_all_same_los) {
            if (h->neighbor_count == h->neighbor_capacity) {
                size_t capacity = h->neighbor_capacity == 0
                    ? 256 : 2*h->neighbor_capacity;
                cb3d_neighbor *replacement;

                if (capacity < h->neighbor_capacity
                    || capacity > SIZE_MAX/sizeof(*replacement)) {
                    snprintf(allocation_error, _ERRORMSGSIZE_,
                             "3D LOS neighbor workspace size overflow");
                    return FAILURE;
                }
                replacement = (cb3d_neighbor *)realloc(
                    h->neighbors, capacity*sizeof(*replacement));
                if (replacement == NULL) {
                    snprintf(allocation_error, _ERRORMSGSIZE_,
                             "not enough memory for 3D LOS neighbor workspace");
                    return FAILURE;
                }
                h->neighbors = replacement;
                h->neighbor_capacity = capacity;
            }
            h->neighbors[h->neighbor_count].body = q;
            h->neighbors[h->neighbor_count].bin = b;
            h->neighbor_count++;
        }
    }

    h->nbbcalcthread++;
    return SUCCESS;
}

local int cb3d_compare_neighbor_los(const void *left, const void *right)
{
    const cb3d_neighbor *a = (const cb3d_neighbor *)left;
    const cb3d_neighbor *b = (const cb3d_neighbor *)right;
    INTEGER a_los = Octree3pcf3dLosId(a->body);
    INTEGER b_los = Octree3pcf3dLosId(b->body);

    return a_los < b_los ? -1 : (a_los > b_los ? 1 : 0);
}

local int cb3d_prepare_same_los_correction(struct global_data* gd,
                                           bodyptr p, cb3d_histptr h,
                                           ErrorMsg allocation_error)
{
    size_t begin;
    int i1_active, i2_active, ell, m;

    (void)allocation_error;
    if (!h->compute_zeta || !h->exclude_all_same_los) return SUCCESS;

    for (i1_active = 0; i1_active < h->nactive; i1_active++) {
        int b1 = h->active_bins[i1_active];
        for (i2_active = 0; i2_active < h->nactive; i2_active++) {
            int b2 = h->active_bins[i2_active];
            h->distinct_los_weight[(size_t)b1*(size_t)h->nbins + (size_t)b2] = 0.0;
            for (ell = 0; ell <= h->lmax; ell++)
                h->distinct_los_power[CB3D_ACC_INDEX(ell, b1, b2, h->nbins)] = 0.0;
        }
    }
    if (h->neighbor_count < 2) return SUCCESS;

    /* Reuse the main moments as a prefix of completed LOS groups. */
    memset(h->alm_re,0,(size_t)h->nbins*h->nlm*sizeof(*h->alm_re));
    memset(h->alm_im,0,(size_t)h->nbins*h->nlm*sizeof(*h->alm_im));
    memset(h->shell_w,0,(size_t)h->nbins*sizeof(*h->shell_w));
    qsort(h->neighbors, h->neighbor_count, sizeof(*h->neighbors),
          cb3d_compare_neighbor_los);
    begin = 0;
    while (begin < h->neighbor_count) {
        size_t end = begin + 1;
        size_t item;

        while (end < h->neighbor_count
               && Octree3pcf3dLosId(h->neighbors[end].body)
                  == Octree3pcf3dLosId(h->neighbors[begin].body))
            end++;
        h->group_nactive = 0;
        memset(h->group_alm_re,0,(size_t)h->nbins*h->nlm*sizeof(*h->group_alm_re));
        memset(h->group_alm_im,0,(size_t)h->nbins*h->nlm*sizeof(*h->group_alm_im));
        memset(h->group_shell_w,0,(size_t)h->nbins*sizeof(*h->group_shell_w));
        for (item = begin; item < end; item++) {
            bodyptr q = h->neighbors[item].body;
            int b = h->neighbors[item].bin;
            REAL dr2, radius;
            REAL field = cballs_product(Weight(q),Kappa(q));
            compute_vector dr;

            if (!h->group_bin_active[b]) {
                size_t offset = (size_t)b*(size_t)h->nlm;
                h->group_bin_active[b] = 1;
                h->group_active_bins[h->group_nactive++] = b;
                memset(h->group_alm_re + offset, 0,
                       (size_t)h->nlm*sizeof(*h->group_alm_re));
                memset(h->group_alm_im + offset, 0,
                       (size_t)h->nlm*sizeof(*h->group_alm_im));
                h->group_shell_w[b] = 0.0;
            }
            DOTPSUBV(dr2, dr, Pos(p), Pos(q));
            if (h->use_periodic) {
                VWrapAll(dr);
                DOTVP(dr2, dr, dr);
            }
            radius = rsqrt(dr2);
            cb3d_accumulate_ylm_to(
                h, h->group_alm_re, h->group_alm_im, b, field,
                (double)dr[0]/(double)radius,
                (double)dr[1]/(double)radius,
                (double)dr[2]/(double)radius);
            h->group_shell_w[b] += Weight(q);
        }

        for (i1_active = 0; i1_active < h->nactive; i1_active++) {
            int b1 = h->active_bins[i1_active];
            for (i2_active = i1_active;
                 i2_active < h->nactive; i2_active++) {
                int b2 = h->active_bins[i2_active];
                REAL group_weight = cballs_product(h->shell_w[b1],h->group_shell_w[b2])
                                  + cballs_product(h->group_shell_w[b1],h->shell_w[b2]);

                h->distinct_los_weight[(size_t)b1*(size_t)h->nbins + (size_t)b2]
                    += group_weight;
                if (b1 != b2)
                    h->distinct_los_weight[(size_t)b2*(size_t)h->nbins + (size_t)b1]
                        += group_weight;
                for (ell = 0; ell <= h->lmax; ell++) {
                    REAL power = 0.0;
                    for (m = 0; m <= ell; m++) {
                        int n = CB3D_LM_INDEX(ell, m);
                        size_t i1 = CB3D_ALM_INDEX(b1, h->nlm, n);
                        size_t i2 = CB3D_ALM_INDEX(b2, h->nlm, n);
                        REAL product = cballs_product(h->alm_re[i1],h->group_alm_re[i2])
                                     + cballs_product(h->alm_im[i1],h->group_alm_im[i2])
                                     + cballs_product(h->group_alm_re[i1],h->alm_re[i2])
                                     + cballs_product(h->group_alm_im[i1],h->alm_im[i2]);
                        power += m == 0 ? product : 2.0*product;
                    }
                    h->distinct_los_power[
                        CB3D_ACC_INDEX(ell, b1, b2, h->nbins)] += power;
                    if (b1 != b2)
                        h->distinct_los_power[
                            CB3D_ACC_INDEX(ell, b2, b1, h->nbins)] += power;
                }
            }
        }
        for (i1_active = 0; i1_active < h->group_nactive; i1_active++) {
            int b=h->group_active_bins[i1_active];
            h->shell_w[b]+=h->group_shell_w[b];
            for (int n=0;n<h->nlm;n++) {
                size_t i=CB3D_ALM_INDEX(b,h->nlm,n);
                h->alm_re[i]+=h->group_alm_re[i]; h->alm_im[i]+=h->group_alm_im[i];
            }
            h->group_bin_active[b]=0;
        }
        begin = end;
    }
    return SUCCESS;
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

                size_t idx = CB3D_ACC_INDEX(ell, b1, b2, h->nbins);
                REAL num_increment;
                REAL den_increment;

                for (m = 0; m <= ell; m++) {
                    int n = CB3D_LM_INDEX(ell, m);
                    size_t i1 = CB3D_ALM_INDEX(b1, h->nlm, n);
                    size_t i2 = CB3D_ALM_INDEX(b2, h->nlm, n);
                    REAL reprod = b1==b2 ? 0.0 : cballs_product(h->alm_re[i1],h->alm_re[i2])
                                + cballs_product(h->alm_im[i1],h->alm_im[i2]);
                    if (b1 == b2 && !h->exclude_all_same_los)
                        reprod = h->alm_distinct[i1];
                    power += (m == 0) ? reprod : 2.0*reprod;
                }

                d = b1==b2 ? 0.0 : cballs_product(h->shell_w[b1],h->shell_w[b2]);
                if (h->exclude_all_same_los) {
                    power = h->distinct_los_power[idx];
                    d = h->distinct_los_weight[
                        (size_t)b1*(size_t)h->nbins + (size_t)b2];
                } else if (b1 == b2) {
                    d = h->shell_weight_distinct[b1];
                }

                num_increment = cballs_product(wp_delta,
                              (REAL)CB3D_FOURPI * power);
                den_increment = cballs_product(wp_count, d);
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
    cb3d_accumulate_ylm_to(h, h->alm_re, h->alm_im, radial_bin, field,
                           xhat, yhat, zhat);
}

local inline void cb3d_accumulate_mode(cb3d_histptr h, REAL *alm_re, REAL *alm_im,
                                        size_t index, double field, double yr, double yi)
{
    const REAL re = cballs_product(field,yr);
    const REAL im = -cballs_product(field,yi);
    if (alm_re == h->alm_re)
        h->alm_distinct[index] += cballs_distinct_complex_cross(
            alm_re[index],alm_im[index],re,im);
    alm_re[index] += re;
    alm_im[index] += im;
}

local void cb3d_accumulate_ylm_to(cb3d_histptr h, REAL *alm_re, REAL *alm_im,
                                  int radial_bin, REAL field, double xhat,
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
        cb3d_accumulate_mode(h,alm_re,alm_im,index,weighted_field,
                                ymm_re,ymm_im);

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
            cb3d_accumulate_mode(h,alm_re,alm_im,index,weighted_field,
                                y_lm1_re,y_lm1_im);

            for (ell = m + 2; ell <= h->lmax; ell++) {
                double ylm_re;
                double ylm_im;

                n = CB3D_LM_INDEX(ell, m);
                ylm_re = (double)h->ylm_a[n] * zhat * y_lm1_re
                       - (double)h->ylm_b[n] * y_lm2_re;
                ylm_im = (double)h->ylm_a[n] * zhat * y_lm1_im
                       - (double)h->ylm_b[n] * y_lm2_im;
                index = offset + (size_t)n;
                cb3d_accumulate_mode(h,alm_re,alm_im,index,weighted_field,
                                ylm_re,ylm_im);
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

local void cb3d_free_correction(cb3d_correction *c)
{
    free(c->zeta);free(c->encore);free(c->valid);free(c->status);
    free(c->rcond);free(c->residual);free(c->xi);free(c->xi_valid);free(c->xi_status);
    memset(c,0,sizeof(*c));
}

local int cb3d_correct_survey(struct cmdline_data *cmd,
    const cb3d_resultptr numerator,const cb3d_resultptr randoms,cb3d_correction *c)
{
    const int nb=cmd->sizeHistN, signal=cb_survey_signal_lmax(cmd);
    const size_t shells=(size_t)nb*nb;
#define CB3D_CORR_ALLOC(field,count) \
    if(cballs_calloc_checked((void **)&c->field,count,sizeof(REAL), \
        "physical survey " #field,cmd->error_message,_ERRORMSGSIZE_)==FAILURE) goto fail
    if(numerator->num) {
        CB3D_CORR_ALLOC(zeta,((size_t)signal+1)*shells);
        CB3D_CORR_ALLOC(encore,((size_t)signal+1)*shells);
        CB3D_CORR_ALLOC(valid,shells);CB3D_CORR_ALLOC(status,shells);
        CB3D_CORR_ALLOC(rcond,shells);CB3D_CORR_ALLOC(residual,shells);
        for(int b1=0;b1<nb;b1++) for(int b2=b1;b2<nb;b2++) {
            double num[CB_SURVEY_MAX_SOLVE+1],ran[CB_SURVEY_MAX_WINDOW+1];
            double x[CB_SURVEY_MAX_SOLVE+1],rcond,residual;
            for(int ell=0;ell<=numerator->lmax;ell++)
                num[ell]=numerator->num[CB3D_ACC_INDEX(ell,b1,b2,nb)];
            for(int ell=0;ell<=randoms->lmax;ell++)
                ran[ell]=randoms->num[CB3D_ACC_INDEX(ell,b1,b2,nb)];
            int status=cb_survey_solve(num,ran,numerator->lmax,randoms->lmax,
                                      cmd->surveyRcondMin,x,&rcond,&residual);
            for(int order=0;order<2;order++) {
                int left=order?b2:b1,right=order?b1:b2;
                size_t pair=(size_t)left*nb+right;
                c->valid[pair]=status==CB_SURVEY_VALID;
                c->status[pair]=status;c->rcond[pair]=rcond;c->residual[pair]=residual;
                for(int ell=0;ell<=signal;ell++) {
                    size_t i=CB3D_ACC_INDEX(ell,left,right,nb);
                    c->zeta[i]=status==CB_SURVEY_VALID ? x[ell] : NAN;
                    c->encore[i]=c->zeta[i]*CB3D_FOURPI*((ell&1)?-1.:1.)/sqrt(2.*ell+1.);
                    if(status==CB_SURVEY_VALID && (!isfinite(c->encore[i])
                        || (c->zeta[i]!=0. && c->encore[i]==0.))) {
                        snprintf(cmd->error_message,_ERRORMSGSIZE_,"physical survey ENCORE conversion is not representable");
                        goto fail;
                    }
                }
            }
        }
    }
    if(numerator->xi_num) {
        CB3D_CORR_ALLOC(xi,(size_t)nb);CB3D_CORR_ALLOC(xi_valid,(size_t)nb);
        CB3D_CORR_ALLOC(xi_status,(size_t)nb);
        for(int b=0;b<nb;b++) {
            double den=randoms->xi_num[b], value=den>0. ? cballs_checked_ratio(numerator->xi_num[b],den) : NAN;
            int status=den==0. ? CB_SURVEY_EMPTY : isfinite(value) ? CB_SURVEY_VALID : CB_SURVEY_NONFINITE;
            c->xi[b]=status==CB_SURVEY_VALID ? value : NAN;
            c->xi_valid[b]=status==CB_SURVEY_VALID;c->xi_status[b]=status;
        }
    }
#undef CB3D_CORR_ALLOC
    return SUCCESS;
fail:
    cb3d_free_correction(c);return FAILURE;
}

local int cb3d_print_survey_zeta(struct cmdline_data* cmd,
    struct global_data* gd,const cb3d_resultptr numerator,
    const cb3d_resultptr randoms,REAL random_scale,const cb3d_correction *c)
{
    char namebuf[MAXLENGTHOFFILES];
    if(format_checked(namebuf,sizeof(namebuf),"3D survey 3PCF output path",
        "%s_3d_survey%s",gd->fpfnamehistZetaMFileName,EXTFILES)!=0) {
        snprintf(cmd->error_message,_ERRORMSGSIZE_,"3D survey 3PCF output path is too long");return FAILURE;
    }
    stream out=fopen(namebuf,"w");
    if(!out) {snprintf(cmd->error_message,_ERRORMSGSIZE_,"cannot open 3D survey output '%s': %s",namebuf,strerror(errno));return FAILURE;}
    fprintf(out,"# cBalls survey 3PCF; computed before serialization\n");
    fprintf(out,"# alpha %.17e signal_lmax %d solve_lmax %d window_lmax %d rcond_min %.17e\n",
        (double)random_scale,cb_survey_signal_lmax(cmd),numerator->lmax,randoms->lmax,(double)cmd->surveyRcondMin);
    fprintf(out,"# N=(D-alpha*R)^3; R=(alpha*R)^3; only bin1<bin2 emitted (API includes diagonal and transpose)\n");
    fprintf(out,"# zeta_legendre multiplies P_ell; zeta_encore multiplies (-1)^ell*sqrt(2ell+1)*P_ell/(4*pi)\n");
    fprintf(out,"# invalid coefficients are NaN; status 0=empty 1=valid 2=singular 3=ill-conditioned 4=residual 5=nonfinite\n");
    fprintf(out,"# condition is infinity-norm condition in standard Legendre basis; relative_residual is backward error\n");
    fprintf(out,"# columns: ell bin1 bin2 r1_center r2_center zeta_legendre zeta_encore N_legendre R_legendre N_encore R_encore R0_encore condition valid rcond relative_residual status\n");
    int nb=cmd->sizeHistN;
    for(int b1=0;b1<nb;b1++) for(int b2=b1+1;b2<nb;b2++) {
        size_t pair=(size_t)b1*nb+b2;
        for(int ell=0;ell<=cb_survey_signal_lmax(cmd);ell++) {
            size_t i=CB3D_ACC_INDEX(ell,b1,b2,nb);
            double conversion=CB3D_FOURPI*((ell&1)?-1.:1.)*sqrt(2.*ell+1.);
            fprintf(out,"%d %d %d %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %d %.17e %.17e %d\n",
                ell,b1+1,b2+1,(double)cb3d_bin_center(cmd,gd,b1),(double)cb3d_bin_center(cmd,gd,b2),
                (double)c->zeta[i],(double)c->encore[i],(double)numerator->num[i],(double)randoms->num[i],
                numerator->num[i]/conversion,randoms->num[i]/conversion,randoms->num[pair]/CB3D_FOURPI,
                c->rcond[pair]>0. ? 1./c->rcond[pair] : HUGE_VAL,(int)c->valid[pair],
                (double)c->rcond[pair],(double)c->residual[pair],(int)c->status[pair]);
        }
    }
    int failed=ferror(out);if(fclose(out)) failed=1;
    if(failed) {snprintf(cmd->error_message,_ERRORMSGSIZE_,"failed writing 3D survey output '%s'",namebuf);return FAILURE;}
    return SUCCESS;
}

local int cb3d_print_survey_xi(struct cmdline_data* cmd,
    struct global_data* gd,const cb3d_resultptr numerator,
    const cb3d_resultptr randoms,REAL random_scale,const cb3d_correction *c)
{
    char namebuf[MAXLENGTHOFFILES];
    if(format_checked(namebuf,sizeof(namebuf),"3D survey 2PCF output path",
        "%s_3d_survey%s",gd->fpfnamehistXi2pcfFileName,EXTFILES)!=0) {
        snprintf(cmd->error_message,_ERRORMSGSIZE_,"3D survey 2PCF output path is too long");return FAILURE;
    }
    stream out=fopen(namebuf,"w");
    if(!out) {snprintf(cmd->error_message,_ERRORMSGSIZE_,"cannot open 3D survey output '%s': %s",namebuf,strerror(errno));return FAILURE;}
    fprintf(out,"# cBalls survey 2PCF; xi=(D-alpha*R)^2/(alpha*R)^2; alpha %.17e\n",(double)random_scale);
    fprintf(out,"# invalid coefficients are NaN; status 0=empty 1=valid 5=nonfinite\n");
    fprintf(out,"# columns: bin r_center xi N_DminusR_squared RR valid status\n");
    for(int b=0;b<cmd->sizeHistN;b++)
        fprintf(out,"%d %.17e %.17e %.17e %.17e %d %d\n",b+1,(double)cb3d_bin_center(cmd,gd,b),
            (double)c->xi[b],(double)numerator->xi_num[b],(double)randoms->xi_num[b],(int)c->xi_valid[b],(int)c->xi_status[b]);
    int failed=ferror(out);if(fclose(out)) failed=1;
    if(failed) {snprintf(cmd->error_message,_ERRORMSGSIZE_,"failed writing 3D survey output '%s'",namebuf);return FAILURE;}
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
