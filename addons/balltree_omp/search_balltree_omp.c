/* ==============================================================================
 MODULE: search_balltree_omp.c        [cTreeBalls]
 Written by: M.A. Rodriguez-Meza.
 Starting date:    april 2023
 Purpose: 2/3-point correlation functions computation
 Language: C
 Use: kd = searchcalc_balltree_omp(cmd, gd, btab, nbody,
                                    ipmin, ipmax, cat1, cat2);
 Major revisions:
 ==============================================================================*/
//        1          2          3          4        ^ 5          6          7

#include <limits.h>

#include "globaldefs.h"

#include "fcfc_balltree.h"
#ifdef BALLTREEMPI
#include "fcfc_balltree_mpi.h"
#endif

static inline bool balltree_intersects(struct cmdline_data *cmd,
                                       struct global_data *gd, bodyptr p,
                                       const fcfc_ballnode *node)
{
#ifdef SINGLEP
    if (!cmd->usePeriodic) {
        const cballs_storage_real dx = Pos(p)[0] - node->center[0];
        const cballs_storage_real dy = Pos(p)[1] - node->center[1];
#if NDIM == 3
        const cballs_storage_real dz = Pos(p)[2] - node->center[2];
#endif
        const cballs_storage_real limit = cballs_store_search_bound(
            gd->Rcut + (real)node->radius);
        const cballs_storage_real distance2 = dx * dx + dy * dy
#if NDIM == 3
            + dz * dz
#endif
            ;

        return distance2 < limit * limit;
    }
#endif
    compute_vector dr;
    real distance2;
    const real limit = gd->Rcut + (real)node->radius;

    DOTPSUBV(distance2, dr, Pos(p), node->center);
    if (cmd->usePeriodic) {
        VWrapAll(dr);
        DOTVP(distance2, dr, dr);
    }
    return distance2 < limit * limit;
}

local void sumnode_sincos(struct  cmdline_data*, struct  global_data*,
                          bodyptr, fcfc_ballnode, fcfc_balltreeptr,
                          INTEGER *, INTEGER *,
                          gdhistptr_sincos_omp);
local void sumnode_sincos_cell(struct  cmdline_data*,
                               struct  global_data*, bodyptr,
                               fcfc_ballnode, bodyptr *,
                               INTEGER *, INTEGER *,
                               gdhistptr_sincos_omp);
local void walk_balltree_exact(struct cmdline_data*, struct global_data*,
                            bodyptr, fcfc_balltreeptr, INTEGER *, INTEGER *,
                            gdhistptr_sincos_omp);
local void walk_balltree_one_ball(struct cmdline_data*, struct global_data*,
                               bodyptr, fcfc_balltreeptr, INTEGER *, INTEGER *,
                               gdhistptr_sincos_omp);
local int print_info(struct cmdline_data* cmd,
                     struct  global_data* gd);

static inline bool balltree_accept_body(struct cmdline_data *cmd,
                                      struct global_data *gd,
                                      bodyptr p,
#ifdef SINGLEP
                                      const fcfc_ballpoint *q,
#else
                                      bodyptr q,
#endif
                                      real *distance, compute_vector dr)
{
    real distance_squared;

#ifdef SINGLEP
    DOTPSUBV(distance_squared, dr, Pos(p), q->pos);
#else
    DOTPSUBV(distance_squared, dr, Pos(p), Pos(q));
#endif
    if (cmd->usePeriodic) {
        VWrapAll(dr);
        DOTVP(distance_squared, dr, dr);
    }

    if (distance_squared >= gd->RcutSq)
        return FALSE;

    *distance = rsqrt(distance_squared);
    return *distance < gd->Rcut;
}

typedef struct {
    INTEGER nbbcalc;
    INTEGER nbccalc;
#ifdef SMOOTHPIVOT
    INTEGER ipfalse;
    INTEGER count_rmin;
    INTEGER count_overlap;
#endif
} balltree_worker_stats;

static void process_balltree_pivot(struct cmdline_data *cmd,
                                   struct global_data *gd, bodyptr p,
                                   bodyptr pivot_base, INTEGER pivot_count,
                                   fcfc_balltreeptr tree, bool use_one_ball,
                                   gdhistptr_sincos_omp hist,
                                   balltree_worker_stats *stats)
{
    int n;

#ifdef SMOOTHPIVOT
    if (Update(p) == FALSE) {
        stats->ipfalse++;
        return;
    }
#endif
    for (n = 1; n <= cmd->sizeHistN; n++) {
        hist->histNNSubthread[n] = 0.0;
        hist->histXi2pcfthreadsub[n] = 0.0;
#ifdef SMOOTHPIVOT
        hist->histNNSubXi2pcfthreadp[n] = 0.0;
#endif
    }
#ifdef TPCF
    CLRM_ext_ext(hist->histXithreadcos, cmd->mChebyshev + 1,
                 cmd->sizeHistN);
    CLRM_ext_ext(hist->histXithreadsin, cmd->mChebyshev + 1,
                 cmd->sizeHistN);
#if NDIM == 3
    dRotation3D(Pos(p), ROTANGLE, ROTANGLE, ROTANGLE, hist->q0);
    DOTPSUBV(hist->drpq2, hist->dr0, Pos(p), hist->q0);
    hist->drpq = rsqrt(hist->drpq2);
#ifdef PTOPIVOTROTATION
    real rtheta = xrandom(0.0, TWOPI);
    compute_vector dr0rot;
    RotationVecAWRtoVecB(dr0rot, hist->dr0, Pos(p), rtheta);
    SETV(hist->dr0, dr0rot);
#endif
#endif
#endif

    if (use_one_ball)
        walk_balltree_one_ball(cmd, gd, p, tree, &stats->nbbcalc,
                               &stats->nbccalc, hist);
    else
        walk_balltree_exact(cmd, gd, p, tree, &stats->nbbcalc,
                            &stats->nbccalc, hist);

#ifdef SMOOTHPIVOT
    for (n = 1; n <= cmd->sizeHistN; n++) {
        hist->histNNSubXi2pcfthreadp[n] =
            (real)NbRmin(p) * hist->histNNSubXi2pcfthreadp[n];
        hist->histNNSubXi2pcfthreadtotal[n] +=
            hist->histNNSubXi2pcfthreadp[n];
        hist->histNNSubthread[n] =
            (real)NbRmin(p) * hist->histNNSubthread[n];
    }
#endif

    computeBodyProperties_sincos(cmd, gd, p, (int)pivot_count, hist);
#ifdef SMOOTHPIVOT
    stats->count_rmin += NbRmin(p);
    stats->count_overlap += NbRminOverlap(p);
#endif

    const INTEGER ip = p - pivot_base + 1;
    if (ip % cmd->stepState == 0)
        verb_log_print(cmd->verbose_log, gd->outlog,
                       " - Completed pivot: %" INTEGER_FMT "\n", ip);
}

#ifdef BALLTREEMPI
static int reduce_balltree_histograms(struct cmdline_data *cmd,
                                      struct global_data *gd,
#ifdef SMOOTHPIVOT
                                      INTEGER *ipfalse,
                                      INTEGER *count_rmin,
                                      INTEGER *count_overlap
#else
                                      INTEGER *unused1,
                                      INTEGER *unused2,
                                      INTEGER *unused3
#endif
                                      )
{
    const size_t bins = (size_t)cmd->sizeHistN;
    size_t count;
    size_t cursor = 0;
    real *packed;
    INTEGER counters[5] = {gd->nbbcalc, gd->nbccalc, gd->ncccalc, 0, 0};

    if (bins > SIZE_MAX / 4) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "balltree-mpi histogram packing size overflow");
        return FAILURE;
    }
    count = 4 * bins;

#ifndef SMOOTHPIVOT
    (void)unused1;
    (void)unused2;
    (void)unused3;
#endif
#ifdef SMOOTHPIVOT
    if (bins > SIZE_MAX - count) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "balltree-mpi histogram packing size overflow");
        return FAILURE;
    }
    count += bins;
    counters[2] = *ipfalse;
    counters[3] = *count_rmin;
    counters[4] = *count_overlap;
#endif
#ifdef TPCF
    const size_t orders = (size_t)cmd->mChebyshev + 1;
    size_t plane;
    size_t zeta_count;
    if ((bins != 0 && bins > SIZE_MAX / bins)
        || (plane = bins * bins, orders != 0 && plane > SIZE_MAX / orders)
        || (zeta_count = orders * plane, zeta_count > (SIZE_MAX - count) / 4)) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "balltree-mpi histogram packing size overflow");
        return FAILURE;
    }
    count += 4 * zeta_count;
#endif
    if (count > SIZE_MAX / sizeof(*packed)) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "balltree-mpi histogram allocation size overflow");
        return FAILURE;
    }
    packed = malloc(count * sizeof(*packed));
    if (fcfc_balltree_mpi_consensus(cmd,
            packed == NULL ? FAILURE : SUCCESS,
            "MPI histogram packing allocation") == FAILURE) {
        free(packed);
        return FAILURE;
    }

    for (int n = 1; n <= cmd->sizeHistN; n++) packed[cursor++] = gd->histNN[n];
    for (int n = 1; n <= cmd->sizeHistN; n++) packed[cursor++] = gd->histNNSub[n];
    for (int n = 1; n <= cmd->sizeHistN; n++) packed[cursor++] = gd->histNNSubXi2pcf[n];
    for (int n = 1; n <= cmd->sizeHistN; n++) packed[cursor++] = gd->histXi2pcf[n];
#ifdef SMOOTHPIVOT
    for (int n = 1; n <= cmd->sizeHistN; n++)
        packed[cursor++] = gd->histNNSubXi2pcftotal[n];
#endif
#ifdef TPCF
#define PACK_ZETA(array)                                                   \
    do {                                                                  \
        for (int m = 1; m <= cmd->mChebyshev + 1; m++)                   \
            for (int n = 1; n <= cmd->sizeHistN; n++)                    \
                for (int l = 1; l <= cmd->sizeHistN; l++)                \
                    packed[cursor++] = (array)[m][n][l];                  \
    } while (0)
    PACK_ZETA(gd->histZetaMcos);
    PACK_ZETA(gd->histZetaMsin);
    PACK_ZETA(gd->histZetaMsincos);
    PACK_ZETA(gd->histZetaMcossin);
#undef PACK_ZETA
#endif
    if (cursor != count
        || fcfc_balltree_mpi_reduce_reals(cmd, packed, count) == FAILURE
        || fcfc_balltree_mpi_reduce_integers(cmd, counters, 5) == FAILURE) {
        free(packed);
        if (cursor != count)
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     "balltree-mpi internal histogram packing mismatch");
        return FAILURE;
    }

    if (fcfc_balltree_mpi_is_root()) {
        cursor = 0;
        for (int n = 1; n <= cmd->sizeHistN; n++) gd->histNN[n] = packed[cursor++];
        for (int n = 1; n <= cmd->sizeHistN; n++) gd->histNNSub[n] = packed[cursor++];
        for (int n = 1; n <= cmd->sizeHistN; n++) gd->histNNSubXi2pcf[n] = packed[cursor++];
        for (int n = 1; n <= cmd->sizeHistN; n++) gd->histXi2pcf[n] = packed[cursor++];
#ifdef SMOOTHPIVOT
        for (int n = 1; n <= cmd->sizeHistN; n++)
            gd->histNNSubXi2pcftotal[n] = packed[cursor++];
#endif
#ifdef TPCF
#define UNPACK_ZETA(array)                                                 \
    do {                                                                  \
        for (int m = 1; m <= cmd->mChebyshev + 1; m++)                   \
            for (int n = 1; n <= cmd->sizeHistN; n++)                    \
                for (int l = 1; l <= cmd->sizeHistN; l++)                \
                    (array)[m][n][l] = packed[cursor++];                  \
    } while (0)
        UNPACK_ZETA(gd->histZetaMcos);
        UNPACK_ZETA(gd->histZetaMsin);
        UNPACK_ZETA(gd->histZetaMsincos);
        UNPACK_ZETA(gd->histZetaMcossin);
#undef UNPACK_ZETA
#endif
        gd->nbbcalc = counters[0];
        gd->nbccalc = counters[1];
#ifdef SMOOTHPIVOT
        *ipfalse = counters[2];
        *count_rmin = counters[3];
        *count_overlap = counters[4];
#else
        gd->ncccalc = counters[2];
#endif
    }
    free(packed);
    return SUCCESS;
}
#endif

static int searchcalc_balltree_driver(struct cmdline_data *cmd,
                                      struct global_data *gd,
                                      bodyptr *btab, INTEGER *nbody,
                                      INTEGER ipmin, INTEGER *ipmax,
                                      int cat1, int cat2, bool distributed)
{
    bodyptr p;
    int n;
    int status = FAILURE;
    double cpustart = CPUTIME;
    fcfc_balltreeptr tree = NULL;
#ifdef BALLTREEMPI
    fcfc_balltreeptr pivot_tree = NULL;
    INTEGER *frontier = NULL;
    INTEGER frontier_count = 0;
    INTEGER task_first = 0;
    INTEGER task_last = 0;
    int scheduler_status = SUCCESS;
    fcfc_balltree_mpi_scheduler scheduler;
    memset(&scheduler, 0, sizeof(scheduler));
#else
    (void)distributed;
#endif
    const int nbucket = cmd->nsmooth;
    const bool use_one_ball = cballs_opt_behavior_ball(cmd)
        && !cballs_opt_no_one_ball(cmd);
    int allocation_failed = FALSE;
#ifdef SMOOTHPIVOT
    INTEGER ipfalse = 0;
    INTEGER icountNbRmin = 0;
    INTEGER icountNbRminOverlap = 0;
#endif

    if (print_info(cmd, gd) == FAILURE) return FAILURE;
#ifdef OPENMPCODE
    ThreadCount(cmd, gd, nbody[cat1], cat1);
#endif
    int hist_status = search_init_gd_hist_sincos(cmd, gd);
#ifdef BALLTREEMPI
    if (distributed)
        hist_status = fcfc_balltree_mpi_consensus(
            cmd, hist_status, "MPI histogram initialization");
#endif
    if (hist_status == FAILURE) return FAILURE;

#ifdef BALLTREEMPI
    if (distributed
        && fcfc_balltree_mpi_broadcast_catalogs(cmd, gd, btab, nbody,
                                                cat1, cat2) == FAILURE)
        goto cleanup;
#endif
    verb_print(cmd->verbose,
               "\nsearchcalc_balltree: Total allocated %g MByte storage so far.\n",
               gd->bytes_tot/(1024.0*1024.0));
    DO_BODY(p, btab[cat1] + ipmin - 1, btab[cat1] + ipmax[cat1])
        Update(p) = TRUE;

    const real cpu_build_balltree = CPUTIME;
    verb_print(cmd->verbose, "\nballtree build: nbucket = %d\n", nbucket);
#ifdef BALLTREEMPI
    if (distributed) {
        if (fcfc_balltree_mpi_build(cmd, gd, btab[cat2], nbody[cat2],
                                    nbucket, &tree) == FAILURE)
            goto cleanup;
    } else
#endif
    if (fcfc_balltree_build(cmd, gd, btab[cat2], nbody[cat2], nbucket,
                            &tree) == FAILURE)
        goto cleanup;
    verb_print(cmd->verbose, "balltree build: CPU time = %lf\n",
               CPUTIME - cpu_build_balltree);
    gd->ncellTable[cat2] = tree->nnode;

#ifdef SMOOTHPIVOT
    int smooth_status = prepare_smooth_pivots(cmd, gd, btab, nbody,
                                               ipmin, ipmax, cat1, cat2);
#ifdef BALLTREEMPI
    if (distributed)
        smooth_status = fcfc_balltree_mpi_consensus(
            cmd, smooth_status, "MPI smooth-pivot preparation");
#endif
    if (smooth_status == FAILURE) goto cleanup;
#endif

#ifdef BALLTREEMPI
    if (distributed) {
        int thread_hint = 1;
#ifdef OPENMPCODE
        thread_hint = omp_get_max_threads();
#endif
        if (cat1 == cat2)
            pivot_tree = tree;
        else if (fcfc_balltree_mpi_build(cmd, gd, btab[cat1], nbody[cat1],
                                         nbucket, &pivot_tree) == FAILURE)
            goto cleanup;
        uintmax_t target = (uintmax_t)fcfc_balltree_mpi_size()
            * (uintmax_t)thread_hint * 8U;
#ifdef LONGINT
        if (target > (uintmax_t)LONG_MAX) target = LONG_MAX;
#else
        if (target > (uintmax_t)INT_MAX) target = INT_MAX;
#endif
        if (fcfc_balltree_mpi_frontier(cmd, pivot_tree, (INTEGER)target,
                                       &frontier, &frontier_count) == FAILURE)
            goto cleanup;
        INTEGER step = thread_hint > 0 ? (INTEGER)thread_hint * 4 : 1;
        if (fcfc_balltree_mpi_scheduler_init(cmd, &scheduler,
                                              frontier_count, step) == FAILURE)
            goto cleanup;
        verb_print(cmd->verbose,
                   "balltree-mpi: %d ranks, %" INTEGER_FMT
                   " balanced frontier tasks\n",
                   fcfc_balltree_mpi_size(), frontier_count);
    }
#endif

#pragma omp parallel
    {
        int hn;
        balltree_worker_stats stats;
        memset(&stats, 0, sizeof(stats));
        gdhist_sincos_omp hist;
        int hist_ready = search_init_sincos_omp(cmd, gd, &hist) == SUCCESS;
        if (!hist_ready) {
#pragma omp atomic write
            allocation_failed = TRUE;
        }
#pragma omp barrier

#ifdef BALLTREEMPI
        if (distributed) {
#pragma omp master
            {
                if (fcfc_balltree_mpi_consensus(
                        cmd, allocation_failed ? FAILURE : SUCCESS,
                        "MPI OpenMP histogram allocation") == FAILURE)
                    allocation_failed = TRUE;
            }
#pragma omp barrier
            while (!allocation_failed) {
#pragma omp master
                {
                    scheduler_status = fcfc_balltree_mpi_scheduler_claim(
                        cmd, &scheduler, &task_first, &task_last);
                    if (scheduler_status == FAILURE)
                        allocation_failed = TRUE;
                }
#pragma omp barrier
                if (allocation_failed || task_first == task_last) break;
#pragma omp for schedule(dynamic,1)
                for (INTEGER itask = task_first; itask < task_last; itask++) {
                    const fcfc_ballnode *task =
                        &pivot_tree->nodes[frontier[itask]];
                    for (INTEGER ip = task->first; ip <= task->last; ip++) {
                        bodyptr pivot = pivot_tree->bptr[ip];
                        const INTEGER original = pivot - btab[cat1] + 1;
                        if (original >= ipmin && original <= ipmax[cat1])
                            process_balltree_pivot(cmd, gd, pivot, btab[cat1],
                                nbody[cat1], tree, use_one_ball, &hist, &stats);
                    }
                }
#pragma omp barrier
            }
        } else
#endif
        {
#pragma omp for nowait schedule(static,1)
            DO_BODY(p, btab[cat1] + ipmin - 1,
                    btab[cat1] + ipmax[cat1]) {
                if (!allocation_failed)
                    process_balltree_pivot(cmd, gd, p, btab[cat1],
                        nbody[cat1], tree, use_one_ball, &hist, &stats);
            }
        }

#ifdef OPENMPCODE
        const int thread_id = omp_get_thread_num();
        const int thread_count = omp_get_num_threads();
#else
        const int thread_id = 0;
        const int thread_count = 1;
#endif
        for (int thread_turn = 0; thread_turn < thread_count; thread_turn++) {
#pragma omp barrier
            if (thread_id == thread_turn && hist_ready && !allocation_failed) {
                for (hn = 1; hn <= cmd->sizeHistN; hn++) {
                    gd->histNN[hn] += hist.histNthread[hn];
                    gd->histNNSub[hn] += hist.histNNSubthread[hn];
                    gd->histNNSubXi2pcf[hn] +=
                        hist.histNNSubXi2pcfthread[hn];
#ifdef SMOOTHPIVOT
                    gd->histNNSubXi2pcftotal[hn] +=
                        hist.histNNSubXi2pcfthreadtotal[hn];
#endif
                    gd->histXi2pcf[hn] += hist.histXi2pcfthread[hn];
                }
#ifdef TPCF
                for (int m = 1; m <= cmd->mChebyshev + 1; m++) {
                    ADDM_ext(gd->histZetaMcos[m], gd->histZetaMcos[m],
                             hist.histZetaMthreadcos[m], cmd->sizeHistN);
                    ADDM_ext(gd->histZetaMsin[m], gd->histZetaMsin[m],
                             hist.histZetaMthreadsin[m], cmd->sizeHistN);
                    ADDM_ext(gd->histZetaMsincos[m], gd->histZetaMsincos[m],
                             hist.histZetaMthreadsincos[m], cmd->sizeHistN);
                    ADDM_ext(gd->histZetaMcossin[m], gd->histZetaMcossin[m],
                             hist.histZetaMthreadcossin[m], cmd->sizeHistN);
                }
#endif
                gd->nbbcalc += stats.nbbcalc;
                gd->nbccalc += stats.nbccalc;
#ifdef SMOOTHPIVOT
                ipfalse += stats.ipfalse;
                icountNbRmin += stats.count_rmin;
                icountNbRminOverlap += stats.count_overlap;
#endif
            }
        }
#pragma omp barrier
        if (hist_ready) search_free_sincos_omp(cmd, gd, &hist);
    }

#ifdef BALLTREEMPI
    if (distributed
        && fcfc_balltree_mpi_consensus(
               cmd, allocation_failed ? FAILURE : SUCCESS,
               "MPI ball-tree workers") == FAILURE)
        allocation_failed = TRUE;
#endif
    if (allocation_failed) {
        if (cmd->error_message[0] == '\0')
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     "searchcalc_balltree: worker histogram or MPI task failure");
        goto cleanup;
    }

#ifdef BALLTREEMPI
    if (distributed) {
        if (fcfc_balltree_mpi_scheduler_destroy(cmd, &scheduler) == FAILURE)
            goto cleanup;
        if (reduce_balltree_histograms(cmd, gd,
#ifdef SMOOTHPIVOT
                &ipfalse, &icountNbRmin, &icountNbRminOverlap
#else
                NULL, NULL, NULL
#endif
                ) == FAILURE)
            goto cleanup;
        if (!fcfc_balltree_mpi_is_root()) {
            status = SUCCESS;
            goto cleanup;
        }
    }
#endif

#ifdef SMOOTHPIVOT
    real xi, den, num;
    int mm;
        num = (real)nbody[cat1];
        den = (real)(nbody[cat1]-ipfalse);
#ifdef NOSTANDARNORMHIST
        xi = 1.0;
#else
        xi = num/den;
#endif // ! NONORMHIST
        verb_print(cmd->verbose,
                   "balltree: p falses found = %" INTEGER_FMT
                   " and %e %e %e\n",
                   ipfalse, num, den, xi);
#ifdef TPCF
            for (mm=1; mm<=cmd->mChebyshev+1; mm++) {
                MULMS_ext(gd->histZetaMcos[mm], gd->histZetaMcos[mm],
                          xi,cmd->sizeHistN);
                MULMS_ext(gd->histZetaMsin[mm], gd->histZetaMsin[mm],
                          xi,cmd->sizeHistN);
                MULMS_ext(gd->histZetaMsincos[mm], gd->histZetaMsincos[mm],
                          xi,cmd->sizeHistN);
                // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
                MULMS_ext(gd->histZetaMcossin[mm], gd->histZetaMcossin[mm],
                          xi,cmd->sizeHistN);
            }
#endif
#endif

    //B Normalization of histograms
        if (!cballs_opt_asymmetric(cmd)) {
            for (n = 1; n <= cmd->sizeHistN; n++) {
#ifdef SMOOTHPIVOT
                if (cmd->verbose>3)
                    printf("%d %e %e\n", n,
                       gd->histNNSubXi2pcf[n], gd->histNNSubXi2pcftotal[n]);
#else
                if (cmd->verbose>3)
                    printf("%d %e\n", n,
                       gd->histNNSubXi2pcf[n]);
#endif
                gd->histXi2pcf[n] /= 2.0;
                gd->histNNSubXi2pcf[n] /= 2.0;
#ifdef SMOOTHPIVOT
                gd->histNNSubXi2pcftotal[n] /= 2.0;
                    gd->histXi2pcf[n] /= MAX(gd->histNNSubXi2pcftotal[n],1.0);
#else
                    gd->histXi2pcf[n] /= MAX(gd->histNNSubXi2pcf[n],1.0);
#endif
    //E
            }
        } else {
            for (n = 1; n <= cmd->sizeHistN; n++) {
#ifdef SMOOTHPIVOT
                if (cmd->verbose>3)
                printf("%d %e %e\n", n,
                       gd->histNNSubXi2pcf[n], gd->histNNSubXi2pcftotal[n]);
#else
                if (cmd->verbose>3)
                printf("%d %e\n", n,
                       gd->histNNSubXi2pcf[n]);
#endif
#ifdef SMOOTHPIVOT
                    gd->histXi2pcf[n] /= MAX(gd->histNNSubXi2pcftotal[n],1.0);
#else
                    gd->histXi2pcf[n] /= MAX(gd->histNNSubXi2pcf[n],1.0);
#endif
            }
        }

    if (cballs_opt_compute_histn(cmd)) {
#ifdef SMOOTHPIVOT
            search_compute_HistN(cmd, gd, nbody[cat1]-ipfalse);
#else
            search_compute_HistN(cmd, gd, nbody[cat1]);
#endif
    }

#ifdef SMOOTHPIVOT
        verb_print(cmd->verbose,
                   "balltree: p falses found = %" INTEGER_FMT "\n", ipfalse);
        //B kappa Avg Rmin
        verb_print(cmd->verbose,
                   "balltree: count NbRmin found = %" INTEGER_FMT "\n",
                   icountNbRmin);
        verb_print(cmd->verbose,
                   "balltree: count overlap found = %" INTEGER_FMT "\n",
                   icountNbRminOverlap);
        
        bodyptr pp;
        INTEGER ifalsecount;
        ifalsecount = 0;
        INTEGER itruecount;
        itruecount = 0;
        for (pp = btab[cat1] + ipmin -1; pp < btab[cat1] + ipmax[cat1]; pp++) {
            if (Update(pp) == FALSE) {
                ifalsecount++;
            } else {
                itruecount++;
            }
        }
        verb_print(cmd->verbose, "balltree: p falses found = %" INTEGER_FMT
                   "\n", ifalsecount);
        verb_print(cmd->verbose, "balltree: p true found = %" INTEGER_FMT
                   "\n", itruecount);
        verb_print(cmd->verbose, "balltree: total = %" INTEGER_FMT "\n",
                   itruecount + ifalsecount);
        //E
#endif

    gd->cpusearch = CPUTIME - cpustart;
    verb_print(cmd->verbose,
               "%s: nbbcalc = %" INTEGER_FMT
               ", nbccalc = %" INTEGER_FMT "\n",
               cmd->searchMethod, gd->nbbcalc, gd->nbccalc);
    verb_print(cmd->verbose, "Going out: CPU time = %lf\n",CPUTIME-cpustart);

    status = SUCCESS;

cleanup:
#ifdef BALLTREEMPI
    if (scheduler.ready
        && fcfc_balltree_mpi_scheduler_destroy(cmd, &scheduler) == FAILURE)
        status = FAILURE;
    free(frontier);
    if (pivot_tree != NULL && pivot_tree != tree)
        fcfc_balltree_free(pivot_tree);
#endif
    fcfc_balltree_free(tree);
    return status;
}

global int searchcalc_balltree_omp(struct cmdline_data *cmd,
                                   struct global_data *gd,
                                   bodyptr *btab, INTEGER *nbody,
                                   INTEGER ipmin, INTEGER *ipmax,
                                   int cat1, int cat2)
{
    return searchcalc_balltree_driver(cmd, gd, btab, nbody, ipmin, ipmax,
                                      cat1, cat2, FALSE);
}

#ifdef BALLTREEMPI
global int searchcalc_balltree_mpi(struct cmdline_data *cmd,
                                   struct global_data *gd,
                                   bodyptr *btab, INTEGER *nbody,
                                   INTEGER ipmin, INTEGER *ipmax,
                                   int cat1, int cat2)
{
    return searchcalc_balltree_driver(cmd, gd, btab, nbody, ipmin, ipmax,
                                      cat1, cat2, TRUE);
}
#endif


local void walk_balltree_exact(struct cmdline_data *cmd,
                             struct global_data *gd, bodyptr p,
                             fcfc_balltreeptr kd, INTEGER *nbbcalcthread,
                             INTEGER *nbccalcthread,
                             gdhistptr_sincos_omp hist)
{
    fcfc_ballnode *ntab = kd->nodes;
    INTEGER stack[sizeof(INTEGER) * CHAR_BIT + 2];
    int top = 0;

    stack[top++] = 0;
    while (top != 0) {
        const INTEGER cp = stack[--top];
        const fcfc_ballnode *node = &ntab[cp];

        if (!balltree_intersects(cmd, gd, p, node)) continue;
        if (node->left >= 0) {
            stack[top++] = node->right;
            stack[top++] = node->left;
        } else {
            sumnode_sincos(cmd, gd, p, *node, kd,
                           nbbcalcthread, nbccalcthread, hist);
        }
    }
}


local void walk_balltree_one_ball(struct cmdline_data *cmd,
                                struct global_data *gd, bodyptr p,
                                fcfc_balltreeptr kd, INTEGER *nbbcalcthread,
                                INTEGER *nbccalcthread,
                                gdhistptr_sincos_omp hist)
{
    fcfc_ballnode *ntab = kd->nodes;
    bodyptr *bptr = kd->bptr;
    INTEGER stack[sizeof(INTEGER) * CHAR_BIT + 2];
    int top = 0;

    stack[top++] = 0;
    while (top != 0) {
        const INTEGER cp = stack[--top];
        const fcfc_ballnode *node = &ntab[cp];
        real dr1;
        real drpq2;
        compute_vector dr;

        if (!balltree_intersects(cmd, gd, p, node)) continue;
        if (node->left < 0) {
            sumnode_sincos(cmd, gd, p, *node, kd,
                           nbbcalcthread, nbccalcthread, hist);
            continue;
        }

        DOTPSUBV(drpq2, dr, Pos(p), node->cmpos);
        if (cmd->usePeriodic) {
            VWrapAll(dr);
            DOTVP(drpq2, dr, dr);
        }
        dr1 = rsqrt(drpq2);
        if (dr1 > 0.0 &&
            (Radius(p) + node->aggregate_radius) / dr1 < gd->deltaR) {
            sumnode_sincos_cell(cmd, gd, p, *node, bptr,
                                nbbcalcthread, nbccalcthread, hist);
        } else {
            stack[top++] = node->right;
            stack[top++] = node->left;
        }
    }
}


local void sumnode_sincos(struct  cmdline_data* cmd,
                          struct  global_data* gd, bodyptr p,
                          fcfc_ballnode ntab, fcfc_balltreeptr tree,
                          INTEGER *nbbcalcthread, INTEGER *nbccalcthread,
                          gdhistptr_sincos_omp hist)
{
#ifdef SINGLEP
    const fcfc_ballpoint *q;
#else
    bodyptr q;
#endif
    real dr1;
    compute_vector dr;
    int n;
    real xi;

    INTEGER pj;

    for (pj = ntab.first; pj <= ntab.last; ++pj) {
#ifdef SINGLEP
        q = &tree->packed_points[pj];
#else
        q = tree->bptr[pj];
#endif
        if (balltree_accept_body(cmd, gd, p, q, &dr1, dr)) {
            if (cmd->useLogHist) {
                if(dr1>cmd->rminHist) {
                    if (cmd->rminHist==0)
                        n = (int)(cmd->logHistBinsPD*(rlog10(dr1)
                                - rlog10(cmd->rangeN)) + cmd->sizeHistN) + 1;
                    else
                        n = (int)(rlog10(dr1/cmd->rminHist) * gd->i_deltaR) + 1;
                    if (n<=cmd->sizeHistN && n>=1) {
                        hist->histNthread[n] = hist->histNthread[n] + 1.;
                        hist->histNNSubXi2pcfthread[n] =
                        hist->histNNSubXi2pcfthread[n] + 1.;
#ifdef SMOOTHPIVOT
                        hist->histNNSubXi2pcfthreadp[n] =
                        hist->histNNSubXi2pcfthreadp[n] + 1.;
#endif
                        hist->histNNSubthread[n] = hist->histNNSubthread[n] + 1.;
#ifdef SINGLEP
                        xi = q->kappa;
#else
                        xi = Kappa(q);
#endif
#ifdef TPCF
                        REAL cosphi,sinphi;
#if NDIM == 3
                        REAL s, sy;
                        compute_vector pr0;
                        DOTVP(s, dr, hist->dr0);
                        cosphi = s/(dr1*hist->drpq);
                        CROSSVP(pr0,hist->dr0,Pos(p));
                        DOTVP(sy, dr, pr0);
#ifdef SINGLEP
                        if (rabs(cosphi)>1.0)
                            sinphi = 0.0;
                        else
                            sinphi = sqrt(1.0 - cosphi*cosphi);
#else
                        sinphi = rsqrt(1.0 - rsqr(cosphi));
#endif
                        if (sy < 0) sinphi *= -1.0;
#else // ! NDIM
                        cosphi = -dr[0]/dr1;
                        sinphi = -dr[1]/dr1;
#endif // ! NDIM
#ifdef SINGLEP
                        if (cosphi>1.0) cosphi = 1.0;
                        if (cosphi<-1.0) cosphi = -1.0;
#else
                        if (rabs(cosphi)>1.0)
                            verb_log_print(cmd->verbose, gd->outlog,
                        "sumenode: Warning!... cossphi must be in (-1,1): %g\n",
                                           cosphi);
#endif
                        CHEBYSHEVTUOMPSINCOS;
#endif
                        hist->histXi2pcfthreadsub[n] += xi;
                        *nbbcalcthread += 1;
                    }
                }
            } else { // ! useLogHist
                if(dr1>cmd->rminHist) {
                    n = (int) ( (dr1-cmd->rminHist) * gd->i_deltaR) + 1;
                    if (n<=cmd->sizeHistN && n>=1) {
                        hist->histNthread[n] = hist->histNthread[n] + 1.;
                        hist->histNNSubXi2pcfthread[n] =
                        hist->histNNSubXi2pcfthread[n] + 1.;
                        hist->histNNSubthread[n] = hist->histNNSubthread[n] + 1.;
#ifdef SINGLEP
                        xi = q->kappa;
#else
                        xi = Kappa(q);
#endif
#ifdef TPCF
                            real cosphi,sinphi;
#if NDIM == 3
                            real s, sy;
                            compute_vector pr0;
                            DOTVP(s, dr, hist->dr0);
                            cosphi = s/(dr1*hist->drpq);
                            CROSSVP(pr0,hist->dr0,Pos(p));
                            DOTVP(sy, dr, pr0);
                            sinphi = rsqrt(1.0 - rsqr(cosphi));;
                            if (sy < 0) sinphi *= -1.0;
#else // ! NDIM
                            cosphi = -dr[0]/dr1;
                            sinphi = -dr[1]/dr1;
#endif // ! NDIM
                            if (rabs(cosphi)>1.0)
                                verb_log_print(cmd->verbose, gd->outlog,
                        "sumenode: Warning!... cossphi must be in (-1,1): %g\n",
                                               cosphi);
                            CHEBYSHEVTUOMPSINCOS;
#endif
                        hist->histXi2pcfthreadsub[n] += xi;
                        *nbbcalcthread += 1;
                    }
                }
            } // ! useLogHist
        } // ! accept_body
    } // ! loop first to last
}


local void sumnode_sincos_cell(struct  cmdline_data* cmd,
                               struct  global_data* gd, bodyptr p,
                               fcfc_ballnode ntab, bodyptr *bptr,
                               INTEGER *nbbcalcthread, INTEGER *nbccalcthread,
                               gdhistptr_sincos_omp hist)
{
    real dr1;
    real drpq2;
    compute_vector dr;
    int n;
    real xi;

    int npoints;
    npoints = ntab.last - ntab.first + 1;
    DOTPSUBV(drpq2, dr, Pos(p), ntab.cmpos);
    if (cmd->usePeriodic) {
        VWrapAll(dr);
        DOTVP(drpq2, dr, dr);
    }
    dr1 = rsqrt(drpq2);
    if (dr1 < cmd->rangeN) {
        if (cmd->useLogHist) {
            if(dr1>cmd->rminHist) {
                if (cmd->rminHist==0)
                    n = (int)(cmd->logHistBinsPD*(rlog10(dr1)
                        - rlog10(cmd->rangeN)) + cmd->sizeHistN) + 1;
                else
                    n = (int)(rlog10(dr1/cmd->rminHist) * gd->i_deltaR) + 1;
                if (n<=cmd->sizeHistN && n>=1) {
                    hist->histNthread[n] = hist->histNthread[n] +  npoints;
                    hist->histNNSubXi2pcfthread[n] =
                    hist->histNNSubXi2pcfthread[n] + npoints;
                    hist->histNNSubthread[n] = hist->histNNSubthread[n] + npoints;
                        
                    xi = npoints*ntab.kappa;

#ifdef TPCF
                    REAL cosphi,sinphi;
#if NDIM == 3
                    REAL s, sy;
                    compute_vector pr0;
                    DOTVP(s, dr, hist->dr0);
                    cosphi = s/(dr1*hist->drpq);
                    CROSSVP(pr0,hist->dr0,Pos(p));
                    DOTVP(sy, dr, pr0);
#ifdef SINGLEP
                    if (rabs(cosphi)>1.0)
                        sinphi = 0.0;
                    else
                        sinphi = sqrt(1.0 - cosphi*cosphi);
#else
                    sinphi = rsqrt(1.0 - rsqr(cosphi));;
#endif
                    if (sy < 0) sinphi *= -1.0;
#else // ! NDIM
                    cosphi = -dr[0]/dr1;
                    sinphi = -dr[1]/dr1;
#endif // ! NDIM
#ifdef SINGLEP
                    if (cosphi>1.0) cosphi = 1.0;
                    if (cosphi<-1.0) cosphi = -1.0;
#else
                    if (rabs(cosphi)>1.0)
                        verb_log_print(cmd->verbose, gd->outlog,
                        "sumenode: Warning!... cossphi must be in (-1,1): %g\n",
                                        cosphi);
#endif
                    CHEBYSHEVTUOMPSINCOS;
#endif
                    hist->histXi2pcfthreadsub[n] += xi;
                    *nbccalcthread += 1;
                } // ! n in (1,sizeHistN)
            } // dr1 > rminHist
        } else {  // ! useLogHist
            if(dr1>cmd->rminHist) {
                n = (int) ( (dr1-cmd->rminHist) * gd->i_deltaR) + 1;
                if (n<=cmd->sizeHistN && n>=1) {
                    hist->histNthread[n] = hist->histNthread[n] +  npoints;
                    hist->histNNSubXi2pcfthread[n] =
                    hist->histNNSubXi2pcfthread[n] + npoints;
                    hist->histNNSubthread[n] = hist->histNNSubthread[n] + npoints;
                    xi = npoints*ntab.kappa;
#ifdef TPCF
                    real cosphi, sinphi;
#if NDIM == 3
                    real s, sy;
                    compute_vector pr0;
                    DOTVP(s, dr, hist->dr0);
                    cosphi = s/(dr1*hist->drpq);
                    CROSSVP(pr0,hist->dr0,Pos(p));
                    DOTVP(sy, dr, pr0);
                    sinphi = rsqrt(1.0 - rsqr(cosphi));;
                    if (sy < 0) sinphi *= -1.0;
#else // ! NDIM
                    cosphi = -dr[0]/dr1;
                    sinphi = -dr[1]/dr1;
#endif // ! NDIM
                    if (rabs(cosphi)>1.0)
                        verb_log_print(cmd->verbose, gd->outlog,
                        "sumenode: Warning!... cossphi must be in (-1,1): %g\n",
                                           cosphi);
                    CHEBYSHEVTUOMPSINCOS;
#endif
                    hist->histXi2pcfthreadsub[n] += xi;
                    *nbccalcthread += 1;
                } // ! n in (1, sizeHistN)
            } // ! dr1 > rminHist
        } // ! useLogHist
    } // ! dr1 < rangeN
}


local int print_info(struct cmdline_data* cmd,
                                  struct  global_data* gd)
{
    const bool behavior_ball = cballs_opt_behavior_ball(cmd);
    const bool no_one_ball = cballs_opt_no_one_ball(cmd);

    verb_print(cmd->verbose, "Search: Running ... (%s) \n",
               cmd->searchMethod);

    if (behavior_ball) {
        verb_print(cmd->verbose, "with option behavior-ball... \n");
        if (!no_one_ball && !cmd->useLogHist) {
//            error("behavior-ball and useLogHist=false are incompatible!");
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     "print_info: behavior-ball and useLogHist=false are incompatible");
            return FAILURE;
        }
    }
    if (no_one_ball) {
        verb_print(cmd->verbose, "with option no-one-ball... \n");
        if (behavior_ball)
            verb_print(cmd->verbose,
                       "no-one-ball disables behavior-ball cell aggregation... \n");
    }
#ifdef SMOOTHPIVOT
        verb_print(cmd->verbose,
                   "with option smooth-pivot... rsmooth=%g\n",gd->rsmooth[0]);
#endif
#ifndef TPCF
        verb_print(cmd->verbose, "computing only 2pcf... \n");
#endif
#ifdef NOSTANDARNORMHIST
    verb_print(cmd->verbose, "warning!! histograms will not be normalized... \n");
#endif

    return SUCCESS;
}
