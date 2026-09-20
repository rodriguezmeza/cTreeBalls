/* ==============================================================================
 MODULE: search_kdtree_omp.c        [cTreeBalls]
 Written by: M.A. Rodriguez-Meza.
 Starting date:    april 2023
 Purpose: 2/3-point correlation functions computation
 Language: C
 Use: kd = searchcalc_kdtree_omp(cmd, gd, btab, nbody,
                                    ipmin, ipmax, cat1, cat2);
 Major revisions:
 ==============================================================================*/
//        1          2          3          4        ^ 5          6          7

#include "globaldefs.h"

#include "kdtree.h"
#include "kdtree_parallel.h"
#include "kdtree_scan_frontier.h"

#ifndef KDTREE_OMP_PIVOT_BLOCK_SIZE
#define KDTREE_OMP_PIVOT_BLOCK_SIZE 64
#endif
#if KDTREE_OMP_PIVOT_BLOCK_SIZE < 1
#error "KDTREE_OMP_PIVOT_BLOCK_SIZE must be positive"
#endif

//B Some macros and definitions
//INTERSECT: macro to determine if node intersects search ball
#ifdef SINGLEP
#define KD_PRUNE_REAL cballs_storage_real
#define KD_PRUNE_LIMIT(r2ball)                                         \
    const cballs_storage_real _radius_limit =                          \
        cballs_store_search_bound(rsqrt((real)(r2ball)));              \
    const cballs_storage_real _r2limit =                               \
        _radius_limit * _radius_limit
#else
#define KD_PRUNE_REAL real
#define KD_PRUNE_LIMIT(r2ball)                                         \
    const real _r2limit = (real)(r2ball)
#endif

#if NDIM == 3
#define Intersect(node, r2ball, pos, done)                  \
{                                                           \
    KD_PRUNE_REAL _dxl, _dxr, _dyl, _dyr, _dzl, _dzr, _dr2; \
    KD_PRUNE_LIMIT(r2ball);                                 \
    _dxl = (KD_PRUNE_REAL)node.bnd.minb[0]                  \
         - (KD_PRUNE_REAL)pos[0];                           \
    _dxr = (KD_PRUNE_REAL)pos[0]                            \
         - (KD_PRUNE_REAL)node.bnd.maxb[0];                 \
    if (_dxl > 0.0) {                                       \
        _dr2 = _dxl*_dxl;                                   \
        if (_dr2 > _r2limit) goto done;                     \
    } else if (_dxr > 0.0) {                                \
        _dr2 = _dxr*_dxr;                                   \
        if (_dr2 > _r2limit) goto done;                     \
    } else                                                  \
        _dr2 = 0.0;                                         \
    _dyl = (KD_PRUNE_REAL)node.bnd.minb[1]                  \
         - (KD_PRUNE_REAL)pos[1];                           \
    _dyr = (KD_PRUNE_REAL)pos[1]                            \
         - (KD_PRUNE_REAL)node.bnd.maxb[1];                 \
    if (_dyl > 0.0) {                                       \
        _dr2 += _dyl*_dyl;                                  \
        if (_dr2 > _r2limit) goto done;                     \
    } else if (_dyr > 0.0) {                                \
        _dr2 += _dyr*_dyr;                                  \
        if (_dr2 > _r2limit) goto done;                     \
    }                                                       \
    _dzl = (KD_PRUNE_REAL)node.bnd.minb[2]                  \
         - (KD_PRUNE_REAL)pos[2];                           \
    _dzr = (KD_PRUNE_REAL)pos[2]                            \
         - (KD_PRUNE_REAL)node.bnd.maxb[2];                 \
    if (_dzl > 0.0) {                                       \
        _dr2 += _dzl*_dzl;                                  \
        if (_dr2 > _r2limit) goto done;                     \
    } else if (_dzr > 0.0) {                                \
        _dr2 += _dzr*_dzr;                                  \
        if (_dr2 > _r2limit) goto done;                     \
    }                                                       \
}

#else
#define Intersect(node, r2ball, pos, done)                  \
{                                                           \
    KD_PRUNE_REAL _dxl, _dxr, _dyl, _dyr, _dr2;             \
    KD_PRUNE_LIMIT(r2ball);                                 \
    _dxl = (KD_PRUNE_REAL)node.bnd.minb[0]                  \
         - (KD_PRUNE_REAL)pos[0];                           \
    _dxr = (KD_PRUNE_REAL)pos[0]                            \
         - (KD_PRUNE_REAL)node.bnd.maxb[0];                 \
    if (_dxl > 0.0) {                                       \
        _dr2 = _dxl*_dxl;                                   \
        if (_dr2 > _r2limit) goto done;                     \
    } else if (_dxr > 0.0) {                                \
        _dr2 = _dxr*_dxr;                                   \
        if (_dr2 > _r2limit) goto done;                     \
    } else                                                  \
        _dr2 = 0.0;                                         \
    _dyl = (KD_PRUNE_REAL)node.bnd.minb[1]                  \
         - (KD_PRUNE_REAL)pos[1];                           \
    _dyr = (KD_PRUNE_REAL)pos[1]                            \
         - (KD_PRUNE_REAL)node.bnd.maxb[1];                 \
    if (_dyl > 0.0) {                                       \
        _dr2 += _dyl*_dyl;                                  \
        if (_dr2 > _r2limit) goto done;                     \
    } else if (_dyr > 0.0) {                                \
        _dr2 += _dyr*_dyr;                                  \
        if (_dr2 > _r2limit) goto done;                     \
    }                                                       \
}

#endif
//E

local void sumnode_sincos(struct  cmdline_data*, struct  global_data*,
                          bodyptr, ballnode, ballxptr,
                          INTEGER *, INTEGER *,
                          gdhistptr_sincos_omp);
local void sumnode_sincos_cell(struct  cmdline_data*,
                               struct  global_data*, bodyptr,
                               ballnode, ballxptr,
                               INTEGER *, INTEGER *,
                               gdhistptr_sincos_omp);
local void walk_kdtree_exact(struct cmdline_data*, struct global_data*,
                            bodyptr, ballxptr, INTEGER *, INTEGER *,
                            gdhistptr_sincos_omp);
local void walk_kdtree_one_ball(struct cmdline_data*, struct global_data*,
                               bodyptr, ballxptr, INTEGER *, INTEGER *,
                               gdhistptr_sincos_omp);
local int print_info(struct cmdline_data* cmd,
                     struct  global_data* gd);

static int kdtree_reduce_results(struct cmdline_data *cmd,
                                 struct global_data *gd)
{
    if (!cballs_opt_only_3pcf(cmd)) {
        real *vectors[] = {gd->histNN, gd->histNNSubXi2pcf,
                           gd->histXi2pcf};
        for (size_t i = 0; i < sizeof(vectors)/sizeof(vectors[0]); i++)
            if (kdtree_reduce(cmd, vectors[i] + 1,
                              (size_t)cmd->sizeHistN) == FAILURE)
                return FAILURE;
#ifdef SMOOTHPIVOT
        if (kdtree_reduce(cmd, gd->histNNSubXi2pcftotal + 1,
                          (size_t)cmd->sizeHistN) == FAILURE)
            return FAILURE;
#endif
    }
#ifdef TPCF
    if (!cballs_opt_only_2pcf(cmd)) {
        if (kdtree_reduce(cmd, gd->histNNSub + 1,
                          (size_t)cmd->sizeHistN) == FAILURE)
            return FAILURE;
        real ***matrices[] = {gd->histZetaMcos, gd->histZetaMsin,
                              gd->histZetaMsincos, gd->histZetaMcossin};
        for (size_t c = 0; c < sizeof(matrices)/sizeof(matrices[0]); c++)
            for (int m = 1; m <= cmd->mChebyshev+1; m++)
                for (int n = 1; n <= cmd->sizeHistN; n++)
                    if (kdtree_reduce(cmd, matrices[c][m][n] + 1,
                                      (size_t)cmd->sizeHistN) == FAILURE)
                        return FAILURE;
    }
#endif
    INTEGER counts[3] = {gd->nbbcalc, gd->nbccalc, gd->ncccalc};
    if (kdtree_reduce_counts(cmd, counts, 3) == FAILURE) return FAILURE;
    if (kdtree_publish(cmd)) {
        gd->nbbcalc = counts[0];
        gd->nbccalc = counts[1];
        gd->ncccalc = counts[2];
    }
    return SUCCESS;
}

static void kdtree_finish_2pcf_pivot(struct cmdline_data *cmd, bodyptr p,
                                    gdhistptr_sincos_omp hist)
{
    real pivot_field = Weight(p)*Kappa(p);
#ifdef SMOOTHPIVOT
    if (cballs_opt_smooth_pivot(cmd)) pivot_field = KappaRmin(p);
#endif
    for (int n = 1; n <= cmd->sizeHistN; n++)
        hist->histXi2pcfthread[n] +=
            pivot_field*hist->histXi2pcfthreadsub[n];
}

static inline bool kdtree_accept_body(struct cmdline_data *cmd,
                                      struct global_data *gd,
                                      bodyptr p,
#ifdef SINGLEP
                                      const kd_leaf_point *q,
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

/*
 Search routine using kdtree method:

 To be called using: search=kdtree-omp

 Arguments:
    * `cmd`: Input: structure cmdline_data pointer
    * `gd`: Input: structure global_data pointer
    * `btable`: Input: point table array
    * `nbody`: Input: number of points in table array
    * `ipmin`: Input: minimum point in table array to analyse
    * `ipmax`: Input: maximum point in table array to analyse
    * `cat1`: Input: catalog tag to act as pivot catalog
    * `cat2`: Input: catalog tag to act as a scanning catalog
    * Global tructures used: gd, cmd
    * Histograms outputs (in global gd): histZetaMcos, histZetaMsin,
    *                                    histZetaMsincos, histN,
    *                                    histNNSubXi2pcf, histNNSubXi2pcftotal,
    *                                    histXi2pcf, histXi,
    * Counting encounters (in global gd): nbbcalc, nbccalc, ncccalc
 Return (the error status):
    int SUCCESS or FAILURE
 */
global int searchcalc_kdtree_omp(struct cmdline_data* cmd,
                                    struct  global_data* gd,
                                    bodyptr *btab, INTEGER *nbody,
                                    INTEGER ipmin, INTEGER *ipmax,
                                    int cat1, int cat2)
{
    bodyptr p;
    int n;
    double cpustart;
    ballxptr kd;
    int nbucket;
    real cpu_build_kdtree;
    const bool use_one_ball = !cballs_opt_no_one_ball(cmd);
    const bool only_2pcf = cballs_opt_only_2pcf(cmd);
    const bool only_3pcf = cballs_opt_only_3pcf(cmd);
#ifdef TWOPCF
    const bool run_2pcf = !only_3pcf;
#else
    const bool run_2pcf = FALSE;
#endif
#ifdef TPCF
    const bool run_3pcf = !only_2pcf;
#else
    const bool run_3pcf = FALSE;
#endif

    if (only_2pcf && only_3pcf)
        cBALLS_FAIL(cmd, "%s: only-2pcf and only-3pcf are mutually exclusive\n",
                    cmd->searchMethod);
    if (only_2pcf && !run_2pcf)
        cBALLS_FAIL(cmd, "%s: only-2pcf requires TWOPCFON=1\n",
                    cmd->searchMethod);
    if (only_3pcf && !run_3pcf)
        cBALLS_FAIL(cmd, "%s: only-3pcf requires TPCFON=1\n",
                    cmd->searchMethod);
    if (!run_2pcf && !run_3pcf)
        cBALLS_FAIL(cmd, "%s: build enables neither 2PCF nor 3PCF\n",
                    cmd->searchMethod);
    if (cballs_opt_edge_corrections(cmd) && only_2pcf)
        cBALLS_FAIL(cmd, "%s: edge-corrections require 3PCF; remove only-2pcf\n",
                    cmd->searchMethod);

    cpustart = CPUTIME;
    if (print_info(cmd, gd) == FAILURE)
        return FAILURE;

#ifdef OPENMPCODE
    ThreadCount(cmd, gd, nbody[cat1], cat1);
#endif

    search_init_gd_hist_sincos(cmd, gd);
    int allocation_failed = FALSE;
    INTEGER ipmask = 0;

#ifdef SMOOTHPIVOT
    INTEGER ipfalse;
    ipfalse=0;
    INTEGER icountNbRmin;
    icountNbRmin=0;
    INTEGER icountNbRminOverlap;
    icountNbRminOverlap=0;
#endif

    verb_print(cmd->verbose,
               "\nsearchcalc_balls: Total allocated %g MByte storage so far.\n",
               gd->bytes_tot/(1024.0*1024.0));


    DO_BODY(p,btab[cat1]+ipmin-1, btab[cat1]+ipmax[cat1])
        Update(p) = TRUE;
//B Building kd-tree
    cpu_build_kdtree = CPUTIME;
    //B version 1.0.1
    nbucket = cmd->nsmooth;
    //E
    verb_print(cmd->verbose, "\nkdtree build: nbucket = %d\n",nbucket);
    kd = init_kdtree(cmd, gd, btab[cat2], nbody[cat2]);
    if (kd == NULL) return FAILURE;
    if (build_kdtree(cmd, gd, kd, nbucket) == FAILURE) {
        finish_kdtree(kd);
        return FAILURE;
    }
    verb_print(cmd->verbose, "kdtree build: CPU time = %lf\n",
               CPUTIME-cpu_build_kdtree);
//E
    gd->ncellTable[cat2] = kd->nnode;               // Equivalent of octree cells

    INTEGER normalization_nbody = nbody[cat1];
    if (cballs_opt_read_mask(cmd)) {
        normalization_nbody = 0;
        for (bodyptr count_body = btab[cat1];
             count_body < btab[cat1] + nbody[cat1]; count_body++)
            if (Mask(count_body) != MASK_NODE_MASKED)
                normalization_nbody++;
    }
    if (normalization_nbody <= 0) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: mask selected no pivots", cmd->searchMethod);
        finish_kdtree(kd);
        return FAILURE;
    }

#ifdef SMOOTHPIVOT
    if (prepare_smooth_pivots(cmd, gd, btab, nbody,
                              ipmin, ipmax, cat1, cat2) == FAILURE) {
        finish_kdtree(kd);
        return FAILURE;
    }
#endif

    if (cballs_opt_edge_corrections(cmd)) {
        int edge_status = kdtree_edge_search(cmd, gd, btab, nbody, ipmin,
                                              ipmax, cat1, cat2, kd);
        finish_kdtree(kd);
        gd->cpusearch = CPUTIME - cpustart;
        return edge_status;
    }

    const INTEGER pivot_count = ipmax[cat1] - ipmin + 1;
#ifdef BALLS4SCANLEV
    const bool complete_auto_catalog = cat1 == cat2 && ipmin == 1
        && ipmax[cat1] == nbody[cat1] && btab[cat1] == kd->body_base;
    const kdtree_scan_frontier pivot_frontier =
        kdtree_scan_frontier_make(kd, pivot_count, complete_auto_catalog);
    const INTEGER block_count = pivot_frontier.count;
    verb_print(cmd->verbose,
               "%s: KD scan-level frontier has %" INTEGER_FMT
               " tasks in %s order\n",
               cmd->searchMethod, block_count,
               pivot_frontier.tree_order ? "tree" : "catalog");
#else
    const INTEGER block_count =
        1 + (pivot_count - 1)/(INTEGER)KDTREE_OMP_PIVOT_BLOCK_SIZE;
    const kdtree_scan_frontier pivot_frontier = {block_count, FALSE};
#endif

#ifdef SMOOTHPIVOT
#pragma omp parallel default(none)   \
    shared(cmd,gd,btab,nbody,roottable,ipmin,ipmax, \
    rootnode, cat1, cat2, kd, ipfalse, icountNbRmin, icountNbRminOverlap, \
    allocation_failed, use_one_ball, run_2pcf, run_3pcf, ipmask, \
    normalization_nbody, block_count, pivot_count, pivot_frontier)
#else
#pragma omp parallel default(none)   \
    shared(cmd,gd,btab,nbody,roottable,ipmin,ipmax, \
    rootnode, cat1, cat2, kd, allocation_failed, use_one_ball, \
    run_2pcf, run_3pcf, ipmask, normalization_nbody, block_count, pivot_count, \
    pivot_frontier)
#endif
    {
        gdhist_sincos_omp hist;
        const int hist_allocated =
            search_init_sincos_omp(cmd, gd, &hist) == SUCCESS;
        int hist_ready = hist_allocated;
        real *histNNSubBlock = NULL;
        if (hist_ready && run_3pcf) {
            histNNSubBlock = calloc((size_t)cmd->sizeHistN + 1,
                                    sizeof(*histNNSubBlock));
            if (histNNSubBlock == NULL)
                hist_ready = FALSE;
        }
        if (!hist_ready) {
#pragma omp atomic write
            allocation_failed = TRUE;
        }

#pragma omp barrier

#ifdef BALLS4SCANLEV
#pragma omp for schedule(dynamic,1) ordered
#else
#pragma omp for schedule(static,1) ordered
#endif
        for (INTEGER block = 0; block < block_count; block++) {
            const bool owned = kdtree_task_owned(cmd, block);
            INTEGER nbbcalcblock = 0;
            INTEGER nbccalcblock = 0;
            INTEGER ipfalseblock = 0;
            INTEGER ipmaskblock = 0;
#ifdef SMOOTHPIVOT
            INTEGER icountNbRminblock = 0;
            INTEGER icountNbRminOverlapblock = 0;
#endif

            if (hist_ready && !allocation_failed && owned) {
                for (int n = 1; n <= cmd->sizeHistN; n++) {
                    hist.histNthread[n] = 0.0;
                    hist.histNNSubthread[n] = 0.0;
                    hist.histNNSubXi2pcfthread[n] = 0.0;
#ifdef SMOOTHPIVOT
                    hist.histNNSubXi2pcfthreadp[n] = 0.0;
                    hist.histNNSubXi2pcfthreadtotal[n] = 0.0;
#endif
                    hist.histXi2pcfthread[n] = 0.0;
                    hist.histXi2pcfthreadsub[n] = 0.0;
                    if (histNNSubBlock != NULL)
                        histNNSubBlock[n] = 0.0;
                }
#ifdef TPCF
                if (run_3pcf) {
                    for (int m = 1; m <= cmd->mChebyshev+1; m++) {
                        CLRM_ext(hist.histZetaMthreadcos[m],
                                 cmd->sizeHistN);
                        CLRM_ext(hist.histZetaMthreadsin[m],
                                 cmd->sizeHistN);
                        CLRM_ext(hist.histZetaMthreadsincos[m],
                                 cmd->sizeHistN);
                        CLRM_ext(hist.histZetaMthreadcossin[m],
                                 cmd->sizeHistN);
                    }
                }
#endif

                INTEGER first;
                INTEGER end;
#ifdef BALLS4SCANLEV
                kdtree_scan_frontier_range(
                    &pivot_frontier, kd, pivot_count, ipmin - 1,
                    block, &first, &end);
#else
                first = ipmin - 1
                      + block*(INTEGER)KDTREE_OMP_PIVOT_BLOCK_SIZE;
                end = MIN(first + (INTEGER)KDTREE_OMP_PIVOT_BLOCK_SIZE,
                          ipmax[cat1]);
#endif
                for (INTEGER pivot_index = first;
                     pivot_index < end; pivot_index++) {
                    bodyptr p = kdtree_scan_frontier_body(
                        &pivot_frontier,
                        kd, btab[cat1], pivot_index);

                    if (cballs_opt_read_mask(cmd)
                        && Mask(p) == MASK_NODE_MASKED) {
                        ipmaskblock++;
                        continue;
                    }
#ifdef SMOOTHPIVOT
                    if (cballs_opt_smooth_pivot(cmd) && Update(p) == FALSE) {
                        ipfalseblock++;
                        continue;
                    }
#endif
                    for (int n = 1; n <= cmd->sizeHistN; n++) {
                        hist.histNNSubthread[n] = 0.0;
                        hist.histXi2pcfthreadsub[n] = 0.0;
#ifdef SMOOTHPIVOT
                        hist.histNNSubXi2pcfthreadp[n] = 0.0;
#endif
                    }
#ifdef TPCF
                    if (run_3pcf) {
                        CLRM_ext_ext(hist.histXithreadcos,
                                     cmd->mChebyshev+1, cmd->sizeHistN);
                        CLRM_ext_ext(hist.histXithreadsin,
                                     cmd->mChebyshev+1, cmd->sizeHistN);
#if NDIM == 3
                        dRotation3D(Pos(p), ROTANGLE, ROTANGLE, ROTANGLE,
                                    hist.q0);
                        DOTPSUBV(hist.drpq2, hist.dr0, Pos(p), hist.q0);
                        hist.drpq = rsqrt(hist.drpq2);
#ifdef PTOPIVOTROTATION
                        real rtheta = xrandom(0.0, TWOPI);
                        compute_vector dr0rot;
                        RotationVecAWRtoVecB(dr0rot, hist.dr0, Pos(p), rtheta);
                        SETV(hist.dr0, dr0rot);
#endif
#endif
                    }
#endif

                    if (use_one_ball)
                        walk_kdtree_one_ball(cmd, gd, p, kd,
                                             &nbbcalcblock, &nbccalcblock,
                                             &hist);
                    else
                        walk_kdtree_exact(cmd, gd, p, kd,
                                          &nbbcalcblock, &nbccalcblock,
                                          &hist);

#ifdef SMOOTHPIVOT
                    for (int n = 1; n <= cmd->sizeHistN; n++) {
                        if (run_2pcf) {
                            hist.histNNSubXi2pcfthreadp[n] =
                                ((real)NbRmin(p))
                                * hist.histNNSubXi2pcfthreadp[n];
                            hist.histNNSubXi2pcfthreadtotal[n] +=
                                hist.histNNSubXi2pcfthreadp[n];
                        }
                        if (run_3pcf)
                            hist.histNNSubthread[n] =
                                ((real)NbRmin(p))*hist.histNNSubthread[n];
                    }
#endif

                    if (run_3pcf) {
                        computeBodyProperties_sincos(
                            cmd, gd, p, normalization_nbody, &hist);
                        for (int n = 1; n <= cmd->sizeHistN; n++)
                            histNNSubBlock[n] += hist.histNNSubthread[n];
                    } else if (run_2pcf) {
                        kdtree_finish_2pcf_pivot(cmd, p, &hist);
                    }

#ifdef SMOOTHPIVOT
                    icountNbRminblock += NbRmin(p);
                    icountNbRminOverlapblock += NbRminOverlap(p);
#endif
                    const INTEGER ip = p - btab[cat1] + 1;
                    if (ip%cmd->stepState == 0)
                        verb_log_print(cmd->verbose_log, gd->outlog,
                                       " - Completed pivot: %ld\n", ip);
                }
            }

#pragma omp ordered
            {
            if (hist_ready && !allocation_failed && owned) {
                for (int n = 1; n <= cmd->sizeHistN; n++) {
                    if (run_2pcf) {
                        gd->histNN[n] += hist.histNthread[n];
                        gd->histNNSubXi2pcf[n] +=
                            hist.histNNSubXi2pcfthread[n];
#ifdef SMOOTHPIVOT
                        gd->histNNSubXi2pcftotal[n] +=
                            hist.histNNSubXi2pcfthreadtotal[n];
#endif
                        gd->histXi2pcf[n] += hist.histXi2pcfthread[n];
                    }
                    if (run_3pcf)
                        gd->histNNSub[n] += histNNSubBlock[n];
                }
#ifdef TPCF
                if (run_3pcf) {
                    for (int m=1; m<=cmd->mChebyshev+1; m++) {
                        ADDM_ext(gd->histZetaMcos[m],gd->histZetaMcos[m],
                                 hist.histZetaMthreadcos[m],cmd->sizeHistN);
                        ADDM_ext(gd->histZetaMsin[m],gd->histZetaMsin[m],
                                 hist.histZetaMthreadsin[m],cmd->sizeHistN);
                        ADDM_ext(gd->histZetaMsincos[m],
                                 gd->histZetaMsincos[m],
                                 hist.histZetaMthreadsincos[m],
                                 cmd->sizeHistN);
                        ADDM_ext(gd->histZetaMcossin[m],
                                 gd->histZetaMcossin[m],
                                 hist.histZetaMthreadcossin[m],
                                 cmd->sizeHistN);
                    }
                }
#endif
                gd->nbbcalc += nbbcalcblock;
                gd->nbccalc += nbccalcblock;
                ipmask += ipmaskblock;
#ifdef SMOOTHPIVOT
                ipfalse += ipfalseblock;
                icountNbRmin += icountNbRminblock;
                icountNbRminOverlap += icountNbRminOverlapblock;
#endif
            }
            }
        }

        free(histNNSubBlock);
        if (hist_allocated)
            search_free_sincos_omp(cmd, gd, &hist);
    } // end pragma omp parallel

    if (allocation_failed) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "searchcalc_kdtree_omp: OpenMP histogram allocation failed");
        finish_kdtree(kd);
        return FAILURE;
    }

    INTEGER selection_counts[4] = {ipmask, 0, 0, 0};
#ifdef SMOOTHPIVOT
    selection_counts[1] = ipfalse;
    selection_counts[2] = icountNbRmin;
    selection_counts[3] = icountNbRminOverlap;
#endif
    if (kdtree_reduce_results(cmd, gd) == FAILURE
        || kdtree_reduce_counts(cmd, selection_counts, 4) == FAILURE) {
        finish_kdtree(kd);
        return FAILURE;
    }
    if (kdtree_publish(cmd)) {
        ipmask = selection_counts[0];
#ifdef SMOOTHPIVOT
        ipfalse = selection_counts[1];
        icountNbRmin = selection_counts[2];
        icountNbRminOverlap = selection_counts[3];
#endif
    }

    int selection_status = SUCCESS;
    if (kdtree_publish(cmd)
        && nbody[cat1] - ipmask
#ifdef SMOOTHPIVOT
           - (cballs_opt_smooth_pivot(cmd) ? ipfalse : 0)
#endif
           <= 0) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: mask/smoothing selected no pivots", cmd->searchMethod);
        selection_status = FAILURE;
    }
    if (kdtree_consensus(cmd, selection_status,
                         "KDTREE pivot selection") == FAILURE) {
        finish_kdtree(kd);
        return FAILURE;
    }

    if (!kdtree_publish(cmd)) {
        finish_kdtree(kd);
        gd->cpusearch = CPUTIME - cpustart;
        return SUCCESS;
    }

#ifdef SMOOTHPIVOT
    real xi, den, num;
    int mm;
        num = (real)normalization_nbody;
        den = (real)(normalization_nbody-ipfalse);
#ifdef NOSTANDARNORMHIST
        xi = 1.0;
#else
        xi = cballs_raw_legacy_multipoles(cmd)
            ? 1.0 : cballs_normalize_or_zero(num, den);
#endif // ! NONORMHIST
        verb_print(cmd->verbose,
                   "%s: p falses found = %" INTEGER_FMT " and %e %e %e\n",
                   cmd->searchMethod,
                   ipfalse, num, den, xi);
#ifdef TPCF
        if (run_3pcf) {
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
        }
#endif
#endif

    //B Normalization of histograms
        if (run_2pcf && !cballs_opt_asymmetric(cmd)) {
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
        } else if (run_2pcf) {
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

    if (run_2pcf && cballs_opt_compute_histn(cmd)) {
#ifdef SMOOTHPIVOT
            if (search_compute_HistN(cmd, gd, nbody[cat1]-ipfalse-ipmask) == FAILURE) {
                finish_kdtree(kd);
                return FAILURE;
            }
#else
            if (search_compute_HistN(cmd, gd, nbody[cat1]-ipmask) == FAILURE) {
                finish_kdtree(kd);
                return FAILURE;
            }
#endif
    }

#ifdef SMOOTHPIVOT
        verb_print(cmd->verbose, "%s: p falses found = %" INTEGER_FMT "\n",
                   cmd->searchMethod, ipfalse);
        //B kappa Avg Rmin
        verb_print(cmd->verbose,
                   "%s: count NbRmin found = %" INTEGER_FMT "\n",
                   cmd->searchMethod, icountNbRmin);
        verb_print(cmd->verbose,
                   "%s: count overlap found = %" INTEGER_FMT "\n",
                   cmd->searchMethod, icountNbRminOverlap);

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
        verb_print(cmd->verbose, "%s: p falses found = %" INTEGER_FMT "\n",
                   cmd->searchMethod, ifalsecount);
        verb_print(cmd->verbose, "%s: p true found = %" INTEGER_FMT "\n",
                   cmd->searchMethod, itruecount);
        verb_print(cmd->verbose, "%s: total = %" INTEGER_FMT "\n",
                   cmd->searchMethod, itruecount+ifalsecount);
        //E
#endif

    gd->cpusearch = CPUTIME - cpustart;
    verb_print(cmd->verbose, "Going out: CPU time = %lf\n",CPUTIME-cpustart);

    finish_kdtree(kd);
    return SUCCESS;
}


local void walk_kdtree_exact(struct cmdline_data *cmd,
                             struct global_data *gd, bodyptr p,
                             ballxptr kd, INTEGER *nbbcalcthread,
                             INTEGER *nbccalcthread,
                             gdhistptr_sincos_omp hist)
{
    ballnode *ntab = kd->ntab;
    INTEGER cp = KDROOT;

    do {
        Intersect(ntab[cp], gd->RcutSq, Pos(p), exact_next_cell);
        if (ntab[cp].valid_count == 0) goto exact_next_cell;
        if (cp < kd->nsplit) {
            cp = Lower(cp);
            continue;
        }

        sumnode_sincos(cmd, gd, p, ntab[cp], kd,
                       nbbcalcthread, nbccalcthread, hist);

exact_next_cell:
        SetNext(cp);
    } while (cp != KDROOT);
}


local void walk_kdtree_one_ball(struct cmdline_data *cmd,
                                struct global_data *gd, bodyptr p,
                                ballxptr kd, INTEGER *nbbcalcthread,
                                INTEGER *nbccalcthread,
                                gdhistptr_sincos_omp hist)
{
    ballnode *ntab = kd->ntab;
    INTEGER cp = KDROOT;

    do {
        real dr1;
        real drpq2;
        compute_vector dr;

        Intersect(ntab[cp], gd->RcutSq, Pos(p), one_ball_next_cell);
        if (ntab[cp].valid_count == 0) goto one_ball_next_cell;
        if (cp < kd->nsplit) {
            DOTPSUBV(drpq2, dr, Pos(p), ntab[cp].cmpos);
            if (cmd->usePeriodic) {
                VWrapAll(dr);
                DOTVP(drpq2, dr, dr);
            }
            dr1 = rsqrt(drpq2);
            if (!kdtree_node_contains_body(kd, &ntab[cp], p)
                && dr1 > 0.0
                && (Radius(p) + ntab[cp].bnd.radius)/dr1 < gd->deltaR
                && cballs_angular_cell_ok(cmd, Pos(p), dr, 0.0,
                                          ntab[cp].bnd.geometric_radius)) {
                sumnode_sincos_cell(cmd, gd, p, ntab[cp], kd,
                                    nbbcalcthread, nbccalcthread, hist);
                SetNext(cp);
                continue;
            }

            cp = Lower(cp);
            continue;
        }

        sumnode_sincos(cmd, gd, p, ntab[cp], kd,
                       nbbcalcthread, nbccalcthread, hist);

one_ball_next_cell:
        SetNext(cp);
    } while (cp != KDROOT);
}


local void sumnode_sincos(struct  cmdline_data* cmd,
                          struct  global_data* gd, bodyptr p,
                          ballnode ntab, ballxptr kd,
                          INTEGER *nbbcalcthread, INTEGER *nbccalcthread,
                          gdhistptr_sincos_omp hist)
{
#ifdef SINGLEP
    const kd_leaf_point *q;
#else
    bodyptr q;
#endif
    real dr1;
    compute_vector dr;
    int n;
    real xi;
#ifdef TWOPCF
    const bool run_2pcf = !cballs_opt_only_3pcf(cmd);
#else
    const bool run_2pcf = FALSE;
#endif
#ifdef TPCF
    const bool run_3pcf = !cballs_opt_only_2pcf(cmd);
#else
    const bool run_3pcf = FALSE;
#endif

    INTEGER pj;

    for (pj = ntab.first; pj <= ntab.last; ++pj) {
#ifdef SINGLEP
        q = &kd->packed_points[pj];
#else
        q = kd->bptr[pj];
#endif
#ifdef SINGLEP
        bodyptr source_q = kd->bptr[pj];
#else
        bodyptr source_q = q;
#endif
        if (source_q == p || (cballs_opt_read_mask(cmd)
            && Mask(source_q) == MASK_NODE_MASKED)) continue;
        if (kdtree_accept_body(cmd, gd, p, q, &dr1, dr)) {
            if (cmd->useLogHist) {
                if(dr1>cmd->rminHist) {
                    if (cmd->rminHist==0)
                        n = (int)(cmd->logHistBinsPD*(rlog10(dr1)
                                - rlog10(cmd->rangeN)) + cmd->sizeHistN) + 1;
                    else
                        n = (int)(rlog10(dr1/cmd->rminHist) * gd->i_deltaR) + 1;
                    if (n<=cmd->sizeHistN && n>=1) {
                        if (run_2pcf) {
                            hist->histNthread[n] += 1.0;
                            hist->histNNSubXi2pcfthread[n] += 1.0;
#ifdef SMOOTHPIVOT
                            hist->histNNSubXi2pcfthreadp[n] += 1.0;
#endif
                        }
                        if (run_3pcf) hist->histNNSubthread[n] += 1.0;
#ifdef SINGLEP
                        xi = cballs_raw_legacy_multipoles(cmd) ? q->weighted_kappa : q->kappa;
#else
                        xi = cballs_raw_legacy_multipoles(cmd) ? Weight(q)*Kappa(q) : Kappa(q);
#endif
#ifdef TPCF
                        if (run_3pcf) {
                        real cosphi, sinphi;
                        if (cballs_angular_phase(Pos(p), dr, &cosphi, &sinphi)) {
                            if (cballs_raw_legacy_multipoles(cmd))
                                cballs_accumulate_raw_moments(cmd, p, hist, n,
                                    xi, xi*xi, cosphi, sinphi);
                            else
                                CHEBYSHEVTUOMPSINCOS;
                        }
                        }
#endif
                        if (run_2pcf)
                            hist->histXi2pcfthreadsub[n] += xi;
                        *nbbcalcthread += 1;
                    }
                }
            } else { // ! useLogHist
                if(dr1>cmd->rminHist) {
                    n = (int) ( (dr1-cmd->rminHist) * gd->i_deltaR) + 1;
                    if (n<=cmd->sizeHistN && n>=1) {
                        if (run_2pcf) {
                            hist->histNthread[n] += 1.0;
                            hist->histNNSubXi2pcfthread[n] += 1.0;
#ifdef SMOOTHPIVOT
                            hist->histNNSubXi2pcfthreadp[n] += 1.0;
#endif
                        }
                        if (run_3pcf) hist->histNNSubthread[n] += 1.0;
#ifdef SINGLEP
                        xi = cballs_raw_legacy_multipoles(cmd) ? q->weighted_kappa : q->kappa;
#else
                        xi = cballs_raw_legacy_multipoles(cmd) ? Weight(q)*Kappa(q) : Kappa(q);
#endif
#ifdef TPCF
                        if (run_3pcf) {
                            real cosphi, sinphi;
                            if (cballs_angular_phase(Pos(p), dr, &cosphi, &sinphi)) {
                                if (cballs_raw_legacy_multipoles(cmd))
                                    cballs_accumulate_raw_moments(cmd, p, hist, n,
                                        xi, xi*xi, cosphi, sinphi);
                                else
                                    CHEBYSHEVTUOMPSINCOS;
                            }
                        }
#endif
                        if (run_2pcf)
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
                               ballnode ntab, ballxptr kd,
                               INTEGER *nbbcalcthread, INTEGER *nbccalcthread,
                               gdhistptr_sincos_omp hist)
{
    real dr1;
    real drpq2;
    compute_vector dr;
    int n;
    real xi;
#ifdef TWOPCF
    const bool run_2pcf = !cballs_opt_only_3pcf(cmd);
#else
    const bool run_2pcf = FALSE;
#endif
#ifdef TPCF
    const bool run_3pcf = !cballs_opt_only_2pcf(cmd);
#else
    const bool run_3pcf = FALSE;
#endif

    const real npoints = (real)ntab.valid_count;
    if (ntab.valid_count == 0 || kdtree_node_contains_body(kd, &ntab, p))
        return;
    DOTPSUBV(drpq2, dr, Pos(p), ntab.cmpos);
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
                    if (run_2pcf) {
                        hist->histNthread[n] += npoints;
                        hist->histNNSubXi2pcfthread[n] += npoints;
#ifdef SMOOTHPIVOT
                        hist->histNNSubXi2pcfthreadp[n] += npoints;
#endif
                    }
                    if (run_3pcf) hist->histNNSubthread[n] += npoints;

                    xi = cballs_raw_legacy_multipoles(cmd) ? ntab.weighted_kappa_sum : npoints*ntab.kappa;

#ifdef TPCF
                    if (run_3pcf) {
                    real cosphi, sinphi;
                    if (cballs_angular_phase(Pos(p), dr, &cosphi, &sinphi)) {
                        if (cballs_raw_legacy_multipoles(cmd))
                            cballs_accumulate_raw_moments(cmd, p, hist, n,
                                xi, ntab.weighted_kappa_sq_sum, cosphi, sinphi);
                        else
                            CHEBYSHEVTUOMPSINCOS;
                    }
                    }
#endif
                    if (run_2pcf)
                        hist->histXi2pcfthreadsub[n] += xi;
                    *nbccalcthread += 1;
                } // ! n in (1,sizeHistN)
            } // dr1 > rminHist
        } else {  // ! useLogHist
            if(dr1>cmd->rminHist) {
                n = (int) ( (dr1-cmd->rminHist) * gd->i_deltaR) + 1;
                if (n<=cmd->sizeHistN && n>=1) {
                    if (run_2pcf) {
                        hist->histNthread[n] += npoints;
                        hist->histNNSubXi2pcfthread[n] += npoints;
#ifdef SMOOTHPIVOT
                        hist->histNNSubXi2pcfthreadp[n] += npoints;
#endif
                    }
                    if (run_3pcf) hist->histNNSubthread[n] += npoints;
                    xi = cballs_raw_legacy_multipoles(cmd) ? ntab.weighted_kappa_sum : npoints*ntab.kappa;
#ifdef TPCF
                    if (run_3pcf) {
                    real cosphi, sinphi;
                    if (cballs_angular_phase(Pos(p), dr, &cosphi, &sinphi)) {
                        if (cballs_raw_legacy_multipoles(cmd))
                            cballs_accumulate_raw_moments(cmd, p, hist, n,
                                xi, ntab.weighted_kappa_sq_sum, cosphi, sinphi);
                        else
                            CHEBYSHEVTUOMPSINCOS;
                    }
                    }
#endif
                    if (run_2pcf)
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

    verb_print(cmd->verbose, "Search: Running ... (%s) \n", cmd->searchMethod);

    if (behavior_ball) {
        verb_print(cmd->verbose,
                   "behavior-ball is explicit; cell aggregation is already the default... \n");
    } else if (!no_one_ball) {
        verb_print(cmd->verbose, "with default one-ball cell aggregation... \n");
    }
    if (no_one_ball) {
        verb_print(cmd->verbose, "with option no-one-ball... \n");
        verb_print(cmd->verbose,
                   "no-one-ball disables default cell aggregation... \n");
    }
    if (cballs_opt_smooth_pivot(cmd))
        verb_print(cmd->verbose,
                   "with compiled-default smooth-pivot... rsmooth=%g\n",gd->rsmooth[0]);
    else
        verb_print(cmd->verbose, "with no-smooth-pivot... \n");
    if (cballs_opt_only_2pcf(cmd))
        verb_print(cmd->verbose, "computing only 2pcf... \n");
    if (cballs_opt_only_3pcf(cmd))
        verb_print(cmd->verbose, "computing only 3pcf... \n");
#ifndef TPCF
        verb_print(cmd->verbose, "computing only 2pcf... \n");
#endif
#ifdef NOSTANDARNORMHIST
    verb_print(cmd->verbose, "warning!! histograms will not be normalized... \n");
#endif

    return SUCCESS;
}
