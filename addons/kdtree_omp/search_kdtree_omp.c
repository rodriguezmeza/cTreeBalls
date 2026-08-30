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
                               ballnode, bodyptr *,
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
    const bool use_one_ball =
        cballs_opt_behavior_ball(cmd)
        && !cballs_opt_no_one_ball(cmd);

    cpustart = CPUTIME;
    if (print_info(cmd, gd) == FAILURE)
        return FAILURE;

#ifdef OPENMPCODE
    ThreadCount(cmd, gd, nbody[cat1], cat1);
#endif

    search_init_gd_hist_sincos(cmd, gd);
    int allocation_failed = FALSE;

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
    if (build_kdtree(cmd, gd, kd, nbucket) == FAILURE) {
        finish_kdtree(kd);
        return FAILURE;
    }
    verb_print(cmd->verbose, "kdtree build: CPU time = %lf\n",
               CPUTIME-cpu_build_kdtree);
//E
    gd->ncellTable[cat1] = kd->nnode;               // Equivalent of octree cells

#ifdef SMOOTHPIVOT
    if (prepare_smooth_pivots(cmd, gd, btab, nbody,
                              ipmin, ipmax, cat1, cat2) == FAILURE) {
        finish_kdtree(kd);
        return FAILURE;
    }
#endif

#ifdef SMOOTHPIVOT
#pragma omp parallel default(none)   \
    shared(cmd,gd,btab,nbody,roottable,ipmin,ipmax, \
    rootnode, cat1, cat2, kd, ipfalse, icountNbRmin, icountNbRminOverlap, \
    allocation_failed, use_one_ball)
#else
#pragma omp parallel default(none)   \
    shared(cmd,gd,btab,nbody,roottable,ipmin,ipmax, \
    rootnode, cat1, cat2, kd, allocation_failed, use_one_ball)
#endif
    {
        bodyptr p;
        int n;
        INTEGER nbbcalcthread = 0;
        INTEGER nbccalcthread = 0;
        
        gdhist_sincos_omp hist;
        int hist_ready =
            search_init_sincos_omp(cmd, gd, &hist) == SUCCESS;
        if (!hist_ready) {
#pragma omp atomic write
            allocation_failed = TRUE;
        }

#pragma omp barrier

        INTEGER ipfalsethreads;
        ipfalsethreads = 0;

#ifdef SMOOTHPIVOT
        INTEGER icountNbRminthread;
        icountNbRminthread=0;
        INTEGER icountNbRminOverlapthread;
        icountNbRminOverlapthread=0;
#endif
#pragma omp for nowait schedule(static,1)
    DO_BODY(p, btab[cat1]+ipmin-1, btab[cat1]+ipmax[cat1]) {
        if (allocation_failed) continue;
// p and q are in differents node structures... cat1!=cat2...
#ifdef SMOOTHPIVOT
        if (Update(p) == FALSE) {
            ipfalsethreads++;
            continue;
        }
#endif
        for (n = 1; n <= cmd->sizeHistN; n++) {
            hist.histNNSubthread[n] = 0.0;          // Affects only 3pcf
            hist.histXi2pcfthreadsub[n] = 0.0;      // Affects only 2pcf
#ifdef SMOOTHPIVOT
            hist.histNNSubXi2pcfthreadp[n] = 0.;    // Affects only 2pcf
#endif
        }
#ifdef TPCF
            CLRM_ext_ext(hist.histXithreadcos, cmd->mChebyshev+1,
                         cmd->sizeHistN);
            CLRM_ext_ext(hist.histXithreadsin, cmd->mChebyshev+1, 
                         cmd->sizeHistN);
#if NDIM == 3
            dRotation3D(Pos(p), ROTANGLE, ROTANGLE, ROTANGLE, hist.q0);
            DOTPSUBV(hist.drpq2, hist.dr0, Pos(p), hist.q0);
            hist.drpq = rsqrt(hist.drpq2);
            //B Random rotation of dr0:
#ifdef PTOPIVOTROTATION
            real rtheta;
            compute_vector dr0rot;
            rtheta = xrandom(0.0, TWOPI);
            RotationVecAWRtoVecB(dr0rot, hist.dr0, Pos(p), rtheta);
            SETV(hist.dr0, dr0rot);
#endif
            //E
#endif // ! NDIM
#endif

        if (use_one_ball)
            walk_kdtree_one_ball(cmd, gd, p, kd,
                                 &nbbcalcthread, &nbccalcthread, &hist);
        else
            walk_kdtree_exact(cmd, gd, p, kd,
                              &nbbcalcthread, &nbccalcthread, &hist);

#ifdef SMOOTHPIVOT
        for (n = 1; n <= cmd->sizeHistN; n++) {
            hist.histNNSubXi2pcfthreadp[n] =
                        ((real)NbRmin(p))*hist.histNNSubXi2pcfthreadp[n];
            hist.histNNSubXi2pcfthreadtotal[n] +=
                        hist.histNNSubXi2pcfthreadp[n];
                hist.histNNSubthread[n] =
                    ((real)NbRmin(p))*hist.histNNSubthread[n];
        }
#endif

//B Normalization of histograms
        computeBodyProperties_sincos(cmd, gd, p, nbody[cat1], &hist);
//E

#ifdef SMOOTHPIVOT
        icountNbRminthread += NbRmin(p);
        icountNbRminOverlapthread += NbRminOverlap(p);
#endif
        INTEGER ip;
        ip = p - btab[cat1] + 1;
        if (ip%cmd->stepState == 0) {
            verb_log_print(cmd->verbose_log, gd->outlog,
                           " - Completed pivot: %ld\n", ip);
        }
    } // end do body p // end pragma omp DO_BODY p

        /* Publish in thread-ID order so floating-point sums are repeatable. */
#ifdef OPENMPCODE
        int thread_id = omp_get_thread_num();
        int thread_count = omp_get_num_threads();
#else
        int thread_id = 0;
        int thread_count = 1;
#endif
        for (int thread_turn = 0; thread_turn < thread_count; thread_turn++) {
#pragma omp barrier
        if (thread_id == thread_turn && hist_ready && !allocation_failed) {
            for (n = 1; n <= cmd->sizeHistN; n++) {
                gd->histNN[n] += hist.histNthread[n];
                gd->histNNSub[n] += hist.histNNSubthread[n];
                gd->histNNSubXi2pcf[n] += hist.histNNSubXi2pcfthread[n];
#ifdef SMOOTHPIVOT
                gd->histNNSubXi2pcftotal[n] += hist.histNNSubXi2pcfthreadtotal[n];
#endif
                gd->histXi2pcf[n] += hist.histXi2pcfthread[n];
            }
#ifdef TPCF
                int m;
                for (m=1; m<=cmd->mChebyshev+1; m++) {
                    ADDM_ext(gd->histZetaMcos[m],gd->histZetaMcos[m],
                             hist.histZetaMthreadcos[m],cmd->sizeHistN);
                    ADDM_ext(gd->histZetaMsin[m],gd->histZetaMsin[m],
                             hist.histZetaMthreadsin[m],cmd->sizeHistN);
                    ADDM_ext(gd->histZetaMsincos[m],gd->histZetaMsincos[m],
                             hist.histZetaMthreadsincos[m],cmd->sizeHistN);
                    // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
                    ADDM_ext(gd->histZetaMcossin[m],gd->histZetaMcossin[m],
                             hist.histZetaMthreadcossin[m],cmd->sizeHistN);
                }
#endif
            gd->nbbcalc += nbbcalcthread;
            gd->nbccalc += nbccalcthread;
#ifdef SMOOTHPIVOT
            ipfalse += ipfalsethreads;
            icountNbRmin += icountNbRminthread;
            icountNbRminOverlap += icountNbRminOverlapthread;
#endif
        }
        }
#pragma omp barrier

        if (hist_ready)
            search_free_sincos_omp(cmd, gd, &hist);
    } // end pragma omp parallel

    if (allocation_failed) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "searchcalc_kdtree_omp: OpenMP histogram allocation failed");
        finish_kdtree(kd);
        return FAILURE;
    }

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
                   "kdtree-omp: p falses found = %ld and %e %e %e\n",
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
        verb_print(cmd->verbose, "kdtree-omp: p falses found = %ld\n",ipfalse);
        //B kappa Avg Rmin
        verb_print(cmd->verbose,
                   "kdtree-omp: count NbRmin found = %ld\n",icountNbRmin);
        verb_print(cmd->verbose,
                   "kdtree-omp: count overlap found = %ld\n",icountNbRminOverlap);
        
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
        verb_print(cmd->verbose, "kdtree-omp: p falses found = %ld\n",ifalsecount);
        verb_print(cmd->verbose, "kdtree-omp: p true found = %ld\n",itruecount);
        verb_print(cmd->verbose, "kdtree-omp: total = %ld\n",itruecount+ifalsecount);
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
    bodyptr *bptr = kd->bptr;
    INTEGER cp = KDROOT;

    do {
        real dr1;
        real drpq2;
        compute_vector dr;

        Intersect(ntab[cp], gd->RcutSq, Pos(p), one_ball_next_cell);
        if (cp < kd->nsplit) {
            DOTPSUBV(drpq2, dr, Pos(p), ntab[cp].cmpos);
            dr1 = rsqrt(drpq2);
            if ((Radius(p) + ntab[cp].bnd.radius)/dr1 < gd->deltaR) {
                sumnode_sincos_cell(cmd, gd, p, ntab[cp], bptr,
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

    INTEGER pj;

    for (pj = ntab.first; pj <= ntab.last; ++pj) {
#ifdef SINGLEP
        q = &kd->packed_points[pj];
#else
        q = kd->bptr[pj];
#endif
        if (kdtree_accept_body(cmd, gd, p, q, &dr1, dr)) {
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
                               ballnode ntab, bodyptr *bptr,
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

    verb_print(cmd->verbose, "Search: Running ... (kdtree-omp) \n");

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
