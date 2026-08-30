/* ==============================================================================
 MODULE: search_kdtree_box_omp.c        [cTreeBalls]
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
#define KD_COORD_DELTA(a, b) ((real)(a) - (real)(b))

//INTERSECT: macro to determine if node intersects search ball
#if NDIM == 3
#define Intersect_bp(node, r2ball, pos, lbox, lbox_h, done, Continue) \
{                                                           \
    real _dl, _dh, _dl2, _dh2;                              \
    _dl2 = 0; _dh2 = 0;                                     \
    if (pos[0] < node.bnd.minb[0]) {                        \
        _dl = KD_COORD_DELTA(node.bnd.minb[0], pos[0]);     \
        _dh = KD_COORD_DELTA(node.bnd.maxb[0], pos[0]);     \
    } else if (pos[0] > node.bnd.maxb[0]) {                 \
        _dl = KD_COORD_DELTA(pos[0], node.bnd.maxb[0]);     \
        _dh = KD_COORD_DELTA(pos[0], node.bnd.minb[0]);     \
    } else {                                                \
        _dl = 0.0;                                          \
        _dh = MAX(KD_COORD_DELTA(node.bnd.maxb[0], pos[0]), \
                  KD_COORD_DELTA(pos[0], node.bnd.minb[0])); \
    }                                                       \
    if(_dl > lbox_h) {                                      \
        _dl = lbox - _dh;                                   \
        _dh = lbox - _dl;                                   \
    } else if(_dh > lbox_h) goto Continue;                  \
    _dl2 += _dl*_dl;                                        \
    _dh2 += _dh*_dh;                                        \
                                                            \
    if (pos[1] < node.bnd.minb[1]) {                        \
        _dl = KD_COORD_DELTA(node.bnd.minb[1], pos[1]);     \
        _dh = KD_COORD_DELTA(node.bnd.maxb[1], pos[1]);     \
    } else if (pos[1] > node.bnd.maxb[1]) {                 \
        _dl = KD_COORD_DELTA(pos[1], node.bnd.maxb[1]);     \
        _dh = KD_COORD_DELTA(pos[1], node.bnd.minb[1]);     \
    } else {                                                \
        _dl = 0.0;                                          \
        _dh = MAX(KD_COORD_DELTA(node.bnd.maxb[1], pos[1]), \
                  KD_COORD_DELTA(pos[1], node.bnd.minb[1])); \
    }                                                       \
    if(_dl > lbox_h) {                                      \
        _dl = lbox - _dh;                                   \
        _dh = lbox - _dl;                                   \
    } else if(_dh > lbox_h) goto Continue;                  \
    _dl2 += _dl*_dl;                                        \
    _dh2 += _dh*_dh;                                        \
                                                            \
    if (pos[2] < node.bnd.minb[2]) {                        \
        _dl = KD_COORD_DELTA(node.bnd.minb[2], pos[2]);     \
        _dh = KD_COORD_DELTA(node.bnd.maxb[2], pos[2]);     \
    } else if (pos[2] > node.bnd.maxb[2]) {                 \
        _dl = KD_COORD_DELTA(pos[2], node.bnd.maxb[2]);     \
        _dh = KD_COORD_DELTA(pos[2], node.bnd.minb[2]);     \
    } else {                                                \
        _dl = 0.0;                                          \
        _dh = MAX(KD_COORD_DELTA(node.bnd.maxb[2], pos[2]), \
                  KD_COORD_DELTA(pos[2], node.bnd.minb[2])); \
    }                                                       \
    if(_dl > lbox_h) {                                      \
        _dl = lbox - _dh;                                   \
        _dh = lbox - _dl;                                   \
    } else if(_dh > lbox_h) goto Continue;                  \
    _dl2 += _dl*_dl;                                        \
    _dh2 += _dh*_dh;                                        \
                                                            \
    if (_dl2 > r2ball) goto done;                           \
}

#define Intersect(node, r2ball, pos, done)                  \
{                                                           \
    real _dxl, _dxr, _dyl, _dyr, _dzl, _dzr, _dr2;          \
    _dxl = KD_COORD_DELTA(node.bnd.minb[0], pos[0]);        \
    _dxr = KD_COORD_DELTA(pos[0], node.bnd.maxb[0]);        \
    if (_dxl > 0.0) {                                       \
        _dr2 = _dxl*_dxl;                                   \
        if (_dr2 > r2ball) goto done;                       \
    } else if (_dxr > 0.0) {                                \
        _dr2 = _dxr*_dxr;                                   \
        if (_dr2 > r2ball) goto done;                       \
    } else                                                  \
        _dr2 = 0.0;                                         \
    _dyl = KD_COORD_DELTA(node.bnd.minb[1], pos[1]);        \
    _dyr = KD_COORD_DELTA(pos[1], node.bnd.maxb[1]);        \
    if (_dyl > 0.0) {                                       \
        _dr2 += _dyl*_dyl;                                  \
        if (_dr2 > r2ball) goto done;                       \
    } else if (_dyr > 0.0) {                                \
        _dr2 += _dyr*_dyr;                                  \
        if (_dr2 > r2ball) goto done;                       \
    }                                                       \
    _dzl = KD_COORD_DELTA(node.bnd.minb[2], pos[2]);        \
    _dzr = KD_COORD_DELTA(pos[2], node.bnd.maxb[2]);        \
    if (_dzl > 0.0) {                                       \
        _dr2 += _dzl*_dzl;                                  \
        if (_dr2 > r2ball) goto done;                       \
    } else if (_dzr > 0.0) {                                \
        _dr2 += _dzr*_dzr;                                  \
        if (_dr2 > r2ball) goto done;                       \
    }                                                       \
}

#else
#define Intersect_bp(node, r2ball, pos, lbox, lbox_h, done, Continue) \
{                                                           \
    real _dl, _dh, _dl2, _dh2;                              \
    _dl2 = 0; _dh2 = 0;                                     \
    if (pos[0] < node.bnd.minb[0]) {                        \
        _dl = KD_COORD_DELTA(node.bnd.minb[0], pos[0]);     \
        _dh = KD_COORD_DELTA(node.bnd.maxb[0], pos[0]);     \
    } else if (pos[0] > node.bnd.maxb[0]) {                 \
        _dl = KD_COORD_DELTA(pos[0], node.bnd.maxb[0]);     \
        _dh = KD_COORD_DELTA(pos[0], node.bnd.minb[0]);     \
    } else {                                                \
        _dl = 0.0;                                          \
        _dh = MAX(KD_COORD_DELTA(node.bnd.maxb[0], pos[0]), \
                  KD_COORD_DELTA(pos[0], node.bnd.minb[0])); \
    }                                                       \
    if(_dl > lbox_h) {                                      \
        _dl = lbox - _dh;                                   \
        _dh = lbox - _dl;                                   \
    } else if(_dh > lbox_h) goto Continue;                  \
    _dl2 += _dl*_dl;                                        \
    _dh2 += _dh*_dh;                                        \
                                                            \
    if (pos[1] < node.bnd.minb[1]) {                        \
        _dl = KD_COORD_DELTA(node.bnd.minb[1], pos[1]);     \
        _dh = KD_COORD_DELTA(node.bnd.maxb[1], pos[1]);     \
    } else if (pos[1] > node.bnd.maxb[1]) {                 \
        _dl = KD_COORD_DELTA(pos[1], node.bnd.maxb[1]);     \
        _dh = KD_COORD_DELTA(pos[1], node.bnd.minb[1]);     \
    } else {                                                \
        _dl = 0.0;                                          \
        _dh = MAX(KD_COORD_DELTA(node.bnd.maxb[1], pos[1]), \
                  KD_COORD_DELTA(pos[1], node.bnd.minb[1])); \
    }                                                       \
    if(_dl > lbox_h) {                                      \
        _dl = lbox - _dh;                                   \
        _dh = lbox - _dl;                                   \
    } else if(_dh > lbox_h) goto Continue;                  \
    _dl2 += _dl*_dl;                                        \
    _dh2 += _dh*_dh;                                        \
                                                            \
    if (_dl2 > r2ball) goto done;                           \
}

#define Intersect(node, r2ball, pos, done)                  \
{                                                           \
    real _dxl, _dxr, _dyl, _dyr, _dr2;                      \
    _dxl = KD_COORD_DELTA(node.bnd.minb[0], pos[0]);        \
    _dxr = KD_COORD_DELTA(pos[0], node.bnd.maxb[0]);        \
    if (_dxl > 0.0) {                                       \
        _dr2 = _dxl*_dxl;                                   \
        if (_dr2 > r2ball) goto done;                       \
    } else if (_dxr > 0.0) {                                \
        _dr2 = _dxr*_dxr;                                   \
        if (_dr2 > r2ball) goto done;                       \
    } else                                                  \
        _dr2 = 0.0;                                         \
    _dyl = KD_COORD_DELTA(node.bnd.minb[1], pos[1]);        \
    _dyr = KD_COORD_DELTA(pos[1], node.bnd.maxb[1]);        \
    if (_dyl > 0.0) {                                       \
        _dr2 += _dyl*_dyl;                                  \
        if (_dr2 > r2ball) goto done;                       \
    } else if (_dyr > 0.0) {                                \
        _dr2 += _dyr*_dyr;                                  \
        if (_dr2 > r2ball) goto done;                       \
    }                                                       \
}
#endif
//E

local void sumnode_sincos(struct  cmdline_data*, struct  global_data*,
                          bodyptr, ballnode, bodyptr *,
                          INTEGER *, INTEGER *,
                          gdhistptr_sincos_omp);
local int print_info(struct cmdline_data* cmd,
                     struct  global_data* gd);

/*
 Search routine using kdtree method:

 To be called using: search=kdtree-box-omp

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
global int searchcalc_kdtree_box_omp(struct cmdline_data* cmd,
                                    struct  global_data* gd,
                                    bodyptr *btab, INTEGER *nbody,
                                    INTEGER ipmin, INTEGER *ipmax,
                                    int cat1, int cat2)
{
    string routineName = "searchcalc_kdtree_box_omp";
    bodyptr p;
    int n;
    double cpustart;
    ballxptr kd;
    int nbucket;
    real cpu_build_kdtree;

    cpustart = CPUTIME;
    if (print_info(cmd, gd) == FAILURE)
        return FAILURE;

#ifdef OPENMPCODE
    ThreadCount(cmd, gd, nbody[cat1], cat1);
#endif

    search_init_gd_hist_sincos(cmd, gd);

    INTEGER ipfalse;
    ipfalse=0;
    int allocation_failed = FALSE;

    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                           "\n%s: Box sizes: %g %g %g and Lbox: %g\n",
                           routineName,
                           gd->Box[0], gd->Box[1], gd->Box[2], cmd->lengthBox);
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                           "\n%s: Total allocated %g MByte storage so far.\n",
                           routineName, gd->bytes_tot*INMB);
    verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                        "\nRunning...\n - Completed pivot node:\n");

    DO_BODY(p,btab[cat1]+ipmin-1, btab[cat1]+ipmax[cat1])
        Update(p) = TRUE;
//B Building kd-tree
    cpu_build_kdtree = CPUTIME;
    //B version 1.0.1
    nbucket = cmd->nsmooth;
    //E
    verb_print(cmd->verbose, "\nkdtree build: nbucket = %d\n",nbucket);
    kd = init_kdtree(cmd, gd, btab[cat2], nbody[cat2]);
    build_kdtree(cmd, gd, kd, nbucket);
    verb_print(cmd->verbose, "kdtree build: CPU time = %lf\n",
               CPUTIME-cpu_build_kdtree);
//E
    gd->ncellTable[cat1] = kd->nnode;               // Equivalent of octree cells

#pragma omp parallel default(none)   \
    shared(cmd,gd,btab,nbody,roottable,ipmin,ipmax, \
    rootnode, cat1, cat2, kd, ipfalse, allocation_failed)
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

        ballnode *ntab = kd->ntab;
        bodyptr *bptr = kd->bptr;
        INTEGER cp;

        real dr1;
        compute_vector dr;
        real drpq2;
        real lbox = cmd->lengthBox;
        real lbox_h = 0.5*lbox;

#pragma omp for nowait schedule(static,1)
    DO_BODY(p, btab[cat1]+ipmin-1, btab[cat1]+ipmax[cat1]) {
        if (allocation_failed) continue;
        for (n = 1; n <= cmd->sizeHistN; n++) {
            hist.histXi2pcfthreadsub[n] = 0.0;      // Affects only 2pcf
        }

        //B Walking the kdtree
        cp = KDROOT;
        do {
            Intersect_bp(ntab[cp], gd->RcutSq, Pos(p),
                         lbox, lbox_h, GetNextCell, Continue);
        Continue:
            if (cp < kd->nsplit) {
                cp = Lower(cp);
                continue;
            } else {
                sumnode_sincos(cmd, gd, p, ntab[cp], bptr,
                               &nbbcalcthread, &nbccalcthread, &hist);
            } // ! cp < nsplit
            GetNextCell:
            SetNext(cp);
        } while (cp != KDROOT);
        //E Walking the kdtree

        computeBodyProperties_sincos(cmd, gd, p, nbody[cat1], &hist);
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
                gd->histNNSubXi2pcf[n] += hist.histNNSubXi2pcfthread[n];
                gd->histXi2pcf[n] += hist.histXi2pcfthread[n];
            }
            gd->nbbcalc += nbbcalcthread;
            gd->nbccalc += nbccalcthread;
        }
        }
#pragma omp barrier

        if (hist_ready)
            search_free_sincos_omp(cmd, gd, &hist);

    } // end pragma omp parallel

    if (allocation_failed) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: OpenMP histogram allocation failed", routineName);
        finish_kdtree(kd);
        return FAILURE;
    }

    //B Normalization of histograms
    for (n = 1; n <= cmd->sizeHistN; n++) {
        gd->histXi2pcf[n] /= 2.0;
        gd->histNNSubXi2pcf[n] /= 2.0;
#ifdef SMOOTHPIVOT
        gd->histNNSubXi2pcftotal[n] /= 2.0;
#endif
        gd->histXi2pcf[n] /= MAX(gd->histNNSubXi2pcf[n],1.0);
    }
    //E

    if (cballs_opt_compute_histn(cmd)) {
        search_compute_HistN(cmd, gd, nbody[cat1]);
    }

    gd->cpusearch = CPUTIME - cpustart;
    verb_print(cmd->verbose, "Going out: CPU time = %lf\n",CPUTIME-cpustart);

    return SUCCESS;
}

local void sumnode_sincos(struct  cmdline_data* cmd,
                          struct  global_data* gd, bodyptr p,
                          ballnode ntab, bodyptr *bptr,
                          INTEGER *nbbcalcthread, INTEGER *nbccalcthread,
                          gdhistptr_sincos_omp hist)
{
    bodyptr q;
    real dr1;
    compute_vector dr;
    int n;
    real xi;

    INTEGER pj;

    for (pj = ntab.first; pj <= ntab.last; ++pj) {
        q = bptr[pj];
        if (accept_body(cmd, gd, p, (nodeptr)q, &dr1, dr)) {
            if(dr1>cmd->rminHist) {
                n = (int) ( (dr1-cmd->rminHist) * gd->i_deltaR) + 1;
                if (n<=cmd->sizeHistN && n>=1) {
                    hist->histNthread[n] = hist->histNthread[n] + 1.;
                    hist->histNNSubXi2pcfthread[n] =
                    hist->histNNSubXi2pcfthread[n] + 1.;
                    xi = Kappa(q);
                    hist->histXi2pcfthreadsub[n] += xi;
                    *nbbcalcthread += 1;
                }
            }
        } // ! accept_body
    } // ! loop first to last
}

local int print_info(struct cmdline_data* cmd,
                                  struct  global_data* gd)
{
    verb_print(cmd->verbose, "Search: Running ... (kdtree-box-omp) \n");

    if (cballs_opt_behavior_ball(cmd)) {
        verb_print(cmd->verbose, "with option behavior-ball... \n");
        if (!cmd->useLogHist) {
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     "print_info: behavior-ball and useLogHist=false are incompatible");
            return FAILURE;
        }
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
