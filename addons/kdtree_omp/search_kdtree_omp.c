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
#if NDIM == 3
#define Intersect(node, r2ball, pos, done)                  \
{                                                           \
    real _dxl, _dxr, _dyl, _dyr, _dzl, _dzr, _dr2;          \
    _dxl = node.bnd.minb[0] - pos[0];                       \
    _dxr = pos[0] - node.bnd.maxb[0];                       \
    if (_dxl > 0.0) {                                       \
        _dr2 = _dxl*_dxl;                                   \
        if (_dr2 > r2ball) goto done;                       \
    } else if (_dxr > 0.0) {                                \
        _dr2 = _dxr*_dxr;                                   \
        if (_dr2 > r2ball) goto done;                       \
    } else                                                  \
        _dr2 = 0.0;                                         \
    _dyl = node.bnd.minb[1] - pos[1];                       \
    _dyr = pos[1] - node.bnd.maxb[1];                       \
    if (_dyl > 0.0) {                                       \
        _dr2 += _dyl*_dyl;                                  \
        if (_dr2 > r2ball) goto done;                       \
    } else if (_dyr > 0.0) {                                \
        _dr2 += _dyr*_dyr;                                  \
        if (_dr2 > r2ball) goto done;                       \
    }                                                       \
    _dzl = node.bnd.minb[2] - pos[2];                       \
    _dzr = pos[2] - node.bnd.maxb[2];                       \
    if (_dzl > 0.0) {                                       \
        _dr2 += _dzl*_dzl;                                  \
        if (_dr2 > r2ball) goto done;                       \
    } else if (_dzr > 0.0) {                                \
        _dr2 += _dzr*_dzr;                                  \
        if (_dr2 > r2ball) goto done;                       \
    }                                                       \
}

#else
#define Intersect(node, r2ball, pos, done)                  \
{                                                           \
    real _dxl, _dxr, _dyl, _dyr, _dr2;                      \
    _dxl = node.bnd.minb[0] - pos[0];                       \
    _dxr = pos[0] - node.bnd.maxb[0];                       \
    if (_dxl > 0.0) {                                       \
        _dr2 = _dxl*_dxl;                                   \
        if (_dr2 > r2ball) goto done;                       \
    } else if (_dxr > 0.0) {                                \
        _dr2 = _dxr*_dxr;                                   \
        if (_dr2 > r2ball) goto done;                       \
    } else                                                  \
        _dr2 = 0.0;                                         \
    _dyl = node.bnd.minb[1] - pos[1];                       \
    _dyr = pos[1] - node.bnd.maxb[1];                       \
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

#ifdef NMultipoles
local void sumnode_sincos_N(struct  cmdline_data*, struct  global_data*,
                          bodyptr, kdnode, bodyptr *,
                          INTEGER *, INTEGER *,
                          gdhistptr_sincos_omp, gdhistNptr_sincos_omp);
local void sumnode_sincos_cell_N(struct  cmdline_data*,
                               struct  global_data*, bodyptr,
                               kdnode, bodyptr *,
                               INTEGER *, INTEGER *,
                               gdhistptr_sincos_omp, gdhistNptr_sincos_omp);
#endif
local void sumnode_sincos(struct  cmdline_data*, struct  global_data*,
                          bodyptr, kdnode, bodyptr *,
                          INTEGER *, INTEGER *,
                          gdhistptr_sincos_omp);
local void sumnode_sincos_cell(struct  cmdline_data*,
                               struct  global_data*, bodyptr,
                               kdnode, bodyptr *,
                               INTEGER *, INTEGER *,
                               gdhistptr_sincos_omp);

// search=kdtree-omp
global int searchcalc_kdtree_omp(struct cmdline_data* cmd,
                                    struct  global_data* gd,
                                    bodyptr *btab, INTEGER *nbody,
                                    INTEGER ipmin, INTEGER *ipmax,
                                    int cat1, int cat2)
{
    bodyptr p;
    int n;
    double cpustart;
    kdxptr kd;
    int nbucket;
    real cpu_build_kdtree;

    cpustart = CPUTIME;
    verb_print(cmd->verbose, "Search: Running ... (kdtree-omp) \n");

    if (scanopt(cmd->options, "behavior-ball")) {
        verb_print(cmd->verbose, "with option behavior-ball... \n");
        if (!cmd->useLogHist)
            error("behavior-ball and useLogHist=false are incompatible!");
    }
    if (scanopt(cmd->options, "smooth-pivot"))
        verb_print(cmd->verbose,
                   "with option smooth-pivot... rsmooth=%g\n",gd->rsmooth[0]);
    if (!cmd->computeTPCF)
        verb_print(cmd->verbose, "computing only 2pcf... \n");
#ifdef NMultipoles
    verb_print(cmd->verbose, "with NMultipoles... \n");
#else
    verb_print(cmd->verbose, "without NMultipoles... \n");
#endif
#ifdef NONORMHIST
    verb_print(cmd->verbose, "with NONORMHIST... \n");
#else
    verb_print(cmd->verbose, "without NONORMHIST... \n");
#endif

#ifdef OPENMPCODE
    ThreadCount(cmd, gd, nbody[cat1], cat1);
#endif

    search_init_gd_hist_sincos(cmd, gd);
#ifdef NMultipoles
    search_init_gd_hist_sincos_N(cmd, gd);
#endif

    INTEGER ipfalse;
    ipfalse=0;
//B kappa Avg Rmin
    INTEGER icountNbRmin;
    icountNbRmin=0;
    INTEGER icountNbRminOverlap;
    icountNbRminOverlap=0;
//E

    verb_print(cmd->verbose,
               "\nsearchcalc_balls: Total allocated %g MByte storage so far.\n",
               gd->bytes_tot/(1024.0*1024.0));


    DO_BODY(p,btab[cat1]+ipmin-1, btab[cat1]+ipmax[cat1])
        Update(p) = TRUE;
//B Building kd-tree
    cpu_build_kdtree = CPUTIME;
    nbucket = gd->nsmooth[0];
    verb_print(cmd->verbose, "\nkdtree build: nbucket = %d\n",nbucket);
    kd = init_kdtree(cmd, gd, btab[cat2], nbody[cat2]);
    build_kdtree(cmd, gd, kd, nbucket);
    verb_print(cmd->verbose, "kdtree build: CPU time = %lf\n",
               CPUTIME-cpu_build_kdtree);
//E
    gd->ncellTable[cat1] = kd->nnode;               // Equivalent of octree cells

#pragma omp parallel default(none)   \
    shared(cmd,gd,btab,nbody,roottable,ipmin,ipmax, \
    rootnode, cat1, cat2, kd, ipfalse, icountNbRmin, icountNbRminOverlap)
    {
        bodyptr p;
        int n;
        INTEGER nbbcalcthread = 0;
        INTEGER nbccalcthread = 0;
        
        gdhist_sincos_omp hist;
        search_init_sincos_omp(cmd, gd, &hist);
#ifdef NMultipoles
        gdhistN_sincos_omp histN;
        search_init_sincos_omp_N(cmd, gd, &histN);
#endif

        INTEGER ipfalsethreads;
        ipfalsethreads = 0;

//B kappa Avg Rmin
        INTEGER icountNbRminthread;
        icountNbRminthread=0;
        INTEGER icountNbRminOverlapthread;
        icountNbRminOverlapthread=0;
//E
        kdnode *ntab = kd->ntab;
        bodyptr *bptr = kd->bptr;
        INTEGER cp;

#pragma omp for nowait schedule(dynamic)
    DO_BODY(p, btab[cat1]+ipmin-1, btab[cat1]+ipmax[cat1]) {
// p and q are in differents node structures... cat1!=cat2...
//B kappa Avg Rmin
        NbRmin(p) = 1;
        NbRminOverlap(p) = 0;
#ifdef NMultipoles
        KappaRminN(p) = 1.0;
        KappaRmin(p) = Kappa(p);
#else
        KappaRmin(p) = Kappa(p);
#endif
        if (scanopt(cmd->options, "smooth-pivot")) {
            if (Update(p) == FALSE) {
                ipfalsethreads++;
                continue;
            }
        }
//E
        for (n = 1; n <= cmd->sizeHistN; n++) {
            hist.histNNSubthread[n] = 0.0;          // Affects only 3pcf
            hist.histXi2pcfthreadsub[n] = 0.0;      // Affects only 2pcf
//B kappa Avg Rmin
            hist.histNNSubXi2pcfthreadp[n] = 0.;    // Affects only 2pcf
//E
        }
        if (cmd->computeTPCF) {
            CLRM_ext_ext(hist.histXithreadcos, cmd->mChebyshev+1, 
                         cmd->sizeHistN);
            CLRM_ext_ext(hist.histXithreadsin, cmd->mChebyshev+1, 
                         cmd->sizeHistN);
#if NDIM == 3
            dRotation3D(Pos(p), ROTANGLE, ROTANGLE, ROTANGLE, hist.q0);
            DOTPSUBV(hist.drpq2, hist.dr0, Pos(p), hist.q0);
#ifdef SINGLEP
            hist.drpq = sqrt(hist.drpq2);
#else
            hist.drpq = rsqrt(hist.drpq2);
#endif
            //B Random rotation of dr0:
#ifdef PTOPIVOTROTATION
            real rtheta;
            vector dr0rot;
            rtheta = xrandom(0.0, TWOPI);
            RotationVecAWRtoVecB(dr0rot, hist.dr0, Pos(p), rtheta);
            SETV(hist.dr0, dr0rot);
#endif
            //E
#endif // ! NDIM
        } // ! computeTPCF

#ifdef NMultipoles
        for (n = 1; n <= cmd->sizeHistN; n++) {
            histN.histNNSubthread[n] = 0.0;         // Affects only 3pcf
            histN.histXi2pcfthreadsub[n] = 0.0;     // Affects only 2pcf
//B kappa Avg Rmin
            histN.histNNSubXi2pcfthreadp[n] = 0.;   // Affects only 2pcf
//E
        }
        if (cmd->computeTPCF) {
            CLRM_ext_ext(histN.histXithreadcos, cmd->mChebyshev+1,
                         cmd->sizeHistN);
            CLRM_ext_ext(histN.histXithreadsin, cmd->mChebyshev+1,
                         cmd->sizeHistN);
#if NDIM == 3
            dRotation3D(Pos(p), ROTANGLE, ROTANGLE, ROTANGLE, histN.q0);
            DOTPSUBV(histN.drpq2, histN.dr0, Pos(p), histN.q0);
#ifdef SINGLEP
            histN.drpq = sqrt(histN.drpq2);
#else
            histN.drpq = rsqrt(histN.drpq2);
#endif
            //B Random rotation of dr0:
#ifdef PTOPIVOTROTATION
            real rtheta;
            vector dr0rot;
            rtheta = xrandom(0.0, TWOPI);
            RotationVecAWRtoVecB(dr0rot, histN.dr0, Pos(p), rtheta);
            SETV(histN.dr0, dr0rot);
#endif
            //E
#endif // ! NDIM
        } // ! computeTPCF
#endif // ! NMultipoles

        real dr1;
        vector dr;
        real drpq2;
        //B Walking the kdtree
        cp = KDROOT;
        do {
            Intersect(ntab[cp], gd->RcutSq, Pos(p), GetNextCell);
            if (scanopt(cmd->options, "behavior-ball")) {
                if (cp < kd->nsplit) {
                    DOTPSUBV(drpq2, dr, Pos(p), ntab[cp].cmpos);
                    dr1 = rsqrt(drpq2);
                    if ( (Radius(p)+ntab[cp].bnd.radius)/(dr1) < gd->deltaR) {
#ifdef NMultipoles
                        sumnode_sincos_cell_N(cmd, gd, p, ntab[cp], bptr,
                                            &nbbcalcthread, &nbccalcthread,
                                            &hist, &histN);
#else
                        sumnode_sincos_cell(cmd, gd, p, ntab[cp], bptr,
                                            &nbbcalcthread, &nbccalcthread,
                                            &hist);
#endif
                        cp = Upper(cp);
                        continue;
                    } else {
//                        cp = Lower(cp);
                        SetNext(cp);
                        continue;
                    }
                } else {
#ifdef NMultipoles
                    sumnode_sincos_N(cmd, gd, p, ntab[cp], bptr,
                                   &nbbcalcthread, &nbccalcthread,
                                   &hist, &histN);
#else
                    sumnode_sincos(cmd, gd, p, ntab[cp], bptr,
                                   &nbbcalcthread, &nbccalcthread, 
                                   &hist);
#endif
                } // ! cp < nsplit
            } else { // ! behavior-ball
                if (cp < kd->nsplit) {
                    cp = Lower(cp);
                    continue;
                } else {
#ifdef NMultipoles
                    sumnode_sincos_N(cmd, gd, p, ntab[cp], bptr,
                                   &nbbcalcthread, &nbccalcthread,
                                   &hist, &histN);
#else
                    sumnode_sincos(cmd, gd, p, ntab[cp], bptr,
                                   &nbbcalcthread, &nbccalcthread,
                                   &hist);
#endif
                } // ! cp < nsplit
            } // ! behavior-ball
            GetNextCell:
            SetNext(cp);
        } while (cp != KDROOT);
        //E Walking the kdtree

//B kappa Avg Rmin
        for (n = 1; n <= cmd->sizeHistN; n++) {
            hist.histNNSubXi2pcfthreadp[n] =
                        ((real)NbRmin(p))*hist.histNNSubXi2pcfthreadp[n];
            hist.histNNSubXi2pcfthreadtotal[n] +=
                        hist.histNNSubXi2pcfthreadp[n];
            if (scanopt(cmd->options, "smooth-pivot"))
                hist.histNNSubthread[n] =
                    ((real)NbRmin(p))*hist.histNNSubthread[n];
        }
//E
#ifdef NMultipoles
        //B kappa Avg Rmin
        for (n = 1; n <= cmd->sizeHistN; n++) {
            histN.histNNSubXi2pcfthreadp[n] =
                        ((real)NbRmin(p))*histN.histNNSubXi2pcfthreadp[n];
            histN.histNNSubXi2pcfthreadtotal[n] +=
                        histN.histNNSubXi2pcfthreadp[n];
            if (scanopt(cmd->options, "smooth-pivot"))
                histN.histNNSubthread[n] =
                        ((real)NbRmin(p))*histN.histNNSubthread[n];
        }
        //E
#endif

//B Normalization of histograms
#ifdef NMultipoles
        computeBodyProperties_sincos(cmd, gd, p, nbody[cat1], &hist);
        computeBodyProperties_sincos_N(cmd, gd, p, nbody[cat1], &histN);
#else
        computeBodyProperties_sincos(cmd, gd, p, nbody[cat1], &hist);
#endif
//E

//B kappa Avg Rmin
        icountNbRminthread += NbRmin(p);
        icountNbRminOverlapthread += NbRminOverlap(p);
//E
        INTEGER ip;
        ip = p - btab[cat1] + 1;
        if (ip%cmd->stepState == 0) {
            verb_log_print(cmd->verbose_log, gd->outlog,
                           " - Completed pivot: %ld\n", ip);
        }
    } // end do body p // end pragma omp DO_BODY p

#pragma omp critical
        {
            for (n = 1; n <= cmd->sizeHistN; n++) {
                gd->histNN[n] += hist.histNthread[n];
                gd->histNNSub[n] += hist.histNNSubthread[n];
                gd->histNNSubXi2pcf[n] += hist.histNNSubXi2pcfthread[n];
//B kappa Avg Rmin
                gd->histNNSubXi2pcftotal[n] += hist.histNNSubXi2pcfthreadtotal[n];
//E
                gd->histXi2pcf[n] += hist.histXi2pcfthread[n];
            }
            if (cmd->computeTPCF) {
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
            }
            gd->nbbcalc += nbbcalcthread;
            gd->nbccalc += nbccalcthread;
#ifdef NMultipoles
            for (n = 1; n <= cmd->sizeHistN; n++) {
                gd->NhistNN[n] += histN.histNthread[n];
                gd->NhistNNSub[n] += histN.histNNSubthread[n];
                gd->NhistNNSubXi2pcf[n] += histN.histNNSubXi2pcfthread[n];
//B kappa Avg Rmin
                gd->NhistNNSubXi2pcftotal[n] += histN.histNNSubXi2pcfthreadtotal[n];
//E
                gd->NhistXi2pcf[n] += histN.histXi2pcfthread[n];
            }
            if (cmd->computeTPCF) {
                int m;
                for (m=1; m<=cmd->mChebyshev+1; m++) {
                    ADDM_ext(gd->NhistZetaMcos[m],gd->NhistZetaMcos[m],
                             histN.histZetaMthreadcos[m],cmd->sizeHistN);
                    ADDM_ext(gd->NhistZetaMsin[m],gd->NhistZetaMsin[m],
                             histN.histZetaMthreadsin[m],cmd->sizeHistN);
                    ADDM_ext(gd->NhistZetaMsincos[m],gd->NhistZetaMsincos[m],
                             histN.histZetaMthreadsincos[m],cmd->sizeHistN);
                    // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
                    ADDM_ext(gd->NhistZetaMcossin[m],gd->NhistZetaMcossin[m],
                             histN.histZetaMthreadcossin[m],cmd->sizeHistN);
                }
            }
#endif
        }

        search_free_sincos_omp(cmd, gd, &hist);
#ifdef NMultipoles
        search_free_sincos_omp_N(cmd, gd, &histN);
#endif

//B kappa Avg Rmin
        ipfalse += ipfalsethreads;
        icountNbRmin += icountNbRminthread;
        icountNbRminOverlap += icountNbRminOverlapthread;
//E
    } // end pragma omp parallel

    real xi, den, num;
    int mm;
    if (scanopt(cmd->options, "smooth-pivot")) {
        num = (real)nbody[cat1];
//B kappa Avg Rmin
        den = (real)(nbody[cat1]-ipfalse);
//E
        xi = num/den;
        verb_print(cmd->verbose,
                   "kdtree-omp: p falses found = %ld and %e %e %e\n",
                   ipfalse, num, den, xi);
        if (cmd->computeTPCF) {
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
#ifdef NMultipoles
        real xiN;
        xiN = num/den;
        if (cmd->computeTPCF) {
            for (mm=1; mm<=cmd->mChebyshev+1; mm++) {
                MULMS_ext(gd->NhistZetaMcos[mm], gd->NhistZetaMcos[mm],
                          xiN,cmd->sizeHistN);
                MULMS_ext(gd->NhistZetaMsin[mm], gd->NhistZetaMsin[mm],
                          xiN,cmd->sizeHistN);
                MULMS_ext(gd->NhistZetaMsincos[mm], gd->NhistZetaMsincos[mm],
                          xiN,cmd->sizeHistN);
                // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
                MULMS_ext(gd->NhistZetaMcossin[mm], gd->NhistZetaMcossin[mm],
                          xiN,cmd->sizeHistN);
            }
        }
#endif
    } // ! smooth-pivot
    
#ifdef NONORMHIST
//B Normalization of histograms
    if (!scanopt(cmd->options, "asymmetric")) {
        for (n = 1; n <= cmd->sizeHistN; n++) {
            if (cmd->verbose>3)
                printf("%d %e %e\n", n,
                   gd->histNNSubXi2pcf[n], gd->histNNSubXi2pcftotal[n]);
            gd->histXi2pcf[n] /= 2.0;
            gd->histNNSubXi2pcf[n] /= 2.0;
//B kappa Avg Rmin
            gd->histNNSubXi2pcftotal[n] /= 2.0;
            if (scanopt(cmd->options, "smooth-pivot")) {
                gd->histXi2pcf[n] /= 1.0;
            } else {
                gd->histXi2pcf[n] /= 1.0;
            }
//E
        }
    } else {
        for (n = 1; n <= cmd->sizeHistN; n++) {
            if (cmd->verbose>3)
            printf("%d %e %e\n", n,
                   gd->histNNSubXi2pcf[n], gd->histNNSubXi2pcftotal[n]);
            if (scanopt(cmd->options, "smooth-pivot")) {
                gd->histXi2pcf[n] /= 1.0;
            } else {
                gd->histXi2pcf[n] /= 1.0;
            }
        }
    }
//E
#ifdef NMultipoles
    //B Normalization of histograms
        if (!scanopt(cmd->options, "asymmetric")) {
            for (n = 1; n <= cmd->sizeHistN; n++) {
                if (cmd->verbose>3)
                    printf("%d %e %e\n", n,
                       gd->NhistNNSubXi2pcf[n], gd->NhistNNSubXi2pcftotal[n]);
                gd->NhistXi2pcf[n] /= 2.0;
                gd->NhistNNSubXi2pcf[n] /= 2.0;
    //B kappa Avg Rmin
                gd->NhistNNSubXi2pcftotal[n] /= 2.0;
                if (scanopt(cmd->options, "smooth-pivot")) {
                    gd->NhistXi2pcf[n] /= 1.0;
                } else {
                    gd->NhistXi2pcf[n] /= 1.0;
                }
    //E
            }
        } else {
            for (n = 1; n <= cmd->sizeHistN; n++) {
                if (cmd->verbose>3)
                printf("%d %e %e\n", n,
                       gd->NhistNNSubXi2pcf[n], gd->NhistNNSubXi2pcftotal[n]);
                if (scanopt(cmd->options, "smooth-pivot")) {
                    gd->NhistXi2pcf[n] /= 1.0;
                } else {
                    gd->NhistXi2pcf[n] /= 1.0;
                }
            }
        }
    //E
#endif // ! NMultipoles
#else // ! NONORMHIST
    //B Normalization of histograms
        if (!scanopt(cmd->options, "asymmetric")) {
            for (n = 1; n <= cmd->sizeHistN; n++) {
                if (cmd->verbose>3)
                    printf("%d %e %e\n", n,
                       gd->histNNSubXi2pcf[n], gd->histNNSubXi2pcftotal[n]);
                gd->histXi2pcf[n] /= 2.0;
                gd->histNNSubXi2pcf[n] /= 2.0;
    //B kappa Avg Rmin
                gd->histNNSubXi2pcftotal[n] /= 2.0;
                if (scanopt(cmd->options, "smooth-pivot")) {
                    gd->histXi2pcf[n] /= MAX(gd->histNNSubXi2pcftotal[n],1.0);
                } else {
                    gd->histXi2pcf[n] /= MAX(gd->histNNSubXi2pcf[n],1.0);
                }
    //E
            }
        } else {
            for (n = 1; n <= cmd->sizeHistN; n++) {
                if (cmd->verbose>3)
                printf("%d %e %e\n", n,
                       gd->histNNSubXi2pcf[n], gd->histNNSubXi2pcftotal[n]);
                if (scanopt(cmd->options, "smooth-pivot")) {
                    gd->histXi2pcf[n] /= MAX(gd->histNNSubXi2pcftotal[n],1.0);
                } else {
                    gd->histXi2pcf[n] /= MAX(gd->histNNSubXi2pcf[n],1.0);
                }
            }
        }
    //E
#ifdef NMultipoles
    //B Normalization of histograms
        if (!scanopt(cmd->options, "asymmetric")) {
            for (n = 1; n <= cmd->sizeHistN; n++) {
                if (cmd->verbose>3)
                    printf("%d %e %e\n", n,
                       gd->NhistNNSubXi2pcf[n], gd->NhistNNSubXi2pcftotal[n]);
                gd->NhistXi2pcf[n] /= 2.0;
                gd->NhistNNSubXi2pcf[n] /= 2.0;
    //B kappa Avg Rmin
                gd->NhistNNSubXi2pcftotal[n] /= 2.0;
                if (scanopt(cmd->options, "smooth-pivot")) {
                    gd->NhistXi2pcf[n] /= MAX(gd->NhistNNSubXi2pcftotal[n],1.0);
                } else {
                    gd->NhistXi2pcf[n] /= MAX(gd->NhistNNSubXi2pcf[n],1.0);
                }
    //E
            }
        } else {
            for (n = 1; n <= cmd->sizeHistN; n++) {
                if (cmd->verbose>3)
                printf("%d %e %e\n", n,
                       gd->NhistNNSubXi2pcf[n], gd->NhistNNSubXi2pcftotal[n]);
                if (scanopt(cmd->options, "smooth-pivot")) {
                    gd->NhistXi2pcf[n] /= MAX(gd->NhistNNSubXi2pcftotal[n],1.0);
                } else {
                    gd->NhistXi2pcf[n] /= MAX(gd->NhistNNSubXi2pcf[n],1.0);
                }
            }
        }
    //E
#endif // ! NMultipoles
#endif // ! NONORMHIST

    if (scanopt(cmd->options, "compute-HistN")) {
        if (scanopt(cmd->options, "smooth-pivot")) {
            search_compute_HistN(cmd, gd, nbody[cat1]-ipfalse);
        } else {
            search_compute_HistN(cmd, gd, nbody[cat1]);
        }
    }

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

    gd->cpusearch = CPUTIME - cpustart;
    verb_print(cmd->verbose, "Going out: CPU time = %lf\n",CPUTIME-cpustart);

    return SUCCESS;
}


#ifdef NMultipoles
local void sumnode_sincos_N(struct  cmdline_data* cmd,
                          struct  global_data* gd, bodyptr p,
                          kdnode ntab, bodyptr *bptr,
                          INTEGER *nbbcalcthread, INTEGER *nbccalcthread,
                          gdhistptr_sincos_omp hist,
                          gdhistNptr_sincos_omp histN)
{
    bodyptr q;
#ifdef SINGLEP
    float dr1;
    float dr[NDIM];
#else
    real dr1;
    vector dr;
#endif
    int n;
    real xiN;                                       // NMultipoles
    real xi;

    INTEGER pj;

//#ifdef OPENMPCODE
//  #pragma omp parallel for
//#endif
    for (pj = ntab.first; pj <= ntab.last; ++pj) {
        q = bptr[pj];
        if (accept_body(cmd, gd, p, (nodeptr)q, &dr1, dr)) {
            if (scanopt(cmd->options, "smooth-pivot")) {
                if (dr1<=gd->rsmooth[0]) {
                    //B kappa Avg Rmin
                    if (Update(q)==TRUE) {
                        Update(q) = FALSE;
                        NbRmin(p) += 1;
                        // Should be another KappaRmin to account to XiMultipoles,
                        //  like the one after "else"
                        KappaRminN(p) += 1;         // NMultipoles
                        KappaRmin(p) += Kappa(q);
                    } else {
                        NbRminOverlap(p) += 1;
                    }
                }
            }
//E
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
                        //B kappa Avg Rmin
                        hist->histNNSubXi2pcfthreadp[n] =
                                            hist->histNNSubXi2pcfthreadp[n] + 1.;
                        //E
                        hist->histNNSubthread[n] = hist->histNNSubthread[n] + 1.;
                        xi = Kappa(q);

                        //B NMultipoles
                        histN->histNthread[n] = histN->histNthread[n] + 1.;
                        histN->histNNSubXi2pcfthread[n] =
                                            histN->histNNSubXi2pcfthread[n] + 1.;
                        //B kappa Avg Rmin
                        histN->histNNSubXi2pcfthreadp[n] =
                                            histN->histNNSubXi2pcfthreadp[n] + 1.;
                        //E
                        histN->histNNSubthread[n] = histN->histNNSubthread[n] + 1.;
                        xiN = 1.0;
                        //E

                        if (cmd->computeTPCF) {
                            REAL cosphi,sinphi;
#if NDIM == 3
                            REAL s, sy;
#ifdef SINGLEP
                            float pr0[NDIM];
#else
                            vector pr0;
#endif
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
                            NCHEBYSHEVTUOMPSINCOS;  // NMultipoles
                            CHEBYSHEVTUOMPSINCOS;
                        } // ! computeTPCF
                        histN->histXi2pcfthreadsub[n] += xiN;   // NMultipoles
                        hist->histXi2pcfthreadsub[n] += xi;
                        *nbbcalcthread += 1;
                    } // ! 1 < n < sizeHistN
                } // ! dr1 > rminHist
            } else { // ! useLogHist
                if(dr1>cmd->rminHist) {
                    n = (int) ( (dr1-cmd->rminHist) * gd->i_deltaR) + 1;
                    if (n<=cmd->sizeHistN && n>=1) {
                        hist->histNthread[n] = hist->histNthread[n] + 1.;
                        hist->histNNSubXi2pcfthread[n] =
                        hist->histNNSubXi2pcfthread[n] + 1.;
                        //B kappa Avg Rmin
                        hist->histNNSubXi2pcfthreadp[n] =
                                            hist->histNNSubXi2pcfthreadp[n] + 1.;
                        //E
                        hist->histNNSubthread[n] = hist->histNNSubthread[n] + 1.;
                        xi = Kappa(q);
                        //B NMultipoles
                        histN->histNthread[n] = histN->histNthread[n] + 1.;
                        histN->histNNSubXi2pcfthread[n] =
                                            histN->histNNSubXi2pcfthread[n] + 1.;
                        //B kappa Avg Rmin
                        histN->histNNSubXi2pcfthreadp[n] =
                                            histN->histNNSubXi2pcfthreadp[n] + 1.;
                        //E
                        histN->histNNSubthread[n] = histN->histNNSubthread[n] + 1.;
                        xiN = 1.0;
                        //E
                        if (cmd->computeTPCF) {
                            real cosphi,sinphi;
#if NDIM == 3
                            real s, sy;
//B Add switch SINGLEP here and below where apply
                            vector pr0;
//E
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
                            NCHEBYSHEVTUOMPSINCOS;  // NMultipoles
                            CHEBYSHEVTUOMPSINCOS;
                        }
                        histN->histXi2pcfthreadsub[n] += xiN;   // NMultipoles
                        hist->histXi2pcfthreadsub[n] += xi;
                        *nbbcalcthread += 1;
                    } // ! 1 < n < sizeHistN
                } // ! dr1 > rminHist
            } // ! useLogHist
        } // ! accept_body
    } // ! loop first to last
}


local void sumnode_sincos_cell_N(struct  cmdline_data* cmd,
                               struct  global_data* gd, bodyptr p,
                               kdnode ntab, bodyptr *bptr,
                               INTEGER *nbbcalcthread, INTEGER *nbccalcthread,
                               gdhistptr_sincos_omp hist,
                               gdhistNptr_sincos_omp histN)
{
//    bodyptr q;
#ifdef SINGLEP
    float dr1;
    float dr[NDIM];
#else
    real dr1;
    real drpq2;
    vector dr;
#endif
    int n;
    real xi;

    int npoints;
    npoints = ntab.last - ntab.first + 1;
//#ifdef OPENMPCODE
//  #pragma omp parallel for
//#endif
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
                                            hist->histNNSubXi2pcfthread[n] + 1.0;
                        hist->histNNSubthread[n] = hist->histNNSubthread[n] + 1.0;
                        //B NMultipoles
                        histN->histNthread[n] = histN->histNthread[n] +  npoints;
                        histN->histNNSubXi2pcfthread[n] =
                                            histN->histNNSubXi2pcfthread[n] + 1.0;
                        histN->histNNSubthread[n] = histN->histNNSubthread[n] + 1.0;
                        //E
#ifdef KappaAvgON
                        real xiN;
                        xiN = 1.0;                  // NMultipoles
                        xi = ntab.kappa;
#else
                        real xiN;                   // NMultipoles
                        xiN = 1.0;
                        xi = ntab.kappa;
#endif

                        if (cmd->computeTPCF) {
                        REAL cosphi,sinphi;
#if NDIM == 3
                        REAL s, sy;
#ifdef SINGLEP
                        float pr0[NDIM];
#else
                        vector pr0;
#endif
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
                        NCHEBYSHEVTUOMPSINCOS;      // NMultipoles
                        CHEBYSHEVTUOMPSINCOS;
                        } // ! computeTPCF
                        histN->histXi2pcfthreadsub[n] += xiN;
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
                                            hist->histNNSubXi2pcfthread[n] + 1.0;
                        hist->histNNSubthread[n] = hist->histNNSubthread[n] + 1.0;
                        //B NMultipoles
                        histN->histNthread[n] = histN->histNthread[n] +  npoints;
                        histN->histNNSubXi2pcfthread[n] =
                                            histN->histNNSubXi2pcfthread[n] + 1.0;
                        histN->histNNSubthread[n] = histN->histNNSubthread[n] + 1.0;
                        //E
#ifdef KappaAvgON
                        real xiN;                   // NMultipoles
                        xiN = 1.0;
                        xi = ntab.kappa;
#else
                        real xiN;                   // NMultipoles
                        xiN = 1.0;
                        xi = ntab.kappa;
#endif
                        if (cmd->computeTPCF) {
                            real cosphi, sinphi;
#if NDIM == 3
                            real s, sy;
                            vector pr0;
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
                            NCHEBYSHEVTUOMPSINCOS;      // NMultipoles
                            CHEBYSHEVTUOMPSINCOS;
                        } // computeTPCF
                        histN->histXi2pcfthreadsub[n] += xiN;       // NMultipoles
                        hist->histXi2pcfthreadsub[n] += xi;
                        *nbccalcthread += 1;
                    } // ! n in (1, sizeHistN)
                } // ! dr1 > rminHist
            } // ! useLogHist
        } // ! dr1 < rangeN
//        } // ! accept_body
//    } // ! loop pj
}

#endif // ! NMultipoles


local void sumnode_sincos(struct  cmdline_data* cmd,
                          struct  global_data* gd, bodyptr p,
                          kdnode ntab, bodyptr *bptr,
                          INTEGER *nbbcalcthread, INTEGER *nbccalcthread,
                          gdhistptr_sincos_omp hist)
{
    bodyptr q;
#ifdef SINGLEP
    float dr1;
    float dr[NDIM];
#else
    real dr1;
    vector dr;
#endif
    int n;
    real xi;

    INTEGER pj;

//#ifdef OPENMPCODE
//  #pragma omp parallel for
//#endif
    for (pj = ntab.first; pj <= ntab.last; ++pj) {
        q = bptr[pj];
        if (accept_body(cmd, gd, p, (nodeptr)q, &dr1, dr)) {
            if (scanopt(cmd->options, "smooth-pivot"))
                if (dr1<=gd->rsmooth[0]) {
//B kappa Avg Rmin
                    if (Update(q)==TRUE) {
                        Update(q) = FALSE;
                        NbRmin(p) += 1;
                        KappaRmin(p) += Kappa(q);
                    } else {
                        NbRminOverlap(p) += 1;
                    }
                }
//E
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
                        //B kappa Avg Rmin
                        hist->histNNSubXi2pcfthreadp[n] =
                        hist->histNNSubXi2pcfthreadp[n] + 1.;
                        //E
                        hist->histNNSubthread[n] = hist->histNNSubthread[n] + 1.;
                        xi = Kappa(q);
                        if (cmd->computeTPCF) {
                        REAL cosphi,sinphi;
#if NDIM == 3
                        REAL s, sy;
#ifdef SINGLEP
                        float pr0[NDIM];
#else
                        vector pr0;
#endif
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
                        }
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
                        xi = Kappa(q);
                        if (cmd->computeTPCF) {
                            real cosphi,sinphi;
#if NDIM == 3
                            real s, sy;
                            vector pr0;
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
                        }
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
                               kdnode ntab, bodyptr *bptr,
                               INTEGER *nbbcalcthread, INTEGER *nbccalcthread,
                               gdhistptr_sincos_omp hist)
{
//    bodyptr q;
#ifdef SINGLEP
    float dr1;
    float dr[NDIM];
#else
    real dr1;
    real drpq2;
    vector dr;
#endif
    int n;
    real xi;

    int npoints;
    npoints = ntab.last - ntab.first + 1;
//#ifdef OPENMPCODE
//  #pragma omp parallel for
//#endif
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
                        hist->histNNSubXi2pcfthread[n] + 1.0;
                        hist->histNNSubthread[n] = hist->histNNSubthread[n] + 1.0;
                        
#ifdef KappaAvgON
                        xi = ntab.kappa;
#else
                        xi = ntab.kappa;
#endif

                        if (cmd->computeTPCF) {
                        REAL cosphi,sinphi;
#if NDIM == 3
                        REAL s, sy;
#ifdef SINGLEP
                        float pr0[NDIM];
#else
                        vector pr0;
#endif
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
                        } // ! computeTPCF
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
                        hist->histNNSubXi2pcfthread[n] + 1.0;
                        hist->histNNSubthread[n] = hist->histNNSubthread[n] + 1.0;
                        xi = ntab.kappa;
                        if (cmd->computeTPCF) {
                        real cosphi, sinphi;
#if NDIM == 3
                        real s, sy;
                        vector pr0;
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
                        } // computeTPCF
                        hist->histXi2pcfthreadsub[n] += xi;
                        *nbccalcthread += 1;
                    } // ! n in (1, sizeHistN)
                } // ! dr1 > rminHist
            } // ! useLogHist
        } // ! dr1 < rangeN
//        } // ! accept_body
//    } // ! loop pj
}


//#endif // ! NMultipoles
