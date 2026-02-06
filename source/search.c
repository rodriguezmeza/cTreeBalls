/* ==============================================================================
 MODULE: search.c		[cTreeBalls]
 Written by: M.A. Rodriguez-Meza.
 Starting date:    april 2023
 Purpose: 2/3-point correlation functions computation
 Language: C
 Use: searchcalc_normal_omp_sincos();
 Major revisions:
 ==============================================================================*/
//        1          2          3          4        ^ 5          6          7

// TODO :: Work to do in order to use with boxes not centered at (0,0,...)

//
// lines where there is a "//B socket:" string are places to include module files
//  that can be found in addons/addons_include folder
//

#include "globaldefs.h"

local void normal_walktree_sincos(struct  cmdline_data* cmd, 
                                  struct  global_data* gd,
                                  bodyptr, nodeptr, real,
                                  INTEGER *, INTEGER *, gdhistptr_sincos_omp);
local void sumnode_sincos(struct  cmdline_data* cmd, 
                          struct  global_data* gd,
                          bodyptr, cellptr, cellptr,
                          INTEGER *, INTEGER *, gdhistptr_sincos_omp);
local void sumnode_sincos_cell(struct  cmdline_data* cmd,
                               struct  global_data* gd,
                               bodyptr, cellptr, cellptr,
                               INTEGER *, INTEGER *, gdhistptr_sincos_omp);

/*
 Search omp parallel/serial routine using normal tree method:

 To be called using: search=octree-sincos-omp (substitute tree-omp-sincos)

 Arguments:
    * `cmd`: Input: structure cmdline_data pointer
    * `gd`: Input: structure global_data pointer
    * `nbody`: Input: number of points in table array
    * `ipmin`: Input: minimum point in table array to analyse
    * `ipmax`: Input: maximum point in table array to analyse
    * `cat1`: Input: catalog tag to act as pivot catalog
    * `cat2`: Input: catalog tag to act as a scanning catalog
    * Global tructures used: gd, cmd
    * Histograms outputs (in global gd): histZetaMcos, histZetaMsin,
    *                                    histZetaMsincos, histNN,
    *                                    histNNSubXi2pcf, histNNSubXi2pcftotal,
    *                                    histXi2pcf, histXi,
    * Counting encounters (in global gd): nbbcalc, nbccalc, ncccalc
 Return (the error status):
    int SUCCESS or FAILURE
 */
global int searchcalc_normal_sincos(struct  cmdline_data* cmd,
                                 struct  global_data* gd,
                                 bodyptr *btable, INTEGER *nbody,
                                 INTEGER ipmin, INTEGER *ipmax,
                                 int cat1, int cat2)
{
    double cpustart;
    int nn;
    
    cpustart = CPUTIME;
    verb_print(cmd->verbose,
               "\nEvalHistograms: Running ... (octree-sincos-omp) \n");
    
    if (scanopt(cmd->options, "behavior-ball"))
        verb_print(cmd->verbose, "with option behavior-ball... \n");
    if (scanopt(cmd->options, "no-one-ball"))
        verb_print(cmd->verbose, "with option no-one-ball... \n");
#ifdef SMOOTHPIVOT
        verb_print(cmd->verbose,
                   "with option smooth-pivot... rsmooth=%g\n",gd->rsmooth[0]);
#endif
#ifndef TPCF
        verb_print(cmd->verbose, "computing only 2pcf... \n");
#endif
    if (scanopt(cmd->options, "kappa-constant-one"))
        verb_print(cmd->verbose, "kappa constant = 1... \n");
#ifdef NOSTANDARNORMHIST
    verb_print(cmd->verbose, "warning!! histograms will not be normalized... \n");
#endif

#ifdef OPENMPCODE
    ThreadCount(cmd, gd, nbody[cat1], cat1);
#endif

    search_init_gd_hist_sincos(cmd, gd);

    INTEGER ipfalse;
    ipfalse=0;
#ifdef SMOOTHPIVOT
    INTEGER icountNbRmin;
    icountNbRmin=0;
    INTEGER icountNbRminOverlap;
    icountNbRminOverlap=0;
#endif

    verb_print(cmd->verbose,
               "\nsearchcalc_balls: Total allocated %g MByte storage so far.\n",
               gd->bytes_tot/(1024.0*1024.0));

#ifdef SMOOTHPIVOT
#pragma omp parallel default(none)   \
    shared(cmd,gd,btable,nbody,roottable,ipmin,ipmax,nodetablescanlev, \
           nodetablescanlev_root,rootnode, cat1, cat2, ipfalse, icountNbRmin, \
           icountNbRminOverlap)
    {
#else
#pragma omp parallel default(none)   \
    shared(cmd,gd,btable,nbody,roottable,ipmin,ipmax,nodetablescanlev, \
           nodetablescanlev_root,rootnode, cat1, cat2, ipfalse)
    {
#endif
        bodyptr p;
        int n;
        INTEGER nbbcalcthread = 0;
        INTEGER nbccalcthread = 0;
        gdhist_sincos_omp hist;

        search_init_sincos_omp(cmd, gd, &hist);

        INTEGER ipfalsethreads;
        ipfalsethreads = 0;

#ifdef SMOOTHPIVOT
        INTEGER icountNbRminthread;
        icountNbRminthread=0;
        INTEGER icountNbRminOverlapthread;
        icountNbRminOverlapthread=0;
#endif
        
#pragma omp for nowait schedule(dynamic)
        for (p = btable[cat1] + ipmin -1; p < btable[cat1] + ipmax[cat1]; p++) {
// p and q are in differents node structures... cat1!=cat2...
#ifdef SMOOTHPIVOT
            NbRmin(p) = 1;
            NbRminOverlap(p) = 0;
            KappaRmin(p) = Kappa(p);
                if (Update(p) == FALSE) {
                    ipfalsethreads++;
                    continue;
                }
#endif
            for (n = 1; n <= cmd->sizeHistN; n++) {
                hist.histNNSubthread[n] = 0.0;      // Affects only 3pcf
                hist.histXi2pcfthreadsub[n] = 0.0;  // Affects only 2pcf
#ifdef SMOOTHPIVOT
                hist.histNNSubXi2pcfthreadp[n] = 0.;// Affects only 2pcf
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
#ifdef SINGLEP
                hist.drpq = sqrt(hist.drpq2);
#else
                hist.drpq = rsqrt(hist.drpq2);
#endif
#endif // ! NDIM
#endif

            normal_walktree_sincos(cmd, gd, p, ((nodeptr) roottable[cat2]),
                    gd->rSizeTable[cat2], &nbbcalcthread, &nbccalcthread, &hist);
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
            computeBodyProperties_sincos(cmd, gd, p, nbody[cat1], &hist);

#ifdef SMOOTHPIVOT
            icountNbRminthread += NbRmin(p);
            icountNbRminOverlapthread += NbRminOverlap(p);
#endif
            INTEGER ip;
            ip = p - btable[cat1] + 1;
            if (ip%cmd->stepState == 0) {
                verb_log_print(cmd->verbose_log, gd->outlog, 
                               " - Completed pivot: %ld\n", ip);
            }
        } // end do body p // end pragma omp DO_BODY p

#pragma omp critical
        {
            for (int n = 1; n <= cmd->sizeHistN; n++) {
                gd->histNN[n] += hist.histNthread[n];
                gd->histNNSub[n] += hist.histNNSubthread[n];
                gd->histNNSubXi2pcf[n] += hist.histNNSubXi2pcfthread[n];
#ifdef SMOOTHPIVOT
                gd->histNNSubXi2pcftotal[n] += hist.histNNSubXi2pcfthreadtotal[n];
#endif
                //B Check this line and the histogram array histXi2pcfthread
                //  and correct if necessary
                gd->histXi2pcf[n] += hist.histXi2pcfthread[n];
                //E
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
        search_free_sincos_omp(cmd, gd, &hist);
    } // end pragma omp parallel

    real xi, den, num;
    int mm;
#ifdef SMOOTHPIVOT
        num = (real)nbody[cat1];
        den = (real)(nbody[cat1]-ipfalse);
#ifdef NOSTANDARNORMHIST
        xi = 1.0;
#else
        xi = num/den;
#endif // ! NONORMHIST
        verb_print(cmd->verbose,
                   "tree-omp-sincos: p falses found = %ld and %e %e %e\n",
                   ipfalse, num, den, xi);
#ifdef TPCF
            for (mm=1; mm<=cmd->mChebyshev+1; mm++) {
                MULMS_ext(gd->histZetaMcos[mm],
                          gd->histZetaMcos[mm],xi,cmd->sizeHistN);
                MULMS_ext(gd->histZetaMsin[mm],
                          gd->histZetaMsin[mm],xi,cmd->sizeHistN);
                MULMS_ext(gd->histZetaMsincos[mm],
                          gd->histZetaMsincos[mm],xi,cmd->sizeHistN);
                // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
                MULMS_ext(gd->histZetaMcossin[mm],
                          gd->histZetaMcossin[mm],xi,cmd->sizeHistN);
            }
#endif
#endif

    if (!scanopt(cmd->options, "asymmetric")) {
        for (nn = 1; nn <= cmd->sizeHistN; nn++) {
#ifdef SMOOTHPIVOT
            if (cmd->verbose>3)
                printf("%d %e %e\n", nn,
                   gd->histNNSubXi2pcf[nn], gd->histNNSubXi2pcftotal[nn]);
#else
            if (cmd->verbose>3)
                printf("%d %e\n", nn,
                   gd->histNNSubXi2pcf[nn]);
#endif
            gd->histXi2pcf[nn] /= 2.0;
            gd->histNNSubXi2pcf[nn] /= 2.0;
#ifdef SMOOTHPIVOT
            gd->histNNSubXi2pcftotal[nn] /= 2.0;
                gd->histXi2pcf[nn] /= MAX(gd->histNNSubXi2pcftotal[nn],1.0);
#else
                gd->histXi2pcf[nn] /= MAX(gd->histNNSubXi2pcf[nn],1.0);
#endif
        }
    } else {
        for (nn = 1; nn <= cmd->sizeHistN; nn++) {
#ifdef SMOOTHPIVOT
            if (cmd->verbose>3)
            printf(0,"%d %e %e\n", nn,
                   gd->histNNSubXi2pcf[nn], gd->histNNSubXi2pcftotal[nn]);
                gd->histXi2pcf[nn] /= MAX(gd->histNNSubXi2pcftotal[nn],1.0);
#else
            if (cmd->verbose>3)
            printf(0,"%d %e\n", nn,
                   gd->histNNSubXi2pcf[nn]);
                gd->histXi2pcf[nn] /= MAX(gd->histNNSubXi2pcf[nn],1.0);
#endif
        }
    }

    if (scanopt(cmd->options, "compute-HistN")) {
#ifdef SMOOTHPIVOT
            search_compute_HistN(cmd, gd, nbody[cat1]-ipfalse);
#else
            search_compute_HistN(cmd, gd, nbody[cat1]);
#endif
    }

#ifdef SMOOTHPIVOT
        verb_print(cmd->verbose, "tree-omp-sincos: p falses found = %ld\n",
                   ipfalse);
        verb_print(cmd->verbose,
                   "tree-omp-sincos: count NbRmin found = %ld\n",icountNbRmin);
        verb_print(cmd->verbose,
                   "tree-omp-sincos: count overlap found = %ld\n",
                   icountNbRminOverlap);
        
        bodyptr pp;
        INTEGER ifalsecount;
        ifalsecount = 0;
        INTEGER itruecount;
        itruecount = 0;
        for (pp = btable[cat1] + ipmin -1; pp < btable[cat1] + ipmax[cat1]; pp++) {
            if (Update(pp) == FALSE) {
                ifalsecount++;
            } else {
                itruecount++;
            }
        }
        verb_print(cmd->verbose, "tree-omp-sincos: p falses found = %ld\n",
                   ifalsecount);
        verb_print(cmd->verbose, "tree-omp-sincos: p true found = %ld\n",
                   itruecount);
        verb_print(cmd->verbose, "tree-omp-sincos: total = %ld\n",
                   itruecount+ifalsecount);
#endif

    gd->cpusearch = CPUTIME - cpustart;
    verb_print(cmd->verbose, "Going out: CPU time = %lf\n",CPUTIME-cpustart);

    return SUCCESS;
}

local void normal_walktree_sincos(struct  cmdline_data* cmd,
                                  struct  global_data* gd,
                                  bodyptr p, nodeptr q, real qsize,
                                  INTEGER *nbbcalcthread, INTEGER *nbccalcthread,
                                  gdhistptr_sincos_omp hist)
{
    nodeptr l;
#ifdef SINGLEP
    float dr1;
    float dr[NDIM];
#else
    real dr1;
    vector dr;
#endif
    int n;

    if (Update(p)) {
        if ( ((nodeptr) p) != q ) {
            if (Type(q) == CELL) {
                if (!reject_cell(cmd, gd, (nodeptr)p, q, qsize)) {
                    if (!scanopt(cmd->options, "no-one-ball")) {
                        accept_body(cmd, gd, p, (nodeptr)q, &dr1, dr);
                        if (cmd->useLogHist) {
                            if (scanopt(cmd->options, "behavior-ball")) {
                                if ( (Radius(p)+Radius(q))/(dr1) < gd->deltaR)
                                    sumnode_sincos_cell(cmd, gd, p, ((cellptr) q),
                                            ( (cellptr) q+1),
                                            nbbcalcthread, nbccalcthread, hist);
                                else
                                    for (l = More(q); l != Next(q); l = Next(l))
                                        normal_walktree_sincos(cmd, gd,
                                                               p,l,qsize/2,
                                            nbbcalcthread, nbccalcthread, hist);
                            } else {
                                //B Segment as original
                                if (cmd->rminHist==0)
                                    n = (int)(cmd->logHistBinsPD*(rlog10(dr1)
                                            - rlog10(cmd->rangeN))
                                            + cmd->sizeHistN);
                                else
                                    n = (int)(rlog10(dr1/cmd->rminHist) 
                                              * gd->i_deltaR);
                                if (n<=cmd->sizeHistN-1 && n>=1) {
                                    if ( gd->deltaRV[n] < dr1 - Radius(q)
                                        && dr1 + Radius(q) < gd->deltaRV[n+1]) {
                                        sumnode_sincos_cell(cmd, gd, p, 
                                            ((cellptr) q), ( (cellptr) q+1),
                                            nbbcalcthread, nbccalcthread, hist);
                                    } else {
                                        for (l = More(q); l != Next(q); 
                                             l = Next(l))
                                            normal_walktree_sincos(cmd, gd, 
                                                p,l,qsize/2,
                                            nbbcalcthread, nbccalcthread, hist);
                                    }
                                } else
                                    for (l = More(q); l != Next(q); l = Next(l))
                                        normal_walktree_sincos(cmd, gd, 
                                            p,l,qsize/2,
                                            nbbcalcthread, nbccalcthread, hist);
                                //E Segment as original
                            }
                        } else {  // ! useLogHist
                            real rBin;
                            n = (int) ((dr1-cmd->rminHist) * gd->i_deltaR) + 1;
                            rBin = cmd->rminHist + ((real)n)*gd->deltaR;
                            if ( rBin-gd->deltaR < dr1 - Radius(q)
                                && dr1 + Radius(q) < rBin ) {
                                sumnode_sincos_cell(cmd, gd, p, ((cellptr) q),
                                                    ( (cellptr) q+1), 
                                                    nbbcalcthread, nbccalcthread,
                                                    hist);
                            } else {
                                for (l = More(q); l != Next(q); l = Next(l))
                                    normal_walktree_sincos(cmd, gd, p,l,qsize/2,
                                            nbbcalcthread, nbccalcthread, hist);
                            }
                        } // ! useLogHist
                    } else {    // ! no-one-ball
                        for (l = More(q); l != Next(q); l = Next(l))
                            normal_walktree_sincos(cmd, gd, p,l,qsize/2,
                                    nbbcalcthread, nbccalcthread, hist);
                    }           // ! no-one-ball
                } // ! reject_cell
            } else { // ! Type(q) == CELL
                sumnode_sincos(cmd, gd, p, ((cellptr) q),
                            ( (cellptr) q+1), nbbcalcthread, nbccalcthread, hist);
            } // ! Type(q) == CELL
        } // ! p != q
    } // ! Update
}

local void sumnode_sincos(struct  cmdline_data* cmd, struct  global_data* gd,
                          bodyptr p, cellptr start, cellptr finish,
                          INTEGER *nbbcalcthread, INTEGER *nbccalcthread,
                          gdhistptr_sincos_omp hist)
{
    cellptr q;
#ifdef SINGLEP
    float dr1;
    float dr[NDIM];
#else
    real dr1;
    vector dr;
#endif
    int n;
    real xi;

    for (q = start; q < finish; q++) {
        if (accept_body(cmd, gd, p, (nodeptr)q, &dr1, dr)) {
#ifdef SMOOTHPIVOT
                if (dr1<=gd->rsmooth[0]) {
                    if (Update(q)==TRUE) {
                        Update(q) = FALSE;
                        NbRmin(p) += 1;
                        KappaRmin(p) += Kappa(q);
                    } else {
                        NbRminOverlap(p) += 1;
                    }
                }
#endif
            if (cmd->useLogHist) {
                if(dr1>cmd->rminHist) {
                    if (cmd->rminHist==0)
                        n = (int)(cmd->logHistBinsPD*(rlog10(dr1)
                                - rlog10(cmd->rangeN))
                                + cmd->sizeHistN) + 1;
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
                        xi = Kappa(q);
#ifdef TPCF
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
                        xi = Kappa(q);
#ifdef TPCF
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
#endif
                        hist->histXi2pcfthreadsub[n] += xi;
                        *nbbcalcthread += 1;
                    }
                }
            } // ! useLogHist
        } // ! accept_body
    }
}

local void sumnode_sincos_cell(struct  cmdline_data* cmd,
                               struct  global_data* gd,
                               bodyptr p, cellptr start, cellptr finish,
                               INTEGER *nbbcalcthread, INTEGER *nbccalcthread,
                               gdhistptr_sincos_omp hist)
{
    cellptr q;
#ifdef SINGLEP
    float dr1;
    float dr[NDIM];
#else
    real dr1;
    vector dr;
#endif
    int n;
    real xi;

    for (q = start; q < finish; q++) {
        if (accept_body(cmd, gd, p, (nodeptr)q, &dr1, dr)) {
            if (cmd->useLogHist) {
                if(dr1>cmd->rminHist) {
                    if (cmd->rminHist==0)
                        n = (int)(cmd->logHistBinsPD*(rlog10(dr1)
                                - rlog10(cmd->rangeN))
                                + cmd->sizeHistN) + 1;
                    else
                        n = (int)(rlog10(dr1/cmd->rminHist) * gd->i_deltaR) + 1;
                    if (n<=cmd->sizeHistN && n>=1) {
                        hist->histNthread[n] = hist->histNthread[n] +  Nb(q);
                        hist->histNNSubXi2pcfthread[n] =
                        hist->histNNSubXi2pcfthread[n] + 1.0;
                        hist->histNNSubthread[n] = hist->histNNSubthread[n] + 1.0;
                        
#ifdef KappaAvgON
                        xi = KappaAvg(q)/Nb(q);
#else
                        xi = Kappa(q);
#endif

#ifdef TPCF
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
#endif
                        hist->histXi2pcfthreadsub[n] += xi;
                        *nbccalcthread += 1;
                    }
                }
            } else {  // ! useLogHist
                if(dr1>cmd->rminHist) {
                    n = (int) ( (dr1-cmd->rminHist) * gd->i_deltaR) + 1;
                    if (n<=cmd->sizeHistN && n>=1) {
                        hist->histNthread[n] = hist->histNthread[n] +  Nb(q);
                        hist->histNNSubXi2pcfthread[n] =
                        hist->histNNSubXi2pcfthread[n] + 1.0;
                        hist->histNNSubthread[n] = hist->histNNSubthread[n] + 1.0;
                        xi = Kappa(q);
#ifdef TPCF
                        real cosphi, sinphi;
#if NDIM == 3
                        real s, sy, phi;
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
#endif
                        hist->histXi2pcfthreadsub[n] += xi;
                        *nbccalcthread += 1;
                    }
                }
            } // ! useLogHist
        } // ! accept_body
    }
}

