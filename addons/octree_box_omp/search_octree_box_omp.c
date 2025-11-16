/* ==============================================================================
 MODULE: search_octree_box_omp.c		[cTreeBalls]
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
local int print_info(struct cmdline_data* cmd,
                     struct  global_data* gd);

#include <pthread.h>

/*
 Search omp parallel/serial routine using normal tree method:

 To be called using: search=octree-box-omp

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
global int searchcalc_octree_box_omp(struct  cmdline_data* cmd,
                                    struct  global_data* gd,
                                    bodyptr *btable, INTEGER *nbody,
                                    INTEGER ipmin, INTEGER *ipmax,
                                    int cat1, int cat2)
{
    string routineName = "searchcalc_octree_box_omp";
    double cpustart;
    int nn;
    
    cpustart = CPUTIME;
    print_info(cmd, gd);

#ifdef OPENMPCODE
    ThreadCount(cmd, gd, nbody[cat1], cat1);
    pthread_t main_thread_id = pthread_self();
#endif

    search_init_gd_hist_sincos(cmd, gd);

    INTEGER ipfalse;
    ipfalse=0;
//B kappa Avg Rmin
    INTEGER icountNbRmin;
    icountNbRmin=0;
    INTEGER icountNbRminOverlap;
    icountNbRminOverlap=0;
//E

    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                           "\n%s: Total allocated %g MByte storage so far.\n",
                           routineName, gd->bytes_tot*INMB);
    verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                        "\nRunning...\n - Completed pivot node:\n");

#ifdef OPENMPCODE
#ifndef BALLS4SCANLEV
#pragma omp parallel default(none)   \
    shared(cmd,gd,btable,nbody,roottable,ipmin,ipmax,nodetablescanlev, \
           nodetablescanlev_root,rootnode, cat1, cat2, ipfalse, icountNbRmin, \
           icountNbRminOverlap, main_thread_id)
#else
#pragma omp parallel default(none)   \
    shared(cmd,gd,btable,nbody,roottable,ipmin,ipmax,nodetablescanlev, \
           nodetablescanlev_root,rootnode, cat1, cat2, ipfalse, icountNbRmin, \
           icountNbRminOverlap,nodetablescanlevB4, main_thread_id)
#endif
    {
        pthread_t current_thread_id = pthread_self();
#endif

#ifndef BALLS4SCANLEV
        bodyptr p;
#else
        nodeptr p;
#endif
        int n;
        int ip;
        INTEGER nbbcalcthread = 0;
        INTEGER nbccalcthread = 0;
        gdhist_sincos_omp hist;

        search_init_sincos_omp(cmd, gd, &hist);

        INTEGER ipfalsethreads;
        ipfalsethreads = 0;

//B kappa Avg Rmin
        INTEGER icountNbRminthread;
        icountNbRminthread=0;
        INTEGER icountNbRminOverlapthread;
        icountNbRminOverlapthread=0;
//E

#ifdef OPENMPCODE
#pragma omp for nowait schedule(dynamic)
#endif
#ifndef BALLS4SCANLEV
        for (p = btable[cat1] + ipmin -1; p < btable[cat1] + ipmax[cat1]; p++) {
#else
        for (INTEGER i=0; i< gd->nnodescanlevTableB4[cat1]; i++) {
            p = nodetablescanlevB4[cat1][i];
#endif
// p and q are in differents node structures... cat1!=cat2...
//B kappa Avg Rmin
            NbRmin(p) = 1;
            NbRminOverlap(p) = 0;
            KappaRmin(p) = Kappa(p);
            if (scanopt(cmd->options, "smooth-pivot")) {
                if (Update(p) == FALSE) {
                    ipfalsethreads++;
                    continue;
                }
            }
//E
            for (n = 1; n <= cmd->sizeHistN; n++) {
                hist.histXi2pcfthreadsub[n] = 0.0;  // Affects only 2pcf
//B kappa Avg Rmin
                hist.histNNSubXi2pcfthreadp[n] = 0.;// Affects only 2pcf
//E
            }

#ifndef BALLS4SCANLEV
            normal_walktree_sincos(cmd, gd, p, ((nodeptr) roottable[cat2]),
                        gd->rSizeTable[cat2], &nbbcalcthread, &nbccalcthread, &hist);
//B kappa Avg Rmin
            for (n = 1; n <= cmd->sizeHistN; n++) {
                hist.histNNSubXi2pcfthreadp[n] =
                            ((real)NbRmin(p))*hist.histNNSubXi2pcfthreadp[n];
                hist.histNNSubXi2pcfthreadtotal[n] +=
                            hist.histNNSubXi2pcfthreadp[n];
            }
//E
            computeBodyProperties_sincos(cmd, gd, p, nbody[cat1], &hist);
#else // ! BALLS4SCANLEV
            normal_walktree_sincos(cmd, gd, (bodyptr)p, ((nodeptr) roottable[cat2]),
                        gd->rSizeTable[cat2], &nbbcalcthread, &nbccalcthread, &hist);
//B kappa Avg Rmin
            for (n = 1; n <= cmd->sizeHistN; n++) {
                hist.histNNSubXi2pcfthreadp[n] =
                            ((real)NbRmin(p))*hist.histNNSubXi2pcfthreadp[n];
                hist.histNNSubXi2pcfthreadtotal[n] +=
                            hist.histNNSubXi2pcfthreadp[n];
            }
//E
            computeBodyProperties_sincos(cmd, gd, (bodyptr)p,
                                         gd->nnodescanlevTableB4[cat1], &hist);
#endif // ! BALLS4SCANLEV

//B kappa Avg Rmin
            icountNbRminthread += NbRmin(p);
            icountNbRminOverlapthread += NbRminOverlap(p);
//E
            INTEGER ip;
#ifndef BALLS4SCANLEV
          ip = p - btable[cat1] + 1;
#else
          ip = i+1;
#endif
//            if (ip%cmd->stepState == 0) {
//                verb_log_print(cmd->verbose_log, gd->outlog,
//                               " - Completed pivot: %ld\n", ip);
//            }
            
            if (ip%gd->stepState == 0)
            verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog, ".");
            if (ip%cmd->stepState == 0)
            verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                                "%d ", ip);
#ifdef OPENMPCODE
        } // end do body p // end pragma omp DO_BODY p

        if (main_thread_id == current_thread_id)
                verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                                "\n", ip);
#else
      verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                                "\n", ip);
#endif

#ifdef OPENMPCODE
#pragma omp critical
        {
#endif
            for (n = 1; n <= cmd->sizeHistN; n++) {
                gd->histNN[n] += hist.histNthread[n];
                gd->histNNSubXi2pcf[n] += hist.histNNSubXi2pcfthread[n];
//B kappa Avg Rmin
                gd->histNNSubXi2pcftotal[n] += hist.histNNSubXi2pcfthreadtotal[n];
//E
                //B Check this line and the histogram array histXi2pcfthread
                //  and correct if necessary
                gd->histXi2pcf[n] += hist.histXi2pcfthread[n];
                //E
            }
            gd->nbbcalc += nbbcalcthread;
            gd->nbccalc += nbccalcthread;
            //B kappa Avg Rmin
            ipfalse += ipfalsethreads;
            icountNbRmin += icountNbRminthread;
            icountNbRminOverlap += icountNbRminOverlapthread;
            //E
#ifdef OPENMPCODE
        }
#endif
        search_free_sincos_omp(cmd, gd, &hist);

#ifdef OPENMPCODE
    } // end pragma omp parallel
#endif

    real xi, den, num;
    int mm;
    if (scanopt(cmd->options, "smooth-pivot")) {
#ifdef BALLS4SCANLEV
        
        num = (real)gd->nnodescanlevTableB4[cat1];
        den = (real)(gd->nnodescanlevTableB4[cat1]-ipfalse);

#ifdef NOSTANDARNORMHIST
        xi = 1.0;
#else
        xi = num/den;
#endif // ! NONORMHIST

#else // ! BALLS4SCANLEV
        num = (real)nbody[cat1];
//B kappa Avg Rmin
        den = (real)(nbody[cat1]-ipfalse);
//E
#ifdef NOSTANDARNORMHIST
        xi = 1.0;
#else
        xi = num/den;
#endif // ! NONORMHIST
#endif // ! BALLS4SCANLEV

        verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                               "%s: p falses found = %ld and %e %e %e\n",
                               routineName, ipfalse, num, den, xi);
    }
    
    if (!scanopt(cmd->options, "asymmetric")) {
        for (nn = 1; nn <= cmd->sizeHistN; nn++) {
            if (cmd->verbose>3)
                printf("%d %e %e\n", nn,
                   gd->histNNSubXi2pcf[nn], gd->histNNSubXi2pcftotal[nn]);
            gd->histXi2pcf[nn] /= 2.0;
            gd->histNNSubXi2pcf[nn] /= 2.0;
//B kappa Avg Rmin
            gd->histNNSubXi2pcftotal[nn] /= 2.0;
            if (scanopt(cmd->options, "smooth-pivot")) {
                gd->histXi2pcf[nn] /= MAX(gd->histNNSubXi2pcftotal[nn],1.0);
            } else {
                gd->histXi2pcf[nn] /= MAX(gd->histNNSubXi2pcf[nn],1.0);
            }
//E
        }
    } else {
        for (nn = 1; nn <= cmd->sizeHistN; nn++) {
            if (cmd->verbose>3)
            printf(0,"%d %e %e\n", nn,
                   gd->histNNSubXi2pcf[nn], gd->histNNSubXi2pcftotal[nn]);
            if (scanopt(cmd->options, "smooth-pivot")) {
                gd->histXi2pcf[nn] /= MAX(gd->histNNSubXi2pcftotal[nn],1.0);
            } else {
                gd->histXi2pcf[nn] /= MAX(gd->histNNSubXi2pcf[nn],1.0);
            }
        }
    }

    if (scanopt(cmd->options, "compute-HistN")) {
        if (scanopt(cmd->options, "smooth-pivot")) {
            search_compute_HistN(cmd, gd, nbody[cat1]-ipfalse);
        } else {
            search_compute_HistN(cmd, gd, nbody[cat1]);
        }
    }

    //B kappa Avg Rmin
    if (scanopt(cmd->options, "smooth-pivot")) {
        verb_print(cmd->verbose, "%s: p falses found = %ld\n",routineName,ipfalse);
        verb_print(cmd->verbose,
                   "%s: count NbRmin found = %ld\n", routineName, icountNbRmin);
        verb_print(cmd->verbose,
                   "%s: count overlap found = %ld\n",routineName, icountNbRminOverlap);
        
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
        verb_print(cmd->verbose, "%s: p falses found = %ld\n",routineName, ifalsecount);
        verb_print(cmd->verbose, "%s: p true found = %ld\n",routineName, itruecount);
        verb_print(cmd->verbose, "%s: total = %ld\n",
                   routineName, itruecount+ifalsecount);
    }
    //E

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

    if (!Update(p)) return;
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
                                                    nbbcalcthread, nbccalcthread,
                                                    hist);
                            else
                                for (l = More(q); l != Next(q); l = Next(l))
                                    normal_walktree_sincos(cmd, gd, p,l,qsize/2,
                                                            nbbcalcthread,
                                                            nbccalcthread, hist);
                        } else {
                            //B Segment as original
                            if (cmd->rminHist==0)
                                n = (int)(cmd->logHistBinsPD*(rlog10(dr1)
                                    - rlog10(cmd->rangeN)) + cmd->sizeHistN);
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
                                    for (l = More(q); l != Next(q); l = Next(l))
                                        normal_walktree_sincos(cmd, gd, p,l,qsize/2,
                                            nbbcalcthread, nbccalcthread, hist);
                                    }
                            } else
                                for (l = More(q); l != Next(q); l = Next(l))
                                    normal_walktree_sincos(cmd, gd, p,l,qsize/2,
                                            nbbcalcthread, nbccalcthread, hist);
                            //E Segment as original
                        }
                    } else {  // ! useLogHist
                        if (scanopt(cmd->options, "behavior-ball")) {
                            if ( (Radius(p)+Radius(q)) < gd->deltaR*THETA)
                                sumnode_sincos_cell(cmd, gd, p, ((cellptr) q),
                                                    ( (cellptr) q+1),
                                                    nbbcalcthread, nbccalcthread,
                                                    hist);
                            else
                                for (l = More(q); l != Next(q); l = Next(l))
                                    normal_walktree_sincos(cmd, gd, p,l,qsize/2,
                                                            nbbcalcthread,
                                                            nbccalcthread, hist);
                        } else {
                            real rBin;
                            n = (int) ((dr1-cmd->rminHist) * gd->i_deltaR) + 1;
                            rBin = cmd->rminHist + ((real)n)*gd->deltaR;
                            if ( rBin-gd->deltaR < dr1 - Radius(q)
                                && dr1 + Radius(q) < rBin ) {
                                sumnode_sincos_cell(cmd, gd, p, ((cellptr) q),
                                                    ( (cellptr) q+1), nbbcalcthread,
                                                    nbccalcthread, hist);
                            } else {
                                for (l = More(q); l != Next(q); l = Next(l))
                                    normal_walktree_sincos(cmd, gd, p,l,qsize/2,
                                                           nbbcalcthread,
                                                           nbccalcthread, hist);
                            }
                        }
                    } // ! useLogHist
                } else {    // ! no-one-ball
                    for (l = More(q); l != Next(q); l = Next(l))
                        normal_walktree_sincos(cmd, gd, p,l,qsize/2,
                                nbbcalcthread, nbccalcthread, hist);
                }           // ! no-one-ball
            } // ! reject_cell
        } else { // ! Type(q) == CELL
            if (Type(q) == BODY)
                sumnode_sincos(cmd, gd, p, ((cellptr) q),
                        ( (cellptr) q+1), nbbcalcthread, nbccalcthread, hist);
        } // ! Type(q) == CELL
    } // ! p != q
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
            //B kappa Avg Rmin
            if (scanopt(cmd->options, "smooth-pivot"))
                if (dr1<=gd->rsmooth[0]) {
                    if (Update(q)==TRUE) {
                        Update(q) = FALSE;
                        NbRmin(p) += 1;
                        KappaRmin(p) += Weight(q)*Kappa(q);
                    } else {
                        NbRminOverlap(p) += 1;
                    }
                }
            //E
            if (cmd->useLogHist) {
                if(dr1>cmd->rminHist) {
                    if (cmd->rminHist==0)
                        n = (int)(cmd->logHistBinsPD*(rlog10(dr1) - rlog10(cmd->rangeN))
                                  + cmd->sizeHistN) + 1;
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
                        xi = Weight(q)*Kappa(q);
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
                        //B kappa Avg Rmin
                        hist->histNNSubXi2pcfthreadp[n] =
                        hist->histNNSubXi2pcfthreadp[n] + 1.;
                        //E
                        xi = Weight(q)*Kappa(q);
                        hist->histXi2pcfthreadsub[n] += xi;
                        *nbbcalcthread += 1;
                    }
                }
            } // ! useLogHist
        } // ! accept_body
    }
}

local void sumnode_sincos_cell(struct  cmdline_data* cmd, struct  global_data* gd,
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
                        n = (int)(cmd->logHistBinsPD*(rlog10(dr1) - rlog10(cmd->rangeN))
                                  + cmd->sizeHistN) + 1;
                    else
                        n = (int)(rlog10(dr1/cmd->rminHist) * gd->i_deltaR) + 1;
                    if (n<=cmd->sizeHistN && n>=1) {
                        hist->histNthread[n] = hist->histNthread[n] +  Nb(q);
                        hist->histNNSubXi2pcfthread[n] =
                        hist->histNNSubXi2pcfthread[n] + 1.0;
                        //B kappa Avg Rmin
                        hist->histNNSubXi2pcfthreadp[n] =
                        hist->histNNSubXi2pcfthreadp[n] + 1.;
                        //E
                        xi = Kappa(q);
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
                        //B kappa Avg Rmin
                        hist->histNNSubXi2pcfthreadp[n] =
                        hist->histNNSubXi2pcfthreadp[n] + 1.;
                        //E
                        xi = Kappa(q);
                        hist->histXi2pcfthreadsub[n] += xi;
                        *nbccalcthread += 1;
                    }
                }
            } // ! useLogHist
        } // ! accept_body
    }
}

local int print_info(struct cmdline_data* cmd,
                                  struct  global_data* gd)
{
    verb_print(cmd->verbose, "\nSearch: Running ... (tree-box-omp) \n");

    if (scanopt(cmd->options, "behavior-ball"))
        verb_print(cmd->verbose, "with option behavior-ball... \n");
    if (scanopt(cmd->options, "no-one-ball"))
        verb_print(cmd->verbose, "with option no-one-ball... \n");
    if (scanopt(cmd->options, "smooth-pivot"))
        verb_print(cmd->verbose,
                   "with option smooth-pivot... rsmooth=%g\n",gd->rsmooth[0]);
    if (!cmd->useLogHist)
        verb_print(cmd->verbose, "using normal histogram scale... \n");
    if (!cmd->computeTPCF)
        verb_print(cmd->verbose, "computing only 2pcf... \n");
    else
        verb_print(cmd->verbose,
                   "Warning: activated computing 3pcf, but won´t do it... \n");
    if (!cmd->usePeriodic)
        verb_print(cmd->verbose, "using non-periodic boundary conditions... \n");
    if (scanopt(cmd->options, "kappa-constant-one"))
        verb_print(cmd->verbose, "kappa constant = 1... \n");
#ifdef NOSTANDARNORMHIST
    verb_print(cmd->verbose, "warning!! histograms will not be normalized... \n");
#endif
#ifdef BALLS4SCANLEV
    verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                        "with BALLS4SCANLEV... \n");
#endif

    return SUCCESS;
}
