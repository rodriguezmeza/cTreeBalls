/* ==============================================================================
 MODULE: search_octree_kk_balls4_omp.c		[cTreeBalls]
 Written by: M.A. Rodriguez-Meza
 Starting date:    april 2023
 Purpose: 2-point correlation function computation
 Language: C
 Use: searchcalc_octree_kk_balls4_omp(cmd, gd, btable, nbody,
                                           ipmin, ipmax, cat1, cat2);
 Major revisions:
 ==============================================================================*/
//        1          2          3          4        ^ 5          6          7


// Work to do in order to use with boxes not centered at (0,0,...)

#include "globaldefs.h"


typedef struct {
    INTEGER actlen;
    nodeptr *active;
    cellptr interact;
} gdhist_omp_balls6_kk, *gdhistptr_omp_balls6_kk;

typedef struct {
    realptr histNNSubthread;

    INTEGER nbbcalcthread;
    INTEGER nbccalcthread;
    INTEGER ncccalcthread;
    INTEGER ibfcountthread;
    INTEGER nsmoothcountthread;
    //B smooth-pivot
    INTEGER ipfalsethreads;
    INTEGER icountNbRminthread;
    INTEGER icountNbRminOverlapthread;
    //E

//    vector q0;
//    real drpq2, drpq;
//    vector dr0;
//    real cosb;
//    real sinb;

    INTEGER ipcount;
} gdhist_sincos_omp_kk, *gdhistptr_sincos_omp_kk;

#ifndef BALLS
local bool nodes_condition_balls(struct cmdline_data* cmd, struct  global_data* gd,
                                 nodeptr p, nodeptr q, real *dr1, vector dr);
#endif

local void walktree_balls6_omp(struct cmdline_data* cmd, struct  global_data* gd,
                               nodeptr *, nodeptr *, cellptr, cellptr,
                               nodeptr, real, vector,
                               gdhistptr_omp_balls6_kk,
                               gdhistptr_sincos_omp_kk, INTEGER, int);
local void walksub6_omp(struct cmdline_data* cmd, struct  global_data* gd,
                        nodeptr *, nodeptr *, cellptr, cellptr,
                        nodeptr, real, vector, gdhistptr_omp_balls6_kk,
                        gdhistptr_sincos_omp_kk, INTEGER, int);
local int sum_balls6_omp(struct cmdline_data* cmd, struct  global_data* gd,
                          cellptr, cellptr, bodyptr, gdhistptr_omp_balls6_kk,
                          gdhistptr_sincos_omp_kk, INTEGER, int);
local void sumnode_balls6_omp(struct cmdline_data* cmd, struct  global_data* gd,
                              cellptr, cellptr, bodyptr,
                              gdhistptr_sincos_omp_kk);
local void sumcell_balls6_omp(struct cmdline_data* cmd, struct  global_data* gd,
                              cellptr, cellptr, bodyptr,
                              gdhistptr_sincos_omp_kk);
local void sumcellcell_balls6_omp(struct cmdline_data* cmd,
                                  struct  global_data* gd,
                                  cellptr, cellptr, nodeptr,
                                  gdhistptr_sincos_omp_kk);

local int search_init_gd_sincos_omp_kk(struct  cmdline_data* cmd,
                                        struct  global_data* gd);
local int search_init_sincos_omp_kk(struct  cmdline_data* cmd,
                                  struct  global_data* gd,
                                     gdhistptr_sincos_omp_kk hist, int);
local int search_free_sincos_omp_kk(struct  cmdline_data* cmd,
                                         struct  global_data* gd,
                                     gdhistptr_sincos_omp_kk hist);
local int computeBodyProperties_sincos_kk(struct  cmdline_data* cmd,
                                            struct  global_data* gd,
                                            bodyptr p, int nbody,
                                       gdhistptr_sincos_omp_kk hist);

local int computeBodyProperties_sincos_kk_sum_balls6_omp(
                                            struct  cmdline_data*,
                                            struct  global_data*,
                                            bodyptr, INTEGER,
                                            gdhistptr_sincos_omp_kk,
                                            gdhistptr_sincos_omp_kk);

local int search_init_omp_balls6_kk(struct cmdline_data* cmd,
                                         struct  global_data* gd,
                                        gdhistptr_omp_balls6_kk hist, int);
local int search_free_omp_balls6_kk(struct cmdline_data* cmd,
                                  struct  global_data* gd,
                                 gdhistptr_omp_balls6_kk hist);

local int print_info(struct cmdline_data* cmd,
                     struct  global_data* gd);

/*
 Search routine using tree brute force direct method:

 To be called using: search=octree-kk-balls4-omp

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
global int searchcalc_octree_kk_balls4_omp(struct cmdline_data* cmd,
                                             struct  global_data* gd,
                                             bodyptr *btable, INTEGER *nbody,
                                             INTEGER ipmin, INTEGER *ipmax, 
                                             int cat1, int cat2)
{
    double cpustart;
#ifdef DEBUG
    INTEGER ibfcount=0;
#endif

    cpustart = CPUTIME;
    
    print_info(cmd, gd);

#ifdef OPENMPCODE
    ThreadCount(cmd, gd, nbody[cat1], cat1);
#endif

    search_init_gd_sincos_omp_kk(cmd, gd);

#ifdef DEBUG
    bodytabbf = (bodyptr) allocate(nbody[cat1] * sizeof(body));
    gd->bytes_tot += nbody[cat1]*sizeof(body);
    fprintf(stdout,"\n\nAllocated %g MByte for particle (found) storage.\n",
            nbody[cat1]*sizeof(body)/(1024.0*1024.0));

    verb_print(cmd->verbose,
            "\nsearchcalc_balls_omp: Total allocated %g MByte storage so far.\n",
            gd->bytes_tot/(1024.0*1024.0));
#endif

    //B smooth-pivot
    INTEGER ipfalse;
    ipfalse=0;
    INTEGER icountNbRmin;
    icountNbRmin=0;
    INTEGER icountNbRminOverlap;
    icountNbRminOverlap=0;
    //E

    verb_print(cmd->verbose,
        "\nsearchcalc_tc_kk_omp: Total allocated %g MByte storage so far.\n",
        gd->bytes_tot/(1024.0*1024.0));

#ifdef DEBUG
#pragma omp parallel default(none)                                          \
    shared(cmd,gd,btable,nbody,roottable,bodytabbf,                         \
           ibfcount,ipfalse,icountNbRmin,icountNbRminOverlap,               \
           ipmin,ipmax,cat1,cat2,nodetablescanlev,nodetablescanlev_root)
#else
#pragma omp parallel default(none)                                          \
    shared(cmd,gd,btable,nbody,roottable,                                   \
           ipfalse,icountNbRmin,icountNbRminOverlap,                        \
           ipmin,ipmax,cat1,cat2,nodetablescanlev,nodetablescanlev_root)
#endif
  {
    nodeptr p;
      int n;
      int m;
      int i;

    //B init:
    gdhist_omp_balls6_kk histb;
    gdhist_sincos_omp_kk hist;
    search_init_omp_balls6_kk(cmd, gd, &histb, cat1);
    search_init_sincos_omp_kk(cmd, gd, &hist, cat1);
    //E

#pragma omp for nowait schedule(dynamic)
      for (i=0; i< gd->nnodescanlevTable[cat1]; i++) {
#ifdef BALLS
          if (cmd->scanLevel==0)
              p = (nodeptr) roottable[cat1];
          else
              p = nodetablescanlev[cat1][i];
#else
          //B if not def BALLS chose one of these:
          p = (nodeptr) roottable[cat1];
//          p = nodetablescanlev[cat1][i];
          //E
#endif
          //B Set histograms to zero for the pivot
          for (n = 1; n <= cmd->sizeHistN; n++)
              hist.histNNSubthread[n] = 0.0;
//          CLRM_ext_ext(hist.histXithreadcos,
//                       cmd->mChebyshev+1, cmd->sizeHistN);
//          CLRM_ext_ext(hist.histXithreadsin,
//                       cmd->mChebyshev+1, cmd->sizeHistN);
          //E
          
          /*
          //B Set a reference axis guess for the pivot
#ifdef POLARAXIS
          hist.q0[0] = 0.0;
          hist.q0[1] = 0.0;
          hist.q0[2] = 1.0;
          DOTPSUBV(hist.drpq2, hist.dr0, Pos(p), hist.q0);
          hist.drpq = rsqrt(hist.drpq2);
          real b = 2.0*rasin(hist.drpq/2.0);
          hist.cosb = rcos(b);
          hist.sinb = rsin(b);
          if (hist.drpq2==0) continue;
#else
          dRotation3D(Pos(p), ROTANGLE, ROTANGLE, ROTANGLE, hist.q0);
          DOTPSUBV(hist.drpq2, hist.dr0, Pos(p), hist.q0);
          hist.drpq = rsqrt(hist.drpq2);
#endif
          */
          
          //E
          histb.active[0] = (nodeptr) (roottable[cat2]);
          walktree_balls6_omp(cmd, gd, histb.active, histb.active + 1,
                  histb.interact, histb.interact + histb.actlen,
                  p, Size(p), Pos(p), &histb, &hist, nbody[cat1], cat1);
          if (!scanopt(cmd->options, "no-two-balls")) {
              computeBodyProperties_sincos_kk(cmd, gd, (bodyptr)p,
                                                  nbody[cat1], &hist);
          }
      } // end do body i

#pragma omp critical
    {
/*        for (m=1; m<=cmd->mChebyshev+1; m++) {
            ADDM_ext(gd->histZetaMcos[m],gd->histZetaMcos[m],
                     hist.histZetaMthreadcos[m],cmd->sizeHistN);
            ADDM_ext(gd->histZetaMsin[m],gd->histZetaMsin[m],
                     hist.histZetaMthreadsin[m],cmd->sizeHistN);
            ADDM_ext(gd->histZetaMsincos[m],gd->histZetaMsincos[m],
                     hist.histZetaMthreadsincos[m],cmd->sizeHistN);
            ADDM_ext(gd->histZetaMcossin[m],gd->histZetaMcossin[m],
                     hist.histZetaMthreadcossin[m],cmd->sizeHistN);
        }
        */
        gd->nsmoothcount += hist.nsmoothcountthread;
        gd->nbbcalc += hist.nbbcalcthread;
        gd->nbccalc += hist.nbccalcthread;
        gd->ncccalc += hist.ncccalcthread;
#ifdef DEBUG
        ibfcount += hist.ibfcountthread;
#endif
    }

      //B smooth-pivot
      ipfalse += hist.ipfalsethreads;
      icountNbRmin += hist.icountNbRminthread;
      icountNbRminOverlap += hist.icountNbRminOverlapthread;
      //E

    search_free_sincos_omp_kk(cmd, gd, &hist);     // free memory
    search_free_omp_balls6_kk(cmd, gd, &histb);

  } // end pragma omp parallel

    //B smooth-pivot
    real xi, den, num;
    int mm;
    if (scanopt(cmd->options, "smooth-pivot")) {
        num = (real)nbody[cat1];
        den = (real)(nbody[cat1]-ipfalse);
        xi = num/den;
        verb_print(cmd->verbose,
                   "octree-kkk-omp: p falses found = %ld and %e %e %e\n",
                   ipfalse, num, den, xi);
/*        for (mm=1; mm<=cmd->mChebyshev+1; mm++) {
            MULMS_ext(gd->histZetaMcos[mm],
                      gd->histZetaMcos[mm],xi,cmd->sizeHistN);
            MULMS_ext(gd->histZetaMsin[mm],
                      gd->histZetaMsin[mm],xi,cmd->sizeHistN);
            MULMS_ext(gd->histZetaMsincos[mm],
                      gd->histZetaMsincos[mm],xi,cmd->sizeHistN);
            MULMS_ext(gd->histZetaMcossin[mm],
                      gd->histZetaMcossin[mm],xi,cmd->sizeHistN);
        } */
    }
    //E

    //B smooth-pivot
    verb_print(cmd->verbose,
               "octree-kkk-balls4-omp: p falses found = %ld\n",ipfalse);
    verb_print(cmd->verbose,
               "octree-kkk-balls4-omp: count NbRmin found = %ld\n",icountNbRmin);
    verb_print(cmd->verbose,
               "octree-kkk-balls4-omp: count overlap found = %ld\n",
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
    verb_print(cmd->verbose, "octree-kkk-balls4-omp: p falses found = %ld\n",
               ifalsecount);
    verb_print(cmd->verbose, "octree-kkk-balls4-omp: p true found = %ld\n",
               itruecount);
    verb_print(cmd->verbose, "octree-kkk-balls4-omp: total = %ld\n",
               itruecount+ifalsecount);
    //E

#ifdef DEBUG
    gd->nbodybf = cmd->ntosave;
    verb_print(cmd->verbose, "\noctree-kkk-balls4-omp: Total bodies found: %ld %ld\n",ibfcount, gd->nbodybf);
#endif

    verb_print(cmd->verbose,
               "octree-kkk-balls4-omp: nsmoothcount = %ld\n",gd->nsmoothcount);

    gd->cpusearch = CPUTIME - cpustart;
    verb_print(cmd->verbose, "Going out: CPU time = %lf\n",CPUTIME-cpustart);

    return SUCCESS;
}

//#ifndef BALLS
//B 2023.11.22
local bool nodes_condition_balls(struct cmdline_data* cmd, struct  global_data* gd,
                                  nodeptr p, nodeptr q, real *dr1, vector dr)
{
//    real drpq, drpq2;
/*
     DOTPSUBV(drpq2, dr, Pos(p), Pos(q));
 #ifdef PERIODIC
     VWrapAll(dr);
     DOTVP(drpq2, dr, dr);
 #endif
     drpq = rsqrt(drpq2);
     *dr1 = drpq;
*/
     int n;

     if ( *dr1 == 0.0)
         return (FALSE);
     else
         if ( (Radius(p)+Radius(q))/(*dr1) < gd->deltaR) {
             if (scanopt(cmd->options, "behavior-tree-omp")) {
//B To behaves as tree-omp
                 if ( (*dr1)<gd->Rcut ) {
                     if((*dr1)>cmd->rminHist) {
                         if (cmd->rminHist==0)
                             n = (int)(cmd->logHistBinsPD*(rlog10(*dr1) -
                                     rlog10(cmd->rangeN)) + cmd->sizeHistN) + 1;
                         else
                             n = (int)(rlog10((*dr1)/cmd->rminHist)
                                       * gd->i_deltaR) + 1;
                         if (n<=cmd->sizeHistN-1 && n>=1) {
                             if ( gd->deltaRV[n] < *dr1 - Radius(q) && *dr1
                                 + Radius(q) < gd->deltaRV[n+1]) {
                                 return (TRUE);
                             } else {
                                 return (FALSE);
                             }
                         } else
                             return (FALSE);
                     } else
                         return (FALSE);
                 } else
                     return (FALSE);
 //E
             } else { // ! behavior-tree-omp
                 return (TRUE);
             }
         } else
             return (FALSE);
}
//E
//#endif

local void walktree_balls6_omp(struct cmdline_data* cmd, struct  global_data* gd,
                               nodeptr *aptr, nodeptr *nptr,
            cellptr cptr, cellptr bptr,
            nodeptr p, real psize, vector pmid, gdhistptr_omp_balls6_kk histb,
            gdhistptr_sincos_omp_kk histsincos,
            INTEGER nbody, int ifile)
{
    nodeptr *np, *ap, q;
    int actsafe;
    real dr1;
    vector dr;

    if (Update(p)) {
        np = nptr;
        actsafe = histb->actlen - NSUB;
        for (ap = aptr; ap < nptr; ap++) {
            if (Type(*ap) == CELL) {
                if (!reject_cell_balls(cmd, gd, p, *ap, &dr1, dr)) {
                    if ( ((Nb(*ap)<=gd->nsmooth[0])
                            || (Size(*ap)<=gd->rminCell[0]))
                            && (dr1 > gd->rminCell[1]) ) {
                        if (np - histb->active >= actsafe) {
                            error("walktree (1): active list overflow\n");
                        }
                        histsincos->nsmoothcountthread += 1;
                        --bptr;
                        Mass(bptr) = Mass(*ap);
                        Kappa(bptr) = Kappa(*ap);
                        SETV(Pos(bptr), Pos(*ap));
                        Id(bptr) = Id(*ap);
                        Type(bptr) = Type(*ap);
                        Nb(bptr) = Nb(*ap);
                     } else { // ! bucket condition
                         if (nodes_condition_balls(cmd, gd, p, *ap, &dr1, dr)) {
                             if (!scanopt(cmd->options, "no-two-balls")
                                 && Type(p) == CELL ) {
                                 sumcellcell_balls6_omp(cmd, gd, (cellptr)(*ap),
                                        (cellptr)*ap+1, p, histsincos);
                             } else {
                                 if (np - histb->active >= actsafe)
                                     error("walktree (2): active list overflow\n");
                                 if (!scanopt(cmd->options, "no-one-ball")) {
                                     Mass(cptr) = Mass(*ap);
                                     Kappa(cptr) = Kappa(*ap);
                                     SETV(Pos(cptr), Pos(*ap));
                                     Id(cptr) = Id(*ap);
                                     Type(cptr) = Type(*ap);
                                     Nb(cptr) = Nb(*ap);
                                     cptr++;
                                 } else // options : ! no-one-ball
                                     for (q = More(*ap); q != Next(*ap);
                                          q = Next(q))
                                         *np++= q;
                             } // meet condition :: no-wo-balls
                         } else // First meet condition
                             for (q = More(*ap); q != Next(*ap); q = Next(q))
                                 *np++= q;
                     } // ! bucket condition
                 } // ! reject_cell
             } else  // ! == CELL
                 if (*ap != p) {
                     --bptr;
                     Mass(bptr) = Mass(*ap);
                     Kappa(bptr) = Kappa(*ap);
                     SETV(Pos(bptr), Pos(*ap));
                     Id(bptr) = Id(*ap);
                     Type(bptr) = Type(*ap);
                     Nb(bptr) = 1;
                 }
        } // ! loop for ap

        gd->actmax = MAX(gd->actmax, np - histb->active);
        if (np != nptr)
            walksub6_omp(cmd, gd, nptr, np, cptr, bptr, p, psize, pmid, histb,
                            histsincos, nbody, ifile);
        else {
            if (Type(p) != BODY)
                error("walktree: recursion terminated with cell\n");

            sum_balls6_omp(cmd, gd, cptr, bptr, (bodyptr) p,
                           histb, histsincos, nbody, ifile);
            Update(p) = FALSE;

#ifdef DEBUG
            pbf = bodytabbf + histsincos->ibfcountthread;
            Mass(pbf) = Mass(p);
            Kappa(pbf) = Kappa(p);
            SETV(Pos(pbf), Pos(p));
            histsincos->ibfcountthread += 1;
            Id(pbf) = histsincos->ibfcountthread;
#endif
        }
    }   // ! update
}

local void walksub6_omp(struct cmdline_data* cmd, struct  global_data* gd,
                        nodeptr *nptr, nodeptr *np, cellptr cptr, cellptr bptr,
        nodeptr p, real psize, vector pmid, gdhistptr_omp_balls6_kk histb,
        gdhistptr_sincos_omp_kk histsincos,
        INTEGER nbody, int ifile)
{
    real poff;
    nodeptr q;
    int k;
    vector nmid;

    poff = psize / 4;
    if (Type(p) == CELL) {
        for (q = More(p); q != Next(p); q = Next(q)) {
            for (k = 0; k < NDIM; k++)
                nmid[k] = pmid[k] + (Pos(q)[k] < pmid[k] ? - poff : poff);
            walktree_balls6_omp(cmd, gd, nptr, np, cptr, bptr, q, psize / 2, nmid,
            histb, histsincos, nbody, ifile);
        }
    } else {
        for (k = 0; k < NDIM; k++)
            nmid[k] = pmid[k] + (Pos(p)[k] < pmid[k] ? - poff : poff);
        walktree_balls6_omp(cmd, gd, nptr, np, cptr, bptr, p, psize / 2, nmid,
            histb, histsincos, nbody, ifile);
    }
}


local int sum_balls6_omp(struct cmdline_data* cmd, struct  global_data* gd,
                          cellptr cptr, cellptr bptr, bodyptr p,
                          gdhistptr_omp_balls6_kk histb,
                          gdhistptr_sincos_omp_kk histsincos,
                          INTEGER nbody, int ifile)
{
    int n;
    int m;
    real xi;
    local INTEGER ip;
    gdhist_sincos_omp_kk hist1sincos;
    
    //B smooth-pivot
    NbRmin(p) = 1;
    NbRminOverlap(p) = 0;
    KappaRmin(p) = Kappa(p);
    if (scanopt(cmd->options, "smooth-pivot")) {
        if (Update(p) == FALSE) {
            histsincos->ipfalsethreads++;
            return SUCCESS;
        }
    }
    //E

    //B init:
    search_init_sincos_omp_kk(cmd, gd, &hist1sincos, ifile);
    for (n = 1; n <= cmd->sizeHistN; n++) {
        hist1sincos.histNNSubthread[n] = 0.0;
    }
    for (n = 1; n <= cmd->sizeHistN; n++) {
        hist1sincos.histNNSubthread[n] = 0.0;
    }
//    CLRM_ext_ext(hist1sincos.histXithreadcos, cmd->mChebyshev+1,
//                 cmd->sizeHistN);
//    CLRM_ext_ext(hist1sincos.histXithreadsin, cmd->mChebyshev+1,
//                 cmd->sizeHistN);
    //E

    /*
    //B Set reference axis for p (pivot)
#ifdef POLARAXIS
    hist1sincos.q0[0] = 0.0;
    hist1sincos.q0[1] = 0.0;
    hist1sincos.q0[2] = 1.0;
    DOTPSUBV(hist1sincos.drpq2, hist1sincos.dr0, Pos(p), hist1sincos.q0);
    hist1sincos.drpq = rsqrt(hist1sincos.drpq2);
    real b = 2.0*rasin(hist1sincos.drpq/2.0);
    hist1sincos.cosb = rcos(b);
    hist1sincos.sinb = rsin(b);
//    if (hist.drpq2==0) continue;
#else
    dRotation3D(Pos(p), ROTANGLE, ROTANGLE, ROTANGLE, hist1sincos.q0);
    DOTPSUBV(hist1sincos.drpq2, hist1sincos.dr0, Pos(p), hist1sincos.q0);
    hist1sincos.drpq = rsqrt(hist1sincos.drpq2);
#endif
    //E
    */

    if (!scanopt(cmd->options, "no-one-ball"))
        sumcell_balls6_omp(cmd, gd, histb->interact, cptr, (bodyptr) p,
                           &hist1sincos);
    sumnode_balls6_omp(cmd, gd, bptr, histb->interact + histb->actlen,
                       (bodyptr) p, &hist1sincos);

    if (scanopt(cmd->options, "smooth-pivot")) {
        for (n = 1; n <= cmd->sizeHistN; n++) {
            hist1sincos.histNNSubthread[n] =
                          ((real)NbRmin(p))*hist1sincos.histNNSubthread[n];
        }
    }
    computeBodyProperties_sincos_kk_sum_balls6_omp(cmd, gd,
                                    p, nbody, histsincos, &hist1sincos);

    histsincos->nbbcalcthread += histb->interact + histb->actlen - bptr;
    histsincos->nbccalcthread += cptr - histb->interact;

    //B smooth-pivote
    histsincos->icountNbRminthread += NbRmin(p);
    histsincos->icountNbRminOverlapthread += NbRminOverlap(p);
    //E

    ip = p - bodytable[gd->iCatalogs[0]] + 1;
    if (ip%cmd->stepState == 0) {
        verb_log_print(cmd->verbose_log, gd->outlog, " - Completed pivot: %ld\n", ip);
    }

    search_free_sincos_omp_kk(cmd, gd, &hist1sincos);
}

local void sumnode_balls6_omp(struct cmdline_data* cmd, struct  global_data* gd,
                              cellptr start, cellptr finish, bodyptr p0,
                              gdhistptr_sincos_omp_kk hist)
{
    cellptr p;
#ifdef SINGLEP
    float dr1;
    float dr[NDIM];
#else
    vector dr;
    real dr1;
#endif
    nodeptr pb;
    INTEGER ibodycount=0;
    int n;
    real xi;
    real cosphi;
    real sinphi;

    for (p = start; p < finish; p++) {
        pb = ((nodeptr) p);
        if (accept_body(cmd, gd, p0, pb, &dr1, dr)) {
            //B smooth-pivot
            if (scanopt(cmd->options, "smooth-pivot")) {
                if (dr1<=gd->rsmooth[0]) {
                    if (Update(pb)==TRUE) {
                        Update(pb) = FALSE;
                        NbRmin(pb) += 1;
                        KappaRmin(pb) += Kappa(pb);
                    } else {
                        NbRminOverlap(pb) += 1;
                    }
                }
            }
            //E
            if(dr1>cmd->rminHist) {
                ibodycount++;
                if (cmd->rminHist==0)
                    n = (int)(cmd->logHistBinsPD*(rlog10(dr1) - rlog10(cmd->rangeN))
                              + cmd->sizeHistN) + 1;
                else
                    n = (int)(rlog10(dr1/cmd->rminHist) * gd->i_deltaR) + 1;
                if (n<=cmd->sizeHistN && n>=1) {
                    hist->histNNSubthread[n] = hist->histNNSubthread[n] + 1.;
                    xi = Kappa(pb);
                    
                    /*
                    //B Component of pb with respect to the axis of reference
#ifdef POLARAXIS
                real a, c2;
                vector vc;
                DOTPSUBV(c2, vc, Pos(pb), hist->q0);
#ifdef NOLIMBER
                a = 2.0*rasin(dr1/2.0);
#else
                a = dr1;
#endif
                real cosc;
                cosc = Pos(pb)[2];
                cosphi = (cosc - (1.0-0.5*rsqr(a))*hist->cosb)
                        /(a*hist->sinb);
                if (rabs(cosphi) <= 1.0)
                    sinphi = rsqrt(1.0 - rsqr(cosphi));
                else
                    sinphi = 0.0;
                if (!crossVecProdSign(Pos(p), hist->q0, Pos(pb)))
                    sinphi *= -1.0;
#else
                    real s, sy;
                    vector pr0;
                    DOTVP(s, dr, hist->dr0);
                    CROSSVP(pr0,hist->dr0,Pos(p));
                    DOTVP(sy, dr, pr0);
                    cosphi = s/((dr1)*hist->drpq);
                    sinphi = rsqrt(1.0 - rsqr(cosphi));;
                    if (sy < 0) sinphi *= -1.0;
                    if (rabs(cosphi)>1.0)
                        verb_log_print(cmd->verbose, gd->outlog,
                            "sumnode: Warning!... cossphi must be in (-1,1): %g\n",
                            cosphi);
#endif
                    //E
                    */
                    
//                    CHEBYSHEVTUOMPKKK;
                } // ! 1 < n < HistN
            } // ! dr1 > rmin
        } // ! accept_body
    } // ! loop p
}

local void sumcell_balls6_omp(struct cmdline_data* cmd, struct  global_data* gd,
                              cellptr start, cellptr finish, bodyptr p0,
                              gdhistptr_sincos_omp_kk hist)
{
    cellptr p;
#ifdef SINGLEP
    float dr1;
    float dr[NDIM];
#else
    vector dr;
    real dr1;
#endif
    nodeptr pb;
    INTEGER ibodycount=0;
    int n;
    real xi;
    real cosphi;
    real sinphi;

    for (p = start; p < finish; p++) {
        pb = ((nodeptr) p);
        if (accept_body(cmd, gd, p0, pb, &dr1, dr)) {
            if(dr1>cmd->rminHist) {
                ibodycount++;
                if (cmd->rminHist==0)
                    n = (int)(cmd->logHistBinsPD*(rlog10(dr1) - rlog10(cmd->rangeN))
                              + cmd->sizeHistN) + 1;
                else
                    n = (int)(rlog10(dr1/cmd->rminHist) * gd->i_deltaR) + 1;
                if (n<=cmd->sizeHistN && n>=1) {
                    hist->histNNSubthread[n] = hist->histNNSubthread[n] + 1.;
                    xi = Kappa(pb);
                    
                    /*
                    //B Component of pb with respect to the axis of reference
#ifdef POLARAXIS
                real a, c2;
                vector vc;
                DOTPSUBV(c2, vc, Pos(pb), hist->q0);
#ifdef NOLIMBER
                a = 2.0*rasin(dr1/2.0);
#else
                a = dr1;
#endif
                real cosc;
                cosc = Pos(pb)[2];
                cosphi = (cosc - (1.0-0.5*rsqr(a))*hist->cosb)
                        /(a*hist->sinb);
                if (rabs(cosphi) <= 1.0)
                    sinphi = rsqrt(1.0 - rsqr(cosphi));
                else
                    sinphi = 0.0;
                if (!crossVecProdSign(Pos(p), hist->q0, Pos(pb)))
                    sinphi *= -1.0;
#else
                    real s, sy;
                    vector pr0;
                    DOTVP(s, dr, hist->dr0);
                    CROSSVP(pr0,hist->dr0,Pos(p));
                    DOTVP(sy, dr, pr0);
                    cosphi = s/((dr1)*hist->drpq);
                    sinphi = rsqrt(1.0 - rsqr(cosphi));;
                    if (sy < 0) sinphi *= -1.0;
                    if (rabs(cosphi)>1.0)
                        verb_log_print(cmd->verbose, gd->outlog,
                            "sumcell: Warning!... cossphi must be in (-1,1): %g\n",
                            cosphi);
#endif
                    //E
                    */
                    
                    // CHEBYSHEVTUOMPKKK;
                } // ! 1 < n < HistN
            } // ! dr1 > rmin
        } // ! accept_body
    } // ! loop p
}

local void sumcellcell_balls6_omp(struct cmdline_data* cmd,
                                  struct  global_data* gd,
                                  cellptr start, cellptr finish, nodeptr p0,
//                                  gdhistptr_omp_balls6_kk histb,
                                  gdhistptr_sincos_omp_kk hist)
{
    cellptr p;
#ifdef SINGLEP
    float dr1;
    float dr[NDIM];
#else
    vector dr;
    real dr1;
#endif
    nodeptr pb;
    INTEGER ibodycount=0;
    int n;
    real xi;
    real cosphi;
    real sinphi;

    for (p = start; p < finish; p++) {
        pb = ((nodeptr) p);
        if (accept_body(cmd, gd, (bodyptr)p0, pb, &dr1, dr)) {
            if(dr1>cmd->rminHist) {
                ibodycount++;
                if (cmd->rminHist==0)
                    n = (int)(cmd->logHistBinsPD*(rlog10(dr1) - rlog10(cmd->rangeN))
                              + cmd->sizeHistN) + 1;
                else
                    n = (int)(rlog10(dr1/cmd->rminHist) * gd->i_deltaR) + 1;
                if (n<=cmd->sizeHistN && n>=1) {
                    hist->histNNSubthread[n] = hist->histNNSubthread[n] + Nb(pb);
                    xi = Kappa(pb);
                    
                    /*
                    //B Component of pb with respect to the axis of reference
#ifdef POLARAXIS
                real a, c2;
                vector vc;
                DOTPSUBV(c2, vc, Pos(pb), hist->q0);
#ifdef NOLIMBER
                a = 2.0*rasin(dr1/2.0);
#else
                a = dr1;
#endif
                real cosc;
                cosc = Pos(pb)[2];
                cosphi = (cosc - (1.0-0.5*rsqr(a))*hist->cosb)
                        /(a*hist->sinb);
                if (rabs(cosphi) <= 1.0)
                    sinphi = rsqrt(1.0 - rsqr(cosphi));
                else
                    sinphi = 0.0;
                if (!crossVecProdSign(Pos(p), hist->q0, Pos(pb)))
                    sinphi *= -1.0;
#else
                    real s, sy;
                    vector pr0;
                    DOTVP(s, dr, hist->dr0);
                    CROSSVP(pr0,hist->dr0,Pos(p));
                    DOTVP(sy, dr, pr0);
                    cosphi = s/((dr1)*hist->drpq);
                    sinphi = rsqrt(1.0 - rsqr(cosphi));;
                    if (sy < 0) sinphi *= -1.0;
                    if (rabs(cosphi)>1.0)
                        verb_log_print(cmd->verbose, gd->outlog,
                        "sumcellcell: Warning!... cossphi must be in (-1,1): %g\n",
                        cosphi);
#endif
                    //E
                    */
                    
                    
                    // CHEBYSHEVTUOMPKKK;
                } // ! 1 < n < HistN
            } // ! dr1 > rmin
        } // ! accept_body
    } // ! loop p

    hist->ncccalcthread += 1;
}


//B Routines as in cballsutils

local int search_init_gd_sincos_omp_kk(struct  cmdline_data* cmd,
                                        struct  global_data* gd)
{
    int n;
    int m;

    for (n = 1; n <= cmd->sizeHistN; n++)
        gd->histNNSub[n] = 0.0;
/*    for (m = 1; m <= cmd->mChebyshev+1; m++) {
        CLRM_ext(gd->histZetaMcos[m], cmd->sizeHistN);
        CLRM_ext(gd->histZetaMsin[m], cmd->sizeHistN);
        CLRM_ext(gd->histZetaMsincos[m], cmd->sizeHistN);
        CLRM_ext(gd->histZetaMcossin[m], cmd->sizeHistN);
        gd->histXi[m][n] = 0.0;
    } */
    gd->nbbcalc = gd->nbccalc = gd->ncccalc = gd->nsmoothcount = 0;

    return SUCCESS;
}

local int search_init_sincos_omp_kk(struct  cmdline_data* cmd,
                                     struct  global_data* gd,
                                     gdhistptr_sincos_omp_kk hist, int ifile)
{
    int n;
    int m;

//    hist->ChebsT = dvector(1,cmd->mChebyshev+1);
//    hist->ChebsU = dvector(1,cmd->mChebyshev+1);
    hist->histNNSubthread = dvector(1,cmd->sizeHistN);
//    hist->histXithreadcos = dmatrix(1,cmd->mChebyshev+1,1,cmd->sizeHistN);
//    hist->histXithreadsin = dmatrix(1,cmd->mChebyshev+1,1,cmd->sizeHistN);
//    hist->histZetaMthreadcos =
//            dmatrix3D(1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
//    hist->histZetaMthreadsin =
//            dmatrix3D(1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
//    hist->histZetaMthreadsincos =
//            dmatrix3D(1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
//    hist->histZetaMthreadcossin =
//            dmatrix3D(1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    
    /*
    hist->xiOUTVPcos = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
    hist->xiOUTVPsin = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
    hist->xiOUTVPsincos = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
    hist->xiOUTVPcossin = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
    hist->histZetaMtmpcos = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
    hist->histZetaMtmpsin = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
    hist->histZetaMtmpsincos = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
    hist->histZetaMtmpcossin = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
*/
    
    for (n = 1; n <= cmd->sizeHistN; n++)
        hist->histNNSubthread[n] = 0.0;
/*    for (m = 1; m <= cmd->mChebyshev+1; m++) {
        CLRM_ext(hist->histZetaMthreadcos[m], cmd->sizeHistN);
        CLRM_ext(hist->histZetaMthreadsin[m], cmd->sizeHistN);
        CLRM_ext(hist->histZetaMthreadsincos[m], cmd->sizeHistN);
        CLRM_ext(hist->histZetaMthreadcossin[m], cmd->sizeHistN);
    } */

    hist->nbbcalcthread = 0;
    hist->nbccalcthread = 0;
    hist->ncccalcthread = 0;
    hist->ibfcountthread = 0;
    hist->nsmoothcountthread = 0;
    //B smooth-pivot
    hist->ipfalsethreads = 0;
    hist->icountNbRminthread = 0;
    hist->icountNbRminOverlapthread = 0;
    //E

    return SUCCESS;
}

local int search_free_sincos_omp_kk(struct  cmdline_data* cmd,
                                         struct  global_data* gd,
                                      gdhistptr_sincos_omp_kk hist)
{
    /*
    free_dmatrix(hist->histZetaMtmpcossin,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix(hist->histZetaMtmpsincos,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix(hist->histZetaMtmpsin,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix(hist->histZetaMtmpcos,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix(hist->xiOUTVPcossin,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix(hist->xiOUTVPsincos,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix(hist->xiOUTVPsin,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix(hist->xiOUTVPcos,1,cmd->sizeHistN,1,cmd->sizeHistN);
    */
/*    free_dmatrix3D(hist->histZetaMthreadcossin,
                   1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix3D(hist->histZetaMthreadsincos,
                   1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix3D(hist->histZetaMthreadsin,
                   1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix3D(hist->histZetaMthreadcos,
                   1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix(hist->histXithreadsin,1,cmd->mChebyshev+1,1,cmd->sizeHistN);
    free_dmatrix(hist->histXithreadcos,1,cmd->mChebyshev+1,1,cmd->sizeHistN); */
    free_dvector(hist->histNNSubthread,1,cmd->sizeHistN);
/*    free_dvector(hist->ChebsU,1,cmd->mChebyshev+1);
    free_dvector(hist->ChebsT,1,cmd->mChebyshev+1);
    */

    return SUCCESS;
}

local int computeBodyProperties_sincos_kk(struct  cmdline_data* cmd,
                                            struct  global_data* gd,
                                            bodyptr p, int nbody,
                                            gdhistptr_sincos_omp_kk hist)
{
    int n;
    int m;
    real xi;

    xi = Weight(p)*Kappa(p)/nbody;
/*
    for (m=1; m<=cmd->mChebyshev+1; m++)
        for (n=1; n<=cmd->sizeHistN; n++) {
            hist->histXithreadcos[m][n] /= MAX(hist->histNNSubthread[n],1.0);
            hist->histXithreadsin[m][n] /= MAX(hist->histNNSubthread[n],1.0);
        }
    */
    /*
    for (m=1; m<=cmd->mChebyshev+1; m++) {
        OUTVP_ext(hist->xiOUTVPcos,
            hist->histXithreadcos[m], hist->histXithreadcos[m], cmd->sizeHistN);
        OUTVP_ext(hist->xiOUTVPsin,
            hist->histXithreadsin[m], hist->histXithreadsin[m],cmd->sizeHistN);
        OUTVP_ext(hist->xiOUTVPsincos,
            hist->histXithreadsin[m], hist->histXithreadcos[m],cmd->sizeHistN);
        OUTVP_ext(hist->xiOUTVPcossin,
            hist->histXithreadcos[m], hist->histXithreadsin[m],cmd->sizeHistN);
        CLRM_ext(hist->histZetaMtmpcos,cmd->sizeHistN);
        CLRM_ext(hist->histZetaMtmpsin,cmd->sizeHistN);
        CLRM_ext(hist->histZetaMtmpsincos,cmd->sizeHistN);
        CLRM_ext(hist->histZetaMtmpcossin,cmd->sizeHistN);
        MULMS_ext(hist->histZetaMtmpcos,hist->xiOUTVPcos,xi,cmd->sizeHistN);
        MULMS_ext(hist->histZetaMtmpsin,hist->xiOUTVPsin,xi,cmd->sizeHistN);
        MULMS_ext(hist->histZetaMtmpsincos,
            hist->xiOUTVPsincos,xi,cmd->sizeHistN);
        MULMS_ext(hist->histZetaMtmpcossin,
            hist->xiOUTVPcossin,xi,cmd->sizeHistN);
        ADDM_ext(hist->histZetaMthreadcos[m],
            hist->histZetaMthreadcos[m],hist->histZetaMtmpcos,cmd->sizeHistN);
        ADDM_ext(hist->histZetaMthreadsin[m],
            hist->histZetaMthreadsin[m],hist->histZetaMtmpsin,cmd->sizeHistN);
        ADDM_ext(hist->histZetaMthreadsincos[m],
            hist->histZetaMthreadsincos[m],
            hist->histZetaMtmpsincos,cmd->sizeHistN);
        ADDM_ext(hist->histZetaMthreadcossin[m],
            hist->histZetaMthreadcossin[m],
            hist->histZetaMtmpcossin,cmd->sizeHistN);
    }
    */

    return SUCCESS;
}

local int computeBodyProperties_sincos_kk_sum_balls6_omp(
                                struct  cmdline_data* cmd,
                                struct  global_data* gd,
                                bodyptr p, INTEGER nbody,
                                gdhistptr_sincos_omp_kk histsincos,
                                gdhistptr_sincos_omp_kk hist1sincos)
{
    int n;
    int m;
    real xi;

//#ifdef NONORMHIST
//    xi = Weight(p)*Kappa(p);
//#else
    xi = Weight(p)*Kappa(p)/nbody;
    if (scanopt(cmd->options, "smooth-pivot")) {
        xi = NbRmin(p)*Weight(p)*KappaRmin(p)/nbody;
    }
//#endif

    /*
    for (m=1; m<=cmd->mChebyshev+1; m++)
        
#ifdef NONORMHIST
        for (n=1; n<=cmd->sizeHistN; n++) {
            hist1sincos->histXithreadcos[m][n] /= 1.0;
            hist1sincos->histXithreadsin[m][n] /= 1.0;
        }
#else
//        for (n=1; n<=cmd->sizeHistN; n++) {
//            hist1sincos->histXithreadcos[m][n] /=
//                    MAX(hist1sincos->histNNSubthread[n],1.0);
//            hist1sincos->histXithreadsin[m][n] /=
//                    MAX(hist1sincos->histNNSubthread[n],1.0);
//        }
#endif
    */

    /*
    for (m=1; m<=cmd->mChebyshev+1; m++){
        OUTVP_ext(hist1sincos->xiOUTVPcos, hist1sincos->histXithreadcos[m],
                  hist1sincos->histXithreadcos[m], cmd->sizeHistN);
        OUTVP_ext(hist1sincos->xiOUTVPsin, hist1sincos->histXithreadsin[m],
                  hist1sincos->histXithreadsin[m],cmd->sizeHistN);
        OUTVP_ext(hist1sincos->xiOUTVPsincos, hist1sincos->histXithreadsin[m],
                  hist1sincos->histXithreadcos[m],cmd->sizeHistN);
        CLRM_ext(hist1sincos->histZetaMtmpcos,cmd->sizeHistN);
        CLRM_ext(hist1sincos->histZetaMtmpsin,cmd->sizeHistN);
        CLRM_ext(hist1sincos->histZetaMtmpsincos,cmd->sizeHistN);
        MULMS_ext(hist1sincos->histZetaMtmpcos,hist1sincos->xiOUTVPcos,
                  xi,cmd->sizeHistN);
        MULMS_ext(hist1sincos->histZetaMtmpsin,hist1sincos->xiOUTVPsin,
                  xi,cmd->sizeHistN);
        MULMS_ext(hist1sincos->histZetaMtmpsincos,hist1sincos->xiOUTVPsincos,
                  xi,cmd->sizeHistN);
        MULMS_ext(hist1sincos->histZetaMtmpcossin,hist1sincos->xiOUTVPcossin,
                  xi,cmd->sizeHistN);
        ADDM_ext(histsincos->histZetaMthreadcos[m],
                 histsincos->histZetaMthreadcos[m],
                 hist1sincos->histZetaMtmpcos,cmd->sizeHistN);
        ADDM_ext(histsincos->histZetaMthreadsin[m],
                 histsincos->histZetaMthreadsin[m],
                 hist1sincos->histZetaMtmpsin,cmd->sizeHistN);
        ADDM_ext(histsincos->histZetaMthreadsincos[m],
                 histsincos->histZetaMthreadsincos[m],
                 hist1sincos->histZetaMtmpsincos,cmd->sizeHistN);
        ADDM_ext(histsincos->histZetaMthreadcossin[m],
                 histsincos->histZetaMthreadcossin[m],
                 hist1sincos->histZetaMtmpcossin,cmd->sizeHistN);
    }
    */

    for (n = 1; n <= cmd->sizeHistN; n++) {
        histsincos->histNNSubthread[n] += hist1sincos->histNNSubthread[n];
    }

    return SUCCESS;
}

//E Routines as in cballsutils


local int search_init_omp_balls6_kk(struct cmdline_data* cmd,
                                 struct  global_data* gd,
                                 gdhistptr_omp_balls6_kk hist, int ifile)
{
#  define FACTIVE  0.75
//#  define FACTOR  1
#  define FACTOR  316
//#  define FACTOR  1024

    hist->actlen = FACTIVE * 216 * FACTOR * gd->tdepthTable[ifile];
    hist->actlen = hist->actlen * rpow(cmd->theta, -2.5);
    verb_log_print(cmd->verbose_log, gd->outlog,
                   "searchcalc_balls6: actlen=%d\n",hist->actlen);
    hist->active = (nodeptr *) allocate(hist->actlen * sizeof(nodeptr));
    gd->bytes_tot += hist->actlen*sizeof(nodeptr);
    verb_log_print(cmd->verbose_log, gd->outlog,
                   "\nAllocated %g MByte for active list storage.\n",
    hist->actlen*sizeof(nodeptr)/(1024.0*1024.0));
    hist->interact = (cellptr) allocate(hist->actlen * sizeof(cell));
    gd->bytes_tot += hist->actlen*sizeof(cell);
    verb_log_print(cmd->verbose_log, gd->outlog,
                   "Allocated %g MByte for interact list storage.\n",
                   hist->actlen*sizeof(cell)/(1024.0*1024.0));

#undef FACTOR
#undef FACTIVE

    return SUCCESS;
}

local int search_free_omp_balls6_kk(struct cmdline_data* cmd,
                                  struct  global_data* gd,
                                  gdhistptr_omp_balls6_kk hist)
{
    free(hist->interact);
    free(hist->active);

    return SUCCESS;
}

local int print_info(struct cmdline_data* cmd,
                                  struct  global_data* gd)
{
    verb_print(cmd->verbose,
               "searchcalc_normal: Running octree-kkk-balls4 ... (tree-omp) \n");
    verb_print(cmd->verbose, "treenode with 2 balls using balls4 method...\n");
    verb_print(cmd->verbose,
            "finding at the same time lists of neighbour cells and bodies...\n");

    if (cmd->usePeriodic==TRUE)
        error("CheckParameters: can´t have periodic boundaries and OCTREEKKKBALLSOMP definition (usePeriodic=%d)\nSet usePeriodic=false\n",
            cmd->usePeriodic);
    if (cmd->useLogHist==FALSE)
        error("CheckParameters: can´t have normal scale hist and OCTREEKKKBALLSOMP definition (useLogHist=%d)\nSet useLogHist=true\n",
            cmd->useLogHist);
    if (cmd->computeTPCF==FALSE)
        error("CheckParameters: can´t have computeTPCF=false and OCTREEKKKBALLSOMP definition (computeTPCF=%d)\nSet computeTPCF=true\n",
            cmd->computeTPCF);
#if NDIM == 2
    error("CheckParameters: OCTREEKKKBALLS4OMP definition works only in a 3D unit sphere");
#endif
#ifdef POLARAXIS
    verb_print(cmd->verbose, "with POLARAXIS... \n");
#endif
#ifdef NOLIMBER
    verb_print(cmd->verbose, "with NOLIMBER (no Limber aproximation)... \n");
#endif
    if (scanopt(cmd->options, "no-one-ball"))
        verb_print(cmd->verbose, "with option no-one-ball... \n");
    if (scanopt(cmd->options, "no-two-balls"))
        verb_print(cmd->verbose, "with option no-two-balls... \n");
    if (scanopt(cmd->options, "behavior-ball"))
        verb_print(cmd->verbose, "with option behavior-ball... \n");
    if (scanopt(cmd->options, "smooth-pivot"))
        verb_print(cmd->verbose,
                   "with option smooth-pivot... rsmooth=%g\n",gd->rsmooth[0]);
    if (scanopt(cmd->options, "bh86"))
        verb_print(cmd->verbose, "with cell opening criterion bh86... \n");
    if (scanopt(cmd->options, "sw94"))
        verb_print(cmd->verbose, "with cell opening criterion sw94... \n");

    return SUCCESS;
}

