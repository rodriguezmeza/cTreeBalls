/* ==============================================================================
 MODULE: search_octree_kkk_balls4_omp_triangles.c		[cTreeBalls]
 Written by: M.A. Rodriguez-Meza
 Starting date:    april 2023
 Purpose: 3-point correlation function computation
 Language: C
 Use: searchcalc_octree_kkk_balls4_omp(cmd, gd, btable, nbody,
                                           ipmin, ipmax, cat1, cat2);
 Major revisions:
 ==============================================================================*/
//        1          2          3          4        ^ 5          6          7


// Work to do in order to use with boxes not centered at (0,0,...)

#include "globaldefs.h"

#define CHEBYSHEVTUOMPANYMMTC                                       \
{real xicosmphi,xisinmphi; int m;                                   \
    hist->ChebsT[1] = 1.0;                                          \
    xicosmphi = xi*xih * hist->ChebsT[1];                           \
    hist->histXithreadcosTC[1] = xicosmphi;                         \
    hist->ChebsT[2] = cosphi;                                       \
    xicosmphi = xi*xih * hist->ChebsT[2];                           \
    hist->histXithreadcosTC[2] = xicosmphi;                         \
    hist->ChebsT[3] = 2.0*(cosphi)*(cosphi) - (1.0);                \
    xicosmphi = xi*xih * hist->ChebsT[3];                           \
    hist->histXithreadcosTC[3] = xicosmphi;                         \
    hist->ChebsU[1] = 0.0;                                          \
    xisinmphi = xi*xih * hist->ChebsU[1] * sinphi;                  \
    hist->histXithreadsinTC[1] = xisinmphi;                         \
    hist->ChebsU[2] = 1.0;                                          \
    xisinmphi = xi*xih * hist->ChebsU[2] * sinphi;                  \
    hist->histXithreadsinTC[2] = xisinmphi;                         \
    hist->ChebsU[3] = 2.0*cosphi;                                   \
    xisinmphi = xi*xih * hist->ChebsU[3] * sinphi;                  \
    hist->histXithreadsinTC[3] = xisinmphi;                         \
    for (m=4; m<=cmd->mChebyshev+1; m++){                           \
        hist->ChebsT[m] = 2.0*(cosphi)*hist->ChebsT[m-1] - hist->ChebsT[m-2]; \
        xicosmphi = xi*xih * hist->ChebsT[m];                       \
        hist->histXithreadcosTC[m] = xicosmphi;                     \
        hist->ChebsU[m] = 2.0*(cosphi)*hist->ChebsU[m-1] - hist->ChebsU[m-2]; \
        xisinmphi = xi*xih * hist->ChebsU[m] * sinphi;              \
        hist->histXithreadsinTC[m] = xisinmphi;                     \
    }}

//B Define structures:
typedef struct {
    realptr histNNSub;
    real ***histZetaMcos;
    real ***histZetaMsin;
    real ***histZetaM;
} gdl_sincos_omp_kkk, *gdlptr_sincos_omp_kkk;

typedef struct {
    INTEGER actlen;
    nodeptr *active;
    cellptr interact;
} gdhist_omp_balls6_kkk, *gdhistptr_omp_balls6_kkk;

typedef struct {
    real **histZetaMtmpcos;
    real **histZetaMtmpsin;
    real *ChebsT;
    real *ChebsU;
    real ***histZetaMthreadcos;
    real ***histZetaMthreadsin;
    
    real *histXithreadcosTC;
    real *histXithreadsinTC;
    real ***xiOUTVPcosTC;
    real ***xiOUTVPsinTC;
    real **histNNSubthreadTC;

    INTEGER nbbcalcthread;
    INTEGER nbccalcthread;
    INTEGER ncccalcthread;
    INTEGER ibfcountthread;
    INTEGER nsmoothcountthread;

    vector q0;
    real drpq2, drpq;
    vector dr0;
    real cosb;
    real sinb;

    INTEGER ipcount;
} gdhist_sincos_omp_kkk, *gdhistptr_sincos_omp_kkk;
//E Define structures:

local bool nodes_condition_balls(struct cmdline_data* cmd, struct  global_data* gd,
                                 nodeptr p, nodeptr q, real *dr1, vector dr);

local void walktree_balls6_omp(struct cmdline_data* cmd, struct  global_data* gd,
                               nodeptr *, nodeptr *, cellptr, cellptr,
                               nodeptr, real, vector,
                               gdhistptr_omp_balls6_kkk,
                               gdhistptr_sincos_omp_kkk, INTEGER, int);
local void walksub6_omp(struct cmdline_data* cmd, struct  global_data* gd,
                        nodeptr *, nodeptr *, cellptr, cellptr,
                        nodeptr, real, vector, gdhistptr_omp_balls6_kkk,
                        gdhistptr_sincos_omp_kkk, INTEGER, int);

local int sum_balls6_omp(struct cmdline_data* cmd, struct  global_data* gd,
                          cellptr, cellptr, bodyptr, gdhistptr_omp_balls6_kkk,
                          gdhistptr_sincos_omp_kkk, INTEGER, int);

local void sumnode_balls6_omp(struct cmdline_data* cmd,
                                        struct  global_data* gd,
                              cellptr start, cellptr finish, bodyptr p0,
                                        gdhistptr_sincos_omp_kkk hist);
local void sumcell_balls6_omp(struct cmdline_data* cmd,
                                        struct  global_data* gd,
                              cellptr, cellptr, bodyptr,
                              gdhistptr_sincos_omp_kkk);
local void sumcellcell_balls6_omp(struct cmdline_data* cmd,
                                  struct  global_data* gd,
                                  cellptr, cellptr, nodeptr,
                                  gdhistptr_sincos_omp_kkk);

local int search_init_gd_sincos_omp_kkk(struct  cmdline_data* cmd,
                                        struct  global_data* gd,
                                        gdlptr_sincos_omp_kkk);
local int search_free_gd_sincos_omp_kkk(struct  cmdline_data* cmd,
                                         struct  global_data* gd,
                                        gdlptr_sincos_omp_kkk);
local int search_init_sincos_omp_kkk(struct  cmdline_data* cmd,
                                  struct  global_data* gd,
                                     gdhistptr_sincos_omp_kkk hist, int);
local int search_free_sincos_omp_kkk(struct  cmdline_data* cmd,
                                         struct  global_data* gd,
                                     gdhistptr_sincos_omp_kkk hist);
local int computeBodyProperties_sincos_kkk(struct  cmdline_data* cmd,
                                            struct  global_data* gd,
                                            bodyptr p, int nbody,
                                       gdhistptr_sincos_omp_kkk hist);
local int computeBodyProperties_sincos_kkk_sum_balls6_omp(
                                            struct  cmdline_data*,
                                            struct  global_data*,
                                            bodyptr, INTEGER,
                                            gdhistptr_sincos_omp_kkk,
                                            gdhistptr_sincos_omp_kkk);

local int search_init_omp_balls6_kkk(struct cmdline_data* cmd,
                                         struct  global_data* gd,
                                        gdhistptr_omp_balls6_kkk hist, int);
local int search_free_omp_balls6_kkk(struct cmdline_data* cmd,
                                  struct  global_data* gd,
                                 gdhistptr_omp_balls6_kkk hist);

local int print_info(struct cmdline_data* cmd,
                     struct  global_data* gd);

//B Saving histograms section: case KKKCORRELATION:
local int PrintHistrBins(struct  cmdline_data* cmd, struct  global_data* gd);
local int PrintHistZetaM_sincos(struct  cmdline_data* cmd,
                                struct  global_data* gd,
                                gdlptr_sincos_omp_kkk);
local int PrintHistZetaMm_sincos(struct  cmdline_data* cmd,
                               struct  global_data* gd,
                                 gdlptr_sincos_omp_kkk);
//E Saving histograms section: case KKKCORRELATION:


/*
 Search routine using tree brute force direct method:

 To be called using: search=octree-kkk-balls4-omp-triangles

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
global int searchcalc_octree_kkk_balls4_omp_triangles(struct cmdline_data* cmd,
                                             struct  global_data* gd,
                                             bodyptr *btable, INTEGER *nbody,
                                             INTEGER ipmin, INTEGER *ipmax, 
                                             int cat1, int cat2)
{
    string routine_name = "searchcalc_octree_kkk_balls4_omp_triangles";
    double cpustart;
    gdl_sincos_omp_kkk gdl;

    cpustart = CPUTIME;
    print_info(cmd, gd);

#ifdef OPENMPCODE
    ThreadCount(cmd, gd, nbody[cat1], cat1);
#endif

    search_init_gd_sincos_omp_kkk(cmd, gd, &gdl);

    verb_print(cmd->verbose,
        "\n%s: Total allocated %g MByte storage so far.\n",
               routine_name, gd->bytes_tot/(1024.0*1024.0));

#pragma omp parallel default(none)                                          \
    shared(cmd,gd,btable,nbody,roottable,                                   \
           ipmin,ipmax,cat1,cat2,nodetablescanlev,nodetablescanlev_root, gdl)
  {
    nodeptr p;
      int n;
      int m;
      int i;

    //B init:
    gdhist_omp_balls6_kkk histb;
    gdhist_sincos_omp_kkk hist;
    search_init_omp_balls6_kkk(cmd, gd, &histb, cat1);
    search_init_sincos_omp_kkk(cmd, gd, &hist, cat1);
    //E

#pragma omp for nowait schedule(dynamic)
      for (i=0; i< gd->nnodescanlevTable[cat1]; i++) {
#ifdef BALLS
          if (cmd->scanLevel==0)
              p = (nodeptr) roottable[cat1];
          else
              p = nodetablescanlev[cat1][i];
#else
#error needs `BALLS`, switch it on addons/Makefile_addons_settings
#endif
          //B Set histograms to zero for the pivot
          CLRM_ext(hist.histNNSubthreadTC, cmd->sizeHistN);
          CLRV_ext(hist.histXithreadcosTC, cmd->mChebyshev+1);
          CLRV_ext(hist.histXithreadsinTC, cmd->mChebyshev+1);
          for (m=1; m<=cmd->mChebyshev+1; m++) {
              CLRM_ext(hist.xiOUTVPcosTC[m], cmd->sizeHistN);
              CLRM_ext(hist.xiOUTVPsinTC[m], cmd->sizeHistN);
          }
          //E

          //B Set a reference axis guess for the pivot
          dRotation3D(Pos(p), ROTANGLE, ROTANGLE, ROTANGLE, hist.q0);
          DOTPSUBV(hist.drpq2, hist.dr0, Pos(p), hist.q0);
          hist.drpq = rsqrt(hist.drpq2);
          //E
          histb.active[0] = (nodeptr) (roottable[cat2]);
          walktree_balls6_omp(cmd, gd, histb.active, histb.active + 1,
                  histb.interact, histb.interact + histb.actlen,
                  p, Size(p), Pos(p), &histb, &hist, nbody[cat1], cat1);
          if (!scanopt(cmd->options, "no-two-balls")) {
              computeBodyProperties_sincos_kkk(cmd, gd, (bodyptr)p,
                                                  nbody[cat1], &hist);
          }
      } // end do body i

#pragma omp critical
    {
        for (m=1; m<=cmd->mChebyshev+1; m++) {
            ADDM_ext(gdl.histZetaMcos[m],gdl.histZetaMcos[m],
                     hist.histZetaMthreadcos[m],cmd->sizeHistN);
            ADDM_ext(gdl.histZetaMsin[m],gdl.histZetaMsin[m],
                     hist.histZetaMthreadsin[m],cmd->sizeHistN);
        }
        gd->nsmoothcount += hist.nsmoothcountthread;
        gd->nbbcalc += hist.nbbcalcthread;
        gd->nbccalc += hist.nbccalcthread;
        gd->ncccalc += hist.ncccalcthread;
    } // ! critical

    search_free_sincos_omp_kkk(cmd, gd, &hist);     // free memory
    search_free_omp_balls6_kkk(cmd, gd, &histb);

  } // end pragma omp parallel

    verb_print(cmd->verbose,
               "octree-kkk-balls4-omp: nsmoothcount = %ld\n",gd->nsmoothcount);

// ===============================================
//B Saving histograms section: case KKKCORRELATION:
// ===============================================
    verb_print(cmd->verbose,
            "\n\tsearch_octree_kkk_omp: printing octree-kkk-omp method...\n\n");
    PrintHistrBins(cmd, gd);
    PrintHistZetaM_sincos(cmd, gd, &gdl);

    if (scanopt(cmd->options, "out-m-HistZeta"))
        PrintHistZetaMm_sincos(cmd, gd, &gdl);

    gd->flagPrint = FALSE;
// ===============================================
//E Saving histograms section: case KKKCORRELATION
// ===============================================

    search_free_gd_sincos_omp_kkk(cmd, gd, &gdl); // free memory

    gd->cpusearch = CPUTIME - cpustart;
    verb_print(cmd->verbose, "\nGoing out: CPU time = %lf %s\n",
               CPUTIME-cpustart, PRNUNITOFTIMEUSED);

    return SUCCESS;
}

local bool nodes_condition_balls(struct cmdline_data* cmd,
                                 struct  global_data* gd,
                                 nodeptr p, nodeptr q, real *dr1, vector dr)
{
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

local void walktree_balls6_omp(struct cmdline_data* cmd, struct  global_data* gd,
                               nodeptr *aptr, nodeptr *nptr,
            cellptr cptr, cellptr bptr,
            nodeptr p, real psize, vector pmid, gdhistptr_omp_balls6_kkk histb,
            gdhistptr_sincos_omp_kkk histsincos,
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
                                 sumcellcell_balls6_omp(cmd, gd,
                                        (cellptr)(*ap),
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
        }
    }   // ! update
}

local void walksub6_omp(struct cmdline_data* cmd, struct  global_data* gd,
                        nodeptr *nptr, nodeptr *np, cellptr cptr, cellptr bptr,
        nodeptr p, real psize, vector pmid, gdhistptr_omp_balls6_kkk histb,
        gdhistptr_sincos_omp_kkk histsincos,
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
                          gdhistptr_omp_balls6_kkk histb,
                          gdhistptr_sincos_omp_kkk histsincos,
                          INTEGER nbody, int ifile)
{
    string routine_name = "sum_balls6_omp";
    int n;
    int m;
    real xi;
    local INTEGER ip;
    gdhist_sincos_omp_kkk hist1sincos;

    if (cmd->verbose_log>=3) {
        verb_log_print(cmd->verbose, gd->outlog,
                       "%s: found %ld bodies to compute...\n",
                       routine_name, histb->interact + histb->actlen - bptr);
        verb_log_print(cmd->verbose, gd->outlog,
                       "%s: found %ld cells to compute...\n",
                       routine_name, cptr - histb->interact);
    }

    //B init:
    search_init_sincos_omp_kkk(cmd, gd, &hist1sincos, ifile);
    CLRM_ext(hist1sincos.histNNSubthreadTC, cmd->sizeHistN);
    CLRV_ext(hist1sincos.histXithreadcosTC, cmd->mChebyshev+1);
    CLRV_ext(hist1sincos.histXithreadsinTC, cmd->mChebyshev+1);
    for (m=1; m<=cmd->mChebyshev+1; m++) {
        CLRM_ext(hist1sincos.xiOUTVPcosTC[m], cmd->sizeHistN);
        CLRM_ext(hist1sincos.xiOUTVPsinTC[m], cmd->sizeHistN);
    }
    //E

    //B Set reference axis for p (pivot)
    dRotation3D(Pos(p), ROTANGLE, ROTANGLE, ROTANGLE, hist1sincos.q0);
    DOTPSUBV(hist1sincos.drpq2, hist1sincos.dr0, Pos(p), hist1sincos.q0);
    hist1sincos.drpq = rsqrt(hist1sincos.drpq2);
    //E

    if (!scanopt(cmd->options, "no-one-ball"))
        sumcell_balls6_omp(cmd, gd, histb->interact, cptr, (bodyptr) p,
                           &hist1sincos);
    sumnode_balls6_omp(cmd, gd, bptr, histb->interact + histb->actlen,
                       (bodyptr) p, &hist1sincos);

    computeBodyProperties_sincos_kkk_sum_balls6_omp(cmd, gd,
                                    p, nbody, histsincos, &hist1sincos);

    histsincos->nbbcalcthread += histb->interact + histb->actlen - bptr;
    histsincos->nbccalcthread += cptr - histb->interact;

    ip = p - bodytable[gd->iCatalogs[0]] + 1;
    if (ip%cmd->stepState == 0) {
        verb_log_print(cmd->verbose_log, gd->outlog, " - Completed pivot: %ld\n", ip);
    }

    search_free_sincos_omp_kkk(cmd, gd, &hist1sincos);
}

//B Routines as in cballsutils

local int search_init_gd_sincos_omp_kkk(struct  cmdline_data* cmd,
                                        struct  global_data* gd,
                                        gdlptr_sincos_omp_kkk gdl)
{
    int n;
    int m;
    INTEGER bytes_tot_local=0;

    gdl->histNNSub = dvector(1,cmd->sizeHistN);
    bytes_tot_local += 1*cmd->sizeHistN*sizeof(real);

    gdl->histZetaMcos =
            dmatrix3D(1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    gdl->histZetaMsin =
            dmatrix3D(1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    gdl->histZetaM =
            dmatrix3D(1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    bytes_tot_local +=
            5*(cmd->mChebyshev+1)*cmd->sizeHistN*cmd->sizeHistN*sizeof(real);

    gd->bytes_tot += bytes_tot_local;

    verb_print(cmd->verbose,
    "\nsearch_init_gd_octree_kkk: Allocated %g MByte for histograms storage.\n",
    bytes_tot_local*INMB);

    for (n = 1; n <= cmd->sizeHistN; n++)
        gd->histNNSub[n] = 0.0;
    for (m = 1; m <= cmd->mChebyshev+1; m++) {
        CLRM_ext(gdl->histZetaMcos[m], cmd->sizeHistN);
        CLRM_ext(gdl->histZetaMsin[m], cmd->sizeHistN);
    }
    gd->nbbcalc = gd->nbccalc = gd->ncccalc = gd->nsmoothcount = 0;

    return SUCCESS;
}

local int search_free_gd_sincos_omp_kkk(struct  cmdline_data* cmd,
                                         struct  global_data* gd,
                                        gdlptr_sincos_omp_kkk gdl)
{
    free_dmatrix3D(gdl->histZetaM,
                   1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix3D(gdl->histZetaMsin,
                   1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix3D(gdl->histZetaMcos,
                   1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dvector(gdl->histNNSub,1,cmd->sizeHistN);

    return SUCCESS;
}

local int search_init_sincos_omp_kkk(struct  cmdline_data* cmd,
                                     struct  global_data* gd,
                                     gdhistptr_sincos_omp_kkk hist, int ifile)
{
    int n;
    int m;

    hist->ChebsT = dvector(1,cmd->mChebyshev+1);
    hist->ChebsU = dvector(1,cmd->mChebyshev+1);

    hist->histNNSubthreadTC = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
    hist->histXithreadcosTC = dvector(1,cmd->mChebyshev+1);
    hist->histXithreadsinTC = dvector(1,cmd->mChebyshev+1);
    hist->xiOUTVPcosTC =
            dmatrix3D(1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    hist->xiOUTVPsinTC =
            dmatrix3D(1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);

    hist->histZetaMthreadcos =
            dmatrix3D(1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    hist->histZetaMthreadsin =
            dmatrix3D(1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);

    hist->histZetaMtmpcos = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
    hist->histZetaMtmpsin = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);

    for (m = 1; m <= cmd->mChebyshev+1; m++) {
        CLRM_ext(hist->histZetaMthreadcos[m], cmd->sizeHistN);
        CLRM_ext(hist->histZetaMthreadsin[m], cmd->sizeHistN);
    }

    hist->nbbcalcthread = 0;
    hist->nbccalcthread = 0;
    hist->ncccalcthread = 0;
    hist->ibfcountthread = 0;
    hist->nsmoothcountthread = 0;

    return SUCCESS;
}

local int search_free_sincos_omp_kkk(struct  cmdline_data* cmd,
                                         struct  global_data* gd,
                                      gdhistptr_sincos_omp_kkk hist)
{
    free_dmatrix(hist->histZetaMtmpsin,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix(hist->histZetaMtmpcos,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix3D(hist->histZetaMthreadsin,
                   1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix3D(hist->histZetaMthreadcos,
                   1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);

    free_dmatrix3D(hist->xiOUTVPsinTC,
                   1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix3D(hist->xiOUTVPcosTC,
                   1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dvector(hist->histXithreadsinTC,1,cmd->mChebyshev+1);
    free_dvector(hist->histXithreadcosTC,1,cmd->mChebyshev+1);

    free_dmatrix(hist->histNNSubthreadTC,1,cmd->sizeHistN,1,cmd->sizeHistN);

    free_dvector(hist->ChebsU,1,cmd->mChebyshev+1);
    free_dvector(hist->ChebsT,1,cmd->mChebyshev+1);

    return SUCCESS;
}

local int computeBodyProperties_sincos_kkk(struct  cmdline_data* cmd,
                                            struct  global_data* gd,
                                            bodyptr p, int nbody,
                                            gdhistptr_sincos_omp_kkk hist)
{
    int n;
    int m;
    real xi;

    xi = Weight(p)*Kappa(p)/nbody;

    for (m=1; m<=cmd->mChebyshev+1; m++) {
        for (n=1; n<=cmd->sizeHistN; n++) {
            for (int l=1; l<=cmd->sizeHistN; l++) {
                hist->xiOUTVPcosTC[m][n][l]
                /= MAX(hist->histNNSubthreadTC[n][l],1.0);
                hist->xiOUTVPsinTC[m][n][l]
                /= MAX(hist->histNNSubthreadTC[n][l],1.0);
            }
        }
    }

    for (m=1; m<=cmd->mChebyshev+1; m++) {
        CLRM_ext(hist->histZetaMtmpcos,cmd->sizeHistN);
        CLRM_ext(hist->histZetaMtmpsin,cmd->sizeHistN);
        MULMS_ext(hist->histZetaMtmpcos,hist->xiOUTVPcosTC[m],xi,cmd->sizeHistN);
        MULMS_ext(hist->histZetaMtmpsin,hist->xiOUTVPsinTC[m],xi,cmd->sizeHistN);
        ADDM_ext(hist->histZetaMthreadcos[m],
            hist->histZetaMthreadcos[m],hist->histZetaMtmpcos,cmd->sizeHistN);
        ADDM_ext(hist->histZetaMthreadsin[m],
            hist->histZetaMthreadsin[m],hist->histZetaMtmpsin,cmd->sizeHistN);
    }

    return SUCCESS;
}

local int computeBodyProperties_sincos_kkk_sum_balls6_omp(
                                struct  cmdline_data* cmd,
                                struct  global_data* gd,
                                bodyptr p, INTEGER nbody,
                                gdhistptr_sincos_omp_kkk histsincos,
                                gdhistptr_sincos_omp_kkk hist1sincos)
{
    int n;
    int m;
    real xi;

    xi = Weight(p)*Kappa(p)/nbody;

    for (m=1; m<=cmd->mChebyshev+1; m++) {
        for (n=1; n<=cmd->sizeHistN; n++) {
            for (int l=1; l<=cmd->sizeHistN; l++) {
                hist1sincos->xiOUTVPcosTC[m][n][l]
                        /= MAX(hist1sincos->histNNSubthreadTC[n][l],1.0);
                hist1sincos->xiOUTVPsinTC[m][n][l]
                        /= MAX(hist1sincos->histNNSubthreadTC[n][l],1.0);
            }
        }
    }

    for (m=1; m<=cmd->mChebyshev+1; m++) {
        CLRM_ext(hist1sincos->histZetaMtmpcos,cmd->sizeHistN);
        CLRM_ext(hist1sincos->histZetaMtmpsin,cmd->sizeHistN);
        MULMS_ext(hist1sincos->histZetaMtmpcos,
                  hist1sincos->xiOUTVPcosTC[m],xi,cmd->sizeHistN);
        MULMS_ext(hist1sincos->histZetaMtmpsin,
                  hist1sincos->xiOUTVPsinTC[m],xi,cmd->sizeHistN);
        ADDM_ext(histsincos->histZetaMthreadcos[m],
                 histsincos->histZetaMthreadcos[m],
                 hist1sincos->histZetaMtmpcos,cmd->sizeHistN);
        ADDM_ext(histsincos->histZetaMthreadsin[m],
                 histsincos->histZetaMthreadsin[m],
                 hist1sincos->histZetaMtmpsin,cmd->sizeHistN);
    }

    for (n = 1; n <= cmd->sizeHistN; n++) {
        for (int l = 1; l <= cmd->sizeHistN; l++) {
            histsincos->histNNSubthreadTC[n][l]
            = hist1sincos->histNNSubthreadTC[n][l];
        }
    }

    return SUCCESS;
}

//E Routines as in cballsutils


local int search_init_omp_balls6_kkk(struct cmdline_data* cmd,
                                 struct  global_data* gd,
                                 gdhistptr_omp_balls6_kkk hist, int ifile)
{
#  define FACTIVE  0.75
#  define FACTOR  316

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

local int search_free_omp_balls6_kkk(struct cmdline_data* cmd,
                                  struct  global_data* gd,
                                  gdhistptr_omp_balls6_kkk hist)
{
    free(hist->interact);
    free(hist->active);

    return SUCCESS;
}

local int print_info(struct cmdline_data* cmd,
                                  struct  global_data* gd)
{
    verb_print(cmd->verbose,
               "searchcalc_normal: Running octree-kkk-balls4... \n");
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
    verb_print(cmd->verbose, "without NMultipoles... \n");
    verb_print(cmd->verbose, "without NONORMHIST... \n");
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

    verb_print(cmd->verbose, "with TRIANGLESTC... \n");

    return SUCCESS;
}

//B Saving histograms section: case KKKCORRELATION:

local int PrintHistrBins(struct  cmdline_data* cmd, struct  global_data* gd)
{
    real rBin, rbinlog;
    int n;
    stream outstr;

    outstr = stropen(gd->fpfnamehistrBinsFileName, "w!");

    verb_print_q(2, cmd->verbose,
               "Printing : to a file %s ...\n",gd->fpfnamehistrBinsFileName);

    for (n=1; n<=cmd->sizeHistN; n++) {
        if (cmd->useLogHist) {
            if (cmd->rminHist==0) {
                rbinlog = ((real)(n-cmd->sizeHistN))/cmd->logHistBinsPD + rlog10(cmd->rangeN);
            } else {
                rbinlog = rlog10(cmd->rminHist) + ((real)(n)-0.5)*gd->deltaR;
            }
            rBin=rpow(10.0,rbinlog);
        } else {
            rBin = cmd->rminHist + ((real)n-0.5)*gd->deltaR;
        }
        fprintf(outstr,"%16.8e\n",rBin);
    }
    fclose(outstr);

    return SUCCESS;
}

// Saves matrix ZetaM for each m multipole
local int PrintHistZetaM_sincos(struct  cmdline_data* cmd,
                                struct  global_data* gd,
                                gdlptr_sincos_omp_kkk gdl)
{
    int n1, n2, m;
    stream outstr;
    char namebuf[256];

    for (m=1; m<=cmd->mChebyshev+1; m++) {
        sprintf(namebuf, "%s%s_%d%s", gd->fpfnamehistZetaMFileName,
                "_cos", m, EXTFILES);
        verb_print_q(2, cmd->verbose, "Printing : to a file %s ...\n",namebuf);
        outstr = stropen(namebuf, "w!");
        for (n1=1; n1<=cmd->sizeHistN; n1++) {
            for (n2=1; n2<=cmd->sizeHistN; n2++) {
                fprintf(outstr,"%16.8e ",gdl->histZetaMcos[m][n1][n2]);
            }
            fprintf(outstr,"\n");
        }
        fclose(outstr);
    }
    for (m=1; m<=cmd->mChebyshev+1; m++) {
        sprintf(namebuf, "%s%s_%d%s", gd->fpfnamehistZetaMFileName,
                "_sin", m, EXTFILES);
        verb_print_q(2, cmd->verbose, "Printing : to a file %s ...\n",namebuf);
        outstr = stropen(namebuf, "w!");
        for (n1=1; n1<=cmd->sizeHistN; n1++) {
            for (n2=1; n2<=cmd->sizeHistN; n2++) {
                fprintf(outstr,"%16.8e ",gdl->histZetaMsin[m][n1][n2]);
            }
            fprintf(outstr,"\n");
        }
        fclose(outstr);
    }

    return SUCCESS;
}

#define MHISTZETA \
"%16.8e %16.8e %16.8e %16.8e %16.8e %16.8e\n"

#define MHISTZETAHEADER \
"# [1] rBins; [2] diagonal; [3] theta2=Nbins/4.0; [4] theta2=2.0*Nbins/4.0; \
[5] theta2=3.0*Nbins/4.0; [6] theta2=4.0*Nbins/4.0 - 1.0\n"


// Saves matrix ZetaM for each m multipole at a set of theta2 angles
local int PrintHistZetaMm_sincos(struct  cmdline_data* cmd,
                                struct  global_data* gd,
                                 gdlptr_sincos_omp_kkk gdl)
{
    real rBin, rbinlog;
    int n1, m;
    stream outstr;
    real Zeta;
    real Zeta2;
    real Zeta3;
    real Zeta4;
    real Zeta5;
    int Nbins;
    char namebuf[256];

    Nbins = cmd->sizeHistN;

    for (m = 1; m <= cmd->mChebyshev+1; m++) {
        sprintf(namebuf, "%s%s_%d%s", gd->fpfnamemhistZetaMFileName,
                "_cos", m, EXTFILES);
        verb_print_q(2, cmd->verbose, "Printing : to a file %s ...\n",namebuf);
        outstr = stropen(namebuf, "w!");
        fprintf(outstr,MHISTZETAHEADER);
        for (n1=1; n1<=cmd->sizeHistN; n1++) {
            if (cmd->useLogHist) {
                if (cmd->rminHist==0) {
                    rbinlog = ((real)(n1-cmd->sizeHistN))/cmd->logHistBinsPD
                    + rlog10(cmd->rangeN);
                } else {
                    rbinlog = rlog10(cmd->rminHist) + ((real)(n1)-0.5)*gd->deltaR;
                }
                rBin=rpow(10.0,rbinlog);
            } else {
                rBin = cmd->rminHist + ((real)n1-0.5)*gd->deltaR;
            }
            Zeta = gdl->histZetaMcos[m][n1][n1];
            Zeta2 = gdl->histZetaMcos[m][n1][(int)(Nbins/4.0)];
            Zeta3 = gdl->histZetaMcos[m][n1][(int)(2.0*Nbins/4.0)];
            Zeta4 = gdl->histZetaMcos[m][n1][(int)(3.0*Nbins/4.0)];
            Zeta5 = gdl->histZetaMcos[m][n1][(int)(4.0*Nbins/4.0 - 1.0)];
            fprintf(outstr,MHISTZETA,rBin,Zeta,Zeta2,Zeta3,Zeta4,Zeta5);
        }
        fclose(outstr);
    }
        
    for (m = 1; m <= cmd->mChebyshev+1; m++) {
        sprintf(namebuf, "%s%s_%d%s", gd->fpfnamemhistZetaMFileName,
                "_sin", m, EXTFILES);
        verb_print_q(2, cmd->verbose, "Printing : to a file %s ...\n",namebuf);
        outstr = stropen(namebuf, "w!");
        fprintf(outstr,MHISTZETAHEADER);
        for (n1=1; n1<=cmd->sizeHistN; n1++) {
            if (cmd->useLogHist) {
                if (cmd->rminHist==0) {
                    rbinlog = ((real)(n1-cmd->sizeHistN))/cmd->logHistBinsPD
                    + rlog10(cmd->rangeN);
                } else {
                    rbinlog = rlog10(cmd->rminHist)
                                + ((real)(n1)-0.5)*gd->deltaR;
                }
                rBin=rpow(10.0,rbinlog);
            } else {
                rBin = cmd->rminHist + ((real)n1-0.5)*gd->deltaR;
            }
                Zeta = gdl->histZetaMsin[m][n1][n1];
                Zeta2 = gdl->histZetaMsin[m][n1][(int)(Nbins/4.0)];
                Zeta3 = gdl->histZetaMsin[m][n1][(int)(2.0*Nbins/4.0)];
                Zeta4 = gdl->histZetaMsin[m][n1][(int)(3.0*Nbins/4.0)];
                Zeta5 = gdl->histZetaMsin[m][n1][(int)(4.0*Nbins/4.0 - 1.0)];
                fprintf(outstr,MHISTZETA,rBin,Zeta,Zeta2,Zeta3,Zeta4,Zeta5);
        }
        fclose(outstr);
    }

    return SUCCESS;
}

#undef MHISTZETAHEADER
#undef MHISTZETA

//E Saving histograms section: case KKKCORRELATION:

local void sumnode_balls6_omp(struct cmdline_data* cmd,
                                        struct  global_data* gd,
                              cellptr start, cellptr finish, bodyptr p0,
                              gdhistptr_sincos_omp_kkk hist)
{
    string routine_name = "sumnode_balls6_omp_triangles";
    cellptr p;
    vector dr;
    real dr1;
    nodeptr pb;
    INTEGER ibodycount=0;
    int n;
    real xi;
    real cosphi;
    real sinphi;
    real theta1;

    for (p = start; p < finish; p++) {
        pb = ((nodeptr) p);
        if (accept_body(cmd, gd, p0, pb, &dr1, dr)) {
            if(dr1>cmd->rminHist) {
                ibodycount++;
                if (cmd->rminHist==0)
                    n = (int)(cmd->logHistBinsPD*(rlog10(dr1)
                            - rlog10(cmd->rangeN)) + cmd->sizeHistN) + 1;
                else
                    n = (int)(rlog10(dr1/cmd->rminHist) * gd->i_deltaR) + 1;
                if (n<=cmd->sizeHistN && n>=1) {
                    xi = Kappa(pb);
                    //B Component of pb with respect to the axis of reference
                    real s, sy;
                    vector pr0;
                    DOTVP(s, dr, hist->dr0);
                    CROSSVP(pr0,hist->dr0,Pos(p));
                    DOTVP(sy, dr, pr0);
                    theta1 = angle_dxdy(s, sy);
                    //E
//B second vertix h
                    cellptr h;
                    nodeptr hb;
                    int l;
                    real dr1_h;
                    vector dr_h;
                    real xih;
                    real theta, theta2;
                    for (h = start; h < finish; h++) {
                        if (h==p) continue;
                        hb = ((nodeptr) h);

                        if (accept_body(cmd, gd, p0, hb, &dr1_h, dr_h)) {
                            if(dr1_h>cmd->rminHist) {
                                ibodycount++;
                                if (cmd->rminHist==0)
                                    l = (int)(cmd->logHistBinsPD*(rlog10(dr1_h)
                                    - rlog10(cmd->rangeN)) + cmd->sizeHistN) + 1;
                                else
                                    l = (int)(rlog10(dr1_h/cmd->rminHist)
                                              * gd->i_deltaR) + 1;
                                if (l<=cmd->sizeHistN && l>=1) {
                                    xih = Kappa(hb);
                //B Component of hb with respect to the axis of reference
                                    DOTVP(s, dr_h, hist->dr0);
                                    CROSSVP(pr0,hist->dr0,Pos(h));
                                    DOTVP(sy, dr_h, pr0);
                                    theta2 = angle_dxdy(s, sy);
                //E
                                    theta = theta2-theta1; // squezed triangle?
                                    cosphi = rcos(theta);
                                    sinphi = rsin(theta);
                                    CHEBYSHEVTUOMPANYMMTC;
                                    for (int m=1; m<=cmd->mChebyshev+1; m++) {
                                        hist->xiOUTVPcosTC[m][n][l]
                                                += hist->histXithreadcosTC[m];
                                        hist->xiOUTVPsinTC[m][n][l]
                                                += hist->histXithreadsinTC[m];
                                    }
                                    hist->histNNSubthreadTC[n][l] =
                                            hist->histNNSubthreadTC[n][l] + 1.;
                                } // ! 1 < l < HistN
                            } // ! dr1_h > rmin
                        } // ! accept_body h
                    } // ! loop h
//E second vertix h
                } // ! 1 < n < HistN
            } // ! dr1 > rmin
        } // ! accept_body
    } // ! loop p

}

local void sumcell_balls6_omp(struct cmdline_data* cmd,
                                        struct  global_data* gd,
                              cellptr start, cellptr finish, bodyptr p0,
                              gdhistptr_sincos_omp_kkk hist)
{
    cellptr p;
    vector dr;
    real dr1;
    nodeptr pb;
    INTEGER ibodycount=0;
    int n;
    real xi;
    real cosphi;
    real sinphi;
    real theta1;

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
                    // needs to multiply xi by Weight(p) of the cell q
                    xi = Kappa(pb);
                    //B Component of pb with respect to the axis of reference
                    real s, sy;
                    vector pr0;
                    DOTVP(s, dr, hist->dr0);
                    CROSSVP(pr0,hist->dr0,Pos(p));
                    DOTVP(sy, dr, pr0);
                    theta1 = angle_dxdy(s, sy);
                    //E
//B second vertix h
                    cellptr h;
                    nodeptr hb;
                    int l;
                    real dr1_h;
                    vector dr_h;
                    real xih;
                    real theta, theta2;
                    for (h = start; h < finish; h++) {
                        if (h==p) continue;
                        hb = ((nodeptr) h);
                        if (accept_body(cmd, gd, p0, hb, &dr1_h, dr_h)) {
                            if(dr1_h>cmd->rminHist) {
                                ibodycount++;
                                if (cmd->rminHist==0)
                                    l = (int)(cmd->logHistBinsPD*(rlog10(dr1)
                                            - rlog10(cmd->rangeN))
                                              + cmd->sizeHistN) + 1;
                                else
                                    l = (int)(rlog10(dr1/cmd->rminHist)
                                              * gd->i_deltaR) + 1;
                                if (l<=cmd->sizeHistN && l>=1) {
                            // needs to multiply xi by Weight(p) of the cell q
                                    xih = Kappa(hb);
//B Component of pb with respect to the axis of reference
                                    DOTVP(s, dr, hist->dr0);
                                    CROSSVP(pr0,hist->dr0,Pos(h));
                                    DOTVP(sy, dr, pr0);
                                    theta2 = angle_dxdy(s, sy);
                                    theta = theta2-theta1; // squezed triangle?
                                    cosphi = rcos(theta);
                                    sinphi = rsin(theta);
//E
                                    CHEBYSHEVTUOMPANYMMTC;
                                    for (int m=1; m<=cmd->mChebyshev+1; m++) {
                                        hist->xiOUTVPcosTC[m][n][l]
                                                += hist->histXithreadcosTC[m];
                                        hist->xiOUTVPsinTC[m][n][l]
                                                += hist->histXithreadsinTC[m];
                                    }
                                    hist->histNNSubthreadTC[n][l] =
                                            hist->histNNSubthreadTC[n][l] + 1.;
                                } // ! 1 < l < HistN
                            } // ! dr1_h > rmin
                        } // ! accept_body
                    } // ! loop h
//E second vertix h
                } // ! 1 < n < HistN
            } // ! dr1 > rmin
        } // ! accept_body
    } // ! loop p
}

local void sumcellcell_balls6_omp(struct cmdline_data* cmd,
                                  struct  global_data* gd,
                                  cellptr start, cellptr finish, nodeptr p0,
                                  gdhistptr_sincos_omp_kkk hist)
{
    cellptr p;
    vector dr;
    real dr1;
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
                    // needs to multiply xi by Weight(p) of the cell q
                    xi = Kappa(pb);
                    //B Component of pb with respect to the axis of reference
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
                    //E
                } // ! 1 < n < HistN
            } // ! dr1 > rmin
        } // ! accept_body
    } // ! loop p

    hist->ncccalcthread += 1;
}

