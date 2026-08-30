/* ==============================================================================
 MODULE: search_octree_ggg_omp_triangles.c		[cTreeBalls]
 Written by: M.A. Rodriguez-Meza
 Starting date:    april 2023
 Purpose: 3-point correlation function computation
 Language: C
 Use: searchcalc_octree_ggg_omp(cmd, gd, btable, nbody,
                                           ipmin, ipmax, cat1, cat2);
 Major revisions:
 ==============================================================================*/
//        1          2          3          4        ^ 5          6          7


// Work to do in order to use with boxes not centered at (0,0,...)

#include "globaldefs.h"

//B Macro for any posible value of mChebyshev
//  for recursivity needs that at least 3 multipoles be evaluated
#define CHEBYSHEVTUOMPANYMMTC                                     \
{real xicosmphi,xisinmphi; int m;                               \
    hist->ChebsT[1] = 1.0;                                      \
    xicosmphi = xi*xih * hist->ChebsT[1];                           \
    hist->histXithreadcosTC[1] = xicosmphi;                   \
    hist->ChebsT[2] = cosphi;                                   \
    xicosmphi = xi*xih * hist->ChebsT[2];                           \
    hist->histXithreadcosTC[2] = xicosmphi;                   \
    hist->ChebsT[3] = 2.0*(cosphi)*(cosphi) - (1.0);            \
    xicosmphi = xi*xih * hist->ChebsT[3];                           \
    hist->histXithreadcosTC[3] = xicosmphi;                   \
    hist->ChebsU[1] = 0.0;                                      \
    xisinmphi = xi*xih * hist->ChebsU[1] * sinphi;                  \
    hist->histXithreadsinTC[1] = xisinmphi;                   \
    hist->ChebsU[2] = 1.0;                                      \
    xisinmphi = xi*xih * hist->ChebsU[2] * sinphi;                  \
    hist->histXithreadsinTC[2] = xisinmphi;                   \
    hist->ChebsU[3] = 2.0*cosphi;                               \
    xisinmphi = xi*xih * hist->ChebsU[3] * sinphi;                  \
    hist->histXithreadsinTC[3] = xisinmphi;                   \
    for (m=4; m<=cmd->mChebyshev+1; m++){                       \
        hist->ChebsT[m] = 2.0*(cosphi)*hist->ChebsT[m-1] - hist->ChebsT[m-2]; \
        xicosmphi = xi*xih * hist->ChebsT[m];                       \
        hist->histXithreadcosTC[m] = xicosmphi;               \
        hist->ChebsU[m] = 2.0*(cosphi)*hist->ChebsU[m-1] - hist->ChebsU[m-2]; \
        xisinmphi = xi*xih * hist->ChebsU[m] * sinphi;              \
        hist->histXithreadsinTC[m] = xisinmphi;               \
    }}
//E

//#if defined(NMultipoles) || defined(NONORMHIST)
//B Macro for any posible value of mChebyshev
//  for recursivity needs that at least 3 multipoles be evaluated
#define CHEBYSHEVTUOMPNANYMMTC                                        \
{real xicosmphi,xisinmphi; int m;                                 \
    histN->ChebsT[1] = 1.0;                                       \
    xicosmphi = xiN*xiNh * histN->ChebsT[1];                           \
    histN->histXithreadcosTC[1] = xicosmphi;                    \
    histN->ChebsT[2] = 1.0;                                       \
    xicosmphi = xiN*xiNh * histN->ChebsT[2];                           \
    histN->histXithreadcosTC[2] = xicosmphi;                    \
    histN->ChebsT[3] = 1.0;                                       \
    xicosmphi = xiN*xiNh * histN->ChebsT[3];                           \
    histN->histXithreadcosTC[3] = xicosmphi;                    \
    histN->ChebsU[1] = 0.0;                                       \
    xisinmphi = xiN*xiNh * histN->ChebsU[1] * sinphi;                  \
    histN->histXithreadsinTC[1] = xisinmphi;                    \
    histN->ChebsU[2] = 0.0;                                       \
    xisinmphi = xiN*xiNh * histN->ChebsU[2] * sinphi;                  \
    histN->histXithreadsinTC[2] = xisinmphi;                    \
    histN->ChebsU[3] = 0.0;                                       \
    xisinmphi = xiN*xiNh * histN->ChebsU[3] * sinphi;                  \
    histN->histXithreadsinTC[3] = xisinmphi;                    \
    for (m=4; m<=cmd->mChebyshev+1; m++){                         \
        histN->ChebsT[m] = 2.0*(cosphi)*histN->ChebsT[m-1]-histN->ChebsT[m-2]; \
        xicosmphi = xiN*xiNh * histN->ChebsT[m];                       \
        histN->histXithreadcosTC[m] = xicosmphi;                \
        histN->ChebsU[m] = 2.0*(cosphi)*histN->ChebsU[m-1]-histN->ChebsU[m-2]; \
        xisinmphi = xiN*xiNh * histN->ChebsU[m] * sinphi;              \
        histN->histXithreadsinTC[m] = xisinmphi;                \
    }}
//E
//#endif // ! NMultipoles


//B Define structures:
typedef struct {
    real **histNNSubTC;
    real ***histZetaMcos;
    real ***histZetaMsin;
} gdl_sincos_omp_ggg, *gdlptr_sincos_omp_ggg;

typedef struct {
    real *ChebsT;
    real *ChebsU;

    real **histZetaMtmpcos;
    real **histZetaMtmpsin;

    real ***histZetaMthreadcos;
    real ***histZetaMthreadsin;

    real *histXithreadcosTC;
    real *histXithreadsinTC;
    real ***xiOUTVPcosTC;
    real ***xiOUTVPsinTC;
    real **histNNSubthreadTC;

    INTEGER nbbcalcthread;
    INTEGER nbccalcthread;

    vector q0;
    real drpq2, drpq;
    vector dr0;
    real cosb;
    real sinb;
} gdhist_sincos_omp_ggg, *gdhistptr_sincos_omp_ggg;

//#if defined(NMultipoles) || defined(NONORMHIST)
typedef struct {
    real **histNNSubTC;
    real ***histZetaMcos;
    real ***histZetaMsin;
} gdl_sincos_omp_ggg_N, *gdlptr_sincos_omp_ggg_N;

// same as gdhist_sincos_omp_ggg above
//  it will necessary to compute only NNN 3pcf
typedef struct {
    real *ChebsT;
    real *ChebsU;

    real **histZetaMtmpcos;
    real **histZetaMtmpsin;

    real ***histZetaMthreadcos;
    real ***histZetaMthreadsin;

    real *histXithreadcosTC;
    real *histXithreadsinTC;
    real ***xiOUTVPcosTC;
    real ***xiOUTVPsinTC;
    real **histNNSubthreadTC;

    INTEGER nbbcalcthread;
    INTEGER nbccalcthread;

    vector q0;
    real drpq2, drpq;
    vector dr0;
    real cosb;
    real sinb;

} gdhist_sincos_omp_ggg_N, *gdhistptr_sincos_omp_ggg_N;
//#endif
//E

local int search_init_gd_sincos_omp_ggg(struct  cmdline_data* cmd,
                                        struct  global_data* gd,
                                        gdlptr_sincos_omp_ggg);
local int search_free_gd_sincos_omp_ggg(struct  cmdline_data* cmd,
                                         struct  global_data* gd,
                                        gdlptr_sincos_omp_ggg);
local int search_init_sincos_omp_ggg(struct  cmdline_data* cmd,
                                  struct  global_data* gd,
                                     gdhistptr_sincos_omp_ggg hist);
local int search_free_sincos_omp_ggg(struct  cmdline_data* cmd,
                                         struct  global_data* gd,
                                     gdhistptr_sincos_omp_ggg hist);
local int computeBodyProperties_sincos_ggg(struct  cmdline_data* cmd,
                                            struct  global_data* gd,
                                            bodyptr p, int nbody,
                                       gdhistptr_sincos_omp_ggg hist);

local int print_info(struct cmdline_data* cmd,
                     struct  global_data* gd);

//#if defined(NMultipoles) || defined(NONORMHIST)
local int search_init_gd_sincos_omp_ggg_N(struct  cmdline_data* cmd,
                                          struct  global_data* gd,
                                          gdlptr_sincos_omp_ggg_N);
local int search_free_gd_sincos_omp_ggg_N(struct  cmdline_data* cmd,
                                         struct  global_data* gd,
                                          gdlptr_sincos_omp_ggg_N);
local int search_init_sincos_omp_ggg_N(struct  cmdline_data* cmd,
                                  struct  global_data* gd,
                                     gdhistptr_sincos_omp_ggg_N hist);
local int search_free_sincos_omp_ggg_N(struct  cmdline_data* cmd,
                                         struct  global_data* gd,
                                     gdhistptr_sincos_omp_ggg_N);
local int computeBodyProperties_sincos_ggg_N(struct  cmdline_data* cmd,
                                            struct  global_data* gd,
                                            bodyptr p, int nbody,
                                       gdhistptr_sincos_omp_ggg_N);
//#endif


//B Saving histograms section: case KKKCORRELATION:
local int PrintHistrBins(struct  cmdline_data* cmd, struct  global_data* gd);
local int PrintHistZetaM_sincos(struct  cmdline_data* cmd,
                                struct  global_data* gd,
                                gdlptr_sincos_omp_ggg);
local int PrintHistZetaMm_sincos(struct  cmdline_data* cmd,
                               struct  global_data* gd,
                                 gdlptr_sincos_omp_ggg);

// NMultipoles precedes NONORMHIST?... Not necessarily.
//  Check consistency!!!
//  (NMultipoles, NONORMHIST):
//      1               1
//      1               0
//      0               1
//      0               0
//#if defined(NMultipoles) || defined(NONORMHIST)
local int PrintHistZetaM_sincos_N(struct  cmdline_data* cmd,
                                struct  global_data* gd,
                                  gdlptr_sincos_omp_ggg_N);
local int PrintHistZetaMm_sincos_N(struct  cmdline_data* cmd,
                               struct  global_data* gd,
                                   gdlptr_sincos_omp_ggg_N);
//#ifdef NONORMHIST
// Saves matrix ZetaM for each m multipole
local int PrintHistZetaM_sincos_normalized(struct  cmdline_data* cmd,
                                           struct  global_data* gd,
                                           gdlptr_sincos_omp_ggg,
                                           gdlptr_sincos_omp_ggg_N);
// Saves matrix ZetaM for each m multipole at a set of theta2 angles
local int PrintHistZetaMm_sincos_normalized(struct  cmdline_data* cmd,
                                            struct  global_data* gd,
                                            gdlptr_sincos_omp_ggg,
                                            gdlptr_sincos_omp_ggg_N);
//B edge effects:
local int PrintHistZetaM_sincos_edge_effects(struct  cmdline_data*,
                                             struct  global_data*,
                                             gdlptr_sincos_omp_ggg,
                                             gdlptr_sincos_omp_ggg_N);

//#endif
//#endif // ! NMultipoles
//E

//B kappa Avg Rmin
#ifdef DEBUG
local char pivotsfilePath[MAXLENGTHOFFILES];
local FILE *outpivots;
#endif

#define FACTIVE  0.75

local int actlen;
local int *activenb;

local void normal_walktree_nblist_omp(struct cmdline_data*,
                                      struct  global_data*,
                                      bodyptr *, int,
                                      bodyptr, nodeptr, real,
                                      gdhistptr_sincos_omp_ggg, int *);

local void find_nblist_omp(struct cmdline_data* cmd,
                           struct  global_data* gd,
                           bodyptr *btable, int cat2,
                           bodyptr p, cellptr start, cellptr finish, int *);

local void sumnode_nblist_omp_triangles(struct cmdline_data*,
                              struct  global_data*,
                              bodyptr *,
                              INTEGER, INTEGER *,
                              int, int,
                              bodyptr,
                              gdhistptr_sincos_omp_ggg, int);

//#if defined(NMultipoles) || defined(NONORMHIST)
local void sumnode_nblist_omp_triangles_N(struct cmdline_data*,
                                          struct  global_data*,
                                          bodyptr *,
                                          INTEGER, INTEGER *,
                                          int, int,
                                          bodyptr,
                                          gdhistptr_sincos_omp_ggg,
                                          gdhistptr_sincos_omp_ggg_N,
                                          int);
//#endif


/*
 Search routine using octree method:

 To be called using: search=octree-ggg-omp-triangles

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
global int searchcalc_octree_ggg_omp_triangles(struct cmdline_data* cmd,
                                     struct  global_data* gd,
                                     bodyptr *btable, INTEGER *nbody,
                                     INTEGER ipmin, INTEGER *ipmax,
                                     int cat1, int cat2)
{
    string routine_name = "searchcalc_octree_ggg_omp_triangles";
    double cpustart;
    gdl_sincos_omp_ggg gdl;
//#if defined(NMultipoles) || defined(NONORMHIST)
    gdl_sincos_omp_ggg_N gdlN;
//#endif

    cpustart = CPUTIME;
    print_info(cmd, gd);

//B kappa Avg Rmin
#ifdef DEBUG
    sprintf(pivotsfilePath,"%s/pivot_info%s.txt",gd->tmpDir,cmd->suffixOutFiles);
    if(!(outpivots=fopen(pivotsfilePath, "w")))
        error("\n%s: error opening file '%s' \n",
              routine_name, pivotsfilePath);
#endif
//E

#ifdef OPENMPCODE
    ThreadCount(cmd, gd, nbody[cat1], cat1);
#else
#error `OPENMPMACHINE` is not defined. Switch it on in Makefile_settings
#endif

    search_init_gd_sincos_omp_ggg(cmd, gd, &gdl);
//#if defined(NMultipoles) || defined(NONORMHIST)
    search_init_gd_sincos_omp_ggg_N(cmd, gd, &gdlN);
//#endif

    //B kappa Avg Rmin
    INTEGER ipfalse;
    ipfalse=0;
    INTEGER icountNbRmin;
    icountNbRmin=0;
    INTEGER icountNbRminOverlap;
    icountNbRminOverlap=0;

    //B Alloc memory for neighbour lists
    actlen = FACTIVE * 216 * 512 * gd->tdepthTable[cat2];
    verb_log_print(cmd->verbose,gd->outlog,"%s: actlen = %ld",
                   routine_name, actlen);
    activenb = (int *) allocate(actlen * sizeof(int)); // change int->INTEGER??
    //E

    verb_print(cmd->verbose,
        "\n%s: Total allocated %g MByte storage so far.\n",
               routine_name, gd->bytes_tot/(1024.0*1024.0));

    if (cmd->verbose >= VERBOSENORMALINFO)
        verb_print(cmd->verbose,
                   "\nRunning...\n - Completed pivot node:\n");

#ifdef DEBUG
#pragma omp parallel default(none)                                          \
    shared(cmd,gd,btable,nbody,roottable,outpivots,                         \
           ipmin,ipmax,cat1,cat2,ipfalse,icountNbRmin,icountNbRminOverlap, gdl)
#else // ! DEBUG
//#if defined(NMultipoles) || defined(NONORMHIST)
#pragma omp parallel default(none)                                          \
    shared(cmd,gd,btable,nbody,roottable, actlen, activenb,                    \
           ipmin,ipmax,cat1,cat2,ipfalse,icountNbRmin,icountNbRminOverlap,  \
           gdl, gdlN)
//#else // ! NMultipoles
//#pragma omp parallel default(none)                                          \
//    shared(cmd,gd,btable,nbody,roottable, actlen, activenb,                    \
//           ipmin,ipmax,cat1,cat2,ipfalse,icountNbRmin,icountNbRminOverlap,gdl)
//#endif // ! NMultipoles
#endif // ! DEBUG
  {
    bodyptr p;
    bodyptr q;
    int n, m, ip;
    int i;

    //B init:
    gdhist_sincos_omp_ggg hist;
    search_init_sincos_omp_ggg(cmd, gd, &hist);
//#if defined(NMultipoles) || defined(NONORMHIST)
    gdhist_sincos_omp_ggg_N histN;
    search_init_sincos_omp_ggg_N(cmd, gd, &histN);
//#endif
    //E

    //B kappa Avg Rmin
    INTEGER ipfalsethreads;
    ipfalsethreads = 0;
    INTEGER icountNbRminthread;
    icountNbRminthread=0;
    INTEGER icountNbRminOverlapthread;
    icountNbRminOverlapthread=0;
      
    int nblist;

#pragma omp for nowait schedule(dynamic)
      for (p = btable[cat1] + ipmin -1; p < btable[cat1] + ipmax[cat1]; p++) {
          //B kappa Avg Rmin
          NbRmin(p) = 1;
          NbRminOverlap(p) = 0;
          KappaRmin(p) = Kappa(p);
          
          nblist=0;

          if (cballs_opt_smooth_pivot(cmd)) {
              if (Update(p) == FALSE) {
                  ipfalsethreads++;
                  continue;
              }
          }
          //E

          //B Set histograms to zero for the pivot
          CLRM_ext(hist.histNNSubthreadTC, cmd->sizeHistN);
          CLRV_ext(hist.histXithreadcosTC, cmd->mChebyshev+1);
          CLRV_ext(hist.histXithreadsinTC, cmd->mChebyshev+1);
          for (m=1; m<=cmd->mChebyshev+1; m++) {
              CLRM_ext(hist.xiOUTVPcosTC[m], cmd->sizeHistN);
              CLRM_ext(hist.xiOUTVPsinTC[m], cmd->sizeHistN);
          }
//#ifdef NMultipoles
          CLRM_ext(histN.histNNSubthreadTC, cmd->sizeHistN);
          CLRV_ext(histN.histXithreadcosTC, cmd->mChebyshev+1);
          CLRV_ext(histN.histXithreadsinTC, cmd->mChebyshev+1);
          for (m=1; m<=cmd->mChebyshev+1; m++) {
              CLRM_ext(histN.xiOUTVPcosTC[m], cmd->sizeHistN);
              CLRM_ext(histN.xiOUTVPsinTC[m], cmd->sizeHistN);
          }
//#endif
          //E
          //B Set reference axis...
          dRotation3D(Pos(p), ROTANGLE, ROTANGLE, ROTANGLE, hist.q0);
          DOTPSUBV(hist.drpq2, hist.dr0, Pos(p), hist.q0);
          hist.drpq = rsqrt(hist.drpq2);
          //E

//#if defined(NMultipoles) || defined(NONORMHIST)
          normal_walktree_nblist_omp(cmd, gd, btable, cat2,
                                     p, ((nodeptr) roottable[cat2]),
                                     gd->rSizeTable[cat2],
                                     &hist, &nblist);
          int_piksrt(nblist, activenb);
          if (cmd->verbose_log>=3)
              verb_log_print(cmd->verbose_log, gd->outlog,
                       " - Summing nblist: %ld\n", nblist);
          sumnode_nblist_omp_triangles_N(cmd, gd,
                             btable, ipmin, ipmax, cat1, cat2,
                             p, &hist, &histN, nblist);
          computeBodyProperties_sincos_ggg(cmd, gd, p, ipmax[cat1]-ipmin+1, &hist);
          computeBodyProperties_sincos_ggg_N(cmd, gd, p, ipmax[cat1]-ipmin+1,
                                           &histN);
/*#else // ! NMultipoles
          normal_walktree_nblist_omp(cmd, gd, btable, cat2,
                                     p, ((nodeptr) roottable[cat2]),
                                     gd->rSizeTable[cat2],
                                     &hist, &nblist);
          int_piksrt(nblist, activenb);
          if (cmd->verbose_log>=3)
              verb_log_print(cmd->verbose_log, gd->outlog,
                       " - Summing nblist: %ld\n", nblist);

          sumnode_nblist_omp_triangles(cmd, gd,
                             btable, ipmin, ipmax, cat1, cat2,
                             p, &hist, nblist);
          computeBodyProperties_sincos_ggg(cmd, gd, p, ipmax[cat1]-ipmin+1, &hist);
#endif // ! NMultipoles
          */

          ip = p - btable[cat1] + 1;
          //B kappa Avg Rmin
          icountNbRminthread += NbRmin(p);
          icountNbRminOverlapthread += NbRminOverlap(p);
#ifdef DEBUG
          fprintf(outpivots,"%ld \t%ld \t%ld \t\t%g\n",
                  ip, NbRmin(p), NbRminOverlap(p),
                  KappaRmin(p)/NbRmin(p));
#endif
          //E
          if (cmd->verbose >= VERBOSENORMALINFO) {
              if (ip%cmd->stepState == 0) {
                  verb_print(cmd->verbose, "%d\n", ip);
              }
          } else
              if (ip%cmd->stepState == 0) {
                  verb_log_print(cmd->verbose_log, gd->outlog,
                                 " - Completed pivot: %ld\n", ip);
              }
      } // end do body p

#pragma omp critical
    {
//#ifndef NONORMHIST
        for (n=1; n<=cmd->sizeHistN; n++) {
            for (int l=1; l<=cmd->sizeHistN; l++) {
                gdl.histNNSubTC[n][l] = hist.histNNSubthreadTC[n][l];
            }
        }
//#endif // ! NNORMHIST

        for (m=1; m<=cmd->mChebyshev+1; m++) {
            ADDM_ext(gdl.histZetaMcos[m],gdl.histZetaMcos[m],
                     hist.histZetaMthreadcos[m],cmd->sizeHistN);
            ADDM_ext(gdl.histZetaMsin[m],gdl.histZetaMsin[m],
                     hist.histZetaMthreadsin[m],cmd->sizeHistN);
        }
/*
#ifndef NONORMHIST
        for (m=1; m<=cmd->mChebyshev+1; m++) {
            for (n=1; n<=cmd->sizeHistN; n++) {
                for (int l=1; l<=cmd->sizeHistN; l++) {
                    gdl.histZetaMcos[m][n][l]
                            /= MAX(gdl.histNNSubTC[n][l],1.0);
                    gdl.histZetaMcos[m][n][l]
                            /= MAX(gdl.histNNSubTC[n][l],1.0);
                }
            }
        }
#endif // ! NNORMHIST
*/
        gd->nbbcalc += hist.nbbcalcthread;
        gd->nbccalc += hist.nbccalcthread;

//#if defined(NMultipoles) || defined(NONORMHIST)
        for (n=1; n<=cmd->sizeHistN; n++) {
            for (int l=1; l<=cmd->sizeHistN; l++) {
                gdlN.histNNSubTC[n][l] = histN.histNNSubthreadTC[n][l];
            }
        }
        for (m=1; m<=cmd->mChebyshev+1; m++) {
            ADDM_ext(gdlN.histZetaMcos[m],gdlN.histZetaMcos[m],
                     histN.histZetaMthreadcos[m],cmd->sizeHistN);
            ADDM_ext(gdlN.histZetaMsin[m],gdlN.histZetaMsin[m],
                     histN.histZetaMthreadsin[m],cmd->sizeHistN);
        }
//#endif
        //B kappa Avg Rmin
        ipfalse += ipfalsethreads;
        icountNbRmin += icountNbRminthread;
        icountNbRminOverlap += icountNbRminOverlapthread;
        //E
    } // ! critical

//#if defined(NMultipoles) || defined(NONORMHIST)
    search_free_sincos_omp_ggg_N(cmd, gd, &histN);  // free memory
//#endif
    search_free_sincos_omp_ggg(cmd, gd, &hist);     // free memory
  } // end pragma omp parallel

    if (cmd->verbose >= VERBOSENORMALINFO)
        verb_print(cmd->verbose, "\n\n");           // end of completed pivot

    //B kappa Avg Rmin
    real xi, den, num;
    int mm;
    if (cballs_opt_smooth_pivot(cmd)) {
        num = (real)nbody[cat1];
        den = (real)(nbody[cat1]-ipfalse);
//#ifdef NONORMHIST
        xi = 1.0;
//#else
//        xi = num/den;
//#endif // ! NONORMHIST
        if (cmd->verbose>=VERBOSENORMALINFO)
            verb_print(cmd->verbose,
                       "%: p falses found = %ld and %e %e %e\n",
                       routine_name, ipfalse, num, den, xi);

        for (mm=1; mm<=cmd->mChebyshev+1; mm++) {
            MULMS_ext(gdl.histZetaMcos[mm],
                      gdl.histZetaMcos[mm],xi,cmd->sizeHistN);
            MULMS_ext(gdl.histZetaMsin[mm],
                      gdl.histZetaMsin[mm],xi,cmd->sizeHistN);
        }

//#if defined(NMultipoles) || defined(NONORMHIST)
        for (mm=1; mm<=cmd->mChebyshev+1; mm++) {
            MULMS_ext(gdlN.histZetaMcos[mm],
                      gdlN.histZetaMcos[mm],xi,cmd->sizeHistN);
            MULMS_ext(gdlN.histZetaMsin[mm],
                      gdlN.histZetaMsin[mm],xi,cmd->sizeHistN);
        }
//#endif
//E
        if (cmd->verbose>=VERBOSENORMALINFO) {
            verb_print(cmd->verbose,
                       "%s: p falses found = %ld\n", routine_name, ipfalse);
            //B kappa Avg Rmin
            verb_print(cmd->verbose,
                       "%s: count NbRmin found = %ld\n",
                       routine_name, icountNbRmin);
            verb_print(cmd->verbose,
                       "%s: count overlap found = %ld\n",
                       routine_name, icountNbRminOverlap);
        }
        
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
        if (cmd->verbose>=VERBOSENORMALINFO) {
            verb_print(cmd->verbose, "%s: p falses found = %ld\n",
                       routine_name, ifalsecount);
            verb_print(cmd->verbose, "%s: p true found = %ld\n",
                       routine_name, itruecount);
            verb_print(cmd->verbose, "%s: total = %ld\n",
                       routine_name, itruecount+ifalsecount);
        }
        //E
    } // ! smooth-pivot

#ifdef DEBUG
    fclose(outpivots);                              // Close file to debug pivots
#endif

// ===============================================
//B Saving histograms section: case GGGCORRELATION:
// ===============================================
    verb_print(cmd->verbose,
        "\n\t%s: printing octree-ggg-omp method...\n\n", routine_name);

    PrintHistrBins(cmd, gd);

//#if defined(NMultipoles)

//#ifdef NONORMHIST
//#ifndef NMultipoles
//    PrintHistZetaM_sincos(cmd, gd, &gdl);
//#else
        if (cballs_opt_no_normalize_histzeta(cmd))
            PrintHistZetaM_sincos(cmd, gd, &gdl);
        else
            PrintHistZetaM_sincos_normalized(cmd, gd, &gdl, &gdlN);
//#endif
//#else
//        PrintHistZetaM_sincos(cmd, gd, &gdl);
//#endif // ! NONORMHIST

        PrintHistZetaM_sincos_N(cmd, gd, &gdlN);

//#else // ! NMultipoles
//    PrintHistZetaM_sincos(cmd, gd, &gdl);
//#endif // ! NMultipoles

        if (cballs_opt_out_m_histzeta(cmd)) {
//#if defined(NMultipoles)

//#ifdef NONORMHIST
//#ifndef NMultipoles
//            PrintHistZetaMm_sincos(cmd, gd, &gdl);
//#else
            if (cballs_opt_no_normalize_histzeta(cmd))
                PrintHistZetaMm_sincos(cmd, gd, &gdl);
            else
                PrintHistZetaMm_sincos_normalized(cmd, gd, &gdl, &gdlN);
//#endif
//#else
//            PrintHistZetaMm_sincos(cmd, gd, &gdl);
//#endif // ! NONORMHIST

            PrintHistZetaMm_sincos_N(cmd, gd, &gdlN);

//#else // ! NMultipoles
//            PrintHistZetaMm_sincos(cmd, gd, &gdl);
//#endif // ! NMultipoles
        }

//#ifdef NMultipoles
//#ifdef NONORMHIST
        if (cballs_opt_no_normalize_histzeta(cmd)) {
            if (cballs_opt_edge_corrections(cmd)) {
                PrintHistZetaM_sincos_edge_effects(cmd, gd, &gdl, &gdlN);
            }
        }
//#endif // ! NONORMHIST
//#endif // ! NMultipoles

    gd->flagPrint = FALSE;
// ===============================================
//E Saving histograms section: case GGGCORRELATION
// ===============================================

//#if defined(NMultipoles) || defined(NONORMHIST)
    search_free_gd_sincos_omp_ggg_N(cmd, gd, &gdlN);// free memory
//#endif
    search_free_gd_sincos_omp_ggg(cmd, gd, &gdl); // free memory

    gd->cpusearch = CPUTIME - cpustart;
    verb_print(cmd->verbose, "\nGoing out: CPU time = %lf %s\n",
               CPUTIME-cpustart, PRNUNITOFTIMEUSED);

    return SUCCESS;
}

#define ASYMMETRICCASE

local void normal_walktree_nblist_omp(struct cmdline_data* cmd,
                                      struct  global_data* gd,
                                      bodyptr *btable, int cat2,
                                      bodyptr p, nodeptr q, real qsize,
                                      gdhistptr_sincos_omp_ggg hist,
                                      int *nblist)
{
    nodeptr l;

    if (Update(p)==FALSE) return;
    if ( ((nodeptr) p) != q ) {
        if (Type(q) == CELL) {
            if (!reject_cell(cmd, gd, (nodeptr)p, q, qsize)) {
                for (l = More(q); l != Next(q); l = Next(l)) {
                    normal_walktree_nblist_omp(cmd, gd,
                                            btable, cat2,
                                            p,l,qsize/2, hist,
                                            nblist);
                }
            }
        } else
            find_nblist_omp(cmd, gd,
                            btable, cat2,
                            p, ((cellptr) q),( (cellptr) q+1), nblist);
    }
}

local void find_nblist_omp(struct cmdline_data* cmd,
                           struct  global_data* gd,
                           bodyptr *btable, int cat2,
                           bodyptr p, cellptr start, cellptr finish, int *nblist)
{
    cellptr q;
    real dr1;
    vector dr;
    INTEGER iq;

    for (q = start; q < finish; q++) {
        if (accept_body(cmd, gd, p, (nodeptr)q, &dr1, dr)) {
            gd->nbbcalc += 1;
            iq = (bodyptr)q-btable[cat2];
            activenb[*nblist]=iq;
            *nblist +=1;
            if (*nblist > actlen)
                error("find_nblist: too many neighbors\n");
        }
    }

}

/*
local void sumnode_nblist_omp_triangles(struct cmdline_data* cmd,
                                        struct  global_data* gd,
                                        bodyptr *btable,
                                        INTEGER ipmin, INTEGER *ipmax,
                                        int cat1, int cat2,
                                        bodyptr p,
                                        gdhistptr_sincos_omp_ggg hist, int nblist)
{
    string routine_name = "sumnode_nblist_omp_triangles";
    bodyptr q;
    real dr1;
    vector dr;
    int i;
    int n;
    real theta1;
    real xi;

    REAL s, sy;
    REAL cosphi,sinphi;
    vector pr0;

    bodyptr h;
    int m, j, l;
    real dr1_h;
    vector dr_h;
    real theta, theta2;

#ifdef ASYMMETRICCASE
    for (i = 0; i < nblist-1; i++) {
#else
    for (i = 0; i < nblist; i++) {
#endif
//B first vertix q
        q = btable[cat2] + activenb[i];
        accept_body(cmd, gd, p, (nodeptr)q, &dr1, dr);

        if(dr1>cmd->rminHist) {
            if (cmd->rminHist==0)
                n = (int)(cmd->logHistBinsPD*(rlog10(dr1)
                                - rlog10(cmd->rangeN)) + cmd->sizeHistN) + 1;
            else
                n = (int)(rlog10(dr1/cmd->rminHist) * gd->i_deltaR) + 1;
            if (n<=cmd->sizeHistN && n>=1) {
                xi = Kappa(q);

                DOTVP(s, dr, hist->dr0);
                CROSSVP(pr0,hist->dr0,Pos(p));
                DOTVP(sy, dr, pr0);
                theta1 = angle_dxdy(s, sy);
//E first vertix q

#ifdef ASYMMETRICCASE
                for (j = i+1; j < nblist; j++) {
#else
                for (j = 0; j < nblist; j++) {
                    if (j==i) continue;
#endif
//B second vertix h
                    h = btable[cat2] + activenb[j];
                    accept_body(cmd, gd, p, (nodeptr)h, &dr1_h, dr_h);
                    if(dr1_h>cmd->rminHist) {
                        if (cmd->rminHist==0)
                            l = (int)(cmd->logHistBinsPD*(rlog10(dr1_h)
                            - rlog10(cmd->rangeN)) + cmd->sizeHistN) + 1;
                        else
                            l = (int)(rlog10(dr1_h/cmd->rminHist)
                                      * gd->i_deltaR) + 1;
                        if (l<=cmd->sizeHistN && l>=1) {
                            DOTVP(s, dr_h, hist->dr0);
                            CROSSVP(pr0,hist->dr0,Pos(p));
                            DOTVP(sy, dr_h, pr0);
                            theta2 = angle_dxdy(s, sy);
                            theta = theta2-theta1; // squezed triangle?

                            cosphi = rcos(theta);
                            sinphi = rsin(theta);
                            real xih = Kappa(h);
                            CHEBYSHEVTUOMPANYMMTC;
                            for (m=1; m<=cmd->mChebyshev+1; m++) {
                                hist->xiOUTVPcosTC[m][n][l]
                                    += hist->histXithreadcosTC[m];
                                hist->xiOUTVPsinTC[m][n][l]
                                    += hist->histXithreadsinTC[m];
                            }
                            hist->histNNSubthreadTC[n][l] =
                                hist->histNNSubthreadTC[n][l] + Nb(q)*Nb(h);
                        } // ! l <= sizeHistN && l >= 1
                    } // ! dr1_h > rminHist
                } // ! end loop j
//E second vertix h
            } // ! n <= sizeHistN && n >= 1
        } // ! dr1 > rminHist (q)
    } // ! end loop i
}
*/

//#if defined(NMultipoles) || defined(NONORMHIST)

local void sumnode_nblist_omp_triangles_N(struct cmdline_data* cmd,
                                          struct  global_data* gd,
                                          bodyptr *btable,
                                          INTEGER ipmin, INTEGER *ipmax,
                                          int cat1, int cat2,
                                          bodyptr p,
                                          gdhistptr_sincos_omp_ggg hist,
                                          gdhistptr_sincos_omp_ggg_N histN,
                                          int nblist)
{
    string routine_name = "sumnode_nblist_omp_triangles_N";
    bodyptr q;
    real dr1;
    vector dr;
    int i;
    int n;
    real theta1;
    real xi;
    real xiN;

    REAL s, sy;
    REAL cosphi,sinphi;
    vector pr0;

    bodyptr h;
    int m, j, l;
    real dr1_h;
    vector dr_h;
    real theta, theta2;

//B first vertix q
#ifdef ASYMMETRICCASE
    for (i = 0; i < nblist-1; i++) {
#else
    for (i = 0; i < nblist; i++) {
#endif
        q = btable[cat2] + activenb[i];
        accept_body(cmd, gd, p, (nodeptr)q, &dr1, dr);
        if(dr1>cmd->rminHist) {
            if (cmd->rminHist==0)
                n = (int)(cmd->logHistBinsPD*(rlog10(dr1)
                    - rlog10(cmd->rangeN)) + cmd->sizeHistN) + 1;
            else
                n = (int)(rlog10(dr1/cmd->rminHist) * gd->i_deltaR) + 1;
            if (n<=cmd->sizeHistN && n>=1) {
#ifndef NOWKAvg
                    xi = Weight(q)*Kappa(q);
                    xiN = Weight(q);
#else
                    xi = Nb(q)*Kappa(q);
                    xiN = Nb(q)*Weight(q);
#endif

                DOTVP(s, dr, hist->dr0);
                CROSSVP(pr0,hist->dr0,Pos(p));
                DOTVP(sy, dr, pr0);
                theta1 = angle_dxdy(s, sy);
//E first vertix q

//B second vertix h
#ifdef ASYMMETRICCASE
                for (j = i+1; j < nblist; j++) {
#else
                for (j = 0; j < nblist; j++) {
                    if (j==i) continue;
#endif
                    h = btable[cat2] + activenb[j];
                    accept_body(cmd, gd, p, (nodeptr)h, &dr1_h, dr_h);
                    if(dr1_h>cmd->rminHist) {
                        if (cmd->rminHist==0)
                            l = (int)(cmd->logHistBinsPD*(rlog10(dr1_h)
                                    - rlog10(cmd->rangeN)) + cmd->sizeHistN) + 1;
                        else
                            l = (int)(rlog10(dr1_h/cmd->rminHist)*gd->i_deltaR)+ 1;
                        if (l<=cmd->sizeHistN && l>=1) {
                            DOTVP(s, dr_h, hist->dr0);
                            CROSSVP(pr0,hist->dr0,Pos(p));
                            DOTVP(sy, dr_h, pr0);
                            theta2 = angle_dxdy(s, sy);
                            theta = theta2-theta1; // squezed triangle?
                            cosphi = rcos(theta);
                            sinphi = rsin(theta);
#ifndef NOWKAvg
                            real xih = Weight(h)*Kappa(h);
                            real xiNh = Weight(h);
#else
                            real xih = Nb(h)*Kappa(h);
                            real xiNh = Nb(h)*Weight(h);
#endif
                            CHEBYSHEVTUOMPNANYMMTC;
                            CHEBYSHEVTUOMPANYMMTC;
                            for (m=1; m<=cmd->mChebyshev+1; m++) {
                                hist->xiOUTVPcosTC[m][n][l]
                                        += hist->histXithreadcosTC[m];
                                hist->xiOUTVPsinTC[m][n][l]
                                        += hist->histXithreadsinTC[m];
                                histN->xiOUTVPcosTC[m][n][l]
                                        += histN->histXithreadcosTC[m];
                                histN->xiOUTVPsinTC[m][n][l]
                                        += histN->histXithreadsinTC[m];
                            }
                            hist->histNNSubthreadTC[n][l] =
                                    hist->histNNSubthreadTC[n][l] + Nb(q)*Nb(h);
                            histN->histNNSubthreadTC[n][l] =
                                    histN->histNNSubthreadTC[n][l] + Nb(q)*Nb(h);
                        } // ! l <= sizeHistN && l >= 1
                    } // ! dr1_h > rminHist
                } // ! end loop j
//E second vertix h
            } // ! n <= sizeHistN && n >= 1
        } // ! dr1 > rminHist (q)
    } // ! end loop i
}

//#endif // ! NMultipoles

#undef ASYMMETRICCASE

local int computeBodyProperties_sincos_ggg(struct  cmdline_data* cmd,
                                            struct  global_data* gd,
                                            bodyptr p, int nbody,
                                            gdhistptr_sincos_omp_ggg hist)
{
    int n;
    int m;
    real xi;

// check Weight factor... must be an average of Weights
//#ifdef NONORMHIST
    xi = Weight(p)*Kappa(p);
    if (cballs_opt_smooth_pivot(cmd)) {
        xi = KappaRmin(p)/NbRmin(p);
    }
/*#else // ! NONORMHIST
    xi = Weight(p)*Kappa(p)/nbody;
//    xi = Weight(p)*Kappa(p);
    //B kappa Avg Rmin
    if (cballs_opt_smooth_pivot(cmd)) {
        xi = (KappaRmin(p)/NbRmin(p))/nbody;
//        xi = (KappaRmin(p)/NbRmin(p));
    }
    //E
#endif // ! NONORMHIST
*/
/*
#ifndef NONORMHIST
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
#endif // ! NNORMHIST
*/
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

//#if defined(NMultipoles) || defined(NONORMHIST)

local int computeBodyProperties_sincos_ggg_N(struct  cmdline_data* cmd,
                                            struct  global_data* gd,
                                            bodyptr p, int nbody,
                                            gdhistptr_sincos_omp_ggg_N hist)
{
    int n;
    int m;
    real xi;
    
//#ifdef NONORMHIST
//    xi = 1.0;
    xi = Weight(p);
    if (cballs_opt_smooth_pivot(cmd)) {
//        xi = KappaRmin(p)/NbRmin(p);
        xi = Weight(p);                 // check this!!!
    }
/*#else
//    xi = 1.0/nbody;
    xi = Weight(p)/nbody;
    //B kappa Avg Rmin
    if (cballs_opt_smooth_pivot(cmd)) {
//        xi = (KappaRmin(p)/NbRmin(p))/nbody;
        xi = Weight(p)/nbody;  // check this!!!
    }
    //E
#endif */

    /*
    for (m=1; m<=cmd->mChebyshev+1; m++) {
#ifndef NONORMHIST
        for (n=1; n<=cmd->sizeHistN; n++) {
            for (int l=1; l<=cmd->sizeHistN; l++) {
                hist->xiOUTVPcosTC[m][n][l]
                        /= MAX(hist->histNNSubthreadTC[n][l],1.0);
                hist->xiOUTVPsinTC[m][n][l]
                        /= MAX(hist->histNNSubthreadTC[n][l],1.0);
            }
        }
#endif // ! NONORMHIST
    }
*/

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

//#endif // ! NMultipoles


//B Routines as in cballsutils

local int search_init_gd_sincos_omp_ggg(struct  cmdline_data* cmd,
                                        struct  global_data* gd,
                                        gdlptr_sincos_omp_ggg gdl)
{
    int n;
    int m;
    INTEGER bytes_tot_local=0;

//    real **histNNSubTC;
//    real ***histZetaMcos;
//    real ***histZetaMsin;

//    gdl->histNNSub = dvector(1,cmd->sizeHistN);
    bytes_tot_local += 1*cmd->sizeHistN*sizeof(real);

    gdl->histNNSubTC = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);

    gdl->histZetaMcos =
            dmatrix3D(1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    gdl->histZetaMsin =
            dmatrix3D(1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
//    gdl->histZetaM =
//            dmatrix3D(1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    bytes_tot_local +=
            5*(cmd->mChebyshev+1)*cmd->sizeHistN*cmd->sizeHistN*sizeof(real);


    gd->bytes_tot += bytes_tot_local;
    verb_print(cmd->verbose,
    "\nsearch_init_gd_octree_ggg: Allocated %g MByte for histograms storage.\n",
    bytes_tot_local*INMB);

//    for (n = 1; n <= cmd->sizeHistN; n++)
//        gdl->histNNSub[n] = 0.0;

    CLRM_ext(gdl->histNNSubTC, cmd->sizeHistN);

    for (m = 1; m <= cmd->mChebyshev+1; m++) {
        CLRM_ext(gdl->histZetaMcos[m], cmd->sizeHistN);
        CLRM_ext(gdl->histZetaMsin[m], cmd->sizeHistN);
//        CLRM_ext(gdl->histZetaM[m], cmd->sizeHistN);
    }

    gd->nbbcalc = gd->nbccalc = gd->ncccalc = 0;

    return SUCCESS;
}

local int search_free_gd_sincos_omp_ggg(struct  cmdline_data* cmd,
                                         struct  global_data* gd,
                                        gdlptr_sincos_omp_ggg gdl)
{
//    free_dmatrix3D(gdl->histZetaM,
//                   1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix3D(gdl->histZetaMsin,
                   1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix3D(gdl->histZetaMcos,
                   1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);

    free_dmatrix(gdl->histNNSubTC,1,cmd->sizeHistN,1,cmd->sizeHistN);

//    free_dvector(gdl->histNNSub,1,cmd->sizeHistN);

    return SUCCESS;
}

local int search_init_sincos_omp_ggg(struct  cmdline_data* cmd,
                                     struct  global_data* gd,
                                     gdhistptr_sincos_omp_ggg hist)
{
    int n;
    int m;
    INTEGER bytes_tot_local=0;

    //B 3pcf convergence & shear
    hist->ChebsT = dvector(1,cmd->mChebyshev+1);
    hist->ChebsU = dvector(1,cmd->mChebyshev+1);
    //E

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

//    for (n = 1; n <= cmd->sizeHistN; n++)
//        hist->histNNSubthread[n] = 0.0;

    for (m = 1; m <= cmd->mChebyshev+1; m++) {
        CLRM_ext(hist->histZetaMthreadcos[m], cmd->sizeHistN);
        CLRM_ext(hist->histZetaMthreadsin[m], cmd->sizeHistN);
        CLRM_ext(hist->histNNSubthreadTC, cmd->sizeHistN);
        CLRM_ext(hist->xiOUTVPcosTC[m], cmd->sizeHistN);
        CLRM_ext(hist->xiOUTVPsinTC[m], cmd->sizeHistN);
    }

    hist->nbbcalcthread = 0;
    hist->nbccalcthread = 0;

    return SUCCESS;
}

local int search_free_sincos_omp_ggg(struct  cmdline_data* cmd,
                                         struct  global_data* gd,
                                      gdhistptr_sincos_omp_ggg hist)
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

    //B 3pcf convergence & shear
    free_dvector(hist->ChebsU,1,cmd->mChebyshev+1);
    free_dvector(hist->ChebsT,1,cmd->mChebyshev+1);
    //E

    return SUCCESS;
}


//#if defined(NMultipoles) || defined(NONORMHIST)
local int search_init_gd_sincos_omp_ggg_N(struct  cmdline_data* cmd,
                                        struct  global_data* gd,
                                          gdlptr_sincos_omp_ggg_N gdl)
{
    int n;
    int m;

    INTEGER bytes_tot_local=0;
/*
//    gdl->histNNSub = dvector(1,cmd->sizeHistN);
//    bytes_tot_local += 1*cmd->sizeHistN*sizeof(real);

    gdl->histZetaMcos =
            dmatrix3D(1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    gdl->histZetaMsin =
            dmatrix3D(1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    bytes_tot_local +=
            4*(cmd->mChebyshev+1)*cmd->sizeHistN*cmd->sizeHistN*sizeof(real);

    gd->bytes_tot += bytes_tot_local;
    verb_print(cmd->verbose,
    "\nsearch_init_gd_octree_ggg: Allocated %g MByte for histograms storage.\n",
    bytes_tot_local*INMB);

//    for (n = 1; n <= cmd->sizeHistN; n++)
//        gdl->histNNSub[n] = 0.0;
    for (m = 1; m <= cmd->mChebyshev+1; m++) {
        CLRM_ext(gdl->histZetaMcos[m], cmd->sizeHistN);
        CLRM_ext(gdl->histZetaMsin[m], cmd->sizeHistN);
    }
*/

    //    real **histNNSubTC;
    //    real ***histZetaMcos;
    //    real ***histZetaMsin;

    //    gdl->histNNSub = dvector(1,cmd->sizeHistN);
        bytes_tot_local += 1*cmd->sizeHistN*sizeof(real);

        gdl->histNNSubTC = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);

        gdl->histZetaMcos =
                dmatrix3D(1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
        gdl->histZetaMsin =
                dmatrix3D(1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    //    gdl->histZetaM =
    //            dmatrix3D(1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
        bytes_tot_local +=
                5*(cmd->mChebyshev+1)*cmd->sizeHistN*cmd->sizeHistN*sizeof(real);


        gd->bytes_tot += bytes_tot_local;
        verb_print(cmd->verbose,
        "\nsearch_init_gd_octree_ggg_N: Allocated %g MByte for histograms storage.\n",
        bytes_tot_local*INMB);

    //    for (n = 1; n <= cmd->sizeHistN; n++)
    //        gdl->histNNSub[n] = 0.0;

        CLRM_ext(gdl->histNNSubTC, cmd->sizeHistN);

        for (m = 1; m <= cmd->mChebyshev+1; m++) {
            CLRM_ext(gdl->histZetaMcos[m], cmd->sizeHistN);
            CLRM_ext(gdl->histZetaMsin[m], cmd->sizeHistN);
    //        CLRM_ext(gdl->histZetaM[m], cmd->sizeHistN);
        }

        gd->nbbcalc = gd->nbccalc = gd->ncccalc = 0;

    return SUCCESS;
}

local int search_free_gd_sincos_omp_ggg_N(struct  cmdline_data* cmd,
                                         struct  global_data* gd,
                                        gdlptr_sincos_omp_ggg_N gdl)
{
    /*
    free_dmatrix3D(gdl->histZetaMsin,
                   1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix3D(gdl->histZetaMcos,
                   1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);

//    free_dvector(gdl->histNNSub,1,cmd->sizeHistN);
*/
    
    //    free_dmatrix3D(gdl->histZetaM,
    //                   1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
        free_dmatrix3D(gdl->histZetaMsin,
                       1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
        free_dmatrix3D(gdl->histZetaMcos,
                       1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);

        free_dmatrix(gdl->histNNSubTC,1,cmd->sizeHistN,1,cmd->sizeHistN);

    //    free_dvector(gdl->histNNSub,1,cmd->sizeHistN);


    return SUCCESS;
}

local int search_init_sincos_omp_ggg_N(struct  cmdline_data* cmd,
                                     struct  global_data* gd,
                                     gdhistptr_sincos_omp_ggg_N hist)
{
    int n;
    int m;
    INTEGER bytes_tot_local=0;

    //B 3pcf convergence & shear
    hist->ChebsT = dvector(1,cmd->mChebyshev+1);
    hist->ChebsU = dvector(1,cmd->mChebyshev+1);
    //E

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
        CLRM_ext(hist->histNNSubthreadTC, cmd->sizeHistN);
        CLRM_ext(hist->xiOUTVPcosTC[m], cmd->sizeHistN);
        CLRM_ext(hist->xiOUTVPsinTC[m], cmd->sizeHistN);
    }

    hist->nbbcalcthread = 0;
    hist->nbccalcthread = 0;

    return SUCCESS;
}

local int search_free_sincos_omp_ggg_N(struct  cmdline_data* cmd,
                                         struct  global_data* gd,
                                      gdhistptr_sincos_omp_ggg_N hist)
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

    //B 3pcf convergence & shear
    free_dvector(hist->ChebsU,1,cmd->mChebyshev+1);
    free_dvector(hist->ChebsT,1,cmd->mChebyshev+1);
    //E


    return SUCCESS;
}

//#endif // ! NMultipoles

//E Routines as in cballsutils

local int print_info(struct cmdline_data* cmd,
                                  struct  global_data* gd)
{
    verb_print(cmd->verbose,
               "searchcalc: Using octree-ggg-omp-triangles... \n");

    if (cballs_opt_ggg_correlation(cmd))
        verb_print(cmd->verbose, "computing GGG correlation... \n");


    if (cmd->usePeriodic==TRUE)
        error("CheckParameters: can´t have periodic boundaries and OCTREEGGGOMP definition (usePeriodic=%d)\nSet usePeriodic=false\n",
            cmd->usePeriodic);
    if (cmd->useLogHist==FALSE)
        error("CheckParameters: can´t have normal scale hist and OCTREEGGGOMP definition (useLogHist=%d)\nSet useLogHist=true\n",
            cmd->useLogHist);
    if (cmd->computeTPCF==FALSE)
        error("CheckParameters: can´t have computeTPCF=false and OCTREEGGGOMP definition (computeTPCF=%d)\nSet computeTPCF=true\n",
            cmd->computeTPCF);

    verb_print(cmd->verbose, "with 3pcf convergence computation... \n");

#if NDIM == 2
    error("CheckParameters: OCTREEGGGOMPTRIANGLES definition works only in a 3D unit sphere");
#endif

//#ifdef NMultipoles
    verb_print(cmd->verbose, "with NMultipoles... \n");
//#else
//    verb_print(cmd->verbose, "without NMultipoles... \n");
//#endif
//#ifdef NONORMHIST
    verb_print(cmd->verbose, "with NONORMHIST... \n");
    if (cballs_opt_no_normalize_histzeta(cmd))
        verb_print(cmd->verbose, "with option no-normalize-HistZeta...\n");
//#else
//    verb_print(cmd->verbose, "without NONORMHIST... \n");
//#endif

//#if defined(NMultipoles) && defined(NONORMHIST)
    if (cballs_opt_no_normalize_histzeta(cmd)) {
        if (cballs_opt_edge_corrections(cmd))
            verb_print(cmd->verbose, "with option edge-corrections... \n");
    } else {
        if (cballs_opt_edge_corrections(cmd)) {
            verb_print(cmd->verbose,
                       "option edge-corrections only works with %s... \n",
                       "no-normalize-HistZeta option added");
            // Check freeing allocated memory...
            error("going out...\n");
        }
    }

/*#else
    if (cballs_opt_edge_corrections(cmd)) {
        verb_print(cmd->verbose,
                   "option edge-corrections only works with %s activated... \n",
                   "NMultipoles && NONORMHIST");
        // Check freeing allocated memory...
        error("going out...");
    }
#endif */

#ifndef USEGSL
    if (cballs_opt_edge_corrections(cmd))
        verb_print(cmd->verbose,
            "option edge-corrections is better computed with %s activated... \n",
            "USEGSL");
#endif

    if (cballs_opt_no_one_ball(cmd))
        verb_print(cmd->verbose, "with option no-one-ball... \n");
    if (cballs_opt_smooth_pivot(cmd))
        verb_print(cmd->verbose,
                   "with option smooth-pivot... rsmooth=%g\n",gd->rsmooth[0]);
    if (cballs_opt_default_rsmooth(cmd))
        verb_print(cmd->verbose, "with option default-rsmooth... \n");
    if (cballs_opt_fix_rsmooth(cmd))
        verb_print(cmd->verbose, "with option fix-rsmooth... \n");

    verb_print(cmd->verbose, "with TRIANGLES... \n");

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

#define MHISTZETA \
"%16.8e %16.8e %16.8e %16.8e %16.8e %16.8e\n"

#define MHISTZETAHEADER \
"# [1] rBins; [2] diagonal; [3] theta2=Nbins/4.0; [4] theta2=2.0*Nbins/4.0; \
[5] theta2=3.0*Nbins/4.0; [6] theta2=4.0*Nbins/4.0 - 1.0\n"


// Saves matrix ZetaM for each m multipole
local int PrintHistZetaM_sincos(struct  cmdline_data* cmd,
                                struct  global_data* gd,
                                gdlptr_sincos_omp_ggg gdl)
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


// Saves matrix ZetaM for each m multipole at a set of theta2 angles
local int PrintHistZetaMm_sincos(struct  cmdline_data* cmd,
                                struct  global_data* gd,
                                 gdlptr_sincos_omp_ggg gdl)
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

//#if defined(NMultipoles) || defined(NONORMHIST)
// Saves matrix ZetaM for each m multipole
local int PrintHistZetaM_sincos_N(struct  cmdline_data* cmd,
                                struct  global_data* gd,
                                  gdlptr_sincos_omp_ggg_N gdlN)
{
    int n1, n2, m;
    stream outstr;
    char namebuf[256];

    for (m=1; m<=cmd->mChebyshev+1; m++) {
        sprintf(namebuf, "%s%s_%d%s", gd->fpfnamehistZetaMFileName,
                "_cos_N", m, EXTFILES);
        verb_print_q(2, cmd->verbose, "Printing : to a file %s ...\n",namebuf);
        outstr = stropen(namebuf, "w!");
        for (n1=1; n1<=cmd->sizeHistN; n1++) {
            for (n2=1; n2<=cmd->sizeHistN; n2++) {
                fprintf(outstr,"%16.8e ",gdlN->histZetaMcos[m][n1][n2]);
            }
            fprintf(outstr,"\n");
        }
        fclose(outstr);
    }
    for (m=1; m<=cmd->mChebyshev+1; m++) {
        sprintf(namebuf, "%s%s_%d%s", gd->fpfnamehistZetaMFileName,
                "_sin_N", m, EXTFILES);
        verb_print_q(2, cmd->verbose, "Printing : to a file %s ...\n",namebuf);
        outstr = stropen(namebuf, "w!");
        for (n1=1; n1<=cmd->sizeHistN; n1++) {
            for (n2=1; n2<=cmd->sizeHistN; n2++) {
                fprintf(outstr,"%16.8e ",gdlN->histZetaMsin[m][n1][n2]);
            }
            fprintf(outstr,"\n");
        }
        fclose(outstr);
    }

    return SUCCESS;
}

// Saves matrix ZetaM for each m multipole at a set of theta2 angles
local int PrintHistZetaMm_sincos_N(struct  cmdline_data* cmd,
                                struct  global_data* gd,
                                   gdlptr_sincos_omp_ggg_N gdlN)
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
                "_cos_N", m, EXTFILES);
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
            Zeta = gdlN->histZetaMcos[m][n1][n1];
            Zeta2 = gdlN->histZetaMcos[m][n1][(int)(Nbins/4.0)];
            Zeta3 = gdlN->histZetaMcos[m][n1][(int)(2.0*Nbins/4.0)];
            Zeta4 = gdlN->histZetaMcos[m][n1][(int)(3.0*Nbins/4.0)];
            Zeta5 = gdlN->histZetaMcos[m][n1][(int)(4.0*Nbins/4.0 - 1.0)];
            fprintf(outstr,MHISTZETA,rBin,Zeta,Zeta2,Zeta3,Zeta4,Zeta5);
        }
        fclose(outstr);
    }
        
    for (m = 1; m <= cmd->mChebyshev+1; m++) {
        sprintf(namebuf, "%s%s_%d%s", gd->fpfnamemhistZetaMFileName,
                "_sin_N", m, EXTFILES);
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
                Zeta = gdlN->histZetaMsin[m][n1][n1];
                Zeta2 = gdlN->histZetaMsin[m][n1][(int)(Nbins/4.0)];
                Zeta3 = gdlN->histZetaMsin[m][n1][(int)(2.0*Nbins/4.0)];
                Zeta4 = gdlN->histZetaMsin[m][n1][(int)(3.0*Nbins/4.0)];
                Zeta5 = gdlN->histZetaMsin[m][n1][(int)(4.0*Nbins/4.0 - 1.0)];
                fprintf(outstr,MHISTZETA,rBin,Zeta,Zeta2,Zeta3,Zeta4,Zeta5);
        }
        fclose(outstr);
    }

    return SUCCESS;
}


//#ifdef NONORMHIST

// Saves matrix ZetaM for each m multipole
local int PrintHistZetaM_sincos_normalized(struct  cmdline_data* cmd,
                                           struct  global_data* gd,
                                           gdlptr_sincos_omp_ggg gdl,
                                           gdlptr_sincos_omp_ggg_N gdlN)
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
                fprintf(outstr,"%16.8e ",
                    gdl->histZetaMcos[m][n1][n2]/gdlN->histZetaMcos[1][n1][n2]);
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
                fprintf(outstr,"%16.8e ",
                    gdl->histZetaMsin[m][n1][n2]/gdlN->histZetaMcos[1][n1][n2]);
            }
            fprintf(outstr,"\n");
        }
        fclose(outstr);
    }

    return SUCCESS;
}

// Saves matrix ZetaM for each m multipole at a set of theta2 angles
//  normalization with N_0 = (N_histZetaMcos[1][n1][n1]
//                              + N_histZetaMsin[1][n1][n1])
//                         = N_histZetaMcos[1][n1][n1]
//  because N_histZetaMsin[1][n1][n1] = 0.
local int PrintHistZetaMm_sincos_normalized(struct  cmdline_data* cmd,
                                            struct  global_data* gd,
                                            gdlptr_sincos_omp_ggg gdl,
                                            gdlptr_sincos_omp_ggg_N gdlN)
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

    real Norm;

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
                    rbinlog = rlog10(cmd->rminHist)
                                + ((real)(n1)-0.5)*gd->deltaR;
                }
                rBin=rpow(10.0,rbinlog);
            } else {
                rBin = cmd->rminHist + ((real)n1-0.5)*gd->deltaR;
            }
            Zeta  = gdl->histZetaMcos[m][n1][n1]
                    / gdlN->histZetaMcos[1][n1][n1];
            Zeta2 = gdl->histZetaMcos[m][n1][(int)(Nbins/4.0)]
                     / gdlN->histZetaMcos[1][n1][(int)(Nbins/4.0)];
            Zeta3 = gdl->histZetaMcos[m][n1][(int)(2.0*Nbins/4.0)]
                    / gdlN->histZetaMcos[1][n1][(int)(2.0*Nbins/4.0)];
            Zeta4 = gdl->histZetaMcos[m][n1][(int)(3.0*Nbins/4.0)]
                    / gdlN->histZetaMcos[1][n1][(int)(3.0*Nbins/4.0)];
            Zeta5 = gdl->histZetaMcos[m][n1][(int)(4.0*Nbins/4.0 - 1.0)]
                    / gdlN->histZetaMcos[1][n1][(int)(4.0*Nbins/4.0 - 1.0)];
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

            Zeta = gdl->histZetaMsin[m][n1][n1]
                 / gdlN->histZetaMcos[1][n1][n1];
            Zeta2 = gdl->histZetaMsin[m][n1][(int)(Nbins/4.0)]
                    / gdlN->histZetaMcos[1][n1][(int)(Nbins/4.0)];
            Zeta3 = gdl->histZetaMsin[m][n1][(int)(2.0*Nbins/4.0)]
                    / gdlN->histZetaMcos[1][n1][(int)(2.0*Nbins/4.0)];
            Zeta4 = gdl->histZetaMsin[m][n1][(int)(3.0*Nbins/4.0)]
                    / gdlN->histZetaMcos[1][n1][(int)(3.0*Nbins/4.0)];
            Zeta5 = gdl->histZetaMsin[m][n1][(int)(4.0*Nbins/4.0 - 1.0)]
                    / gdlN->histZetaMcos[1][n1][(int)(4.0*Nbins/4.0 - 1.0)];
                fprintf(outstr,MHISTZETA,rBin,Zeta,Zeta2,Zeta3,Zeta4,Zeta5);
        }
        fclose(outstr);
    }

    return SUCCESS;
}

// Saves matrix ZetaM for each m multipole
local int PrintHistZetaM_sincos_edge_effects(struct  cmdline_data* cmd,
                                              struct  global_data* gd,
                                              gdlptr_sincos_omp_ggg gdl,
                                              gdlptr_sincos_omp_ggg_N gdlN)
{
    int n1, n2, m;
    stream outstr;
    char namebuf[256];
    real rBin, rbinlog;

    real ***mat3;
    real ***mat4;
    real ***mat5;
    mat3 = dmatrix3D(1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    mat4 = dmatrix3D(1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    mat5 = dmatrix3D(1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    for (m = 1; m <= cmd->mChebyshev+1; m++) {
        CLRM_ext(mat3[m], cmd->sizeHistN);
        CLRM_ext(mat4[m], cmd->sizeHistN);
        CLRM_ext(mat5[m], cmd->sizeHistN);
    }

    for (n1=1; n1<=cmd->sizeHistN; n1++) {
        for (n2=1; n2<=cmd->sizeHistN; n2++) {
            for (m=1; m<=cmd->mChebyshev+1; m++) {
                mat3[m][n1][n2] = gdl->histZetaMcos[m][n1][n2]
                                    + gdl->histZetaMsin[m][n1][n2];
                mat4[m][n1][n2] = gdlN->histZetaMcos[m][n1][n2]
                                    + gdlN->histZetaMsin[m][n1][n2];
            }
        }
    }

    for (n1=1; n1<=cmd->sizeHistN; n1++) {
        for (n2=1; n2<=cmd->sizeHistN; n2++) {
            matrixClm(cmd, gd, mat3, mat4, n1, n2, mat5);
            if (cmd->verbose_log>=3)  {
                verb_log_print(cmd->verbose_log, gd->outlog,
                               "\n\nhistZetaM elements again (%d, %d):\n\n",
                               n1, n2);
                for (m=1; m<=cmd->mChebyshev+1; m++) {
                        verb_log_print(cmd->verbose_log, gd->outlog,
                                       "%g\n",
                                       mat5[m][n1][n2]);
                }
            }

        }
    }

    for (m=1; m<=cmd->mChebyshev+1; m++) {
        sprintf(namebuf, "%s%s_%d%s", gd->fpfnamehistZetaMFileName,
                "_EE", m, EXTFILES);
        verb_print_q(2, cmd->verbose, "Printing : to a file %s ...\n",namebuf);
        outstr = stropen(namebuf, "w!");
        for (n1=1; n1<=cmd->sizeHistN; n1++) {
            for (n2=1; n2<=cmd->sizeHistN; n2++) {
                fprintf(outstr,"%16.8e ",mat5[m][n1][n2]);
            }
            fprintf(outstr,"\n");
        }
        fclose(outstr);
    }

    //B  and saves matrix ZetaM for each m multipole at a set of theta2 angles
    if (cballs_opt_out_m_histzeta(cmd)) {
        real Zeta;
        real Zeta2;
        real Zeta3;
        real Zeta4;
        real Zeta5;
        int Nbins;
        
        Nbins = cmd->sizeHistN;
        for (m = 1; m <= cmd->mChebyshev+1; m++) {
            sprintf(namebuf, "%s%s_%d%s", gd->fpfnamemhistZetaMFileName,
                    "_EE", m, EXTFILES);
            verb_print_q(2, cmd->verbose, "Printing : to a file %s ...\n",
                         namebuf);
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
                Zeta = mat5[m][n1][n1];
                Zeta2 = mat5[m][n1][(int)(Nbins/4.0)];
                Zeta3 = mat5[m][n1][(int)(2.0*Nbins/4.0)];
                Zeta4 = mat5[m][n1][(int)(3.0*Nbins/4.0)];
                Zeta5 = mat5[m][n1][(int)(4.0*Nbins/4.0 - 1.0)];
                fprintf(outstr,MHISTZETA,rBin,Zeta,Zeta2,Zeta3,Zeta4,Zeta5);
            }
            fclose(outstr);
        }
    }
    //E

    free_dmatrix3D(mat5,
                   1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix3D(mat4,
                   1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix3D(mat3,
                   1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);

    return SUCCESS;
}

//#endif // ! NONORMHIST

//#endif // ! NMultipoles

#undef MHISTZETAHEADER
#undef MHISTZETA

//E Saving histograms section: case GGGCORRELATION:

