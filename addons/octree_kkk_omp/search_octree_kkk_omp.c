/* ==============================================================================
 MODULE: search_octree_kkk_omp.c		[cTreeBalls]
 Written by: M.A. Rodriguez-Meza
 Starting date:    april 2023
 Purpose: 3-point correlation function computation
 Language: C
 Use: searchcalc_octree_kkk_omp(cmd, gd, btable, nbody,
                                           ipmin, ipmax, cat1, cat2);
 Major revisions:
 ==============================================================================*/
//        1          2          3          4        ^ 5          6          7


// Work to do in order to use with boxes not centered at (0,0,...)

#include "globaldefs.h"

#ifdef MANUALCHEBYSHEV
#define CHEBYSHEVTUOMP                                            \
{REAL xicosmphi,xisinmphi; int m;                                 \
    REAL cosphi2, cosphi3, cosphi4;                               \
    cosphi2 = cosphi*cosphi; cosphi3 = cosphi2*cosphi;            \
    cosphi4 = cosphi3*cosphi;                                     \
    hist->ChebsT[1] = 1.0;                                        \
    xicosmphi = xi * hist->ChebsT[1];                             \
    hist->histXithreadcos[1][n] += xicosmphi;                     \
    hist->ChebsT[2] = cosphi;                                     \
    xicosmphi = xi * hist->ChebsT[2];                             \
    hist->histXithreadcos[2][n] += xicosmphi;                     \
    hist->ChebsT[3] = 2.0*cosphi2 - (1.0);                        \
    xicosmphi = xi * hist->ChebsT[3];                             \
    hist->histXithreadcos[3][n] += xicosmphi;                     \
    hist->ChebsT[4] = -3.0*cosphi + 4.0*cosphi3;                  \
    xicosmphi = xi * hist->ChebsT[4];                             \
    hist->histXithreadcos[4][n] += xicosmphi;                     \
    hist->ChebsT[5] = 1.0 - 8.0*cosphi2 + 8.0*cosphi4;            \
    xicosmphi = xi * hist->ChebsT[5];                             \
    hist->histXithreadcos[5][n] += xicosmphi;                     \
    hist->ChebsT[6] = 5.0*cosphi - 20.0*cosphi3 + 16.0*cosphi4*cosphi; \
    xicosmphi = xi * hist->ChebsT[6];                             \
    hist->histXithreadcos[6][n] += xicosmphi;                     \
    hist->ChebsU[1] = 0.0;                                        \
    xisinmphi = xi * hist->ChebsU[1] * sinphi;                    \
    hist->histXithreadsin[1][n] += xisinmphi;                     \
    hist->ChebsU[2] = 1.0;                                        \
    xisinmphi = xi * hist->ChebsU[2] * sinphi;                    \
    hist->histXithreadsin[2][n] += xisinmphi;                     \
    hist->ChebsU[3] = 2.0*cosphi;                                 \
    xisinmphi = xi * hist->ChebsU[3] * sinphi;                    \
    hist->histXithreadsin[3][n] += xisinmphi;                     \
    hist->ChebsU[4] = -1.0 + 4.0*cosphi2;                         \
    xisinmphi = xi * hist->ChebsU[4] * sinphi;                    \
    hist->histXithreadsin[4][n] += xisinmphi;                     \
    hist->ChebsU[5] = -4.0*cosphi + 8.0*cosphi3;                  \
    xisinmphi = xi * hist->ChebsU[5] * sinphi;                    \
    hist->histXithreadsin[5][n] += xisinmphi;                     \
    hist->ChebsU[6] = 1.0 -12.0*cosphi2 + 16.0*cosphi4;           \
    xisinmphi = xi * hist->ChebsU[6] * sinphi;                    \
    hist->histXithreadsin[6][n] += xisinmphi;                     \
    for (m=7; m<=cmd->mChebyshev+1; m++){                         \
        hist->ChebsT[m] = 2.0*(cosphi)*hist->ChebsT[m-1] - hist->ChebsT[m-2]; \
        xicosmphi = xi * hist->ChebsT[m];                         \
        hist->histXithreadcos[m][n] += xicosmphi;                 \
        hist->ChebsU[m] = 2.0*(cosphi)*hist->ChebsU[m-1] - hist->ChebsU[m-2]; \
        xisinmphi = xi * hist->ChebsU[m] * sinphi;                \
        hist->histXithreadsin[m][n] += xisinmphi;                 \
    }}
#else
#define CHEBYSHEVTUOMP                                            \
{real xicosmphi,xisinmphi; int m;                                 \
    hist->ChebsT[1] = 1.0;                                        \
    xicosmphi = xi * hist->ChebsT[1];                             \
    hist->histXithreadcos[1][n] += xicosmphi;                     \
    hist->ChebsT[2] = cosphi;                                     \
    xicosmphi = xi * hist->ChebsT[2];                             \
    hist->histXithreadcos[2][n] += xicosmphi;                     \
    hist->ChebsT[3] = 2.0*(cosphi)*(cosphi) - (1.0);              \
    xicosmphi = xi * hist->ChebsT[3];                             \
    hist->histXithreadcos[3][n] += xicosmphi;                     \
    hist->ChebsU[1] = 0.0;                                        \
    xisinmphi = xi * hist->ChebsU[1] * sinphi;                    \
    hist->histXithreadsin[1][n] += xisinmphi;                     \
    hist->ChebsU[2] = 1.0;                                        \
    xisinmphi = xi * hist->ChebsU[2] * sinphi;                    \
    hist->histXithreadsin[2][n] += xisinmphi;                     \
    hist->ChebsU[3] = 2.0*cosphi;                                 \
    xisinmphi = xi * hist->ChebsU[3] * sinphi;                    \
    hist->histXithreadsin[3][n] += xisinmphi;                     \
    for (m=4; m<=cmd->mChebyshev+1; m++){                         \
        hist->ChebsT[m] = 2.0*(cosphi)*hist->ChebsT[m-1] - hist->ChebsT[m-2]; \
        xicosmphi = xi * hist->ChebsT[m];                         \
        hist->histXithreadcos[m][n] += xicosmphi;                 \
        hist->ChebsU[m] = 2.0*(cosphi)*hist->ChebsU[m-1] - hist->ChebsU[m-2]; \
        xisinmphi = xi * hist->ChebsU[m] * sinphi;                \
        hist->histXithreadsin[m][n] += xisinmphi;                 \
    }}
#endif

#ifdef NMultipoles
#ifdef MANUALCHEBYSHEV
#define CHEBYSHEVTUOMPN                                           \
{REAL xicosmphi,xisinmphi; int m;                                 \
    REAL cosphi2, cosphi3, cosphi4;                               \
    cosphi2 = cosphi*cosphi; cosphi3 = cosphi2*cosphi;            \
    cosphi4 = cosphi3*cosphi;                                     \
    histN->ChebsT[1] = 1.0;                                       \
    xicosmphi = xiN * histN->ChebsT[1];                           \
    histN->histXithreadcos[1][n] += xicosmphi;                    \
    histN->ChebsT[2] = cosphi;                                    \
    xicosmphi = xiN * histN->ChebsT[2];                           \
    histN->histXithreadcos[2][n] += xicosmphi;                    \
    histN->ChebsT[3] = 2.0*cosphi2 - (1.0);                       \
    xicosmphi = xiN * histN->ChebsT[3];                           \
    histN->histXithreadcos[3][n] += xicosmphi;                    \
    histN->ChebsT[4] = -3.0*cosphi + 4.0*cosphi3;                 \
    xicosmphi = xiN * histN->ChebsT[4];                           \
    histN->histXithreadcos[4][n] += xicosmphi;                    \
    histN->ChebsT[5] = 1.0 - 8.0*cosphi2 + 8.0*cosphi4;           \
    xicosmphi = xiN * histN->ChebsT[5];                           \
    histN->histXithreadcos[5][n] += xicosmphi;                    \
    histN->ChebsT[6] = 5.0*cosphi - 20.0*cosphi3 + 16.0*cosphi4*cosphi; \
    xicosmphi = xiN * histN->ChebsT[6];                           \
    histN->histXithreadcos[6][n] += xicosmphi;                    \
    histN->ChebsU[1] = 0.0;                                       \
    xisinmphi = xiN * histN->ChebsU[1] * sinphi;                  \
    histN->histXithreadsin[1][n] += xisinmphi;                    \
    histN->ChebsU[2] = 1.0;                                       \
    xisinmphi = xiN * histN->ChebsU[2] * sinphi;                  \
    histN->histXithreadsin[2][n] += xisinmphi;                    \
    histN->ChebsU[3] = 2.0*cosphi;                                \
    xisinmphi = xiN * histN->ChebsU[3] * sinphi;                  \
    histN->histXithreadsin[3][n] += xisinmphi;                    \
    histN->ChebsU[4] = -1.0 + 4.0*cosphi2;                        \
    xisinmphi = xiN * histN->ChebsU[4] * sinphi;                  \
    histN->histXithreadsin[4][n] += xisinmphi;                    \
    histN->ChebsU[5] = -4.0*cosphi + 8.0*cosphi3;                 \
    xisinmphi = xiN * histN->ChebsU[5] * sinphi;                  \
    histN->histXithreadsin[5][n] += xisinmphi;                    \
    histN->ChebsU[6] = 1.0 -12.0*cosphi2 + 16.0*cosphi4;          \
    xisinmphi = xiN * histN->ChebsU[6] * sinphi;                  \
    histN->histXithreadsin[6][n] += xisinmphi;                    \
    for (m=7; m<=cmd->mChebyshev+1; m++){                         \
        histN->ChebsT[m] = 2.0*(cosphi)*histN->ChebsT[m-1]-histN->ChebsT[m-2]; \
        xicosmphi = xiN * histN->ChebsT[m];                       \
        histN->histXithreadcos[m][n] += xicosmphi;                \
        histN->ChebsU[m] = 2.0*(cosphi)*histN->ChebsU[m-1]-histN->ChebsU[m-2]; \
        xisinmphi = xiN * histN->ChebsU[m] * sinphi;              \
        histN->histXithreadsin[m][n] += xisinmphi;                \
    }}
#else
#define CHEBYSHEVTUOMPN                                           \
{real xicosmphi,xisinmphi; int m;                                 \
    histN->ChebsT[1] = 1.0;                                       \
    xicosmphi = xiN * histN->ChebsT[1];                           \
    histN->histXithreadcos[1][n] += xicosmphi;                    \
    histN->ChebsT[2] = cosphi;                                    \
    xicosmphi = xiN * histN->ChebsT[2];                           \
    histN->histXithreadcos[2][n] += xicosmphi;                    \
    histN->ChebsT[3] = 2.0*(cosphi)*(cosphi) - (1.0);             \
    xicosmphi = xiN * histN->ChebsT[3];                           \
    histN->histXithreadcos[3][n] += xicosmphi;                    \
    histN->ChebsU[1] = 0.0;                                       \
    xisinmphi = xiN * histN->ChebsU[1] * sinphi;                  \
    histN->histXithreadsin[1][n] += xisinmphi;                    \
    histN->ChebsU[2] = 1.0;                                       \
    xisinmphi = xiN * histN->ChebsU[2] * sinphi;                  \
    histN->histXithreadsin[2][n] += xisinmphi;                    \
    histN->ChebsU[3] = 2.0*cosphi;                                \
    xisinmphi = xiN * histN->ChebsU[3] * sinphi;                  \
    histN->histXithreadsin[3][n] += xisinmphi;                    \
    for (m=4; m<=cmd->mChebyshev+1; m++){                         \
        histN->ChebsT[m] = 2.0*(cosphi)*histN->ChebsT[m-1]-histN->ChebsT[m-2]; \
        xicosmphi = xiN * histN->ChebsT[m];                       \
        histN->histXithreadcos[m][n] += xicosmphi;                \
        histN->ChebsU[m] = 2.0*(cosphi)*histN->ChebsU[m-1]-histN->ChebsU[m-2]; \
        xisinmphi = xiN * histN->ChebsU[m] * sinphi;              \
        histN->histXithreadsin[m][n] += xisinmphi;                \
    }}
#endif
#endif

typedef struct {
    real **xiOUTVPcos;
    real **xiOUTVPsin;
    real **xiOUTVPsincos;
    real **xiOUTVPcossin;
    real **histZetaMtmpcos;
    real **histZetaMtmpsin;
    real **histZetaMtmpsincos;
    real **histZetaMtmpcossin;
    real *ChebsT;
    real *ChebsU;

    real ***histZetaMthreadcos;
    real ***histZetaMthreadsin;
    real ***histZetaMthreadsincos;
    real ***histZetaMthreadcossin;
    realptr histNthread;
    realptr histNNSubthread;

    real **histXithreadcos;
    real **histXithreadsin;

    INTEGER nbbcalcthread;
    INTEGER nbccalcthread;

    vector q0;
    real drpq2, drpq;
    vector dr0;
    real cosb;
    real sinb;
} gdhist_sincos_omp_kkk, *gdhistptr_sincos_omp_kkk;

#ifdef NMultipoles
typedef struct {
    real **xiOUTVPcos;
    real **xiOUTVPsin;
    real **xiOUTVPsincos;
    real **xiOUTVPcossin;
    real **histZetaMtmpcos;
    real **histZetaMtmpsin;
    real **histZetaMtmpsincos;
    real **histZetaMtmpcossin;
    real *ChebsT;
    real *ChebsU;

    real ***histZetaMthreadcos;
    real ***histZetaMthreadsin;
    real ***histZetaMthreadsincos;
    real ***histZetaMthreadcossin;
    realptr histNthread;
    realptr histNNSubthread;

    real **histXithreadcos;
    real **histXithreadsin;

    INTEGER nbbcalcthread;
    INTEGER nbccalcthread;

    vector q0;
    real drpq2, drpq;
    vector dr0;
} gdhist_sincos_omp_kkk_N, *gdhistptr_sincos_omp_kkk_N;
#endif

local void normal_walktree_sincos(struct  cmdline_data* cmd,
                                  struct  global_data* gd,
                                  bodyptr *btable, int cat2,
                                  bodyptr, nodeptr, real,
                                  gdhistptr_sincos_omp_kkk, int *, int *);
local void sumnode_sincos(struct  cmdline_data* cmd,
                          struct  global_data* gd,
                          bodyptr *btable, int cat2,
                          bodyptr, cellptr, cellptr,
                          gdhistptr_sincos_omp_kkk, int *, int *);
local void sumnode_sincos_cell(struct  cmdline_data* cmd,
                               struct  global_data* gd,
                               bodyptr *btable, int cat2,
                               bodyptr p, cellptr start, cellptr finish,
                               gdhistptr_sincos_omp_kkk hist,
                               int *nbList, int *intList);

local int search_init_gd_sincos_omp_kkk(struct  cmdline_data* cmd,
                                        struct  global_data* gd);
local int search_init_sincos_omp_kkk(struct  cmdline_data* cmd,
                                  struct  global_data* gd,
                                     gdhistptr_sincos_omp_kkk hist);
local int search_free_sincos_omp_kkk(struct  cmdline_data* cmd,
                                         struct  global_data* gd,
                                     gdhistptr_sincos_omp_kkk hist);
local int computeBodyProperties_sincos_kkk(struct  cmdline_data* cmd,
                                            struct  global_data* gd,
                                            bodyptr p, int nbody,
                                       gdhistptr_sincos_omp_kkk hist);
local int print_info(struct cmdline_data* cmd,
                     struct  global_data* gd);

#ifdef NMultipoles
local void normal_walktree_sincos_N(struct  cmdline_data* cmd,
                                    struct  global_data* gd,
                                    bodyptr *btable, int cat2,
                                    bodyptr, nodeptr, real,
                                    gdhistptr_sincos_omp_kkk,
                                    gdhistptr_sincos_omp_kkk_N,  int *, int *);
local void sumnode_sincos_N(struct  cmdline_data* cmd,
                            struct  global_data* gd,
                            bodyptr *btable, int cat2,
                            bodyptr, cellptr, cellptr,
                            gdhistptr_sincos_omp_kkk,
                            gdhistptr_sincos_omp_kkk_N, int *, int *);
local void sumnode_sincos_cell_N(struct  cmdline_data*,
                                 struct  global_data*,
                                 bodyptr *btable, int cat2,
                                 bodyptr, cellptr, cellptr,
                                 gdhistptr_sincos_omp_kkk,
                                 gdhistptr_sincos_omp_kkk_N,
                                 int *nbList, int *intList);
local int search_init_gd_sincos_omp_kkk_N(struct  cmdline_data* cmd,
                                        struct  global_data* gd);
local int search_init_sincos_omp_kkk_N(struct  cmdline_data* cmd,
                                  struct  global_data* gd,
                                     gdhistptr_sincos_omp_kkk_N hist);
local int search_free_sincos_omp_kkk_N(struct  cmdline_data* cmd,
                                         struct  global_data* gd,
                                     gdhistptr_sincos_omp_kkk_N);
local int computeBodyProperties_sincos_kkk_N(struct  cmdline_data* cmd,
                                            struct  global_data* gd,
                                            bodyptr p, int nbody,
                                       gdhistptr_sincos_omp_kkk_N);
#endif

//B kappa Avg Rmin
#ifdef DEBUG
local char pivotsfilePath[MAXLENGTHOFFILES];
local FILE *outpivots;
#endif

#ifdef ADDPIVOTNEIGHBOURS
#define FACTIVENB  10
#define FACTIVEINT  10
 
local int actlenNb;
local int *activeNb;

local int actlenInt;
local int *activeInt;

local void sumnode_nblist_omp(struct cmdline_data* cmd,
                              struct  global_data* gd,
                              bodyptr *btable,
                              INTEGER ipmin, INTEGER *ipmax,
                              int cat1, int cat2,
                              bodyptr p,
                              gdhistptr_sincos_omp_kkk hist, int);

#endif

//E

/*
 Search routine using octree method:

 To be called using: search=octree-kkk-omp

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
global int searchcalc_octree_kkk_omp(struct cmdline_data* cmd,
                                             struct  global_data* gd,
                                             bodyptr *btable, INTEGER *nbody,
                                             INTEGER ipmin, INTEGER *ipmax, 
                                             int cat1, int cat2)
{
    double cpustart;

    cpustart = CPUTIME;
    print_info(cmd, gd);

//B kappa Avg Rmin
#ifdef DEBUG
    sprintf(pivotsfilePath,"%s/pivot_info%s.txt",gd->tmpDir,cmd->suffixOutFiles);
    if(!(outpivots=fopen(pivotsfilePath, "w")))
        error("\nsearchcalc_tc_kkk_omp: error opening file '%s' \n",
              pivotsfilePath);
#endif
//E

#ifdef OPENMPCODE
    ThreadCount(cmd, gd, nbody[cat1], cat1);
#endif

    search_init_gd_sincos_omp_kkk(cmd, gd);
//    verb_print_debug(1, "\nAqui voy (1)\n");
#ifdef NMultipoles
    search_init_gd_sincos_omp_kkk_N(cmd, gd);
#endif

    //B kappa Avg Rmin
    INTEGER ipfalse;
    ipfalse=0;
    INTEGER icountNbRmin;
    icountNbRmin=0;
    INTEGER icountNbRminOverlap;
    icountNbRminOverlap=0;

#ifdef ADDPIVOTNEIGHBOURS
    if (scanopt(cmd->options, "smooth-pivot")) {
        //B Alloc memory for neighbour lists
        int Nrsmooth;
        Nrsmooth = 0.25*nbody[cat2]*rsqr(gd->rsmooth[0]);
        actlenNb = FACTIVENB * Nrsmooth;
        verb_print(cmd->verbose, "\n- Nrsmooth and actlenNb: %d %d\n",
                   Nrsmooth, actlenNb);
        verb_log_print(cmd->verbose,gd->outlog,
                       "searchcalc: actlenNb = %ld\n",actlenNb);
        activeNb = (int *) allocate(actlenNb * sizeof(int));
        //E
        //B Alloc memory for intercation lists
        int NrangeN;
        NrangeN = 0.25*nbody[cat2]*rsqr(cmd->rangeN);
        actlenInt = FACTIVEINT * NrangeN;
        verb_print(cmd->verbose, "- NrangeN and actlenInt: %d %d\n",
                   NrangeN, actlenInt);
        verb_log_print(cmd->verbose,gd->outlog,
                       "searchcalc: actlenInt = %ld\n",actlenInt);
        activeInt = (int *) allocate(actlenInt * sizeof(int));
        //E
    }
#endif
    //E

    verb_print(cmd->verbose,
        "\nsearchcalc_tc_kkk_omp: Total allocated %g MByte storage so far.\n",
        gd->bytes_tot/(1024.0*1024.0));

    if (cmd->verbose >= VERBOSENORMALINFO)
        verb_print(cmd->verbose,
                   "\nRunning...\n - Completed pivot node:\n");

#ifdef DEBUG
#ifdef ADDPIVOTNEIGHBOURS
#pragma omp parallel default(none)                                          \
    shared(cmd,gd,btable,nbody,roottable,outpivots,                         \
           actlenNb,activeNb,actlenInt,activeInt,                           \
           ipmin,ipmax,cat1,cat2,ipfalse,icountNbRmin,icountNbRminOverlap)
#else
#pragma omp parallel default(none)                                          \
    shared(cmd,gd,btable,nbody,roottable,outpivots,                         \
           ipmin,ipmax,cat1,cat2,ipfalse,icountNbRmin,icountNbRminOverlap)
#endif
#else
#ifdef ADDPIVOTNEIGHBOURS
#pragma omp parallel default(none)                                          \
    shared(cmd,gd,btable,nbody,roottable,                         \
           actlenNb,activeNb,actlenInt,activeInt,                           \
           ipmin,ipmax,cat1,cat2,ipfalse,icountNbRmin,icountNbRminOverlap)
#else
#pragma omp parallel default(none)                                          \
    shared(cmd,gd,btable,nbody,roottable,                         \
           ipmin,ipmax,cat1,cat2,ipfalse,icountNbRmin,icountNbRminOverlap)
#endif
#endif
  {
    bodyptr p;
    bodyptr q;
    int n, m, ip;
    int i;

    //B init:
    gdhist_sincos_omp_kkk hist;
    search_init_sincos_omp_kkk(cmd, gd, &hist);
#ifdef NMultipoles
    gdhist_sincos_omp_kkk_N histN;
    search_init_sincos_omp_kkk_N(cmd, gd, &histN);
#endif
    //E

    //B kappa Avg Rmin
    INTEGER ipfalsethreads;
    ipfalsethreads = 0;
    INTEGER icountNbRminthread;
    icountNbRminthread=0;
    INTEGER icountNbRminOverlapthread;
    icountNbRminOverlapthread=0;

//#ifdef ADDPIVOTNEIGHBOURS
    int nbList;
    int intList;
//#endif
    //E

#pragma omp for nowait schedule(dynamic)
      for (p = btable[cat1] + ipmin -1; p < btable[cat1] + ipmax[cat1]; p++) {
          //B kappa Avg Rmin
          NbRmin(p) = 1;
          NbRminOverlap(p) = 0;
          KappaRmin(p) = Kappa(p);
//#ifdef ADDPIVOTNEIGHBOURS
          nbList=0;
          intList=0;
//#endif
          if (scanopt(cmd->options, "smooth-pivot")) {
              if (Update(p) == FALSE) {
                  ipfalsethreads++;
                  continue;
              }
          }
          //E

//B segment to be included below...
          //B Set histograms to zero for the pivot
          for (n = 1; n <= cmd->sizeHistN; n++)
              hist.histNNSubthread[n] = 0.0;
          CLRM_ext_ext(hist.histXithreadcos,
                       cmd->mChebyshev+1, cmd->sizeHistN);
          CLRM_ext_ext(hist.histXithreadsin,
                       cmd->mChebyshev+1, cmd->sizeHistN);
#ifdef NMultipoles
          for (n = 1; n <= cmd->sizeHistN; n++)
              histN.histNNSubthread[n] = 0.0;
          CLRM_ext_ext(histN.histXithreadcos,
                       cmd->mChebyshev+1, cmd->sizeHistN);
          CLRM_ext_ext(histN.histXithreadsin,
                       cmd->mChebyshev+1, cmd->sizeHistN);
#endif
          //E
          //B Set reference axis...
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
          //E
//E

#ifdef NMultipoles
          normal_walktree_sincos_N(cmd, gd, btable, cat2,
                                   p, ((nodeptr) roottable[cat2]),
                                   gd->rSizeTable[cat2], &hist, &histN,
                                   &nbList, &intList);
          computeBodyProperties_sincos_kkk(cmd, gd, p, nbody[cat1], &hist);
          computeBodyProperties_sincos_kkk_N(cmd, gd, p, nbody[cat1], &histN);

#ifdef ADDPIVOTNEIGHBOURS
          if (scanopt(cmd->options, "smooth-pivot")) {
              if (cmd->verbose_log>=3)
                  verb_log_print(cmd->verbose_log, gd->outlog,
                                 " - Summing nbList: %ld\n", nbList);
              for (i = 0; i < nbList; i++) {        // loop over neighbours
                  if (cmd->verbose_log>=3)
                      verb_log_print(cmd->verbose_log, gd->outlog,
                                     " - Summing intList: %ld\n", intList);
                  q = btable[cat1] + activeNb[i];
                  //B Set histograms to zero for the pivot
                  for (n = 1; n <= cmd->sizeHistN; n++)
                      hist.histNNSubthread[n] = 0.0;
                  CLRM_ext_ext(hist.histXithreadcos,
                               cmd->mChebyshev+1, cmd->sizeHistN);
                  CLRM_ext_ext(hist.histXithreadsin,
                               cmd->mChebyshev+1, cmd->sizeHistN);
#ifdef NMultipoles
                  for (n = 1; n <= cmd->sizeHistN; n++)
                      histN.histNNSubthread[n] = 0.0;
                  CLRM_ext_ext(histN.histXithreadcos,
                               cmd->mChebyshev+1, cmd->sizeHistN);
                  CLRM_ext_ext(histN.histXithreadsin,
                               cmd->mChebyshev+1, cmd->sizeHistN);
#endif
                  //E
//B Set reference axis...
#ifdef POLARAXIS
                  hist.q0[0] = 0.0;
                  hist.q0[1] = 0.0;
                  hist.q0[2] = 1.0;
                  DOTPSUBV(hist.drpq2, hist.dr0, Pos(q), hist.q0);
                  hist.drpq = rsqrt(hist.drpq2);
                  real b = 2.0*rasin(hist.drpq/2.0);
                  hist.cosb = rcos(b);
                  hist.sinb = rsin(b);
                  if (hist.drpq2==0) continue;
#else
                  dRotation3D(Pos(q), ROTANGLE, ROTANGLE, ROTANGLE, hist.q0);
                  DOTPSUBV(hist.drpq2, hist.dr0, Pos(q), hist.q0);
                  hist.drpq = rsqrt(hist.drpq2);
#endif
//E
                  sumnode_nblist_omp(cmd, gd, btable, ipmin, ipmax, cat1, cat2,
                                     q, &hist, intList);
                  computeBodyProperties_sincos_kkk(cmd, gd, q,
                                                   ipmax[cat1]-ipmin+1, &hist);
                  computeBodyProperties_sincos_kkk_N(cmd, gd, q,
                                                     ipmax[cat1]-ipmin+1, &histN);
              } // ! end i loop
          } // ! scanoption smooth-pivot
#endif // ! ADDPIVOTNEIGHBOURS


#else // ! NMultipoles
          normal_walktree_sincos(cmd, gd, btable, cat2,
                                 p, ((nodeptr) roottable[cat2]),
                                 gd->rSizeTable[cat2], &hist, &nbList, &intList);
          computeBodyProperties_sincos_kkk(cmd, gd, p,
                                           ipmax[cat1]-ipmin+1, &hist);
#ifdef ADDPIVOTNEIGHBOURS
          if (scanopt(cmd->options, "smooth-pivot")) {
              if (cmd->verbose_log>=3)
                  verb_log_print(cmd->verbose_log, gd->outlog,
                                 " - Summing nbList: %ld\n", nbList);
              for (i = 0; i < nbList; i++) {        // loop over neighbours
                  if (cmd->verbose_log>=3)
                      verb_log_print(cmd->verbose_log, gd->outlog,
                                     " - Summing intList: %ld\n", intList);
                  q = btable[cat1] + activeNb[i];
                  //B Set histograms to zero for the pivot
                  for (n = 1; n <= cmd->sizeHistN; n++)
                      hist.histNNSubthread[n] = 0.0;
                  CLRM_ext_ext(hist.histXithreadcos,
                               cmd->mChebyshev+1, cmd->sizeHistN);
                  CLRM_ext_ext(hist.histXithreadsin,
                               cmd->mChebyshev+1, cmd->sizeHistN);
#ifdef NMultipoles
                  for (n = 1; n <= cmd->sizeHistN; n++)
                      histN.histNNSubthread[n] = 0.0;
                  CLRM_ext_ext(histN.histXithreadcos,
                               cmd->mChebyshev+1, cmd->sizeHistN);
                  CLRM_ext_ext(histN.histXithreadsin,
                               cmd->mChebyshev+1, cmd->sizeHistN);
#endif
                  //E
//B Set reference axis...
#ifdef POLARAXIS
                  hist.q0[0] = 0.0;
                  hist.q0[1] = 0.0;
                  hist.q0[2] = 1.0;
                  DOTPSUBV(hist.drpq2, hist.dr0, Pos(q), hist.q0);
                  hist.drpq = rsqrt(hist.drpq2);
                  real b = 2.0*rasin(hist.drpq/2.0);
                  hist.cosb = rcos(b);
                  hist.sinb = rsin(b);
                  if (hist.drpq2==0) continue;
#else
                  dRotation3D(Pos(q), ROTANGLE, ROTANGLE, ROTANGLE, hist.q0);
                  DOTPSUBV(hist.drpq2, hist.dr0, Pos(q), hist.q0);
                  hist.drpq = rsqrt(hist.drpq2);
#endif
//E
                  sumnode_nblist_omp(cmd, gd, btable, ipmin, ipmax, cat1, cat2,
                                     q, &hist, intList);
                  computeBodyProperties_sincos_kkk(cmd, gd, q,
                                                   ipmax[cat1]-ipmin+1, &hist);
              } // ! end i loop
          } // ! scanoption smooth-pivot
#endif // ! ADDPIVOTNEIGHBOURS
#endif // ! NMultipoles

          ip = p - btable[cat1] + 1;
          //B kappa Avg Rmin
          icountNbRminthread += NbRmin(p);
          icountNbRminOverlapthread += NbRminOverlap(p);
#ifdef DEBUG
#ifdef ADDPIVOTNEIGHBOURS
          fprintf(outpivots,"%ld \t%ld \t%ld \t%ld \t\t%g\n",
                  ip, NbRmin(p), NbRminOverlap(p), intList,
                  KappaRmin(p)/NbRmin(p));
#else
          fprintf(outpivots,"%ld \t%ld \t%ld \t\t%g\n",
                  ip, NbRmin(p), NbRminOverlap(p),
                  KappaRmin(p)/NbRmin(p));
#endif
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
        for (m=1; m<=cmd->mChebyshev+1; m++) {
            ADDM_ext(gd->histZetaMcos[m],gd->histZetaMcos[m],
                     hist.histZetaMthreadcos[m],cmd->sizeHistN);
            ADDM_ext(gd->histZetaMsin[m],gd->histZetaMsin[m],
                     hist.histZetaMthreadsin[m],cmd->sizeHistN);
            ADDM_ext(gd->histZetaMsincos[m],gd->histZetaMsincos[m],
                     hist.histZetaMthreadsincos[m],cmd->sizeHistN);
            ADDM_ext(gd->histZetaMcossin[m],gd->histZetaMcossin[m],
                     hist.histZetaMthreadcossin[m],cmd->sizeHistN);
        }
        gd->nbbcalc += hist.nbbcalcthread;
        gd->nbccalc += hist.nbccalcthread;
#ifdef NMultipoles
        for (m=1; m<=cmd->mChebyshev+1; m++) {
            ADDM_ext(gd->NhistZetaMcos[m],gd->NhistZetaMcos[m],
                     histN.histZetaMthreadcos[m],cmd->sizeHistN);
            ADDM_ext(gd->NhistZetaMsin[m],gd->NhistZetaMsin[m],
                     histN.histZetaMthreadsin[m],cmd->sizeHistN);
            ADDM_ext(gd->NhistZetaMsincos[m],gd->NhistZetaMsincos[m],
                     histN.histZetaMthreadsincos[m],cmd->sizeHistN);
            ADDM_ext(gd->NhistZetaMcossin[m],gd->NhistZetaMcossin[m],
                     histN.histZetaMthreadcossin[m],cmd->sizeHistN);
        }
#endif
        //B kappa Avg Rmin
        ipfalse += ipfalsethreads;
        icountNbRmin += icountNbRminthread;
        icountNbRminOverlap += icountNbRminOverlapthread;
        //E
    } // ! critical

#ifdef NMultipoles
    search_free_sincos_omp_kkk_N(cmd, gd, &histN);  // free memory
#endif
    search_free_sincos_omp_kkk(cmd, gd, &hist);     // free memory
  } // end pragma omp parallel

    if (cmd->verbose >= VERBOSENORMALINFO)
        verb_print(cmd->verbose, "\n\n");             // end of completed pivot

    //B kappa Avg Rmin
    real xi, den, num;
    int mm;
    if (scanopt(cmd->options, "smooth-pivot")) {
        num = (real)nbody[cat1];
        den = (real)(nbody[cat1]-ipfalse);
#ifdef NONORMHIST
        xi = 1.0;
#else
#ifdef ADDPIVOTNEIGHBOURS
        xi = 1.0;
#else
        xi = num/den;
#endif
#endif
        if (cmd->verbose>=VERBOSENORMALINFO)
            verb_print(cmd->verbose,
                       "octree-kkk-omp: p falses found = %ld and %e %e %e\n",
                       ipfalse, num, den, xi);
        for (mm=1; mm<=cmd->mChebyshev+1; mm++) {
            MULMS_ext(gd->histZetaMcos[mm],
                      gd->histZetaMcos[mm],xi,cmd->sizeHistN);
            MULMS_ext(gd->histZetaMsin[mm],
                      gd->histZetaMsin[mm],xi,cmd->sizeHistN);
            MULMS_ext(gd->histZetaMsincos[mm],
                      gd->histZetaMsincos[mm],xi,cmd->sizeHistN);
            MULMS_ext(gd->histZetaMcossin[mm],
                      gd->histZetaMcossin[mm],xi,cmd->sizeHistN);
        }
    }
#ifdef NMultipoles
    if (scanopt(cmd->options, "smooth-pivot")) {
        for (mm=1; mm<=cmd->mChebyshev+1; mm++) {
            MULMS_ext(gd->NhistZetaMcos[mm],
                      gd->NhistZetaMcos[mm],xi,cmd->sizeHistN);
            MULMS_ext(gd->NhistZetaMsin[mm],
                      gd->NhistZetaMsin[mm],xi,cmd->sizeHistN);
            MULMS_ext(gd->NhistZetaMsincos[mm],
                      gd->NhistZetaMsincos[mm],xi,cmd->sizeHistN);
            MULMS_ext(gd->NhistZetaMcossin[mm],
                      gd->NhistZetaMcossin[mm],xi,cmd->sizeHistN);
        }
    }
#endif
    //E

    if (scanopt(cmd->options, "smooth-pivot")) {
        if (cmd->verbose>=VERBOSENORMALINFO) {
            verb_print(cmd->verbose,
                       "octree-kkk-omp: p falses found = %ld\n",ipfalse);
            //B kappa Avg Rmin
            verb_print(cmd->verbose,
                       "octree-kkk-omp: count NbRmin found = %ld\n",
                       icountNbRmin);
            verb_print(cmd->verbose,
                       "octree-kkk-omp: count overlap found = %ld\n",
                       icountNbRminOverlap);
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
            verb_print(cmd->verbose, "octree-kkk-omp: p falses found = %ld\n",
                       ifalsecount);
            verb_print(cmd->verbose, "octree-kkk-omp: p true found = %ld\n",itruecount);
            verb_print(cmd->verbose, "octree-kkk-omp: total = %ld\n",
                       itruecount+ifalsecount);
        }
        //E
    }

#ifdef DEBUG
    fclose(outpivots);                              // Close file to debug pivots
#endif

    gd->cpusearch = CPUTIME - cpustart;
    verb_print(cmd->verbose, "\nGoing out: CPU time = %lf %s\n",
               CPUTIME-cpustart, PRNUNITOFTIMEUSED);

    return SUCCESS;
}

local void normal_walktree_sincos(struct  cmdline_data* cmd, 
                                  struct  global_data* gd,
                                  bodyptr *btable, int cat2,
                                  bodyptr p, nodeptr q, real qsize,
                                  gdhistptr_sincos_omp_kkk hist,
                                  int *nbList, int *intList)
{
    nodeptr l;
    real dr1;
    vector dr;

    if (Update(p)) {
        if ( ((nodeptr) p) != q ) {
            if (Type(q) == CELL) {
                if (!reject_cell(cmd, gd, (nodeptr)p, q, qsize)) {
                    if (!scanopt(cmd->options, "no-one-ball")) {
                        accept_body(cmd, gd,
                                    p, (nodeptr)q, &dr1, dr);
                        if ( (Radius(p)+Radius(q))/(dr1) < gd->deltaR)
                            sumnode_sincos_cell(cmd, gd, btable, cat2, p,
                                                ((cellptr) q),
                                                ((cellptr) q+1), hist,
                                                nbList, intList);
                        else
                            for (l = More(q); l != Next(q); l = Next(l))
                                normal_walktree_sincos(cmd, gd, btable, cat2,
                                                       p,l,qsize/2, hist,
                                                       nbList, intList);
                    } else {
                        for (l = More(q); l != Next(q); l = Next(l))
                            normal_walktree_sincos(cmd, gd, btable, cat2,
                                                   p,l,qsize/2, hist,
                                                   nbList, intList);
                    }
                }
            } else { // ! Type(q) == CELL
                sumnode_sincos(cmd, gd, btable, cat2,
                               p, ((cellptr)q), ((cellptr)q+1), hist,
                               nbList, intList);
            } // ! Type(q) == CELL
        } // ! p != q
    }
}

#ifdef ADDPIVOTNEIGHBOURS
//B kappa Avg Rmin
local void sumnode_nblist_omp(struct cmdline_data* cmd,
                              struct  global_data* gd,
                              bodyptr *btable,
                              INTEGER ipmin, INTEGER *ipmax,
                              int cat1, int cat2,
                              bodyptr p,
                              gdhistptr_sincos_omp_kkk hist, int intList)
{
    bodyptr q;
    real dr1;
    vector dr;
    int i;
    int n;
    real xi;
    REAL cosphi,sinphi;

    REAL s, sy;
    vector pr0;

    for (i = 0; i < intList; i++) {
        q = btable[cat2] + activeInt[i];
        accept_body(cmd, gd, p, (nodeptr)q, &dr1, dr);
        if(dr1>cmd->rminHist) {
            if (cmd->rminHist==0)
                n = (int)(cmd->logHistBinsPD*(rlog10(dr1) - rlog10(cmd->rangeN))
                                  + cmd->sizeHistN) + 1;
            else
                n = (int)(rlog10(dr1/cmd->rminHist) * gd->i_deltaR) + 1;
            if (n<=cmd->sizeHistN && n>=1) {
                hist->histNNSubthread[n] = hist->histNNSubthread[n] + 1.0;
                xi = Kappa(q);
                DOTVP(s, dr, hist->dr0);
                cosphi = s/(dr1*hist->drpq);
                CROSSVP(pr0,hist->dr0,Pos(p));
                DOTVP(sy, dr, pr0);
                sinphi = rsqrt(1.0 - rsqr(cosphi));
                if (sy < 0) sinphi *= -1.0;
                if (rabs(cosphi)>1.0)
                    verb_log_print(cmd->verbose, gd->outlog,
                        "sumenode: Warning!... cossphi must be in (-1,1): %g\n",
                        cosphi);
                CHEBYSHEVTUOMP;
                hist->nbccalcthread += 1;
            } // ! 1 < n < sizeHistN
        } // ! dr1 > rminHist
    } // ! end loop i
}
//E
#endif

local void sumnode_sincos(struct  cmdline_data* cmd,
                          struct  global_data* gd,
                          bodyptr *btable, int cat2,
                          bodyptr p, cellptr start, cellptr finish,
                          gdhistptr_sincos_omp_kkk hist,
                          int *nbList, int *intList)
{
    cellptr q;
    real dr1;
    vector dr;
    int n;
    real xi;
    REAL cosphi,sinphi;
    int iq;

    q = start;
    if (accept_body(cmd, gd, p, (nodeptr)q, &dr1, dr)) {
        //B kappa Avg Rmin
        if (scanopt(cmd->options, "smooth-pivot")) {
            if (dr1<=gd->rsmooth[0]) {
                if (Update(q)==TRUE) {
#ifdef ADDPIVOTNEIGHBOURS
                    iq = (bodyptr)q-btable[cat2];
                    activeNb[*nbList]=iq;
                    *nbList +=1;
                    if (*nbList > actlenNb)
                        error("nbList: too many neighbors, %d %d\n",
                              *nbList, actlenNb);
#endif
                    Update(q) = FALSE;
                    NbRmin(p) += 1;
                    KappaRmin(p) += Weight(q)*Kappa(q);
                } else {
                    NbRminOverlap(p) += 1;
                }
            }
        }
        //E
        if(dr1>cmd->rminHist) {
#ifdef ADDPIVOTNEIGHBOURS
            if (scanopt(cmd->options, "smooth-pivot")) {
                iq = (bodyptr)q-btable[cat2];
                activeInt[*intList]=iq;
                *intList +=1;
                if (*intList > actlenInt)
                    error("intList: too many neighbors\n");
            }
#endif
            if (cmd->rminHist==0)
                n = (int)(cmd->logHistBinsPD*(rlog10(dr1)
                    - rlog10(cmd->rangeN)) + cmd->sizeHistN) + 1;
            else
                n = (int)(rlog10(dr1/cmd->rminHist) * gd->i_deltaR) + 1;
            if (n<=cmd->sizeHistN && n>=1) {
                hist->histNNSubthread[n] = hist->histNNSubthread[n] + 1.;
                xi = Weight(q)*Kappa(q);
#ifdef POLARAXIS
                real a, c, c2;
                vector vc;
                DOTPSUBV(c2, vc, Pos(q), hist->q0);
#ifdef NOLIMBER
                a = 2.0*rasin(dr1/2.0);
#else
                a = dr1;
#endif
                real cosc;
                cosc = Pos(q)[2];
                cosphi = (cosc - (1.0-0.5*rsqr(a))*hist->cosb)
                        /(a*hist->sinb);
                if (rabs(cosphi) <= 1.0)
                    sinphi = rsqrt(1.0 - rsqr(cosphi));
                else
                    sinphi = 0.0;
                if (!crossVecProdSign(Pos(p), hist->q0, Pos(q)))
                    sinphi *= -1.0;
#else // ! POLARAXIS
                REAL s, sy;
                vector pr0;
                DOTVP(s, dr, hist->dr0);
                cosphi = s/(dr1*hist->drpq);
                CROSSVP(pr0,hist->dr0,Pos(p));
                DOTVP(sy, dr, pr0);
                sinphi = rsqrt(1.0 - rsqr(cosphi));
                if (sy < 0) sinphi *= -1.0;
                if (rabs(cosphi)>1.0)
                    verb_log_print(cmd->verbose, gd->outlog,
                    "sumenode: Warning!... cossphi must be in (-1,1): %g\n",
                                        cosphi);
#endif // ! POLARAXIS
                CHEBYSHEVTUOMP;
                hist->nbbcalcthread += 1;
            } // ! 1 < n < sizeHistN
        } // ! dr1>cmd->rminHist
    } // ! accept_body
}

local void sumnode_sincos_cell(struct  cmdline_data* cmd,
                               struct  global_data* gd,
                               bodyptr *btable, int cat2,
                               bodyptr p, cellptr start, cellptr finish,
                               gdhistptr_sincos_omp_kkk hist,
                               int *nbList, int *intList)
{
    cellptr q;
    real dr1;
    vector dr;
    int n;
    real xi;
    REAL cosphi,sinphi;

    q = start;
    if (accept_body(cmd, gd, p, (nodeptr)q, &dr1, dr)) {
        if(dr1>cmd->rminHist) {
#ifdef ADDPIVOTNEIGHBOURS
            INTEGER iq;
            iq = (bodyptr)q-btable[cat2];
            activeInt[*intList]=iq;
            *intList +=1;
            if (*intList > actlenInt)
                error("intList: too many neighbors\n");
#endif
            if (cmd->rminHist==0)
                n = (int)(cmd->logHistBinsPD*(rlog10(dr1)
                    - rlog10(cmd->rangeN)) + cmd->sizeHistN) + 1;
            else
                n = (int)(rlog10(dr1/cmd->rminHist) * gd->i_deltaR) + 1;
            if (n<=cmd->sizeHistN && n>=1) {
                hist->histNNSubthread[n] = hist->histNNSubthread[n] + 1.0;
                xi = Kappa(q);
#ifdef POLARAXIS
                real a, c, c2;
                vector vc;
                DOTPSUBV(c2, vc, Pos(q), hist->q0);
#ifdef NOLIMBER
                a = 2.0*rasin(dr1/2.0);
#else
                a = dr1;
#endif
                real cosc;
                cosc = Pos(q)[2];
                cosphi = (cosc - (1.0-0.5*rsqr(a))*hist->cosb)
                        /(a*hist->sinb);
                if (rabs(cosphi) <= 1.0)
                    sinphi = rsqrt(1.0 - rsqr(cosphi));
                else
                    sinphi = 0.0;
                if (!crossVecProdSign(Pos(p), hist->q0, Pos(q)))
                    sinphi *= -1.0;
#else
                REAL s, sy;
                vector pr0;
                DOTVP(s, dr, hist->dr0);
                cosphi = s/(dr1*hist->drpq);
                CROSSVP(pr0,hist->dr0,Pos(p));
                DOTVP(sy, dr, pr0);
                sinphi = rsqrt(1.0 - rsqr(cosphi));;
                if (sy < 0) sinphi *= -1.0;
                if (rabs(cosphi)>1.0)
                    verb_log_print(cmd->verbose, gd->outlog,
                        "sumenode: Warning!... cossphi must be in (-1,1): %g\n",
                                           cosphi);
#endif
                CHEBYSHEVTUOMP;
                hist->nbccalcthread += 1;
            } // ! 1 < n < sizeHistN
        } // ! dr1 > rminHist
    } // ! accept_body
}

#ifdef NMultipoles
local void normal_walktree_sincos_N(struct  cmdline_data* cmd,
                                    struct  global_data* gd,
                                    bodyptr *btable, int cat2,
                                    bodyptr p, nodeptr q, real qsize,
                                    gdhistptr_sincos_omp_kkk hist,
                                    gdhistptr_sincos_omp_kkk_N histN,
                                    int *nbList, int *intList)
{
    nodeptr l;
    real dr1;
    vector dr;

    if (Update(p)) {
        if ( ((nodeptr) p) != q ) {
            if (Type(q) == CELL) {
                if (!reject_cell(cmd, gd, (nodeptr)p, q, qsize)) {
                    if (!scanopt(cmd->options, "no-one-ball")) {
                        accept_body(cmd, gd, p, (nodeptr)q, &dr1, dr);
                        if ( (Radius(p)+Radius(q))/(dr1) < gd->deltaR)
                            sumnode_sincos_cell_N(cmd, gd, btable, cat2,
                                                  p, ((cellptr) q),
                                                  ((cellptr) q+1), hist, histN,
                                                  nbList, intList);
                        else
                            for (l = More(q); l != Next(q); l = Next(l))
                                normal_walktree_sincos_N(cmd, gd, btable, cat2,
                                                         p,l,qsize/2,
                                                         hist, histN,
                                                         nbList, intList);
                    } else {
                        for (l = More(q); l != Next(q); l = Next(l))
                            normal_walktree_sincos_N(cmd, gd, btable, cat2,
                                                     p,l,qsize/2,
                                                     hist, histN,
                                                     nbList, intList);
                    }
                }
            } else { // ! Type(q) == CELL
                sumnode_sincos_N(cmd, gd, btable, cat2,
                                 p, ((cellptr)q), ((cellptr)q+1),
                                 hist, histN, nbList, intList);
            } // ! Type(q) == CELL
        } // ! p != q
    } // ! Update
}

local void sumnode_sincos_N(struct  cmdline_data* cmd,
                            struct  global_data* gd,
                            bodyptr *btable, int cat2,
                            bodyptr p, cellptr start, cellptr finish,
                            gdhistptr_sincos_omp_kkk hist,
                            gdhistptr_sincos_omp_kkk_N histN,
                            int *nbList, int *intList)
{
    cellptr q;
    real dr1;
    vector dr;
    int n;
    real xi;
    real xiN;
    REAL cosphi,sinphi;
    int iq;

    q = start;
    if (accept_body(cmd, gd, p, (nodeptr)q, &dr1, dr)) {
        //B kappa Avg Rmin
        if (scanopt(cmd->options, "smooth-pivot")) {
            if (dr1<=gd->rsmooth[0]) {
                if (Update(q)==TRUE) {
#ifdef ADDPIVOTNEIGHBOURS
                    iq = (bodyptr)q-btable[cat2];
                    activeNb[*nbList]=iq;
                    *nbList +=1;
                    if (*nbList > actlenNb)
                        error("nbList: too many neighbors, %d %d\n",
                              *nbList, actlenNb);
#endif
                    Update(q) = FALSE;
                    NbRmin(p) += 1;
                    KappaRmin(p) += Kappa(q);
                } else {
                    NbRminOverlap(p) += 1;
                }
            }
        }
        //E
        if(dr1>cmd->rminHist) {
#ifdef ADDPIVOTNEIGHBOURS
            iq = (bodyptr)q-btable[cat2];
            activeInt[*intList]=iq;
            *intList +=1;
            if (*intList > actlenInt)
                error("intList: too many neighbors\n");
#endif
            if (cmd->rminHist==0)
                n = (int)(cmd->logHistBinsPD*(rlog10(dr1)
                    - rlog10(cmd->rangeN))
                    + cmd->sizeHistN) + 1;
            else
                n = (int)(rlog10(dr1/cmd->rminHist) * gd->i_deltaR) + 1;
            if (n<=cmd->sizeHistN && n>=1) {
                hist->histNNSubthread[n] = hist->histNNSubthread[n] + 1.;
                histN->histNNSubthread[n] = histN->histNNSubthread[n] + 1.;
                xi = Kappa(q);
                xiN = 1.0;
#ifdef POLARAXIS
                real a, c, c2;
                vector vc;
                DOTPSUBV(c2, vc, Pos(q), hist->q0);
#ifdef NOLIMBER
                a = 2.0*rasin(dr1/2.0);
#else
                a = dr1;
#endif
                real cosc;
                cosc = Pos(q)[2];
                cosphi = (cosc - (1.0-0.5*rsqr(a))*hist->cosb)
                        /(a*hist->sinb);
                if (rabs(cosphi) <= 1.0)
                    sinphi = rsqrt(1.0 - rsqr(cosphi));
                else
                    sinphi = 0.0;
                if (!crossVecProdSign(Pos(p), hist->q0, Pos(q)))
                    sinphi *= -1.0;
#else // ! POLARAXIS
                REAL s, sy;
                vector pr0;
                DOTVP(s, dr, hist->dr0);
                cosphi = s/(dr1*hist->drpq);
                CROSSVP(pr0,hist->dr0,Pos(p));
                DOTVP(sy, dr, pr0);
                sinphi = rsqrt(1.0 - rsqr(cosphi));
                if (sy < 0) sinphi *= -1.0;
                if (rabs(cosphi)>1.0)
                    verb_log_print(cmd->verbose, gd->outlog,
                    "sumenode: Warning!... cossphi must be in (-1,1): %g\n",
                                    cosphi);
#endif // ! POLARAXIS
                CHEBYSHEVTUOMPN;
                CHEBYSHEVTUOMP;
                hist->nbbcalcthread += 1;
            }
        }
    } // ! accept_body
}

local void sumnode_sincos_cell_N(struct  cmdline_data* cmd,
                                 struct  global_data* gd,
                                 bodyptr *btable, int cat2,
                                 bodyptr p, cellptr start, cellptr finish,
                                 gdhistptr_sincos_omp_kkk hist,
                                 gdhistptr_sincos_omp_kkk_N histN,
                                 int *nbList, int *intList)
{
    cellptr q;
    real dr1;
    vector dr;
    int n;
    real xi;
    real xiN;
    REAL cosphi,sinphi;

    q = start;
    if (accept_body(cmd, gd, p, (nodeptr)q, &dr1, dr)) {
        if(dr1>cmd->rminHist) {
#ifdef ADDPIVOTNEIGHBOURS
            INTEGER iq;
            iq = (bodyptr)q-btable[cat2];
            activeInt[*intList]=iq;
            *intList +=1;
            if (*intList > actlenInt)
                error("intList: too many neighbors\n");
#endif
            if (cmd->rminHist==0)
                n = (int)(cmd->logHistBinsPD*(rlog10(dr1) - rlog10(cmd->rangeN))
                                  + cmd->sizeHistN) + 1;
            else
                n = (int)(rlog10(dr1/cmd->rminHist) * gd->i_deltaR) + 1;
            if (n<=cmd->sizeHistN && n>=1) {
                hist->histNNSubthread[n] = hist->histNNSubthread[n] + 1.0;
                histN->histNNSubthread[n] = histN->histNNSubthread[n] + 1.0;
                xi = Kappa(q);
                xiN = 1.0;
#ifdef POLARAXIS
                real a, c, c2;
                vector vc;
                DOTPSUBV(c2, vc, Pos(q), hist->q0);
#ifdef NOLIMBER
                a = 2.0*rasin(dr1/2.0);
#else
                a = dr1;
#endif
                real cosc;
                cosc = Pos(q)[2];
                cosphi = (cosc - (1.0-0.5*rsqr(a))*hist->cosb)
                        /(a*hist->sinb);
                if (rabs(cosphi) <= 1.0)
                    sinphi = rsqrt(1.0 - rsqr(cosphi));
                else
                    sinphi = 0.0;
                if (!crossVecProdSign(Pos(p), hist->q0, Pos(q)))
                    sinphi *= -1.0;
#else // ! POLARAXIS
                REAL s, sy;
                vector pr0;
                DOTVP(s, dr, hist->dr0);
                cosphi = s/(dr1*hist->drpq);
                CROSSVP(pr0,hist->dr0,Pos(p));
                DOTVP(sy, dr, pr0);
                sinphi = rsqrt(1.0 - rsqr(cosphi));;
                if (sy < 0) sinphi *= -1.0;
                if (rabs(cosphi)>1.0)
                    verb_log_print(cmd->verbose, gd->outlog,
                        "sumenode: Warning!... cossphi must be in (-1,1): %g\n",
                                           cosphi);
#endif // ! POLARAXIS
                CHEBYSHEVTUOMPN;
                CHEBYSHEVTUOMP;
                hist->nbccalcthread += 1;
            }
        }
    } // ! accept_body
}
#endif


//B Routines as in cballsutils

local int search_init_gd_sincos_omp_kkk(struct  cmdline_data* cmd,
                                        struct  global_data* gd)
{
    int n;
    int m;

    for (n = 1; n <= cmd->sizeHistN; n++)
        gd->histNNSub[n] = 0.0;
    for (m = 1; m <= cmd->mChebyshev+1; m++) {
        CLRM_ext(gd->histZetaMcos[m], cmd->sizeHistN);
        CLRM_ext(gd->histZetaMsin[m], cmd->sizeHistN);
        CLRM_ext(gd->histZetaMsincos[m], cmd->sizeHistN);
        CLRM_ext(gd->histZetaMcossin[m], cmd->sizeHistN);
        gd->histXi[m][n] = 0.0;
    }
    gd->nbbcalc = gd->nbccalc = gd->ncccalc = 0;

    return SUCCESS;
}

local int search_init_sincos_omp_kkk(struct  cmdline_data* cmd,
                                     struct  global_data* gd,
                                     gdhistptr_sincos_omp_kkk hist)
{
    int n;
    int m;

    hist->ChebsT = dvector(1,cmd->mChebyshev+1);
    hist->ChebsU = dvector(1,cmd->mChebyshev+1);
    hist->histNNSubthread = dvector(1,cmd->sizeHistN);
    hist->histXithreadcos = dmatrix(1,cmd->mChebyshev+1,1,cmd->sizeHistN);
    hist->histXithreadsin = dmatrix(1,cmd->mChebyshev+1,1,cmd->sizeHistN);
    hist->histZetaMthreadcos =
            dmatrix3D(1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    hist->histZetaMthreadsin =
            dmatrix3D(1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    hist->histZetaMthreadsincos =
            dmatrix3D(1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    hist->histZetaMthreadcossin =
            dmatrix3D(1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    hist->xiOUTVPcos = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
    hist->xiOUTVPsin = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
    hist->xiOUTVPsincos = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
    hist->xiOUTVPcossin = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
    hist->histZetaMtmpcos = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
    hist->histZetaMtmpsin = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
    hist->histZetaMtmpsincos = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
    hist->histZetaMtmpcossin = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);

    for (n = 1; n <= cmd->sizeHistN; n++)
        hist->histNNSubthread[n] = 0.0;
    for (m = 1; m <= cmd->mChebyshev+1; m++) {
        CLRM_ext(hist->histZetaMthreadcos[m], cmd->sizeHistN);
        CLRM_ext(hist->histZetaMthreadsin[m], cmd->sizeHistN);
        CLRM_ext(hist->histZetaMthreadsincos[m], cmd->sizeHistN);
        CLRM_ext(hist->histZetaMthreadcossin[m], cmd->sizeHistN);
    }

    hist->nbbcalcthread = 0;
    hist->nbccalcthread = 0;

    return SUCCESS;
}

local int search_free_sincos_omp_kkk(struct  cmdline_data* cmd,
                                         struct  global_data* gd,
                                      gdhistptr_sincos_omp_kkk hist)
{
    free_dmatrix(hist->histZetaMtmpcossin,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix(hist->histZetaMtmpsincos,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix(hist->histZetaMtmpsin,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix(hist->histZetaMtmpcos,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix(hist->xiOUTVPcossin,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix(hist->xiOUTVPsincos,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix(hist->xiOUTVPsin,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix(hist->xiOUTVPcos,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix3D(hist->histZetaMthreadcossin,
                   1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix3D(hist->histZetaMthreadsincos,
                   1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix3D(hist->histZetaMthreadsin,
                   1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix3D(hist->histZetaMthreadcos,
                   1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix(hist->histXithreadsin,1,cmd->mChebyshev+1,1,cmd->sizeHistN);
    free_dmatrix(hist->histXithreadcos,1,cmd->mChebyshev+1,1,cmd->sizeHistN);
    free_dvector(hist->histNNSubthread,1,cmd->sizeHistN);
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

// check Weight factor... must be an average of Weights
#ifdef NONORMHIST
    xi = Weight(p)*Kappa(p);
    if (scanopt(cmd->options, "smooth-pivot")) {
        xi = KappaRmin(p)/NbRmin(p);
    }
#else
#ifdef ADDPIVOTNEIGHBOURS
    xi = Weight(p)*Kappa(p)/nbody;
#else
    xi = Weight(p)*Kappa(p)/nbody;
    //B kappa Avg Rmin
    if (scanopt(cmd->options, "smooth-pivot")) {
        xi = (KappaRmin(p)/NbRmin(p))/nbody;
    }
    //E
#endif
#endif

    for (m=1; m<=cmd->mChebyshev+1; m++)
#ifdef NONORMHIST
        for (n=1; n<=cmd->sizeHistN; n++) {
            hist->histXithreadcos[m][n] /= 1.0;
            hist->histXithreadsin[m][n] /= 1.0;
        }
#else
        for (n=1; n<=cmd->sizeHistN; n++) {
            hist->histXithreadcos[m][n] /= MAX(hist->histNNSubthread[n],1.0);
            hist->histXithreadsin[m][n] /= MAX(hist->histNNSubthread[n],1.0);
        }
#endif
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
        MULMS_ext(hist->histZetaMtmpsincos,hist->xiOUTVPsincos,xi,cmd->sizeHistN);
        MULMS_ext(hist->histZetaMtmpcossin,hist->xiOUTVPcossin,xi,cmd->sizeHistN);
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

    return SUCCESS;
}

#ifdef NMultipoles
local int search_init_gd_sincos_omp_kkk_N(struct  cmdline_data* cmd,
                                        struct  global_data* gd)
{
    int n;
    int m;

//    verb_print_debug(1, "\nAqui voy (0)\n");
    for (n = 1; n <= cmd->sizeHistN; n++)
        gd->NhistNNSub[n] = 0.0;
//    verb_print_debug(1, "\nAqui voy (1)\n");
    for (m = 1; m <= cmd->mChebyshev+1; m++) {
        CLRM_ext(gd->NhistZetaMcos[m], cmd->sizeHistN);
        CLRM_ext(gd->NhistZetaMsin[m], cmd->sizeHistN);
        CLRM_ext(gd->NhistZetaMsincos[m], cmd->sizeHistN);
        CLRM_ext(gd->NhistZetaMcossin[m], cmd->sizeHistN);
        gd->NhistXi[m][n] = 0.0;
    }
//    verb_print_debug(1, "\nAqui voy (2)\n");

    return SUCCESS;
}

local int search_init_sincos_omp_kkk_N(struct  cmdline_data* cmd,
                                     struct  global_data* gd,
                                     gdhistptr_sincos_omp_kkk_N hist)
{
    int n;
    int m;

    hist->ChebsT = dvector(1,cmd->mChebyshev+1);
    hist->ChebsU = dvector(1,cmd->mChebyshev+1);
    hist->histNNSubthread = dvector(1,cmd->sizeHistN);
    hist->histXithreadcos = dmatrix(1,cmd->mChebyshev+1,1,cmd->sizeHistN);
    hist->histXithreadsin = dmatrix(1,cmd->mChebyshev+1,1,cmd->sizeHistN);
    hist->histZetaMthreadcos =
            dmatrix3D(1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    hist->histZetaMthreadsin =
            dmatrix3D(1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    hist->histZetaMthreadsincos =
            dmatrix3D(1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    hist->histZetaMthreadcossin =
            dmatrix3D(1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    hist->xiOUTVPcos = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
    hist->xiOUTVPsin = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
    hist->xiOUTVPsincos = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
    hist->xiOUTVPcossin = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
    hist->histZetaMtmpcos = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
    hist->histZetaMtmpsin = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
    hist->histZetaMtmpsincos = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
    hist->histZetaMtmpcossin = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);

    for (n = 1; n <= cmd->sizeHistN; n++)
        hist->histNNSubthread[n] = 0.0;
    for (m = 1; m <= cmd->mChebyshev+1; m++) {
        CLRM_ext(hist->histZetaMthreadcos[m], cmd->sizeHistN);
        CLRM_ext(hist->histZetaMthreadsin[m], cmd->sizeHistN);
        CLRM_ext(hist->histZetaMthreadsincos[m], cmd->sizeHistN);
        CLRM_ext(hist->histZetaMthreadcossin[m], cmd->sizeHistN);
    }

    return SUCCESS;
}

local int search_free_sincos_omp_kkk_N(struct  cmdline_data* cmd,
                                         struct  global_data* gd,
                                      gdhistptr_sincos_omp_kkk_N hist)
{
    free_dmatrix(hist->histZetaMtmpcossin,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix(hist->histZetaMtmpsincos,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix(hist->histZetaMtmpsin,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix(hist->histZetaMtmpcos,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix(hist->xiOUTVPcossin,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix(hist->xiOUTVPsincos,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix(hist->xiOUTVPsin,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix(hist->xiOUTVPcos,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix3D(hist->histZetaMthreadcossin,
                   1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix3D(hist->histZetaMthreadsincos,
                   1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix3D(hist->histZetaMthreadsin,
                   1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix3D(hist->histZetaMthreadcos,
                   1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix(hist->histXithreadsin,1,cmd->mChebyshev+1,1,cmd->sizeHistN);
    free_dmatrix(hist->histXithreadcos,1,cmd->mChebyshev+1,1,cmd->sizeHistN);
    free_dvector(hist->histNNSubthread,1,cmd->sizeHistN);
    free_dvector(hist->ChebsU,1,cmd->mChebyshev+1);
    free_dvector(hist->ChebsT,1,cmd->mChebyshev+1);

    return SUCCESS;
}

local int computeBodyProperties_sincos_kkk_N(struct  cmdline_data* cmd,
                                            struct  global_data* gd,
                                            bodyptr p, int nbody,
                                            gdhistptr_sincos_omp_kkk_N hist)
{
    int n;
    int m;
    real xi;

#ifdef NONORMHIST
    xi = 1.0;
#else
    xi = 1.0/nbody;
#endif

    for (m=1; m<=cmd->mChebyshev+1; m++)
#ifdef NONORMHIST
        for (n=1; n<=cmd->sizeHistN; n++) {
            hist->histXithreadcos[m][n] /= 1.0;
            hist->histXithreadsin[m][n] /= 1.0;
        }
#else
        for (n=1; n<=cmd->sizeHistN; n++) {
            hist->histXithreadcos[m][n] /= MAX(hist->histNNSubthread[n],1.0);
            hist->histXithreadsin[m][n] /= MAX(hist->histNNSubthread[n],1.0);
        }
#endif
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
        MULMS_ext(hist->histZetaMtmpsincos,hist->xiOUTVPsincos,xi,cmd->sizeHistN);
        MULMS_ext(hist->histZetaMtmpcossin,hist->xiOUTVPcossin,xi,cmd->sizeHistN);
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

    return SUCCESS;
}
#endif

//E Routines as in cballsutils

local int print_info(struct cmdline_data* cmd,
                                  struct  global_data* gd)
{
    verb_print(cmd->verbose,
               "searchcalc: Using octree-kkk-omp... \n");

    if (cmd->usePeriodic==TRUE)
        error("CheckParameters: can´t have periodic boundaries and OCTREEKKKOMP definition (usePeriodic=%d)\nSet usePeriodic=false\n",
            cmd->usePeriodic);
    if (cmd->useLogHist==FALSE)
        error("CheckParameters: can´t have normal scale hist and OCTREEKKKOMP definition (useLogHist=%d)\nSet useLogHist=true\n",
            cmd->useLogHist);
    if (cmd->computeTPCF==FALSE)
        error("CheckParameters: can´t have computeTPCF=false and OCTREEKKKOMP definition (computeTPCF=%d)\nSet computeTPCF=true\n",
            cmd->computeTPCF);
#if NDIM == 2
error("CheckParameters: OCTREEKKKOMP definition works only in a 3D unit sphere")
#endif
#ifdef NMultipoles
    verb_print(cmd->verbose, "with NMultipoles... \n");
#else
    verb_print(cmd->verbose, "without NMultipoles... \n");
#endif
#ifdef NONORMHIST
    verb_print(cmd->verbose, "with NONORMHIST... \n");
    if (scanopt(cmd->options, "no-normalize-HistZeta"))
        verb_print(cmd->verbose, "with option no-normalize-HistZeta...\n");
#else
    verb_print(cmd->verbose, "without NONORMHIST... \n");
#endif
#ifdef POLARAXIS
    verb_print(cmd->verbose, "with POLARAXIS... \n");
#endif
#ifdef ADDPIVOTNEIGHBOURS
    if (scanopt(cmd->options, "smooth-pivot"))
        verb_print(cmd->verbose, "with ADDPIVOTNEIGHBOURS... \n");
#endif
    if (scanopt(cmd->options, "no-one-ball"))
        verb_print(cmd->verbose, "with option no-one-ball... \n");
    if (scanopt(cmd->options, "smooth-pivot"))
        verb_print(cmd->verbose,
                   "with option smooth-pivot... rsmooth=%g\n",gd->rsmooth[0]);
    if (scanopt(cmd->options, "default-rsmooth"))
        verb_print(cmd->verbose, "with option default-rsmooth... \n");
    if (scanopt(cmd->options, "fix-rsmooth"))
        verb_print(cmd->verbose, "with option fix-rsmooth... \n");

    return SUCCESS;
}


#ifdef ADDPIVOTNEIGHBOURS
#undef FACTIVENB
#undef FACTIVEINT
#endif
