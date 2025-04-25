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

//B correction 2025-04-06
//B What if mChebyshev is less than 7... correct!
//#ifdef MANUALCHEBYSHEV
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
/*
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
*/

#ifdef NMultipoles
//#ifdef MANUALCHEBYSHEV
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
/*
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
*/

//B Macro for any posible value of mChebyshev
//  for recursivity needs that at least 3 multipoles be evaluated
#define CHEBYSHEVTUOMPNANY                                        \
{real xicosmphi,xisinmphi; int m;                                 \
    histN->ChebsT[1] = 1.0;                                       \
    xicosmphi = xiN * histN->ChebsT[1];                           \
    histN->histXithreadcos[1][n] += xicosmphi;                    \
    histN->ChebsT[2] = 1.0;                                       \
    xicosmphi = xiN * histN->ChebsT[2];                           \
    histN->histXithreadcos[2][n] += xicosmphi;                    \
    histN->ChebsT[3] = 1.0;                                       \
    xicosmphi = xiN * histN->ChebsT[3];                           \
    histN->histXithreadcos[3][n] += xicosmphi;                    \
    histN->ChebsU[1] = 0.0;                                       \
    xisinmphi = xiN * histN->ChebsU[1] * sinphi;                  \
    histN->histXithreadsin[1][n] += xisinmphi;                    \
    histN->ChebsU[2] = 0.0;                                       \
    xisinmphi = xiN * histN->ChebsU[2] * sinphi;                  \
    histN->histXithreadsin[2][n] += xisinmphi;                    \
    histN->ChebsU[3] = 0.0;                                       \
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
//E
#endif // ! NMultipoles

//E correction 2025-04-06


//B Define structures:
typedef struct {
    realptr histNNSub;
    //B TPCF
    real ***histZetaMcos;
    real ***histZetaMsin;
    real ***histZetaMsincos;
    // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
    real ***histZetaMcossin;
    //E
    //B Used to compute ZetaG using FFT
    real ***histZetaGmRe;
    real ***histZetaGmIm;
    real ***histXi3pcf;
    //E
    real ***histZetaM;
} gdl_sincos_omp_kkk, *gdlptr_sincos_omp_kkk;

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
    realptr histNNSub;
    //B TPCF
    real ***histZetaMcos;
    real ***histZetaMsin;
    real ***histZetaMsincos;
    // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
    real ***histZetaMcossin;
    //E
} gdl_sincos_omp_kkk_N, *gdlptr_sincos_omp_kkk_N;

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
//E

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
                                        struct  global_data* gd,
                                        gdlptr_sincos_omp_kkk);
local int search_free_gd_sincos_omp_kkk(struct  cmdline_data* cmd,
                                         struct  global_data* gd,
                                        gdlptr_sincos_omp_kkk);
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
                                          struct  global_data* gd,
                                          gdlptr_sincos_omp_kkk_N);
local int search_free_gd_sincos_omp_kkk_N(struct  cmdline_data* cmd,
                                         struct  global_data* gd,
                                          gdlptr_sincos_omp_kkk_N);
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

#ifdef NMultipoles
local int matrixClm(struct cmdline_data* cmd, struct  global_data* gd,
                    gdlptr_sincos_omp_kkk, gdlptr_sincos_omp_kkk_N,
                    int, int);
//B correction 2025-04-06
local int matrixClm_normalized_HistZeta(struct cmdline_data* cmd,
                                        struct  global_data* gd,
                    gdlptr_sincos_omp_kkk, gdlptr_sincos_omp_kkk_N,
                    int, int);
//E correction 2025-04-06
#endif


//B Saving histograms section: case KKKCORRELATION:
local int PrintHistrBins(struct  cmdline_data* cmd, struct  global_data* gd);
local int PrintHistZetaM_sincos(struct  cmdline_data* cmd,
                                struct  global_data* gd,
                                gdlptr_sincos_omp_kkk);
local int PrintHistZetaMm_sincos(struct  cmdline_data* cmd,
                               struct  global_data* gd,
                                 gdlptr_sincos_omp_kkk);
local int PrintHistZetaG(struct  cmdline_data* cmd,
                         struct  global_data* gd,
                         gdlptr_sincos_omp_kkk);
local int PrintHistZetaGm_sincos(struct  cmdline_data* cmd,
                                 struct  global_data* gd,
                                 gdlptr_sincos_omp_kkk);
local int PrintHistZetaMZetaGm_sincos(struct  cmdline_data* cmd,
                                      struct  global_data* gd,
                                      gdlptr_sincos_omp_kkk);

// NMultipoles precedes NONORMHIST?... Not necessarily.
//  Check consistency!!!
//  (NMultipoles, NONORMHIST):
//      1               1
//      1               0
//      0               1
//      0               0
#ifdef NMultipoles
local int PrintHistZetaM_sincos_N(struct  cmdline_data* cmd,
                                struct  global_data* gd,
                                  gdlptr_sincos_omp_kkk_N);
local int PrintHistZetaMm_sincos_N(struct  cmdline_data* cmd,
                               struct  global_data* gd,
                                   gdlptr_sincos_omp_kkk_N);
#ifdef NONORMHIST
// Saves matrix ZetaM for each m multipole
local int PrintHistZetaM_sincos_normalized(struct  cmdline_data* cmd,
                                           struct  global_data* gd,
                                           gdlptr_sincos_omp_kkk,
                                           gdlptr_sincos_omp_kkk_N);
// Saves matrix ZetaM for each m multipole at a set of theta2 angles
local int PrintHistZetaMm_sincos_normalized(struct  cmdline_data* cmd,
                                            struct  global_data* gd,
                                            gdlptr_sincos_omp_kkk,
                                            gdlptr_sincos_omp_kkk_N);
//B edge effects:
local int PrintHistZetaM_sincos_edge_effects(struct  cmdline_data*,
                                             struct  global_data*,
                                             gdlptr_sincos_omp_kkk,
                                             gdlptr_sincos_omp_kkk_N);
//B correction 2025-04-06
local int PrintHistZetaM_sincos_edge_effects_normalized_HistZeta(
                                            struct  cmdline_data*,
                                            struct  global_data*,
                                            gdlptr_sincos_omp_kkk,
                                            gdlptr_sincos_omp_kkk_N);
//E correction 2025-04-06

//B correction 2025-04-06
// Saves matrix ZetaM for each m multipole at a set of theta2 angles
//local int PrintHistZetaMm_sincos_edge_effects(struct  cmdline_data* cmd,
//                                              struct  global_data* gd,
//                                              gdlptr_sincos_omp_kkk,
//                                              gdlptr_sincos_omp_kkk_N);
//E correction 2025-04-06
//E
#endif
#endif // ! NMultipoles
//E

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
    gdl_sincos_omp_kkk gdl;
#ifdef NMultipoles
    gdl_sincos_omp_kkk_N gdlN;
#endif

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

    search_init_gd_sincos_omp_kkk(cmd, gd, &gdl);
#ifdef NMultipoles
    search_init_gd_sincos_omp_kkk_N(cmd, gd, &gdlN);
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
        //B Alloc memory for interaction lists
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

//
// Check that all posibilities are take in to account...
//
#ifdef DEBUG
//B In this segment there must be NMultipoles acting...
#ifdef ADDPIVOTNEIGHBOURS
#pragma omp parallel default(none)                                          \
    shared(cmd,gd,btable,nbody,roottable,outpivots,                         \
           actlenNb,activeNb,actlenInt,activeInt,                           \
           ipmin,ipmax,cat1,cat2,ipfalse,icountNbRmin,icountNbRminOverlap, gdl)
#else // ! ADDPIVOTNEIGHBOURS
#pragma omp parallel default(none)                                          \
    shared(cmd,gd,btable,nbody,roottable,outpivots,                         \
           ipmin,ipmax,cat1,cat2,ipfalse,icountNbRmin,icountNbRminOverlap, gdl)
#endif // ! ADDPIVOTNEIGHBOURS
//E
#else // ! DEBUG

#ifdef ADDPIVOTNEIGHBOURS
//B In this segment there must be NMultipoles acting...
#pragma omp parallel default(none)                                          \
    shared(cmd,gd,btable,nbody,roottable,                         \
           actlenNb,activeNb,actlenInt,activeInt,                           \
           ipmin,ipmax,cat1,cat2,ipfalse,icountNbRmin,icountNbRminOverlap, gdl)
//E
#else // ! ADDPIVOTNEIGHBOURS

#ifdef NMultipoles
#pragma omp parallel default(none)                                          \
    shared(cmd,gd,btable,nbody,roottable,                                   \
           ipmin,ipmax,cat1,cat2,ipfalse,icountNbRmin,icountNbRminOverlap,  \
           gdl, gdlN)
#else // ! NMultipoles
#pragma omp parallel default(none)                                          \
    shared(cmd,gd,btable,nbody,roottable,                         \
           ipmin,ipmax,cat1,cat2,ipfalse,icountNbRmin,icountNbRminOverlap,gdl)
#endif // ! NMultipoles

#endif // ! ADDPIVOTNEIGHBOURS

#endif // ! DEBUG
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
//B correction 2025-04-06
//#ifdef NMultipoles
                  for (n = 1; n <= cmd->sizeHistN; n++)
                      histN.histNNSubthread[n] = 0.0;
                  CLRM_ext_ext(histN.histXithreadcos,
                               cmd->mChebyshev+1, cmd->sizeHistN);
                  CLRM_ext_ext(histN.histXithreadsin,
                               cmd->mChebyshev+1, cmd->sizeHistN);
//B correction 2025-04-06
//#endif
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
//B correction 2025-04-06
/*#ifdef NMultipoles
                  for (n = 1; n <= cmd->sizeHistN; n++)
                      histN.histNNSubthread[n] = 0.0;
                  CLRM_ext_ext(histN.histXithreadcos,
                               cmd->mChebyshev+1, cmd->sizeHistN);
                  CLRM_ext_ext(histN.histXithreadsin,
                               cmd->mChebyshev+1, cmd->sizeHistN);
#endif */
//B correction 2025-04-06
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
            ADDM_ext(gdl.histZetaMcos[m],gdl.histZetaMcos[m],
                     hist.histZetaMthreadcos[m],cmd->sizeHistN);
            ADDM_ext(gdl.histZetaMsin[m],gdl.histZetaMsin[m],
                     hist.histZetaMthreadsin[m],cmd->sizeHistN);
            ADDM_ext(gdl.histZetaMsincos[m],gdl.histZetaMsincos[m],
                     hist.histZetaMthreadsincos[m],cmd->sizeHistN);
            ADDM_ext(gdl.histZetaMcossin[m],gdl.histZetaMcossin[m],
                     hist.histZetaMthreadcossin[m],cmd->sizeHistN);
        }
        gd->nbbcalc += hist.nbbcalcthread;
        gd->nbccalc += hist.nbccalcthread;
#ifdef NMultipoles
        for (m=1; m<=cmd->mChebyshev+1; m++) {
            ADDM_ext(gdlN.histZetaMcos[m],gdlN.histZetaMcos[m],
                     histN.histZetaMthreadcos[m],cmd->sizeHistN);
            ADDM_ext(gdlN.histZetaMsin[m],gdlN.histZetaMsin[m],
                     histN.histZetaMthreadsin[m],cmd->sizeHistN);
            ADDM_ext(gdlN.histZetaMsincos[m],gdlN.histZetaMsincos[m],
                     histN.histZetaMthreadsincos[m],cmd->sizeHistN);
            ADDM_ext(gdlN.histZetaMcossin[m],gdlN.histZetaMcossin[m],
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
#endif // ! NONORMHIST
        if (cmd->verbose>=VERBOSENORMALINFO)
            verb_print(cmd->verbose,
                       "octree-kkk-omp: p falses found = %ld and %e %e %e\n",
                       ipfalse, num, den, xi);
        for (mm=1; mm<=cmd->mChebyshev+1; mm++) {
            MULMS_ext(gdl.histZetaMcos[mm],
                      gdl.histZetaMcos[mm],xi,cmd->sizeHistN);
            MULMS_ext(gdl.histZetaMsin[mm],
                      gdl.histZetaMsin[mm],xi,cmd->sizeHistN);
            MULMS_ext(gdl.histZetaMsincos[mm],
                      gdl.histZetaMsincos[mm],xi,cmd->sizeHistN);
            MULMS_ext(gdl.histZetaMcossin[mm],
                      gdl.histZetaMcossin[mm],xi,cmd->sizeHistN);
        }
//B correction 2025-04-06
//    }
#ifdef NMultipoles
//    if (scanopt(cmd->options, "smooth-pivot")) {
        for (mm=1; mm<=cmd->mChebyshev+1; mm++) {
            MULMS_ext(gdlN.histZetaMcos[mm],
                      gdlN.histZetaMcos[mm],xi,cmd->sizeHistN);
            MULMS_ext(gdlN.histZetaMsin[mm],
                      gdlN.histZetaMsin[mm],xi,cmd->sizeHistN);
            MULMS_ext(gdlN.histZetaMsincos[mm],
                      gdlN.histZetaMsincos[mm],xi,cmd->sizeHistN);
            MULMS_ext(gdlN.histZetaMcossin[mm],
                      gdlN.histZetaMcossin[mm],xi,cmd->sizeHistN);
        }
//    }
#endif
    //E

//    if (scanopt(cmd->options, "smooth-pivot")) {
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
    } // ! smooth-pivot
//E correction 2025-04-06


#ifdef DEBUG
    fclose(outpivots);                              // Close file to debug pivots
#endif

//B Saving histograms section: case KKKCORRELATION:
    verb_print(cmd->verbose,
        "\n\tsearch_octree_kkk_omp: printing octree-kkk-omp method...\n\n");
    PrintHistrBins(cmd, gd);
#ifdef NMultipoles

#ifdef NONORMHIST
        if (scanopt(cmd->options, "no-normalize-HistZeta"))
            PrintHistZetaM_sincos(cmd, gd, &gdl);
        else
            PrintHistZetaM_sincos_normalized(cmd, gd, &gdl, &gdlN);
#else
        PrintHistZetaM_sincos(cmd, gd, &gdl);
#endif // ! NONORMHIST

//B correction 2025-04-06
        PrintHistZetaM_sincos_N(cmd, gd, &gdlN);
//E correction 2025-04-06

#else // ! NMultipoles
    PrintHistZetaM_sincos(cmd, gd, &gdl);
#endif // ! NMultipoles


//B correction 2025-04-06
//#ifdef NMultipoles
//        PrintHistZetaM_sincos_N(cmd, gd, &gdlN);
//#endif
//E correction 2025-04-06

        if (scanopt(cmd->options, "out-m-HistZeta")) {
#ifdef NMultipoles

#ifdef NONORMHIST
            if (scanopt(cmd->options, "no-normalize-HistZeta"))
                PrintHistZetaMm_sincos(cmd, gd, &gdl);
            else
                PrintHistZetaMm_sincos_normalized(cmd, gd, &gdl, &gdlN);
#else
            PrintHistZetaMm_sincos(cmd, gd, &gdl);
#endif // ! NONORMHIST

//B correction 2025-04-06
            PrintHistZetaMm_sincos_N(cmd, gd, &gdlN);
//E correction 2025-04-06

#else // ! NMultipoles
            PrintHistZetaMm_sincos(cmd, gd, &gdl);
#endif // ! NMultipoles

//B correction 2025-04-06
//#ifdef NMultipoles
//            PrintHistZetaMm_sincos_N(cmd, gd, &gdlN);
//#endif
//E correction 2025-04-06
        }

        if (scanopt(cmd->options, "out-HistZetaG")) {
            PrintHistZetaGm_sincos(cmd, gd, &gdl);
            PrintHistZetaG(cmd, gd, &gdl);
            PrintHistZetaMZetaGm_sincos(cmd, gd, &gdl);
        }

#ifdef NMultipoles
#ifdef NONORMHIST
    //B correction 2025-04-06
        if (scanopt(cmd->options, "no-normalize-HistZeta")) {
            if (scanopt(cmd->options, "edge-corrections")) {
                PrintHistZetaM_sincos_edge_effects(cmd, gd, &gdl, &gdlN);
            }
        } else {
            if (scanopt(cmd->options, "edge-corrections")) {
                PrintHistZetaM_sincos_edge_effects_normalized_HistZeta(cmd,
                                                                gd, &gdl, &gdlN);
            }
        }
    //E correction 2025-04-06
#endif // ! NONORMHIST
#endif // ! NMultipoles
    gd->flagPrint = FALSE;
//E Saving histograms section: case KKKCORRELATION


#ifdef NMultipoles
    search_free_gd_sincos_omp_kkk_N(cmd, gd, &gdlN);// free memory
#endif
    search_free_gd_sincos_omp_kkk(cmd, gd, &gdl); // free memory

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
                if (cmd->mChebyshev<7) {
                    CHEBYSHEVTUOMPSINCOSANY
                } else {
                    CHEBYSHEVTUOMP;
                }
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
                if (cmd->mChebyshev<7) {
                    CHEBYSHEVTUOMPSINCOSANY
                } else {
                    CHEBYSHEVTUOMP;
                }
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

//B correction 2025-04-06
#ifdef NONORMHIST
                if (scanopt(cmd->options, "no-normalize-HistZeta")) {
                    xi = Nb(q)*Kappa(q);
                } else {
                    xi = Kappa(q);
                }
#else
                xi = Kappa(q);
#endif
//E

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
                if (cmd->mChebyshev<7) {
                    CHEBYSHEVTUOMPSINCOSANY
                } else {
                    CHEBYSHEVTUOMP;
                }
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
                if (cmd->mChebyshev<7) {
                    CHEBYSHEVTUOMPNANY;
                    CHEBYSHEVTUOMPSINCOSANY
                } else {
                    CHEBYSHEVTUOMPN;
                    CHEBYSHEVTUOMP;
                }
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

//B correction 2025-04-06
//                xi = Kappa(q);
//                xiN = 1.0;
#ifdef NONORMHIST
                if (scanopt(cmd->options, "no-normalize-HistZeta")) {
                    xi = Nb(q)*Kappa(q);
                    xiN = 1.0;
                } else {
                    xi = Kappa(q);
                    xiN = 1.0;
                }
#else
                xi = Kappa(q);
                xiN = 1.0;
#endif
//E

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
                if (cmd->mChebyshev<7) {
                    CHEBYSHEVTUOMPNANY;
                    CHEBYSHEVTUOMPSINCOSANY
                } else {
                    CHEBYSHEVTUOMPN;
                    CHEBYSHEVTUOMP;
                }
                hist->nbccalcthread += 1;
            }
        }
    } // ! accept_body
}
#endif


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
    gdl->histZetaMsincos =
            dmatrix3D(1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
    gdl->histZetaMcossin =
            dmatrix3D(1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    gdl->histZetaM =
            dmatrix3D(1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    bytes_tot_local +=
            5*(cmd->mChebyshev+1)*cmd->sizeHistN*cmd->sizeHistN*sizeof(real);

    gdl->histZetaGmRe =
                dmatrix3D(1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    gdl->histZetaGmIm =
                dmatrix3D(1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    bytes_tot_local +=
            2*(cmd->mChebyshev+1)*cmd->sizeHistN*cmd->sizeHistN*sizeof(real);
    gdl->histXi3pcf = dmatrix3D(1,cmd->sizeHistPhi,1,cmd->sizeHistN,1,cmd->sizeHistN);
    bytes_tot_local +=
                (cmd->sizeHistN*cmd->sizeHistN*cmd->sizeHistPhi)*sizeof(real);

    gd->bytes_tot += bytes_tot_local;
    verb_print(cmd->verbose,
    "\nsearch_init_gd_octree_kkk: Allocated %g MByte for histograms storage.\n",
    bytes_tot_local*INMB);

    for (n = 1; n <= cmd->sizeHistN; n++)
        gdl->histNNSub[n] = 0.0;
    for (m = 1; m <= cmd->mChebyshev+1; m++) {
        CLRM_ext(gdl->histZetaMcos[m], cmd->sizeHistN);
        CLRM_ext(gdl->histZetaMsin[m], cmd->sizeHistN);
        CLRM_ext(gdl->histZetaMsincos[m], cmd->sizeHistN);
        CLRM_ext(gdl->histZetaMcossin[m], cmd->sizeHistN);
        CLRM_ext(gdl->histZetaM[m], cmd->sizeHistN);
    }
    
    gd->nbbcalc = gd->nbccalc = gd->ncccalc = 0;

    return SUCCESS;
}

local int search_free_gd_sincos_omp_kkk(struct  cmdline_data* cmd,
                                         struct  global_data* gd,
                                        gdlptr_sincos_omp_kkk gdl)
{
    free_dmatrix3D(gdl->histXi3pcf,1,cmd->sizeHistPhi,
                   1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix3D(gdl->histZetaGmIm,
                   1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix3D(gdl->histZetaGmRe,
                   1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);

    free_dmatrix3D(gdl->histZetaM,
                   1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
    free_dmatrix3D(gdl->histZetaMcossin,
                   1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix3D(gdl->histZetaMsincos,
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
#else // ! NONORMHIST
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
#endif // ! NONORMHIST

//B correction 2025-04-06
    for (m=1; m<=cmd->mChebyshev+1; m++) {
        /*
#ifdef NONORMHIST
        for (n=1; n<=cmd->sizeHistN; n++) {
            hist->histXithreadcos[m][n] /= 1.0;
            hist->histXithreadsin[m][n] /= 1.0;
        }
#else
        */
#ifndef NONORMHIST
        for (n=1; n<=cmd->sizeHistN; n++) {
            hist->histXithreadcos[m][n] /= MAX(hist->histNNSubthread[n],1.0);
            hist->histXithreadsin[m][n] /= MAX(hist->histNNSubthread[n],1.0);
        }
#endif
    }
//E correction 2025-04-06
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
                                        struct  global_data* gd,
                                          gdlptr_sincos_omp_kkk_N gdl)
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
    gdl->histZetaMsincos =
            dmatrix3D(1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
    gdl->histZetaMcossin =
            dmatrix3D(1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    bytes_tot_local +=
            4*(cmd->mChebyshev+1)*cmd->sizeHistN*cmd->sizeHistN*sizeof(real);

    gd->bytes_tot += bytes_tot_local;
    verb_print(cmd->verbose,
    "\nsearch_init_gd_octree_kkk: Allocated %g MByte for histograms storage.\n",
    bytes_tot_local*INMB);

    for (n = 1; n <= cmd->sizeHistN; n++)
        gdl->histNNSub[n] = 0.0;
    for (m = 1; m <= cmd->mChebyshev+1; m++) {
        CLRM_ext(gdl->histZetaMcos[m], cmd->sizeHistN);
        CLRM_ext(gdl->histZetaMsin[m], cmd->sizeHistN);
        CLRM_ext(gdl->histZetaMsincos[m], cmd->sizeHistN);
        CLRM_ext(gdl->histZetaMcossin[m], cmd->sizeHistN);
    }

    return SUCCESS;
}

local int search_free_gd_sincos_omp_kkk_N(struct  cmdline_data* cmd,
                                         struct  global_data* gd,
                                        gdlptr_sincos_omp_kkk_N gdl)
{
    // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
    free_dmatrix3D(gdl->histZetaMcossin,
                   1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix3D(gdl->histZetaMsincos,
                   1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix3D(gdl->histZetaMsin,
                   1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix3D(gdl->histZetaMcos,
                   1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);

    free_dvector(gdl->histNNSub,1,cmd->sizeHistN);

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

//B correction 2025-04-06
    for (m=1; m<=cmd->mChebyshev+1; m++) {
/*#ifdef NONORMHIST
        for (n=1; n<=cmd->sizeHistN; n++) {
            hist->histXithreadcos[m][n] /= 1.0;
            hist->histXithreadsin[m][n] /= 1.0;
        }
#else */
#ifndef NONORMHIST
        for (n=1; n<=cmd->sizeHistN; n++) {
            hist->histXithreadcos[m][n] /= MAX(hist->histNNSubthread[n],1.0);
            hist->histXithreadsin[m][n] /= MAX(hist->histNNSubthread[n],1.0);
        }
#endif
    }
//B correction 2025-04-06

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


#ifdef NMultipoles
#ifdef USEGSL
local int matrixClm(struct cmdline_data* cmd, struct  global_data* gd,
                    gdlptr_sincos_omp_kkk gdl, gdlptr_sincos_omp_kkk_N gdlN,
                    int n1, int n2)
{
    // 1 <= l, m <= 2*mChebyshev + 1
    // 1 <= n1, n2 <= sizeHistN
    int l, m;
    int lmx;
    int neqs=2*cmd->mChebyshev+1;
    int mx=cmd->mChebyshev;
    real C1;
    int s;

    gsl_vector * bl = gsl_vector_alloc (neqs);
    gsl_matrix * Clm = gsl_matrix_alloc (neqs, neqs);
    gsl_matrix * ClmChk = gsl_matrix_alloc (neqs, neqs);
    gsl_vector *x = gsl_vector_alloc (neqs);
    gsl_permutation * p = gsl_permutation_alloc (neqs);

    gsl_vector *t = gsl_vector_alloc (neqs);
    gsl_matrix * u = gsl_matrix_alloc (neqs, neqs);
    real v;

    for (l=0; l<neqs; l++) {
        gsl_vector_set(bl, l, 0.0);
        for (m=0; m<neqs; m++) {
            gsl_matrix_set(Clm, l, m, 0.0);
            gsl_matrix_set(ClmChk, l, m, 0.0);
        }
    }

    if (cmd->verbose_log>=3)
        verb_log_print(cmd->verbose_log, gd->outlog,"\n\nMatrix and b elements:\n\n");
    for (l=0; l<neqs; l++) {
        if (l<=mx)
            lmx = mx-(l+1)+2;
        else
            lmx = (l-1)+2-mx;

        gsl_vector_set(bl, l,
                       (gdl->histZetaMcos[lmx][n1][n2]
                        + gdl->histZetaMsin[lmx][n1][n2])
                       /(gdlN->histZetaMcos[1][n1][n2]
                         + gdlN->histZetaMsin[1][n1][n2])
                       );
        if (l<=mx) {
            if (cmd->verbose_log>=3)
                verb_log_print(cmd->verbose_log, gd->outlog,
                               "b%d : %g :: %d %d\n",
                               -(lmx-1), gsl_vector_get(bl, l), l, lmx);
        } else {
            if (cmd->verbose_log>=3)
                verb_log_print(cmd->verbose_log, gd->outlog,
                               "b%d : %g :: %d %d\n",
                               lmx-1, gsl_vector_get(bl, l), l, lmx);
        }
        for (m=0; m<neqs; m++) {
            if (l-m>=-mx && l-m<0) {
                C1 = (gdlN->histZetaMcos[m-l+1][n1][n2]
                      + gdlN->histZetaMsin[m-l+1][n1][n2]
                      )/(gdlN->histZetaMcos[1][n1][n2]
                         +gdlN->histZetaMsin[1][n1][n2]);
                gsl_matrix_set(Clm, l, m, C1);
                gsl_matrix_set( ClmChk, l, m, gsl_matrix_get(Clm, l, m) );
                if (cmd->verbose_log>=3)
                    verb_log_print(cmd->verbose_log, gd->outlog,
                                   "%g ", gsl_matrix_get(Clm, l, m));
                continue;
            }
            if (l-m>=0 && l-m<=mx) {
                C1 = (gdlN->histZetaMcos[l-m+1][n1][n2]
                      +gdlN->histZetaMsin[l-m+1][n1][n2])
                      /(gdlN->histZetaMcos[1][n1][n2]+gdlN->histZetaMsin[1][n1][n2]);
                gsl_matrix_set(Clm, l, m, C1);
                gsl_matrix_set( ClmChk, l, m, gsl_matrix_get(Clm, l, m) );
                if (cmd->verbose_log>=3)
                    verb_log_print(cmd->verbose_log, gd->outlog,
                                   "%g ", gsl_matrix_get(Clm, l, m));
                continue;
            }
            if (cmd->verbose_log>=3)
                verb_log_print(cmd->verbose_log, gd->outlog,
                               "%g ", gsl_matrix_get(Clm, l, m));
        }
        if (cmd->verbose_log>=3)
            verb_log_print(cmd->verbose_log, gd->outlog,"\n\n");
    }

    gsl_linalg_LU_decomp (Clm, p, &s);
    gsl_linalg_LU_solve (Clm, p, bl, x);
    if (cmd->verbose_log>=3) {
        verb_log_print(cmd->verbose_log, gd->outlog,"x = \n");
        gsl_vector_fprintf (gd->outlog, x, "%g");
    }

    // check A x = b
    if (cmd->verbose_log>=3) {
        verb_log_print(cmd->verbose_log, gd->outlog,"\nA x = b:\n");
        for (l=0; l<neqs; l++) {
            v = 0.0;
            for (m=0; m<neqs; m++) {
                v += ( gsl_matrix_get(ClmChk,l,m)*gsl_vector_get(x,m) );
            }
            gsl_vector_set(t, l, v);
            verb_log_print(cmd->verbose_log, gd->outlog,
                           "%8s %g %g\n"," ",
                           gsl_vector_get(bl,l),gsl_vector_get(t,l));
        }
    }

    for (l=0; l<neqs; l++) {
        if (l<=mx)
            lmx = mx-(l+1)+2;
        else
            lmx = (l-1)+2-mx;

        gdl->histZetaM[lmx][n1][n2] = gsl_vector_get(x, l);
    }
    
    gsl_matrix_free (u);
    gsl_vector_free (t);
    gsl_permutation_free (p);
    gsl_vector_free (x);
    gsl_matrix_free (ClmChk);
    gsl_matrix_free (Clm);
    gsl_vector_free (bl);

    return SUCCESS;
}

//B correction 2025-04-06
local int matrixClm_normalized_HistZeta(struct cmdline_data* cmd,
                                        struct  global_data* gd,
                    gdlptr_sincos_omp_kkk gdl, gdlptr_sincos_omp_kkk_N gdlN,
                    int n1, int n2)
{
    // 1 <= l, m <= 2*mChebyshev + 1
    // 1 <= n1, n2 <= sizeHistN
    int l, m;
    int lmx;
    int neqs=2*cmd->mChebyshev+1;
    int mx=cmd->mChebyshev;
    real C1;
    int s;

    gsl_vector * bl = gsl_vector_alloc (neqs);
    gsl_matrix * Clm = gsl_matrix_alloc (neqs, neqs);
    gsl_matrix * ClmChk = gsl_matrix_alloc (neqs, neqs);
    gsl_vector *x = gsl_vector_alloc (neqs);
    gsl_permutation * p = gsl_permutation_alloc (neqs);

    gsl_vector *t = gsl_vector_alloc (neqs);
    gsl_matrix * u = gsl_matrix_alloc (neqs, neqs);
    real v;

    for (l=0; l<neqs; l++) {
        gsl_vector_set(bl, l, 0.0);
        for (m=0; m<neqs; m++) {
            gsl_matrix_set(Clm, l, m, 0.0);
            gsl_matrix_set(ClmChk, l, m, 0.0);
        }
    }

    if (cmd->verbose_log>=3)
        verb_log_print(cmd->verbose_log, gd->outlog,"\n\nMatrix and b elements:\n\n");
    for (l=0; l<neqs; l++) {
        if (l<=mx)
            lmx = mx-(l+1)+2;
        else
            lmx = (l-1)+2-mx;

        gsl_vector_set(bl, l,
                       (gdl->histZetaMcos[lmx][n1][n2]
                        + gdl->histZetaMsin[lmx][n1][n2])
                       /(gdlN->histZetaMcos[1][n1][n2]
                         + gdlN->histZetaMsin[1][n1][n2])
                       );
        if (l<=mx) {
            if (cmd->verbose_log>=3)
                verb_log_print(cmd->verbose_log, gd->outlog,
                               "b%d : %g :: %d %d\n",
                               -(lmx-1), gsl_vector_get(bl, l), l, lmx);
        } else {
            if (cmd->verbose_log>=3)
                verb_log_print(cmd->verbose_log, gd->outlog,
                               "b%d : %g :: %d %d\n",
                               lmx-1, gsl_vector_get(bl, l), l, lmx);
        }
        for (m=0; m<neqs; m++) {
            if (l-m>=-mx && l-m<0) {
                C1 = (gdlN->histZetaMcos[m-l+1][n1][n2]
                      + gdlN->histZetaMsin[m-l+1][n1][n2]
                      )/(gdlN->histZetaMcos[1][n1][n2]
                         +gdlN->histZetaMsin[1][n1][n2]);
                gsl_matrix_set(Clm, l, m, C1);
                gsl_matrix_set( ClmChk, l, m, gsl_matrix_get(Clm, l, m) );
                if (cmd->verbose_log>=3)
                    verb_log_print(cmd->verbose_log, gd->outlog,
                                   "%g ", gsl_matrix_get(Clm, l, m));
                continue;
            }
            if (l-m>=0 && l-m<=mx) {
                C1 = (gdlN->histZetaMcos[l-m+1][n1][n2]
                      +gdlN->histZetaMsin[l-m+1][n1][n2])
                      /(gdlN->histZetaMcos[1][n1][n2]+gdlN->histZetaMsin[1][n1][n2]);
                gsl_matrix_set(Clm, l, m, C1);
                gsl_matrix_set( ClmChk, l, m, gsl_matrix_get(Clm, l, m) );
                if (cmd->verbose_log>=3)
                    verb_log_print(cmd->verbose_log, gd->outlog,
                                   "%g ", gsl_matrix_get(Clm, l, m));
                continue;
            }
            if (cmd->verbose_log>=3)
                verb_log_print(cmd->verbose_log, gd->outlog,
                               "%g ", gsl_matrix_get(Clm, l, m));
        }
        if (cmd->verbose_log>=3)
            verb_log_print(cmd->verbose_log, gd->outlog,"\n\n");
    }

    gsl_linalg_LU_decomp (Clm, p, &s);
    gsl_linalg_LU_solve (Clm, p, bl, x);
    if (cmd->verbose_log>=3) {
        verb_log_print(cmd->verbose_log, gd->outlog,"x = \n");
        gsl_vector_fprintf (gd->outlog, x, "%g");
    }

    // check A x = b
    if (cmd->verbose_log>=3) {
        verb_log_print(cmd->verbose_log, gd->outlog,"\nA x = b:\n");
        for (l=0; l<neqs; l++) {
            v = 0.0;
            for (m=0; m<neqs; m++) {
                v += ( gsl_matrix_get(ClmChk,l,m)*gsl_vector_get(x,m) );
            }
            gsl_vector_set(t, l, v);
            verb_log_print(cmd->verbose_log, gd->outlog,
                           "%8s %g %g\n"," ",
                           gsl_vector_get(bl,l),gsl_vector_get(t,l));
        }
    }

    for (l=0; l<neqs; l++) {
        if (l<=mx)
            lmx = mx-(l+1)+2;
        else
            lmx = (l-1)+2-mx;

        gdl->histZetaM[lmx][n1][n2] = gsl_vector_get(x, l);
    }
    
    gsl_matrix_free (u);
    gsl_vector_free (t);
    gsl_permutation_free (p);
    gsl_vector_free (x);
    gsl_matrix_free (ClmChk);
    gsl_matrix_free (Clm);
    gsl_vector_free (bl);

    return SUCCESS;
}
//E correction 2025-04-06

#else // ! USEGSL
// routine to compute matriz elements Clm and its inverse
//  Be careful with the use of NONORMHISTON and NMultipoles switches:
//  NMultipolesON = 1 and NONORMHISTON = 1
local int matrixClm(struct cmdline_data* cmd, struct  global_data* gd,
                    gdlptr_sincos_omp_kkk gdl, gdlptr_sincos_omp_kkk_N gdlN,
                                         int n1, int n2)
{
// 1 <= l, m <= 2*mChebyshev + 1
// 1 <= n1, n2 <= sizeHistN
    int l, m;
    int j;
    int *indx;
    real p;
    real **Clm;
    real **ClmChk;
    real **u;
    real *bl;
    real *blChk;
    real *t;
    real C1;

    int lm, lp;
    int neqs=2*cmd->mChebyshev+1;
    int mx=cmd->mChebyshev;
    int lmx;

    Clm = dmatrix(1,neqs,1,neqs);
    ClmChk = dmatrix(1,neqs,1,neqs);
    u = dmatrix(1,neqs,1,neqs);
    bl = dvector(1,neqs);
    blChk = dvector(1,neqs);
    t = dvector(1,neqs);
    indx = ivector(1,neqs);

    CLRM_ext(Clm,neqs);
    CLRM_ext(ClmChk,neqs);
    CLRV_ext(bl,neqs);

    if (cmd->verbose_log>=3)
        verb_log_print(cmd->verbose_log, gd->outlog,"\n\nMatrix and b elements:\n\n");
    for (l=1; l<=neqs; l++) {
        if (l<=mx+1)
            lmx = (mx-(l+1)+2) +1;
        else
            lmx = (l-2)+2-mx;

        bl[l] = (gdl->histZetaMcos[lmx][n1][n2] + gdl->histZetaMsin[lmx][n1][n2])
                       /gdlN->histZetaMcos[1][n1][n2];
        blChk[l] = bl[l];
        if (l<=mx+1) {
            if (cmd->verbose_log>=3)
                verb_log_print(cmd->verbose_log, gd->outlog,
                               "b%d : %g :: %d %d\n",
                               -(lmx-1), bl[l], l, lmx);
        } else {
            if (cmd->verbose_log>=3)
                verb_log_print(cmd->verbose_log, gd->outlog,
                               "b%d : %g :: %d %d\n",
                               lmx-1, bl[l], l, lmx);
        }
        for (m=1; m<=neqs; m++) {
            if (l-m>=-mx && l-m<0) {
                C1 = (gdlN->histZetaMcos[m-l+1][n1][n2]
                      + gdlN->histZetaMsin[m-l+1][n1][n2]
                      )/gdlN->histZetaMcos[1][n1][n2];
                Clm[l][m] = C1;
                ClmChk[l][m] = Clm[l][m];
                if (cmd->verbose_log>=3)
                    verb_log_print(cmd->verbose_log, gd->outlog,
                                   "%g ", Clm[l][m]);
                continue;
            }
            if (l-m>=0 && l-m<=mx) {
                C1 = (gdlN->histZetaMcos[l-m+1][n1][n2]
                      +gdlN->histZetaMsin[l-m+1][n1][n2])
                      /gdlN->histZetaMcos[1][n1][n2];
                Clm[l][m] = C1;
                ClmChk[l][m] = Clm[l][m];
                if (cmd->verbose_log>=3)
                    verb_log_print(cmd->verbose_log, gd->outlog,
                                   "%g ", Clm[l][m]);
                continue;
            }
            if (cmd->verbose_log>=3)
                verb_log_print(cmd->verbose_log, gd->outlog,
                               "%g ", Clm[l][m]);
        }
        if (cmd->verbose_log>=3)
            verb_log_print(cmd->verbose_log, gd->outlog,"\n\n");
    }

    ludcmp(Clm,neqs,indx,&p);
    lubksb(Clm,neqs,indx,bl);

    // vector solutions
    if (cmd->verbose_log>=3) {
        verb_log_print(cmd->verbose_log, gd->outlog,
                       "\nVector solution:\n");
        for (l=1;l<=neqs;l++) {
            verb_log_print(cmd->verbose_log, gd->outlog,"%8s %g\n"," ", bl[l]);
        }
    }

    // check A x = b
    real v;
    if (cmd->verbose_log>=3) {
        verb_log_print(cmd->verbose_log, gd->outlog,"\nA x = b:\n");
        for (l=1; l<=neqs; l++) {
            v = 0.0;
            for (m=0; m<neqs; m++) {
                v += ClmChk[l][m]*bl[m];
            }
            t[l] = v;
            verb_log_print(cmd->verbose_log, gd->outlog,
                           "%8s %g %g\n"," ",blChk[l],t[l]);
        }
    }

//B correction 2025-04-06
    for (l=1; l<=neqs; l++) {
        if (l<=mx+1)
            lmx = (mx-(l+1)+2) +1;
        else
            lmx = (l-2)+2-mx;

        gdl->histZetaM[lmx][n1][n2] = bl[l];
    }
//E correction 2025-04-06

    free_ivector(indx,1,cmd->mChebyshev+1);
    free_dvector(t,1,cmd->mChebyshev+1);
    free_dvector(blChk,1,cmd->mChebyshev+1);
    free_dvector(bl,1,cmd->mChebyshev+1);
    free_dmatrix(u,1,cmd->mChebyshev+1,1,cmd->mChebyshev+1);
    free_dmatrix(ClmChk,1,cmd->mChebyshev+1,1,cmd->mChebyshev+1);
    free_dmatrix(Clm,1,cmd->mChebyshev+1,1,cmd->mChebyshev+1);

    return SUCCESS;
}

//B correction 2025-04-06
local int matrixClm_normalized_HistZeta(struct cmdline_data* cmd,
                                        struct  global_data* gd,
                    gdlptr_sincos_omp_kkk gdl, gdlptr_sincos_omp_kkk_N gdlN,
                                         int n1, int n2)
{
// 1 <= l, m <= 2*mChebyshev + 1
// 1 <= n1, n2 <= sizeHistN
    int l, m;
    int j;
    int *indx;
    real p;
    real **Clm;
    real **ClmChk;
    real **u;
    real *bl;
    real *blChk;
    real *t;
    real C1;

    int lm, lp;
    int neqs=2*cmd->mChebyshev+1;
    int mx=cmd->mChebyshev;
    int lmx;

    Clm = dmatrix(1,neqs,1,neqs);
    ClmChk = dmatrix(1,neqs,1,neqs);
    u = dmatrix(1,neqs,1,neqs);
    bl = dvector(1,neqs);
    blChk = dvector(1,neqs);
    t = dvector(1,neqs);
    indx = ivector(1,neqs);

    CLRM_ext(Clm,neqs);
    CLRM_ext(ClmChk,neqs);
    CLRV_ext(bl,neqs);

    if (cmd->verbose_log>=3)
        verb_log_print(cmd->verbose_log, gd->outlog,"\n\nMatrix and b elements:\n\n");
    for (l=1; l<=neqs; l++) {
        if (l<=mx+1)
            lmx = (mx-(l+1)+2) +1;
        else
            lmx = (l-2)+2-mx;

        bl[l] = (gdl->histZetaMcos[lmx][n1][n2] + gdl->histZetaMsin[lmx][n1][n2])
                       /gdlN->histZetaMcos[1][n1][n2];
        blChk[l] = bl[l];
        if (l<=mx+1) {
            if (cmd->verbose_log>=3)
                verb_log_print(cmd->verbose_log, gd->outlog,
                               "b%d : %g :: %d %d\n",
                               -(lmx-1), bl[l], l, lmx);
        } else {
            if (cmd->verbose_log>=3)
                verb_log_print(cmd->verbose_log, gd->outlog,
                               "b%d : %g :: %d %d\n",
                               lmx-1, bl[l], l, lmx);
        }
        for (m=1; m<=neqs; m++) {
            if (l-m>=-mx && l-m<0) {
                C1 = (gdlN->histZetaMcos[m-l+1][n1][n2]
                      + gdlN->histZetaMsin[m-l+1][n1][n2]
                      )/gdlN->histZetaMcos[1][n1][n2];
                Clm[l][m] = C1;
                ClmChk[l][m] = Clm[l][m];
                if (cmd->verbose_log>=3)
                    verb_log_print(cmd->verbose_log, gd->outlog,
                                   "%g ", Clm[l][m]);
                continue;
            }
            if (l-m>=0 && l-m<=mx) {
                C1 = (gdlN->histZetaMcos[l-m+1][n1][n2]
                      +gdlN->histZetaMsin[l-m+1][n1][n2])
                      /gdlN->histZetaMcos[1][n1][n2];
                Clm[l][m] = C1;
                ClmChk[l][m] = Clm[l][m];
                if (cmd->verbose_log>=3)
                    verb_log_print(cmd->verbose_log, gd->outlog,
                                   "%g ", Clm[l][m]);
                continue;
            }
            if (cmd->verbose_log>=3)
                verb_log_print(cmd->verbose_log, gd->outlog,
                               "%g ", Clm[l][m]);
        }
        if (cmd->verbose_log>=3)
            verb_log_print(cmd->verbose_log, gd->outlog,"\n\n");
    }

    ludcmp(Clm,neqs,indx,&p);
    lubksb(Clm,neqs,indx,bl);

    // vector solutions
    if (cmd->verbose_log>=3) {
        verb_log_print(cmd->verbose_log, gd->outlog,
                       "\nVector solution:\n");
        for (l=1;l<=neqs;l++) {
            verb_log_print(cmd->verbose_log, gd->outlog,"%8s %g\n"," ", bl[l]);
        }
    }

    // check A x = b
    real v;
    if (cmd->verbose_log>=3) {
        verb_log_print(cmd->verbose_log, gd->outlog,"\nA x = b:\n");
        for (l=1; l<=neqs; l++) {
            v = 0.0;
            for (m=0; m<neqs; m++) {
                v += ClmChk[l][m]*bl[m];
            }
            t[l] = v;
            verb_log_print(cmd->verbose_log, gd->outlog,
                           "%8s %g %g\n"," ",blChk[l],t[l]);
        }
    }

//B correction 2025-04-06
    for (l=1; l<=neqs; l++) {
        if (l<=mx+1)
            lmx = (mx-(l+1)+2) +1;
        else
            lmx = (l-2)+2-mx;

        gdl->histZetaM[lmx][n1][n2] = bl[l];
    }
//E correction 2025-04-06

    free_ivector(indx,1,cmd->mChebyshev+1);
    free_dvector(t,1,cmd->mChebyshev+1);
    free_dvector(blChk,1,cmd->mChebyshev+1);
    free_dvector(bl,1,cmd->mChebyshev+1);
    free_dmatrix(u,1,cmd->mChebyshev+1,1,cmd->mChebyshev+1);
    free_dmatrix(ClmChk,1,cmd->mChebyshev+1,1,cmd->mChebyshev+1);
    free_dmatrix(Clm,1,cmd->mChebyshev+1,1,cmd->mChebyshev+1);

    return SUCCESS;
}
//E correction 2025-04-06

#endif // ! USEGSL
#endif

//E Routines as in cballsutils

local int print_info(struct cmdline_data* cmd,
                                  struct  global_data* gd)
{
    verb_print(cmd->verbose,
               "searchcalc: Using octree-kkk-omp... \n");

    if (scanopt(cmd->options, "KKKCorrelation"))
        verb_print(cmd->verbose, "computing KKK correlation... \n");


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

#if defined(NMultipoles) && defined(NONORMHIST)
    if (scanopt(cmd->options, "no-normalize-HistZeta")) {
        if (scanopt(cmd->options, "edge-corrections"))
            verb_print(cmd->verbose, "with option edge-corrections... \n");
    } else {
        if (scanopt(cmd->options, "edge-corrections")) {
            verb_print(cmd->verbose,
                       "option edge-corrections only works with %s... \n",
                       "no-normalize-HistZeta option added");
            // Check freeing allocated memory...
            error("going out...\n");
        }
    }
#else
    if (scanopt(cmd->options, "edge-corrections")) {
        verb_print(cmd->verbose,
                   "option edge-corrections only works with %s activated... \n",
                   "NMultipoles && NONORMHIST");
        // Check freeing allocated memory...
        error("going out...");
    }
#endif

#ifndef USEGSL
    if (scanopt(cmd->options, "edge-corrections"))
        verb_print(cmd->verbose,
                   "option edge-corrections is better computed with %s activated... \n",
                   "USEGSL");
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
    for (m=1; m<=cmd->mChebyshev+1; m++) {
        sprintf(namebuf, "%s%s_%d%s", gd->fpfnamehistZetaMFileName,
                "_sincos", m, EXTFILES);
        verb_print_q(2, cmd->verbose, "Printing : to a file %s ...\n",namebuf);
        outstr = stropen(namebuf, "w!");
        for (n1=1; n1<=cmd->sizeHistN; n1++) {
            for (n2=1; n2<=cmd->sizeHistN; n2++) {
                fprintf(outstr,"%16.8e ",gdl->histZetaMsincos[m][n1][n2]);
            }
            fprintf(outstr,"\n");
        }
        fclose(outstr);
    }
    // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
    for (m=1; m<=cmd->mChebyshev+1; m++) {
        sprintf(namebuf, "%s%s_%d%s", gd->fpfnamehistZetaMFileName,
                "_cossin", m, EXTFILES);
        verb_print_q(2, cmd->verbose, "Printing : to a file %s ...\n",namebuf);
        outstr = stropen(namebuf, "w!");
        for (n1=1; n1<=cmd->sizeHistN; n1++) {
            for (n2=1; n2<=cmd->sizeHistN; n2++) {
                fprintf(outstr,"%16.8e ",gdl->histZetaMcossin[m][n1][n2]);
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

    for (m = 1; m <= cmd->mChebyshev+1; m++) {
        sprintf(namebuf, "%s%s_%d%s", gd->fpfnamemhistZetaMFileName,
                "_sincos", m, EXTFILES);
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
            Zeta = gdl->histZetaMsincos[m][n1][n1];
            Zeta2 = gdl->histZetaMsincos[m][n1][(int)(Nbins/4.0)];
            Zeta3 = gdl->histZetaMsincos[m][n1][(int)(2.0*Nbins/4.0)];
            Zeta4 = gdl->histZetaMsincos[m][n1][(int)(3.0*Nbins/4.0)];
            Zeta5 = gdl->histZetaMsincos[m][n1][(int)(4.0*Nbins/4.0 - 1.0)];
            fprintf(outstr,MHISTZETA,rBin,Zeta,Zeta2,Zeta3,Zeta4,Zeta5);
        }
        fclose(outstr);
    }

    // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
    for (m = 1; m <= cmd->mChebyshev+1; m++) {
        sprintf(namebuf, "%s%s_%d%s", gd->fpfnamemhistZetaMFileName,
                "_cossin", m, EXTFILES);
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
            Zeta = gdl->histZetaMcossin[m][n1][n1];
            Zeta2 = gdl->histZetaMcossin[m][n1][(int)(Nbins/4.0)];
            Zeta3 = gdl->histZetaMcossin[m][n1][(int)(2.0*Nbins/4.0)];
            Zeta4 = gdl->histZetaMcossin[m][n1][(int)(3.0*Nbins/4.0)];
            Zeta5 = gdl->histZetaMcossin[m][n1][(int)(4.0*Nbins/4.0 - 1.0)];
            fprintf(outstr,MHISTZETA,rBin,Zeta,Zeta2,Zeta3,Zeta4,Zeta5);
        }
        fclose(outstr);
    }

    return SUCCESS;
}


// Saves matrix ZetaG, full correlation function at each phi bins
local int PrintHistZetaG(struct  cmdline_data* cmd,
                        struct  global_data* gd, gdlptr_sincos_omp_kkk gdl)
{
    int n1, n2, l, m;
    stream outstr;
    char namebuf[256];
    real theta;
    real ***histXi3pcfIm;

    histXi3pcfIm =
        dmatrix3D(1,cmd->sizeHistPhi,1,cmd->sizeHistN,1,cmd->sizeHistN);

    for (l=1; l<=cmd->sizeHistPhi; l++) {
        sprintf(namebuf, "%s%s_%d%s", gd->fpfnamehistZetaGFileName,
                "_Xi3pcf_",l, EXTFILES);
        theta = (real)l * gd->deltaPhi;
        verb_print_q(2, cmd->verbose,
                    "Printing : to a file %s with theta %g...\n",namebuf, theta);
        outstr = stropen(namebuf, "w!");
        for (n1=1; n1<=cmd->sizeHistN; n1++) {
            for (n2=1; n2<=cmd->sizeHistN; n2++) {
                gdl->histXi3pcf[l][n1][n2] = gdl->histZetaMcos[1][n1][n2]
                                            + gdl->histZetaMsin[1][n1][n2];
                histXi3pcfIm[l][n1][n2] = 0.0;
                for (m=2; m<=cmd->mChebyshev+1; m++) {
                    gdl->histXi3pcf[l][n1][n2] += 2.0*(
                            gdl->histZetaMcos[m][n1][n2]
                            + gdl->histZetaMsin[m][n1][n2]
                                                       )*rcos((real)m*theta);
                    //B This will be a measurement error
                    histXi3pcfIm[l][n1][n2] +=
                                            2.0*(
                            gdl->histZetaMsincos[m][n1][n2]
                    // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
                            - gdl->histZetaMcossin[m][n1][n2]
                                                 )*rcos((real)m*theta);
                    //E
                }
                fprintf(outstr,"%16.8e ",gdl->histXi3pcf[l][n1][n2]);
            }
            fprintf(outstr,"\n");
        }
        fclose(outstr);
    }

    free_dmatrix3D(histXi3pcfIm,1,cmd->sizeHistPhi,
                   1,cmd->sizeHistN,1,cmd->sizeHistN);

    return SUCCESS;
}



// Saves matrix ZetaM as obtained from ZetaG, for each m multipole
//  It seems this is doing in routine below... check and delete this if the case
local int PrintHistZetaGm_sincos(struct  cmdline_data* cmd,
                                 struct  global_data* gd,
                                 gdlptr_sincos_omp_kkk gdl)
{
    int n1, n2, m;
    stream outstr;
    char namebuf[256];

    for (m=1; m<=cmd->mChebyshev+1; m++) {
        sprintf(namebuf, "%s%s_%d%s", gd->fpfnamehistZetaGmFileName,
                "_Re", m, EXTFILES);
        verb_print_q(2, cmd->verbose, "Printing : to a file %s ...\n",namebuf);
        outstr = stropen(namebuf, "w!");
        for (n1=1; n1<=cmd->sizeHistN; n1++) {
            for (n2=1; n2<=cmd->sizeHistN; n2++) {
                fprintf(outstr,"%16.8e ",gdl->histZetaGmRe[m][n1][n2]);
            }
            fprintf(outstr,"\n");
        }
        fclose(outstr);
    }
    for (m=1; m<=cmd->mChebyshev+1; m++) {
        sprintf(namebuf, "%s%s_%d%s", gd->fpfnamehistZetaGmFileName,
                "_Im", m, EXTFILES);
        verb_print_q(2, cmd->verbose, "Printing : to a file %s ...\n",namebuf);
        outstr = stropen(namebuf, "w!");
        for (n1=1; n1<=cmd->sizeHistN; n1++) {
            for (n2=1; n2<=cmd->sizeHistN; n2++) {
                fprintf(outstr,"%16.8e ",gdl->histZetaGmIm[m][n1][n2]);
            }
            fprintf(outstr,"\n");
        }
        fclose(outstr);
    }

    return SUCCESS;
}


// Saves matrix ZetaG, real and imaginary parts, obtained from ZetaM multipoles
//  also saves full 3pcf ZetaG matrix for each phi bins obtained from inverse FFT
local int PrintHistZetaMZetaGm_sincos(struct  cmdline_data* cmd,
                                     struct  global_data* gd,
                                      gdlptr_sincos_omp_kkk gdl)
{
    int n1, n2, m, l;
    stream outstr;
    char namebuf[256];

    int NP = 2*(cmd->mChebyshev+1);
    double ***histZetaG;
    double ***histZetaG_Im;
    histZetaG = dmatrix3D(1,NP,1,cmd->sizeHistN,1,cmd->sizeHistN);
    histZetaG_Im = dmatrix3D(1,NP,1,cmd->sizeHistN,1,cmd->sizeHistN);

#ifdef USEGSL
    double *data;
    gsl_fft_real_wavetable * real;
    gsl_fft_real_workspace * work;
    gsl_fft_halfcomplex_wavetable * hc;

    //B Test and check this allocation of memory...
    //data=dvector(0,NP-1);
    data=(double *)allocate(NP*sizeof(double));
    //E

    work = gsl_fft_real_workspace_alloc (NP);
    real = gsl_fft_real_wavetable_alloc (NP);
    hc = gsl_fft_halfcomplex_wavetable_alloc (NP);
#else
    double *data;
    data=dvector(1,NP);
#endif

    //B Sum cos^2 + sin^2 and sincos - sincos
    // mchebyshev + 1 < sizeHistPhi/2
    // and mchebyshev + 1 must be a power of 2 also
    for (n1=1; n1<=cmd->sizeHistN; n1++) {
        for (n2=1; n2<=cmd->sizeHistN; n2++) {
            for (m=1; m<=cmd->mChebyshev+1; m++) {
                gdl->histZetaGmRe[m][n1][n2] =
                        gdl->histZetaMcos[m][n1][n2]
                        + gdl->histZetaMsin[m][n1][n2];
                gdl->histZetaGmIm[m][n1][n2] =
                        gdl->histZetaMsincos[m][n1][n2]
                // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
                        - gdl->histZetaMcossin[m][n1][n2];
            }
#ifdef USEGSL
            for (m=0; m<cmd->mChebyshev+1; m++) {
                data[2*m] = gdl->histZetaGmRe[m+1][n1][n2];
                data[2*m+1] = gdl->histZetaGmIm[m+1][n1][n2];
            }

            gsl_fft_complex_radix2_inverse (data, 1, NP/2);
            for (l=0; l<NP; l++) {                 // l denote angular separation
                histZetaG[l+1][n1][n2] = data[l];
            }
#else
            int isign = -1;                         // sign in imaginary unit
            for (m=1; m<=cmd->mChebyshev+1; m++) {
                data[2*m-1] = gdl->histZetaGmRe[m][n1][n2];
                data[2*m] = gdl->histZetaGmIm[m][n1][n2];
            }
            dfour1(data,NP/2,isign);                // Inverse Fourier transform
                                                    // data has Re and Im parts
            for (l=1; l<=NP; l++) {                 // l denote angular
                                                    //  separation
                histZetaG[l][n1][n2] = (2.0/(double)NP)*data[l];
            }
#endif
        }
    }
    //E

    for (m=1; m<=cmd->mChebyshev+1; m++) {
        sprintf(namebuf, "%s%s_%d%s", gd->fpfnamehistZetaGmFileName,
                "_Re", m, EXTFILES);
        verb_print_q(2, cmd->verbose,
                    "Printing : to a file %s ...\n",namebuf);
        outstr = stropen(namebuf, "w!");
        for (n1=1; n1<=cmd->sizeHistN; n1++) {
            for (n2=1; n2<=cmd->sizeHistN; n2++) {
                fprintf(outstr,"%16.8e ",gdl->histZetaGmRe[m][n1][n2]);
            }
            fprintf(outstr,"\n");
        }
        fclose(outstr);
    }

    for (m=1; m<=cmd->mChebyshev+1; m++) {
        sprintf(namebuf, "%s%s_%d%s", gd->fpfnamehistZetaGmFileName,
                "_Im", m, EXTFILES);
        verb_print_q(2, cmd->verbose,
                    "Printing : to a file %s ...\n",namebuf);
        outstr = stropen(namebuf, "w!");
        for (n1=1; n1<=cmd->sizeHistN; n1++) {
            for (n2=1; n2<=cmd->sizeHistN; n2++) {
                fprintf(outstr,"%16.8e ",gdl->histZetaGmIm[m][n1][n2]);
            }
            fprintf(outstr,"\n");
        }
        fclose(outstr);
    }

    for (l=1; l<=cmd->mChebyshev+1; l++) {
        sprintf(namebuf, "%s%s_%d%s", gd->fpfnamehistZetaGFileName,
                "_fftinv_Re",l, EXTFILES);
        verb_print_q(2, cmd->verbose,
                    "Printing : to a file %s ...\n",namebuf);
        outstr = stropen(namebuf, "w!");
        for (n1=1; n1<=cmd->sizeHistN; n1++) {
            for (n2=1; n2<=cmd->sizeHistN; n2++) {
                fprintf(outstr,"%16.8e ",histZetaG[2*l-1][n1][n2]);
            }
            fprintf(outstr,"\n");
        }
        fclose(outstr);
    }
    for (l=1; l<=cmd->mChebyshev+1; l++) {
        sprintf(namebuf, "%s%s_%d%s", gd->fpfnamehistZetaGFileName,
                "_fftinv_Im",l, EXTFILES);
        verb_print_q(2, cmd->verbose,
                    "Printing : to a file %s ...\n",namebuf);
        outstr = stropen(namebuf, "w!");
        for (n1=1; n1<=cmd->sizeHistN; n1++) {
            for (n2=1; n2<=cmd->sizeHistN; n2++) {
                fprintf(outstr,"%16.8e ",histZetaG[2*l][n1][n2]);
            }
            fprintf(outstr,"\n");
        }
        fclose(outstr);
    }

    //B Sum cos^2 + sin^2 and sincos - sincos
    // mchebyshev + 1 < sizeHistPhi/2
    // and mchebyshev + 1 must be a power of 2 also
    double deltaPhi = TWOPI/((double)NP);
    for (n1=1; n1<=cmd->sizeHistN; n1++) {
        for (n2=1; n2<=cmd->sizeHistN; n2++) {
            for (l=1; l<=NP; l++) {              // l denote angular sep.
                histZetaG[l][n1][n2] = gdl->histZetaMcos[1][n1][n2]
                                        + gdl->histZetaMsin[1][n1][n2];
                histZetaG_Im[l][n1][n2] = gdl->histZetaMsincos[1][n1][n2]
                // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
                                            - gdl->histZetaMcossin[1][n1][n2];
                for (m=2; m<=cmd->mChebyshev+1; m++) {
                    histZetaG[l][n1][n2] += 2.0*(gdl->histZetaMcos[m][n1][n2]
                                            + gdl->histZetaMsin[m][n1][n2])
                                            *rcos(((double)(m*l))*deltaPhi);
                    histZetaG_Im[l][n1][n2] +=
                                    2.0*(gdl->histZetaMsincos[m][n1][n2]
                    // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
                                    - gdl->histZetaMcossin[m][n1][n2])
                                    *rcos(((double)(m*l))*deltaPhi);
                }
            }
        }
    }
    //E
    for (l=1; l<=NP; l++) {
        sprintf(namebuf, "%s_%s_%d%s",
                gd->fpfnamehistZetaGFileName, "Re", l, EXTFILES);
        verb_print_q(2, cmd->verbose,
                    "Printing : to a file %s ...\n",namebuf);
        outstr = stropen(namebuf, "w!");
        for (n1=1; n1<=cmd->sizeHistN; n1++) {
            for (n2=1; n2<=cmd->sizeHistN; n2++) {
                fprintf(outstr,"%16.8e ",histZetaG[l][n1][n2]);
            }
            fprintf(outstr,"\n");
        }
        fclose(outstr);
    }
    for (l=1; l<=NP; l++) {
        sprintf(namebuf, "%s_%s_%d%s",
                gd->fpfnamehistZetaGFileName, "Im", l, EXTFILES);
        verb_print_q(2, cmd->verbose,
                    "Printing : to a file %s ...\n",namebuf);
        outstr = stropen(namebuf, "w!");
        for (n1=1; n1<=cmd->sizeHistN; n1++) {
            for (n2=1; n2<=cmd->sizeHistN; n2++) {
                fprintf(outstr,"%16.8e ",histZetaG_Im[l][n1][n2]);
            }
            fprintf(outstr,"\n");
        }
        fclose(outstr);
    }

#ifdef USEGSL
    gsl_fft_halfcomplex_wavetable_free (hc);
    gsl_fft_real_wavetable_free (real);
    gsl_fft_real_workspace_free (work);
//    free_dvector(data,0,NP-1);
    free(data);
#else
    free_dvector(data,1,NP);
#endif

    free_dmatrix3D(histZetaG_Im,1,NP,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix3D(histZetaG,1,NP,1,cmd->sizeHistN,1,cmd->sizeHistN);

    return SUCCESS;
}


//#ifdef ADDONS
//#include "cballs_include_05.h"
//#endif


//#undef MHISTZETAHEADER
//#undef MHISTZETA



#ifdef NMultipoles
// Saves matrix ZetaM for each m multipole
local int PrintHistZetaM_sincos_N(struct  cmdline_data* cmd,
                                struct  global_data* gd,
                                  gdlptr_sincos_omp_kkk_N gdlN)
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
    for (m=1; m<=cmd->mChebyshev+1; m++) {
        sprintf(namebuf, "%s%s_%d%s", gd->fpfnamehistZetaMFileName,
                "_sincos_N", m, EXTFILES);
        verb_print_q(2, cmd->verbose, "Printing : to a file %s ...\n",namebuf);
        outstr = stropen(namebuf, "w!");
        for (n1=1; n1<=cmd->sizeHistN; n1++) {
            for (n2=1; n2<=cmd->sizeHistN; n2++) {
                fprintf(outstr,"%16.8e ",gdlN->histZetaMsincos[m][n1][n2]);
            }
            fprintf(outstr,"\n");
        }
        fclose(outstr);
    }
    // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
    for (m=1; m<=cmd->mChebyshev+1; m++) {
        sprintf(namebuf, "%s%s_%d%s", gd->fpfnamehistZetaMFileName,
                "_cossin_N", m, EXTFILES);
        verb_print_q(2, cmd->verbose, "Printing : to a file %s ...\n",namebuf);
        outstr = stropen(namebuf, "w!");
        for (n1=1; n1<=cmd->sizeHistN; n1++) {
            for (n2=1; n2<=cmd->sizeHistN; n2++) {
                fprintf(outstr,"%16.8e ",gdlN->histZetaMcossin[m][n1][n2]);
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
                                   gdlptr_sincos_omp_kkk_N gdlN)
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

    for (m = 1; m <= cmd->mChebyshev+1; m++) {
        sprintf(namebuf, "%s%s_%d%s", gd->fpfnamemhistZetaMFileName,
                "_sincos_N", m, EXTFILES);
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
            Zeta = gdlN->histZetaMsincos[m][n1][n1];
            Zeta2 = gdlN->histZetaMsincos[m][n1][(int)(Nbins/4.0)];
            Zeta3 = gdlN->histZetaMsincos[m][n1][(int)(2.0*Nbins/4.0)];
            Zeta4 = gdlN->histZetaMsincos[m][n1][(int)(3.0*Nbins/4.0)];
            Zeta5 = gdlN->histZetaMsincos[m][n1][(int)(4.0*Nbins/4.0 - 1.0)];
            fprintf(outstr,MHISTZETA,rBin,Zeta,Zeta2,Zeta3,Zeta4,Zeta5);
        }
        fclose(outstr);
    }

    // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
    for (m = 1; m <= cmd->mChebyshev+1; m++) {
        sprintf(namebuf, "%s%s_%d%s", gd->fpfnamemhistZetaMFileName,
                "_cossin_", m, EXTFILES);
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
            Zeta = gdlN->histZetaMcossin[m][n1][n1];
            Zeta2 = gdlN->histZetaMcossin[m][n1][(int)(Nbins/4.0)];
            Zeta3 = gdlN->histZetaMcossin[m][n1][(int)(2.0*Nbins/4.0)];
            Zeta4 = gdlN->histZetaMcossin[m][n1][(int)(3.0*Nbins/4.0)];
            Zeta5 = gdlN->histZetaMcossin[m][n1][(int)(4.0*Nbins/4.0 - 1.0)];
            fprintf(outstr,MHISTZETA,rBin,Zeta,Zeta2,Zeta3,Zeta4,Zeta5);
        }
        fclose(outstr);
    }

    return SUCCESS;
}
//#endif


//B Additional definitions... delete ASAP

/*
//#ifdef NMultipoles
local int PrintHistXi2pcf_N(struct  cmdline_data* cmd, struct  global_data* gd)
{
    real rBin, rbinlog;
    int n;
    stream outstr;
    char namebuf[256];

    sprintf(namebuf, "%s%s%s%s", gd->fpfnamehistXi2pcfFileName,
            "_N", cmd->suffixOutFiles, EXTFILES);
    verb_print_q(2, cmd->verbose, "Printing : to a file %s ...\n",namebuf);
    outstr = stropen(namebuf, "w!");

    for (n=1; n<=cmd->sizeHistN; n++) {
        if (cmd->useLogHist) {
            if (cmd->rminHist==0) {
                rbinlog = ((real)(n-cmd->sizeHistN))/cmd->logHistBinsPD
                            + rlog10(cmd->rangeN);
            } else {
                rbinlog = rlog10(cmd->rminHist) + ((real)(n)-0.5)*gd->deltaR;
            }
            rBin=rpow(10.0,rbinlog);
        } else {
            rBin = cmd->rminHist + ((real)n-0.5)*gd->deltaR;
        }
        fprintf(outstr,"%16.8e %16.8e\n",rBin,gd->NhistXi2pcf[n]);
    }
    fclose(outstr);

    return SUCCESS;
}
*/

#ifdef NONORMHIST

// Saves matrix ZetaM for each m multipole
local int PrintHistZetaM_sincos_normalized(struct  cmdline_data* cmd,
                                           struct  global_data* gd,
                                           gdlptr_sincos_omp_kkk gdl,
                                           gdlptr_sincos_omp_kkk_N gdlN)
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
    for (m=1; m<=cmd->mChebyshev+1; m++) {
        sprintf(namebuf, "%s%s_%d%s", gd->fpfnamehistZetaMFileName,
                "_sincos", m, EXTFILES);
        verb_print_q(2, cmd->verbose, "Printing : to a file %s ...\n",namebuf);
        outstr = stropen(namebuf, "w!");
        for (n1=1; n1<=cmd->sizeHistN; n1++) {
            for (n2=1; n2<=cmd->sizeHistN; n2++) {
                fprintf(outstr,"%16.8e ",
                gdl->histZetaMsincos[m][n1][n2]/gdlN->histZetaMcos[1][n1][n2]);
            }
            fprintf(outstr,"\n");
        }
        fclose(outstr);
    }
    // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
    for (m=1; m<=cmd->mChebyshev+1; m++) {
        sprintf(namebuf, "%s%s_%d%s", gd->fpfnamehistZetaMFileName,
                "_cossin", m, EXTFILES);
        verb_print_q(2, cmd->verbose, "Printing : to a file %s ...\n",namebuf);
        outstr = stropen(namebuf, "w!");
        for (n1=1; n1<=cmd->sizeHistN; n1++) {
            for (n2=1; n2<=cmd->sizeHistN; n2++) {
                fprintf(outstr,"%16.8e ",
                gdl->histZetaMcossin[m][n1][n2]/gdlN->histZetaMcos[1][n1][n2]);
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
                                            gdlptr_sincos_omp_kkk gdl,
                                            gdlptr_sincos_omp_kkk_N gdlN)
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
                    rbinlog = rlog10(cmd->rminHist) + ((real)(n1)-0.5)*gd->deltaR;
                }
                rBin=rpow(10.0,rbinlog);
            } else {
                rBin = cmd->rminHist + ((real)n1-0.5)*gd->deltaR;
            }
            Zeta  = gdl->histZetaMcos[m][n1][n1]
                  /
                    (
                     gdlN->histZetaMcos[1][n1][n1]
                     + gdlN->histZetaMsin[1][n1][n1]
                     );
            Zeta2 = gdl->histZetaMcos[m][n1][(int)(Nbins/4.0)]
                  /
                    (
                     gdlN->histZetaMcos[1][n1][(int)(Nbins/4.0)]
                     + gdlN->histZetaMsin[1][n1][(int)(Nbins/4.0)]
                     );
            Zeta3 = gdl->histZetaMcos[m][n1][(int)(2.0*Nbins/4.0)]
                  /
                    (
                    gdlN->histZetaMcos[1][n1][(int)(2.0*Nbins/4.0)]
                     + gdlN->histZetaMsin[1][n1][(int)(2.0*Nbins/4.0)]
                     );
            Zeta4 = gdl->histZetaMcos[m][n1][(int)(3.0*Nbins/4.0)]
                  /
                    (
                    gdlN->histZetaMcos[1][n1][(int)(3.0*Nbins/4.0)]
                     + gdlN->histZetaMsin[1][n1][(int)(3.0*Nbins/4.0)]
                     );
            Zeta5 = gdl->histZetaMcos[m][n1][(int)(4.0*Nbins/4.0 - 1.0)]
                  /
                    (
                    gdlN->histZetaMcos[1][n1][(int)(4.0*Nbins/4.0 - 1.0)]
                     + gdlN->histZetaMsin[1][n1][(int)(4.0*Nbins/4.0 - 1.0)]
                     );
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
                    rbinlog = rlog10(cmd->rminHist) + ((real)(n1)-0.5)*gd->deltaR;
                }
                rBin=rpow(10.0,rbinlog);
            } else {
                rBin = cmd->rminHist + ((real)n1-0.5)*gd->deltaR;
            }
                Zeta = gdl->histZetaMsin[m][n1][n1]
                     /
                        (
                        gdlN->histZetaMcos[1][n1][n1]
                         + gdlN->histZetaMsin[1][n1][n1]
                         );
                Zeta2 = gdl->histZetaMsin[m][n1][(int)(Nbins/4.0)]
                        /
                        (
                        gdlN->histZetaMcos[1][n1][(int)(Nbins/4.0)]
                         + gdlN->histZetaMsin[1][n1][(int)(Nbins/4.0)]
                         );
                Zeta3 = gdl->histZetaMsin[m][n1][(int)(2.0*Nbins/4.0)]
                        /
                        (
                         gdlN->histZetaMcos[1][n1][(int)(2.0*Nbins/4.0)]
                         + gdlN->histZetaMsin[1][n1][(int)(2.0*Nbins/4.0)]
                         );
                Zeta4 = gdl->histZetaMsin[m][n1][(int)(3.0*Nbins/4.0)]
                        /
                        (
                         gdlN->histZetaMcos[1][n1][(int)(3.0*Nbins/4.0)]
                         + gdlN->histZetaMsin[1][n1][(int)(3.0*Nbins/4.0)]
                         );
                Zeta5 = gdl->histZetaMsin[m][n1][(int)(4.0*Nbins/4.0 - 1.0)]
                        /
                        (
                         gdlN->histZetaMcos[1][n1][(int)(4.0*Nbins/4.0 - 1.0)]
                         + gdlN->histZetaMsin[1][n1][(int)(4.0*Nbins/4.0 - 1.0)]
                         );
                fprintf(outstr,MHISTZETA,rBin,Zeta,Zeta2,Zeta3,Zeta4,Zeta5);
        }
        fclose(outstr);
    }

    for (m = 1; m <= cmd->mChebyshev+1; m++) {
        sprintf(namebuf, "%s%s_%d%s", gd->fpfnamemhistZetaMFileName,
                "_sincos", m, EXTFILES);
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
            Zeta = gdl->histZetaMsincos[m][n1][n1]
                    /
                    (
                    gdlN->histZetaMcos[1][n1][n1]
                     + gdlN->histZetaMsin[1][n1][n1]
                     );
            Zeta2 = gdl->histZetaMsincos[m][n1][(int)(Nbins/4.0)]
                    /
                    (
                    gdlN->histZetaMcos[1][n1][(int)(Nbins/4.0)]
                     + gdlN->histZetaMsin[1][n1][(int)(Nbins/4.0)]
                     );
            Zeta3 = gdl->histZetaMsincos[m][n1][(int)(2.0*Nbins/4.0)]
                    /
                    (
                     gdlN->histZetaMcos[1][n1][(int)(2.0*Nbins/4.0)]
                     + gdlN->histZetaMsin[1][n1][(int)(2.0*Nbins/4.0)]
                     );
            Zeta4 = gdl->histZetaMsincos[m][n1][(int)(3.0*Nbins/4.0)]
                    /
                    (
                     gdlN->histZetaMcos[1][n1][(int)(3.0*Nbins/4.0)]
                     + gdlN->histZetaMsin[1][n1][(int)(3.0*Nbins/4.0)]
                     );
            Zeta5 = gdl->histZetaMsincos[m][n1][(int)(4.0*Nbins/4.0 - 1.0)]
                    /
                    (
                     gdlN->histZetaMcos[1][n1][(int)(4.0*Nbins/4.0 - 1.0)]
                     + gdlN->histZetaMsin[1][n1][(int)(4.0*Nbins/4.0 - 1.0)]
                     );

            fprintf(outstr,MHISTZETA,rBin,Zeta,Zeta2,Zeta3,Zeta4,Zeta5);
        }
        fclose(outstr);
    }

    // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
    for (m = 1; m <= cmd->mChebyshev+1; m++) {
        sprintf(namebuf, "%s%s_%d%s", gd->fpfnamemhistZetaMFileName,
                "_cossin", m, EXTFILES);
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
            Zeta = gdl->histZetaMcossin[m][n1][n1]
                    /
                    (
                    gdlN->histZetaMcos[1][n1][n1]
                     + gdlN->histZetaMsin[1][n1][n1]
                     );
            Zeta2 = gdl->histZetaMcossin[m][n1][(int)(Nbins/4.0)]
                    /
                    (
                    gdlN->histZetaMcos[1][n1][(int)(Nbins/4.0)]
                     + gdlN->histZetaMsin[1][n1][(int)(Nbins/4.0)]
                     );
            Zeta3 = gdl->histZetaMcossin[m][n1][(int)(2.0*Nbins/4.0)]
                    /
                    (
                     gdlN->histZetaMcos[1][n1][(int)(2.0*Nbins/4.0)]
                     + gdlN->histZetaMsin[1][n1][(int)(2.0*Nbins/4.0)]
                     );
            Zeta4 = gdl->histZetaMcossin[m][n1][(int)(3.0*Nbins/4.0)]
                    /
                    (
                     gdlN->histZetaMcos[1][n1][(int)(3.0*Nbins/4.0)]
                     + gdlN->histZetaMsin[1][n1][(int)(3.0*Nbins/4.0)]
                     );
            Zeta5 = gdl->histZetaMcossin[m][n1][(int)(4.0*Nbins/4.0 - 1.0)]
                    /
                    (
                     gdlN->histZetaMcos[1][n1][(int)(4.0*Nbins/4.0 - 1.0)]
                     + gdlN->histZetaMsin[1][n1][(int)(4.0*Nbins/4.0 - 1.0)]
                     );

            fprintf(outstr,MHISTZETA,rBin,Zeta,Zeta2,Zeta3,Zeta4,Zeta5);
        }
        fclose(outstr);
    }

    return SUCCESS;
}


// Saves matrix ZetaM for each m multipole
local int PrintHistZetaM_sincos_edge_effects(struct  cmdline_data* cmd,
                                              struct  global_data* gd,
                                              gdlptr_sincos_omp_kkk gdl,
                                              gdlptr_sincos_omp_kkk_N gdlN)
{
    int n1, n2, m;
    stream outstr;
    char namebuf[256];
    real rBin, rbinlog;

    for (m=1; m<=cmd->mChebyshev+1; m++) {
        sprintf(namebuf, "%s%s_%d%s", gd->fpfnamehistZetaMFileName,
                "_EE", m, EXTFILES);
        verb_print_q(2, cmd->verbose, "Printing : to a file %s ...\n",namebuf);
        outstr = stropen(namebuf, "w!");
        for (n1=1; n1<=cmd->sizeHistN; n1++) {
            for (n2=1; n2<=cmd->sizeHistN; n2++) {
                matrixClm(cmd, gd, gdl, gdlN, n1, n2);
                fprintf(outstr,"%16.8e ",gdl->histZetaM[m][n1][n2]);
            }
            fprintf(outstr,"\n");
        }
        fclose(outstr);
    }

    //B  and saves matrix ZetaM for each m multipole at a set of theta2 angles
    if (scanopt(cmd->options, "out-m-HistZeta")) {
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
                Zeta = gdl->histZetaM[m][n1][n1];
                Zeta2 = gdl->histZetaM[m][n1][(int)(Nbins/4.0)];
                Zeta3 = gdl->histZetaM[m][n1][(int)(2.0*Nbins/4.0)];
                Zeta4 = gdl->histZetaM[m][n1][(int)(3.0*Nbins/4.0)];
                Zeta5 = gdl->histZetaM[m][n1][(int)(4.0*Nbins/4.0 - 1.0)];
                fprintf(outstr,MHISTZETA,rBin,Zeta,Zeta2,Zeta3,Zeta4,Zeta5);
            }
            fclose(outstr);
        }
    }
    //E

    return SUCCESS;
}

//B correction 2025-04-06
// Saves matrix ZetaM for each m multipole
//  Normalized matrices:
//  ZetaM_cos, ZetaM_sin, Zeta_cossin, Zeta_sincos
//  NZetaM_cos, NZetaM_sin, NZeta_cossin, NZeta_sincos
local int PrintHistZetaM_sincos_edge_effects_normalized_HistZeta(
                                            struct  cmdline_data* cmd,
                                            struct  global_data* gd,
                                            gdlptr_sincos_omp_kkk gdl,
                                            gdlptr_sincos_omp_kkk_N gdlN)
{
    int n1, n2, m;
    stream outstr;
    char namebuf[256];
    real rBin, rbinlog;

    for (m=1; m<=cmd->mChebyshev+1; m++) {
        sprintf(namebuf, "%s%s_%d%s", gd->fpfnamehistZetaMFileName,
                "_EE", m, EXTFILES);
        verb_print_q(2, cmd->verbose, "Printing : to a file %s ...\n",namebuf);
        outstr = stropen(namebuf, "w!");
        for (n1=1; n1<=cmd->sizeHistN; n1++) {
            for (n2=1; n2<=cmd->sizeHistN; n2++) {
                matrixClm_normalized_HistZeta(cmd, gd, gdl, gdlN, n1, n2);
                fprintf(outstr,"%16.8e ",gdl->histZetaM[m][n1][n2]);
            }
            fprintf(outstr,"\n");
        }
        fclose(outstr);
    }

    //B  and saves matrix ZetaM for each m multipole at a set of theta2 angles
    if (scanopt(cmd->options, "out-m-HistZeta")) {
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
                Zeta = gdl->histZetaM[m][n1][n1];
                Zeta2 = gdl->histZetaM[m][n1][(int)(Nbins/4.0)];
                Zeta3 = gdl->histZetaM[m][n1][(int)(2.0*Nbins/4.0)];
                Zeta4 = gdl->histZetaM[m][n1][(int)(3.0*Nbins/4.0)];
                Zeta5 = gdl->histZetaM[m][n1][(int)(4.0*Nbins/4.0 - 1.0)];
                fprintf(outstr,MHISTZETA,rBin,Zeta,Zeta2,Zeta3,Zeta4,Zeta5);
            }
            fclose(outstr);
        }
    }
    //E

    return SUCCESS;
}
//E correction 2025-04-06


//B correction 2025-04-06
/*
// Saves matrix ZetaM for each m multipole at a set of theta2 angles
local int PrintHistZetaMm_sincos_edge_effects(struct  cmdline_data* cmd,
                                              struct  global_data* gd,
                                              gdlptr_sincos_omp_kkk gdl,
                                              gdlptr_sincos_omp_kkk_N gdlN)
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
                "_EE", m, EXTFILES);
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
            matrixClm(cmd, gd, gdl, gdlN, n1, n1);
            Zeta = gdl->histZetaM[m][n1][n1];
            matrixClm(cmd, gd, gdl, gdlN, n1, (int)(Nbins/4.0));
            Zeta2 = gdl->histZetaM[m][n1][(int)(Nbins/4.0)];
            matrixClm(cmd, gd, gdl, gdlN, n1, (int)(2.0*Nbins/4.0));
            Zeta3 = gdl->histZetaM[m][n1][(int)(2.0*Nbins/4.0)];
            matrixClm(cmd, gd, gdl, gdlN, n1, (int)(3.0*Nbins/4.0));
            Zeta4 = gdl->histZetaM[m][n1][(int)(3.0*Nbins/4.0)];
            matrixClm(cmd, gd, gdl, gdlN, n1, (int)(4.0*Nbins/4.0 - 1.0));
            Zeta5 = gdl->histZetaM[m][n1][(int)(4.0*Nbins/4.0 - 1.0)];
            fprintf(outstr,MHISTZETA,rBin,Zeta,Zeta2,Zeta3,Zeta4,Zeta5);
        }
        fclose(outstr);
    }

    return SUCCESS;
}
*/
//E correction 2025-04-06

#endif // ! NONORMHIST

#endif // ! NMultipoles

#undef MHISTZETAHEADER
#undef MHISTZETA


//E Saving histograms section: case KKKCORRELATION:
