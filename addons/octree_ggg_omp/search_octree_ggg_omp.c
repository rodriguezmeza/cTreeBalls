/* ==============================================================================
 MODULE: search_octree_ggg_omp.c		[cTreeBalls]
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

// edge-corrections in the pivot loop
#define PIVOTLOOP
#undef PIVOTLOOP

#ifdef THREEPCFSHEAR
#define mpOffSet        3
int local mCheb;                                    // mCheb =
                                                    // cmd->mChebyshev + mpOffSet
#endif

#ifdef THREEPCFCONVERGENCE
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
#endif

#ifdef NMultipoles
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

#ifdef THREEPCFSHEAR
//B shear macro definition
//B Macro for any posible value of mChebyshev
//  for recursivity needs that at least 3 multipoles be evaluated
#define CHEBYSHEVTUOMPGGGANY                                        \
{real g1cosmphi,g1sinmphi,g2cosmphi,g2sinmphi; int m;               \
    hist->ChebsT[1] = 1.0;                                          \
    g1cosmphi = gamma1 * hist->ChebsT[1];                           \
    g2cosmphi = gamma2 * hist->ChebsT[1];                           \
    hist->histg1threadcos[1][n] += g1cosmphi;                       \
    hist->histg2threadcos[1][n] += g2cosmphi;                       \
    hist->ChebsT[2] = 1.0;                                          \
    g1cosmphi = gamma1 * hist->ChebsT[2];                           \
    g2cosmphi = gamma2 * hist->ChebsT[2];                           \
    hist->histg1threadcos[2][n] += g1cosmphi;                       \
    hist->histg2threadcos[2][n] += g2cosmphi;                       \
    hist->ChebsT[3] = 1.0;                                          \
    g1cosmphi = gamma1 * hist->ChebsT[3];                           \
    g2cosmphi = gamma2 * hist->ChebsT[3];                           \
    hist->histg1threadcos[3][n] += g1cosmphi;                       \
    hist->histg2threadcos[3][n] += g2cosmphi;                       \
    hist->ChebsU[1] = 0.0;                                          \
    g1sinmphi = gamma1 * hist->ChebsU[1] * sinphi;                  \
    g2sinmphi = gamma2 * hist->ChebsU[1] * sinphi;                  \
    hist->histg1threadsin[1][n] += g1sinmphi;                       \
    hist->histg2threadsin[1][n] += g2sinmphi;                       \
    hist->ChebsU[2] = 0.0;                                          \
    g1sinmphi = gamma1 * hist->ChebsU[2] * sinphi;                  \
    g2sinmphi = gamma2 * hist->ChebsU[2] * sinphi;                  \
    hist->histg1threadsin[2][n] += g1sinmphi;                       \
    hist->histg2threadsin[2][n] += g2sinmphi;                       \
    hist->ChebsU[3] = 0.0;                                          \
    g1sinmphi = gamma1 * hist->ChebsU[3] * sinphi;                  \
    g2sinmphi = gamma2 * hist->ChebsU[3] * sinphi;                  \
    hist->histg1threadsin[3][n] += g1sinmphi;                       \
    hist->histg2threadsin[3][n] += g2sinmphi;                      \
    for (m=4; m<=cmd->mChebyshev+1; m++){                           \
        hist->ChebsT[m] = 2.0*(cosphi)*hist->ChebsT[m-1]-hist->ChebsT[m-2]; \
        g1cosmphi = gamma1 * hist->ChebsT[m];                       \
        g2cosmphi = gamma2 * hist->ChebsT[m];                       \
        hist->histg1threadcos[m][n] += g1cosmphi;                   \
        hist->histg2threadcos[m][n] += g2cosmphi;                   \
        hist->ChebsU[m] = 2.0*(cosphi)*hist->ChebsU[m-1]-hist->ChebsU[m-2]; \
        g1sinmphi = gamma1 * hist->ChebsU[m] * sinphi;              \
        g2sinmphi = gamma2 * hist->ChebsU[m] * sinphi;              \
        hist->histg1threadsin[m][n] += g1sinmphi;                   \
        hist->histg2threadsin[m][n] += g2sinmphi;                   \
    }}
//E
//E shear macro definition
#endif // ! THREEPCFSHEAR

//B Define structures:
typedef struct {
#ifdef TWOPCF
    real *histNN;                   // used
    real *histWW;                   // used
    real *histCF;                   // used
    real *histNNSubXi2pcf;          // used
    //B kappa Avg Rmin
    realptr histNNSubXi2pcftotal;   // used
    //E
    real *histXi2pcf;               // used
#endif

    realptr histNNSub;

#ifdef THREEPCFCONVERGENCE
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
#endif

#ifdef THREEPCFSHEAR
#endif

} gdl_sincos_omp_ggg, *gdlptr_sincos_omp_ggg;

typedef struct {
#ifdef TWOPCF
    realptr histNthread;            // used
    realptr histWthread;            // used
    realptr histWWthread;            // used
    realptr histNNSubXi2pcfthread;  // used
    //B kappa Avg Rmin
    realptr histNNSubXi2pcfthreadp; // used
    realptr histNNSubXi2pcfthreadtotal; // used
    //E
    real *histXi2pcfthread;         // used
    real *histXi2pcfthreadsub;      // used
#endif

    real *ChebsT;
    real *ChebsU;

    realptr histNNSubthread;

#ifdef THREEPCFCONVERGENCE
    real **xiOUTVPcos;
    real **xiOUTVPsin;
    real **xiOUTVPsincos;
    real **xiOUTVPcossin;
#ifdef PIVOTLOOP
    real ***histZetaMtmpcos;
    real ***histZetaMtmpsin;
    real ***histZetaMtmpsincos;
    real ***histZetaMtmpcossin;
#else
    real **histZetaMtmpcos;
    real **histZetaMtmpsin;
    real **histZetaMtmpsincos;
    real **histZetaMtmpcossin;
#endif
    real ***histZetaMthreadcos;
    real ***histZetaMthreadsin;
    real ***histZetaMthreadsincos;
    real ***histZetaMthreadcossin;

    real **histXithreadcos;
    real **histXithreadsin;
#ifdef PIVOTLOOP
    real ***histZetaMtmp;
#endif
#endif

#ifdef THREEPCFSHEAR
    real **histg1threadcos;
    real **histg1threadsin;
    real **histg2threadcos;
    real **histg2threadsin;
    real **histReGthread;
    real **histImGthread;
    real **histReGNthread;
    real **histImGNthread;
    real **ReUpsilonOUTVP0;
    real **ImUpsilonOUTVP0;
#endif

    INTEGER nbbcalcthread;
    INTEGER nbccalcthread;

    vector q0;
    real drpq2, drpq;
    vector dr0;
    real cosb;
    real sinb;
} gdhist_sincos_omp_ggg, *gdhistptr_sincos_omp_ggg;

#ifdef NMultipoles
//B estructures for NMultipoles should be the same as KKK multipoles
//  check the memory allocation and freeing is updated...
typedef struct {
 #ifdef TWOPCF
     real *histNN;                   // used
    real *histWW;                   // used
     real *histCF;                   // used
     real *histNNSubXi2pcf;          // used
     //B kappa Avg Rmin
     realptr histNNSubXi2pcftotal;   // used
     //E
     real *histXi2pcf;               // used
 #endif

     realptr histNNSub;

 #ifdef THREEPCFCONVERGENCE
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
 #endif

 #ifdef THREEPCFSHEAR
 #endif
} gdl_sincos_omp_ggg_N, *gdlptr_sincos_omp_ggg_N;

typedef struct {
#ifdef TWOPCF
    realptr histNthread;            // used
    realptr histWthread;            // used
    realptr histWWthread;            // used
    realptr histNNSubXi2pcfthread;  // used
    //B kappa Avg Rmin
    realptr histNNSubXi2pcfthreadp; // used
    realptr histNNSubXi2pcfthreadtotal; // used
    //E
    real *histXi2pcfthread;         // used
    real *histXi2pcfthreadsub;      // used
#endif

    real *ChebsT;
    real *ChebsU;

    realptr histNNSubthread;

#ifdef THREEPCFCONVERGENCE
    real **xiOUTVPcos;
    real **xiOUTVPsin;
    real **xiOUTVPsincos;
    real **xiOUTVPcossin;
#ifdef PIVOTLOOP
    real ***histZetaMtmpcos;
    real ***histZetaMtmpsin;
    real ***histZetaMtmpsincos;
    real ***histZetaMtmpcossin;
#else
    real **histZetaMtmpcos;
    real **histZetaMtmpsin;
    real **histZetaMtmpsincos;
    real **histZetaMtmpcossin;
#endif
    real ***histZetaMthreadcos;
    real ***histZetaMthreadsin;
    real ***histZetaMthreadsincos;
    real ***histZetaMthreadcossin;

    real **histXithreadcos;
    real **histXithreadsin;
#ifdef PIVOTLOOP
    real ***histZetaMtmp;
#endif
#endif

#ifdef THREEPCFSHEAR
    real **histg1threadcos;
    real **histg1threadsin;
    real **histg2threadcos;
    real **histg2threadsin;
    real **histReGthread;
    real **histImGthread;
    real **histReGNthread;
    real **histImGNthread;
#endif

    INTEGER nbbcalcthread;
    INTEGER nbccalcthread;

    vector q0;
    real drpq2, drpq;
    vector dr0;
    real cosb;
    real sinb;
} gdhist_sincos_omp_ggg_N, *gdhistptr_sincos_omp_ggg_N;
//E estructures for NMultipoles should be the same as KKK multipoles
#endif // ! NMultipoles
//E Define structures

local void normal_walktree_sincos(struct  cmdline_data* cmd,
                                  struct  global_data* gd,
                                  bodyptr *btable, int cat2,
                                  bodyptr, nodeptr, real,
                                  gdhistptr_sincos_omp_ggg, int *, int *);
local void sumnode_sincos(struct  cmdline_data* cmd,
                          struct  global_data* gd,
                          bodyptr *btable, int cat2,
                          bodyptr, cellptr, cellptr,
                          gdhistptr_sincos_omp_ggg, int *, int *);
local void sumnode_sincos_cell(struct  cmdline_data* cmd,
                               struct  global_data* gd,
                               bodyptr *btable, int cat2,
                               bodyptr p, cellptr start, cellptr finish,
                               gdhistptr_sincos_omp_ggg hist,
                               int *nbList, int *intList);

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

#ifdef TWOPCF
local int search_compute_Xi_ggg(struct  cmdline_data* cmd,
                                struct  global_data* gd, int nbody,
                                gdlptr_sincos_omp_ggg gdl);
local int search_compute_HistN_ggg(struct  cmdline_data* cmd,
                               struct  global_data* gd, int nbody,
                                   gdlptr_sincos_omp_ggg gdl);
local int PrintHistNN(struct cmdline_data* cmd, struct  global_data* gd,
                      gdlptr_sincos_omp_ggg gdl);
local int PrintHistCF(struct  cmdline_data* cmd, struct  global_data* gd,
                      gdlptr_sincos_omp_ggg gdl);
local int PrintHistXi2pcf(struct  cmdline_data* cmd, struct  global_data* gd,
                          gdlptr_sincos_omp_ggg gdl);
#endif

local int print_info(struct cmdline_data* cmd,
                     struct  global_data* gd);

#ifdef NMultipoles
local void normal_walktree_sincos_N(struct  cmdline_data* cmd,
                                    struct  global_data* gd,
                                    bodyptr *btable, int cat2,
                                    bodyptr, nodeptr, real,
                                    gdhistptr_sincos_omp_ggg,
                                    gdhistptr_sincos_omp_ggg_N,  int *, int *);
local void sumnode_sincos_N(struct  cmdline_data* cmd,
                            struct  global_data* gd,
                            bodyptr *btable, int cat2,
                            bodyptr, cellptr, cellptr,
                            gdhistptr_sincos_omp_ggg,
                            gdhistptr_sincos_omp_ggg_N, int *, int *);
local void sumnode_sincos_cell_N(struct  cmdline_data*,
                                 struct  global_data*,
                                 bodyptr *btable, int cat2,
                                 bodyptr, cellptr, cellptr,
                                 gdhistptr_sincos_omp_ggg,
                                 gdhistptr_sincos_omp_ggg_N,
                                 int *nbList, int *intList);
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
#endif


//B Saving histograms section: case KKKCORRELATION:
local int PrintHistrBins(struct  cmdline_data* cmd, struct  global_data* gd);
local int PrintHistZetaM_sincos(struct  cmdline_data* cmd,
                                struct  global_data* gd,
                                gdlptr_sincos_omp_ggg);
local int PrintHistZetaMm_sincos(struct  cmdline_data* cmd,
                               struct  global_data* gd,
                                 gdlptr_sincos_omp_ggg);
local int PrintHistZetaG(struct  cmdline_data* cmd,
                         struct  global_data* gd,
                         gdlptr_sincos_omp_ggg);
local int PrintHistZetaGm_sincos(struct  cmdline_data* cmd,
                                 struct  global_data* gd,
                                 gdlptr_sincos_omp_ggg);
local int PrintHistZetaMZetaGm_sincos(struct  cmdline_data* cmd,
                                      struct  global_data* gd,
                                      gdlptr_sincos_omp_ggg);
//B POLARAXIS
//      check NDIM and periodic...
local int polarix_init(struct  cmdline_data* cmd,
                               struct  global_data* gd,
                               bodyptr p, gdhistptr_sincos_omp_ggg hist);

local int polarix_init(struct  cmdline_data* cmd,
                               struct  global_data* gd,
                               bodyptr p, gdhistptr_sincos_omp_ggg hist)
{
    //B Set reference axis...
#ifdef POLARAXIS
    //B check NDIM and periodic...
    hist->q0[0] = 0.0;
    hist->q0[1] = 0.0;
    hist->q0[2] = 1.0;
    DOTPSUBV(hist->drpq2, hist->dr0, Pos(p), hist->q0);
    hist->drpq = rsqrt(hist->drpq2);
    real b = 2.0*rasin(hist->drpq/2.0);
    hist->cosb = rcos(b);
    hist->sinb = rsin(b);
    if (hist.drpq2==0) continue;
    //E
#else // ! POLARAXIS
    //B check NDIM and periodic...
    dRotation3D(Pos(p), ROTANGLE, ROTANGLE, ROTANGLE, hist->q0);
    DOTPSUBV(hist->drpq2, hist->dr0, Pos(p), hist->q0);
#ifdef SINGLEP
    hist->drpq = sqrt(hist->drpq2);
#else
    hist->drpq = rsqrt(hist->drpq2);
#endif
    //E
#endif // ! POLARAXIS
    //E

    return SUCCESS;
}

#ifdef POLARAXIS
#if NDIM == 3
#ifdef NOLIMBER
#define POLARAXIS_MAIN                              \
{                                                   \
    real a, c, c2;                                  \
    vector vc;                                      \
    DOTPSUBV(c2, vc, Pos(q), hist->q0);             \
    a = 2.0*rasin(dr1/2.0);                         \
    real cosc;                                      \
    cosc = Pos(q)[2];                               \
    cosphi = (cosc - (1.0-0.5*rsqr(a))*hist->cosb)  \
              /(a*hist->sinb);                      \
    if (rabs(cosphi) <= 1.0)                        \
        sinphi = rsqrt(1.0 - rsqr(cosphi));         \
    else                                            \
        sinphi = 0.0;                               \
    if (!crossVecProdSign(Pos(p), hist->q0, Pos(q)))\
        sinphi *= -1.0;                             \
}
#else // ! NOLIMBER
#define POLARAXIS_MAIN                              \
{                                                   \
    real a, c, c2;                                  \
    vector vc;                                      \
    DOTPSUBV(c2, vc, Pos(q), hist->q0);             \
    a = dr1;                                        \
    real cosc;                                      \
    cosc = Pos(q)[2];                               \
    cosphi = (cosc - (1.0-0.5*rsqr(a))*hist->cosb)  \
              /(a*hist->sinb);                      \
    if (rabs(cosphi) <= 1.0)                        \
        sinphi = rsqrt(1.0 - rsqr(cosphi));         \
    else                                            \
        sinphi = 0.0;                               \
    if (!crossVecProdSign(Pos(p), hist->q0, Pos(q)))\
        sinphi *= -1.0;                             \
}
#endif // ! NOLIMBER
#else // ! NDIM == 3
// work to do in 2D....
#endif // ! NDIM == 3

#else // ! POLARAXIS

#if NDIM == 3
#ifdef SINGLEP
#define POLARAXIS_MAIN                              \
{                                                   \
    float s, sy; float pr0[NDIM];                   \
    DOTVP(s, dr, hist->dr0);                        \
    cosphi = s/(dr1*hist->drpq);                    \
    CROSSVP(pr0,hist->dr0,Pos(p));                  \
    DOTVP(sy, dr, pr0);                             \
    if (rabs(cosphi)>1.0)                           \
        sinphi = 0.0;                               \
    else                                            \
        sinphi = sqrt(1.0 - cosphi*cosphi);         \
    if (sy < 0) sinphi *= -1.0;                     \
    if (cosphi>1.0) cosphi = 1.0;                   \
    if (cosphi<-1.0) cosphi = -1.0;                 \
}
#else // ! SINGLEP
//B DEFINITION USE BY DEFAULT...
#define POLARAXIS_MAIN                              \
{                                                   \
    real s, sy; vector pr0;                         \
    DOTVP(s, dr, hist->dr0);                        \
    cosphi = s/(dr1*hist->drpq);                    \
    CROSSVP(pr0,hist->dr0,Pos(p));                  \
    DOTVP(sy, dr, pr0);                             \
    sinphi = rsqrt(1.0 - rsqr(cosphi));             \
    if (sy < 0) sinphi *= -1.0;                     \
    if (rabs(cosphi)>1.0)                           \
    verb_log_print(cmd->verbose, gd->outlog,        \
    "sumenode: Warning!... cossphi must be in (-1,1): %g\n", \
                    cosphi);                        \
}
//E
#endif // ! SINGLEP
#else // ! NDIM == 3 ... 2 dimensions...
#ifdef SINGLEP
#define POLARAXIS_MAIN                              \
{                                                   \
    cosphi = -dr[0]/dr1;                            \
    sinphi = -dr[1]/dr1;                            \
    if (cosphi>1.0) cosphi = 1.0;                   \
    if (cosphi<-1.0) cosphi = -1.0;                 \
}
#else // ! SINGLEP
#define POLARAXIS_MAIN                              \
{                                                   \
    cosphi = -dr[0]/dr1;                            \
    sinphi = -dr[1]/dr1;                            \
    if (rabs(cosphi)>1.0)                           \
    verb_log_print(cmd->verbose, gd->outlog,        \
    "sumenode: Warning!... cossphi must be in (-1,1): %g\n", \
                    cosphi);                        \
}
#endif // ! SINGLEP
#endif // ! NDIM == 3 ... 2 dimensions...

#endif // ! POLARAXIS
//E POLARAXIS


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
                                  gdlptr_sincos_omp_ggg_N);
local int PrintHistZetaMm_sincos_N(struct  cmdline_data* cmd,
                               struct  global_data* gd,
                                   gdlptr_sincos_omp_ggg_N);
#ifdef NONORMHIST
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
#ifdef PIVOTLOOP
local int HistZetaM_sincos_edge_effects(struct  cmdline_data*,
                                        struct  global_data*,
                                        gdhistptr_sincos_omp_ggg,
                                        gdhistptr_sincos_omp_ggg_N);
#endif
local int PrintHistZetaM_sincos_edge_effects(struct  cmdline_data*,
                                             struct  global_data*,
                                             gdlptr_sincos_omp_ggg,
                                             gdlptr_sincos_omp_ggg_N);

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
                              gdhistptr_sincos_omp_ggg hist, int);

#endif
//E

#include <pthread.h>

/*
 Search routine using octree method:

 To be called using: search=octree-ggg-omp

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
global int searchcalc_octree_ggg_omp(struct cmdline_data* cmd,
                                     struct  global_data* gd,
                                     bodyptr *btable, INTEGER *nbody,
                                     INTEGER ipmin, INTEGER *ipmax,
                                     int cat1, int cat2)
{
    string routineName = "searchcalc_octree_ggg_omp";
    double cpustart;
    gdl_sincos_omp_ggg gdl;
#ifdef NMultipoles
    gdl_sincos_omp_ggg_N gdlN;
#endif

    cpustart = CPUTIME;
    print_info(cmd, gd);
    if (cmd->useLogHist==FALSE &&
        (strcmp(cmd->searchMethod,"octree-ggg-omp") == 0))
//        error("%s: can´t have loghist false and octree-ggg-omp (%d %s)\n",
//              routineName, cmd->useLogHist, cmd->searchMethod);
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                           "%s: can´t have loghist false and octree-ggg-omp (%d %s)\n",
                            routineName, cmd->useLogHist, cmd->searchMethod);

#ifdef THREEPCFSHEAR
    mCheb =cmd->mChebyshev + mpOffSet;
#endif

//B kappa Avg Rmin
#ifdef DEBUG
    sprintf(pivotsfilePath,"%s/pivot_info%s.txt",gd->tmpDir,cmd->suffixOutFiles);
    if(!(outpivots=fopen(pivotsfilePath, "w")))
        error("\n%s: error opening file '%s' \n",
              routineName, pivotsfilePath);
#endif
//E

#ifdef OPENMPCODE
    ThreadCount(cmd, gd, nbody[cat1], cat1);
    pthread_t main_thread_id = pthread_self();
#else
#error `OPENMPMACHINE` is not defined. Switch it on in Makefile_settings
#endif

    search_init_gd_sincos_omp_ggg(cmd, gd, &gdl);
#ifdef NMultipoles
    search_init_gd_sincos_omp_ggg_N(cmd, gd, &gdlN);
#endif

    //B Mask correction
    INTEGER ipmask;
    if (scanopt(cmd->options, "read-mask"))
        ipmask=0;
    //E

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
                       "%s: actlenNb = %ld\n", routineName, actlenNb);
        activeNb = (int *) allocate(actlenNb * sizeof(int));
        //E
        //B Alloc memory for interaction lists
        int NrangeN;
        NrangeN = 0.25*nbody[cat2]*rsqr(cmd->rangeN);
        actlenInt = FACTIVEINT * NrangeN;
        verb_print(cmd->verbose, "- NrangeN and actlenInt: %d %d\n",
                   NrangeN, actlenInt);
        verb_log_print(cmd->verbose,gd->outlog,
                       "%s: actlenInt = %ld\n", routineName, actlenInt);
        activeInt = (int *) allocate(actlenInt * sizeof(int));
        //E
    }
#endif
    //E

#if defined(NMultipoles) && defined(NONORMHIST)
    if (scanopt(cmd->options, "patch-with-all")) {
        verb_print(cmd->verbose,
            "\n%s: total number of pixels to be pivots: %ld\n",
                   routineName, gd->pivotCount);
    }
#endif

    if (scanopt(cmd->options, "pivot-number")) {
        ipmax[cat1] = gd->pivotNumber;
        verb_print(cmd->verbose, "\n%s: pivot number: %ld\n",
                   routineName, gd->pivotNumber);
    }

    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                           "\n%s: Total allocated %g MByte storage so far.\n",
                           routineName, gd->bytes_tot*INMB);

//B sometimes when running happens a "Floating point exception: 8" erreor
//      indicates a division by zero error within a program.
//          numeric code "8" specifically points to integer division
//          by zero as the root cause
//    verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
//                        "\n%s: stepState: %ld\n", routineName, gd->stepState);
//E

    verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                        "\nRunning...\n - Completed pivot node:\n");

//
// Check that all posibilities are taken in to account...
//
#ifdef DEBUG
//B In this segment there must be NMultipoles acting...
#ifdef ADDPIVOTNEIGHBOURS
#pragma omp parallel default(none)                                          \
    shared(cmd,gd,btable,nbody,roottable,outpivots,                         \
           actlenNb,activeNb,actlenInt,activeInt,                           \
           ipmin,ipmax,cat1,cat2,ipfalse,ipmask,                            \
           icountNbRmin,icountNbRminOverlap, gdl)
#else // ! ADDPIVOTNEIGHBOURS
#pragma omp parallel default(none)                                          \
    shared(cmd,gd,btable,nbody,roottable,outpivots,                         \
           ipmin,ipmax,cat1,cat2,ipfalse,ipmask,                            \
           icountNbRmin,icountNbRminOverlap, gdl)
#endif // ! ADDPIVOTNEIGHBOURS
//E
#else // ! DEBUG

#ifdef ADDPIVOTNEIGHBOURS
//B In this segment there must be NMultipoles acting...
#pragma omp parallel default(none)                                          \
    shared(cmd,gd,btable,nbody,roottable,                                   \
           actlenNb,activeNb,actlenInt,activeInt,                           \
           ipmin,ipmax,cat1,cat2,ipfalse,ipmask,                            \
           icountNbRmin,icountNbRminOverlap, gdl)
//E
#else // ! ADDPIVOTNEIGHBOURS

#ifdef NMultipoles
#ifndef BALLS4SCANLEV
#pragma omp parallel default(none)                                          \
    shared(cmd,gd,btable,nbody,roottable,                                   \
           ipmin,ipmax,cat1,cat2,ipfalse,ipmask,                            \
           icountNbRmin,icountNbRminOverlap,                                \
           gdl, gdlN, main_thread_id)
#else // ! BALLS4SCANLEV
#pragma omp parallel default(none)                                          \
    shared(cmd,gd,btable,nbody,roottable,nodetablescanlevB4,                \
           ipmin,ipmax,cat1,cat2,ipfalse,ipmask,                            \
           icountNbRmin,icountNbRminOverlap,                                \
           gdl, gdlN, main_thread_id)
#endif // ! BALLS4SCANLEV
#else // ! NMultipoles

#ifndef BALLS4SCANLEV
#pragma omp parallel default(none)                                          \
    shared(cmd,gd,btable,nbody,roottable,                                   \
           ipmin,ipmax,cat1,cat2,ipfalse,ipmask,                            \
           icountNbRmin,icountNbRminOverlap,gdl, main_thread_id)
#else // ! BALLS4SCANLEV
#pragma omp parallel default(none)                                          \
    shared(cmd,gd,btable,nbody,roottable,nodetablescanlevB4,                \
           ipmin,ipmax,cat1,cat2,ipfalse,ipmask,                            \
           icountNbRmin,icountNbRminOverlap,gdl, main_thread_id)
#endif // ! BALLS4SCANLEV

#endif // ! NMultipoles
#endif // ! ADDPIVOTNEIGHBOURS
#endif // ! DEBUG
  {
      pthread_t current_thread_id = pthread_self();

#ifndef BALLS4SCANLEV
    bodyptr p;
#else
    nodeptr p;
#endif
    bodyptr q;
    INTEGER n, m, ip;
    INTEGER i;

    //B init:
    gdhist_sincos_omp_ggg hist;
    search_init_sincos_omp_ggg(cmd, gd, &hist);
#ifdef NMultipoles
    gdhist_sincos_omp_ggg_N histN;
    search_init_sincos_omp_ggg_N(cmd, gd, &histN);
#endif
      if (main_thread_id == current_thread_id)
          verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                "\n\t Total allocated %g MByte storage so far (threads).\n",
                gd->bytes_tot*INMB);
    //E

    INTEGER ipmaskthreads;
      if (scanopt(cmd->options, "read-mask"))
          ipmaskthreads = 0;

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
#ifndef BALLS4SCANLEV
      for (p = btable[cat1] + ipmin -1; p < btable[cat1] + ipmax[cat1]; p++) {
#else
      for (INTEGER i=0; i< gd->nnodescanlevTableB4[cat1]; i++) {
          p = nodetablescanlevB4[cat1][i];
#endif
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

          if (scanopt(cmd->options, "read-mask")) {
            if (Mask(p) == FALSE) {
                ipmaskthreads++;
                continue;
            }
          }

#if defined(NMultipoles) && defined(NONORMHIST)
          if (scanopt(cmd->options, "patch-with-all")) {
              if (UpdatePivot(p) == FALSE) {
                  continue;
              }
          }
#endif

//B segment to be included below...
#ifdef TWOPCF
          for (n = 1; n <= cmd->sizeHistN; n++) {
              hist.histXi2pcfthreadsub[n] = 0.0;
              //B kappa Avg Rmin
              hist.histNNSubXi2pcfthreadp[n] = 0.;
              //E
          }
#endif

          for (n = 1; n <= cmd->sizeHistN; n++)
              hist.histNNSubthread[n] = 0.0;

          //B Set histograms to zero for the pivot
#ifdef THREEPCFCONVERGENCE
          CLRM_ext_ext(hist.histXithreadcos,
                       cmd->mChebyshev+1, cmd->sizeHistN);
          CLRM_ext_ext(hist.histXithreadsin,
                       cmd->mChebyshev+1, cmd->sizeHistN);
#endif

#ifdef THREEPCFCONVERGENCE
#ifdef NMultipoles
          //B 3pcf convergence & shear counting
          for (n = 1; n <= cmd->sizeHistN; n++)
              histN.histNNSubthread[n] = 0.0;
          CLRM_ext_ext(histN.histXithreadcos,
                       cmd->mChebyshev+1, cmd->sizeHistN);
          CLRM_ext_ext(histN.histXithreadsin,
                       cmd->mChebyshev+1, cmd->sizeHistN);
          //E
#endif
#endif

#ifdef THREEPCFSHEAR
          CLRM_ext_ext(hist.histg1threadcos,
                       cmd->mChebyshev+1, cmd->sizeHistN);
          CLRM_ext_ext(hist.histg1threadsin,
                       cmd->mChebyshev+1, cmd->sizeHistN);
          CLRM_ext_ext(hist.histg2threadcos,
                       cmd->mChebyshev+1, cmd->sizeHistN);
          CLRM_ext_ext(hist.histg2threadsin,
                       cmd->mChebyshev+1, cmd->sizeHistN);
#endif
          //E

#if defined(THREEPCFCONVERGENCE) || defined(THREEPCFSHEAR)
          //B 3pcf convergence & shear
#ifndef BALLS4SCANLEV
          polarix_init(cmd, gd, p, &hist);
#else
          polarix_init(cmd, gd, (bodyptr)p, &hist);
#endif
          //E 3pcf convergence & shear
#endif
//E segment to be included below...

#ifdef NMultipoles
#ifndef BALLS4SCANLEV
          normal_walktree_sincos_N(cmd, gd, btable, cat2,
                                   p, ((nodeptr) roottable[cat2]),
                                   gd->rSizeTable[cat2], &hist, &histN,
                                   &nbList, &intList);

#ifdef TWOPCF
          //B kappa Avg Rmin
          for (n = 1; n <= cmd->sizeHistN; n++) {
              hist.histNNSubXi2pcfthreadp[n] =
                        ((real)NbRmin(p))*hist.histNNSubXi2pcfthreadp[n];
              hist.histNNSubXi2pcfthreadtotal[n] +=
                        hist.histNNSubXi2pcfthreadp[n];
// Check if these lines are needed!!!
              if (scanopt(cmd->options, "smooth-pivot"))
                  hist.histNNSubthread[n] =
                                    ((real)NbRmin(p))*hist.histNNSubthread[n];
          }
          //E
#endif

          computeBodyProperties_sincos_ggg(cmd, gd, p,
                                           ipmax[cat1]-ipmin+1, &hist);
          computeBodyProperties_sincos_ggg_N(cmd, gd, p,
                                             ipmax[cat1]-ipmin+1, &histN);
#else // ! BALLS4SCANLEV
          normal_walktree_sincos_N(cmd, gd, btable, cat2,
                                   (bodyptr)p, ((nodeptr) roottable[cat2]),
                                   gd->rSizeTable[cat2], &hist, &histN,
                                   &nbList, &intList);

#ifdef TWOPCF
          //B kappa Avg Rmin
          for (n = 1; n <= cmd->sizeHistN; n++) {
              hist.histNNSubXi2pcfthreadp[n] =
                        ((real)NbRmin(p))*hist.histNNSubXi2pcfthreadp[n];
              hist.histNNSubXi2pcfthreadtotal[n] +=
                        hist.histNNSubXi2pcfthreadp[n];
// Check if these lines are needed!!!
              if (scanopt(cmd->options, "smooth-pivot"))
                  hist.histNNSubthread[n] =
                                    ((real)NbRmin(p))*hist.histNNSubthread[n];
          }
          //E
#endif

//          computeBodyProperties_sincos_ggg(cmd, gd, (bodyptr)p,
//                                           ipmax[cat1]-ipmin+1, &hist);
//          computeBodyProperties_sincos_ggg_N(cmd, gd, (bodyptr)p,
//                                             ipmax[cat1]-ipmin+1, &histN);
          computeBodyProperties_sincos_ggg(cmd, gd, (bodyptr)p,
                                           gd->nnodescanlevTableB4[cat1], &hist);
          computeBodyProperties_sincos_ggg_N(cmd, gd, (bodyptr)p,
                                             gd->nnodescanlevTableB4[cat1], &histN);
#endif // ! BALLS4SCANLEV

#ifdef PIVOTLOOP
#if defined(NMultipoles) && defined(NONORMHIST)
    //  make this two questions more precise...
        if (scanopt(cmd->options, "no-normalize-HistZeta")) {
            if (scanopt(cmd->options, "edge-corrections")) {
                HistZetaM_sincos_edge_effects(cmd, gd, &hist, &histN);
            }
        }
#endif
#endif

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

                  for (n = 1; n <= cmd->sizeHistN; n++)
                      histN.histNNSubthread[n] = 0.0;
                  CLRM_ext_ext(histN.histXithreadcos,
                               cmd->mChebyshev+1, cmd->sizeHistN);
                  CLRM_ext_ext(histN.histXithreadsin,
                               cmd->mChebyshev+1, cmd->sizeHistN);
                  //E
//B Set reference axis...
                  //E 3pcf convergence & shear
#ifndef BALLS4SCANLEV
          polarix_init(cmd, gd, p, &hist);
#else
          polarix_init(cmd, gd, (bodyptr)p, &hist);
                  //E 3pcf convergence & shear
#endif
//E
                  sumnode_nblist_omp(cmd, gd, btable, ipmin, ipmax, cat1, cat2,
                                     q, &hist, intList);
                  computeBodyProperties_sincos_ggg(cmd, gd, q,
                                                   ipmax[cat1]-ipmin+1, &hist);
                  computeBodyProperties_sincos_ggg_N(cmd, gd, q,
                                                     ipmax[cat1]-ipmin+1, &histN);
              } // ! end i loop
          } // ! scanoption smooth-pivot
#endif // ! ADDPIVOTNEIGHBOURS

#else // ! NMultipoles
#ifndef BALLS4SCANLEV
          normal_walktree_sincos(cmd, gd, btable, cat2,
                                 p, ((nodeptr) roottable[cat2]),
                                 gd->rSizeTable[cat2], &hist, &nbList, &intList);

#ifdef TWOPCF
          //B kappa Avg Rmin
          for (n = 1; n <= cmd->sizeHistN; n++) {
              hist.histNNSubXi2pcfthreadp[n] =
                        ((real)NbRmin(p))*hist.histNNSubXi2pcfthreadp[n];
              hist.histNNSubXi2pcfthreadtotal[n] +=
                        hist.histNNSubXi2pcfthreadp[n];
// Check if these lines are needed!!!
              if (scanopt(cmd->options, "smooth-pivot"))
                  hist.histNNSubthread[n] =
                                    ((real)NbRmin(p))*hist.histNNSubthread[n];
          }
          //E
#endif

          computeBodyProperties_sincos_ggg(cmd, gd, p,
                                           ipmax[cat1]-ipmin+1, &hist);
#else // ! BALLS4SCANLEV
          normal_walktree_sincos(cmd, gd, btable, cat2,
                                 (bodyptr)p, ((nodeptr) roottable[cat2]),
                                 gd->rSizeTable[cat2], &hist, &nbList, &intList);

#ifdef TWOPCF
          //B kappa Avg Rmin
          for (n = 1; n <= cmd->sizeHistN; n++) {
              hist.histNNSubXi2pcfthreadp[n] =
                        ((real)NbRmin(p))*hist.histNNSubXi2pcfthreadp[n];
              hist.histNNSubXi2pcfthreadtotal[n] +=
                        hist.histNNSubXi2pcfthreadp[n];
// Check if these lines are needed!!!
              if (scanopt(cmd->options, "smooth-pivot"))
                  hist.histNNSubthread[n] =
                                    ((real)NbRmin(p))*hist.histNNSubthread[n];
          }
          //E
#endif

//          computeBodyProperties_sincos_ggg(cmd, gd, (bodyptr)p,
//                                           ipmax[cat1]-ipmin+1, &hist);
          computeBodyProperties_sincos_ggg(cmd, gd, (bodyptr)p,
                                           gd->nnodescanlevTableB4[cat1], &hist);
#endif // ! BALLS4SCANLEV

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
                  //E
//B Set reference axis...
                  //E 3pcf convergence & shear
#ifndef BALLS4SCANLEV
          polarix_init(cmd, gd, p, &hist);
#else
          polarix_init(cmd, gd, (bodyptr)p, &hist);
                  //E 3pcf convergence & shear
#endif
//E
                  sumnode_nblist_omp(cmd, gd, btable, ipmin, ipmax, cat1, cat2,
                                     q, &hist, intList);
                  computeBodyProperties_sincos_ggg(cmd, gd, q,
                                                   ipmax[cat1]-ipmin+1, &hist);
              } // ! end i loop
          } // ! scanoption smooth-pivot
#endif // ! ADDPIVOTNEIGHBOURS

#endif // ! NMultipoles

#ifndef BALLS4SCANLEV
          ip = p - btable[cat1] + 1;
#else
          ip = i+1;
#endif
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
          if (ip%gd->stepState == 0)
          verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog, ".");
          if (ip%cmd->stepState == 0)
          verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                              "%d ", ip);
      } // end do body p

          if (main_thread_id == current_thread_id)
    verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                              "\n", ip);

#pragma omp critical
    {
#ifdef TWOPCF
        for (n = 1; n <= cmd->sizeHistN; n++) {
            gdl.histNN[n] += hist.histNthread[n];
            gdl.histWW[n] += hist.histWWthread[n];
            gdl.histNNSub[n] += hist.histNNSubthread[n];
            gdl.histNNSubXi2pcf[n] += hist.histNNSubXi2pcfthread[n];
//B kappa Avg Rmin
            gdl.histNNSubXi2pcftotal[n] += hist.histNNSubXi2pcfthreadtotal[n];
//E

// Check this line and the histogram array histXi2pcfthread
//  in source search.c and correct if necessary
            gdl.histXi2pcf[n] += hist.histXi2pcfthread[n];
//            gdl.histXi2pcf[n] += hist.histXi2pcfthreadsub[n];
        }
#endif

#ifdef THREEPCFCONVERGENCE
        for (m=1; m<=cmd->mChebyshev+1; m++) {
            ADDM_ext(gdl.histZetaMcos[m],gdl.histZetaMcos[m],
                     hist.histZetaMthreadcos[m],cmd->sizeHistN);
            ADDM_ext(gdl.histZetaMsin[m],gdl.histZetaMsin[m],
                     hist.histZetaMthreadsin[m],cmd->sizeHistN);
            ADDM_ext(gdl.histZetaMsincos[m],gdl.histZetaMsincos[m],
                     hist.histZetaMthreadsincos[m],cmd->sizeHistN);
            ADDM_ext(gdl.histZetaMcossin[m],gdl.histZetaMcossin[m],
                     hist.histZetaMthreadcossin[m],cmd->sizeHistN);
#ifdef PIVOTLOOP
            ADDM_ext(gdl.histZetaM[m],gdl.histZetaM[m],
                     hist.histZetaMtmp[m],cmd->sizeHistN);
#endif
        }
#endif

        gd->nbbcalc += hist.nbbcalcthread;
        gd->nbccalc += hist.nbccalcthread;

#ifdef THREEPCFCONVERGENCE
#ifdef NMultipoles
        //B 3pcf convergence & shear
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
        //E 3pcf convergence & shear
#endif
#endif

        if (scanopt(cmd->options, "read-mask"))
            ipmask += ipmaskthreads;

        //B kappa Avg Rmin
        ipfalse += ipfalsethreads;
        icountNbRmin += icountNbRminthread;
        icountNbRminOverlap += icountNbRminOverlapthread;
        //E
    } // ! critical

#ifdef NMultipoles
    search_free_sincos_omp_ggg_N(cmd, gd, &histN);  // free memory
#endif
    search_free_sincos_omp_ggg(cmd, gd, &hist);     // free memory
  } // end pragma omp parallel

    //B end of completed pivot
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog, "\n\n");
    //E

    //B kappa Avg Rmin
    real xi, den, num;
    int mm;
    if (scanopt(cmd->options, "smooth-pivot")) {
#ifdef BALLS4SCANLEV
        num = (real)gd->nnodescanlevTableB4[cat1];
        den = (real)(gd->nnodescanlevTableB4[cat1]-ipfalse);
#else
        num = (real)nbody[cat1];
        den = (real)(nbody[cat1]-ipfalse);
#endif

#ifdef NONORMHIST
        xi = 1.0;
#else
#ifdef ADDPIVOTNEIGHBOURS
        xi = 1.0;
#else
        xi = num/den;
#endif
#endif // ! NONORMHIST
        verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                               "%s: p falses found = %ld and %e %e %e\n",
                               routineName, ipfalse, num, den, xi);
#ifdef THREEPCFCONVERGENCE
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
#endif

#ifdef THREEPCFCONVERGENCE
#ifdef NMultipoles
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
#endif
#endif

        verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                               "%s: p falses found = %ld\n", routineName, ipfalse);
        //B kappa Avg Rmin
        verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                               "%s: count NbRmin found = %ld\n",
                               routineName, icountNbRmin);
        verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                               "%s: count overlap found = %ld\n",
                               routineName, icountNbRminOverlap);
        //E

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
        verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                               "%s: p falses found = %ld\n",
                               routineName, ifalsecount);
        verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                               "%s: p true found = %ld\n",
                               routineName, itruecount);
        verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                               "%s: total = %ld\n",
                               routineName, itruecount+ifalsecount);

        //E
    } // ! smooth-pivot

    if (scanopt(cmd->options, "read-mask"))
        verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                               "%s: p masked found = %ld\n", routineName, ipmask);

#ifdef DEBUG
    fclose(outpivots);                              // Close file to debug pivots
#endif

#ifdef TWOPCF
    int nn;
    if (!scanopt(cmd->options, "asymmetric")) {
        for (nn = 1; nn <= cmd->sizeHistN; nn++) {
            if (cmd->verbose>3)
                printf("%d %e %e %e\n", nn,
                   gdl.histNNSubXi2pcf[nn], gdl.histNNSubXi2pcftotal[nn],
                       gdl.histNN[nn]);
            gdl.histXi2pcf[nn] /= 2.0;
            gdl.histNNSubXi2pcf[nn] /= 2.0;
            if (scanopt(cmd->options, "weights-norm")) {
//                gdl.histXi2pcf[nn] /= MAX(gdl.histNN[nn],1.0);
                gdl.histXi2pcf[nn] /= MAX(gdl.histWW[nn],1.0);
            } else {
                //B kappa Avg Rmin
                gdl.histNNSubXi2pcftotal[nn] /= 2.0;
                if (scanopt(cmd->options, "smooth-pivot")) {
                    gdl.histXi2pcf[nn] /= MAX(gdl.histNNSubXi2pcftotal[nn],1.0);
                } else {
                    gdl.histXi2pcf[nn] /= MAX(gdl.histNNSubXi2pcf[nn],1.0);
                }
                //E
            }
        }
    } else {
        for (nn = 1; nn <= cmd->sizeHistN; nn++) {
            if (cmd->verbose>3)
            printf(0,"%d %e %e\n", nn,
                   gdl.histNNSubXi2pcf[nn], gdl.histNNSubXi2pcftotal[nn]);
            if (scanopt(cmd->options, "smooth-pivot")) {
                gdl.histXi2pcf[nn] /= MAX(gdl.histNNSubXi2pcftotal[nn],1.0);
            } else {
                gdl.histXi2pcf[nn] /= MAX(gdl.histNNSubXi2pcf[nn],1.0);
            }
        }
    }

    if (scanopt(cmd->options, "compute-HistN")) {
        if (scanopt(cmd->options, "smooth-pivot")) {
            search_compute_HistN_ggg(cmd, gd, nbody[cat1]-ipfalse, &gdl);
        } else {
            search_compute_HistN_ggg(cmd, gd, nbody[cat1], &gdl);
        }
    }
#endif

// ===============================================
//B Saving histograms section: case GGGCORRELATION:
// ===============================================
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
            "\n\t%s: printing octree-ggg-omp method...\n\n", routineName);

#ifdef TWOPCF
    if (scanopt(cmd->options, "compute-HistN")) PrintHistNN(cmd, gd, &gdl);
    PrintHistXi2pcf(cmd, gd, &gdl);
#endif

    PrintHistrBins(cmd, gd);

#ifdef THREEPCFCONVERGENCE
#ifdef NMultipoles

#ifdef NONORMHIST
        if (scanopt(cmd->options, "no-normalize-HistZeta"))
            PrintHistZetaM_sincos(cmd, gd, &gdl);
        else
            PrintHistZetaM_sincos_normalized(cmd, gd, &gdl, &gdlN);
#else
        PrintHistZetaM_sincos(cmd, gd, &gdl);
#endif // ! NONORMHIST

        PrintHistZetaM_sincos_N(cmd, gd, &gdlN);

#else // ! NMultipoles
    PrintHistZetaM_sincos(cmd, gd, &gdl);
#endif // ! NMultipoles

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

            PrintHistZetaMm_sincos_N(cmd, gd, &gdlN);

#else // ! NMultipoles
            PrintHistZetaMm_sincos(cmd, gd, &gdl);
#endif // ! NMultipoles
        }

        if (scanopt(cmd->options, "out-HistZetaG")) {
            PrintHistZetaGm_sincos(cmd, gd, &gdl);
            PrintHistZetaG(cmd, gd, &gdl);
            PrintHistZetaMZetaGm_sincos(cmd, gd, &gdl);
        }

#ifdef NMultipoles
#ifdef NONORMHIST
    //  make this two questions more precise...
        if (scanopt(cmd->options, "no-normalize-HistZeta")) {
            if (scanopt(cmd->options, "edge-corrections")) {
                PrintHistZetaM_sincos_edge_effects(cmd, gd, &gdl, &gdlN);
            }
        }
#endif // ! NONORMHIST
#endif // ! NMultipoles
#endif // ! THREEPCFCONVERGENCE

    gd->flagPrint = FALSE;
// ===============================================
//E Saving histograms section: case GGGCORRELATION
// ===============================================

// ===============================================
//B Making histograms public (cballys PXD) section
// ===============================================
#ifdef PXD
    int m, n, n1, n2;
    for (n=1; n<=cmd->sizeHistN; n++) {
        gd->histNN[n] = gdl.histNN[n];
        gd->histCF[n] = gdl.histCF[n];
        gd->histXi2pcf[n] = gdl.histXi2pcf[n];
    }

#ifdef THREEPCFCONVERGENCE
    for (m=1; m<=cmd->mChebyshev+1; m++) {
        for (n1=1; n1<=cmd->sizeHistN; n1++) {
            for (n2=1; n2<=cmd->sizeHistN; n2++) {
                gd->histZetaMcos[m][n1][n2] = gdl.histZetaMcos[m][n1][n2];
                gd->histZetaMsin[m][n1][n2] = gdl.histZetaMsin[m][n1][n2];
                gd->histZetaMsincos[m][n1][n2] = gdl.histZetaMsincos[m][n1][n2];
                gd->histZetaMcossin[m][n1][n2] = gdl.histZetaMcossin[m][n1][n2];
            }
        }
    }
#endif
#endif
// ===============================================
//E Making histograms public (cballys PXD) section
// ===============================================

//B free memory
#ifdef NMultipoles
    search_free_gd_sincos_omp_ggg_N(cmd, gd, &gdlN);
#endif
    search_free_gd_sincos_omp_ggg(cmd, gd, &gdl);
//E

    gd->cpusearch = CPUTIME - cpustart;
    verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                        "\nGoing out: CPU time = %lf %s\n",
                        CPUTIME-cpustart, PRNUNITOFTIMEUSED);

    return SUCCESS;
}

local void normal_walktree_sincos(struct  cmdline_data* cmd, 
                                  struct  global_data* gd,
                                  bodyptr *btable, int cat2,
                                  bodyptr p, nodeptr q, real qsize,
                                  gdhistptr_sincos_omp_ggg hist,
                                  int *nbList, int *intList)
{
    nodeptr l;
#ifdef SINGLEP
    float dr1;
    float dr[NDIM];
#else
    real dr1;
    vector dr;
#endif

    if (Update(p)==FALSE) return;
    if ( ((nodeptr) p) != q ) {
        if (Type(q) == CELL) {
            if (!reject_cell(cmd, gd, (nodeptr)p, q, qsize)) {
                if (!scanopt(cmd->options, "no-one-ball")) {
                    accept_body(cmd, gd, p, (nodeptr)q, &dr1, dr);

#ifndef NORMALHISTSCALE
//B useLogHist section
                    if ( (Radius(p)+Radius(q))/(dr1) < gd->deltaR)
                        sumnode_sincos_cell(cmd, gd, btable, cat2, p,
                                            ((cellptr) q), ((cellptr) q+1), hist,
                                            nbList, intList);
                    else
                        for (l = More(q); l != Next(q); l = Next(l))
                            normal_walktree_sincos(cmd, gd, btable, cat2,
                                                   p,l,qsize/2, hist,
                                                   nbList, intList);
//E useLogHist section
#else // ! NORMALHISTSCALE
                    if ( (Radius(p)+Radius(q)) < gd->deltaR*THETA)
                        sumnode_sincos_cell(cmd, gd, btable, cat2, p,
                                            ((cellptr) q), ((cellptr) q+1), hist,
                                            nbList, intList);
                    else
                        for (l = More(q); l != Next(q); l = Next(l))
                            normal_walktree_sincos(cmd, gd, btable, cat2,
                                                   p,l,qsize/2, hist,
                                                   nbList, intList);
#endif // ! NORMALHISTSCALE

                } else { // ! no-one-ball
                    for (l = More(q); l != Next(q); l = Next(l))
                        normal_walktree_sincos(cmd, gd, btable, cat2,
                                               p,l,qsize/2, hist,
                                               nbList, intList);
                } // ! no-one-ball
            } // ! !reject_cell
        } else { // ! Type(q) == CELL
            sumnode_sincos(cmd, gd, btable, cat2,
                           p, ((cellptr)q), ((cellptr)q+1), hist,
                           nbList, intList);
        } // ! Type(q) == CELL
    } // ! p != q
}

#ifdef ADDPIVOTNEIGHBOURS
//B kappa Avg Rmin
local void sumnode_nblist_omp(struct cmdline_data* cmd,
                              struct  global_data* gd,
                              bodyptr *btable,
                              INTEGER ipmin, INTEGER *ipmax,
                              int cat1, int cat2,
                              bodyptr p,
                              gdhistptr_sincos_omp_ggg hist, int intList)
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

#ifdef THREEPCFSHEAR
#endif
#ifdef TWOPCF
#endif
                //B check NDIM and periodic...
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
                //E
                if (cmd->mChebyshev<7) {
                    CHEBYSHEVTUOMPSINCOSANY
                } else {
                    CHEBYSHEVTUOMP;
                }

#ifdef THREEPCFSHEAR
#endif
#ifdef TWOPCF
#endif

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
                          gdhistptr_sincos_omp_ggg hist,
                          int *nbList, int *intList)
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
    REAL cosphi,sinphi;
    int iq;
#ifdef THREEPCFSHEAR
    real gamma1, gamma2;
#endif

    q = start;
    if (scanopt(cmd->options, "read-mask"))
        if (Mask(q)==FALSE) return;

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


#ifndef NORMALHISTSCALE
//B useLogHist section
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

#ifdef TWOPCF
                hist->histNthread[n] = hist->histNthread[n] + 1.;
                hist->histWthread[n] = hist->histWthread[n] + 1.;
                hist->histNNSubXi2pcfthread[n] =
                                            hist->histNNSubXi2pcfthread[n] + 1.;
                //B kappa Avg Rmin
                hist->histNNSubXi2pcfthreadp[n] =
                                            hist->histNNSubXi2pcfthreadp[n] + 1.;
                //E
#endif

                hist->histNNSubthread[n] = hist->histNNSubthread[n] + 1.;

                // needs to multiply xi by Weight(p)
                xi = Weight(q)*Kappa(q);

#ifdef THREEPCFSHEAR
                gamma1 = Gamma1(q);
                gamma2 = Gamma2(q);
#endif
#if defined(THREEPCFCONVERGENCE) || defined(THREEPCFSHEAR)
                POLARAXIS_MAIN;
#endif
#ifdef THREEPCFCONVERGENCE
                if (cmd->mChebyshev<7) {
                    CHEBYSHEVTUOMPSINCOSANY;
                } else {
                    CHEBYSHEVTUOMP;
                }
#endif

#ifdef THREEPCFSHEAR
                CHEBYSHEVTUOMPGGGANY;
#endif
#ifdef TWOPCF
                hist->histXi2pcfthreadsub[n] += xi;
#endif
                hist->nbbcalcthread += 1;
            } // ! 1 < n < sizeHistN
        } // ! dr1>cmd->rminHist
//E useLogHist section
#else // ! NORMALHISTSCALE

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
            n = (int) ( (dr1-cmd->rminHist) * gd->i_deltaR) + 1;
            if (n<=cmd->sizeHistN && n>=1) {

#ifdef TWOPCF
                hist->histNthread[n] = hist->histNthread[n] + 1.;
                hist->histWthread[n] = hist->histWthread[n] + 1.;
                hist->histNNSubXi2pcfthread[n] =
                                            hist->histNNSubXi2pcfthread[n] + 1.;
                //B kappa Avg Rmin
                hist->histNNSubXi2pcfthreadp[n] =
                                            hist->histNNSubXi2pcfthreadp[n] + 1.;
                //E
#endif

                hist->histNNSubthread[n] = hist->histNNSubthread[n] + 1.;

                // needs to multiply xi by Weight(p)
                xi = Weight(q)*Kappa(q);

#ifdef THREEPCFSHEAR
                gamma1 = Gamma1(q);
                gamma2 = Gamma2(q);
#endif
#if defined(THREEPCFCONVERGENCE) || defined(THREEPCFSHEAR)
                POLARAXIS_MAIN;
#endif
#ifdef THREEPCFCONVERGENCE
                if (cmd->mChebyshev<7) {
                    CHEBYSHEVTUOMPSINCOSANY;
                } else {
                    CHEBYSHEVTUOMP;
                }
#endif

#ifdef THREEPCFSHEAR
                CHEBYSHEVTUOMPGGGANY;
#endif
#ifdef TWOPCF
                hist->histXi2pcfthreadsub[n] += xi;
#endif
                hist->nbbcalcthread += 1;
            } // ! 1 < n < sizeHistN
        } // ! dr1>cmd->rminHist

#endif // ! NORMALHISTSCALE

    } // ! accept_body
}

local void sumnode_sincos_cell(struct  cmdline_data* cmd,
                               struct  global_data* gd,
                               bodyptr *btable, int cat2,
                               bodyptr p, cellptr start, cellptr finish,
                               gdhistptr_sincos_omp_ggg hist,
                               int *nbList, int *intList)
{
    cellptr q;
//B implement these memory optimization and check its results
#ifdef SINGLEP
    float dr1;
    float dr[NDIM];
#else
    real dr1;
    vector dr;
#endif
//E
    int n;
    real xi;
    REAL cosphi,sinphi;
#ifdef THREEPCFSHEAR
    real gamma1, gamma2;
#endif

    q = start;

    if (scanopt(cmd->options, "read-mask"))
        if (Mask(q)==FALSE) return;

    if (accept_body(cmd, gd, p, (nodeptr)q, &dr1, dr)) {
        

#ifndef NORMALHISTSCALE
//B useLogHist section
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

#ifdef TWOPCF
                hist->histNthread[n] = hist->histNthread[n] +  Nb(q);
                hist->histWthread[n] = hist->histWthread[n] +  Nb(q);
                hist->histNNSubXi2pcfthread[n] =
                                            hist->histNNSubXi2pcfthread[n] + 1.0;
                //B kappa Avg Rmin
                hist->histNNSubXi2pcfthreadp[n] =
                                            hist->histNNSubXi2pcfthreadp[n] + 1.;
                //E
#endif

                hist->histNNSubthread[n] = hist->histNNSubthread[n] + 1.0;

//B needs to multiply xi by Weight(p) of the cell q
#ifdef NONORMHIST
                if (scanopt(cmd->options, "no-normalize-HistZeta")) {
                    xi = Nb(q)*Kappa(q);
                } else {
#ifdef KappaAvgON
                    xi = KappaAvg(q)/Nb(q);
#else
#ifdef NMultipole
                    xi = Kappa(q);
#else
                    xi = Nb(q)*Kappa(q);
#endif
#endif
                }
#else // ! NONORMHIST
#ifdef KappaAvgON
                xi = KappaAvg(q)/Nb(q);
#else
                xi = Kappa(q);
#endif
#endif // ! NONORMHIST
//E

#ifdef THREEPCFSHEAR
                gamma1 = Gamma1(q);
                gamma2 = Gamma2(q);
#endif
#if defined(THREEPCFCONVERGENCE) || defined(THREEPCFSHEAR)
                POLARAXIS_MAIN;
#endif
#ifdef THREEPCFCONVERGENCE
                if (cmd->mChebyshev<7) {
                    CHEBYSHEVTUOMPSINCOSANY
                } else {
                    CHEBYSHEVTUOMP;
                }
#endif

#ifdef THREEPCFSHEAR
                CHEBYSHEVTUOMPGGGANY;
#endif
#ifdef TWOPCF
                hist->histXi2pcfthreadsub[n] += xi;
#endif

                hist->nbccalcthread += 1;
            } // ! 1 < n < sizeHistN
        } // ! dr1 > rminHist
//E useLogHist section
#else // ! NORMALHISTSCALE

        if(dr1>cmd->rminHist) {
#ifdef ADDPIVOTNEIGHBOURS
            INTEGER iq;
            iq = (bodyptr)q-btable[cat2];
            activeInt[*intList]=iq;
            *intList +=1;
            if (*intList > actlenInt)
                error("intList: too many neighbors\n");
#endif
            n = (int) ( (dr1-cmd->rminHist) * gd->i_deltaR) + 1;
            if (n<=cmd->sizeHistN && n>=1) {

#ifdef TWOPCF
                hist->histNthread[n] = hist->histNthread[n] +  Nb(q);
                hist->histWthread[n] = hist->histWthread[n] +  Nb(q);
                hist->histNNSubXi2pcfthread[n] =
                                            hist->histNNSubXi2pcfthread[n] + 1.0;
                //B kappa Avg Rmin
                hist->histNNSubXi2pcfthreadp[n] =
                                            hist->histNNSubXi2pcfthreadp[n] + 1.;
                //E
#endif

                hist->histNNSubthread[n] = hist->histNNSubthread[n] + 1.0;

//B needs to multiply xi by Weight(p) of the cell q
#ifdef NONORMHIST
                if (scanopt(cmd->options, "no-normalize-HistZeta")) {
                    xi = Nb(q)*Kappa(q);
                } else {
#ifdef KappaAvgON
                    xi = KappaAvg(q)/Nb(q);
#else
#ifdef NMultipole
                    xi = Kappa(q);
#else
                    xi = Nb(q)*Kappa(q);
#endif
#endif
                }
#else // ! NONORMHIST
#ifdef KappaAvgON
                xi = KappaAvg(q)/Nb(q);
#else
                xi = Kappa(q);
#endif
#endif // ! NONORMHIST
//E

#ifdef THREEPCFSHEAR
                gamma1 = Gamma1(q);
                gamma2 = Gamma2(q);
#endif
#if defined(THREEPCFCONVERGENCE) || defined(THREEPCFSHEAR)
                POLARAXIS_MAIN;
#endif
#ifdef THREEPCFCONVERGENCE
                if (cmd->mChebyshev<7) {
                    CHEBYSHEVTUOMPSINCOSANY
                } else {
                    CHEBYSHEVTUOMP;
                }
#endif

#ifdef THREEPCFSHEAR
                CHEBYSHEVTUOMPGGGANY;
#endif
#ifdef TWOPCF
                hist->histXi2pcfthreadsub[n] += xi;
#endif

                hist->nbccalcthread += 1;
            } // ! 1 < n < sizeHistN
        } // ! dr1 > rminHist

#endif // ! NORMALHISTSCALE

    } // ! accept_body
}

local int computeBodyProperties_sincos_ggg(struct  cmdline_data* cmd,
                                            struct  global_data* gd,
                                            bodyptr p, int nbody,
                                            gdhistptr_sincos_omp_ggg hist)
{
    int n;
    int m;
    real xi;
#ifdef TWOPCF
    real xi_2p;
    real wi;
#endif

// check Weight factor... must be an average of Weights
#ifdef NONORMHIST

#ifdef BALLS4SCANLEV
    xi = Weight(p)*Kappa(p);                    // equiv to Nb*(Weight/Nb)*Kappa
    if (scanopt(cmd->options, "smooth-pivot")) {
        xi = Nb(p)*KappaRmin(p)/NbRmin(p);
    }
#else
    xi = Weight(p)*Kappa(p);
    if (scanopt(cmd->options, "smooth-pivot")) {
        xi = KappaRmin(p)/NbRmin(p);
    }
#endif

#ifdef TWOPCF
#ifdef BALLS4SCANLEV
    wi = Weight(p);
    xi_2p = (Weight(p)/Nb(p))*Kappa(p);
    if (scanopt(cmd->options, "smooth-pivot")) {
        xi_2p = KappaRmin(p);
    }
#else
    wi = Weight(p);
    xi_2p = Weight(p)*Kappa(p);
    if (scanopt(cmd->options, "smooth-pivot")) {
        xi_2p = KappaRmin(p);
    }
#endif // ! BALLS4SCANLEV
#endif

#else // ! NONORMHIST

#ifdef ADDPIVOTNEIGHBOURS
    //B not finished yet
    xi = Weight(p)*Kappa(p)/nbody;
    //E
#else

#ifdef BALLS4SCANLEV
//    xi = Weight(p)*Kappa(p)/nbody;               // equiv to Nb*(Weight/Nb)*Kappa
    xi = (Weight(p)/Nb(p))*Kappa(p)/nbody;               // equiv to Nb*(Weight/Nb)*Kappa
    if (scanopt(cmd->options, "smooth-pivot")) {
//        xi = Nb(p)*KappaRmin(p)/NbRmin(p)/nbody;
//        xi = Nb(p)*KappaRmin(p)/nbody;
        xi = KappaRmin(p)/NbRmin(p)/nbody;
    }
#else
    //B Not working yet
    xi = Weight(p)*Kappa(p)/nbody;
    if (scanopt(cmd->options, "smooth-pivot")) {
        xi = (KappaRmin(p)/NbRmin(p))/nbody;
    }
    //E
#endif

#ifdef TWOPCF
#ifdef BALLS4SCANLEV
    wi = Weight(p);
    xi_2p = (Weight(p)/Nb(p))*Kappa(p);
    if (scanopt(cmd->options, "smooth-pivot")) {
        xi_2p = KappaRmin(p);
    }
#else
    wi = Weight(p);
    xi_2p = Weight(p)*Kappa(p);
    if (scanopt(cmd->options, "smooth-pivot")) {
        xi_2p = (KappaRmin(p));
    }
#endif
#endif

#endif // ! ADDPIVOTNEIGHBOURS

#endif // ! NONORMHIST


#ifndef NONORMHIST
#ifdef THREEPCFCONVERGENCE
    for (m=1; m<=cmd->mChebyshev+1; m++) {
        for (n=1; n<=cmd->sizeHistN; n++) {
            hist->histXithreadcos[m][n] /= MAX(hist->histNNSubthread[n],1.0);
            hist->histXithreadsin[m][n] /= MAX(hist->histNNSubthread[n],1.0);
        }
    }
#endif
#endif

#ifdef THREEPCFCONVERGENCE
    for (m=1; m<=cmd->mChebyshev+1; m++) {
        OUTVP_ext(hist->xiOUTVPcos,
            hist->histXithreadcos[m], hist->histXithreadcos[m], cmd->sizeHistN);
        OUTVP_ext(hist->xiOUTVPsin,
            hist->histXithreadsin[m], hist->histXithreadsin[m],cmd->sizeHistN);
        OUTVP_ext(hist->xiOUTVPsincos,
            hist->histXithreadsin[m], hist->histXithreadcos[m],cmd->sizeHistN);
        OUTVP_ext(hist->xiOUTVPcossin,
            hist->histXithreadcos[m], hist->histXithreadsin[m],cmd->sizeHistN);
#ifdef PIVOTLOOP
        CLRM_ext(hist->histZetaMtmpcos[m], cmd->sizeHistN);
        CLRM_ext(hist->histZetaMtmpsin[m], cmd->sizeHistN);
        CLRM_ext(hist->histZetaMtmpsincos[m], cmd->sizeHistN);
        CLRM_ext(hist->histZetaMtmpcossin[m], cmd->sizeHistN);
        MULMS_ext(hist->histZetaMtmpcos[m],hist->xiOUTVPcos,xi,cmd->sizeHistN);
        MULMS_ext(hist->histZetaMtmpsin[m],hist->xiOUTVPsin,xi,cmd->sizeHistN);
        MULMS_ext(hist->histZetaMtmpsincos[m],hist->xiOUTVPsincos,xi,cmd->sizeHistN);
        MULMS_ext(hist->histZetaMtmpcossin[m],hist->xiOUTVPcossin,xi,cmd->sizeHistN);
        ADDM_ext(hist->histZetaMthreadcos[m],
            hist->histZetaMthreadcos[m],hist->histZetaMtmpcos[m],cmd->sizeHistN);
        ADDM_ext(hist->histZetaMthreadsin[m],
            hist->histZetaMthreadsin[m],hist->histZetaMtmpsin[m],cmd->sizeHistN);
        ADDM_ext(hist->histZetaMthreadsincos[m],
            hist->histZetaMthreadsincos[m],
            hist->histZetaMtmpsincos[m],cmd->sizeHistN);
        ADDM_ext(hist->histZetaMthreadcossin[m],
            hist->histZetaMthreadcossin[m],
            hist->histZetaMtmpcossin[m],cmd->sizeHistN);
#else
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
#endif
    }
#endif

#ifdef THREEPCFSHEAR
    for (m=1; m<=cmd->mChebyshev+1; m++) {
        MULVS_ext(hist->histg2threadsin[m],hist->histg2threadsin[m],
                  -1.0,cmd->sizeHistN);
        ADDV_ext(hist->histReGthread[m], hist->histg1threadcos[m],
            hist->histg2threadsin[m], cmd->sizeHistN);
        ADDV_ext(hist->histImGthread[m], hist->histg1threadsin[m],
            hist->histg2threadcos[m], cmd->sizeHistN);
    }
    for (m=2; m<=cmd->mChebyshev+1; m++) {
        MULVS_ext(hist->histg2threadsin[m],hist->histg2threadsin[m],
                  -1.0,cmd->sizeHistN);
        MULVS_ext(hist->histg1threadsin[m],hist->histg1threadsin[m],
                  -1.0,cmd->sizeHistN);
        ADDV_ext(hist->histReGNthread[m-1], hist->histg1threadcos[m],
            hist->histg2threadsin[m], cmd->sizeHistN);
        ADDV_ext(hist->histImGNthread[m-1], hist->histg1threadsin[m],
            hist->histg2threadcos[m], cmd->sizeHistN);
    }
#endif
#ifdef TWOPCF
    for (n=1; n<=cmd->sizeHistN; n++) {
        hist->histXi2pcfthread[n] += xi_2p*hist->histXi2pcfthreadsub[n];
        hist->histWWthread[n] += wi*hist->histWthread[n];
    }
#endif

    return SUCCESS;
}

#ifdef NMultipoles
local void normal_walktree_sincos_N(struct  cmdline_data* cmd,
                                    struct  global_data* gd,
                                    bodyptr *btable, int cat2,
                                    bodyptr p, nodeptr q, real qsize,
                                    gdhistptr_sincos_omp_ggg hist,
                                    gdhistptr_sincos_omp_ggg_N histN,
                                    int *nbList, int *intList)
{
    nodeptr l;
#ifdef SINGLEP
    float dr1;
    float dr[NDIM];
#else
    real dr1;
    vector dr;
#endif

    if (Update(p)==FALSE) return;
    if ( ((nodeptr) p) != q ) {
        if (Type(q) == CELL) {
            if (!reject_cell(cmd, gd, (nodeptr)p, q, qsize)) {
                if (!scanopt(cmd->options, "no-one-ball")) {
                    accept_body(cmd, gd, p, (nodeptr)q, &dr1, dr);

#ifndef NORMALHISTSCALE
//B useLogHist section
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
//E useLogHist section
#else // ! NORMALHISTSCALE
                    if ( (Radius(p)+Radius(q)) < gd->deltaR*THETA)
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
#endif // ! NORMALHISTSCALE

                } else { // ! no-one-ball
                    for (l = More(q); l != Next(q); l = Next(l))
                        normal_walktree_sincos_N(cmd, gd, btable, cat2,
                                                 p,l,qsize/2,
                                                 hist, histN,
                                                 nbList, intList);
                } // ! no-one-ball
            }
        } else { // ! Type(q) == CELL
            sumnode_sincos_N(cmd, gd, btable, cat2,
                             p, ((cellptr)q), ((cellptr)q+1),
                             hist, histN, nbList, intList);
        } // ! Type(q) == CELL
    } // ! p != q
}

local void sumnode_sincos_N(struct  cmdline_data* cmd,
                            struct  global_data* gd,
                            bodyptr *btable, int cat2,
                            bodyptr p, cellptr start, cellptr finish,
                            gdhistptr_sincos_omp_ggg hist,
                            gdhistptr_sincos_omp_ggg_N histN,
                            int *nbList, int *intList)
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
    real xiN;
    REAL cosphi,sinphi;
    int iq;
#ifdef THREEPCFSHEAR
    real gamma1, gamma2;
#endif

    q = start;
    if (scanopt(cmd->options, "read-mask"))
        if (Mask(q)==FALSE) return;

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
        } // ! smooth-pivot
        //E


#ifndef NORMALHISTSCALE
//B useLogHist section
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

#ifdef TWOPCF
                hist->histNthread[n] = hist->histNthread[n] + 1.;
                hist->histWthread[n] = hist->histWthread[n] + 1.;
                hist->histNNSubXi2pcfthread[n] =
                                            hist->histNNSubXi2pcfthread[n] + 1.;
                //B kappa Avg Rmin
                hist->histNNSubXi2pcfthreadp[n] =
                                            hist->histNNSubXi2pcfthreadp[n] + 1.;
                //E
#endif

                hist->histNNSubthread[n] = hist->histNNSubthread[n] + 1.;
                histN->histNNSubthread[n] = histN->histNNSubthread[n] + 1.;

                xi = Weight(q)*Kappa(q);
                xiN = Weight(q);

#ifdef THREEPCFSHEAR
                gamma1 = Gamma1(q);
                gamma2 = Gamma2(q);
#endif
#if defined(THREEPCFCONVERGENCE) || defined(THREEPCFSHEAR)
                POLARAXIS_MAIN;
#endif
#ifdef THREEPCFCONVERGENCE
                if (cmd->mChebyshev<7) {
                    CHEBYSHEVTUOMPNANY;
                    CHEBYSHEVTUOMPSINCOSANY
                } else {
                    CHEBYSHEVTUOMPN;
                    CHEBYSHEVTUOMP;
                }
#endif

#ifdef THREEPCFSHEAR
                CHEBYSHEVTUOMPGGGANY;
#endif
#ifdef TWOPCF
                hist->histXi2pcfthreadsub[n] += xi;
#endif

                hist->nbbcalcthread += 1;
            } // ! n <= sizeHistN && n >= 1
        } // ! dr1 > rminHist

//E useLogHist section
#else // ! NORMALHISTSCALE

        if(dr1>cmd->rminHist) {
#ifdef ADDPIVOTNEIGHBOURS
            iq = (bodyptr)q-btable[cat2];
            activeInt[*intList]=iq;
            *intList +=1;
            if (*intList > actlenInt)
                error("intList: too many neighbors\n");
#endif
            n = (int) ( (dr1-cmd->rminHist) * gd->i_deltaR) + 1;
            if (n<=cmd->sizeHistN && n>=1) {

#ifdef TWOPCF
                hist->histNthread[n] = hist->histNthread[n] + 1.;
                hist->histWthread[n] = hist->histWthread[n] + 1.;
                hist->histNNSubXi2pcfthread[n] =
                                            hist->histNNSubXi2pcfthread[n] + 1.;
                //B kappa Avg Rmin
                hist->histNNSubXi2pcfthreadp[n] =
                                            hist->histNNSubXi2pcfthreadp[n] + 1.;
                //E
#endif

                hist->histNNSubthread[n] = hist->histNNSubthread[n] + 1.;
                histN->histNNSubthread[n] = histN->histNNSubthread[n] + 1.;

                xi = Weight(q)*Kappa(q);
                xiN = Weight(q);

#ifdef THREEPCFSHEAR
                gamma1 = Gamma1(q);
                gamma2 = Gamma2(q);
#endif
#if defined(THREEPCFCONVERGENCE) || defined(THREEPCFSHEAR)
                POLARAXIS_MAIN;
#endif
#ifdef THREEPCFCONVERGENCE
                if (cmd->mChebyshev<7) {
                    CHEBYSHEVTUOMPNANY;
                    CHEBYSHEVTUOMPSINCOSANY
                } else {
                    CHEBYSHEVTUOMPN;
                    CHEBYSHEVTUOMP;
                }
#endif

#ifdef THREEPCFSHEAR
                CHEBYSHEVTUOMPGGGANY;
#endif
#ifdef TWOPCF
                hist->histXi2pcfthreadsub[n] += xi;
#endif

                hist->nbbcalcthread += 1;
            } // ! n <= sizeHistN && n >= 1
        } // ! dr1 > rminHist

#endif // ! NORMALHISTSCALE

    } // ! accept_body
}

local void sumnode_sincos_cell_N(struct  cmdline_data* cmd,
                                 struct  global_data* gd,
                                 bodyptr *btable, int cat2,
                                 bodyptr p, cellptr start, cellptr finish,
                                 gdhistptr_sincos_omp_ggg hist,
                                 gdhistptr_sincos_omp_ggg_N histN,
                                 int *nbList, int *intList)
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
    real xiN;
    REAL cosphi,sinphi;
#ifdef THREEPCFSHEAR
    real gamma1, gamma2;
#endif

    q = start;
    if (scanopt(cmd->options, "read-mask"))
        if (Mask(q)==FALSE) return;

    if (accept_body(cmd, gd, p, (nodeptr)q, &dr1, dr)) {

#ifndef NORMALHISTSCALE
//B useLogHist section
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

#ifdef TWOPCF
                hist->histNthread[n] = hist->histNthread[n] +  Nb(q);
                hist->histWthread[n] = hist->histWthread[n] +  Nb(q);
                hist->histNNSubXi2pcfthread[n] =
                                            hist->histNNSubXi2pcfthread[n] + 1.0;
                //B kappa Avg Rmin
                hist->histNNSubXi2pcfthreadp[n] =
                                        hist->histNNSubXi2pcfthreadp[n] + 1.0;
                //E
#endif

                //B 3pcf convergence
                hist->histNNSubthread[n] = hist->histNNSubthread[n] + 1.0;
                histN->histNNSubthread[n] = histN->histNNSubthread[n] + 1.0;
                //E

#ifdef NONORMHIST
                if (scanopt(cmd->options, "no-normalize-HistZeta")) {
//B begin corrections
                    xi = Nb(q)*Kappa(q);            // Kappa is cell average
                    xiN = Weight(q);                // Weight sum in cell
                } else {
#ifdef KappaAvgON
                    xi = KappaAvg(q)/Nb(q);         // Weight(q)?
#else
//                    xi = Weight(q)*Kappa(q);
                    xi = Kappa(q);                  // original line
#endif
                    xiN = 1.0;                      // original line
//                    xiN = Weight(q);              // average of weights
                }
#else // ! NONORMHIST
#ifdef KappaAvgON
                xi = KappaAvg(q)/Nb(q);             // Weight(q)
#else
                xi = Kappa(q);                      // original line
//                xi = Weight(q)*Kappa(q);
#endif
                xiN = 1.0;                          // original line
//                xiN = Weight(q);                  // average of weights
#endif // ! NONORMHIST
//E

#ifdef THREEPCFSHEAR
                gamma1 = Gamma1(q);
                gamma2 = Gamma2(q);
#endif
#if defined(THREEPCFCONVERGENCE) || defined(THREEPCFSHEAR)
                POLARAXIS_MAIN;
#endif
#ifdef THREEPCFCONVERGENCE
                if (cmd->mChebyshev<7) {
                    CHEBYSHEVTUOMPNANY;
                    CHEBYSHEVTUOMPSINCOSANY
                } else {
                    CHEBYSHEVTUOMPN;
                    CHEBYSHEVTUOMP;
                }
#endif

#ifdef THREEPCFSHEAR
                CHEBYSHEVTUOMPGGGANY;
#endif
#ifdef TWOPCF
#ifdef NONORMHIST
                if (scanopt(cmd->options, "no-normalize-HistZeta")) {
                    hist->histXi2pcfthreadsub[n] += xi/Nb(q);
                } else {
                    hist->histXi2pcfthreadsub[n] += xi;
                }
#else
                hist->histXi2pcfthreadsub[n] += xi;
#endif
#endif

                hist->nbccalcthread += 1;
            } // ! n <= sizeHistN && n >= 1
        } // ! dr1 > rminHist
//E useLogHist section

#else // ! NORMALHISTSCALE

        if(dr1>cmd->rminHist) {
#ifdef ADDPIVOTNEIGHBOURS
            INTEGER iq;
            iq = (bodyptr)q-btable[cat2];
            activeInt[*intList]=iq;
            *intList +=1;
            if (*intList > actlenInt)
                error("intList: too many neighbors\n");
#endif
            n = (int) ( (dr1-cmd->rminHist) * gd->i_deltaR) + 1;
            if (n<=cmd->sizeHistN && n>=1) {

#ifdef TWOPCF
                hist->histNthread[n] = hist->histNthread[n] +  Nb(q);
                hist->histWthread[n] = hist->histWthread[n] +  Nb(q);
                hist->histNNSubXi2pcfthread[n] =
                                            hist->histNNSubXi2pcfthread[n] + 1.0;
                //B kappa Avg Rmin
                hist->histNNSubXi2pcfthreadp[n] =
                                        hist->histNNSubXi2pcfthreadp[n] + 1.0;
                //E
#endif

                //B 3pcf convergence
                hist->histNNSubthread[n] = hist->histNNSubthread[n] + 1.0;
                histN->histNNSubthread[n] = histN->histNNSubthread[n] + 1.0;
                //E

#ifdef NONORMHIST
                if (scanopt(cmd->options, "no-normalize-HistZeta")) {
//B begin corrections
                    xi = Nb(q)*Kappa(q);            // Kappa is cell average
                    xiN = Weight(q);                // Weight sum in cell
                } else {
#ifdef KappaAvgON
                    xi = KappaAvg(q)/Nb(q);         // Weight(q)?
#else
//                    xi = Weight(q)*Kappa(q);
                    xi = Kappa(q);                  // original line
#endif
                    xiN = 1.0;                      // original line
//                    xiN = Weight(q);              // average of weights
                }
#else // ! NONORMHIST
#ifdef KappaAvgON
                xi = KappaAvg(q)/Nb(q);             // Weight(q)
#else
                xi = Kappa(q);                      // original line
//                xi = Weight(q)*Kappa(q);
#endif
                xiN = 1.0;                          // original line
//                xiN = Weight(q);                  // average of weights
#endif // ! NONORMHIST
//E

#ifdef THREEPCFSHEAR
                gamma1 = Gamma1(q);
                gamma2 = Gamma2(q);
#endif
#if defined(THREEPCFCONVERGENCE) || defined(THREEPCFSHEAR)
                POLARAXIS_MAIN;
#endif
#ifdef THREEPCFCONVERGENCE
                if (cmd->mChebyshev<7) {
                    CHEBYSHEVTUOMPNANY;
                    CHEBYSHEVTUOMPSINCOSANY
                } else {
                    CHEBYSHEVTUOMPN;
                    CHEBYSHEVTUOMP;
                }
#endif

#ifdef THREEPCFSHEAR
                CHEBYSHEVTUOMPGGGANY;
#endif
#ifdef TWOPCF
#ifdef NONORMHIST
                if (scanopt(cmd->options, "no-normalize-HistZeta")) {
                    hist->histXi2pcfthreadsub[n] += xi/Nb(q);
                } else {
                    hist->histXi2pcfthreadsub[n] += xi;
                }
#else
                hist->histXi2pcfthreadsub[n] += xi;
#endif
#endif

                hist->nbccalcthread += 1;
            } // ! n <= sizeHistN && n >= 1
        } // ! dr1 > rminHist

#endif // ! NORMALHISTSCALE

    } // ! accept_body
}

local int computeBodyProperties_sincos_ggg_N(struct  cmdline_data* cmd,
                                            struct  global_data* gd,
                                            bodyptr p, int nbody,
                                            gdhistptr_sincos_omp_ggg_N hist)
{
    int n;
    int m;
    real xi;
    
//B add smooth-pivot corrections
#ifdef NONORMHIST
    xi = Weight(p);
#else
    xi = Weight(p)/nbody;
#endif
//E

#ifndef NONORMHIST
    for (m=1; m<=cmd->mChebyshev+1; m++) {
        for (n=1; n<=cmd->sizeHistN; n++) {
            hist->histXithreadcos[m][n] /= MAX(hist->histNNSubthread[n],1.0);
            hist->histXithreadsin[m][n] /= MAX(hist->histNNSubthread[n],1.0);
        }
    }
#endif


#ifdef THREEPCFCONVERGENCE
    for (m=1; m<=cmd->mChebyshev+1; m++) {
        OUTVP_ext(hist->xiOUTVPcos,
            hist->histXithreadcos[m], hist->histXithreadcos[m], cmd->sizeHistN);
        OUTVP_ext(hist->xiOUTVPsin,
            hist->histXithreadsin[m], hist->histXithreadsin[m],cmd->sizeHistN);
        OUTVP_ext(hist->xiOUTVPsincos,
            hist->histXithreadsin[m], hist->histXithreadcos[m],cmd->sizeHistN);
        OUTVP_ext(hist->xiOUTVPcossin,
            hist->histXithreadcos[m], hist->histXithreadsin[m],cmd->sizeHistN);
#ifdef PIVOTLOOP
        CLRM_ext(hist->histZetaMtmpcos[m], cmd->sizeHistN);
        CLRM_ext(hist->histZetaMtmpsin[m], cmd->sizeHistN);
        CLRM_ext(hist->histZetaMtmpsincos[m], cmd->sizeHistN);
        CLRM_ext(hist->histZetaMtmpcossin[m], cmd->sizeHistN);
        MULMS_ext(hist->histZetaMtmpcos[m],hist->xiOUTVPcos,xi,cmd->sizeHistN);
        MULMS_ext(hist->histZetaMtmpsin[m],hist->xiOUTVPsin,xi,cmd->sizeHistN);
        MULMS_ext(hist->histZetaMtmpsincos[m],
                  hist->xiOUTVPsincos,xi,cmd->sizeHistN);
        MULMS_ext(hist->histZetaMtmpcossin[m],
                  hist->xiOUTVPcossin,xi,cmd->sizeHistN);
        ADDM_ext(hist->histZetaMthreadcos[m],
            hist->histZetaMthreadcos[m],hist->histZetaMtmpcos[m],cmd->sizeHistN);
        ADDM_ext(hist->histZetaMthreadsin[m],
            hist->histZetaMthreadsin[m],hist->histZetaMtmpsin[m],cmd->sizeHistN);
        ADDM_ext(hist->histZetaMthreadsincos[m],
            hist->histZetaMthreadsincos[m],
            hist->histZetaMtmpsincos[m],cmd->sizeHistN);
        ADDM_ext(hist->histZetaMthreadcossin[m],
            hist->histZetaMthreadcossin[m],
            hist->histZetaMtmpcossin[m],cmd->sizeHistN);
#else
        CLRM_ext(hist->histZetaMtmpcos, cmd->sizeHistN);
        CLRM_ext(hist->histZetaMtmpsin, cmd->sizeHistN);
        CLRM_ext(hist->histZetaMtmpsincos, cmd->sizeHistN);
        CLRM_ext(hist->histZetaMtmpcossin, cmd->sizeHistN);
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
#endif
    }
#endif

    return SUCCESS;
}

#endif // ! NMultipoles


//B Routines as in cballsutils

local int search_init_gd_sincos_omp_ggg(struct  cmdline_data* cmd,
                                        struct  global_data* gd,
                                        gdlptr_sincos_omp_ggg gdl)
{
    int n;
    int m;
    INTEGER bytes_tot_local=0;

#ifdef TWOPCF
    gdl->histNN = dvector(1,cmd->sizeHistN);
    gdl->histWW = dvector(1,cmd->sizeHistN);
    gdl->histCF = dvector(1,cmd->sizeHistN);
    gdl->histNNSubXi2pcf = dvector(1,cmd->sizeHistN);
    //B kappa Avg Rmin
    gdl->histNNSubXi2pcftotal = dvector(1,cmd->sizeHistN);
    //E
    gdl->histXi2pcf = dvector(1,cmd->sizeHistN);
    bytes_tot_local += 5*cmd->sizeHistN*sizeof(real);
#endif

    gdl->histNNSub = dvector(1,cmd->sizeHistN);
    bytes_tot_local += 1*cmd->sizeHistN*sizeof(real);

#ifdef THREEPCFCONVERGENCE
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
    gdl->histXi3pcf = dmatrix3D(1,cmd->sizeHistPhi,
                                1,cmd->sizeHistN,1,cmd->sizeHistN);
    bytes_tot_local +=
                (cmd->sizeHistN*cmd->sizeHistN*cmd->sizeHistPhi)*sizeof(real);
#endif

    gd->bytes_tot += bytes_tot_local;
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
    "\nsearch_init_gd_octree_ggg: Allocated %g MByte for histograms storage.\n",
    bytes_tot_local*INMB);


#ifdef TWOPCF
    for (n = 1; n <= cmd->sizeHistN; n++) {
        gdl->histNN[n] = 0.0;
        gdl->histWW[n] = 0.0;
        gdl->histCF[n] = 0.0;
        gdl->histNNSubXi2pcf[n] = 0.0;
//B kappa Avg Rmin
        gdl->histNNSubXi2pcftotal[n] = 0.0;
//E
        gdl->histXi2pcf[n] = 0.0;
    }
#endif

    for (n = 1; n <= cmd->sizeHistN; n++)
        gdl->histNNSub[n] = 0.0;

#ifdef THREEPCFCONVERGENCE
    for (m = 1; m <= cmd->mChebyshev+1; m++) {
        CLRM_ext(gdl->histZetaMcos[m], cmd->sizeHistN);
        CLRM_ext(gdl->histZetaMsin[m], cmd->sizeHistN);
        CLRM_ext(gdl->histZetaMsincos[m], cmd->sizeHistN);
        CLRM_ext(gdl->histZetaMcossin[m], cmd->sizeHistN);
//#ifdef PIVOTLOOP
        CLRM_ext(gdl->histZetaM[m], cmd->sizeHistN);
//#endif
    }
#endif

    gd->nbbcalc = gd->nbccalc = gd->ncccalc = 0;

    return SUCCESS;
}

local int search_free_gd_sincos_omp_ggg(struct  cmdline_data* cmd,
                                         struct  global_data* gd,
                                        gdlptr_sincos_omp_ggg gdl)
{
#ifdef THREEPCFCONVERGENCE
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
#endif

    free_dvector(gdl->histNNSub,1,cmd->sizeHistN);

#ifdef TWOPCF
    free_dvector(gdl->histXi2pcf,1,cmd->sizeHistN);
    //B kappa Avg Rmin
    free_dvector(gdl->histNNSubXi2pcftotal,1,cmd->sizeHistN);
    //E
    free_dvector(gdl->histNNSubXi2pcf,1,cmd->sizeHistN);
    //
    free_dvector(gdl->histCF,1,cmd->sizeHistN);
    free_dvector(gdl->histNN,1,cmd->sizeHistN);
#endif

    return SUCCESS;
}

local int search_init_sincos_omp_ggg(struct  cmdline_data* cmd,
                                     struct  global_data* gd,
                                     gdhistptr_sincos_omp_ggg hist)
{
    int n;
    int m;
    INTEGER bytes_tot_local=0;

#ifdef TWOPCF
    hist->histNthread = dvector(1,cmd->sizeHistN);
    hist->histWthread = dvector(1,cmd->sizeHistN);
    hist->histWWthread = dvector(1,cmd->sizeHistN);
    hist->histNNSubXi2pcfthread = dvector(1,cmd->sizeHistN);
    //B kappa Avg Rmin
    hist->histNNSubXi2pcfthreadp = dvector(1,cmd->sizeHistN);
    hist->histNNSubXi2pcfthreadtotal = dvector(1,cmd->sizeHistN);
    //E
    hist->histXi2pcfthread = dvector(1,cmd->sizeHistN);
    hist->histXi2pcfthreadsub = dvector(1,cmd->sizeHistN);
    bytes_tot_local += 6*cmd->sizeHistN*sizeof(real);
#endif

    //B 3pcf convergence & shear
    hist->ChebsT = dvector(1,cmd->mChebyshev+1);
    hist->ChebsU = dvector(1,cmd->mChebyshev+1);
    //E 3pcf convergence & shear
    bytes_tot_local += 2.0*(cmd->mChebyshev+1)*sizeof(real);

    hist->histNNSubthread = dvector(1,cmd->sizeHistN);
    bytes_tot_local += cmd->sizeHistN*sizeof(real);

#ifdef THREEPCFCONVERGENCE
    hist->histXithreadcos = dmatrix(1,cmd->mChebyshev+1,1,cmd->sizeHistN);
    hist->histXithreadsin = dmatrix(1,cmd->mChebyshev+1,1,cmd->sizeHistN);
    bytes_tot_local += 2.0*(cmd->mChebyshev+1)*cmd->sizeHistN*sizeof(real);
    hist->histZetaMthreadcos =
            dmatrix3D(1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    hist->histZetaMthreadsin =
            dmatrix3D(1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    hist->histZetaMthreadsincos =
            dmatrix3D(1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    hist->histZetaMthreadcossin =
            dmatrix3D(1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    bytes_tot_local +=
            4.0*(cmd->mChebyshev+1)*cmd->sizeHistN*cmd->sizeHistN*sizeof(real);
    hist->xiOUTVPcos = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
    hist->xiOUTVPsin = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
    hist->xiOUTVPsincos = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
    hist->xiOUTVPcossin = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
    bytes_tot_local += 4.0*cmd->sizeHistN*cmd->sizeHistN*sizeof(real);
#ifdef PIVOTLOOP
    hist->histZetaMtmpcos = dmatrix3D(1,cmd->mChebyshev+1,
                                      1,cmd->sizeHistN,1,cmd->sizeHistN);
    hist->histZetaMtmpsin = dmatrix3D(1,cmd->mChebyshev+1,
                                      1,cmd->sizeHistN,1,cmd->sizeHistN);
    hist->histZetaMtmpsincos = dmatrix3D(1,cmd->mChebyshev+1,
                                         1,cmd->sizeHistN,1,cmd->sizeHistN);
    hist->histZetaMtmpcossin = dmatrix3D(1,cmd->mChebyshev+1,
                                         1,cmd->sizeHistN,1,cmd->sizeHistN);
    hist->histZetaMtmp = dmatrix3D(1,cmd->mChebyshev+1,
                                      1,cmd->sizeHistN,1,cmd->sizeHistN);
    bytes_tot_local +=
            5.0*(cmd->mChebyshev+1)*cmd->sizeHistN*cmd->sizeHistN*sizeof(real);
#else
    hist->histZetaMtmpcos = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
    hist->histZetaMtmpsin = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
    hist->histZetaMtmpsincos = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
    hist->histZetaMtmpcossin = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
    bytes_tot_local += 4.0*cmd->sizeHistN*cmd->sizeHistN*sizeof(real);
#endif
#endif

#ifdef THREEPCFSHEAR
    hist->histg1threadcos = dmatrix(1,cmd->mChebyshev+1,1,cmd->sizeHistN);
    hist->histg1threadsin = dmatrix(1,cmd->mChebyshev+1,1,cmd->sizeHistN);
    hist->histg2threadcos = dmatrix(1,cmd->mChebyshev+1,1,cmd->sizeHistN);
    hist->histg2threadsin = dmatrix(1,cmd->mChebyshev+1,1,cmd->sizeHistN);
    hist->histReGthread = dmatrix(1,cmd->mChebyshev+1,1,cmd->sizeHistN);
    hist->histImGthread = dmatrix(1,cmd->mChebyshev+1,1,cmd->sizeHistN);
    bytes_tot_local += 6.0*(cmd->mChebyshev+1)*cmd->sizeHistN*sizeof(real);
    hist->histReGNthread = dmatrix(1,cmd->mChebyshev,1,cmd->sizeHistN);
    hist->histImGNthread = dmatrix(1,cmd->mChebyshev,1,cmd->sizeHistN);
    bytes_tot_local += 2.0*(cmd->mChebyshev)*cmd->sizeHistN*sizeof(real);
#endif

#ifdef TWOPCF
    for (n = 1; n <= cmd->sizeHistN; n++) {
        hist->histNthread[n] = 0.0;
        hist->histWthread[n] = 0.0;
        hist->histWWthread[n] = 0.0;
        hist->histNNSubXi2pcfthread[n] = 0.0;
//B kappa Avg Rmin
        hist->histNNSubXi2pcfthreadp[n] = 0.0;
        hist->histNNSubXi2pcfthreadtotal[n] = 0.0;
//E
        hist->histXi2pcfthread[n] = 0.0;
        hist->histXi2pcfthreadsub[n] = 0.0;
    }
#endif

    for (n = 1; n <= cmd->sizeHistN; n++)
        hist->histNNSubthread[n] = 0.0;

#ifdef THREEPCFCONVERGENCE
    for (m = 1; m <= cmd->mChebyshev+1; m++) {
        CLRM_ext(hist->histZetaMthreadcos[m], cmd->sizeHistN);
        CLRM_ext(hist->histZetaMthreadsin[m], cmd->sizeHistN);
        CLRM_ext(hist->histZetaMthreadsincos[m], cmd->sizeHistN);
        CLRM_ext(hist->histZetaMthreadcossin[m], cmd->sizeHistN);
#ifdef PIVOTLOOP
        CLRM_ext(hist->histZetaMtmp[m], cmd->sizeHistN);
#endif
    }
#endif

#ifdef THREEPCFSHEAR
    for (m = 1; m <= cmd->mChebyshev+1; m++) {
        CLRV_ext(hist->histg1threadcos[m], cmd->sizeHistN);
        CLRV_ext(hist->histg1threadsin[m], cmd->sizeHistN);
        CLRV_ext(hist->histg2threadcos[m], cmd->sizeHistN);
        CLRV_ext(hist->histg2threadsin[m], cmd->sizeHistN);
        CLRV_ext(hist->histReGthread[m], cmd->sizeHistN);
        CLRV_ext(hist->histImGthread[m], cmd->sizeHistN);
    }
    for (m = 1; m <= cmd->mChebyshev; m++) {
        CLRV_ext(hist->histReGNthread[m], cmd->sizeHistN);
        CLRV_ext(hist->histImGNthread[m], cmd->sizeHistN);
    }
#endif

    hist->nbbcalcthread = 0;
    hist->nbccalcthread = 0;

    gd->bytes_tot += (cmd->numthreads)*bytes_tot_local;

    return SUCCESS;
}

local int search_free_sincos_omp_ggg(struct  cmdline_data* cmd,
                                         struct  global_data* gd,
                                      gdhistptr_sincos_omp_ggg hist)
{
#ifdef THREEPCFSHEAR
    free_dmatrix(hist->histImGNthread,1,cmd->mChebyshev,1,cmd->sizeHistN);
    free_dmatrix(hist->histReGNthread,1,cmd->mChebyshev,1,cmd->sizeHistN);
    free_dmatrix(hist->histImGthread,1,cmd->mChebyshev+1,1,cmd->sizeHistN);
    free_dmatrix(hist->histReGthread,1,cmd->mChebyshev+1,1,cmd->sizeHistN);
    free_dmatrix(hist->histg2threadsin,1,cmd->mChebyshev+1,1,cmd->sizeHistN);
    free_dmatrix(hist->histg2threadcos,1,cmd->mChebyshev+1,1,cmd->sizeHistN);
    free_dmatrix(hist->histg1threadsin,1,cmd->mChebyshev+1,1,cmd->sizeHistN);
    free_dmatrix(hist->histg1threadcos,1,cmd->mChebyshev+1,1,cmd->sizeHistN);
#endif

#ifdef THREEPCFCONVERGENCE
#ifdef PIVOTLOOP
    free_dmatrix3D(hist->histZetaMtmpcossin,1,cmd->mChebyshev+1,
                   1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix3D(hist->histZetaMtmpsincos,1,cmd->mChebyshev+1,
                   1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix3D(hist->histZetaMtmpsin,1,cmd->mChebyshev+1,
                   1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix3D(hist->histZetaMtmpcos,1,cmd->mChebyshev+1,
                   1,cmd->sizeHistN,1,cmd->sizeHistN);
#else
    free_dmatrix(hist->histZetaMtmpcossin,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix(hist->histZetaMtmpsincos,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix(hist->histZetaMtmpsin,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix(hist->histZetaMtmpcos,1,cmd->sizeHistN,1,cmd->sizeHistN);
#endif
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
#endif

    free_dvector(hist->histNNSubthread,1,cmd->sizeHistN);

    //B 3pcf convergence & shear
    free_dvector(hist->ChebsU,1,cmd->mChebyshev+1);
    free_dvector(hist->ChebsT,1,cmd->mChebyshev+1);
    //E 3pcf convergence & shear

#ifdef TWOPCF
    free_dvector(hist->histXi2pcfthread,1,cmd->sizeHistN);
    free_dvector(hist->histXi2pcfthreadsub,1,cmd->sizeHistN);
    //B kappa Avg Rmin
    free_dvector(hist->histNNSubXi2pcfthreadtotal,1,cmd->sizeHistN);
    free_dvector(hist->histNNSubXi2pcfthreadp,1,cmd->sizeHistN);
    //E
    free_dvector(hist->histNNSubXi2pcfthread,1,cmd->sizeHistN);
    //
    free_dvector(hist->histWWthread,1,cmd->sizeHistN);
    free_dvector(hist->histWthread,1,cmd->sizeHistN);
    free_dvector(hist->histNthread,1,cmd->sizeHistN);
#endif

    return SUCCESS;
}


#ifdef NMultipoles
local int search_init_gd_sincos_omp_ggg_N(struct  cmdline_data* cmd,
                                        struct  global_data* gd,
                                          gdlptr_sincos_omp_ggg_N gdl)
{
    int n;
    int m;
    INTEGER bytes_tot_local=0;

#ifdef TWOPCF
    gdl->histNN = dvector(1,cmd->sizeHistN);
    gdl->histWW = dvector(1,cmd->sizeHistN);
    gdl->histCF = dvector(1,cmd->sizeHistN);
    gdl->histNNSubXi2pcf = dvector(1,cmd->sizeHistN);
    //B kappa Avg Rmin
    gdl->histNNSubXi2pcftotal = dvector(1,cmd->sizeHistN);
    //E
    gdl->histXi2pcf = dvector(1,cmd->sizeHistN);
    bytes_tot_local += 5*cmd->sizeHistN*sizeof(real);
#endif

    gdl->histNNSub = dvector(1,cmd->sizeHistN);
    bytes_tot_local += 1*cmd->sizeHistN*sizeof(real);

#ifdef THREEPCFCONVERGENCE
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
    gdl->histXi3pcf = dmatrix3D(1,cmd->sizeHistPhi,
                                1,cmd->sizeHistN,1,cmd->sizeHistN);
    bytes_tot_local +=
                (cmd->sizeHistN*cmd->sizeHistN*cmd->sizeHistPhi)*sizeof(real);
#endif

    gd->bytes_tot += bytes_tot_local;
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
    "search_init_gd_octree_ggg_N: Allocated %g MByte for histograms storage.\n",
    bytes_tot_local*INMB);


#ifdef TWOPCF
    for (n = 1; n <= cmd->sizeHistN; n++) {
        gdl->histNN[n] = 0.0;
        gdl->histWW[n] = 0.0;
        gdl->histCF[n] = 0.0;
        gdl->histNNSubXi2pcf[n] = 0.0;
//B kappa Avg Rmin
        gdl->histNNSubXi2pcftotal[n] = 0.0;
//E
        gdl->histXi2pcf[n] = 0.0;
    }
#endif

    for (n = 1; n <= cmd->sizeHistN; n++)
        gdl->histNNSub[n] = 0.0;

#ifdef THREEPCFCONVERGENCE
    for (m = 1; m <= cmd->mChebyshev+1; m++) {
        CLRM_ext(gdl->histZetaMcos[m], cmd->sizeHistN);
        CLRM_ext(gdl->histZetaMsin[m], cmd->sizeHistN);
        CLRM_ext(gdl->histZetaMsincos[m], cmd->sizeHistN);
        CLRM_ext(gdl->histZetaMcossin[m], cmd->sizeHistN);
//#ifdef PIVOTLOOP
        CLRM_ext(gdl->histZetaM[m], cmd->sizeHistN);
//#endif
    }
#endif

    gd->nbbcalc = gd->nbccalc = gd->ncccalc = 0;

    return SUCCESS;
}

local int search_free_gd_sincos_omp_ggg_N(struct  cmdline_data* cmd,
                                         struct  global_data* gd,
                                        gdlptr_sincos_omp_ggg_N gdl)
{
#ifdef THREEPCFCONVERGENCE
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
#endif

    free_dvector(gdl->histNNSub,1,cmd->sizeHistN);

#ifdef TWOPCF
    free_dvector(gdl->histXi2pcf,1,cmd->sizeHistN);
    //B kappa Avg Rmin
    free_dvector(gdl->histNNSubXi2pcftotal,1,cmd->sizeHistN);
    //E
    free_dvector(gdl->histNNSubXi2pcf,1,cmd->sizeHistN);
    //
    free_dvector(gdl->histCF,1,cmd->sizeHistN);
    free_dvector(gdl->histNN,1,cmd->sizeHistN);
#endif

    return SUCCESS;
}

local int search_init_sincos_omp_ggg_N(struct  cmdline_data* cmd,
                                     struct  global_data* gd,
                                     gdhistptr_sincos_omp_ggg_N hist)
{
    int n;
    int m;
    INTEGER bytes_tot_local=0;

#ifdef TWOPCF
    hist->histNthread = dvector(1,cmd->sizeHistN);
    hist->histWthread = dvector(1,cmd->sizeHistN);
    hist->histWWthread = dvector(1,cmd->sizeHistN);
    hist->histNNSubXi2pcfthread = dvector(1,cmd->sizeHistN);
    //B kappa Avg Rmin
    hist->histNNSubXi2pcfthreadp = dvector(1,cmd->sizeHistN);
    hist->histNNSubXi2pcfthreadtotal = dvector(1,cmd->sizeHistN);
    //E
    hist->histXi2pcfthread = dvector(1,cmd->sizeHistN);
    hist->histXi2pcfthreadsub = dvector(1,cmd->sizeHistN);
    bytes_tot_local += 6*cmd->sizeHistN*sizeof(real);
#endif

    //B 3pcf convergence & shear
    hist->ChebsT = dvector(1,cmd->mChebyshev+1);
    hist->ChebsU = dvector(1,cmd->mChebyshev+1);
    //E 3pcf convergence & shear
    bytes_tot_local += 2.0*(cmd->mChebyshev+1)*sizeof(real);

    hist->histNNSubthread = dvector(1,cmd->sizeHistN);
    bytes_tot_local += cmd->sizeHistN*sizeof(real);

#ifdef THREEPCFCONVERGENCE
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
#ifdef PIVOTLOOP
    hist->histZetaMtmpcos = dmatrix3D(1,cmd->mChebyshev+1,
                                      1,cmd->sizeHistN,1,cmd->sizeHistN);
    hist->histZetaMtmpsin = dmatrix3D(1,cmd->mChebyshev+1,
                                      1,cmd->sizeHistN,1,cmd->sizeHistN);
    hist->histZetaMtmpsincos = dmatrix3D(1,cmd->mChebyshev+1,
                                         1,cmd->sizeHistN,1,cmd->sizeHistN);
    hist->histZetaMtmpcossin = dmatrix3D(1,cmd->mChebyshev+1,
                                         1,cmd->sizeHistN,1,cmd->sizeHistN);
    hist->histZetaMtmp = dmatrix3D(1,cmd->mChebyshev+1,
                                      1,cmd->sizeHistN,1,cmd->sizeHistN);
#else
    hist->histZetaMtmpcos = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
    hist->histZetaMtmpsin = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
    hist->histZetaMtmpsincos = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
    hist->histZetaMtmpcossin = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
#endif
#endif

#ifdef THREEPCFSHEAR
    hist->histg1threadcos = dmatrix(1,cmd->mChebyshev+1,1,cmd->sizeHistN);
    hist->histg1threadsin = dmatrix(1,cmd->mChebyshev+1,1,cmd->sizeHistN);
    hist->histg2threadcos = dmatrix(1,cmd->mChebyshev+1,1,cmd->sizeHistN);
    hist->histg2threadsin = dmatrix(1,cmd->mChebyshev+1,1,cmd->sizeHistN);
    hist->histReGthread = dmatrix(1,cmd->mChebyshev+1,1,cmd->sizeHistN);
    hist->histImGthread = dmatrix(1,cmd->mChebyshev+1,1,cmd->sizeHistN);
    hist->histReGNthread = dmatrix(1,cmd->mChebyshev,1,cmd->sizeHistN);
    hist->histImGNthread = dmatrix(1,cmd->mChebyshev,1,cmd->sizeHistN);
#endif

#ifdef TWOPCF
    for (n = 1; n <= cmd->sizeHistN; n++) {
        hist->histNthread[n] = 0.0;
        hist->histWthread[n] = 0.0;
        hist->histWWthread[n] = 0.0;
        hist->histNNSubXi2pcfthread[n] = 0.0;
//B kappa Avg Rmin
        hist->histNNSubXi2pcfthreadp[n] = 0.0;
        hist->histNNSubXi2pcfthreadtotal[n] = 0.0;
//E
        hist->histXi2pcfthread[n] = 0.0;
        hist->histXi2pcfthreadsub[n] = 0.0;
    }
#endif

    for (n = 1; n <= cmd->sizeHistN; n++)
        hist->histNNSubthread[n] = 0.0;

#ifdef THREEPCFCONVERGENCE
    for (m = 1; m <= cmd->mChebyshev+1; m++) {
        CLRM_ext(hist->histZetaMthreadcos[m], cmd->sizeHistN);
        CLRM_ext(hist->histZetaMthreadsin[m], cmd->sizeHistN);
        CLRM_ext(hist->histZetaMthreadsincos[m], cmd->sizeHistN);
        CLRM_ext(hist->histZetaMthreadcossin[m], cmd->sizeHistN);
#ifdef PIVOTLOOP
        CLRM_ext(hist->histZetaMtmp[m], cmd->sizeHistN);
#endif
    }
#endif

#ifdef THREEPCFSHEAR
    for (m = 1; m <= cmd->mChebyshev+1; m++) {
        CLRV_ext(hist->histg1threadcos[m], cmd->sizeHistN);
        CLRV_ext(hist->histg1threadsin[m], cmd->sizeHistN);
        CLRV_ext(hist->histg2threadcos[m], cmd->sizeHistN);
        CLRV_ext(hist->histg2threadsin[m], cmd->sizeHistN);
        CLRV_ext(hist->histReGthread[m], cmd->sizeHistN);
        CLRV_ext(hist->histImGthread[m], cmd->sizeHistN);
    }
    for (m = 1; m <= cmd->mChebyshev; m++) {
        CLRV_ext(hist->histReGNthread[m], cmd->sizeHistN);
        CLRV_ext(hist->histImGNthread[m], cmd->sizeHistN);
    }
#endif

    hist->nbbcalcthread = 0;
    hist->nbccalcthread = 0;

    return SUCCESS;
}

local int search_free_sincos_omp_ggg_N(struct  cmdline_data* cmd,
                                         struct  global_data* gd,
                                      gdhistptr_sincos_omp_ggg_N hist)
{
#ifdef THREEPCFSHEAR
    free_dmatrix(hist->histImGNthread,1,cmd->mChebyshev,1,cmd->sizeHistN);
    free_dmatrix(hist->histReGNthread,1,cmd->mChebyshev,1,cmd->sizeHistN);
    free_dmatrix(hist->histImGthread,1,cmd->mChebyshev+1,1,cmd->sizeHistN);
    free_dmatrix(hist->histReGthread,1,cmd->mChebyshev+1,1,cmd->sizeHistN);
    free_dmatrix(hist->histg2threadsin,1,cmd->mChebyshev+1,1,cmd->sizeHistN);
    free_dmatrix(hist->histg2threadcos,1,cmd->mChebyshev+1,1,cmd->sizeHistN);
    free_dmatrix(hist->histg1threadsin,1,cmd->mChebyshev+1,1,cmd->sizeHistN);
    free_dmatrix(hist->histg1threadcos,1,cmd->mChebyshev+1,1,cmd->sizeHistN);
#endif

#ifdef THREEPCFCONVERGENCE
#ifdef PIVOTLOOP
    free_dmatrix3D(hist->histZetaMtmpcossin,1,cmd->mChebyshev+1,
                   1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix3D(hist->histZetaMtmpsincos,1,cmd->mChebyshev+1,
                   1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix3D(hist->histZetaMtmpsin,1,cmd->mChebyshev+1,
                   1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix3D(hist->histZetaMtmpcos,1,cmd->mChebyshev+1,
                   1,cmd->sizeHistN,1,cmd->sizeHistN);
#else
    free_dmatrix(hist->histZetaMtmpcossin,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix(hist->histZetaMtmpsincos,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix(hist->histZetaMtmpsin,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix(hist->histZetaMtmpcos,1,cmd->sizeHistN,1,cmd->sizeHistN);
#endif
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
#endif

    free_dvector(hist->histNNSubthread,1,cmd->sizeHistN);

    //B 3pcf convergence & shear
    free_dvector(hist->ChebsU,1,cmd->mChebyshev+1);
    free_dvector(hist->ChebsT,1,cmd->mChebyshev+1);
    //E 3pcf convergence & shear

#ifdef TWOPCF
    free_dvector(hist->histXi2pcfthread,1,cmd->sizeHistN);
    free_dvector(hist->histXi2pcfthreadsub,1,cmd->sizeHistN);
    //B kappa Avg Rmin
    free_dvector(hist->histNNSubXi2pcfthreadtotal,1,cmd->sizeHistN);
    free_dvector(hist->histNNSubXi2pcfthreadp,1,cmd->sizeHistN);
    //E
    free_dvector(hist->histNNSubXi2pcfthread,1,cmd->sizeHistN);
    //
    free_dvector(hist->histWWthread,1,cmd->sizeHistN);
    free_dvector(hist->histWthread,1,cmd->sizeHistN);
    free_dvector(hist->histNthread,1,cmd->sizeHistN);
#endif

    return SUCCESS;
}

#endif // ! NMultipoles


#ifdef TWOPCF
//B Computation of histogram of all B-B encounters
// The correlation function is estimated as:
//    xi=(V/v(r))*(DD(r)/N^2)
// where v(r)=4*pi*((r+dr/2)^3-(r-dr/2)^3)/3, V=box_size^3 and N is the
// total # particles.
local int search_compute_Xi_ggg(struct  cmdline_data* cmd,
                                struct  global_data* gd, int nbody,
                                gdlptr_sincos_omp_ggg gdl)
{
    int k;
    int n;
    real normFac;
    real Vol;

    Vol = 1.0;
    DO_COORD(k)
        Vol = Vol*gd->Box[k];

if (!cmd->useLogHist) {
    if ((scanopt(cmd->options, "cute-box"))) {
        gdl->histNN[1]-=nbody;
    }
}
    real *edd;
    real *corr;
    real *ercorr;
    edd = dvector(1,cmd->sizeHistN);
    corr = dvector(1,cmd->sizeHistN);
    ercorr = dvector(1,cmd->sizeHistN);
    real rho_av=(real)nbody/Vol;

    for (n = 1; n <= cmd->sizeHistN; n++)
        edd[n] = 1./rsqrt(gdl->histNN[n]);

    for (n = 1; n <= cmd->sizeHistN; n++) {
        if(gdl->histNN[n]==0) {
            corr[n]=0;
            ercorr[n]=0;
        } else {
            double r0,r1,vr,rho_r;
            if (cmd->useLogHist) {
                if (cmd->rminHist==0) {
                    r0 = rpow(10.0, ((real)(n-cmd->sizeHistN))/cmd->logHistBinsPD
                              + rlog10(cmd->rangeN) );
                    r1 = rpow(10.0, ((real)(n+1-cmd->sizeHistN))/cmd->logHistBinsPD
                              + rlog10(cmd->rangeN) );
                } else {
                    r0 = rpow(10.0, rlog10(cmd->rminHist) + ((real)(n))*gd->deltaR );
                    r1 = rpow(10.0, rlog10(cmd->rminHist) + ((real)(n+1))*gd->deltaR );
                }
            } else {
                r0=(real)n*gd->deltaR;
                r1=(real)(n+1)*gd->deltaR;
            }

#if (NDIM==3)
            if (scanopt(cmd->options, "cute-box")) {
                //B this version does not give same results as CB
                //      although the programming is the same...
                vr=4.0*PI*(r1*r1*r1-r0*r0*r0)/3.0;
                rho_r=gdl->histNN[n]/((real)nbody*vr);
                corr[n]=rho_r/rho_av-1;             // Correlation function
                ercorr[n]=(1+corr[n])*edd[n];       // Poisson errors
                gdl->histCF[n] = corr[n];            // Original line
                //E
            } else {
                normFac = Vol/(2.0*PI*rpow(gd->deltaR,3.0)*nbody*nbody);
// This line gives results for rdf (radial distribution function):
//                gd->histCF[n] = gd->histNN[n] * normFac / rsqr((int)n-0.5);
// This line gives results in agreement with CB:
                gdl->histCF[n] = gdl->histNN[n] * normFac / rsqr((int)n-0.5) -1.0;
            }
#else
            if (scanopt(cmd->options, "cute-box")) {
                // This should be CB version...
                normFac = Vol/(PI*rpow(gd->deltaR,2.0)*nbody*nbody);
                gdl->histCF[n] = gdl->histNN[n] * normFac / ((int)n-0.5) - 1.0;
            } else {
                normFac = Vol/(PI*rpow(gd->deltaR,2.0)*nbody*nbody);
// This line gives results for rdf (radial distribution function):
//                gd->histCF[n] = gd->histNN[n] * normFac / ((int)n-0.5);
// This line gives results in agreement with CB:
                gdl->histCF[n] = gdl->histNN[n] * normFac / ((int)n-0.5) - 1.0;
            }
#endif // ! NDIM
        }
    }

    free_dvector(ercorr,1,cmd->sizeHistN);
    free_dvector(corr,1,cmd->sizeHistN);
    free_dvector(edd,1,cmd->sizeHistN);

    return SUCCESS;
}


local int search_compute_HistN_ggg(struct  cmdline_data* cmd,
                                struct  global_data* gd, int nbody,
                                   gdlptr_sincos_omp_ggg gdl)
{
    int n;
    real normFac;

// Check this factor is correct
    normFac = 0.5;

    for (n = 1; n <= cmd->sizeHistN; n++)
        gdl->histNN[n] *= normFac;

    if (scanopt(cmd->options, "and-CF"))
        search_compute_Xi_ggg(cmd, gd, nbody, gdl);

    return SUCCESS;
}
#endif // ! TWOPCF


//E Routines as in cballsutils

local int print_info(struct cmdline_data* cmd,
                                  struct  global_data* gd)
{
    string routineName = "print_info";

    verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                           "searchcalc: Using octree-ggg-omp... \n");

    if (scanopt(cmd->options, "GGGCorrelation"))
        verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                            "computing using GGG routine... \n");


    if (cmd->usePeriodic==TRUE)
//        error("CheckParameters: can´t have periodic boundaries and OCTREEGGGOMP definition (usePeriodic=%d)\nSet usePeriodic=false\n",
//            cmd->usePeriodic);
        verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
    "%s: warning!! can´t have periodic boundaries and OCTREEGGGOMP definition (usePeriodic=%d)\nSet usePeriodic=false\n",
                            routineName,cmd->usePeriodic);

    if (cmd->useLogHist==FALSE)
//        error("CheckParameters: can´t have normal scale hist and OCTREEGGGOMP definition (useLogHist=%d)\nSet useLogHist=true\n",
//            cmd->useLogHist);
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
"%s: warning!! can´t have normal scale hist and OCTREEGGGOMP definition (useLogHist=%d)\nSet useLogHist=true\n",
                        routineName,cmd->useLogHist);

    if (cmd->computeTPCF==FALSE)
//        error("CheckParameters: can´t have computeTPCF=false and OCTREEGGGOMP definition (computeTPCF=%d)\nSet computeTPCF=true\n",
//            cmd->computeTPCF);
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
"%s: warning!! can´t have computeTPCF=false and OCTREEGGGOMP definition (computeTPCF=%d)\nSet computeTPCF=true\n",
                        routineName,cmd->computeTPCF);

#ifdef TWOPCF
    verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                        "with 2pcf computation... \n");
#endif
#ifdef THREEPCFCONVERGENCE
    if (cmd->computeTPCF==FALSE)
    error("\ncan´t have computeTPCF=false and THREEPCFCONVERGENCE definition\n",
          routineName);
    verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                        "with 3pcf convergence computation... \n");
#endif
#ifdef THREEPCFSHEAR
    verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                        "with 3pcf shear computation... \n");
#endif

#if NDIM == 2
//error("CheckParameters: OCTREEGGGOMP definition works only in a 3D unit sphere");
    verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
            "OCTREEGGGOMP definition working in a 2D box... \n");
#endif

#ifdef NMultipoles
    verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                        "with NMultipoles... \n");
#else
    verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                        "without NMultipoles... \n");
#endif
#ifdef NONORMHIST
    verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                        "with NONORMHIST... \n");
    if (scanopt(cmd->options, "no-normalize-HistZeta"))
        verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                            "with option no-normalize-HistZeta...\n");
#else
    verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                        "without NONORMHIST... \n");
#endif

#if defined(NMultipoles) && defined(NONORMHIST)
    if (scanopt(cmd->options, "no-normalize-HistZeta")) {
        if (scanopt(cmd->options, "edge-corrections"))
            verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                                "with option edge-corrections... \n");
    } else {
        if (scanopt(cmd->options, "edge-corrections")) {
            verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                                "option edge-corrections only works with %s... \n",
                                "no-normalize-HistZeta option added");
            // Check freeing allocated memory...
            error("going out...\n");
        }
    }
#else
    if (scanopt(cmd->options, "edge-corrections")) {
        verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                    "option edge-corrections only works with %s activated... \n",
                    "NMultipoles && NONORMHIST");
        // Check freeing allocated memory...
        error("going out...");
    }
#endif

#ifndef USEGSL
    if (scanopt(cmd->options, "edge-corrections"))
        verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
            "option edge-corrections is better computed with %s activated... \n",
            "USEGSL");
#endif

#ifdef POLARAXIS
    verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                        "with POLARAXIS... \n");
#endif
#ifdef ADDPIVOTNEIGHBOURS
    if (scanopt(cmd->options, "smooth-pivot"))
        verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                            "with ADDPIVOTNEIGHBOURS... \n");
#endif
    if (scanopt(cmd->options, "no-one-ball"))
        verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                            "with option no-one-ball... \n");
    if (scanopt(cmd->options, "smooth-pivot"))
        verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                            "with option smooth-pivot... rsmooth=%g\n",
                            gd->rsmooth[0]);
    if (scanopt(cmd->options, "default-rsmooth"))
        verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                            "with option default-rsmooth... \n");
    if (scanopt(cmd->options, "fix-rsmooth"))
        verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                            "with option fix-rsmooth... \n");

#ifdef PIVOTLOOP
    verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                        "with PIVOTLOOP... \n");
#endif

#ifdef BALLS4SCANLEV
    verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                        "with BALLS4SCANLEV... \n");
#endif

#ifdef NORMALHISTSCALE
    verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                        "with NORMALHISTSCALE... \n");
#endif

    if (scanopt(cmd->options, "read-mask"))
        verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                            "with option read-mask... \n");

    if (scanopt(cmd->options, "kappa-constant"))
        verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                            "with option kappa-constant... \n");
    if (scanopt(cmd->options, "celestial"))
        verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                            "with option celestial... \n");
    if (scanopt(cmd->options, "ra-reversed"))
        verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                            "with option ra-reversed... \n");

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

#define MHISTZETA \
"%16.8e %16.8e %16.8e %16.8e %16.8e %16.8e\n"

#define MHISTZETAHEADER \
"# [1] rBins; [2] diagonal; [3] theta2=Nbins/4.0; [4] theta2=2.0*Nbins/4.0; \
[5] theta2=3.0*Nbins/4.0; [6] theta2=4.0*Nbins/4.0 - 1.0\n"


#ifdef THREEPCFCONVERGENCE
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
                        struct  global_data* gd, gdlptr_sincos_omp_ggg gdl)
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
                                 gdlptr_sincos_omp_ggg gdl)
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
                                      gdlptr_sincos_omp_ggg gdl)
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
    free(data);
#else
    free_dvector(data,1,NP);
#endif

    free_dmatrix3D(histZetaG_Im,1,NP,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix3D(histZetaG,1,NP,1,cmd->sizeHistN,1,cmd->sizeHistN);

    return SUCCESS;
}

#endif // ! THREEPCFCONVERGENCE


#ifdef NMultipoles

#ifdef THREEPCFCONVERGENCE
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
                "_cossin_N", m, EXTFILES);
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
#endif // ! THREEPCFCONVERGENCE


#ifdef NONORMHIST

#ifdef THREEPCFCONVERGENCE
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
                    rbinlog = rlog10(cmd->rminHist)
                                + ((real)(n1)-0.5)*gd->deltaR;
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

#ifdef PIVOTLOOP
// Edge-correct each matrix ZetaM for each m multipole
local int HistZetaM_sincos_edge_effects(struct  cmdline_data* cmd,
                                        struct  global_data* gd,
                                        gdhistptr_sincos_omp_ggg hist,
                                        gdhistptr_sincos_omp_ggg_N histN)
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
                mat3[m][n1][n2] = hist->histZetaMtmpcos[m][n1][n2]
                                    + hist->histZetaMtmpsin[m][n1][n2];
                mat4[m][n1][n2] = histN->histZetaMtmpcos[m][n1][n2]
                                    + histN->histZetaMtmpsin[m][n1][n2];
            }
        }
    }

    for (n1=1; n1<=cmd->sizeHistN; n1++) {
        for (n2=1; n2<=cmd->sizeHistN; n2++) {
            matrixClm(cmd, gd, mat3, mat4, n1, n2, mat5);
            for (m=1; m<=cmd->mChebyshev+1; m++) {
                if (cmd->verbose_log==4)  {
                    verb_log_print(cmd->verbose_log, gd->outlog,
                                    "\n\nhistZetaM elements again (%d, %d):\n\n",
                                    n1, n2);
                    verb_log_print(cmd->verbose_log, gd->outlog,
                                    "%g\n",
                                    mat5[m][n1][n2]);
                }
                hist->histZetaMtmp[m][n1][n2] += mat5[m][n1][n2];
            }
        }
    }

    free_dmatrix3D(mat5,
                   1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix3D(mat4,
                   1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix3D(mat3,
                   1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);

    return SUCCESS;
}
#endif // ! PIVOTLOOP

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

#ifndef PIVOTLOOP
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
#else

    for (m=1; m<=cmd->mChebyshev+1; m++) {
        sprintf(namebuf, "%s%s_%d%s", gd->fpfnamehistZetaMFileName,
                "_EE", m, EXTFILES);
        verb_print_q(2, cmd->verbose, "Printing : to a file %s ...\n",namebuf);
        outstr = stropen(namebuf, "w!");
        for (n1=1; n1<=cmd->sizeHistN; n1++) {
            for (n2=1; n2<=cmd->sizeHistN; n2++) {
                fprintf(outstr,"%16.8e ",
                        mat5[m][n1][n2]
                        );
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

#endif

    return SUCCESS;
}
#endif // ! THREEPCFCONVERGENCE

#endif // ! NONORMHIST

#endif // ! NMultipoles

#undef MHISTZETAHEADER
#undef MHISTZETA


#ifdef TWOPCF
local int PrintHistNN(struct  cmdline_data* cmd, struct  global_data* gd,
                      gdlptr_sincos_omp_ggg gdl)
{
    real rBin, rbinlog;
    int n;
    stream outstr;

    outstr = stropen(gd->fpfnamehistNNFileName, "w!");

    verb_print_q(2, cmd->verbose,
               "Printing : to a file %s ...\n",gd->fpfnamehistNNFileName);

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
        fprintf(outstr,"%16.8e %16.8e\n",rBin,gdl->histNN[n]);
    }
    fclose(outstr);

    if (scanopt(cmd->options, "and-CF"))
        PrintHistCF(cmd, gd, gdl);

    return SUCCESS;
}

local int PrintHistCF(struct  cmdline_data* cmd, struct  global_data* gd,
                      gdlptr_sincos_omp_ggg gdl)
{
    real rBin, rbinlog;
    int n;
    stream outstr;

    outstr = stropen(gd->fpfnamehistCFFileName, "w!");

    verb_print_q(2, cmd->verbose,
               "Printing : to a file %s ...\n",gd->fpfnamehistCFFileName);

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
        fprintf(outstr,"%16.8e %16.8e\n",rBin,gdl->histCF[n]);
    }
    fclose(outstr);

    return SUCCESS;
}

// 180*60/Pi
#define RADTOARCMIN   3437.74677
local int PrintHistXi2pcf(struct  cmdline_data* cmd, struct  global_data* gd,
                          gdlptr_sincos_omp_ggg gdl)
{
    real rBin, rbinlog;
    int n;
    stream outstr;
    char namebuf[256];

    sprintf(namebuf, "%s%s%s", gd->fpfnamehistXi2pcfFileName,
            cmd->suffixOutFiles, EXTFILES);
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
        if (scanopt(cmd->options, "rbin-arcmin"))
            rBin = rBin*RADTOARCMIN;
        else
            if (scanopt(cmd->options, "rbin-degree"))
                rBin = rBin*RADTOARCMIN/60.0;
        fprintf(outstr,"%16.8e %16.8e\n",rBin,gdl->histXi2pcf[n]);
    }
    fclose(outstr);

    return SUCCESS;
}
#undef RADTOARCMIN

#endif // ! TWOPCF


//E Saving histograms section: case GGGCORRELATION:





