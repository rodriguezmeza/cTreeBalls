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
#include <float.h>
#ifdef OCTREEGGGMPI
#include "fcfc_octree_ggg_mpi.h"
#endif

#ifdef TPCFSHEAR
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

#ifdef TPCFSHEAR
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
#endif // ! TPCFSHEAR

//B Define structures:
typedef struct {
    bool threepcf_enabled;
#ifdef TWOPCF
    real *histNN;                                   // used
    real *histWW;                                   // used
    real *histCF;                                   // used
    real *histNNSubXi2pcf;                          // used
#ifdef SMOOTHPIVOT
    realptr histNNSubXi2pcftotal;                   // used
#endif
    real *histXi2pcf;                               // used
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
    // (EE) edge_effects
    real ***histZetaM_EE;
    real ***histZetaM_EE_Im;
#endif

#ifdef TPCFSHEAR
#endif

} gdl_sincos_omp_ggg, *gdlptr_sincos_omp_ggg;

typedef struct {
    bool threepcf_enabled;
#ifdef TWOPCF
    realptr histNthread;                            // used
    realptr histWthread;                            // used
    realptr histWWthread;                           // used
    realptr histNNSubXi2pcfthread;                  // used
#ifdef SMOOTHPIVOT
    realptr histNNSubXi2pcfthreadp;                 // used
    realptr histNNSubXi2pcfthreadtotal;             // used
#endif
    real *histXi2pcfthread;                         // used
    real *histXi2pcfthreadsub;                      // used
#endif

    real *ChebsT;
    real *ChebsU;

    realptr histNNSubthread;

#ifdef THREEPCFCONVERGENCE
    real **xiOUTVPcos;
    real **xiOUTVPsin;
    real **xiOUTVPsincos;
    real **xiOUTVPcossin;
    real **histZetaMtmpcos;
    real **histZetaMtmpsin;
    real **histZetaMtmpsincos;
    real **histZetaMtmpcossin;
    real ***histZetaMthreadcos;
    real ***histZetaMthreadsin;
    real ***histZetaMthreadsincos;
    real ***histZetaMthreadcossin;

    real **histXithreadcos;
    real **histXithreadsin;
    real **histXithreaddiagcos;
    real **histXithreaddiagsin;
    real **histXithreaddiagsincos;
#endif

#ifdef TPCFSHEAR
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

    compute_vector q0;
    real drpq2, drpq;
    compute_vector dr0;
    real cosb;
    real sinb;
} gdhist_sincos_omp_ggg, *gdhistptr_sincos_omp_ggg;

#ifdef NMultipoles
//B estructures for NMultipoles should be the same as KKK multipoles
//  check the memory allocation and freeing is updated...
typedef struct {
 #ifdef TWOPCF
     real *histNN;                                  // used
    real *histWW;                                   // used
     real *histCF;                                  // used
     real *histNNSubXi2pcf;                         // used
#ifdef SMOOTHPIVOT
     realptr histNNSubXi2pcftotal;                  // used
#endif
     real *histXi2pcf;                              // used
 #endif

     realptr histNNSub;

 #ifdef THREEPCFCONVERGENCE
     real ***histZetaMcos;
     real ***histZetaMsin;
     real ***histZetaMsincos;
     // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
     real ***histZetaMcossin;
     //E
 #endif

 #ifdef TPCFSHEAR
 #endif
} gdl_sincos_omp_ggg_N, *gdlptr_sincos_omp_ggg_N;

typedef struct {
#ifdef TWOPCF
    realptr histNthread;                            // used
    realptr histWthread;                            // used
    realptr histWWthread;                           // used
    realptr histNNSubXi2pcfthread;                  // used
#ifdef SMOOTHPIVOT
    realptr histNNSubXi2pcfthreadp;                 // used
    realptr histNNSubXi2pcfthreadtotal;             // used
#endif
    real *histXi2pcfthread;                         // used
    real *histXi2pcfthreadsub;                      // used
#endif

    real *ChebsT;
    real *ChebsU;

    realptr histNNSubthread;

#ifdef THREEPCFCONVERGENCE
    real **xiOUTVPcos;
    real **xiOUTVPsin;
    real **xiOUTVPsincos;
    real **xiOUTVPcossin;
    real **histZetaMtmpcos;
    real **histZetaMtmpsin;
    real **histZetaMtmpsincos;
    real **histZetaMtmpcossin;
    real ***histZetaMthreadcos;
    real ***histZetaMthreadsin;
    real ***histZetaMthreadsincos;
    real ***histZetaMthreadcossin;

    real **histXithreadcos;
    real **histXithreadsin;
    real **histXithreaddiagcos;
    real **histXithreaddiagsin;
    real **histXithreaddiagsincos;
#endif

#ifdef TPCFSHEAR
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

    compute_vector q0;
    real drpq2, drpq;
    compute_vector dr0;
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
                                  gdhistptr_sincos_omp_ggg);
local void sumnode_sincos(struct  cmdline_data* cmd,
                          struct  global_data* gd,
                          bodyptr *btable, int cat2,
                          bodyptr, cellptr, cellptr,
                          gdhistptr_sincos_omp_ggg);
local void sumnode_sincos_cell(struct  cmdline_data* cmd,
                               struct  global_data* gd,
                               bodyptr *btable, int cat2,
                               bodyptr p, cellptr start, cellptr finish,
                               gdhistptr_sincos_omp_ggg hist);

local int search_init_gd_sincos_omp_ggg(struct  cmdline_data* cmd,
                                        struct  global_data* gd,
                                        gdlptr_sincos_omp_ggg);
local int search_init_gd_sincos_omp_ggg_unguarded(struct cmdline_data *cmd,
                                        struct global_data *gd,
                                        gdlptr_sincos_omp_ggg);
local int search_free_gd_sincos_omp_ggg(struct  cmdline_data* cmd,
                                         struct  global_data* gd,
                                        gdlptr_sincos_omp_ggg);
local int search_init_sincos_omp_ggg(struct  cmdline_data* cmd,
                                  struct  global_data* gd,
                                     gdhistptr_sincos_omp_ggg hist);
local int search_init_sincos_omp_ggg_unguarded(struct cmdline_data *cmd,
                                  struct global_data *gd,
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
                                    gdhistptr_sincos_omp_ggg_N);
local void sumnode_sincos_N(struct  cmdline_data* cmd,
                            struct  global_data* gd,
                            bodyptr *btable, int cat2,
                            bodyptr, cellptr, cellptr,
                            gdhistptr_sincos_omp_ggg,
                            gdhistptr_sincos_omp_ggg_N);
local void sumnode_sincos_cell_N(struct  cmdline_data*,
                                 struct  global_data*,
                                 bodyptr *btable, int cat2,
                                 bodyptr, cellptr, cellptr,
                                 gdhistptr_sincos_omp_ggg,
                                 gdhistptr_sincos_omp_ggg_N);
local int search_init_gd_sincos_omp_ggg_N(struct  cmdline_data* cmd,
                                          struct  global_data* gd,
                                          gdlptr_sincos_omp_ggg_N);
local int search_init_gd_sincos_omp_ggg_N_unguarded(struct cmdline_data *cmd,
                                          struct global_data *gd,
                                          gdlptr_sincos_omp_ggg_N);
local int search_free_gd_sincos_omp_ggg_N(struct  cmdline_data* cmd,
                                         struct  global_data* gd,
                                          gdlptr_sincos_omp_ggg_N);
local int search_init_sincos_omp_ggg_N(struct  cmdline_data* cmd,
                                  struct  global_data* gd,
                                     gdhistptr_sincos_omp_ggg_N hist);
local int search_init_sincos_omp_ggg_N_unguarded(struct cmdline_data *cmd,
                                  struct global_data *gd,
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
    if (hist->drpq2==0) continue;
    //E
#else // ! POLARAXIS
    //B check NDIM and periodic...
    dRotation3D(Pos(p), ROTANGLE, ROTANGLE, ROTANGLE, hist->q0);
    DOTPSUBV(hist->drpq2, hist->dr0, Pos(p), hist->q0);
    hist->drpq = rsqrt(hist->drpq2);
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
    compute_vector vc;                              \
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
    compute_vector vc;                              \
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
    real s, sy; compute_vector pr0;                 \
    DOTVP(s, dr, hist->dr0);                        \
    cosphi = s/(dr1*hist->drpq);                    \
    CROSSVP(pr0,hist->dr0,Pos(p));                  \
    DOTVP(sy, dr, pr0);                             \
    if (rabs(cosphi)>1.0)                           \
        sinphi = 0.0;                               \
    else                                            \
        sinphi = rsqrt(1.0 - rsqr(cosphi));         \
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
local int PrintHistZetaM_sincos_edge_effects(struct  cmdline_data*,
                                             struct  global_data*,
                                             gdlptr_sincos_omp_ggg,
                                             gdlptr_sincos_omp_ggg_N);

#endif
#endif // ! NMultipoles

#ifdef NMultipoles
#define GGG_WINDOW_ORDERS(cmd) (2*(cmd)->mChebyshev + 1)
#endif

#if defined(NMultipoles) && defined(THREEPCFCONVERGENCE)
/* N_{ell-n} requires window modes through |ell-n| = 2*mChebyshev. */
#define GGG_ACCUMULATE_WINDOW_MODES                                      \
do {                                                                     \
    int _ggg_index;                                                       \
    int _ggg_orders = GGG_WINDOW_ORDERS(cmd);                            \
    real _ggg_cos_mode;                                                   \
    real _ggg_sin_mode;                                                   \
    real _ggg_xi_cos;                                                     \
    real _ggg_xi_sin;                                                     \
    histN->ChebsT[1] = 1.0;                                              \
    histN->ChebsU[1] = 0.0;                                              \
    for (_ggg_index = 1; _ggg_index <= _ggg_orders; _ggg_index++) {     \
        if (_ggg_index == 2) {                                           \
            histN->ChebsT[2] = cosphi;                                   \
            histN->ChebsU[2] = 1.0;                                      \
        } else if (_ggg_index > 2) {                                     \
            histN->ChebsT[_ggg_index] =                                  \
                2.0*cosphi*histN->ChebsT[_ggg_index-1]                   \
                - histN->ChebsT[_ggg_index-2];                           \
            histN->ChebsU[_ggg_index] =                                  \
                2.0*cosphi*histN->ChebsU[_ggg_index-1]                   \
                - histN->ChebsU[_ggg_index-2];                           \
        }                                                                \
        _ggg_cos_mode = histN->ChebsT[_ggg_index];                       \
        _ggg_sin_mode = histN->ChebsU[_ggg_index]*sinphi;                \
        _ggg_xi_cos = xiN*_ggg_cos_mode;                                 \
        _ggg_xi_sin = xiN*_ggg_sin_mode;                                 \
        histN->histXithreadcos[_ggg_index][n] += _ggg_xi_cos;            \
        histN->histXithreadsin[_ggg_index][n] += _ggg_xi_sin;            \
        histN->histXithreaddiagcos[_ggg_index][n] +=                     \
            _ggg_xi_cos*_ggg_xi_cos;                                     \
        histN->histXithreaddiagsin[_ggg_index][n] +=                     \
            _ggg_xi_sin*_ggg_xi_sin;                                     \
        histN->histXithreaddiagsincos[_ggg_index][n] +=                  \
            _ggg_xi_sin*_ggg_xi_cos;                                     \
    }                                                                    \
} while (0)

#define GGG_ACCUMULATE_NUMERATOR_DIAGONAL                               \
do {                                                                     \
    int _ggg_index;                                                       \
    for (_ggg_index = 1; _ggg_index <= cmd->mChebyshev + 1;            \
         _ggg_index++) {                                                 \
        real _ggg_xi_cos = xi*hist->ChebsT[_ggg_index];                 \
        real _ggg_xi_sin =                                               \
            xi*hist->ChebsU[_ggg_index]*sinphi;                          \
        hist->histXithreaddiagcos[_ggg_index][n] +=                     \
            _ggg_xi_cos*_ggg_xi_cos;                                     \
        hist->histXithreaddiagsin[_ggg_index][n] +=                     \
            _ggg_xi_sin*_ggg_xi_sin;                                     \
        hist->histXithreaddiagsincos[_ggg_index][n] +=                  \
            _ggg_xi_sin*_ggg_xi_cos;                                     \
    }                                                                    \
} while (0)

#endif

/* Radius(q) is an opening radius and is divided by theta in treeload.c.
 * Cell rejection needs the undivided geometric bound or it can drop valid
 * neighbors when theta > 1. */
local bool reject_cell_ggg(struct cmdline_data *cmd,
                           struct global_data *gd,
                           nodeptr p, nodeptr q)
{
    real distance_squared;
    real distance;
    real bounding_radius = Radius(q);
    compute_vector dr;

    DOTPSUBV(distance_squared, dr, Pos(p), Pos(q));
    if (cmd->usePeriodic) {
        VWrapAll(dr);
        DOTVP(distance_squared, dr, dr);
    }
    distance = rsqrt(distance_squared);
    if (cmd->theta > 0.0)
        bounding_radius *= cmd->theta;

    return distance >= gd->Rcut + bounding_radius;
}
//E

#if defined(NMultipoles) && defined(NONORMHIST)
local real normalize_zeta_ggg(real numerator,
                              gdlptr_sincos_omp_ggg_N gdlN,
                              int n1, int n2)
{
    real denominator = gdlN->histZetaMcos[1][n1][n2]
                     + gdlN->histZetaMsin[1][n1][n2];

    return cballs_normalize_or_zero(numerator, denominator);
}
#endif

#if defined(NMultipoles) && defined(NONORMHIST) \
    && defined(THREEPCFCONVERGENCE)
typedef struct {
    real re;
    real im;
} ggg_complex;

local ggg_complex ggg_complex_make(real re, real im)
{
    ggg_complex value = {re, im};
    return value;
}

local ggg_complex ggg_complex_sub(ggg_complex a, ggg_complex b)
{
    return ggg_complex_make(a.re-b.re, a.im-b.im);
}

local ggg_complex ggg_complex_mul(ggg_complex a, ggg_complex b)
{
    return ggg_complex_make(a.re*b.re-a.im*b.im,
                            a.re*b.im+a.im*b.re);
}

local real ggg_complex_abs2(ggg_complex value)
{
    return value.re*value.re + value.im*value.im;
}

local ggg_complex ggg_complex_div(ggg_complex numerator,
                                  ggg_complex denominator)
{
    real norm = ggg_complex_abs2(denominator);
    return ggg_complex_make(
        (numerator.re*denominator.re + numerator.im*denominator.im)/norm,
        (numerator.im*denominator.re - numerator.re*denominator.im)/norm);
}

local ggg_complex ggg_histogram_mode(real ***coscos, real ***sinsin,
                                     real ***sincos, real ***cossin,
                                     int order, int n1, int n2)
{
    int index = abs(order) + 1;
    ggg_complex value = ggg_complex_make(
        coscos[index][n1][n2] + sinsin[index][n1][n2],
        sincos[index][n1][n2] - cossin[index][n1][n2]);
    if (order < 0)
        value.im = -value.im;
    return value;
}

local int compute_edge_corrections_ggg(struct cmdline_data *cmd,
                                       struct global_data *gd,
                                       gdlptr_sincos_omp_ggg gdl,
                                       gdlptr_sincos_omp_ggg_N gdlN)
{
    int nmax = cmd->mChebyshev;
    int multipoles = 2*nmax + 1;
    size_t matrix_count;
    ggg_complex *matrix = NULL;
    ggg_complex *rhs = NULL;
    int n1, n2;
    int singular_bins = 0;

    if ((size_t)multipoles > SIZE_MAX/(size_t)multipoles) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "octree-ggg-omp: edge-correction workspace overflow");
        return FAILURE;
    }
    matrix_count = (size_t)multipoles*(size_t)multipoles;
    matrix = calloc(matrix_count, sizeof(*matrix));
    rhs = calloc((size_t)multipoles, sizeof(*rhs));
    if (matrix == NULL || rhs == NULL) {
        free(matrix);
        free(rhs);
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "octree-ggg-omp: cannot allocate edge-correction workspace");
        return FAILURE;
    }

    for (n1=1; n1<=cmd->sizeHistN; n1++) {
        for (n2=1; n2<=cmd->sizeHistN; n2++) {
            ggg_complex nzero = ggg_histogram_mode(
                gdlN->histZetaMcos, gdlN->histZetaMsin,
                gdlN->histZetaMsincos, gdlN->histZetaMcossin,
                0, n1, n2);
            real tolerance;
            real matrix_scale = 0.0;
            bool singular = FALSE;
            int row, column;

            memset(matrix, 0, matrix_count*sizeof(*matrix));
            memset(rhs, 0, (size_t)multipoles*sizeof(*rhs));
            if (ggg_complex_abs2(nzero) <= DBL_MIN)
                continue;

            for (row=0; row<multipoles; row++) {
                int ell = row-nmax;
                for (column=0; column<multipoles; column++) {
                    int order = ell-(column-nmax);
                    ggg_complex value = ggg_histogram_mode(
                        gdlN->histZetaMcos, gdlN->histZetaMsin,
                        gdlN->histZetaMsincos, gdlN->histZetaMcossin,
                        order, n1, n2);
                    matrix[(size_t)row*multipoles+column] =
                        ggg_complex_div(value, nzero);
                    if (ggg_complex_abs2(
                            matrix[(size_t)row*multipoles+column])
                        > matrix_scale)
                        matrix_scale = ggg_complex_abs2(
                            matrix[(size_t)row*multipoles+column]);
                }
                rhs[row] = ggg_complex_div(ggg_histogram_mode(
                    gdl->histZetaMcos, gdl->histZetaMsin,
                    gdl->histZetaMsincos, gdl->histZetaMcossin,
                    ell, n1, n2), nzero);
            }
            tolerance = 128.0*DBL_EPSILON*(1.0+sqrt(matrix_scale));

            for (column=0; column<multipoles; column++) {
                int pivot = column;
                real pivot_norm = ggg_complex_abs2(
                    matrix[(size_t)column*multipoles+column]);
                int candidate;
                for (candidate=column+1; candidate<multipoles; candidate++) {
                    real candidate_norm = ggg_complex_abs2(
                        matrix[(size_t)candidate*multipoles+column]);
                    if (candidate_norm > pivot_norm) {
                        pivot = candidate;
                        pivot_norm = candidate_norm;
                    }
                }
                if (pivot_norm <= tolerance*tolerance) {
                    singular = TRUE;
                    break;
                }
                if (pivot != column) {
                    int k;
                    for (k=0; k<multipoles; k++) {
                        ggg_complex swap = matrix[
                            (size_t)column*multipoles+k];
                        matrix[(size_t)column*multipoles+k] = matrix[
                            (size_t)pivot*multipoles+k];
                        matrix[(size_t)pivot*multipoles+k] = swap;
                    }
                    {
                        ggg_complex swap = rhs[column];
                        rhs[column] = rhs[pivot];
                        rhs[pivot] = swap;
                    }
                }
                {
                    ggg_complex divisor = matrix[
                        (size_t)column*multipoles+column];
                    int k;
                    for (k=0; k<multipoles; k++)
                        matrix[(size_t)column*multipoles+k] =
                            ggg_complex_div(
                                matrix[(size_t)column*multipoles+k], divisor);
                    rhs[column] = ggg_complex_div(rhs[column], divisor);
                }
                for (row=0; row<multipoles; row++) {
                    ggg_complex factor;
                    int k;
                    if (row == column)
                        continue;
                    factor = matrix[(size_t)row*multipoles+column];
                    if (ggg_complex_abs2(factor) == 0.0)
                        continue;
                    for (k=0; k<multipoles; k++)
                        matrix[(size_t)row*multipoles+k] = ggg_complex_sub(
                            matrix[(size_t)row*multipoles+k],
                            ggg_complex_mul(factor,
                                matrix[(size_t)column*multipoles+k]));
                    rhs[row] = ggg_complex_sub(
                        rhs[row], ggg_complex_mul(factor, rhs[column]));
                }
            }
            if (singular) {
                singular_bins++;
                continue;
            }
            for (row=nmax; row<multipoles; row++) {
                int output_index = row-nmax+1;
                gdl->histZetaM_EE[output_index][n1][n2] = rhs[row].re;
                gdl->histZetaM_EE_Im[output_index][n1][n2] = rhs[row].im;
            }
        }
    }

    if (singular_bins > 0)
        verb_print_normal_info(
            cmd->verbose, cmd->verbose_log, gd->outlog,
            "octree-ggg-omp: %d edge-correction radial-bin pairs "
            "were singular and set to zero\n", singular_bins);
    free(matrix);
    free(rhs);
    return SUCCESS;
}
#endif

#ifdef SMOOTHPIVOT
#ifdef DEBUG
local char pivotsfilePath[MAXLENGTHOFFILES];
local FILE *outpivots;
#endif
#endif

#include <pthread.h>

/* Keep floating-point reduction groups independent of the OpenMP team size.
 * Tree walks within a chunk run sequentially on one worker, while different
 * chunks run concurrently and are published in ascending pivot order.
 *
 * Both values can be overridden with CPPFLAGS for focused benchmarking.  The
 * MPI claim is deliberately independent of the deterministic OpenMP chunk. */
#ifndef GGG_OMP_PIVOT_CHUNK_SIZE
#define GGG_OMP_PIVOT_CHUNK_SIZE 4096
#endif
#ifndef GGG_MPI_PIVOT_CLAIM_SIZE
#define GGG_MPI_PIVOT_CLAIM_SIZE 65536
#endif

#if GGG_OMP_PIVOT_CHUNK_SIZE < 1
#error GGG_OMP_PIVOT_CHUNK_SIZE must be positive
#endif
#if GGG_MPI_PIVOT_CLAIM_SIZE < 1
#error GGG_MPI_PIVOT_CLAIM_SIZE must be positive
#endif

#ifdef OCTREEGGGMPI
local int reduce_octree_ggg_histograms(
    struct cmdline_data *cmd, struct global_data *gd,
    gdlptr_sincos_omp_ggg gdl, bool run_threepcf,
#ifdef NMultipoles
    gdlptr_sincos_omp_ggg_N gdlN,
#endif
    INTEGER *ipmask,
#ifdef SMOOTHPIVOT
    INTEGER *ipfalse, INTEGER *count_rmin, INTEGER *count_overlap
#else
    INTEGER *unused1, INTEGER *unused2, INTEGER *unused3
#endif
    )
{
    const size_t bins = (size_t)cmd->sizeHistN;
    size_t count = bins;
    size_t cursor = 0;
    real *packed;
    INTEGER counters[7] = {gd->nbbcalc, gd->nbccalc, gd->ncccalc,
                           0, 0, 0, 0};

#ifdef TWOPCF
    if (bins > (SIZE_MAX - count) / 4) goto overflow;
    count += 4 * bins;
#ifdef SMOOTHPIVOT
    if (bins > SIZE_MAX - count) goto overflow;
    count += bins;
#endif
#endif
#ifdef THREEPCFCONVERGENCE
    if (run_threepcf) {
        const size_t numerator_orders = (size_t)cmd->mChebyshev + 1;
        size_t plane;
        size_t numerator_count;
#ifdef NMultipoles
        const size_t window_orders = (size_t)GGG_WINDOW_ORDERS(cmd);
        size_t window_count;
#endif
        if ((bins != 0 && bins > SIZE_MAX / bins)
            || (plane = bins * bins,
                numerator_orders != 0 && plane > SIZE_MAX / numerator_orders)
            || (numerator_count = numerator_orders * plane,
                numerator_count > (SIZE_MAX - count) / 4))
            goto overflow;
        count += 4 * numerator_count;
#ifdef NMultipoles
        if ((window_orders != 0 && plane > SIZE_MAX / window_orders)
            || (window_count = window_orders * plane,
                window_count > (SIZE_MAX - count) / 4))
            goto overflow;
        count += 4 * window_count;
#endif
    }
#endif
    if (count > SIZE_MAX / sizeof(*packed)) goto overflow;

    if (cballs_opt_read_mask(cmd)) counters[3] = *ipmask;
#ifndef SMOOTHPIVOT
    (void)unused1;
    (void)unused2;
    (void)unused3;
#else
    counters[4] = *ipfalse;
    counters[5] = *count_rmin;
    counters[6] = *count_overlap;
#endif

    packed = malloc(count * sizeof(*packed));
    if (fcfc_octree_ggg_mpi_consensus(
            cmd, packed == NULL ? FAILURE : SUCCESS,
            "MPI octree-GGG histogram packing allocation") == FAILURE) {
        free(packed);
        return FAILURE;
    }

#define PACK_VECTOR(array)                                                 \
    do {                                                                  \
        for (int n = 1; n <= cmd->sizeHistN; n++)                        \
            packed[cursor++] = (array)[n];                               \
    } while (0)
#define PACK_ZETA(array, orders)                                           \
    do {                                                                  \
        for (int m = 1; m <= (orders); m++)                              \
            for (int n = 1; n <= cmd->sizeHistN; n++)                    \
                for (int l = 1; l <= cmd->sizeHistN; l++)                \
                    packed[cursor++] = (array)[m][n][l];                  \
    } while (0)

    PACK_VECTOR(gdl->histNNSub);
#ifdef TWOPCF
    PACK_VECTOR(gdl->histNN);
    PACK_VECTOR(gdl->histWW);
    PACK_VECTOR(gdl->histNNSubXi2pcf);
    PACK_VECTOR(gdl->histXi2pcf);
#ifdef SMOOTHPIVOT
    PACK_VECTOR(gdl->histNNSubXi2pcftotal);
#endif
#endif
#ifdef THREEPCFCONVERGENCE
    if (run_threepcf) {
    PACK_ZETA(gdl->histZetaMcos, cmd->mChebyshev + 1);
    PACK_ZETA(gdl->histZetaMsin, cmd->mChebyshev + 1);
    PACK_ZETA(gdl->histZetaMsincos, cmd->mChebyshev + 1);
    PACK_ZETA(gdl->histZetaMcossin, cmd->mChebyshev + 1);
#ifdef NMultipoles
    PACK_ZETA(gdlN->histZetaMcos, GGG_WINDOW_ORDERS(cmd));
    PACK_ZETA(gdlN->histZetaMsin, GGG_WINDOW_ORDERS(cmd));
    PACK_ZETA(gdlN->histZetaMsincos, GGG_WINDOW_ORDERS(cmd));
    PACK_ZETA(gdlN->histZetaMcossin, GGG_WINDOW_ORDERS(cmd));
#endif
    }
#endif

    if (cursor != count
        || fcfc_octree_ggg_mpi_reduce_reals(cmd, packed, count) == FAILURE
        || fcfc_octree_ggg_mpi_reduce_integers(cmd, counters, 7) == FAILURE) {
        free(packed);
        if (cursor != count)
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     "octree-ggg-mpi internal histogram packing mismatch");
        return FAILURE;
    }

    if (fcfc_octree_ggg_mpi_is_root()) {
        cursor = 0;
#define UNPACK_VECTOR(array)                                               \
        do {                                                              \
            for (int n = 1; n <= cmd->sizeHistN; n++)                    \
                (array)[n] = packed[cursor++];                            \
        } while (0)
#define UNPACK_ZETA(array, orders)                                         \
        do {                                                              \
            for (int m = 1; m <= (orders); m++)                          \
                for (int n = 1; n <= cmd->sizeHistN; n++)                \
                    for (int l = 1; l <= cmd->sizeHistN; l++)            \
                        (array)[m][n][l] = packed[cursor++];              \
        } while (0)
        UNPACK_VECTOR(gdl->histNNSub);
#ifdef TWOPCF
        UNPACK_VECTOR(gdl->histNN);
        UNPACK_VECTOR(gdl->histWW);
        UNPACK_VECTOR(gdl->histNNSubXi2pcf);
        UNPACK_VECTOR(gdl->histXi2pcf);
#ifdef SMOOTHPIVOT
        UNPACK_VECTOR(gdl->histNNSubXi2pcftotal);
#endif
#endif
#ifdef THREEPCFCONVERGENCE
        if (run_threepcf) {
        UNPACK_ZETA(gdl->histZetaMcos, cmd->mChebyshev + 1);
        UNPACK_ZETA(gdl->histZetaMsin, cmd->mChebyshev + 1);
        UNPACK_ZETA(gdl->histZetaMsincos, cmd->mChebyshev + 1);
        UNPACK_ZETA(gdl->histZetaMcossin, cmd->mChebyshev + 1);
#ifdef NMultipoles
        UNPACK_ZETA(gdlN->histZetaMcos, GGG_WINDOW_ORDERS(cmd));
        UNPACK_ZETA(gdlN->histZetaMsin, GGG_WINDOW_ORDERS(cmd));
        UNPACK_ZETA(gdlN->histZetaMsincos, GGG_WINDOW_ORDERS(cmd));
        UNPACK_ZETA(gdlN->histZetaMcossin, GGG_WINDOW_ORDERS(cmd));
#endif
        }
#endif
        gd->nbbcalc = counters[0];
        gd->nbccalc = counters[1];
        gd->ncccalc = counters[2];
        if (cballs_opt_read_mask(cmd)) *ipmask = counters[3];
#ifdef SMOOTHPIVOT
        *ipfalse = counters[4];
        *count_rmin = counters[5];
        *count_overlap = counters[6];
#endif
#undef UNPACK_VECTOR
#undef UNPACK_ZETA
    }
    free(packed);
#undef PACK_VECTOR
#undef PACK_ZETA
    return SUCCESS;

overflow:
    snprintf(cmd->error_message, _ERRORMSGSIZE_,
             "octree-ggg-mpi histogram packing size overflow");
    return FAILURE;
}
#endif

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
local int searchcalc_octree_ggg_driver(struct cmdline_data* cmd,
                                       struct  global_data* gd,
                                       bodyptr *btable, INTEGER *nbody,
                                       INTEGER ipmin, INTEGER *ipmax,
                                       int cat1, int cat2, bool distributed)
{
    string routineName = distributed ? "searchcalc_octree_ggg_mpi"
                                     : "searchcalc_octree_ggg_omp";
    double cpustart;
    gdl_sincos_omp_ggg gdl;
#ifdef NMultipoles
    gdl_sincos_omp_ggg_N gdlN;
#endif
#ifndef OCTREEGGGMPI
    (void)distributed;
#endif

    cpustart = CPUTIME;
    gd->cpusearch = 0.0;
    const bool run_threepcf = !cballs_opt_only_2pcf(cmd);
#ifndef TWOPCF
    if (!run_threepcf) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: only-2pcf requires TWOPCFON=1", routineName);
        return FAILURE;
    }
#endif
#ifdef OCTREEGGGMPI
    if (distributed && !fcfc_octree_ggg_mpi_active()) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: MPI runtime is not initialized", routineName);
        return FAILURE;
    }
#ifdef DEBUG
    if (distributed) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: DEBUG pivot-file output is not MPI-safe", routineName);
        return FAILURE;
    }
#endif
#endif
    if (print_info(cmd, gd) == FAILURE)
        return FAILURE;

    if (cmd->useLogHist==FALSE &&
        (strcmp(cmd->searchMethod,"octree-ggg-omp") == 0
         || strcmp(cmd->searchMethod,"octree-ggg-mpi") == 0))
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
            "%s: can´t have loghist false and octree-ggg-omp (%d %s)\n",
            routineName, cmd->useLogHist, cmd->searchMethod);

#ifdef TPCFSHEAR
    mCheb =cmd->mChebyshev + mpOffSet;
#endif

#ifdef SMOOTHPIVOT
#ifdef DEBUG
    if (format_checked(pivotsfilePath, sizeof(pivotsfilePath),
        "pivotsfilePath", "%s/pivot_info%s.txt",
                       gd->tmpDir,cmd->suffixOutFiles) != 0)
        return FAILURE;
    if(!(outpivots=fopen(pivotsfilePath, "w"))) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "\n%s: error opening file '%s' \n",
                       routineName, pivotsfilePath);
        return FAILURE;
    }
#endif
#endif

#ifdef OPENMPCODE
    ThreadCount(cmd, gd, nbody[cat1], cat1);
    pthread_t main_thread_id = pthread_self();
#else
#error `OPENMPMACHINE` is not defined. Switch it on in Makefile_settings
#endif

    if (cballs_opt_pivot_number(cmd)) {
        ipmax[cat1] = gd->pivotNumber;
        verb_print(cmd->verbose, "\n%s: pivot number: %ld\n",
                   routineName, gd->pivotNumber);
    }

#ifdef SMOOTHPIVOT
    int smooth_status = prepare_smooth_pivots(cmd, gd, btable, nbody,
                                               ipmin, ipmax, cat1, cat2);
#ifdef OCTREEGGGMPI
    if (distributed)
        smooth_status = fcfc_octree_ggg_mpi_consensus(
            cmd, smooth_status, "MPI smooth-pivot preparation");
#endif
    if (smooth_status == FAILURE)
        return FAILURE;
#endif

    int gdl_local_status = search_init_gd_sincos_omp_ggg(cmd, gd, &gdl);
    int gdl_status = gdl_local_status;
#ifdef OCTREEGGGMPI
    if (distributed)
        gdl_status = fcfc_octree_ggg_mpi_consensus(
            cmd, gdl_status, "MPI octree-GGG histogram initialization");
#endif
    if (gdl_status == FAILURE) {
        if (gdl_local_status == SUCCESS)
            search_free_gd_sincos_omp_ggg(cmd, gd, &gdl);
#ifdef DEBUG
        fclose(outpivots);
#endif
        return FAILURE;
    }
#ifdef NMultipoles
    int gdln_local_status = SUCCESS;
    if (run_threepcf)
        gdln_local_status =
            search_init_gd_sincos_omp_ggg_N(cmd, gd, &gdlN);
    int gdln_status = gdln_local_status;
#ifdef OCTREEGGGMPI
    if (distributed)
        gdln_status = fcfc_octree_ggg_mpi_consensus(
            cmd, gdln_status,
            "MPI octree-GGG N-multipole histogram initialization");
#endif
    if (gdln_status == FAILURE) {
        if (run_threepcf && gdln_local_status == SUCCESS)
            search_free_gd_sincos_omp_ggg_N(cmd, gd, &gdlN);
        search_free_gd_sincos_omp_ggg(cmd, gd, &gdl);
#ifdef DEBUG
        fclose(outpivots);
#endif
        return FAILURE;
    }
#endif

    //B Mask correction
    INTEGER ipmask = 0;
    //E

#ifdef SMOOTHPIVOT
    INTEGER ipfalse;
    ipfalse=0;
    INTEGER icountNbRmin;
    icountNbRmin=0;
    INTEGER icountNbRminOverlap;
    icountNbRminOverlap=0;
#endif

#ifndef BALLS4SCANLEV
    INTEGER pivot_task_count = ipmax[cat1] - ipmin + 1;
#else
    INTEGER pivot_task_count = gd->nnodescanlevTableB4[cat1];
#endif
    INTEGER task_first = 0;
    INTEGER task_last = pivot_task_count;
    int scheduler_status = SUCCESS;
    int work_done = FALSE;
#ifdef OCTREEGGGMPI
    fcfc_octree_ggg_mpi_scheduler scheduler;
    memset(&scheduler, 0, sizeof(scheduler));
    if (distributed) {
        int thread_hint = omp_get_max_threads();
        INTEGER step = (INTEGER)GGG_MPI_PIVOT_CLAIM_SIZE;
        if (fcfc_octree_ggg_mpi_scheduler_init(
                cmd, &scheduler, pivot_task_count, step) == FAILURE) {
#ifdef NMultipoles
            if (run_threepcf)
                search_free_gd_sincos_omp_ggg_N(cmd, gd, &gdlN);
#endif
            search_free_gd_sincos_omp_ggg(cmd, gd, &gdl);
#ifdef DEBUG
            fclose(outpivots);
#endif
            return FAILURE;
        }
        verb_print(cmd->verbose,
                   "octree-ggg-mpi: %d ranks, %d OpenMP threads, "
                   "%" INTEGER_FMT " pivot tasks\n",
                   fcfc_octree_ggg_mpi_size(), thread_hint,
                   pivot_task_count);
    }
#endif

#if defined(NMultipoles) && defined(NONORMHIST)
    if (cballs_opt_patch_with_all(cmd)) {
        verb_print(cmd->verbose,
            "\n%s: total number of pixels to be pivots: %ld\n",
                   routineName, gd->pivotCount);
    }
#endif

    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                           "\n%s: Total allocated %g MByte storage so far.\n",
                           routineName, gd->bytes_tot*INMB);

//B sometimes when running happens a "Floating point exception: 8" error
//      indicates a division by zero error within a program.
//          numeric code "8" specifically points to integer division
//          by zero as the root cause
//    verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
//                        "\n%s: stepState: %ld\n", routineName, gd->stepState);
//E

    verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                        "\nRunning...\n - Completed pivot node:\n");
    int allocation_failed = FALSE;

//
// Check that all posibilities are taken in to account...
//
#ifdef DEBUG
#ifndef SMOOTHPIVOT
    snprintf(cmd->error_message, _ERRORMSGSIZE_,
             "%s: DEBUG can not be used with SMOOTHPIVOT turned OFF.",
             routineName);
    return FAILURE;
#else
#pragma omp parallel default(shared)                                        \
    shared(cmd,gd,btable,nbody,roottable,outpivots,                         \
           ipmin,ipmax,cat1,cat2,ipfalse,ipmask,                            \
           icountNbRmin,icountNbRminOverlap, gdl, main_thread_id, \
           allocation_failed)
#endif // ! SMOOTHPIVOT
#else // ! DEBUG
#ifdef NMultipoles
#ifndef BALLS4SCANLEV
#ifdef SMOOTHPIVOT
#pragma omp parallel default(shared)                                        \
    shared(cmd,gd,btable,nbody,roottable,                                   \
           ipmin,ipmax,cat1,cat2,ipfalse,ipmask,                            \
           icountNbRmin,icountNbRminOverlap,                                \
           gdl, gdlN, main_thread_id, allocation_failed)
#else // ! SMOOTHPIVOT
#pragma omp parallel default(shared)                                        \
    shared(cmd,gd,btable,nbody,roottable,                                   \
           ipmin,ipmax,cat1,cat2,ipmask,                                    \
           gdl, gdlN, main_thread_id, allocation_failed)
#endif // ! SMOOTHPIVOT
#else // ! BALLS4SCANLEV
#ifdef SMOOTHPIVOT
#pragma omp parallel default(shared)                                        \
    shared(cmd,gd,btable,nbody,roottable,nodetablescanlevB4,                \
           ipmin,ipmax,cat1,cat2,ipfalse,ipmask,                            \
           icountNbRmin,icountNbRminOverlap, gdl, gdlN, main_thread_id, \
           allocation_failed)
#else // ! SMOOTHPIVOT
#pragma omp parallel default(shared)                                        \
    shared(cmd,gd,btable,nbody,roottable,nodetablescanlevB4,                \
           ipmin,ipmax,cat1,cat2,ipmask, gdl, gdlN, main_thread_id, \
           allocation_failed)
#endif // ! SMOOTHPIVOT
#endif // ! BALLS4SCANLEV
#else // ! NMultipoles
#ifndef BALLS4SCANLEV
#ifdef SMOOTHPIVOT
#pragma omp parallel default(shared)                                        \
    shared(cmd,gd,btable,nbody,roottable,                                   \
           ipmin,ipmax,cat1,cat2,ipfalse,ipmask,                            \
           icountNbRmin,icountNbRminOverlap, gdl, main_thread_id, \
           allocation_failed)
#else // ! SMOOTHPIVOT
#pragma omp parallel default(shared)                                        \
    shared(cmd,gd,btable,nbody,roottable,                                   \
           ipmin,ipmax,cat1,cat2,ipmask, gdl, main_thread_id, \
           allocation_failed)
#endif // ! SMOOTHPIVOT
#else // ! BALLS4SCANLEV
#ifdef SMOOTHPIVOT
#pragma omp parallel default(shared)                                        \
    shared(cmd,gd,btable,nbody,roottable,nodetablescanlevB4,                \
           ipmin,ipmax,cat1,cat2,ipfalse,ipmask,                            \
           icountNbRmin,icountNbRminOverlap, gdl, main_thread_id, \
           allocation_failed)
#else // ! SMOOTHPIVOT
#pragma omp parallel default(shared)                                        \
    shared(cmd,gd,btable,nbody,roottable,nodetablescanlevB4,                \
           ipmin,ipmax,cat1,cat2,ipmask, gdl, main_thread_id, \
           allocation_failed)
#endif // ! SMOOTHPIVOT
#endif // ! BALLS4SCANLEV
#endif // ! NMultipoles
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
    int hist_ready =
        search_init_sincos_omp_ggg(cmd, gd, &hist) == SUCCESS;
    int worker_ready = hist_ready;
#ifdef NMultipoles
    gdhist_sincos_omp_ggg_N histN;
    int hist_n_ready = !run_threepcf;
    if (run_threepcf)
        hist_n_ready = hist_ready &&
            search_init_sincos_omp_ggg_N(cmd, gd, &histN) == SUCCESS;
    worker_ready = worker_ready && hist_n_ready;
#endif
    if (!worker_ready) {
#pragma omp atomic write
        allocation_failed = TRUE;
    }

#pragma omp barrier
      if (main_thread_id == current_thread_id)
          verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                "\n\t Total allocated %g MByte storage so far (threads).\n",
                gd->bytes_tot*INMB);
    //E

    INTEGER ipmaskthreads = 0;

#ifdef SMOOTHPIVOT
    INTEGER ipfalsethreads;
    ipfalsethreads = 0;
    INTEGER icountNbRminthread;
    icountNbRminthread=0;
    INTEGER icountNbRminOverlapthread;
    icountNbRminOverlapthread=0;

#endif

#ifdef OCTREEGGGMPI
      if (distributed) {
#pragma omp master
        {
          if (fcfc_octree_ggg_mpi_consensus(
                  cmd, allocation_failed ? FAILURE : SUCCESS,
                  "MPI OpenMP histogram allocation") == FAILURE)
              allocation_failed = TRUE;
        }
#pragma omp barrier
      }
#endif

      for (;;) {
#pragma omp master
        {
          if (!distributed) {
              task_first = 0;
              task_last = pivot_task_count;
              scheduler_status = SUCCESS;
          }
#ifdef OCTREEGGGMPI
          else {
              scheduler_status = fcfc_octree_ggg_mpi_scheduler_claim(
                  cmd, &scheduler, &task_first, &task_last);
          }
#endif
          if (scheduler_status == FAILURE) allocation_failed = TRUE;
          work_done = task_first == task_last;
        }
#pragma omp barrier
        if (allocation_failed || work_done) break;

      INTEGER chunk_count =
          (task_last - task_first + GGG_OMP_PIVOT_CHUNK_SIZE - 1)
          / GGG_OMP_PIVOT_CHUNK_SIZE;
#pragma omp for schedule(dynamic, 1) ordered
      for (INTEGER ichunk = 0; ichunk < chunk_count; ichunk++) {
          INTEGER chunk_first = task_first
              + ichunk * GGG_OMP_PIVOT_CHUNK_SIZE;
          INTEGER chunk_last = MIN(
              chunk_first + GGG_OMP_PIVOT_CHUNK_SIZE, task_last);

          /* These are chunk accumulators.  Clearing them here gives every
           * run exactly the same floating-point grouping. */
#ifdef TWOPCF
          for (n = 1; n <= cmd->sizeHistN; n++) {
              hist.histNthread[n] = 0.0;
              hist.histWWthread[n] = 0.0;
              hist.histNNSubXi2pcfthread[n] = 0.0;
#ifdef SMOOTHPIVOT
              hist.histNNSubXi2pcfthreadtotal[n] = 0.0;
#endif
              hist.histXi2pcfthread[n] = 0.0;
          }
#endif
          for (n = 1; n <= cmd->sizeHistN; n++)
              hist.histNNSubthread[n] = 0.0;
#ifdef THREEPCFCONVERGENCE
          if (run_threepcf) {
          for (m = 1; m <= cmd->mChebyshev+1; m++) {
              CLRM_ext(hist.histZetaMthreadcos[m], cmd->sizeHistN);
              CLRM_ext(hist.histZetaMthreadsin[m], cmd->sizeHistN);
              CLRM_ext(hist.histZetaMthreadsincos[m], cmd->sizeHistN);
              CLRM_ext(hist.histZetaMthreadcossin[m], cmd->sizeHistN);
          }
#ifdef NMultipoles
          for (m = 1; m <= GGG_WINDOW_ORDERS(cmd); m++) {
              CLRM_ext(histN.histZetaMthreadcos[m], cmd->sizeHistN);
              CLRM_ext(histN.histZetaMthreadsin[m], cmd->sizeHistN);
              CLRM_ext(histN.histZetaMthreadsincos[m], cmd->sizeHistN);
              CLRM_ext(histN.histZetaMthreadcossin[m], cmd->sizeHistN);
          }
#endif
          }
#endif
          hist.nbbcalcthread = 0;
          hist.nbccalcthread = 0;
          ipmaskthreads = 0;
#ifdef SMOOTHPIVOT
          ipfalsethreads = 0;
          icountNbRminthread = 0;
          icountNbRminOverlapthread = 0;
#endif

          for (INTEGER itask = chunk_first; itask < chunk_last; itask++) {
#ifndef BALLS4SCANLEV
              p = btable[cat1] + ipmin - 1 + itask;
#else
              i = itask;
              p = nodetablescanlevB4[cat1][i];
#endif
          if (allocation_failed) continue;
#ifdef SMOOTHPIVOT
          if (Update(p) == FALSE) {
              ipfalsethreads++;
              continue;
          }
#endif

          if (cballs_opt_read_mask(cmd)) {
            if (Mask(p) == FALSE) {
                ipmaskthreads++;
                continue;
            }
          }

#if defined(NMultipoles) && defined(NONORMHIST)
          if (run_threepcf && cballs_opt_patch_with_all(cmd)) {
              if (UpdatePivot(p) == FALSE) {
                  continue;
              }
          }
#endif

//B segment to be included below...
#ifdef TWOPCF
          for (n = 1; n <= cmd->sizeHistN; n++) {
              hist.histWthread[n] = 0.0;
              hist.histXi2pcfthreadsub[n] = 0.0;
#ifdef SMOOTHPIVOT
              hist.histNNSubXi2pcfthreadp[n] = 0.;
#endif
          }
#endif

          for (n = 1; n <= cmd->sizeHistN; n++)
              hist.histNNSubthread[n] = 0.0;
          //B 3pcf convergence & shear counting
#ifdef NMultipoles
          if (run_threepcf)
              for (n = 1; n <= cmd->sizeHistN; n++)
                  histN.histNNSubthread[n] = 0.0;
#endif

          //B Set histograms to zero for the pivot
#ifdef THREEPCFCONVERGENCE
          if (run_threepcf) {
          CLRM_ext_ext(hist.histXithreadcos,
                       cmd->mChebyshev+1, cmd->sizeHistN);
          CLRM_ext_ext(hist.histXithreadsin,
                       cmd->mChebyshev+1, cmd->sizeHistN);
          CLRM_ext_ext(hist.histXithreaddiagcos,
                       cmd->mChebyshev+1, cmd->sizeHistN);
          CLRM_ext_ext(hist.histXithreaddiagsin,
                       cmd->mChebyshev+1, cmd->sizeHistN);
          CLRM_ext_ext(hist.histXithreaddiagsincos,
                       cmd->mChebyshev+1, cmd->sizeHistN);
#endif

#ifdef THREEPCFCONVERGENCE
#ifdef NMultipoles
          //B 3pcf convergence & shear counting
          CLRM_ext_ext(histN.histXithreadcos,
                       GGG_WINDOW_ORDERS(cmd), cmd->sizeHistN);
          CLRM_ext_ext(histN.histXithreadsin,
                       GGG_WINDOW_ORDERS(cmd), cmd->sizeHistN);
          CLRM_ext_ext(histN.histXithreaddiagcos,
                       GGG_WINDOW_ORDERS(cmd), cmd->sizeHistN);
          CLRM_ext_ext(histN.histXithreaddiagsin,
                       GGG_WINDOW_ORDERS(cmd), cmd->sizeHistN);
          CLRM_ext_ext(histN.histXithreaddiagsincos,
                       GGG_WINDOW_ORDERS(cmd), cmd->sizeHistN);
          //E
#endif
          }
#endif

#ifdef TPCFSHEAR
          if (run_threepcf) {
          CLRM_ext_ext(hist.histg1threadcos,
                       cmd->mChebyshev+1, cmd->sizeHistN);
          CLRM_ext_ext(hist.histg1threadsin,
                       cmd->mChebyshev+1, cmd->sizeHistN);
          CLRM_ext_ext(hist.histg2threadcos,
                       cmd->mChebyshev+1, cmd->sizeHistN);
          CLRM_ext_ext(hist.histg2threadsin,
                       cmd->mChebyshev+1, cmd->sizeHistN);
          }
#endif
          //E

#if defined(THREEPCFCONVERGENCE) || defined(TPCFSHEAR)
          //B 3pcf convergence & shear
          if (run_threepcf) {
#ifndef BALLS4SCANLEV
          polarix_init(cmd, gd, p, &hist);
#else
          polarix_init(cmd, gd, (bodyptr)p, &hist);
#endif
          }
          //E 3pcf convergence & shear
#endif
//E segment to be included below...

//================
#ifdef NMultipoles
//================

#ifndef BALLS4SCANLEV
          if (run_threepcf)
              normal_walktree_sincos_N(cmd, gd, btable, cat2,
                                       p, ((nodeptr) roottable[cat2]),
                                       gd->rSizeTable[cat2], &hist, &histN);
          else
              normal_walktree_sincos(cmd, gd, btable, cat2,
                                     p, ((nodeptr) roottable[cat2]),
                                     gd->rSizeTable[cat2], &hist);

#ifdef TWOPCF
          //B kappa Avg Rmin
#ifdef SMOOTHPIVOT
          for (n = 1; n <= cmd->sizeHistN; n++) {
              hist.histNNSubXi2pcfthreadp[n] =
                        ((real)NbRmin(p))*hist.histNNSubXi2pcfthreadp[n];
              hist.histNNSubXi2pcfthreadtotal[n] +=
                        hist.histNNSubXi2pcfthreadp[n];
// Check if these lines are needed!!! they are in not a BALLS4SCANLEV domain
//                  hist.histNNSubthread[n] =
//                                    ((real)NbRmin(p))*hist.histNNSubthread[n];
          }
#endif
          //E
#endif

#ifdef SMOOTHPIVOT
          for (n = 1; n <= cmd->sizeHistN; n++)
              // Check if these lines are needed!!!
              hist.histNNSubthread[n] =
                    ((real)NbRmin(p))*hist.histNNSubthread[n];
#endif

          computeBodyProperties_sincos_ggg(cmd, gd, p,
                                           ipmax[cat1]-ipmin+1, &hist);
          if (run_threepcf)
              computeBodyProperties_sincos_ggg_N(
                  cmd, gd, p, ipmax[cat1]-ipmin+1, &histN);
#else // ! BALLS4SCANLEV
          if (run_threepcf)
              normal_walktree_sincos_N(cmd, gd, btable, cat2,
                                       (bodyptr)p,
                                       ((nodeptr) roottable[cat2]),
                                       gd->rSizeTable[cat2], &hist, &histN);
          else
              normal_walktree_sincos(cmd, gd, btable, cat2,
                                     (bodyptr)p,
                                     ((nodeptr) roottable[cat2]),
                                     gd->rSizeTable[cat2], &hist);

#ifdef TWOPCF
#ifdef SMOOTHPIVOT
          for (n = 1; n <= cmd->sizeHistN; n++) {
              hist.histNNSubXi2pcfthreadp[n] =
                        ((real)NbRmin(p))*hist.histNNSubXi2pcfthreadp[n];
              hist.histNNSubXi2pcfthreadtotal[n] +=
                        hist.histNNSubXi2pcfthreadp[n];
// Check if these lines are needed!!! they are in a BALLS4SCANLEV domain
//                  hist.histNNSubthread[n] =
//                                    ((real)NbRmin(p))*hist.histNNSubthread[n];
          }
#endif
          //E
#endif

#ifdef SMOOTHPIVOT
          for (n = 1; n <= cmd->sizeHistN; n++)
        // Check if these lines are needed!!! they are in a BALLS4SCANLEV domain
                hist.histNNSubthread[n] =
                                ((real)NbRmin(p))*hist.histNNSubthread[n];
#endif

          computeBodyProperties_sincos_ggg(cmd, gd, (bodyptr)p,
                                           gd->nnodescanlevTableB4[cat1], &hist);
          if (run_threepcf)
              computeBodyProperties_sincos_ggg_N(
                  cmd, gd, (bodyptr)p,
                  gd->nnodescanlevTableB4[cat1], &histN);
#endif // ! BALLS4SCANLEV

//================
#else // ! NMultipoles
//================

#ifndef BALLS4SCANLEV
          normal_walktree_sincos(cmd, gd, btable, cat2,
                                 p, ((nodeptr) roottable[cat2]),
                                 gd->rSizeTable[cat2], &hist);
#ifdef TWOPCF
#ifdef SMOOTHPIVOT
          for (n = 1; n <= cmd->sizeHistN; n++) {
              hist.histNNSubXi2pcfthreadp[n] =
                        ((real)NbRmin(p))*hist.histNNSubXi2pcfthreadp[n];
              hist.histNNSubXi2pcfthreadtotal[n] +=
                        hist.histNNSubXi2pcfthreadp[n];
// Check if these lines are needed!!! they are in not a BALLS4SCANLEV domain
//                  hist.histNNSubthread[n] =
//                                    ((real)NbRmin(p))*hist.histNNSubthread[n];
          }
#endif
#endif

#ifdef SMOOTHPIVOT
          for (n = 1; n <= cmd->sizeHistN; n++)
              // Check if these lines are needed!!!
                hist.histNNSubthread[n] =
                            ((real)NbRmin(p))*hist.histNNSubthread[n];
#endif

          computeBodyProperties_sincos_ggg(cmd, gd, p,
                                           ipmax[cat1]-ipmin+1, &hist);
#else // ! BALLS4SCANLEV
          normal_walktree_sincos(cmd, gd, btable, cat2,
                                 (bodyptr)p, ((nodeptr) roottable[cat2]),
                                 gd->rSizeTable[cat2], &hist);
#ifdef TWOPCF
#ifdef SMOOTHPIVOT
          for (n = 1; n <= cmd->sizeHistN; n++) {
              hist.histNNSubXi2pcfthreadp[n] =
                        ((real)NbRmin(p))*hist.histNNSubXi2pcfthreadp[n];
              hist.histNNSubXi2pcfthreadtotal[n] +=
                        hist.histNNSubXi2pcfthreadp[n];
// Check if these lines are needed!!! they are in a BALLS4SCANLEV domain
//                  hist.histNNSubthread[n] =
//                                    ((real)NbRmin(p))*hist.histNNSubthread[n];
          }
#endif
#endif

#ifdef SMOOTHPIVOT
          for (n = 1; n <= cmd->sizeHistN; n++)
              // Check if these lines are needed!!!
              hist.histNNSubthread[n] =
                        ((real)NbRmin(p))*hist.histNNSubthread[n];
#endif

          computeBodyProperties_sincos_ggg(cmd, gd, (bodyptr)p,
                                           gd->nnodescanlevTableB4[cat1], &hist);
#endif // ! BALLS4SCANLEV
          
//================
#endif // ! NMultipoles
//================

#ifndef BALLS4SCANLEV
          ip = p - btable[cat1] + 1;
#else
          ip = i+1;
#endif
#ifdef SMOOTHPIVOT
          icountNbRminthread += NbRmin(p);
          icountNbRminOverlapthread += NbRminOverlap(p);
#ifdef DEBUG
          fprintf(outpivots,"%ld \t%ld \t%ld \t\t%g \t\t%g\n",
                  ip, NbRmin(p), NbRminOverlap(p),
                  KappaRmin(p)/NbRmin(p), WeightRmin(p)/NbRmin(p));
#endif
#endif
          if (ip%gd->stepState == 0)
          verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog, ".");
          if (ip%cmd->stepState == 0)
          verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                              "%d ", ip);
          } // end do body p

          /* OpenMP's ordered region follows logical chunk order, regardless
           * of which worker computed the chunk. */
#pragma omp ordered
          {
#ifdef TWOPCF
            for (n = 1; n <= cmd->sizeHistN; n++) {
                gdl.histNN[n] += hist.histNthread[n];
                gdl.histWW[n] += hist.histWWthread[n];
                gdl.histNNSubXi2pcf[n] += hist.histNNSubXi2pcfthread[n];
#ifdef SMOOTHPIVOT
                gdl.histNNSubXi2pcftotal[n] +=
                    hist.histNNSubXi2pcfthreadtotal[n];
#endif
                gdl.histXi2pcf[n] += hist.histXi2pcfthread[n];
            }
#endif
            for (n = 1; n <= cmd->sizeHistN; n++)
                gdl.histNNSub[n] += hist.histNNSubthread[n];

#ifdef THREEPCFCONVERGENCE
            if (run_threepcf) {
            for (m = 1; m <= cmd->mChebyshev+1; m++) {
                ADDM_ext(gdl.histZetaMcos[m], gdl.histZetaMcos[m],
                         hist.histZetaMthreadcos[m], cmd->sizeHistN);
                ADDM_ext(gdl.histZetaMsin[m], gdl.histZetaMsin[m],
                         hist.histZetaMthreadsin[m], cmd->sizeHistN);
                ADDM_ext(gdl.histZetaMsincos[m], gdl.histZetaMsincos[m],
                         hist.histZetaMthreadsincos[m], cmd->sizeHistN);
                ADDM_ext(gdl.histZetaMcossin[m], gdl.histZetaMcossin[m],
                         hist.histZetaMthreadcossin[m], cmd->sizeHistN);
            }
#ifdef NMultipoles
            for (m = 1; m <= GGG_WINDOW_ORDERS(cmd); m++) {
                ADDM_ext(gdlN.histZetaMcos[m], gdlN.histZetaMcos[m],
                         histN.histZetaMthreadcos[m], cmd->sizeHistN);
                ADDM_ext(gdlN.histZetaMsin[m], gdlN.histZetaMsin[m],
                         histN.histZetaMthreadsin[m], cmd->sizeHistN);
                ADDM_ext(gdlN.histZetaMsincos[m], gdlN.histZetaMsincos[m],
                         histN.histZetaMthreadsincos[m], cmd->sizeHistN);
                ADDM_ext(gdlN.histZetaMcossin[m], gdlN.histZetaMcossin[m],
                         histN.histZetaMthreadcossin[m], cmd->sizeHistN);
            }
#endif
            }
#endif
            gd->nbbcalc += hist.nbbcalcthread;
            gd->nbccalc += hist.nbccalcthread;
            if (cballs_opt_read_mask(cmd))
                ipmask += ipmaskthreads;
#ifdef SMOOTHPIVOT
            ipfalse += ipfalsethreads;
            icountNbRmin += icountNbRminthread;
            icountNbRminOverlap += icountNbRminOverlapthread;
#endif
          }
      } // end deterministic chunk loop
        if (!distributed) break;
      } // end FCFC MPI claim/process loop

      if (main_thread_id == current_thread_id)
          verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog, "\n");

#ifdef NMultipoles
    if (run_threepcf && hist_n_ready)
        search_free_sincos_omp_ggg_N(cmd, gd, &histN);  // free memory
#endif
    if (hist_ready)
        search_free_sincos_omp_ggg(cmd, gd, &hist);     // free memory
  } // end pragma omp parallel

#ifdef OCTREEGGGMPI
    if (distributed) {
        if (fcfc_octree_ggg_mpi_consensus(
                cmd, allocation_failed ? FAILURE : SUCCESS,
                "MPI octree-GGG workers") == FAILURE)
            allocation_failed = TRUE;
        if (scheduler.ready
            && fcfc_octree_ggg_mpi_scheduler_destroy(
                   cmd, &scheduler) == FAILURE)
            allocation_failed = TRUE;
        if (fcfc_octree_ggg_mpi_consensus(
                cmd, allocation_failed ? FAILURE : SUCCESS,
                "MPI octree-GGG scheduler cleanup") == FAILURE)
            allocation_failed = TRUE;
    }
#endif

    if (allocation_failed) {
#ifdef DEBUG
        fclose(outpivots);
#endif
#ifdef NMultipoles
        if (run_threepcf)
            search_free_gd_sincos_omp_ggg_N(cmd, gd, &gdlN);
#endif
        search_free_gd_sincos_omp_ggg(cmd, gd, &gdl);
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: OpenMP histogram allocation failed", routineName);
        return FAILURE;
    }

#ifdef OCTREEGGGMPI
    if (distributed) {
        if (reduce_octree_ggg_histograms(
                cmd, gd, &gdl, run_threepcf,
#ifdef NMultipoles
                &gdlN,
#endif
                &ipmask,
#ifdef SMOOTHPIVOT
                &ipfalse, &icountNbRmin, &icountNbRminOverlap
#else
                NULL, NULL, NULL
#endif
                ) == FAILURE) {
#ifdef NMultipoles
            if (run_threepcf)
                search_free_gd_sincos_omp_ggg_N(cmd, gd, &gdlN);
#endif
            search_free_gd_sincos_omp_ggg(cmd, gd, &gdl);
            return FAILURE;
        }
        if (!fcfc_octree_ggg_mpi_is_root()) {
#ifdef NMultipoles
            if (run_threepcf)
                search_free_gd_sincos_omp_ggg_N(cmd, gd, &gdlN);
#endif
            search_free_gd_sincos_omp_ggg(cmd, gd, &gdl);
            return SUCCESS;
        }
    }
#endif

    //B end of completed pivot
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog, "\n\n");
    //E

#ifdef SMOOTHPIVOT
    real xi, den, num;
    int mm;
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
        xi = num/den;
#endif // ! NONORMHIST
        verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                               "%s: p falses found = %ld and %e %e %e\n",
                               routineName, ipfalse, num, den, xi);
#ifdef THREEPCFCONVERGENCE
        if (run_threepcf) {
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
        }
#endif

#ifdef THREEPCFCONVERGENCE
#ifdef NMultipoles
        if (run_threepcf) {
        for (mm=1; mm<=GGG_WINDOW_ORDERS(cmd); mm++) {
            MULMS_ext(gdlN.histZetaMcos[mm],
                      gdlN.histZetaMcos[mm],xi,cmd->sizeHistN);
            MULMS_ext(gdlN.histZetaMsin[mm],
                      gdlN.histZetaMsin[mm],xi,cmd->sizeHistN);
            MULMS_ext(gdlN.histZetaMsincos[mm],
                      gdlN.histZetaMsincos[mm],xi,cmd->sizeHistN);
            MULMS_ext(gdlN.histZetaMcossin[mm],
                      gdlN.histZetaMcossin[mm],xi,cmd->sizeHistN);
        }
        }
#endif
#endif

        verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                               "%s: p falses found = %ld\n", routineName, ipfalse);
        verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                               "%s: count NbRmin found = %ld\n",
                               routineName, icountNbRmin);
        verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                               "%s: count overlap found = %ld\n",
                               routineName, icountNbRminOverlap);

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
#endif

    if (cballs_opt_read_mask(cmd))
        verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                               "%s: p masked found = %ld\n", routineName, ipmask);

#ifdef DEBUG
    fclose(outpivots);                              // Close file to debug pivots
#endif

#ifdef TWOPCF
    //B asymmetric,weights-norm options give same results without them...
    int nn;
    if (!cballs_opt_asymmetric(cmd)) {
        for (nn = 1; nn <= cmd->sizeHistN; nn++) {
#ifdef SMOOTHPIVOT
            if (cmd->verbose>3)
                printf("%d %e %e %e\n",
                       nn, gdl.histNNSubXi2pcf[nn],
                       gdl.histNNSubXi2pcftotal[nn], gdl.histNN[nn]);
#else
            if (cmd->verbose>3)
                printf("%d %e %e\n",
                       nn, gdl.histNNSubXi2pcf[nn], gdl.histNN[nn]);
#endif
            if (cballs_opt_weights_norm(cmd)) {
//           gdl.histXi2pcf[nn] /= MAX(gdl.histNN[nn],1.0);// gives same as below
//                gdl.histXi2pcf[nn] /= MAX(gdl.histWW[nn],1.0);
                gdl.histXi2pcf[nn] /= gdl.histWW[nn];   // gives same as above
            } else {
                gdl.histXi2pcf[nn] /= 2.0;
                gdl.histNNSubXi2pcf[nn] /= 2.0;
#ifdef SMOOTHPIVOT
                gdl.histNNSubXi2pcftotal[nn] /= 2.0;
                    gdl.histXi2pcf[nn] /= MAX(gdl.histNNSubXi2pcftotal[nn],1.0);
#else
                gdl.histXi2pcf[nn] /= MAX(gdl.histNNSubXi2pcf[nn],1.0);
#endif
            }
        }
    } else {
        for (nn = 1; nn <= cmd->sizeHistN; nn++) {
//B verbose should be 0--3.
//  But there is special place (debuging) where is useful to have > 3.
#ifdef SMOOTHPIVOT
            if (cmd->verbose>3)
                verb_print_debug(cmd->verbose,"%d %e %e\n", nn,
                   gdl.histNNSubXi2pcf[nn], gdl.histNNSubXi2pcftotal[nn]);
            gdl.histXi2pcf[nn] /= MAX(gdl.histNNSubXi2pcftotal[nn],1.0);
#else
            if (cmd->verbose>3)
                verb_print_debug(cmd->verbose,"%d %e\n", nn,
                                 gdl.histNNSubXi2pcf[nn]);
            if (cballs_opt_weights_norm(cmd))
                gdl.histXi2pcf[nn] /= gdl.histWW[nn];
            else
                gdl.histXi2pcf[nn] /= MAX(gdl.histNNSubXi2pcf[nn],1.0);
#endif
//E
        }
    }
    //E

    if (cballs_opt_compute_histn(cmd)) {
#ifdef SMOOTHPIVOT
            search_compute_HistN_ggg(cmd, gd, nbody[cat1]-ipfalse, &gdl);
#else
            search_compute_HistN_ggg(cmd, gd, nbody[cat1], &gdl);
#endif
    }
#endif // ! TWOPCF

#if defined(THREEPCFCONVERGENCE) && defined(NMultipoles) \
    && defined(NONORMHIST)
    if (run_threepcf
        && cballs_opt_no_normalize_histzeta(cmd)
        && cballs_opt_edge_corrections(cmd)) {
        if (compute_edge_corrections_ggg(cmd, gd, &gdl, &gdlN)
            == FAILURE) {
            search_free_gd_sincos_omp_ggg_N(cmd, gd, &gdlN);
            search_free_gd_sincos_omp_ggg(cmd, gd, &gdl);
            return FAILURE;
        }
    }
#endif

// ===============================================
//B Saving histograms section: case GGGCORRELATION:
// ===============================================

      if (gd->rootDirFlag == TRUE) {
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
            "\n\t%s: printing %s method...\n\n",
            routineName, cmd->searchMethod);

#ifdef TWOPCF
    if (cballs_opt_compute_histn(cmd)) PRINT_OR_FAIL(PrintHistNN(cmd, gd, &gdl));
    PRINT_OR_FAIL(PrintHistXi2pcf(cmd, gd, &gdl));
#endif

    PRINT_OR_FAIL(PrintHistrBins(cmd, gd));

#ifdef THREEPCFCONVERGENCE
    if (run_threepcf) {
#ifdef NMultipoles

#ifdef NONORMHIST
        if (cballs_opt_no_normalize_histzeta(cmd))
            PRINT_OR_FAIL(PrintHistZetaM_sincos(cmd, gd, &gdl));
        else
            PRINT_OR_FAIL(PrintHistZetaM_sincos_normalized(cmd, gd, &gdl, &gdlN));
#else
        PRINT_OR_FAIL(PrintHistZetaM_sincos(cmd, gd, &gdl));
#endif // ! NONORMHIST

          PRINT_OR_FAIL(PrintHistZetaM_sincos_N(cmd, gd, &gdlN));

#else // ! NMultipoles
    PRINT_OR_FAIL(PrintHistZetaM_sincos(cmd, gd, &gdl));
#endif // ! NMultipoles

        if (cballs_opt_out_m_histzeta(cmd)) {
#ifdef NMultipoles

#ifdef NONORMHIST
            if (cballs_opt_no_normalize_histzeta(cmd))
                PRINT_OR_FAIL(PrintHistZetaMm_sincos(cmd, gd, &gdl));
            else
                PRINT_OR_FAIL(PrintHistZetaMm_sincos_normalized(cmd, gd, &gdl, &gdlN));
#else
            PRINT_OR_FAIL(PrintHistZetaMm_sincos(cmd, gd, &gdl));
#endif // ! NONORMHIST

            PRINT_OR_FAIL(PrintHistZetaMm_sincos_N(cmd, gd, &gdlN));

#else // ! NMultipoles
            PRINT_OR_FAIL(PrintHistZetaMm_sincos(cmd, gd, &gdl));
#endif // ! NMultipoles
        }

        if (cballs_opt_out_histzetag(cmd)) {
            PRINT_OR_FAIL(PrintHistZetaGm_sincos(cmd, gd, &gdl));
            PRINT_OR_FAIL(PrintHistZetaG(cmd, gd, &gdl));
            PRINT_OR_FAIL(PrintHistZetaMZetaGm_sincos(cmd, gd, &gdl));
        }

#ifdef NMultipoles
#ifdef NONORMHIST
    //  make this two questions more precise...
        if (cballs_opt_no_normalize_histzeta(cmd)) {
            if (cballs_opt_edge_corrections(cmd)) {
                PRINT_OR_FAIL(PrintHistZetaM_sincos_edge_effects(cmd, gd, &gdl, &gdlN));
            }
        }
#endif // ! NONORMHIST
#endif // ! NMultipoles
    }
#endif // ! THREEPCFCONVERGENCE

          gd->flagPrint = FALSE;
      }

// ===============================================
//E Saving histograms section: case GGGCORRELATION
// ===============================================

// ===============================================
//B Making histograms public (cballys PXD) section
// ===============================================
#ifdef PXD
    int m, n, n1, n2;
#ifdef TWOPCF
    for (n=1; n<=cmd->sizeHistN; n++) {
        gd->histNN[n] = gdl.histNN[n];
        gd->histCF[n] = gdl.histCF[n];
        gd->histXi2pcf[n] = gdl.histXi2pcf[n];
    }
#endif

//  gdl->histZetaMcos[m][n1][n2]/gdlN->histZetaMcos[1][n1][n2]);

#ifdef THREEPCFCONVERGENCE
    if (run_threepcf) {
    for (m=1; m<=cmd->mChebyshev+1; m++) {
        for (n1=1; n1<=cmd->sizeHistN; n1++) {
            for (n2=1; n2<=cmd->sizeHistN; n2++) {
                //B Working OK in cyballs
                if (cballs_opt_no_normalize_histzeta(cmd)) {
                    gd->histZetaMcos[m][n1][n2] = gdl.histZetaMcos[m][n1][n2];
                    gd->histZetaMsin[m][n1][n2] = gdl.histZetaMsin[m][n1][n2];
                    gd->histZetaMsincos[m][n1][n2] = gdl.histZetaMsincos[m][n1][n2];
                    gd->histZetaMcossin[m][n1][n2] = gdl.histZetaMcossin[m][n1][n2];
                } else {
#if defined(NMultipoles) && defined(NONORMHIST)
                    gd->histZetaMcos[m][n1][n2] =
                    normalize_zeta_ggg(gdl.histZetaMcos[m][n1][n2],
                                       &gdlN, n1, n2);
                    gd->histZetaMsin[m][n1][n2] =
                    normalize_zeta_ggg(gdl.histZetaMsin[m][n1][n2],
                                       &gdlN, n1, n2);
                    gd->histZetaMsincos[m][n1][n2] =
                    normalize_zeta_ggg(gdl.histZetaMsincos[m][n1][n2],
                                       &gdlN, n1, n2);
                    gd->histZetaMcossin[m][n1][n2] =
                    normalize_zeta_ggg(gdl.histZetaMcossin[m][n1][n2],
                                       &gdlN, n1, n2);
#else
                    gd->histZetaMcos[m][n1][n2] = gdl.histZetaMcos[m][n1][n2];
                    gd->histZetaMsin[m][n1][n2] = gdl.histZetaMsin[m][n1][n2];
                    gd->histZetaMsincos[m][n1][n2] = gdl.histZetaMsincos[m][n1][n2];
                    gd->histZetaMcossin[m][n1][n2] = gdl.histZetaMcossin[m][n1][n2];
#endif
                }
                // (EE) edge_effects
                if (cballs_opt_no_normalize_histzeta(cmd))
                    if (cballs_opt_edge_corrections(cmd))
                {
                    gd->histZetaM_EE[m][n1][n2] =
                        gdl.histZetaM_EE[m][n1][n2];
                    gd->histZetaM_EE_Im[m][n1][n2] =
                        gdl.histZetaM_EE_Im[m][n1][n2];
                }
            }
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
    if (run_threepcf)
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

global int searchcalc_octree_ggg_omp(struct cmdline_data *cmd,
                                     struct global_data *gd,
                                     bodyptr *btable, INTEGER *nbody,
                                     INTEGER ipmin, INTEGER *ipmax,
                                     int cat1, int cat2)
{
    return searchcalc_octree_ggg_driver(cmd, gd, btable, nbody,
                                        ipmin, ipmax, cat1, cat2, FALSE);
}

#ifdef OCTREEGGGMPI
global int searchcalc_octree_ggg_mpi(struct cmdline_data *cmd,
                                     struct global_data *gd,
                                     bodyptr *btable, INTEGER *nbody,
                                     INTEGER ipmin, INTEGER *ipmax,
                                     int cat1, int cat2)
{
    return searchcalc_octree_ggg_driver(cmd, gd, btable, nbody,
                                        ipmin, ipmax, cat1, cat2, TRUE);
}
#endif

local void normal_walktree_sincos(struct  cmdline_data* cmd, 
                                  struct  global_data* gd,
                                  bodyptr *btable, int cat2,
                                  bodyptr p, nodeptr q, real qsize,
                                  gdhistptr_sincos_omp_ggg hist)
{
    nodeptr l;
    real dr1;
    compute_vector dr;

    if (Update(p)==FALSE) return;
    if ( ((nodeptr) p) != q ) {
        if (Type(q) == CELL) {
            if (cballs_opt_read_mask(cmd)
                && Mask(q) == MASK_NODE_MASKED)
                return;
            if (!reject_cell_ggg(cmd, gd, (nodeptr)p, q)) {
                if (!cballs_opt_no_one_ball(cmd)
                    && (!cballs_opt_read_mask(cmd)
                        || mask_node_can_approximate(Mask(q)))) {
                    accept_body(cmd, gd, p, (nodeptr)q, &dr1, dr);

#ifndef NORMALHISTSCALE
//B useLogHist section
                    if ( (Radius(p)+Radius(q))/(dr1) < gd->deltaR)
                        sumnode_sincos_cell(cmd, gd, btable, cat2, p,
                                            ((cellptr) q), ((cellptr) q+1), hist);
                    else
                        for (l = More(q); l != Next(q); l = Next(l))
                            normal_walktree_sincos(cmd, gd, btable, cat2,
                                                   p,l,qsize/2, hist);
//E
#else // ! NORMALHISTSCALE
                    if ( (Radius(p)+Radius(q)) < gd->deltaR*THETA)
                        sumnode_sincos_cell(cmd, gd, btable, cat2, p,
                                            ((cellptr) q), ((cellptr) q+1), hist);
                    else
                        for (l = More(q); l != Next(q); l = Next(l))
                            normal_walktree_sincos(cmd, gd, btable, cat2,
                                                   p,l,qsize/2, hist);
#endif // ! NORMALHISTSCALE

                } else { // ! no-one-ball
                    for (l = More(q); l != Next(q); l = Next(l))
                        normal_walktree_sincos(cmd, gd, btable, cat2,
                                               p,l,qsize/2, hist);
                } // ! no-one-ball
            } // ! !reject_cell
        } else { // ! Type(q) == CELL
            sumnode_sincos(cmd, gd, btable, cat2,
                           p, ((cellptr)q), ((cellptr)q+1), hist);
        } // ! Type(q) == CELL
    } // ! p != q
}

local void sumnode_sincos(struct  cmdline_data* cmd,
                          struct  global_data* gd,
                          bodyptr *btable, int cat2,
                          bodyptr p, cellptr start, cellptr finish,
                          gdhistptr_sincos_omp_ggg hist)
{
    cellptr q;
    real dr1;
    compute_vector dr;
    int n;
    real xi;
    REAL cosphi,sinphi;
    int iq;
#ifdef TPCFSHEAR
    real gamma1, gamma2;
#endif

    q = start;
    if (cballs_opt_read_mask(cmd))
        if (Mask(q) != MASK_NODE_VALID) return;

    if (accept_body(cmd, gd, p, (nodeptr)q, &dr1, dr)) {
#ifndef NORMALHISTSCALE
//B useLogHist section
        if(dr1>cmd->rminHist) {
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
#ifdef SMOOTHPIVOT
                hist->histNNSubXi2pcfthreadp[n] =
                                            hist->histNNSubXi2pcfthreadp[n] + 1.;
#endif
#endif

                hist->histNNSubthread[n] = hist->histNNSubthread[n] + 1.;

                xi = Weight(q)*Kappa(q);

                if (hist->threepcf_enabled) {
#ifdef TPCFSHEAR
                gamma1 = Gamma1(q);
                gamma2 = Gamma2(q);
#endif
#if defined(THREEPCFCONVERGENCE) || defined(TPCFSHEAR)
                POLARAXIS_MAIN;
#endif
#ifdef THREEPCFCONVERGENCE
                if (cmd->mChebyshev<7) {
                    CHEBYSHEVTUOMPSINCOSANY;
                } else {
                    CHEBYSHEVTUOMP;
                }
#endif

#ifdef TPCFSHEAR
                CHEBYSHEVTUOMPGGGANY;
#endif
                }
#ifdef TWOPCF
                hist->histXi2pcfthreadsub[n] += xi;
#endif
                hist->nbbcalcthread += 1;
            } // ! 1 < n < sizeHistN
        } // ! dr1>cmd->rminHist
//E useLogHist section
#else // ! NORMALHISTSCALE

        if(dr1>cmd->rminHist) {
            n = (int) ( (dr1-cmd->rminHist) * gd->i_deltaR) + 1;
            if (n<=cmd->sizeHistN && n>=1) {

#ifdef TWOPCF
                hist->histNthread[n] = hist->histNthread[n] + 1.;
                hist->histWthread[n] = hist->histWthread[n] + 1.;
                hist->histNNSubXi2pcfthread[n] =
                                            hist->histNNSubXi2pcfthread[n] + 1.;
#ifdef SMOOTHPIVOT
                hist->histNNSubXi2pcfthreadp[n] =
                                            hist->histNNSubXi2pcfthreadp[n] + 1.;
#endif
#endif

                hist->histNNSubthread[n] = hist->histNNSubthread[n] + 1.;

                xi = Weight(q)*Kappa(q);

                if (hist->threepcf_enabled) {
#ifdef TPCFSHEAR
                gamma1 = Gamma1(q);
                gamma2 = Gamma2(q);
#endif
#if defined(THREEPCFCONVERGENCE) || defined(TPCFSHEAR)
                POLARAXIS_MAIN;
#endif
#ifdef THREEPCFCONVERGENCE
                if (cmd->mChebyshev<7) {
                    CHEBYSHEVTUOMPSINCOSANY;
                } else {
                    CHEBYSHEVTUOMP;
                }
#endif

#ifdef TPCFSHEAR
                CHEBYSHEVTUOMPGGGANY;
#endif
                }
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
                               gdhistptr_sincos_omp_ggg hist)
{
    cellptr q;
//B implement these memory optimization and check its results
    real dr1;
    compute_vector dr;
//E
    int n;
    real xi;
    REAL cosphi,sinphi;
#ifdef TPCFSHEAR
    real gamma1, gamma2;
#endif

    q = start;

    if (cballs_opt_read_mask(cmd))
        if (Mask(q) != MASK_NODE_VALID) return;

    if (accept_body(cmd, gd, p, (nodeptr)q, &dr1, dr)) {
        

#ifndef NORMALHISTSCALE
//B useLogHist section
        if(dr1>cmd->rminHist) {
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
                                            hist->histNNSubXi2pcfthread[n] + Nb(q);
#ifdef SMOOTHPIVOT
                hist->histNNSubXi2pcfthreadp[n] =
                                            hist->histNNSubXi2pcfthreadp[n] + Nb(q);
#endif
#endif

                hist->histNNSubthread[n] = hist->histNNSubthread[n] + Nb(q);

#ifndef NOWKAvg
                xi = Weight(q)*Kappa(q);
#else
                xi = Nb(q)*Kappa(q);
#endif
//E

                if (hist->threepcf_enabled) {
#ifdef TPCFSHEAR
                gamma1 = Gamma1(q);
                gamma2 = Gamma2(q);
#endif
#if defined(THREEPCFCONVERGENCE) || defined(TPCFSHEAR)
                POLARAXIS_MAIN;
#endif
#ifdef THREEPCFCONVERGENCE
                if (cmd->mChebyshev<7) {
                    CHEBYSHEVTUOMPSINCOSANY
                } else {
                    CHEBYSHEVTUOMP;
                }
#endif

#ifdef TPCFSHEAR
                CHEBYSHEVTUOMPGGGANY;
#endif
                }
#ifdef TWOPCF
                hist->histXi2pcfthreadsub[n] += xi;
#endif
                hist->nbccalcthread += 1;
            } // ! 1 < n < sizeHistN
        } // ! dr1 > rminHist
//E useLogHist section
#else // ! NORMALHISTSCALE

        if(dr1>cmd->rminHist) {
            n = (int) ( (dr1-cmd->rminHist) * gd->i_deltaR) + 1;
            if (n<=cmd->sizeHistN && n>=1) {

#ifdef TWOPCF
                hist->histNthread[n] = hist->histNthread[n] +  Nb(q);
                hist->histWthread[n] = hist->histWthread[n] +  Nb(q);
                hist->histNNSubXi2pcfthread[n] =
                                            hist->histNNSubXi2pcfthread[n] + Nb(q);
#ifdef SMOOTHPIVOT
                hist->histNNSubXi2pcfthreadp[n] =
                                            hist->histNNSubXi2pcfthreadp[n] + Nb(q);
#endif
#endif

                hist->histNNSubthread[n] = hist->histNNSubthread[n] + Nb(q);

#ifndef NOWKAvg
                //B ... this definition will give tests OK
                xi = Weight(q)*Kappa(q);
                //E
#else
                xi = Nb(q)*Kappa(q);
#endif

                if (hist->threepcf_enabled) {
#ifdef TPCFSHEAR
                gamma1 = Gamma1(q);
                gamma2 = Gamma2(q);
#endif
#if defined(THREEPCFCONVERGENCE) || defined(TPCFSHEAR)
                POLARAXIS_MAIN;
#endif
#ifdef THREEPCFCONVERGENCE
                if (cmd->mChebyshev<7) {
                    CHEBYSHEVTUOMPSINCOSANY
                } else {
                    CHEBYSHEVTUOMP;
                }
#endif

#ifdef TPCFSHEAR
                CHEBYSHEVTUOMPGGGANY;
#endif
                }
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

//B Weight: check Weight factor... must be an average of Weights
//==============
#ifndef NOWKAvg     // this definition will give tests OK
//==============

#ifdef NONORMHIST

#ifdef BALLS4SCANLEV
    xi = Weight(p)*Kappa(p);                    // equiv to Nb*(Weight/Nb)*Kappa
#ifdef SMOOTHPIVOT
        xi = Nb(p)*KappaRmin(p)/NbRmin(p);
#endif
#else
    xi = Weight(p)*Kappa(p);
#ifdef SMOOTHPIVOT
        xi = KappaRmin(p)/NbRmin(p);
#endif
#endif

#ifdef TWOPCF
#ifdef BALLS4SCANLEV
    wi = Weight(p);
    xi_2p = (Weight(p)/Nb(p))*Kappa(p);
#ifdef SMOOTHPIVOT
        xi_2p = KappaRmin(p);
#endif
#else
    wi = Weight(p);
    xi_2p = Weight(p)*Kappa(p);
#ifdef SMOOTHPIVOT
        xi_2p = KappaRmin(p);
#endif
#endif // ! BALLS4SCANLEV
#endif

#else // ! NONORMHIST

#ifdef BALLS4SCANLEV
    xi = (Weight(p)/Nb(p))*Kappa(p)/nbody;    // equiv to Nb*(Weight/Nb)*Kappa
#ifdef SMOOTHPIVOT                            // not working inside BALLS4SCANLEV
    xi = KappaRmin(p)/nbody;
#endif
#else // ! BALLS4SCANLEV
    xi = Weight(p)*Kappa(p)/nbody;
#ifdef SMOOTHPIVOT
    //B works if added a line for WeightRmin at the begining of
    //      sumnode_sincos above (first SMOOTHPIVOT segment)
    xi = WeightRmin(p)*KappaRmin(p)/nbody;
    //E
#endif
    //E
#endif // ! BALLS4SCANLEV

#ifdef TWOPCF
#ifdef BALLS4SCANLEV
    wi = Weight(p);
    xi_2p = (Weight(p)/Nb(p))*Kappa(p);
#ifdef SMOOTHPIVOT                            // not working inside BALLS4SCANLEV
        xi_2p = KappaRmin(p);
#endif
#else
    wi = Weight(p);
    xi_2p = Weight(p)*Kappa(p);
#ifdef SMOOTHPIVOT
        xi_2p = (KappaRmin(p));
#endif
#endif
#endif

#endif // ! NONORMHIST

//==============
#else // ! NOWKAvg
//==============

//B Weight:: not corrected yet this segment... use Weight(p) or not...
#ifdef NONORMHIST

    xi = Weight(p)*Kappa(p);                    // equiv to Nb*(Weight/Nb)*Kappa
#ifdef BALLS4SCANLEV
#ifdef SMOOTHPIVOT
        xi = Nb(p)*Weight(p)*KappaRmin(p)/NbRmin(p)/NbRmin(p);
#endif
#else // ! BALLS4SCANLEV
#ifdef SMOOTHPIVOT
    // check this line thoroughly... Weight must be smoothed also...
    xi = (WeightRmin(p)/NbRmin(p))*KappaRmin(p)/NbRmin(p);
#endif
#endif // ! BALLS4SCANLEV

#ifdef TWOPCF
    wi = Weight(p);
#ifdef BALLS4SCANLEV
    // check this line thoroughly against treeload.c
    xi_2p = (Weight(p)/Nb(p))*(Kappa(p)/Nb(p));
#ifdef SMOOTHPIVOT
        xi_2p = Nb(p)*Weight(p)*KappaRmin(p)/NbRmin(p)/NbRmin(p);
#endif
#else // ! BALLS4SCANLEV
    xi_2p = Weight(p)*Kappa(p);
#ifdef SMOOTHPIVOT
        xi_2p = (WeightRmin(p)/NbRmin(p))*KappaRmin(p)/NbRmin(p);
#endif
#endif // ! BALLS4SCANLEV
#endif // ! TWOPCF

#else // ! NONORMHIST

#ifdef BALLS4SCANLEV
    xi = (Weight(p)/Nb(p))*Kappa(p)/nbody;      // equiv to Nb*(Weight/Nb)*Kappa
#ifdef SMOOTHPIVOT
        xi = KappaRmin(p)/NbRmin(p)/nbody;
#endif
#else // ! BALLS4SCANLEV
    xi = Weight(p)*Kappa(p)/nbody;
#ifdef SMOOTHPIVOT
        xi = (WeightRmin(p)/NbRmin(p))*(KappaRmin(p)/NbRmin(p))/nbody;
#endif
#endif // ! BALLS4SCANLEV

#ifdef TWOPCF
    wi = Weight(p);
#ifdef BALLS4SCANLEV
    //B check if Weight and Kappa are averaged or not in treeload.c
    xi_2p = (Weight(p)/Nb(p))*Kappa(p)/Nb(p);
#ifdef SMOOTHPIVOT
        xi_2p = WeightRmin(p)*KappaRmin(p)/NbRmin(p)/NbRmin(p);
#endif
#else // ! BALLS4SCANLEV
    xi_2p = Weight(p)*Kappa(p);
#ifdef SMOOTHPIVOT
    //B divide by Nb(p)
        xi_2p = WeightRmin(p)*KappaRmin(p)/NbRmin(p)/NbRmin(p);
#endif
#endif // ! BALLS4SCANLEV
#endif // ! TWOPCF

#endif // ! NONORMHIST
//E Weight:: not corrected yet this segment... use Weight(p) or not...
#endif // ! NOWKAvg
//E

    if (hist->threepcf_enabled) {
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
        for (n=1; n<=cmd->sizeHistN; n++) {
            hist->xiOUTVPcos[n][n] -= hist->histXithreaddiagcos[m][n];
            hist->xiOUTVPsin[n][n] -= hist->histXithreaddiagsin[m][n];
            hist->xiOUTVPsincos[n][n] -=
                hist->histXithreaddiagsincos[m][n];
            hist->xiOUTVPcossin[n][n] -=
                hist->histXithreaddiagsincos[m][n];
        }
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
#endif

#ifdef TPCFSHEAR
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
    }
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
                                    gdhistptr_sincos_omp_ggg_N histN)
{
    nodeptr l;
    real dr1;
    compute_vector dr;

    if (Update(p)==FALSE) return;
    if ( ((nodeptr) p) != q ) {
        if (Type(q) == CELL) {
            if (cballs_opt_read_mask(cmd)
                && Mask(q) == MASK_NODE_MASKED)
                return;
            if (!reject_cell_ggg(cmd, gd, (nodeptr)p, q)) {
                if (!cballs_opt_no_one_ball(cmd)
                    && (!cballs_opt_read_mask(cmd)
                        || mask_node_can_approximate(Mask(q)))) {
                    accept_body(cmd, gd, p, (nodeptr)q, &dr1, dr);

#ifndef NORMALHISTSCALE
//B useLogHist section
                    if ( (Radius(p)+Radius(q))/(dr1) < gd->deltaR)
                        sumnode_sincos_cell_N(cmd, gd, btable, cat2,
                                              p, ((cellptr) q),
                                              ((cellptr) q+1), hist, histN);
                    else
                        for (l = More(q); l != Next(q); l = Next(l))
                            normal_walktree_sincos_N(cmd, gd, btable, cat2,
                                                     p,l,qsize/2,
                                                     hist, histN);
//E
#else // ! NORMALHISTSCALE
                    if ( (Radius(p)+Radius(q)) < gd->deltaR*THETA)
                        sumnode_sincos_cell_N(cmd, gd, btable, cat2,
                                              p, ((cellptr) q),
                                              ((cellptr) q+1), hist, histN);
                    else
                        for (l = More(q); l != Next(q); l = Next(l))
                            normal_walktree_sincos_N(cmd, gd, btable, cat2,
                                                     p,l,qsize/2,
                                                     hist, histN);
#endif // ! NORMALHISTSCALE

                } else { // ! no-one-ball
                    for (l = More(q); l != Next(q); l = Next(l))
                        normal_walktree_sincos_N(cmd, gd, btable, cat2,
                                                 p,l,qsize/2,
                                                 hist, histN);
                } // ! no-one-ball
            }
        } else { // ! Type(q) == CELL
            sumnode_sincos_N(cmd, gd, btable, cat2,
                             p, ((cellptr)q), ((cellptr)q+1),
                             hist, histN);
            //nbList, intList);
        } // ! Type(q) == CELL
    } // ! p != q
}

local void sumnode_sincos_N(struct  cmdline_data* cmd,
                            struct  global_data* gd,
                            bodyptr *btable, int cat2,
                            bodyptr p, cellptr start, cellptr finish,
                            gdhistptr_sincos_omp_ggg hist,
                            gdhistptr_sincos_omp_ggg_N histN)
{
    cellptr q;
    real dr1;
    compute_vector dr;
    int n;
    real xi;
    real xiN;
    REAL cosphi,sinphi;
    int iq;
#ifdef TPCFSHEAR
    real gamma1, gamma2;
#endif

    q = start;
    if (cballs_opt_read_mask(cmd))
        if (Mask(q) != MASK_NODE_VALID) return;

    if (accept_body(cmd, gd, p, (nodeptr)q, &dr1, dr)) {
#ifndef NORMALHISTSCALE
//B useLogHist section
        if(dr1>cmd->rminHist) {
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
#ifdef SMOOTHPIVOT
                hist->histNNSubXi2pcfthreadp[n] =
                                            hist->histNNSubXi2pcfthreadp[n] + 1.;
#endif
#endif

                hist->histNNSubthread[n] = hist->histNNSubthread[n] + 1.;
                histN->histNNSubthread[n] = histN->histNNSubthread[n] + 1.;

                xi = Weight(q)*Kappa(q);
                xiN = Weight(q);

#ifdef TPCFSHEAR
                gamma1 = Gamma1(q);
                gamma2 = Gamma2(q);
#endif
#if defined(THREEPCFCONVERGENCE) || defined(TPCFSHEAR)
                POLARAXIS_MAIN;
#endif
#ifdef THREEPCFCONVERGENCE
                GGG_ACCUMULATE_WINDOW_MODES;
                CHEBYSHEVTUOMPSINCOSANY;
                GGG_ACCUMULATE_NUMERATOR_DIAGONAL;
#endif

#ifdef TPCFSHEAR
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
            n = (int) ( (dr1-cmd->rminHist) * gd->i_deltaR) + 1;
            if (n<=cmd->sizeHistN && n>=1) {

#ifdef TWOPCF
                hist->histNthread[n] = hist->histNthread[n] + 1.;
                hist->histWthread[n] = hist->histWthread[n] + 1.;
                hist->histNNSubXi2pcfthread[n] =
                                            hist->histNNSubXi2pcfthread[n] + 1.;
#ifdef SMOOTHPIVOT
                hist->histNNSubXi2pcfthreadp[n] =
                                            hist->histNNSubXi2pcfthreadp[n] + 1.;
#endif
#endif

                hist->histNNSubthread[n] = hist->histNNSubthread[n] + 1.;
                histN->histNNSubthread[n] = histN->histNNSubthread[n] + 1.;

                xi = Weight(q)*Kappa(q);
                xiN = Weight(q);

#ifdef TPCFSHEAR
                gamma1 = Gamma1(q);
                gamma2 = Gamma2(q);
#endif
#if defined(THREEPCFCONVERGENCE) || defined(TPCFSHEAR)
                POLARAXIS_MAIN;
#endif
#ifdef THREEPCFCONVERGENCE
                GGG_ACCUMULATE_WINDOW_MODES;
                CHEBYSHEVTUOMPSINCOSANY;
                GGG_ACCUMULATE_NUMERATOR_DIAGONAL;
#endif

#ifdef TPCFSHEAR
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
                                 gdhistptr_sincos_omp_ggg_N histN)
{
    cellptr q;
    real dr1;
    compute_vector dr;
    int n;
    real xi;
    real xiN;
    REAL cosphi,sinphi;
#ifdef TPCFSHEAR
    real gamma1, gamma2;
#endif

    q = start;
    if (cballs_opt_read_mask(cmd))
        if (Mask(q) != MASK_NODE_VALID) return;

    if (accept_body(cmd, gd, p, (nodeptr)q, &dr1, dr)) {

#ifndef NORMALHISTSCALE
//B useLogHist section
        if(dr1>cmd->rminHist) {
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
                                            hist->histNNSubXi2pcfthread[n] + Nb(q);
#ifdef SMOOTHPIVOT
                hist->histNNSubXi2pcfthreadp[n] =
                                        hist->histNNSubXi2pcfthreadp[n] + Nb(q);
#endif
#endif

                //B 3pcf convergence
                hist->histNNSubthread[n] = hist->histNNSubthread[n] + Nb(q);
                histN->histNNSubthread[n] = histN->histNNSubthread[n] + Nb(q);
                //E

#ifndef NOWKAvg
                xi = Weight(q)*Kappa(q);
                xiN = Weight(q);
#else
                xi = Nb(q)*Kappa(q);
                xiN = Nb(q)*Weight(q);
#endif

#ifdef TPCFSHEAR
                gamma1 = Gamma1(q);
                gamma2 = Gamma2(q);
#endif
#if defined(THREEPCFCONVERGENCE) || defined(TPCFSHEAR)
                POLARAXIS_MAIN;
#endif
#ifdef THREEPCFCONVERGENCE
                GGG_ACCUMULATE_WINDOW_MODES;
                CHEBYSHEVTUOMPSINCOSANY;
                GGG_ACCUMULATE_NUMERATOR_DIAGONAL;
#endif

#ifdef TPCFSHEAR
                CHEBYSHEVTUOMPGGGANY;
#endif
#ifdef TWOPCF
                hist->histXi2pcfthreadsub[n] += xi;
#endif

                hist->nbccalcthread += 1;
            } // ! n <= sizeHistN && n >= 1
        } // ! dr1 > rminHist
//E useLogHist section

#else // ! NORMALHISTSCALE

        if(dr1>cmd->rminHist) {
            n = (int) ( (dr1-cmd->rminHist) * gd->i_deltaR) + 1;
            if (n<=cmd->sizeHistN && n>=1) {

#ifdef TWOPCF
                hist->histNthread[n] = hist->histNthread[n] +  Nb(q);
                hist->histWthread[n] = hist->histWthread[n] +  Nb(q);
                hist->histNNSubXi2pcfthread[n] =
                                            hist->histNNSubXi2pcfthread[n] + Nb(q);
#ifdef SMOOTHPIVOT
                hist->histNNSubXi2pcfthreadp[n] =
                                        hist->histNNSubXi2pcfthreadp[n] + Nb(q);
#endif
#endif

                //B 3pcf convergence
                hist->histNNSubthread[n] = hist->histNNSubthread[n] + Nb(q);
                histN->histNNSubthread[n] = histN->histNNSubthread[n] + Nb(q);
                //E

#ifndef NOWKAvg
                xi = Weight(q)*Kappa(q);
                xiN = Weight(q);
#else
                xi = Nb(q)*Kappa(q);
                xiN = Nb(q)*Weight(q);
#endif

#ifdef TPCFSHEAR
                gamma1 = Gamma1(q);
                gamma2 = Gamma2(q);
#endif
#if defined(THREEPCFCONVERGENCE) || defined(TPCFSHEAR)
                POLARAXIS_MAIN;
#endif
#ifdef THREEPCFCONVERGENCE
                GGG_ACCUMULATE_WINDOW_MODES;
                CHEBYSHEVTUOMPSINCOSANY;
                GGG_ACCUMULATE_NUMERATOR_DIAGONAL;
#endif

#ifdef TPCFSHEAR
                CHEBYSHEVTUOMPGGGANY;
#endif
#ifdef TWOPCF
                hist->histXi2pcfthreadsub[n] += xi;
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
    for (m=1; m<=GGG_WINDOW_ORDERS(cmd); m++) {
        for (n=1; n<=cmd->sizeHistN; n++) {
            hist->histXithreadcos[m][n] /= MAX(hist->histNNSubthread[n],1.0);
            hist->histXithreadsin[m][n] /= MAX(hist->histNNSubthread[n],1.0);
        }
    }
#endif


#ifdef THREEPCFCONVERGENCE
    for (m=1; m<=GGG_WINDOW_ORDERS(cmd); m++) {
        OUTVP_ext(hist->xiOUTVPcos,
            hist->histXithreadcos[m], hist->histXithreadcos[m], cmd->sizeHistN);
        OUTVP_ext(hist->xiOUTVPsin,
            hist->histXithreadsin[m], hist->histXithreadsin[m],cmd->sizeHistN);
        OUTVP_ext(hist->xiOUTVPsincos,
            hist->histXithreadsin[m], hist->histXithreadcos[m],cmd->sizeHistN);
        OUTVP_ext(hist->xiOUTVPcossin,
            hist->histXithreadcos[m], hist->histXithreadsin[m],cmd->sizeHistN);
        for (n=1; n<=cmd->sizeHistN; n++) {
            hist->xiOUTVPcos[n][n] -= hist->histXithreaddiagcos[m][n];
            hist->xiOUTVPsin[n][n] -= hist->histXithreaddiagsin[m][n];
            hist->xiOUTVPsincos[n][n] -=
                hist->histXithreaddiagsincos[m][n];
            hist->xiOUTVPcossin[n][n] -=
                hist->histXithreaddiagsincos[m][n];
        }
        CLRM_ext(hist->histZetaMtmpcos, cmd->sizeHistN);
        CLRM_ext(hist->histZetaMtmpsin, cmd->sizeHistN);
        CLRM_ext(hist->histZetaMtmpsincos, cmd->sizeHistN);
        CLRM_ext(hist->histZetaMtmpcossin, cmd->sizeHistN);
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
#endif

    return SUCCESS;
}

#endif // ! NMultipoles


//B Routines as in cballsutils

typedef struct {
    struct cmdline_data *cmd;
    struct global_data *gd;
    gdlptr_sincos_omp_ggg gdl;
} ggg_global_hist_init_context;

local int search_init_gd_sincos_omp_ggg_callback(void *argument)
{
    ggg_global_hist_init_context *context = argument;
    return search_init_gd_sincos_omp_ggg_unguarded(
        context->cmd, context->gd, context->gdl);
}

local int search_init_gd_sincos_omp_ggg(struct cmdline_data *cmd,
                                        struct global_data *gd,
                                        gdlptr_sincos_omp_ggg gdl)
{
    ggg_global_hist_init_context context = {cmd, gd, gdl};

    memset(gdl, 0, sizeof(*gdl));
    if (cballs_allocation_guard(search_init_gd_sincos_omp_ggg_callback,
                                &context, cmd->error_message,
                                sizeof(cmd->error_message)) == FAILURE) {
        search_free_gd_sincos_omp_ggg(cmd, gd, gdl);
        return FAILURE;
    }
    return SUCCESS;
}

local int search_init_gd_sincos_omp_ggg_unguarded(struct cmdline_data *cmd,
                                        struct global_data *gd,
                                        gdlptr_sincos_omp_ggg gdl)
{
    int n;
    int m;
    INTEGER bytes_tot_local=0;

    gdl->threepcf_enabled = !cballs_opt_only_2pcf(cmd);

#ifdef TWOPCF
    gdl->histNN = dvector(1,cmd->sizeHistN);
    gdl->histWW = dvector(1,cmd->sizeHistN);
    gdl->histCF = dvector(1,cmd->sizeHistN);
    gdl->histNNSubXi2pcf = dvector(1,cmd->sizeHistN);
#ifdef SMOOTHPIVOT
    gdl->histNNSubXi2pcftotal = dvector(1,cmd->sizeHistN);
#endif
    gdl->histXi2pcf = dvector(1,cmd->sizeHistN);
    bytes_tot_local += 5*cmd->sizeHistN*sizeof(real);
#endif

    gdl->histNNSub = dvector(1,cmd->sizeHistN);
    bytes_tot_local += 1*cmd->sizeHistN*sizeof(real);

#ifdef THREEPCFCONVERGENCE
    if (gdl->threepcf_enabled) {
    gdl->histZetaMcos = dmatrix3D(
            1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    gdl->histZetaMsin = dmatrix3D(
            1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    gdl->histZetaMsincos = dmatrix3D(
            1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
    gdl->histZetaMcossin = dmatrix3D(
            1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    gdl->histZetaM_EE = dmatrix3D(
            1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    gdl->histZetaM_EE_Im = dmatrix3D(
            1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    bytes_tot_local +=
            7*(cmd->mChebyshev+1)*cmd->sizeHistN*cmd->sizeHistN*sizeof(real);

    gdl->histZetaGmRe = dmatrix3D(
            1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    gdl->histZetaGmIm = dmatrix3D(
            1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    bytes_tot_local +=
            2*(cmd->mChebyshev+1)*cmd->sizeHistN*cmd->sizeHistN*sizeof(real);
    gdl->histXi3pcf = dmatrix3D(
            1,cmd->sizeHistPhi,1,cmd->sizeHistN,1,cmd->sizeHistN);
    bytes_tot_local +=
            cmd->sizeHistN*cmd->sizeHistN*cmd->sizeHistPhi*sizeof(real);
    }
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
#ifdef SMOOTHPIVOT
        gdl->histNNSubXi2pcftotal[n] = 0.0;
#endif
        gdl->histXi2pcf[n] = 0.0;
    }
#endif

    for (n = 1; n <= cmd->sizeHistN; n++)
        gdl->histNNSub[n] = 0.0;

#ifdef THREEPCFCONVERGENCE
    if (gdl->threepcf_enabled) {
    for (m = 1; m <= cmd->mChebyshev+1; m++) {
        CLRM_ext(gdl->histZetaMcos[m], cmd->sizeHistN);
        CLRM_ext(gdl->histZetaMsin[m], cmd->sizeHistN);
        CLRM_ext(gdl->histZetaMsincos[m], cmd->sizeHistN);
        CLRM_ext(gdl->histZetaMcossin[m], cmd->sizeHistN);
        // (EE) edge_effects
        CLRM_ext(gdl->histZetaM_EE[m], cmd->sizeHistN);
        CLRM_ext(gdl->histZetaM_EE_Im[m], cmd->sizeHistN);
    }
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
    if (gdl->threepcf_enabled) {
    free_dmatrix3D(gdl->histXi3pcf,1,cmd->sizeHistPhi,
                   1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix3D(gdl->histZetaGmIm,
                   1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix3D(gdl->histZetaGmRe,
                   1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);

    // (EE) edge_effects
    free_dmatrix3D(gdl->histZetaM_EE_Im,
                   1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix3D(gdl->histZetaM_EE,
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
    }
#endif

    free_dvector(gdl->histNNSub,1,cmd->sizeHistN);

#ifdef TWOPCF
    free_dvector(gdl->histXi2pcf,1,cmd->sizeHistN);
#ifdef SMOOTHPIVOT
    free_dvector(gdl->histNNSubXi2pcftotal,1,cmd->sizeHistN);
#endif
    free_dvector(gdl->histNNSubXi2pcf,1,cmd->sizeHistN);
    //
    free_dvector(gdl->histCF,1,cmd->sizeHistN);
    free_dvector(gdl->histWW,1,cmd->sizeHistN);
    free_dvector(gdl->histNN,1,cmd->sizeHistN);
#endif

    return SUCCESS;
}

typedef struct {
    struct cmdline_data *cmd;
    struct global_data *gd;
    gdhistptr_sincos_omp_ggg hist;
} ggg_hist_init_context;

local int search_init_sincos_omp_ggg_callback(void *argument)
{
    ggg_hist_init_context *context = argument;
    return search_init_sincos_omp_ggg_unguarded(
        context->cmd, context->gd, context->hist);
}

local int search_init_sincos_omp_ggg(struct cmdline_data *cmd,
                                     struct global_data *gd,
                                     gdhistptr_sincos_omp_ggg hist)
{
    ggg_hist_init_context context = {cmd, gd, hist};
    ErrorMsg allocation_error;

    memset(hist, 0, sizeof(*hist));
    if (cballs_allocation_guard(search_init_sincos_omp_ggg_callback,
                                &context, allocation_error,
                                sizeof(allocation_error)) == FAILURE) {
        search_free_sincos_omp_ggg(cmd, gd, hist);
        return FAILURE;
    }
    return SUCCESS;
}

local int search_init_sincos_omp_ggg_unguarded(struct cmdline_data *cmd,
                                     struct global_data *gd,
                                     gdhistptr_sincos_omp_ggg hist)
{
    int n;
    int m;
    INTEGER bytes_tot_local=0;

    hist->threepcf_enabled = !cballs_opt_only_2pcf(cmd);

#ifdef TWOPCF
    hist->histNthread = dvector(1,cmd->sizeHistN);
    hist->histWthread = dvector(1,cmd->sizeHistN);
    hist->histWWthread = dvector(1,cmd->sizeHistN);
    hist->histNNSubXi2pcfthread = dvector(1,cmd->sizeHistN);
    //B kappa Avg Rmin
#ifdef SMOOTHPIVOT
    hist->histNNSubXi2pcfthreadp = dvector(1,cmd->sizeHistN);
    hist->histNNSubXi2pcfthreadtotal = dvector(1,cmd->sizeHistN);
#endif
    //E
    hist->histXi2pcfthread = dvector(1,cmd->sizeHistN);
    hist->histXi2pcfthreadsub = dvector(1,cmd->sizeHistN);
    bytes_tot_local += 6*cmd->sizeHistN*sizeof(real);
#endif

    if (hist->threepcf_enabled) {
        hist->ChebsT = dvector(1,cmd->mChebyshev+1);
        hist->ChebsU = dvector(1,cmd->mChebyshev+1);
        bytes_tot_local += 2.0*(cmd->mChebyshev+1)*sizeof(real);
    }

    hist->histNNSubthread = dvector(1,cmd->sizeHistN);
    bytes_tot_local += cmd->sizeHistN*sizeof(real);

#ifdef THREEPCFCONVERGENCE
    if (hist->threepcf_enabled) {
    hist->histXithreadcos = dmatrix(1,cmd->mChebyshev+1,1,cmd->sizeHistN);
    hist->histXithreadsin = dmatrix(1,cmd->mChebyshev+1,1,cmd->sizeHistN);
    hist->histXithreaddiagcos =
            dmatrix(1,cmd->mChebyshev+1,1,cmd->sizeHistN);
    hist->histXithreaddiagsin =
            dmatrix(1,cmd->mChebyshev+1,1,cmd->sizeHistN);
    hist->histXithreaddiagsincos =
            dmatrix(1,cmd->mChebyshev+1,1,cmd->sizeHistN);
    bytes_tot_local += 5.0*(cmd->mChebyshev+1)*cmd->sizeHistN*sizeof(real);
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
    hist->histZetaMtmpcos = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
    hist->histZetaMtmpsin = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
    hist->histZetaMtmpsincos = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
    hist->histZetaMtmpcossin = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
    bytes_tot_local += 4.0*cmd->sizeHistN*cmd->sizeHistN*sizeof(real);
    }
#endif

#ifdef TPCFSHEAR
    if (hist->threepcf_enabled) {
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
    }
#endif

#ifdef TWOPCF
    for (n = 1; n <= cmd->sizeHistN; n++) {
        hist->histNthread[n] = 0.0;
        hist->histWthread[n] = 0.0;
        hist->histWWthread[n] = 0.0;
        hist->histNNSubXi2pcfthread[n] = 0.0;
#ifdef SMOOTHPIVOT
        hist->histNNSubXi2pcfthreadp[n] = 0.0;
        hist->histNNSubXi2pcfthreadtotal[n] = 0.0;
#endif
        hist->histXi2pcfthread[n] = 0.0;
        hist->histXi2pcfthreadsub[n] = 0.0;
    }
#endif

    for (n = 1; n <= cmd->sizeHistN; n++)
        hist->histNNSubthread[n] = 0.0;

#ifdef THREEPCFCONVERGENCE
    if (hist->threepcf_enabled) {
    for (m = 1; m <= cmd->mChebyshev+1; m++) {
        CLRM_ext(hist->histZetaMthreadcos[m], cmd->sizeHistN);
        CLRM_ext(hist->histZetaMthreadsin[m], cmd->sizeHistN);
        CLRM_ext(hist->histZetaMthreadsincos[m], cmd->sizeHistN);
        CLRM_ext(hist->histZetaMthreadcossin[m], cmd->sizeHistN);
        CLRV_ext(hist->histXithreaddiagcos[m], cmd->sizeHistN);
        CLRV_ext(hist->histXithreaddiagsin[m], cmd->sizeHistN);
        CLRV_ext(hist->histXithreaddiagsincos[m], cmd->sizeHistN);
    }
    }
#endif

#ifdef TPCFSHEAR
    if (hist->threepcf_enabled) {
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
#ifdef TPCFSHEAR
    if (hist->threepcf_enabled) {
    free_dmatrix(hist->histImGNthread,1,cmd->mChebyshev,1,cmd->sizeHistN);
    free_dmatrix(hist->histReGNthread,1,cmd->mChebyshev,1,cmd->sizeHistN);
    free_dmatrix(hist->histImGthread,1,cmd->mChebyshev+1,1,cmd->sizeHistN);
    free_dmatrix(hist->histReGthread,1,cmd->mChebyshev+1,1,cmd->sizeHistN);
    free_dmatrix(hist->histg2threadsin,1,cmd->mChebyshev+1,1,cmd->sizeHistN);
    free_dmatrix(hist->histg2threadcos,1,cmd->mChebyshev+1,1,cmd->sizeHistN);
    free_dmatrix(hist->histg1threadsin,1,cmd->mChebyshev+1,1,cmd->sizeHistN);
    free_dmatrix(hist->histg1threadcos,1,cmd->mChebyshev+1,1,cmd->sizeHistN);
    }
#endif

#ifdef THREEPCFCONVERGENCE
    if (hist->threepcf_enabled) {
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
    free_dmatrix(hist->histXithreaddiagsincos,
                 1,cmd->mChebyshev+1,1,cmd->sizeHistN);
    free_dmatrix(hist->histXithreaddiagsin,
                 1,cmd->mChebyshev+1,1,cmd->sizeHistN);
    free_dmatrix(hist->histXithreaddiagcos,
                 1,cmd->mChebyshev+1,1,cmd->sizeHistN);
    free_dmatrix(hist->histXithreadsin,1,cmd->mChebyshev+1,1,cmd->sizeHistN);
    free_dmatrix(hist->histXithreadcos,1,cmd->mChebyshev+1,1,cmd->sizeHistN);
    }
#endif

    free_dvector(hist->histNNSubthread,1,cmd->sizeHistN);

    if (hist->threepcf_enabled) {
        free_dvector(hist->ChebsU,1,cmd->mChebyshev+1);
        free_dvector(hist->ChebsT,1,cmd->mChebyshev+1);
    }

#ifdef TWOPCF
    free_dvector(hist->histXi2pcfthread,1,cmd->sizeHistN);
    free_dvector(hist->histXi2pcfthreadsub,1,cmd->sizeHistN);
#ifdef SMOOTHPIVOT
    free_dvector(hist->histNNSubXi2pcfthreadtotal,1,cmd->sizeHistN);
    free_dvector(hist->histNNSubXi2pcfthreadp,1,cmd->sizeHistN);
#endif
    free_dvector(hist->histNNSubXi2pcfthread,1,cmd->sizeHistN);
    //
    free_dvector(hist->histWWthread,1,cmd->sizeHistN);
    free_dvector(hist->histWthread,1,cmd->sizeHistN);
    free_dvector(hist->histNthread,1,cmd->sizeHistN);
#endif

    return SUCCESS;
}


#ifdef NMultipoles
typedef struct {
    struct cmdline_data *cmd;
    struct global_data *gd;
    gdlptr_sincos_omp_ggg_N gdl;
} ggg_global_hist_n_init_context;

local int search_init_gd_sincos_omp_ggg_N_callback(void *argument)
{
    ggg_global_hist_n_init_context *context = argument;
    return search_init_gd_sincos_omp_ggg_N_unguarded(
        context->cmd, context->gd, context->gdl);
}

local int search_init_gd_sincos_omp_ggg_N(struct cmdline_data *cmd,
                                        struct global_data *gd,
                                        gdlptr_sincos_omp_ggg_N gdl)
{
    ggg_global_hist_n_init_context context = {cmd, gd, gdl};

    memset(gdl, 0, sizeof(*gdl));
    if (cballs_allocation_guard(search_init_gd_sincos_omp_ggg_N_callback,
                                &context, cmd->error_message,
                                sizeof(cmd->error_message)) == FAILURE) {
        search_free_gd_sincos_omp_ggg_N(cmd, gd, gdl);
        return FAILURE;
    }
    return SUCCESS;
}

local int search_init_gd_sincos_omp_ggg_N_unguarded(struct cmdline_data *cmd,
                                        struct global_data *gd,
                                        gdlptr_sincos_omp_ggg_N gdl)
{
    int n;
    int m;
    int orders = GGG_WINDOW_ORDERS(cmd);
    INTEGER bytes_tot_local=0;

#ifdef TWOPCF
    gdl->histNN = dvector(1,cmd->sizeHistN);
    gdl->histWW = dvector(1,cmd->sizeHistN);
    gdl->histCF = dvector(1,cmd->sizeHistN);
    gdl->histNNSubXi2pcf = dvector(1,cmd->sizeHistN);
#ifdef SMOOTHPIVOT
    gdl->histNNSubXi2pcftotal = dvector(1,cmd->sizeHistN);
#endif
    gdl->histXi2pcf = dvector(1,cmd->sizeHistN);
    bytes_tot_local += 5*cmd->sizeHistN*sizeof(real);
#endif

    gdl->histNNSub = dvector(1,cmd->sizeHistN);
    bytes_tot_local += 1*cmd->sizeHistN*sizeof(real);

#ifdef THREEPCFCONVERGENCE
    gdl->histZetaMcos =
            dmatrix3D(1,orders,1,cmd->sizeHistN,1,cmd->sizeHistN);
    gdl->histZetaMsin =
            dmatrix3D(1,orders,1,cmd->sizeHistN,1,cmd->sizeHistN);
    gdl->histZetaMsincos =
            dmatrix3D(1,orders,1,cmd->sizeHistN,1,cmd->sizeHistN);
    // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
    gdl->histZetaMcossin =
            dmatrix3D(1,orders,1,cmd->sizeHistN,1,cmd->sizeHistN);
    bytes_tot_local +=
            4*orders*cmd->sizeHistN*cmd->sizeHistN*sizeof(real);
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
#ifdef SMOOTHPIVOT
        gdl->histNNSubXi2pcftotal[n] = 0.0;
#endif
        gdl->histXi2pcf[n] = 0.0;
    }
#endif

    for (n = 1; n <= cmd->sizeHistN; n++)
        gdl->histNNSub[n] = 0.0;

#ifdef THREEPCFCONVERGENCE
    for (m = 1; m <= orders; m++) {
        CLRM_ext(gdl->histZetaMcos[m], cmd->sizeHistN);
        CLRM_ext(gdl->histZetaMsin[m], cmd->sizeHistN);
        CLRM_ext(gdl->histZetaMsincos[m], cmd->sizeHistN);
        CLRM_ext(gdl->histZetaMcossin[m], cmd->sizeHistN);
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
    // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
    free_dmatrix3D(gdl->histZetaMcossin,
                   1,GGG_WINDOW_ORDERS(cmd),1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix3D(gdl->histZetaMsincos,
                   1,GGG_WINDOW_ORDERS(cmd),1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix3D(gdl->histZetaMsin,
                   1,GGG_WINDOW_ORDERS(cmd),1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix3D(gdl->histZetaMcos,
                   1,GGG_WINDOW_ORDERS(cmd),1,cmd->sizeHistN,1,cmd->sizeHistN);
#endif

    free_dvector(gdl->histNNSub,1,cmd->sizeHistN);

#ifdef TWOPCF
    free_dvector(gdl->histXi2pcf,1,cmd->sizeHistN);
#ifdef SMOOTHPIVOT
    free_dvector(gdl->histNNSubXi2pcftotal,1,cmd->sizeHistN);
#endif
    free_dvector(gdl->histNNSubXi2pcf,1,cmd->sizeHistN);
    //
    free_dvector(gdl->histCF,1,cmd->sizeHistN);
    free_dvector(gdl->histWW,1,cmd->sizeHistN);
    free_dvector(gdl->histNN,1,cmd->sizeHistN);
#endif

    return SUCCESS;
}

typedef struct {
    struct cmdline_data *cmd;
    struct global_data *gd;
    gdhistptr_sincos_omp_ggg_N hist;
} ggg_hist_n_init_context;

local int search_init_sincos_omp_ggg_N_callback(void *argument)
{
    ggg_hist_n_init_context *context = argument;
    return search_init_sincos_omp_ggg_N_unguarded(
        context->cmd, context->gd, context->hist);
}

local int search_init_sincos_omp_ggg_N(struct cmdline_data *cmd,
                                     struct global_data *gd,
                                     gdhistptr_sincos_omp_ggg_N hist)
{
    ggg_hist_n_init_context context = {cmd, gd, hist};
    ErrorMsg allocation_error;

    memset(hist, 0, sizeof(*hist));
    if (cballs_allocation_guard(search_init_sincos_omp_ggg_N_callback,
                                &context, allocation_error,
                                sizeof(allocation_error)) == FAILURE) {
        search_free_sincos_omp_ggg_N(cmd, gd, hist);
        return FAILURE;
    }
    return SUCCESS;
}

local int search_init_sincos_omp_ggg_N_unguarded(struct cmdline_data *cmd,
                                     struct global_data *gd,
                                     gdhistptr_sincos_omp_ggg_N hist)
{
    int n;
    int m;
    int orders = GGG_WINDOW_ORDERS(cmd);
    INTEGER bytes_tot_local=0;

#ifdef TWOPCF
    hist->histNthread = dvector(1,cmd->sizeHistN);
    hist->histWthread = dvector(1,cmd->sizeHistN);
    hist->histWWthread = dvector(1,cmd->sizeHistN);
    hist->histNNSubXi2pcfthread = dvector(1,cmd->sizeHistN);
#ifdef SMOOTHPIVOT
    hist->histNNSubXi2pcfthreadp = dvector(1,cmd->sizeHistN);
    hist->histNNSubXi2pcfthreadtotal = dvector(1,cmd->sizeHistN);
#endif
    hist->histXi2pcfthread = dvector(1,cmd->sizeHistN);
    hist->histXi2pcfthreadsub = dvector(1,cmd->sizeHistN);
    bytes_tot_local += 6*cmd->sizeHistN*sizeof(real);
#endif

    //B 3pcf convergence & shear
    hist->ChebsT = dvector(1,orders);
    hist->ChebsU = dvector(1,orders);
    //E 3pcf convergence & shear
    bytes_tot_local += 2.0*orders*sizeof(real);

    hist->histNNSubthread = dvector(1,cmd->sizeHistN);
    bytes_tot_local += cmd->sizeHistN*sizeof(real);

#ifdef THREEPCFCONVERGENCE
    hist->histXithreadcos = dmatrix(1,orders,1,cmd->sizeHistN);
    hist->histXithreadsin = dmatrix(1,orders,1,cmd->sizeHistN);
    hist->histXithreaddiagcos = dmatrix(1,orders,1,cmd->sizeHistN);
    hist->histXithreaddiagsin = dmatrix(1,orders,1,cmd->sizeHistN);
    hist->histXithreaddiagsincos = dmatrix(1,orders,1,cmd->sizeHistN);
    bytes_tot_local += 5.0*orders*cmd->sizeHistN*sizeof(real);
    hist->histZetaMthreadcos =
            dmatrix3D(1,orders,1,cmd->sizeHistN,1,cmd->sizeHistN);
    hist->histZetaMthreadsin =
            dmatrix3D(1,orders,1,cmd->sizeHistN,1,cmd->sizeHistN);
    hist->histZetaMthreadsincos =
            dmatrix3D(1,orders,1,cmd->sizeHistN,1,cmd->sizeHistN);
    hist->histZetaMthreadcossin =
            dmatrix3D(1,orders,1,cmd->sizeHistN,1,cmd->sizeHistN);
    bytes_tot_local +=
        4.0*orders*cmd->sizeHistN*cmd->sizeHistN*sizeof(real);
    hist->xiOUTVPcos = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
    hist->xiOUTVPsin = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
    hist->xiOUTVPsincos = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
    hist->xiOUTVPcossin = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
    hist->histZetaMtmpcos = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
    hist->histZetaMtmpsin = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
    hist->histZetaMtmpsincos = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
    hist->histZetaMtmpcossin = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
    bytes_tot_local += 8.0*cmd->sizeHistN*cmd->sizeHistN*sizeof(real);
#endif

#ifdef TPCFSHEAR
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
#ifdef SMOOTHPIVOT
        hist->histNNSubXi2pcfthreadp[n] = 0.0;
        hist->histNNSubXi2pcfthreadtotal[n] = 0.0;
#endif
        hist->histXi2pcfthread[n] = 0.0;
        hist->histXi2pcfthreadsub[n] = 0.0;
    }
#endif

    for (n = 1; n <= cmd->sizeHistN; n++)
        hist->histNNSubthread[n] = 0.0;

#ifdef THREEPCFCONVERGENCE
    for (m = 1; m <= orders; m++) {
        CLRM_ext(hist->histZetaMthreadcos[m], cmd->sizeHistN);
        CLRM_ext(hist->histZetaMthreadsin[m], cmd->sizeHistN);
        CLRM_ext(hist->histZetaMthreadsincos[m], cmd->sizeHistN);
        CLRM_ext(hist->histZetaMthreadcossin[m], cmd->sizeHistN);
        CLRV_ext(hist->histXithreaddiagcos[m], cmd->sizeHistN);
        CLRV_ext(hist->histXithreaddiagsin[m], cmd->sizeHistN);
        CLRV_ext(hist->histXithreaddiagsincos[m], cmd->sizeHistN);
    }
#endif

#ifdef TPCFSHEAR
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
#ifdef TPCFSHEAR
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
    free_dmatrix(hist->histZetaMtmpcossin,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix(hist->histZetaMtmpsincos,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix(hist->histZetaMtmpsin,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix(hist->histZetaMtmpcos,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix(hist->xiOUTVPcossin,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix(hist->xiOUTVPsincos,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix(hist->xiOUTVPsin,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix(hist->xiOUTVPcos,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix3D(hist->histZetaMthreadcossin,
                   1,GGG_WINDOW_ORDERS(cmd),1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix3D(hist->histZetaMthreadsincos,
                   1,GGG_WINDOW_ORDERS(cmd),1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix3D(hist->histZetaMthreadsin,
                   1,GGG_WINDOW_ORDERS(cmd),1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix3D(hist->histZetaMthreadcos,
                   1,GGG_WINDOW_ORDERS(cmd),1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix(hist->histXithreaddiagsincos,
                 1,GGG_WINDOW_ORDERS(cmd),1,cmd->sizeHistN);
    free_dmatrix(hist->histXithreaddiagsin,
                 1,GGG_WINDOW_ORDERS(cmd),1,cmd->sizeHistN);
    free_dmatrix(hist->histXithreaddiagcos,
                 1,GGG_WINDOW_ORDERS(cmd),1,cmd->sizeHistN);
    free_dmatrix(hist->histXithreadsin,
                 1,GGG_WINDOW_ORDERS(cmd),1,cmd->sizeHistN);
    free_dmatrix(hist->histXithreadcos,
                 1,GGG_WINDOW_ORDERS(cmd),1,cmd->sizeHistN);
#endif

    free_dvector(hist->histNNSubthread,1,cmd->sizeHistN);

    //B 3pcf convergence & shear
    free_dvector(hist->ChebsU,1,GGG_WINDOW_ORDERS(cmd));
    free_dvector(hist->ChebsT,1,GGG_WINDOW_ORDERS(cmd));
    //E 3pcf convergence & shear

#ifdef TWOPCF
    free_dvector(hist->histXi2pcfthread,1,cmd->sizeHistN);
    free_dvector(hist->histXi2pcfthreadsub,1,cmd->sizeHistN);
#ifdef SMOOTHPIVOT
    free_dvector(hist->histNNSubXi2pcfthreadtotal,1,cmd->sizeHistN);
    free_dvector(hist->histNNSubXi2pcfthreadp,1,cmd->sizeHistN);
#endif
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
    if ((cballs_opt_cute_box(cmd))) {
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
            if (cballs_opt_cute_box(cmd)) {
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
            if (cballs_opt_cute_box(cmd)) {
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

    if (cballs_opt_and_cf(cmd))
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
                        "searchcalc: Using %s... \n", cmd->searchMethod);

    if (cballs_opt_ggg_correlation(cmd))
        verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                            "computing using GGG routine... \n");


    if (cmd->usePeriodic==TRUE)
        verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
    "%s: warning!! can´t have periodic boundaries and OCTREEGGGOMP definition (usePeriodic=%d)\nSet usePeriodic=false\n",
                            routineName,cmd->usePeriodic);

    if (cmd->useLogHist==FALSE)
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
"%s: warning!! can´t have normal scale hist and OCTREEGGGOMP definition (useLogHist=%d)\nSet useLogHist=true\n",
                        routineName,cmd->useLogHist);
#ifdef TWOPCF
    verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                        "with 2pcf convergence computation... \n");
#endif
#ifdef THREEPCFCONVERGENCE
    verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                        "with 3pcf convergence computation... \n");
#endif
#ifdef THREEPCFCONVERGENCE
#ifndef TPCF
    snprintf(cmd->error_message, _ERRORMSGSIZE_,
             "\n%s: can´t have computeTPCF=false and THREEPCFCONVERGENCE definition\n",
                   routineName);
    return FAILURE;
#endif
#endif
#ifdef TPCFSHEAR
    verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                        "with 3pcf shear computation... \n");
#endif

#if NDIM == 2
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
    if (cballs_opt_no_normalize_histzeta(cmd))
        verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                            "with option no-normalize-HistZeta...\n");
#else
    verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                        "without NONORMHIST... \n");
#endif

#if defined(NMultipoles) && defined(NONORMHIST)
    if (cballs_opt_no_normalize_histzeta(cmd)) {
        if (cballs_opt_edge_corrections(cmd))
            verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                                "with option edge-corrections... \n");
    } else {
        if (cballs_opt_edge_corrections(cmd)) {
            verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                            "option edge-corrections only works with %s... \n",
                                "no-normalize-HistZeta option added");
            // Check freeing allocated memory...
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     "option edge-corrections only works with %s... \n",
                         "no-normalize-HistZeta option added");
            return FAILURE;
        }
    }
#else
    if (cballs_opt_edge_corrections(cmd)) {
        verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                    "option edge-corrections only works with %s activated... \n",
                    "NMultipoles && NONORMHIST");
        // Check freeing allocated memory...
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "option edge-corrections only works with %s activated... \n",
                 "NMultipoles && NONORMHIST");
        return FAILURE;
    }
#endif

#ifndef USEGSL
    if (cballs_opt_edge_corrections(cmd))
        verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
            "option edge-corrections is better computed with %s activated... \n",
            "USEGSL");
#endif

#ifdef POLARAXIS
    verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                        "with POLARAXIS... \n");
#endif
    if (cballs_opt_no_one_ball(cmd))
        verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                            "with option no-one-ball... \n");
#ifdef SMOOTHPIVOT
        verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                            "with option smooth-pivot... rsmooth=%g\n",
                            gd->rsmooth[0]);
#endif
    if (cballs_opt_default_rsmooth(cmd))
        verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                            "with option default-rsmooth... \n");
    if (cballs_opt_fix_rsmooth(cmd))
        verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                            "with option fix-rsmooth... \n");
#ifdef BALLS4SCANLEV
    verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                        "with BALLS4SCANLEV... \n");
#endif

#ifdef NORMALHISTSCALE
    verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                        "with NORMALHISTSCALE... \n");
#endif

#ifdef NOWKAvg
    verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                        "with NOWKAvg to correct WK product... \n");
#endif

    if (cballs_opt_read_mask(cmd))
        verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                            "with option read-mask... \n");

    if (cballs_opt_kappa_constant(cmd))
        verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                            "with option kappa-constant... \n");
    if (cballs_opt_celestial(cmd))
        verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                            "with option celestial... \n");
    if (cballs_opt_ra_reversed(cmd))
        verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                            "with option ra-reversed... \n");

    return SUCCESS;
}

//B Saving histograms section: case KKKCORRELATION:

local int PrintHistrBins(struct  cmdline_data* cmd, struct  global_data* gd)
{
    string routineName = "PrintHistrBins";
    real rBin, rbinlog;
    int n;
    stream outstr;

    OPEN_OUTPUT_OR_FAIL(outstr, gd->fpfnamehistrBinsFileName, "w!");

    verb_print_q(2, cmd->verbose,
               "Printing : to a file %s ...\n",gd->fpfnamehistrBinsFileName);

    for (n=1; n<=cmd->sizeHistN; n++) {
        if (cmd->useLogHist) {
            if (cmd->rminHist==0) {
                rbinlog =
                    ((real)(n-cmd->sizeHistN))/cmd->logHistBinsPD
                    + rlog10(cmd->rangeN);
            } else {
                rbinlog = rlog10(cmd->rminHist) + ((real)(n)-0.5)*gd->deltaR;
            }
            rBin=rpow(10.0,rbinlog);
        } else {
            rBin = cmd->rminHist + ((real)n-0.5)*gd->deltaR;
        }
        WRITE_OUTPUT_OR_FAIL(outstr, gd->fpfnamehistrBinsFileName,
                             "%16.8e\n",rBin);
    }
    CLOSE_OUTPUT_OR_FAIL(outstr, gd->fpfnamehistrBinsFileName);

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
    string routineName = "PrintHistZetaM_sincos";
    int n1, n2, m;
    stream outstr;
    char namebuf[256];

    for (m=1; m<=cmd->mChebyshev+1; m++) {
        if (format_checked(namebuf, sizeof(namebuf),
            "namebuf", "%s%s_%d%s", gd->fpfnamehistZetaMFileName,
                           "_cos", m, EXTFILES) != 0)
            return FAILURE;

        verb_print_q(2, cmd->verbose, "Printing : to a file %s ...\n",namebuf);
        OPEN_OUTPUT_OR_FAIL(outstr, namebuf, "w!");
        for (n1=1; n1<=cmd->sizeHistN; n1++) {
            for (n2=1; n2<=cmd->sizeHistN; n2++) {
                WRITE_OUTPUT_OR_FAIL(outstr, namebuf, "%16.8e ",
                                     gd->histZetaMcos[m][n1][n2]);
            }
            WRITE_OUTPUT_OR_FAIL(outstr, namebuf, "\n");
        }
        CLOSE_OUTPUT_OR_FAIL(outstr, namebuf);
    }
    for (m=1; m<=cmd->mChebyshev+1; m++) {
        if (format_checked(namebuf, sizeof(namebuf),
            "namebuf", "%s%s_%d%s", gd->fpfnamehistZetaMFileName,
                           "_sin", m, EXTFILES) != 0)
            return FAILURE;

        verb_print_q(2, cmd->verbose, "Printing : to a file %s ...\n",namebuf);
        OPEN_OUTPUT_OR_FAIL(outstr, namebuf, "w!");
        for (n1=1; n1<=cmd->sizeHistN; n1++) {
            for (n2=1; n2<=cmd->sizeHistN; n2++) {
                WRITE_OUTPUT_OR_FAIL(outstr, namebuf, "%16.8e ",
                                     gd->histZetaMsin[m][n1][n2]);
            }
            WRITE_OUTPUT_OR_FAIL(outstr, namebuf, "\n");
        }
        CLOSE_OUTPUT_OR_FAIL(outstr, namebuf);
    }
    for (m=1; m<=cmd->mChebyshev+1; m++) {
        if (format_checked(namebuf, sizeof(namebuf),
            "namebuf", "%s%s_%d%s", gd->fpfnamehistZetaMFileName,
                           "_sincos", m, EXTFILES) != 0)
            return FAILURE;

        verb_print_q(2, cmd->verbose, "Printing : to a file %s ...\n",namebuf);
        OPEN_OUTPUT_OR_FAIL(outstr, namebuf, "w!");
        for (n1=1; n1<=cmd->sizeHistN; n1++) {
            for (n2=1; n2<=cmd->sizeHistN; n2++) {
                WRITE_OUTPUT_OR_FAIL(outstr, namebuf, "%16.8e ",
                                     gd->histZetaMsincos[m][n1][n2]);
            }
            WRITE_OUTPUT_OR_FAIL(outstr, namebuf, "\n");
        }
        CLOSE_OUTPUT_OR_FAIL(outstr, namebuf);
    }
    // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
    for (m=1; m<=cmd->mChebyshev+1; m++) {
        if (format_checked(namebuf, sizeof(namebuf),
            "namebuf", "%s%s_%d%s", gd->fpfnamehistZetaMFileName,
                           "_cossin", m, EXTFILES) != 0)
            return FAILURE;
        verb_print_q(2, cmd->verbose, "Printing : to a file %s ...\n",namebuf);
        OPEN_OUTPUT_OR_FAIL(outstr, namebuf, "w!");
        for (n1=1; n1<=cmd->sizeHistN; n1++) {
            for (n2=1; n2<=cmd->sizeHistN; n2++) {
                WRITE_OUTPUT_OR_FAIL(outstr, namebuf, "%16.8e ",
                                     gd->histZetaMcossin[m][n1][n2]);
            }
            WRITE_OUTPUT_OR_FAIL(outstr, namebuf, "\n");
        }
        CLOSE_OUTPUT_OR_FAIL(outstr, namebuf);
    }

    return SUCCESS;
}


// Saves matrix ZetaM for each m multipole at a set of theta2 angles
local int PrintHistZetaMm_sincos(struct  cmdline_data* cmd,
                                struct  global_data* gd,
                                 gdlptr_sincos_omp_ggg gdl)
{
    string routineName = "PrintHistZetaMm_sincos";
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
        if (format_checked(namebuf, sizeof(namebuf),
            "namebuf", "%s%s_%d%s", gd->fpfnamemhistZetaMFileName,
                           "_cos", m, EXTFILES) != 0)
            return FAILURE;
        verb_print_q(2, cmd->verbose, "Printing : to a file %s ...\n",namebuf);
        OPEN_OUTPUT_OR_FAIL(outstr, namebuf, "w!");
        WRITE_OUTPUT_OR_FAIL(outstr, namebuf, MHISTZETAHEADER);
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
            WRITE_OUTPUT_OR_FAIL(outstr, namebuf,
                                 MHISTZETA,rBin,Zeta,Zeta2,Zeta3,Zeta4,Zeta5);
        }
        CLOSE_OUTPUT_OR_FAIL(outstr, namebuf);
    }
        
    for (m = 1; m <= cmd->mChebyshev+1; m++) {
        if (format_checked(namebuf, sizeof(namebuf),
            "namebuf", "%s%s_%d%s", gd->fpfnamemhistZetaMFileName,
                           "_sin", m, EXTFILES) != 0)
            return FAILURE;
        verb_print_q(2, cmd->verbose, "Printing : to a file %s ...\n",namebuf);
        OPEN_OUTPUT_OR_FAIL(outstr, namebuf, "w!");
        WRITE_OUTPUT_OR_FAIL(outstr, namebuf, MHISTZETAHEADER);
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
            WRITE_OUTPUT_OR_FAIL(outstr, namebuf,
                             MHISTZETA,rBin,Zeta,Zeta2,Zeta3,Zeta4,Zeta5);
        }
        CLOSE_OUTPUT_OR_FAIL(outstr, namebuf);
    }

    for (m = 1; m <= cmd->mChebyshev+1; m++) {
        if (format_checked(namebuf, sizeof(namebuf),
            "namebuf", "%s%s_%d%s", gd->fpfnamemhistZetaMFileName,
                           "_sincos", m, EXTFILES) != 0)
            return FAILURE;
        verb_print_q(2, cmd->verbose, "Printing : to a file %s ...\n",namebuf);
        OPEN_OUTPUT_OR_FAIL(outstr, namebuf, "w!");
        WRITE_OUTPUT_OR_FAIL(outstr, namebuf, MHISTZETAHEADER);
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
            WRITE_OUTPUT_OR_FAIL(outstr, namebuf,
                                 MHISTZETA,rBin,Zeta,Zeta2,Zeta3,Zeta4,Zeta5);
        }
        CLOSE_OUTPUT_OR_FAIL(outstr, namebuf);
    }

    // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
    for (m = 1; m <= cmd->mChebyshev+1; m++) {
        if (format_checked(namebuf, sizeof(namebuf),
            "namebuf", "%s%s_%d%s", gd->fpfnamemhistZetaMFileName,
                           "_cossin", m, EXTFILES) != 0)
            return FAILURE;
        verb_print_q(2, cmd->verbose, "Printing : to a file %s ...\n",namebuf);
        OPEN_OUTPUT_OR_FAIL(outstr, namebuf, "w!");
        WRITE_OUTPUT_OR_FAIL(outstr, namebuf, MHISTZETAHEADER);
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
            WRITE_OUTPUT_OR_FAIL(outstr, namebuf,
                                 MHISTZETA,rBin,Zeta,Zeta2,Zeta3,Zeta4,Zeta5);
        }
        CLOSE_OUTPUT_OR_FAIL(outstr, namebuf);
    }

    return SUCCESS;
}


// Saves matrix ZetaG, full correlation function at each phi bins
local int PrintHistZetaG(struct  cmdline_data* cmd,
                        struct  global_data* gd, gdlptr_sincos_omp_ggg gdl)
{
    string routineName = "PrintHistZetaG";
    int n1, n2, l, m;
    stream outstr;
    char namebuf[256];
    real theta;
    real ***histXi3pcfIm;

    histXi3pcfIm =
        dmatrix3D(1,cmd->sizeHistPhi,1,cmd->sizeHistN,1,cmd->sizeHistN);

    for (l=1; l<=cmd->sizeHistPhi; l++) {
        if (format_checked(namebuf, sizeof(namebuf),
            "namebuf", "%s%s_%d%s", gd->fpfnamehistZetaGFileName,
                           "_Xi3pcf_",l, EXTFILES) != 0)
            return FAILURE;
        theta = (real)l * gd->deltaPhi;
        verb_print_q(2, cmd->verbose,
                    "Printing : to a file %s with theta %g...\n",namebuf, theta);
        OPEN_OUTPUT_OR_FAIL(outstr, namebuf, "w!");
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
                WRITE_OUTPUT_OR_FAIL(outstr, namebuf,
                                     "%16.8e ",gd->histXi3pcf[l][n1][n2]);
            }
            WRITE_OUTPUT_OR_FAIL(outstr, namebuf, "\n");
        }
        CLOSE_OUTPUT_OR_FAIL(outstr, namebuf);
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
    string routineName = "PrintHistZetaGm_sincos";
    int n1, n2, m;
    stream outstr;
    char namebuf[256];

    for (m=1; m<=cmd->mChebyshev+1; m++) {
        if (format_checked(namebuf, sizeof(namebuf),
            "namebuf", "%s%s_%d%s", gd->fpfnamehistZetaGmFileName,
                           "_Re", m, EXTFILES) != 0)
            return FAILURE;
        verb_print_q(2, cmd->verbose, "Printing : to a file %s ...\n",namebuf);
        OPEN_OUTPUT_OR_FAIL(outstr, namebuf, "w!");
        for (n1=1; n1<=cmd->sizeHistN; n1++) {
            for (n2=1; n2<=cmd->sizeHistN; n2++) {
                WRITE_OUTPUT_OR_FAIL(outstr, namebuf,
                                     "%16.8e ",gd->histZetaGmRe[m][n1][n2]);
            }
            WRITE_OUTPUT_OR_FAIL(outstr, namebuf, "\n");
        }
        CLOSE_OUTPUT_OR_FAIL(outstr, namebuf);
    }
    for (m=1; m<=cmd->mChebyshev+1; m++) {
        if (format_checked(namebuf, sizeof(namebuf),
            "namebuf", "%s%s_%d%s", gd->fpfnamehistZetaGmFileName,
                           "_Im", m, EXTFILES) != 0)
            return FAILURE;
        verb_print_q(2, cmd->verbose, "Printing : to a file %s ...\n",namebuf);
        OPEN_OUTPUT_OR_FAIL(outstr, namebuf, "w!");
        for (n1=1; n1<=cmd->sizeHistN; n1++) {
            for (n2=1; n2<=cmd->sizeHistN; n2++) {
                WRITE_OUTPUT_OR_FAIL(outstr, namebuf,
                                     "%16.8e ",gd->histZetaGmIm[m][n1][n2]);
            }
            WRITE_OUTPUT_OR_FAIL(outstr, namebuf, "\n");
        }
        CLOSE_OUTPUT_OR_FAIL(outstr, namebuf);
    }

    return SUCCESS;
}


// Saves matrix ZetaG, real and imaginary parts, obtained from ZetaM multipoles
//  also saves full 3pcf ZetaG matrix for each phi bins obtained from inverse FFT
local int PrintHistZetaMZetaGm_sincos(struct  cmdline_data* cmd,
                                     struct  global_data* gd,
                                      gdlptr_sincos_omp_ggg gdl)
{
#define OPEN_OUTPUT_OR_FAIL_local(outstr, filename, mode)              \
    do {                                                               \
        if (stropen_checked((filename), (mode), &(outstr),             \
                            cmd->error_message, _ERRORMSGSIZE_)        \
            == FAILURE)                                                \
            goto cleanup;                                            \
    } while (0)

    string routineName = "PrintHistZetaMZetaGm_sincos";
    int n1, n2, m, l;
    stream outstr;
    char namebuf[256];

    int NP = 2*(cmd->mChebyshev+1);

    //B
    int status = FAILURE;
    double ***histZetaG = NULL;
    double ***histZetaG_Im = NULL;
    double *data = NULL;

    #ifdef USEGSL
    gsl_fft_real_wavetable *real = NULL;
    gsl_fft_real_workspace *work = NULL;
    gsl_fft_halfcomplex_wavetable *hc = NULL;
    #endif
    //E
    
    histZetaG = dmatrix3D(1,NP,1,cmd->sizeHistN,1,cmd->sizeHistN);
    histZetaG_Im = dmatrix3D(1,NP,1,cmd->sizeHistN,1,cmd->sizeHistN);

#ifdef USEGSL
    //B Test and check this allocation of memory...
    data=(double *)allocate(NP*sizeof(double));
    //E

    work = gsl_fft_real_workspace_alloc (NP);
    real = gsl_fft_real_wavetable_alloc (NP);
    hc = gsl_fft_halfcomplex_wavetable_alloc (NP);
#else
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
        if (format_checked(namebuf, sizeof(namebuf),
            "namebuf", "%s%s_%d%s", gd->fpfnamehistZetaGmFileName,
                           "_Re", m, EXTFILES) != 0)
            return FAILURE;
        verb_print_q(2, cmd->verbose,
                    "Printing : to a file %s ...\n",namebuf);
        OPEN_OUTPUT_OR_FAIL_local(outstr, namebuf, "w!");
        for (n1=1; n1<=cmd->sizeHistN; n1++) {
            for (n2=1; n2<=cmd->sizeHistN; n2++) {
                WRITE_OUTPUT_OR_FAIL(outstr, namebuf,
                                     "%16.8e ",gd->histZetaGmRe[m][n1][n2]);
            }
            WRITE_OUTPUT_OR_FAIL(outstr, namebuf, "\n");
        }
        CLOSE_OUTPUT_OR_FAIL(outstr, namebuf);
    }

    for (m=1; m<=cmd->mChebyshev+1; m++) {
        if (format_checked(namebuf, sizeof(namebuf),
            "namebuf", "%s%s_%d%s", gd->fpfnamehistZetaGmFileName,
                           "_Im", m, EXTFILES) != 0)
            return FAILURE;
        verb_print_q(2, cmd->verbose,
                    "Printing : to a file %s ...\n",namebuf);
        OPEN_OUTPUT_OR_FAIL_local(outstr, namebuf, "w!");
        for (n1=1; n1<=cmd->sizeHistN; n1++) {
            for (n2=1; n2<=cmd->sizeHistN; n2++) {
                WRITE_OUTPUT_OR_FAIL(outstr, namebuf,
                                     "%16.8e ",gd->histZetaGmIm[m][n1][n2]);
            }
            WRITE_OUTPUT_OR_FAIL(outstr, namebuf, "\n");
        }
        CLOSE_OUTPUT_OR_FAIL(outstr, namebuf);
    }

    for (l=1; l<=cmd->mChebyshev+1; l++) {
        if (format_checked(namebuf, sizeof(namebuf),
            "namebuf", "%s%s_%d%s", gd->fpfnamehistZetaGFileName,
                           "_fftinv_Re",l, EXTFILES) != 0)
            return FAILURE;
        verb_print_q(2, cmd->verbose,
                    "Printing : to a file %s ...\n",namebuf);
        OPEN_OUTPUT_OR_FAIL_local(outstr, namebuf, "w!");
        for (n1=1; n1<=cmd->sizeHistN; n1++) {
            for (n2=1; n2<=cmd->sizeHistN; n2++) {
                WRITE_OUTPUT_OR_FAIL(outstr, namebuf,
                                     "%16.8e ",histZetaG[2*l-1][n1][n2]);
            }
            WRITE_OUTPUT_OR_FAIL(outstr, namebuf, "\n");
        }
        CLOSE_OUTPUT_OR_FAIL(outstr, namebuf);
    }
    for (l=1; l<=cmd->mChebyshev+1; l++) {
        if (format_checked(namebuf, sizeof(namebuf),
            "namebuf", "%s%s_%d%s", gd->fpfnamehistZetaGFileName,
                           "_fftinv_Im",l, EXTFILES) != 0)
            return FAILURE;
        verb_print_q(2, cmd->verbose,
                    "Printing : to a file %s ...\n",namebuf);
        OPEN_OUTPUT_OR_FAIL_local(outstr, namebuf, "w!");
        for (n1=1; n1<=cmd->sizeHistN; n1++) {
            for (n2=1; n2<=cmd->sizeHistN; n2++) {
                WRITE_OUTPUT_OR_FAIL(outstr, namebuf,
                                     "%16.8e ",histZetaG[2*l][n1][n2]);
            }
            WRITE_OUTPUT_OR_FAIL(outstr, namebuf, "\n");
        }
        CLOSE_OUTPUT_OR_FAIL(outstr, namebuf);
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
        if (format_checked(namebuf, sizeof(namebuf),
            "namebuf", "%s_%s_%d%s",
                           gd->fpfnamehistZetaGFileName, "Re", l, EXTFILES) != 0)
            return FAILURE;
        verb_print_q(2, cmd->verbose,
                    "Printing : to a file %s ...\n",namebuf);
        OPEN_OUTPUT_OR_FAIL_local(outstr, namebuf, "w!");
        for (n1=1; n1<=cmd->sizeHistN; n1++) {
            for (n2=1; n2<=cmd->sizeHistN; n2++) {
                WRITE_OUTPUT_OR_FAIL(outstr, namebuf,
                                     "%16.8e ",histZetaG[l][n1][n2]);
            }
            WRITE_OUTPUT_OR_FAIL(outstr, namebuf, "\n");
        }
        CLOSE_OUTPUT_OR_FAIL(outstr, namebuf);
    }
    for (l=1; l<=NP; l++) {
        if (format_checked(namebuf, sizeof(namebuf),
            "namebuf", "%s_%s_%d%s",
                           gd->fpfnamehistZetaGFileName, "Im", l, EXTFILES) != 0)
            return FAILURE;
        verb_print_q(2, cmd->verbose,
                    "Printing : to a file %s ...\n",namebuf);
        OPEN_OUTPUT_OR_FAIL_local(outstr, namebuf, "w!");
        for (n1=1; n1<=cmd->sizeHistN; n1++) {
            for (n2=1; n2<=cmd->sizeHistN; n2++) {
                WRITE_OUTPUT_OR_FAIL(outstr, namebuf,
                                     "%16.8e ",histZetaG_Im[l][n1][n2]);
            }
            WRITE_OUTPUT_OR_FAIL(outstr, namebuf, "\n");
        }
        CLOSE_OUTPUT_OR_FAIL(outstr, namebuf);
    }
    
    status = SUCCESS;

cleanup:
    #ifdef USEGSL
    if (hc) gsl_fft_halfcomplex_wavetable_free(hc);
    if (real) gsl_fft_real_wavetable_free(real);
    if (work) gsl_fft_real_workspace_free(work);
    if (data) free(data);
    #else
    if (data) free_dvector(data, 1, NP);
    #endif

    if (histZetaG_Im)
        free_dmatrix3D(histZetaG_Im, 1, NP, 1, cmd->sizeHistN, 1, cmd->sizeHistN);
    if (histZetaG)
        free_dmatrix3D(histZetaG, 1, NP, 1, cmd->sizeHistN, 1, cmd->sizeHistN);

#undef OPEN_OUTPUT_OR_FAIL_local

    return status;
}

#endif // ! THREEPCFCONVERGENCE


#ifdef NMultipoles

#ifdef THREEPCFCONVERGENCE
// Saves matrix ZetaM for each m multipole
local int PrintHistZetaM_sincos_N(struct  cmdline_data* cmd,
                                struct  global_data* gd,
                                  gdlptr_sincos_omp_ggg_N gdlN)
{
    string routineName = "PrintHistZetaM_sincos_N";
    int n1, n2, m;
    stream outstr;
    char namebuf[256];

    for (m=1; m<=GGG_WINDOW_ORDERS(cmd); m++) {
        if (format_checked(namebuf, sizeof(namebuf),
            "namebuf", "%s%s_%d%s", gd->fpfnamehistZetaMFileName,
                           "_cos_N", m, EXTFILES) != 0)
            return FAILURE;
        verb_print_q(2, cmd->verbose, "Printing : to a file %s ...\n",namebuf);
        OPEN_OUTPUT_OR_FAIL(outstr, namebuf, "w!");
        for (n1=1; n1<=cmd->sizeHistN; n1++) {
            for (n2=1; n2<=cmd->sizeHistN; n2++) {
                WRITE_OUTPUT_OR_FAIL(outstr, namebuf,
                                     "%16.8e ",gdlN->histZetaMcos[m][n1][n2]);
            }
            WRITE_OUTPUT_OR_FAIL(outstr, namebuf, "\n");
        }
        CLOSE_OUTPUT_OR_FAIL(outstr, namebuf);
    }
    for (m=1; m<=GGG_WINDOW_ORDERS(cmd); m++) {
        if (format_checked(namebuf, sizeof(namebuf),
            "namebuf", "%s%s_%d%s", gd->fpfnamehistZetaMFileName,
                           "_sin_N", m, EXTFILES) != 0)
            return FAILURE;
        verb_print_q(2, cmd->verbose, "Printing : to a file %s ...\n",namebuf);
        OPEN_OUTPUT_OR_FAIL(outstr, namebuf, "w!");
       for (n1=1; n1<=cmd->sizeHistN; n1++) {
            for (n2=1; n2<=cmd->sizeHistN; n2++) {
                WRITE_OUTPUT_OR_FAIL(outstr, namebuf,
                                     "%16.8e ",gdlN->histZetaMsin[m][n1][n2]);
            }
           WRITE_OUTPUT_OR_FAIL(outstr, namebuf, "\n");
        }
        CLOSE_OUTPUT_OR_FAIL(outstr, namebuf);
    }
    for (m=1; m<=GGG_WINDOW_ORDERS(cmd); m++) {
        if (format_checked(namebuf, sizeof(namebuf),
            "namebuf", "%s%s_%d%s", gd->fpfnamehistZetaMFileName,
                           "_sincos_N", m, EXTFILES) != 0)
            return FAILURE;
        verb_print_q(2, cmd->verbose, "Printing : to a file %s ...\n",namebuf);
        OPEN_OUTPUT_OR_FAIL(outstr, namebuf, "w!");
        for (n1=1; n1<=cmd->sizeHistN; n1++) {
            for (n2=1; n2<=cmd->sizeHistN; n2++) {
                WRITE_OUTPUT_OR_FAIL(outstr, namebuf,
                                     "%16.8e ",gdlN->histZetaMsincos[m][n1][n2]);
            }
            WRITE_OUTPUT_OR_FAIL(outstr, namebuf, "\n");
        }
        CLOSE_OUTPUT_OR_FAIL(outstr, namebuf);
    }
    // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
    for (m=1; m<=GGG_WINDOW_ORDERS(cmd); m++) {
        if (format_checked(namebuf, sizeof(namebuf),
            "namebuf", "%s%s_%d%s", gd->fpfnamehistZetaMFileName,
                           "_cossin_N", m, EXTFILES) != 0)
            return FAILURE;
        verb_print_q(2, cmd->verbose, "Printing : to a file %s ...\n",namebuf);
        OPEN_OUTPUT_OR_FAIL(outstr, namebuf, "w!");
       for (n1=1; n1<=cmd->sizeHistN; n1++) {
            for (n2=1; n2<=cmd->sizeHistN; n2++) {
                WRITE_OUTPUT_OR_FAIL(outstr, namebuf,
                                     "%16.8e ",gdlN->histZetaMcossin[m][n1][n2]);
            }
           WRITE_OUTPUT_OR_FAIL(outstr, namebuf, "\n");
        }
        CLOSE_OUTPUT_OR_FAIL(outstr, namebuf);
    }

    return SUCCESS;
}

// Saves matrix ZetaM for each m multipole at a set of theta2 angles
local int PrintHistZetaMm_sincos_N(struct  cmdline_data* cmd,
                                struct  global_data* gd,
                                   gdlptr_sincos_omp_ggg_N gdlN)
{
    string routineName = "PrintHistZetaMm_sincos_N";
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

    for (m = 1; m <= GGG_WINDOW_ORDERS(cmd); m++) {
        if (format_checked(namebuf, sizeof(namebuf),
            "namebuf", "%s%s_%d%s", gd->fpfnamemhistZetaMFileName,
                           "_cos_N", m, EXTFILES) != 0)
            return FAILURE;
        verb_print_q(2, cmd->verbose, "Printing : to a file %s ...\n",namebuf);
        OPEN_OUTPUT_OR_FAIL(outstr, namebuf, "w!");
        WRITE_OUTPUT_OR_FAIL(outstr, namebuf, MHISTZETAHEADER);
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
            WRITE_OUTPUT_OR_FAIL(outstr, namebuf,
                                 MHISTZETA,rBin,Zeta,Zeta2,Zeta3,Zeta4,Zeta5);
        }
        CLOSE_OUTPUT_OR_FAIL(outstr, namebuf);
    }
        
    for (m = 1; m <= GGG_WINDOW_ORDERS(cmd); m++) {
        if (format_checked(namebuf, sizeof(namebuf),
            "namebuf", "%s%s_%d%s", gd->fpfnamemhistZetaMFileName,
                           "_sin_N", m, EXTFILES) != 0)
            return FAILURE;
        verb_print_q(2, cmd->verbose, "Printing : to a file %s ...\n",namebuf);
        OPEN_OUTPUT_OR_FAIL(outstr, namebuf, "w!");
        WRITE_OUTPUT_OR_FAIL(outstr, namebuf, MHISTZETAHEADER);
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
                WRITE_OUTPUT_OR_FAIL(outstr, namebuf,
                             MHISTZETA,rBin,Zeta,Zeta2,Zeta3,Zeta4,Zeta5);
        }
        CLOSE_OUTPUT_OR_FAIL(outstr, namebuf);
    }

    for (m = 1; m <= GGG_WINDOW_ORDERS(cmd); m++) {
        if (format_checked(namebuf, sizeof(namebuf),
            "namebuf", "%s%s_%d%s", gd->fpfnamemhistZetaMFileName,
                           "_sincos_N", m, EXTFILES) != 0)
            return FAILURE;
        verb_print_q(2, cmd->verbose, "Printing : to a file %s ...\n",namebuf);
        OPEN_OUTPUT_OR_FAIL(outstr, namebuf, "w!");
        WRITE_OUTPUT_OR_FAIL(outstr, namebuf, MHISTZETAHEADER);
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
            WRITE_OUTPUT_OR_FAIL(outstr, namebuf,
                                 MHISTZETA,rBin,Zeta,Zeta2,Zeta3,Zeta4,Zeta5);
        }
        CLOSE_OUTPUT_OR_FAIL(outstr, namebuf);
    }

    // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
    for (m = 1; m <= GGG_WINDOW_ORDERS(cmd); m++) {
        if (format_checked(namebuf, sizeof(namebuf),
            "namebuf", "%s%s_%d%s", gd->fpfnamemhistZetaMFileName,
                           "_cossin_N", m, EXTFILES) != 0)
            return FAILURE;
        verb_print_q(2, cmd->verbose, "Printing : to a file %s ...\n",namebuf);
        OPEN_OUTPUT_OR_FAIL(outstr, namebuf, "w!");
        WRITE_OUTPUT_OR_FAIL(outstr, namebuf, MHISTZETAHEADER);
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
            WRITE_OUTPUT_OR_FAIL(outstr, namebuf,
                                 MHISTZETA,rBin,Zeta,Zeta2,Zeta3,Zeta4,Zeta5);
        }
        CLOSE_OUTPUT_OR_FAIL(outstr, namebuf);
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
    string routineName = "PrintHistZetaM_sincos_normalized";
    int n1, n2, m;
    stream outstr;
    char namebuf[256];

    for (m=1; m<=cmd->mChebyshev+1; m++) {
        if (format_checked(namebuf, sizeof(namebuf),
            "namebuf", "%s%s_%d%s", gd->fpfnamehistZetaMFileName,
                           "_cos", m, EXTFILES) != 0)
            return FAILURE;
        verb_print_q(2, cmd->verbose, "Printing : to a file %s ...\n",namebuf);
        OPEN_OUTPUT_OR_FAIL(outstr, namebuf, "w!");
        for (n1=1; n1<=cmd->sizeHistN; n1++) {
            for (n2=1; n2<=cmd->sizeHistN; n2++) {
                WRITE_OUTPUT_OR_FAIL(outstr, namebuf, "%16.8e ",
                    normalize_zeta_ggg(gdl->histZetaMcos[m][n1][n2],
                                       gdlN, n1, n2));
            }
            WRITE_OUTPUT_OR_FAIL(outstr, namebuf, "\n");
        }
        CLOSE_OUTPUT_OR_FAIL(outstr, namebuf);
    }
    for (m=1; m<=cmd->mChebyshev+1; m++) {
        if (format_checked(namebuf, sizeof(namebuf),
            "namebuf", "%s%s_%d%s", gd->fpfnamehistZetaMFileName,
                           "_sin", m, EXTFILES) != 0)
            return FAILURE;
        verb_print_q(2, cmd->verbose, "Printing : to a file %s ...\n",namebuf);
        OPEN_OUTPUT_OR_FAIL(outstr, namebuf, "w!");
        for (n1=1; n1<=cmd->sizeHistN; n1++) {
            for (n2=1; n2<=cmd->sizeHistN; n2++) {
                WRITE_OUTPUT_OR_FAIL(outstr, namebuf, "%16.8e ",
                    normalize_zeta_ggg(gdl->histZetaMsin[m][n1][n2],
                                       gdlN, n1, n2));
            }
            WRITE_OUTPUT_OR_FAIL(outstr, namebuf, "\n");
        }
        CLOSE_OUTPUT_OR_FAIL(outstr, namebuf);
    }
    for (m=1; m<=cmd->mChebyshev+1; m++) {
        if (format_checked(namebuf, sizeof(namebuf),
            "namebuf", "%s%s_%d%s", gd->fpfnamehistZetaMFileName,
                           "_sincos", m, EXTFILES) != 0)
            return FAILURE;
        verb_print_q(2, cmd->verbose, "Printing : to a file %s ...\n",namebuf);
        OPEN_OUTPUT_OR_FAIL(outstr, namebuf, "w!");
        for (n1=1; n1<=cmd->sizeHistN; n1++) {
            for (n2=1; n2<=cmd->sizeHistN; n2++) {
                WRITE_OUTPUT_OR_FAIL(outstr, namebuf, "%16.8e ",
                    normalize_zeta_ggg(gdl->histZetaMsincos[m][n1][n2],
                                       gdlN, n1, n2));
            }
            WRITE_OUTPUT_OR_FAIL(outstr, namebuf, "\n");
        }
        CLOSE_OUTPUT_OR_FAIL(outstr, namebuf);
    }
    // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
    for (m=1; m<=cmd->mChebyshev+1; m++) {
        if (format_checked(namebuf, sizeof(namebuf),
            "namebuf", "%s%s_%d%s", gd->fpfnamehistZetaMFileName,
                           "_cossin", m, EXTFILES) != 0)
            return FAILURE;
        verb_print_q(2, cmd->verbose, "Printing : to a file %s ...\n",namebuf);
        OPEN_OUTPUT_OR_FAIL(outstr, namebuf, "w!");
        for (n1=1; n1<=cmd->sizeHistN; n1++) {
            for (n2=1; n2<=cmd->sizeHistN; n2++) {
                WRITE_OUTPUT_OR_FAIL(outstr, namebuf, "%16.8e ",
                    normalize_zeta_ggg(gdl->histZetaMcossin[m][n1][n2],
                                       gdlN, n1, n2));
            }
            WRITE_OUTPUT_OR_FAIL(outstr, namebuf, "\n");
        }
        CLOSE_OUTPUT_OR_FAIL(outstr, namebuf);
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
    string routineName = "PrintHistZetaMm_sincos_normalized";
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
        if (format_checked(namebuf, sizeof(namebuf),
            "namebuf", "%s%s_%d%s", gd->fpfnamemhistZetaMFileName,
                           "_cos", m, EXTFILES) != 0)
            return FAILURE;
        verb_print_q(2, cmd->verbose, "Printing : to a file %s ...\n",namebuf);
        OPEN_OUTPUT_OR_FAIL(outstr, namebuf, "w!");
        WRITE_OUTPUT_OR_FAIL(outstr, namebuf, MHISTZETAHEADER);
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
            Zeta = normalize_zeta_ggg(gdl->histZetaMcos[m][n1][n1],
                                      gdlN, n1, n1);
            Zeta2 = normalize_zeta_ggg(
                gdl->histZetaMcos[m][n1][(int)(Nbins/4.0)],
                gdlN, n1, (int)(Nbins/4.0));
            Zeta3 = normalize_zeta_ggg(
                gdl->histZetaMcos[m][n1][(int)(2.0*Nbins/4.0)],
                gdlN, n1, (int)(2.0*Nbins/4.0));
            Zeta4 = normalize_zeta_ggg(
                gdl->histZetaMcos[m][n1][(int)(3.0*Nbins/4.0)],
                gdlN, n1, (int)(3.0*Nbins/4.0));
            Zeta5 = normalize_zeta_ggg(
                gdl->histZetaMcos[m][n1][(int)(4.0*Nbins/4.0 - 1.0)],
                gdlN, n1, (int)(4.0*Nbins/4.0 - 1.0));
            WRITE_OUTPUT_OR_FAIL(outstr, namebuf,
                                 MHISTZETA,rBin,Zeta,Zeta2,Zeta3,Zeta4,Zeta5);
        }
        CLOSE_OUTPUT_OR_FAIL(outstr, namebuf);
    }
        
    for (m = 1; m <= cmd->mChebyshev+1; m++) {
        if (format_checked(namebuf, sizeof(namebuf),
            "namebuf", "%s%s_%d%s", gd->fpfnamemhistZetaMFileName,
                           "_sin", m, EXTFILES) != 0)
            return FAILURE;
        verb_print_q(2, cmd->verbose, "Printing : to a file %s ...\n",namebuf);
        OPEN_OUTPUT_OR_FAIL(outstr, namebuf, "w!");
        WRITE_OUTPUT_OR_FAIL(outstr, namebuf, MHISTZETAHEADER);
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
            Zeta = normalize_zeta_ggg(gdl->histZetaMsin[m][n1][n1],
                                      gdlN, n1, n1);
            Zeta2 = normalize_zeta_ggg(
                gdl->histZetaMsin[m][n1][(int)(Nbins/4.0)],
                gdlN, n1, (int)(Nbins/4.0));
            Zeta3 = normalize_zeta_ggg(
                gdl->histZetaMsin[m][n1][(int)(2.0*Nbins/4.0)],
                gdlN, n1, (int)(2.0*Nbins/4.0));
            Zeta4 = normalize_zeta_ggg(
                gdl->histZetaMsin[m][n1][(int)(3.0*Nbins/4.0)],
                gdlN, n1, (int)(3.0*Nbins/4.0));
            Zeta5 = normalize_zeta_ggg(
                gdl->histZetaMsin[m][n1][(int)(4.0*Nbins/4.0 - 1.0)],
                gdlN, n1, (int)(4.0*Nbins/4.0 - 1.0));
            WRITE_OUTPUT_OR_FAIL(outstr, namebuf,
                                 MHISTZETA,rBin,Zeta,Zeta2,Zeta3,Zeta4,Zeta5);
        }
        CLOSE_OUTPUT_OR_FAIL(outstr, namebuf);
    }

    for (m = 1; m <= cmd->mChebyshev+1; m++) {
        if (format_checked(namebuf, sizeof(namebuf),
            "namebuf", "%s%s_%d%s", gd->fpfnamemhistZetaMFileName,
                           "_sincos", m, EXTFILES) != 0)
            return FAILURE;
        verb_print_q(2, cmd->verbose, "Printing : to a file %s ...\n",namebuf);
        OPEN_OUTPUT_OR_FAIL(outstr, namebuf, "w!");
        WRITE_OUTPUT_OR_FAIL(outstr, namebuf, MHISTZETAHEADER);
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
            Zeta = normalize_zeta_ggg(gdl->histZetaMsincos[m][n1][n1],
                                      gdlN, n1, n1);
            Zeta2 = normalize_zeta_ggg(
                gdl->histZetaMsincos[m][n1][(int)(Nbins/4.0)],
                gdlN, n1, (int)(Nbins/4.0));
            Zeta3 = normalize_zeta_ggg(
                gdl->histZetaMsincos[m][n1][(int)(2.0*Nbins/4.0)],
                gdlN, n1, (int)(2.0*Nbins/4.0));
            Zeta4 = normalize_zeta_ggg(
                gdl->histZetaMsincos[m][n1][(int)(3.0*Nbins/4.0)],
                gdlN, n1, (int)(3.0*Nbins/4.0));
            Zeta5 = normalize_zeta_ggg(
                gdl->histZetaMsincos[m][n1][(int)(4.0*Nbins/4.0 - 1.0)],
                gdlN, n1, (int)(4.0*Nbins/4.0 - 1.0));

            WRITE_OUTPUT_OR_FAIL(outstr, namebuf,
                                 MHISTZETA,rBin,Zeta,Zeta2,Zeta3,Zeta4,Zeta5);
        }
        CLOSE_OUTPUT_OR_FAIL(outstr, namebuf);
    }

    // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
    for (m = 1; m <= cmd->mChebyshev+1; m++) {
        if (format_checked(namebuf, sizeof(namebuf),
            "namebuf", "%s%s_%d%s", gd->fpfnamemhistZetaMFileName,
                           "_cossin", m, EXTFILES) != 0)
            return FAILURE;
        verb_print_q(2, cmd->verbose, "Printing : to a file %s ...\n",namebuf);
        OPEN_OUTPUT_OR_FAIL(outstr, namebuf, "w!");
        WRITE_OUTPUT_OR_FAIL(outstr, namebuf, MHISTZETAHEADER);
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
            Zeta = normalize_zeta_ggg(gdl->histZetaMcossin[m][n1][n1],
                                      gdlN, n1, n1);
            Zeta2 = normalize_zeta_ggg(
                gdl->histZetaMcossin[m][n1][(int)(Nbins/4.0)],
                gdlN, n1, (int)(Nbins/4.0));
            Zeta3 = normalize_zeta_ggg(
                gdl->histZetaMcossin[m][n1][(int)(2.0*Nbins/4.0)],
                gdlN, n1, (int)(2.0*Nbins/4.0));
            Zeta4 = normalize_zeta_ggg(
                gdl->histZetaMcossin[m][n1][(int)(3.0*Nbins/4.0)],
                gdlN, n1, (int)(3.0*Nbins/4.0));
            Zeta5 = normalize_zeta_ggg(
                gdl->histZetaMcossin[m][n1][(int)(4.0*Nbins/4.0 - 1.0)],
                gdlN, n1, (int)(4.0*Nbins/4.0 - 1.0));

            WRITE_OUTPUT_OR_FAIL(outstr, namebuf,
                                 MHISTZETA,rBin,Zeta,Zeta2,Zeta3,Zeta4,Zeta5);
        }
        CLOSE_OUTPUT_OR_FAIL(outstr, namebuf);
    }

    return SUCCESS;
}

// Saves matrix ZetaM for each m multipole
local int PrintHistZetaM_sincos_edge_effects(struct  cmdline_data* cmd,
                                              struct  global_data* gd,
                                              gdlptr_sincos_omp_ggg gdl,
                                              gdlptr_sincos_omp_ggg_N gdlN)
{
    string routineName = "PrintHistZetaM_sincos_edge_effects";
    int n1, n2, m;
    stream outstr;
    char namebuf[256];
    real rBin, rbinlog;
    (void)routineName;
    (void)gdlN;

    for (m=1; m<=cmd->mChebyshev+1; m++) {
        if (format_checked(namebuf, sizeof(namebuf),
            "namebuf", "%s%s_%d%s", gd->fpfnamehistZetaMFileName,
                           "_EE", m, EXTFILES) != 0)
            return FAILURE;
        verb_print_q(2, cmd->verbose, "Printing : to a file %s ...\n",namebuf);
        OPEN_OUTPUT_OR_FAIL(outstr, namebuf, "w!");
        for (n1=1; n1<=cmd->sizeHistN; n1++) {
            for (n2=1; n2<=cmd->sizeHistN; n2++) {
                WRITE_OUTPUT_OR_FAIL(outstr, namebuf, "%16.8e ",
                                     gdl->histZetaM_EE[m][n1][n2]);
            }
            WRITE_OUTPUT_OR_FAIL(outstr, namebuf, "\n");
        }
        CLOSE_OUTPUT_OR_FAIL(outstr, namebuf);
    }

    for (m=1; m<=cmd->mChebyshev+1; m++) {
        if (format_checked(namebuf, sizeof(namebuf), "namebuf",
                           "%s%s_%d%s", gd->fpfnamehistZetaMFileName,
                           "_EE_Im", m, EXTFILES) != 0)
            return FAILURE;
        OPEN_OUTPUT_OR_FAIL(outstr, namebuf, "w!");
        for (n1=1; n1<=cmd->sizeHistN; n1++) {
            for (n2=1; n2<=cmd->sizeHistN; n2++)
                WRITE_OUTPUT_OR_FAIL(outstr, namebuf, "%16.8e ",
                                     gdl->histZetaM_EE_Im[m][n1][n2]);
            WRITE_OUTPUT_OR_FAIL(outstr, namebuf, "\n");
        }
        CLOSE_OUTPUT_OR_FAIL(outstr, namebuf);
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
            if (format_checked(namebuf, sizeof(namebuf),
                "namebuf", "%s%s_%d%s", gd->fpfnamemhistZetaMFileName,
                               "_EE", m, EXTFILES) != 0)
                return FAILURE;
            verb_print_q(2, cmd->verbose, "Printing : to a file %s ...\n",
                         namebuf);
            OPEN_OUTPUT_OR_FAIL(outstr, namebuf, "w!");
            WRITE_OUTPUT_OR_FAIL(outstr, namebuf, MHISTZETAHEADER);
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
                Zeta = gdl->histZetaM_EE[m][n1][n1];
                Zeta2 = gdl->histZetaM_EE[m][n1][(int)(Nbins/4.0)];
                Zeta3 = gdl->histZetaM_EE[m][n1][(int)(2.0*Nbins/4.0)];
                Zeta4 = gdl->histZetaM_EE[m][n1][(int)(3.0*Nbins/4.0)];
                Zeta5 = gdl->histZetaM_EE[m][n1][(int)(4.0*Nbins/4.0 - 1.0)];
                WRITE_OUTPUT_OR_FAIL(outstr, namebuf,
                    MHISTZETA,rBin,Zeta,Zeta2,Zeta3,Zeta4,Zeta5);
            }
            CLOSE_OUTPUT_OR_FAIL(outstr, namebuf);
        }
    }
    //E

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
    string routineName = "PrintHistNN";
    real rBin, rbinlog;
    int n;
    stream outstr;

    OPEN_OUTPUT_OR_FAIL(outstr, gd->fpfnamehistNNFileName, "w!");

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
        WRITE_OUTPUT_OR_FAIL(outstr, gd->fpfnamehistNNFileName,
                             "%16.8e %16.8e\n",rBin,gdl->histNN[n]);
    }
    CLOSE_OUTPUT_OR_FAIL(outstr, gd->fpfnamehistNNFileName);

    if (cballs_opt_and_cf(cmd))
        PRINT_OR_FAIL(PrintHistCF(cmd, gd, gdl));

    return SUCCESS;
}

local int PrintHistCF(struct  cmdline_data* cmd, struct  global_data* gd,
                      gdlptr_sincos_omp_ggg gdl)
{
    string routineName = "PrintHistCF";
    real rBin, rbinlog;
    int n;
    stream outstr;

    OPEN_OUTPUT_OR_FAIL(outstr, gd->fpfnamehistCFFileName, "w!");

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
        WRITE_OUTPUT_OR_FAIL(outstr, gd->fpfnamehistCFFileName,
                             "%16.8e %16.8e\n",rBin,gdl->histCF[n]);
    }
    CLOSE_OUTPUT_OR_FAIL(outstr, gd->fpfnamehistCFFileName);

    return SUCCESS;
}

// in common_defs.h
// 180*60/Pi
//#define RADTOARCMIN   3437.74677
local int PrintHistXi2pcf(struct  cmdline_data* cmd, struct  global_data* gd,
                          gdlptr_sincos_omp_ggg gdl)
{
    string routineName = "PrintHistXi2pcf";
    real rBin, rbinlog;
    int n;
    stream outstr;
    char namebuf[256];

    if (format_checked(namebuf, sizeof(namebuf),
        "namebuf", "%s%s%s", gd->fpfnamehistXi2pcfFileName,
                       cmd->suffixOutFiles, EXTFILES) != 0)
        return FAILURE;
    verb_print_q(2, cmd->verbose, "Printing : to a file %s ...\n",namebuf);
    OPEN_OUTPUT_OR_FAIL(outstr, namebuf, "w!");

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
        if (cballs_opt_rbin_arcmin(cmd))
            rBin = rBin*RADTOARCMIN;
        else
            if (cballs_opt_rbin_degree(cmd))
                rBin = rBin*RADTOARCMIN/60.0;
        WRITE_OUTPUT_OR_FAIL(outstr, namebuf,
                             "%16.8e %16.8e\n",rBin,gdl->histXi2pcf[n]);
    }
    CLOSE_OUTPUT_OR_FAIL(outstr, namebuf);

    return SUCCESS;
}
//#undef RADTOARCMIN

#endif // ! TWOPCF


//E Saving histograms section: case GGGCORRELATION:

