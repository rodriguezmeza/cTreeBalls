/* ==============================================================================
 MODULE: search_octree_balls4_omp.c		[cTreeBalls]
 Written by: M.A. Rodriguez-Meza
 Starting date:    april 2023
 Purpose: 2- and 3-point correlation function computation
 Language: C
 Use: searchcalc_octree_balls4_omp(cmd, gd, btable, nbody,
                                           ipmin, ipmax, cat1, cat2);
 Major revisions:
 ==============================================================================*/
//        1          2          3          4        ^ 5          6          7


// Work to do in order to use with boxes not centered at (0,0,...)

#include "globaldefs.h"
#include <stdint.h>
#include "balls4_parallel.h"

#ifdef DEBUG
local bodyptr bodytabbf_balls4;
#endif

//B experimental, do not activate
#define PRAGMABCINTLIST
#undef PRAGMABCINTLIST
//E

//B What if mChebyshev is less than 7... correct!
#ifdef MANUALCHEBYSHEV
#define CHEBYSHEVTUOMPKKK                                         \
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
#define CHEBYSHEVTUOMPKKK                                         \
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
#define CHEBYSHEVTUOMPKKKN                                           \
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
#define CHEBYSHEVTUOMPKKKN                                           \
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

//B Define structures:
typedef struct {
#ifdef TWOPCF
    realptr histNN;
    realptr histWW;
    realptr histCF;
    realptr histNNSubXi2pcf;
    realptr histXi2pcf;
#endif
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
    INTEGER actlen;
    nodeptr *active;
    cellptr interact;
} gdhist_omp_balls6_kkk, *gdhistptr_omp_balls6_kkk;

typedef struct {
#ifdef TWOPCF
    realptr histNthread;
    realptr histWthread;
    realptr histWWthread;
    realptr histNNSubXi2pcfthread;
    realptr histXi2pcfthread;
    realptr histXi2pcfthreadsub;
#endif
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
    realptr histNNSubthread;

    real **histXithreadcos;
    real **histXithreadsin;

    INTEGER nbbcalcthread;
    INTEGER nbccalcthread;
    INTEGER ncccalcthread;
    INTEGER ibfcountthread;
    INTEGER nsmoothcountthread;

    compute_vector q0;
    real drpq2, drpq;
    compute_vector dr0;
    real cosb;
    real sinb;

    int thrID;

    INTEGER ipcount;
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
    realptr histNNSubthread;

    real **histXithreadcos;
    real **histXithreadsin;

    INTEGER nbbcalcthread;
    INTEGER nbccalcthread;
    INTEGER ncccalcthread;
    INTEGER ibfcountthread;
    INTEGER nsmoothcountthread;

    compute_vector q0;
    real drpq2, drpq;
    compute_vector dr0;
} gdhist_sincos_omp_kkk_N, *gdhistptr_sincos_omp_kkk_N;
#endif

//E Define structures:

local bool balls4_node_is_masked(struct cmdline_data *cmd, nodeptr q)
{
    return cballs_opt_read_mask(cmd) && Mask(q) == MASK_NODE_MASKED;
}

local bool balls4_cell_requires_open(struct cmdline_data *cmd, nodeptr q)
{
    return cballs_opt_read_mask(cmd) && Type(q) == CELL
        && !mask_node_can_approximate(Mask(q));
}

#ifdef TWOPCF
local real balls4_node_count(nodeptr q)
{
    return Type(q) == CELL ? (real)Nb(q) : 1.0;
}

local real balls4_node_weight_sum(nodeptr q)
{
#ifndef NOWKAvg
    return Weight(q);
#else
    return balls4_node_count(q)*Weight(q);
#endif
}

local real balls4_node_field_sum(nodeptr q)
{
#ifndef NOWKAvg
    return Weight(q)*Kappa(q);
#else
    return balls4_node_count(q)*Weight(q)*Kappa(q);
#endif
}

local void balls4_accumulate_neighbor_2pcf(
        gdhistptr_sincos_omp_kkk hist, nodeptr q, int n)
{
    real count = balls4_node_count(q);

    hist->histNthread[n] += count;
    hist->histWthread[n] += balls4_node_weight_sum(q);
    hist->histNNSubXi2pcfthread[n] += count;
    hist->histXi2pcfthreadsub[n] += balls4_node_field_sum(q);
}

local void balls4_accumulate_cell_pair_2pcf(
        gdhistptr_sincos_omp_kkk hist, nodeptr p, nodeptr q, int n)
{
    real pair_count = balls4_node_count(p)*balls4_node_count(q);

    hist->histNthread[n] += pair_count;
    hist->histWWthread[n] +=
        balls4_node_weight_sum(p)*balls4_node_weight_sum(q);
    hist->histNNSubXi2pcfthread[n] += pair_count;
    hist->histXi2pcfthread[n] +=
        balls4_node_field_sum(p)*balls4_node_field_sum(q);
}
#endif

//#ifndef BALLS
local bool nodes_condition_balls(struct cmdline_data* cmd, struct  global_data* gd,
                                 nodeptr p, nodeptr q, real *dr1,
                                 compute_vector dr);
local bool collapse_small_cell_balls4(struct cmdline_data* cmd,
                                      struct global_data* gd,
                                      nodeptr q, real distance);
//#endif
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
local void sumnode_balls6_omp(struct cmdline_data* cmd, struct  global_data* gd,
                              cellptr, cellptr, bodyptr,
                              gdhistptr_sincos_omp_kkk);
local void sumcell_balls6_omp(struct cmdline_data* cmd, struct  global_data* gd,
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

#ifdef NMultipoles
local void walktree_balls6_omp_N(struct cmdline_data* cmd, struct  global_data* gd,
                                 nodeptr *, nodeptr *, cellptr, cellptr,
                                 nodeptr, real, vector,
                                 gdhistptr_omp_balls6_kkk,
                                 gdhistptr_sincos_omp_kkk,
                                 gdhistptr_sincos_omp_kkk_N,
                                 INTEGER, int);
local void walksub6_omp_N(struct cmdline_data* cmd, struct  global_data* gd,
                        nodeptr *, nodeptr *, cellptr, cellptr,
                        nodeptr, real, vector,
                          gdhistptr_omp_balls6_kkk,
                          gdhistptr_sincos_omp_kkk,
                          gdhistptr_sincos_omp_kkk_N,
                          INTEGER, int);
local int sum_balls6_omp_N(struct cmdline_data* cmd, struct  global_data* gd,
                            cellptr, cellptr, bodyptr,
                            gdhistptr_omp_balls6_kkk,
                            gdhistptr_sincos_omp_kkk,
                            gdhistptr_sincos_omp_kkk_N,
                            INTEGER, int);
local void sumnode_balls6_omp_N(struct cmdline_data* cmd, struct  global_data* gd,
                                cellptr, cellptr, bodyptr,
                                gdhistptr_sincos_omp_kkk,
                                gdhistptr_sincos_omp_kkk_N);
local void sumcell_balls6_omp_N(struct cmdline_data* cmd, struct  global_data* gd,
                                cellptr, cellptr, bodyptr,
                                gdhistptr_sincos_omp_kkk,
                                gdhistptr_sincos_omp_kkk_N);
local void sumcellcell_balls6_omp_N(struct cmdline_data* cmd,
                                    struct  global_data* gd,
                                    cellptr, cellptr, nodeptr,
                                    gdhistptr_sincos_omp_kkk,
                                    gdhistptr_sincos_omp_kkk_N);

local int search_init_gd_sincos_omp_kkk_N(struct  cmdline_data* cmd,
                                        struct  global_data* gd,
                                          gdlptr_sincos_omp_kkk_N);
local int search_free_gd_sincos_omp_kkk_N(struct  cmdline_data* cmd,
                                         struct  global_data* gd,
                                          gdlptr_sincos_omp_kkk_N gdl);
local int search_init_sincos_omp_kkk_N(struct  cmdline_data* cmd,
                                       struct  global_data* gd,
                                       gdhistptr_sincos_omp_kkk_N, int);
local int search_free_sincos_omp_kkk_N(struct  cmdline_data* cmd,
                                         struct  global_data* gd,
                                     gdhistptr_sincos_omp_kkk_N);
local int computeBodyProperties_sincos_kkk_N(struct  cmdline_data* cmd,
                                            struct  global_data* gd,
                                            bodyptr p, int nbody,
                                       gdhistptr_sincos_omp_kkk_N);
local int computeBodyProperties_sincos_kkk_sum_balls6_omp_N(
                                            struct  cmdline_data*,
                                            struct  global_data*,
                                            bodyptr, INTEGER,
                                            gdhistptr_sincos_omp_kkk_N,
                                            gdhistptr_sincos_omp_kkk_N);
local int search_init_omp_balls6_kkk_N(struct cmdline_data* cmd,
                                         struct  global_data* gd,
                                        gdhistptr_sincos_omp_kkk_N, int);
local int search_free_omp_balls6_kkk_N(struct cmdline_data* cmd,
                                  struct  global_data* gd,
                                       gdhistptr_sincos_omp_kkk_N);
#endif // ! NMultipoles

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
#ifdef TWOPCF
local int search_compute_Xi_kkk_balls4(struct cmdline_data* cmd,
                                      struct global_data* gd,
                                      INTEGER nbody,
                                      gdlptr_sincos_omp_kkk gdl);
local int search_compute_HistN_kkk_balls4(struct cmdline_data* cmd,
                                         struct global_data* gd,
                                         INTEGER nbody,
                                         gdlptr_sincos_omp_kkk gdl);
local int PrintHistNN_kkk_balls4(struct cmdline_data* cmd,
                                struct global_data* gd,
                                gdlptr_sincos_omp_kkk gdl);
local int PrintHistCF_kkk_balls4(struct cmdline_data* cmd,
                                struct global_data* gd,
                                gdlptr_sincos_omp_kkk gdl);
local int PrintHistXi2pcf_kkk_balls4(struct cmdline_data* cmd,
                                    struct global_data* gd,
                                    gdlptr_sincos_omp_kkk gdl);
#endif

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
//local int HistZetaM_sincos_edge_effects(struct  cmdline_data*,
//                                        struct  global_data*,
//                                        gdhistptr_sincos_omp_kkk,
//                                        gdhistptr_sincos_omp_kkk_N);

#endif
#endif // ! NMultipoles
//E
//E Saving histograms section: case KKKCORRELATION:

/*
 Search routine using tree brute force direct method:

 To be called using: search=octree-balls4-omp

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
    *                                    histZetaMsincos, histNN, histCF,
    *                                    histNNSubXi2pcf, histXi2pcf,
    * Counting encounters (in global gd): nbbcalc, nbccalc, ncccalc
 Return (the error status):
    int SUCCESS or FAILURE
 */
static int balls4_reduce_legacy(struct cmdline_data *cmd, struct global_data *gd,
                                gdlptr_sincos_omp_kkk hist
#ifdef NMultipoles
                                , gdlptr_sincos_omp_kkk_N window
#endif
                                )
{
    if (!balls4_distributed(cmd)) return SUCCESS;
    size_t bins = (size_t)cmd->sizeHistN, orders = (size_t)cmd->mChebyshev+1;
    size_t planes = 4;
#ifdef NMultipoles
    planes += 4;
#endif
    if (bins > SIZE_MAX/bins || orders > SIZE_MAX/planes/(bins*bins)
        || orders*planes*bins*bins > SIZE_MAX-4*bins)
        cBALLS_FAIL(cmd, "BALLS4 MPI histogram size overflow\n");
    const size_t count = orders*planes*bins*bins + 4*bins;
    real *packed = count <= SIZE_MAX/sizeof(real) ? calloc(count, sizeof(real)) : NULL;
    if (balls4_consensus(cmd, packed ? SUCCESS : FAILURE,
            "BALLS4 MPI reduction allocation") == FAILURE) {
        free(packed); return FAILURE;
    }
    size_t offset = 0;
#define B4_PACK_VECTOR(v) do { for (size_t n=1; n<=bins; n++) packed[offset++] = (v)[n]; } while (0)
#define B4_PACK_MATRIX(v) do { for (size_t m=1; m<=orders; m++) for (size_t i=1; i<=bins; i++) for (size_t j=1; j<=bins; j++) packed[offset++] = (v)[m][i][j]; } while (0)
#ifdef TWOPCF
    B4_PACK_VECTOR(hist->histNN); B4_PACK_VECTOR(hist->histWW);
    B4_PACK_VECTOR(hist->histNNSubXi2pcf); B4_PACK_VECTOR(hist->histXi2pcf);
#endif
#ifdef THREEPCFCONVERGENCE
    B4_PACK_MATRIX(hist->histZetaMcos); B4_PACK_MATRIX(hist->histZetaMsin);
    B4_PACK_MATRIX(hist->histZetaMsincos); B4_PACK_MATRIX(hist->histZetaMcossin);
#ifdef NMultipoles
    B4_PACK_MATRIX(window->histZetaMcos); B4_PACK_MATRIX(window->histZetaMsin);
    B4_PACK_MATRIX(window->histZetaMsincos); B4_PACK_MATRIX(window->histZetaMcossin);
#endif
#endif
    INTEGER counters[4] = {gd->nbbcalc, gd->nbccalc, gd->ncccalc, gd->nsmoothcount};
    int status = balls4_reduce(cmd, packed, offset);
    if (status == SUCCESS) status = balls4_reduce_counts(cmd, counters, 4);
    if (status == SUCCESS && balls4_publish(cmd)) {
        offset = 0;
#define B4_UNPACK_VECTOR(v) do { for (size_t n=1; n<=bins; n++) (v)[n] = packed[offset++]; } while (0)
#define B4_UNPACK_MATRIX(v) do { for (size_t m=1; m<=orders; m++) for (size_t i=1; i<=bins; i++) for (size_t j=1; j<=bins; j++) (v)[m][i][j] = packed[offset++]; } while (0)
#ifdef TWOPCF
        B4_UNPACK_VECTOR(hist->histNN); B4_UNPACK_VECTOR(hist->histWW);
        B4_UNPACK_VECTOR(hist->histNNSubXi2pcf); B4_UNPACK_VECTOR(hist->histXi2pcf);
#endif
#ifdef THREEPCFCONVERGENCE
        B4_UNPACK_MATRIX(hist->histZetaMcos); B4_UNPACK_MATRIX(hist->histZetaMsin);
        B4_UNPACK_MATRIX(hist->histZetaMsincos); B4_UNPACK_MATRIX(hist->histZetaMcossin);
#ifdef NMultipoles
        B4_UNPACK_MATRIX(window->histZetaMcos); B4_UNPACK_MATRIX(window->histZetaMsin);
        B4_UNPACK_MATRIX(window->histZetaMsincos); B4_UNPACK_MATRIX(window->histZetaMcossin);
#endif
#endif
        gd->nbbcalc = counters[0]; gd->nbccalc = counters[1];
        gd->ncccalc = counters[2]; gd->nsmoothcount = counters[3];
    }
#undef B4_PACK_VECTOR
#undef B4_PACK_MATRIX
#undef B4_UNPACK_VECTOR
#undef B4_UNPACK_MATRIX
    free(packed);
    return status;
}

global int searchcalc_octree_balls4_omp(struct cmdline_data* cmd,
                                             struct  global_data* gd,
                                             bodyptr *btable, INTEGER *nbody,
                                             INTEGER ipmin, INTEGER *ipmax,
                                             int cat1, int cat2)
{
    string routine_name = "searchcalc_octree_balls4_omp";
    double cpustart;
    int result_status = SUCCESS;
    gdl_sincos_omp_kkk gdl;
#ifdef NMultipoles
    gdl_sincos_omp_kkk_N gdlN;
#endif
#ifdef DEBUG
    INTEGER ibfcount=0;
#endif

    cpustart = CPUTIME;
    if (print_info(cmd, gd) == FAILURE)
        return FAILURE;

#ifdef OPENMPCODE
    ThreadCount(cmd, gd, nbody[cat1], cat1);
#else
#error `OPENMPMACHINE` is not defined. Switch it on in Makefile_settings
#endif
    if (cballs_opt_edge_corrections(cmd) || cballs_opt_no_normalize_histzeta(cmd))
        return balls4_edge_search(cmd, gd, btable, nbody, ipmin, ipmax, cat1, cat2);

//    search_init_gd_sincos_omp_kkk(cmd, gd);
//#ifdef NMultipoles
//    search_init_gd_sincos_omp_kkk_N(cmd, gd);
//#endif

    search_init_gd_sincos_omp_kkk(cmd, gd, &gdl);
#ifdef NMultipoles
    search_init_gd_sincos_omp_kkk_N(cmd, gd, &gdlN);
#endif

#ifdef DEBUG
    bodytabbf_balls4 = (bodyptr) allocate(nbody[cat1] * sizeof(body));
    if (bodytabbf_balls4 == NULL) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: unable to allocate the debug pivot table\n",
                 routine_name);
        return FAILURE;
    }
    gd->bytes_tot += nbody[cat1]*sizeof(body);
    fprintf(stdout,"\n\nAllocated %g MByte for particle (found) storage.\n",
            nbody[cat1]*sizeof(body)/(1024.0*1024.0));

    verb_print(cmd->verbose,
            "\n%s: Total allocated %g MByte storage so far.\n",
               routine_name, gd->bytes_tot/(1024.0*1024.0));
#endif

    verb_print(cmd->verbose,
        "\n%s: Total allocated %g MByte storage so far.\n",
               routine_name, gd->bytes_tot/(1024.0*1024.0));

#ifdef DEBUG
#pragma omp parallel default(none)                                          \
    shared(cmd,gd,btable,nbody,roottable,                                   \
           ibfcount,                                                     \
           ipmin,ipmax,cat1,cat2,nodetablescanlevB4,gdl)
#else // ! DEBUG

#ifdef NMultipoles
#pragma omp parallel default(none)                                          \
    shared(cmd,gd,btable,nbody,roottable,                                   \
           ipmin,ipmax,cat1,cat2,nodetablescanlevB4,gdl,gdlN)
#else
#pragma omp parallel default(none)                                          \
    shared(cmd,gd,btable,nbody,roottable,                                   \
           ipmin,ipmax,cat1,cat2,nodetablescanlevB4,gdl)
#endif

#endif // ! DEBUG
  {
    nodeptr p;
    int n;
    int m;
    INTEGER i;

    //B init:
    gdhist_omp_balls6_kkk histb;
    gdhist_sincos_omp_kkk hist;
    search_init_omp_balls6_kkk(cmd, gd, &histb, cat1);
    search_init_sincos_omp_kkk(cmd, gd, &hist, cat1);
#ifdef NMultipoles
    gdhist_sincos_omp_kkk_N histN;
    search_init_sincos_omp_kkk_N(cmd, gd, &histN, cat1);
#endif
    hist.thrID = omp_get_thread_num();
    //E

#ifndef PRAGMABCINTLIST
#pragma omp for nowait schedule(dynamic)
#endif
      for (i=0; i<gd->nnodescanlevTableB4[cat1]; i++) {
          if (!balls4_task_owned(cmd, i)) continue;
          p = nodetablescanlevB4[cat1][i];
          if (cballs_opt_read_mask(cmd)
              && !mask_node_can_approximate(Mask(p)))
              continue;

          //B Set histograms to zero for the pivot
#ifdef THREEPCFCONVERGENCE
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
#endif
          //E
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
          //E
          histb.active[0] = (nodeptr) (roottable[cat2]);
#ifdef NMultipoles
          walktree_balls6_omp_N(cmd, gd, histb.active, histb.active + 1,
                  histb.interact, histb.interact + histb.actlen,
                  p, Size(p), Pos(p), &histb, &hist, &histN, nbody[cat1], cat1);
#ifdef THREEPCFCONVERGENCE
          if (!cballs_opt_no_two_balls(cmd)) {
              computeBodyProperties_sincos_kkk(cmd, gd, (bodyptr)p,
                                                  nbody[cat1], &hist);
              computeBodyProperties_sincos_kkk_N(cmd, gd, (bodyptr)p,
                                                  nbody[cat1], &histN);
          }
#endif
#else
          walktree_balls6_omp(cmd, gd, histb.active, histb.active + 1,
                  histb.interact, histb.interact + histb.actlen,
                  p, Size(p), Pos(p), &histb, &hist, nbody[cat1], cat1);
#ifdef THREEPCFCONVERGENCE
          if (!cballs_opt_no_two_balls(cmd)) {
              computeBodyProperties_sincos_kkk(cmd, gd, (bodyptr)p,
                                                  nbody[cat1], &hist);
          }
#endif
#endif
      } // end do body i

#pragma omp critical
    {
#ifdef TWOPCF
        for (n = 1; n <= cmd->sizeHistN; n++) {
            gdl.histNN[n] += hist.histNthread[n];
            gdl.histWW[n] += hist.histWWthread[n];
            gdl.histNNSubXi2pcf[n] += hist.histNNSubXi2pcfthread[n];
            gdl.histXi2pcf[n] += hist.histXi2pcfthread[n];
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
        }
#endif
        gd->nsmoothcount += hist.nsmoothcountthread;
        gd->nbbcalc += hist.nbbcalcthread;
        gd->nbccalc += hist.nbccalcthread;
        gd->ncccalc += hist.ncccalcthread;
#ifdef DEBUG
        ibfcount += hist.ibfcountthread;
#endif
#ifdef NMultipoles
#ifdef THREEPCFCONVERGENCE
/*        for (m=1; m<=cmd->mChebyshev+1; m++) {
            ADDM_ext(gd->NhistZetaMcos[m],gd->NhistZetaMcos[m],
                     histN.histZetaMthreadcos[m],cmd->sizeHistN);
            ADDM_ext(gd->NhistZetaMsin[m],gd->NhistZetaMsin[m],
                     histN.histZetaMthreadsin[m],cmd->sizeHistN);
            ADDM_ext(gd->NhistZetaMsincos[m],gd->NhistZetaMsincos[m],
                     histN.histZetaMthreadsincos[m],cmd->sizeHistN);
            ADDM_ext(gd->NhistZetaMcossin[m],gd->NhistZetaMcossin[m],
                     histN.histZetaMthreadcossin[m],cmd->sizeHistN);
        } */
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
#endif
    } // ! critical

#ifdef NMultipoles
    search_free_sincos_omp_kkk_N(cmd, gd, &histN);  // free memory
#endif
    search_free_sincos_omp_kkk(cmd, gd, &hist);     // free memory
    search_free_omp_balls6_kkk(cmd, gd, &histb);

  } // end pragma omp parallel

    result_status = balls4_reduce_legacy(cmd, gd, &gdl
#ifdef NMultipoles
                                         , &gdlN
#endif
                                         );
    if (result_status == FAILURE || !balls4_publish(cmd))
        goto balls4_cleanup;

//#ifdef DEBUG
//    gd->nbodybf = cmd->ntosave;
//    verb_print(cmd->verbose, "\noctree-balls4-omp: Total bodies found: %ld %ld\n",ibfcount, gd->nbodybf);
//#endif

    verb_print(cmd->verbose,
               "octree-balls4-omp: nsmoothcount = %ld\n",gd->nsmoothcount);

#ifdef TWOPCF
    {
        int nn;
        bool symmetric = !cballs_opt_asymmetric(cmd);

        for (nn = 1; nn <= cmd->sizeHistN; nn++) {
            if (cballs_opt_weights_norm(cmd)) {
                gdl.histXi2pcf[nn] = cballs_normalize_or_zero(
                    gdl.histXi2pcf[nn], gdl.histWW[nn]);
            } else {
                real numerator = gdl.histXi2pcf[nn];
                real denominator = gdl.histNNSubXi2pcf[nn];

                if (symmetric) {
                    numerator *= 0.5;
                    denominator *= 0.5;
                }
                gdl.histXi2pcf[nn] =
                    cballs_normalize_or_zero(numerator, denominator);
                gdl.histNNSubXi2pcf[nn] = denominator;
            }
        }

        if (cballs_opt_compute_histn(cmd)) {
            INTEGER valid_bodies = nbody[cat1];
            if (cballs_opt_read_mask(cmd)) {
                valid_bodies = 0;
                for (INTEGER b=0; b<nbody[cat1]; b++)
                    if (Mask(btable[cat1]+b) != MASK_NODE_MASKED) valid_bodies++;
            }
            if (search_compute_HistN_kkk_balls4(
                    cmd, gd, valid_bodies, &gdl) == FAILURE) {
#ifdef NMultipoles
                search_free_gd_sincos_omp_kkk_N(cmd, gd, &gdlN);
#endif
                search_free_gd_sincos_omp_kkk(cmd, gd, &gdl);
                return FAILURE;
            }
        }
    }
#endif

// ===============================================
//B Saving histograms section: case KKKCORRELATION:
// ===============================================
    if (!cballs_opt_no_out_hist(cmd)) {
    verb_print(cmd->verbose,
            "\n\tsearch_octree_kkk_omp: printing octree-kkk-omp method...\n\n");
#ifdef TWOPCF
    if (cballs_opt_compute_histn(cmd))
        PRINT_OR_FAIL(PrintHistNN_kkk_balls4(cmd, gd, &gdl));
    PRINT_OR_FAIL(PrintHistXi2pcf_kkk_balls4(cmd, gd, &gdl));
#endif
    PrintHistrBins(cmd, gd);
#ifdef THREEPCFCONVERGENCE
#ifdef NMultipoles

#ifdef NONORMHIST
    if (cballs_opt_no_normalize_histzeta(cmd))
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

    if (cballs_opt_out_m_histzeta(cmd)) {
#ifdef NMultipoles

#ifdef NONORMHIST
        if (cballs_opt_no_normalize_histzeta(cmd))
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
    }

    if (cballs_opt_out_histzetag(cmd)) {
        PrintHistZetaGm_sincos(cmd, gd, &gdl);
        PrintHistZetaG(cmd, gd, &gdl);
        PrintHistZetaMZetaGm_sincos(cmd, gd, &gdl);
    }

#ifdef NMultipoles
#ifdef NONORMHIST
#endif // ! NONORMHIST
#endif // ! NMultipoles
#endif // ! THREEPCFCONVERGENCE
    }
    gd->flagPrint = FALSE;
// ===============================================
//E Saving histograms section: case KKKCORRELATION
// ===============================================

#ifdef PXD
    {
        int publish_m, publish_n1, publish_n2;
#ifdef TWOPCF
        for (publish_n1=1; publish_n1<=cmd->sizeHistN; publish_n1++) {
            gd->histNN[publish_n1] = gdl.histNN[publish_n1];
            gd->histCF[publish_n1] = gdl.histCF[publish_n1];
            gd->histXi2pcf[publish_n1] = gdl.histXi2pcf[publish_n1];
        }
#endif
#ifdef THREEPCFCONVERGENCE
        for (publish_m=1; publish_m<=cmd->mChebyshev+1; publish_m++) {
            for (publish_n1=1; publish_n1<=cmd->sizeHistN; publish_n1++) {
                for (publish_n2=1; publish_n2<=cmd->sizeHistN; publish_n2++) {
                    real scale = 1.0;
#if defined(NMultipoles) && defined(NONORMHIST)
                    if (!cballs_opt_no_normalize_histzeta(cmd)) {
                        real denominator =
                            gdlN.histZetaMcos[1][publish_n1][publish_n2];
                        scale = denominator != 0.0 ? 1.0/denominator : 0.0;
                    }
#endif
                    gd->histZetaMcos[publish_m][publish_n1][publish_n2] =
                        scale*gdl.histZetaMcos[publish_m][publish_n1][publish_n2];
                    gd->histZetaMsin[publish_m][publish_n1][publish_n2] =
                        scale*gdl.histZetaMsin[publish_m][publish_n1][publish_n2];
                    gd->histZetaMsincos[publish_m][publish_n1][publish_n2] =
                        scale*gdl.histZetaMsincos[publish_m][publish_n1][publish_n2];
                    gd->histZetaMcossin[publish_m][publish_n1][publish_n2] =
                        scale*gdl.histZetaMcossin[publish_m][publish_n1][publish_n2];
                }
            }
        }
#endif
    }
#endif

balls4_cleanup:
    gd->flagPrint = FALSE;
#ifdef NMultipoles
    search_free_gd_sincos_omp_kkk_N(cmd, gd, &gdlN);// free memory
#endif
    search_free_gd_sincos_omp_kkk(cmd, gd, &gdl); // free memory

#ifdef DEBUG
    free(bodytabbf_balls4);
    bodytabbf_balls4 = NULL;
#endif

    gd->cpusearch = CPUTIME - cpustart;
    verb_print(cmd->verbose, "\nGoing out: CPU time = %lf %s\n",
               CPUTIME-cpustart, PRNUNITOFTIMEUSED);

    return result_status;
}

//#ifndef BALLS
//B 2023.11.22
local bool nodes_condition_balls(struct cmdline_data* cmd,
                                 struct  global_data* gd,
                                 nodeptr p, nodeptr q, real *dr1,
                                 compute_vector dr)
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

     if (cballs_opt_read_mask(cmd)
         && (!mask_node_can_approximate(Mask(p))
             || !mask_node_can_approximate(Mask(q))))
         return FALSE;

     if ( *dr1 == 0.0)
         return (FALSE);
     else
         if ( (Radius(p)+Radius(q))/(*dr1) < gd->deltaR) {
             if (cballs_opt_behavior_tree_omp(cmd)) {
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

local bool collapse_small_cell_balls4(struct cmdline_data* cmd,
                                      struct global_data* gd,
                                      nodeptr q, real distance)
{
    if (balls4_cell_requires_open(cmd, q))
        return FALSE;
#ifdef OCTREESMOOTHING
    return ((Nb(q) <= gd->nsmooth[0] || Size(q) <= gd->rminCell[0])
            && distance > gd->rminCell[1]);
#else
    (void)gd;
    (void)distance;
    return FALSE;
#endif
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
    compute_vector dr;

    if (balls4_node_is_masked(cmd, p))
        return;

    if (Update(p)) {
        np = nptr;
        actsafe = histb->actlen - NSUB;
        for (ap = aptr; ap < nptr; ap++) {
            if (balls4_node_is_masked(cmd, *ap))
                continue;
            if (Type(*ap) == CELL) {
                if (!reject_cell_balls(cmd, gd, p, *ap, &dr1, dr)) {
                    if (balls4_cell_requires_open(cmd, *ap)) {
                        if (np - histb->active >= actsafe)
                            error("walktree (mask): active list overflow\n");
                        for (q = More(*ap); q != Next(*ap); q = Next(q))
                            *np++ = q;
                        continue;
                    }
                    if (collapse_small_cell_balls4(cmd, gd, *ap, dr1)) {
                        if (np - histb->active >= actsafe) {
                            error("walktree (1): active list overflow\n");
                        }
                        histsincos->nsmoothcountthread += 1;
                        --bptr;
                        Mass(bptr) = Mass(*ap);
                        Weight(bptr) = Weight(*ap);
                        Kappa(bptr) = Kappa(*ap);
                        SETV(Pos(bptr), Pos(*ap));
                        Id(bptr) = Id(*ap);
                        Type(bptr) = Type(*ap);
                        Nb(bptr) = Nb(*ap);
                     } else { // ! bucket condition
                         if (nodes_condition_balls(cmd, gd, p, *ap, &dr1, dr)) {
                             if (!cballs_opt_no_two_balls(cmd)
                                 && Type(p) == CELL ) {
                                 sumcellcell_balls6_omp(cmd, gd, (cellptr)(*ap),
                                        (cellptr)*ap+1, p, histsincos);
                             } else {
                                 if (np - histb->active >= actsafe)
                                     error("walktree (2): active list overflow\n");
                                 if (!cballs_opt_no_one_ball(cmd)) {
                                     Mass(cptr) = Mass(*ap);
                                     Weight(cptr) = Weight(*ap);
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
                     Weight(bptr) = Weight(*ap);
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
            bodyptr pbf = bodytabbf_balls4 + histsincos->ibfcountthread;
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
            if (balls4_node_is_masked(cmd, q))
                continue;
            for (k = 0; k < NDIM; k++)
                nmid[k] = pmid[k] + (Pos(q)[k] < pmid[k] ? - poff : poff);
            walktree_balls6_omp(cmd, gd, nptr, np, cptr, bptr, q, psize / 2, nmid,
            histb, histsincos, nbody, ifile);
        }
    } else {
        if (balls4_node_is_masked(cmd, p))
            return;
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
    int n;
    INTEGER ip;
    gdhist_sincos_omp_kkk hist1sincos;


    //B init:
    search_init_sincos_omp_kkk(cmd, gd, &hist1sincos, ifile);
    for (n = 1; n <= cmd->sizeHistN; n++) {
        hist1sincos.histNNSubthread[n] = 0.0;
    }
    for (n = 1; n <= cmd->sizeHistN; n++) {
        hist1sincos.histNNSubthread[n] = 0.0;
    }
    CLRM_ext_ext(hist1sincos.histXithreadcos, cmd->mChebyshev+1,
                 cmd->sizeHistN);
    CLRM_ext_ext(hist1sincos.histXithreadsin, cmd->mChebyshev+1,
                 cmd->sizeHistN);
    //E

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

    if (!cballs_opt_no_one_ball(cmd))
        sumcell_balls6_omp(cmd, gd, histb->interact, cptr, (bodyptr) p,
                           &hist1sincos);
    sumnode_balls6_omp(cmd, gd, bptr, histb->interact + histb->actlen,
                       (bodyptr) p, &hist1sincos);

    computeBodyProperties_sincos_kkk_sum_balls6_omp(cmd, gd,
                                    p, nbody, histsincos, &hist1sincos);

    histsincos->nbbcalcthread += histb->interact + histb->actlen - bptr;
    histsincos->nbccalcthread += cptr - histb->interact;

    ip = p - bodytable[gd->iCatalogs[0]] + 1;
//    if (histsincos->thrID==0)
    if (ip%cmd->stepState == 0) {
        verb_print(cmd->verbose_log, " - Completed pivot: %ld\n", ip);
    }

    search_free_sincos_omp_kkk(cmd, gd, &hist1sincos);
    return SUCCESS;
}

local void sumnode_balls6_omp(struct cmdline_data* cmd, struct  global_data* gd,
                              cellptr start, cellptr finish, bodyptr p0,
                              gdhistptr_sincos_omp_kkk hist)
{
    cellptr p;
    compute_vector dr;
    real dr1;
    nodeptr pb;
    INTEGER ibodycount=0;
    int n;
    real xi;
    real cosphi;
    real sinphi;

#ifdef PRAGMABCINTLIST
#pragma omp for nowait schedule(dynamic)
#endif
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
#ifdef TWOPCF
                    balls4_accumulate_neighbor_2pcf(hist, pb, n);
#endif
                    xi = Kappa(pb);
                    //B Component of pb with respect to the axis of reference
                    if (!cballs_angular_phase(Pos(p0), dr, &cosphi, &sinphi))
                        continue;
                    CHEBYSHEVTUOMPKKK;
                } // ! 1 < n < HistN
            } // ! dr1 > rmin
        } // ! accept_body
    } // ! loop p
}

local void sumcell_balls6_omp(struct cmdline_data* cmd, struct  global_data* gd,
                              cellptr start, cellptr finish, bodyptr p0,
                              gdhistptr_sincos_omp_kkk hist)
{
    cellptr p;
    compute_vector dr;
    real dr1;
    nodeptr pb;
    INTEGER ibodycount=0;
    int n;
    real xi;
    real cosphi;
    real sinphi;

#ifdef PRAGMABCINTLIST
#pragma omp for nowait schedule(dynamic)
#endif
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
#ifdef TWOPCF
                    balls4_accumulate_neighbor_2pcf(hist, pb, n);
#endif
                    //B correction 2025-04-06
                    // needs to multiply xi by Weight(p) of the cell q
#ifdef NONORMHIST
                    if (cballs_opt_no_normalize_histzeta(cmd)) {
                        xi = Nb(pb)*Kappa(pb);
                    } else {
                        xi = Kappa(pb);
                        }
#else
                    xi = Kappa(pb);
#endif
                    //E
//                    xi = Kappa(pb);
                    //B Component of pb with respect to the axis of reference
                    if (!cballs_angular_phase(Pos(p0), dr, &cosphi, &sinphi))
                        continue;
                    CHEBYSHEVTUOMPKKK;
                } // ! 1 < n < HistN
            } // ! dr1 > rmin
        } // ! accept_body
    } // ! loop p
}

local void sumcellcell_balls6_omp(struct cmdline_data* cmd,
                                  struct  global_data* gd,
                                  cellptr start, cellptr finish, nodeptr p0,
//                                  gdhistptr_omp_balls6_kkk histb,
                                  gdhistptr_sincos_omp_kkk hist)
{
    cellptr p;
    compute_vector dr;
    real dr1;
    nodeptr pb;
    INTEGER ibodycount=0;
    int n;
    real xi;
    real cosphi;
    real sinphi;

//    verb_print_debug(1, "\nAqui voy (0): inside cell-cell counting");

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
#ifdef TWOPCF
                    balls4_accumulate_cell_pair_2pcf(hist, p0, pb, n);
#endif
                    //B correction 2025-04-06
                    // needs to multiply xi by Weight(p) of the cell q
#ifdef NONORMHIST
                    if (cballs_opt_no_normalize_histzeta(cmd)) {
                        xi = Nb(pb)*Kappa(pb);
                    } else {
                        xi = Kappa(pb);
                        }
#else
                    xi = Kappa(pb);
#endif
                    //E
//                    xi = Kappa(pb);
                    //B Component of pb with respect to the axis of reference
                    if (!cballs_angular_phase(Pos(p0), dr, &cosphi, &sinphi))
                        continue;
                    CHEBYSHEVTUOMPKKK;
                } // ! 1 < n < HistN
            } // ! dr1 > rmin
        } // ! accept_body
    } // ! loop p

    hist->ncccalcthread += 1;
//    verb_print_debug(1, "Aqui voy (1)\n");
}


//B Routines as in cballsutils

local int search_init_gd_sincos_omp_kkk(struct  cmdline_data* cmd,
                                        struct  global_data* gd,
                                        gdlptr_sincos_omp_kkk gdl)
{
    int n;
    int m;
    INTEGER bytes_tot_local=0;

#ifdef TWOPCF
    gdl->histNN = dvector(1,cmd->sizeHistN);
    gdl->histWW = dvector(1,cmd->sizeHistN);
    gdl->histCF = dvector(1,cmd->sizeHistN);
    gdl->histNNSubXi2pcf = dvector(1,cmd->sizeHistN);
    gdl->histXi2pcf = dvector(1,cmd->sizeHistN);
    bytes_tot_local += 5*cmd->sizeHistN*sizeof(real);
#endif

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
    gdl->histXi3pcf = dmatrix3D(1,cmd->sizeHistPhi,
                                1,cmd->sizeHistN,1,cmd->sizeHistN);
    bytes_tot_local +=
                (cmd->sizeHistN*cmd->sizeHistN*cmd->sizeHistPhi)*sizeof(real);

    gd->bytes_tot += bytes_tot_local;
    verb_print(cmd->verbose,
    "\nsearch_init_gd_octree_kkk: Allocated %g MByte for histograms storage.\n",
    bytes_tot_local*INMB);

#ifdef TWOPCF
    for (n = 1; n <= cmd->sizeHistN; n++) {
        gdl->histNN[n] = 0.0;
        gdl->histWW[n] = 0.0;
        gdl->histCF[n] = 0.0;
        gdl->histNNSubXi2pcf[n] = 0.0;
        gdl->histXi2pcf[n] = 0.0;
    }
#endif
    for (n = 1; n <= cmd->sizeHistN; n++)
        gdl->histNNSub[n] = 0.0;
    for (m = 1; m <= cmd->mChebyshev+1; m++) {
        CLRM_ext(gdl->histZetaMcos[m], cmd->sizeHistN);
        CLRM_ext(gdl->histZetaMsin[m], cmd->sizeHistN);
        CLRM_ext(gdl->histZetaMsincos[m], cmd->sizeHistN);
        CLRM_ext(gdl->histZetaMcossin[m], cmd->sizeHistN);
//        gd->histXi[m][n] = 0.0;
    }
    gd->nbbcalc = gd->nbccalc = gd->ncccalc = gd->nsmoothcount = 0;

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

#ifdef TWOPCF
    free_dvector(gdl->histXi2pcf,1,cmd->sizeHistN);
    free_dvector(gdl->histNNSubXi2pcf,1,cmd->sizeHistN);
    free_dvector(gdl->histCF,1,cmd->sizeHistN);
    free_dvector(gdl->histWW,1,cmd->sizeHistN);
    free_dvector(gdl->histNN,1,cmd->sizeHistN);
#endif

    return SUCCESS;
}

local int search_init_sincos_omp_kkk(struct  cmdline_data* cmd,
                                     struct  global_data* gd,
                                     gdhistptr_sincos_omp_kkk hist, int ifile)
{
    int n;
    int m;

#ifdef TWOPCF
    hist->histNthread = dvector(1,cmd->sizeHistN);
    hist->histWthread = dvector(1,cmd->sizeHistN);
    hist->histWWthread = dvector(1,cmd->sizeHistN);
    hist->histNNSubXi2pcfthread = dvector(1,cmd->sizeHistN);
    hist->histXi2pcfthread = dvector(1,cmd->sizeHistN);
    hist->histXi2pcfthreadsub = dvector(1,cmd->sizeHistN);
#endif

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

    for (n = 1; n <= cmd->sizeHistN; n++) {
        hist->histNNSubthread[n] = 0.0;
#ifdef TWOPCF
        hist->histNthread[n] = 0.0;
        hist->histWthread[n] = 0.0;
        hist->histWWthread[n] = 0.0;
        hist->histNNSubXi2pcfthread[n] = 0.0;
        hist->histXi2pcfthread[n] = 0.0;
        hist->histXi2pcfthreadsub[n] = 0.0;
#endif
    }
    for (m = 1; m <= cmd->mChebyshev+1; m++) {
        CLRM_ext(hist->histZetaMthreadcos[m], cmd->sizeHistN);
        CLRM_ext(hist->histZetaMthreadsin[m], cmd->sizeHistN);
        CLRM_ext(hist->histZetaMthreadsincos[m], cmd->sizeHistN);
        CLRM_ext(hist->histZetaMthreadcossin[m], cmd->sizeHistN);
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

#ifdef TWOPCF
    free_dvector(hist->histXi2pcfthreadsub,1,cmd->sizeHistN);
    free_dvector(hist->histXi2pcfthread,1,cmd->sizeHistN);
    free_dvector(hist->histNNSubXi2pcfthread,1,cmd->sizeHistN);
    free_dvector(hist->histWWthread,1,cmd->sizeHistN);
    free_dvector(hist->histWthread,1,cmd->sizeHistN);
    free_dvector(hist->histNthread,1,cmd->sizeHistN);
#endif

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

#ifdef NONORMHIST
    xi = Weight(p)*Kappa(p);
#else
    xi = Weight(p)*Kappa(p)/nbody;
#endif

#ifndef NONORMHIST
    for (m=1; m<=cmd->mChebyshev+1; m++)
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
#ifdef TWOPCF
    real pivot_weight_sum;
    real pivot_field_sum;
#endif

#ifdef NONORMHIST
    xi = Weight(p)*Kappa(p);
#else
    xi = Weight(p)*Kappa(p)/nbody;
#endif

#ifdef TWOPCF
    pivot_weight_sum = balls4_node_weight_sum((nodeptr)p);
    pivot_field_sum = balls4_node_field_sum((nodeptr)p);
    for (n = 1; n <= cmd->sizeHistN; n++) {
        hist1sincos->histXi2pcfthread[n] +=
            pivot_field_sum*hist1sincos->histXi2pcfthreadsub[n];
        hist1sincos->histWWthread[n] +=
            pivot_weight_sum*hist1sincos->histWthread[n];
    }
#endif

    for (m=1; m<=cmd->mChebyshev+1; m++) {
#ifndef NONORMHIST
        for (n=1; n<=cmd->sizeHistN; n++) {
            hist1sincos->histXithreadcos[m][n] /=
            MAX(hist1sincos->histNNSubthread[n],1.0);
            hist1sincos->histXithreadsin[m][n] /=
            MAX(hist1sincos->histNNSubthread[n],1.0);
        }
#endif
    }

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

    for (n = 1; n <= cmd->sizeHistN; n++) {
        histsincos->histNNSubthread[n] += hist1sincos->histNNSubthread[n];
#ifdef TWOPCF
        histsincos->histNthread[n] += hist1sincos->histNthread[n];
        histsincos->histWWthread[n] += hist1sincos->histWWthread[n];
        histsincos->histNNSubXi2pcfthread[n] +=
            hist1sincos->histNNSubXi2pcfthread[n];
        histsincos->histXi2pcfthread[n] +=
            hist1sincos->histXi2pcfthread[n];
#endif
    }

    return SUCCESS;
}

//E Routines as in cballsutils


local int search_init_omp_balls6_kkk(struct cmdline_data* cmd,
                                 struct  global_data* gd,
                                 gdhistptr_omp_balls6_kkk hist, int ifile)
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
        cBALLS_FAIL(cmd,
            "CheckParameters: can't have periodic boundaries and "
            "OCTREEBALLS4OMP (usePeriodic=%d); set usePeriodic=false\n",
            cmd->usePeriodic);
    if (cmd->useLogHist==FALSE)
        cBALLS_FAIL(cmd,
            "CheckParameters: can't have normal-scale histograms and "
            "OCTREEBALLS4OMP (useLogHist=%d); set useLogHist=true\n",
            cmd->useLogHist);
#if !defined(THREEPCFCONVERGENCE) && !defined(TWOPCF)
    cBALLS_FAIL(cmd,
        "CheckParameters: OCTREEBALLS4OMP requires TWOPCF or TPCF support\n");
#endif
#if NDIM == 2
    cBALLS_FAIL(cmd,
        "CheckParameters: OCTREEBALLS4OMP works only in 3D\n");
#endif
#ifdef NMultipoles
    verb_print(cmd->verbose, "with NMultipoles... \n");
#else
    verb_print(cmd->verbose, "without NMultipoles... \n");
#endif
#ifdef TWOPCF
    verb_print(cmd->verbose, "with 2pcf convergence computation... \n");
#endif
#ifdef THREEPCFCONVERGENCE
    verb_print(cmd->verbose, "with 3pcf convergence computation... \n");
#endif
#ifdef NONORMHIST
    verb_print(cmd->verbose, "with NONORMHIST... \n");
    if (cballs_opt_no_normalize_histzeta(cmd))
        verb_print(cmd->verbose, "with option no-normalize-HistZeta...\n");
#else
    verb_print(cmd->verbose, "without NONORMHIST... \n");
#endif
#ifdef POLARAXIS
    verb_print(cmd->verbose, "with POLARAXIS... \n");
#endif
#ifdef NOLIMBER
    verb_print(cmd->verbose, "with NOLIMBER (no Limber aproximation)... \n");
#endif
    if (cballs_opt_no_one_ball(cmd))
        verb_print(cmd->verbose, "with option no-one-ball... \n");
    if (cballs_opt_no_two_balls(cmd))
        verb_print(cmd->verbose, "with option no-two-balls... \n");
    if (cballs_opt_behavior_ball(cmd))
        verb_print(cmd->verbose, "with option behavior-ball... \n");
    if (cballs_opt_smooth_pivot(cmd))
        cBALLS_FAIL(cmd, "octree-balls4-omp does not support options=smooth-pivot\n");
    if (cballs_opt_bh86(cmd))
        verb_print(cmd->verbose, "with cell opening criterion bh86... \n");
    if (cballs_opt_sw94(cmd))
        verb_print(cmd->verbose, "with cell opening criterion sw94... \n");

    verb_print(cmd->verbose, "with the dedicated BALLS4 scan partition...\n");

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
//        gd->NhistXi[m][n] = 0.0;
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
                                     gdhistptr_sincos_omp_kkk_N hist, int ifile)
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
    //B kappa Avg Rmin
    //E
#endif

    for (m=1; m<=cmd->mChebyshev+1; m++) {
#ifndef NONORMHIST
        for (n=1; n<=cmd->sizeHistN; n++) {
            hist->histXithreadcos[m][n] /= MAX(hist->histNNSubthread[n],1.0);
            hist->histXithreadsin[m][n] /= MAX(hist->histNNSubthread[n],1.0);
        }
#endif
    }
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
#endif // ! NMultipoles


#ifdef NMultipoles

local void walktree_balls6_omp_N(struct cmdline_data* cmd,
                                 struct  global_data* gd,
                                 nodeptr *aptr, nodeptr *nptr,
                                 cellptr cptr, cellptr bptr,
                                 nodeptr p, real psize, vector pmid,
                                 gdhistptr_omp_balls6_kkk histb,
                                 gdhistptr_sincos_omp_kkk histsincos,
                                 gdhistptr_sincos_omp_kkk_N histsincos_N,
                                 INTEGER nbody, int ifile)
{
    nodeptr *np, *ap, q;
    int actsafe;
    real dr1;
    compute_vector dr;

    if (balls4_node_is_masked(cmd, p))
        return;

    if (Update(p)) {
        np = nptr;
        actsafe = histb->actlen - NSUB;
        for (ap = aptr; ap < nptr; ap++) {
            if (balls4_node_is_masked(cmd, *ap))
                continue;
            if (Type(*ap) == CELL) {
                if (!reject_cell_balls(cmd, gd, p, *ap, &dr1, dr)) {
                    if (balls4_cell_requires_open(cmd, *ap)) {
                        if (np - histb->active >= actsafe)
                            error("walktree (mask): active list overflow\n");
                        for (q = More(*ap); q != Next(*ap); q = Next(q))
                            *np++ = q;
                        continue;
                    }
                    if (collapse_small_cell_balls4(cmd, gd, *ap, dr1)) {
                        if (np - histb->active >= actsafe) {
                            error("walktree (1): active list overflow\n");
                        }
                        histsincos->nsmoothcountthread += 1;
                        --bptr;
                        Mass(bptr) = Mass(*ap);
                        Weight(bptr) = Weight(*ap);
                        Kappa(bptr) = Kappa(*ap);
                        SETV(Pos(bptr), Pos(*ap));
                        Id(bptr) = Id(*ap);
                        Type(bptr) = Type(*ap);
                        Nb(bptr) = Nb(*ap);
                     } else { // ! bucket condition
                         if (nodes_condition_balls(cmd, gd, p, *ap, &dr1, dr)) {
                             if (!cballs_opt_no_two_balls(cmd)
                                 && Type(p) == CELL ) {
                                 sumcellcell_balls6_omp_N(cmd, gd, (cellptr)(*ap),
                                        (cellptr)*ap+1, p, histsincos, histsincos_N);
                             } else {
                                 if (np - histb->active >= actsafe)
                                     error("walktree (2): active list overflow\n");
                                 if (!cballs_opt_no_one_ball(cmd)) {
                                     Mass(cptr) = Mass(*ap);
                                     Weight(cptr) = Weight(*ap);
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
                     Weight(bptr) = Weight(*ap);
                     Kappa(bptr) = Kappa(*ap);
                     SETV(Pos(bptr), Pos(*ap));
                     Id(bptr) = Id(*ap);
                     Type(bptr) = Type(*ap);
                     Nb(bptr) = 1;
                 }
        } // ! loop for ap

        gd->actmax = MAX(gd->actmax, np - histb->active);
        if (np != nptr)
            walksub6_omp_N(cmd, gd, nptr, np, cptr, bptr, p, psize, pmid, histb,
                         histsincos, histsincos_N, nbody, ifile);
        else {
            if (Type(p) != BODY)
                error("walktree: recursion terminated with cell\n");

            sum_balls6_omp_N(cmd, gd, cptr, bptr, (bodyptr) p,
                             histb, histsincos, histsincos_N, nbody, ifile);
            Update(p) = FALSE;

#ifdef DEBUG
            bodyptr pbf = bodytabbf_balls4 + histsincos->ibfcountthread;
            Mass(pbf) = Mass(p);
            Kappa(pbf) = Kappa(p);
            SETV(Pos(pbf), Pos(p));
            histsincos->ibfcountthread += 1;
            Id(pbf) = histsincos->ibfcountthread;
#endif
        }
    }   // ! update
}

local void walksub6_omp_N(struct cmdline_data* cmd,
                          struct  global_data* gd,
                          nodeptr *nptr, nodeptr *np, cellptr cptr, cellptr bptr,
                          nodeptr p, real psize, vector pmid,
                          gdhistptr_omp_balls6_kkk histb,
                          gdhistptr_sincos_omp_kkk histsincos,
                          gdhistptr_sincos_omp_kkk_N histsincos_N,
                          INTEGER nbody, int ifile)
{
    real poff;
    nodeptr q;
    int k;
    vector nmid;

    poff = psize / 4;
    if (Type(p) == CELL) {
        for (q = More(p); q != Next(p); q = Next(q)) {
            if (balls4_node_is_masked(cmd, q))
                continue;
            for (k = 0; k < NDIM; k++)
                nmid[k] = pmid[k] + (Pos(q)[k] < pmid[k] ? - poff : poff);
            walktree_balls6_omp_N(cmd, gd,
                                  nptr, np, cptr, bptr, q, psize / 2, nmid,
                                  histb, histsincos, histsincos_N, nbody, ifile);
        }
    } else {
        if (balls4_node_is_masked(cmd, p))
            return;
        for (k = 0; k < NDIM; k++)
            nmid[k] = pmid[k] + (Pos(p)[k] < pmid[k] ? - poff : poff);
        walktree_balls6_omp_N(cmd, gd,
                              nptr, np, cptr, bptr, p, psize / 2, nmid,
                              histb, histsincos, histsincos_N, nbody, ifile);
    }
}


local int sum_balls6_omp_N(struct cmdline_data* cmd,
                         struct  global_data* gd,
                         cellptr cptr, cellptr bptr, bodyptr p,
                         gdhistptr_omp_balls6_kkk histb,
                         gdhistptr_sincos_omp_kkk histsincos,
                         gdhistptr_sincos_omp_kkk_N histsincos_N,
                         INTEGER nbody, int ifile)
{
    int n;
    INTEGER ip;
    gdhist_sincos_omp_kkk hist1sincos;
    gdhist_sincos_omp_kkk_N hist1sincos_N;

    //B init:
    search_init_sincos_omp_kkk(cmd, gd, &hist1sincos, ifile);
    search_init_sincos_omp_kkk_N(cmd, gd, &hist1sincos_N, ifile);
    for (n = 1; n <= cmd->sizeHistN; n++) {
        hist1sincos.histNNSubthread[n] = 0.0;
    }
    CLRM_ext_ext(hist1sincos.histXithreadcos, cmd->mChebyshev+1,
                 cmd->sizeHistN);
    CLRM_ext_ext(hist1sincos.histXithreadsin, cmd->mChebyshev+1,
                 cmd->sizeHistN);
    CLRM_ext_ext(hist1sincos_N.histXithreadcos, cmd->mChebyshev+1,
                 cmd->sizeHistN);
    CLRM_ext_ext(hist1sincos_N.histXithreadsin, cmd->mChebyshev+1,
                 cmd->sizeHistN);
    for (n = 1; n <= cmd->sizeHistN; n++) {
        hist1sincos_N.histNNSubthread[n] = 0.0;
    }
    //E

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

    if (!cballs_opt_no_one_ball(cmd))
        sumcell_balls6_omp_N(cmd, gd, histb->interact, cptr, (bodyptr) p,
                           &hist1sincos, &hist1sincos_N);
    sumnode_balls6_omp_N(cmd, gd, bptr, histb->interact + histb->actlen,
                       (bodyptr) p, &hist1sincos, &hist1sincos_N);

    computeBodyProperties_sincos_kkk_sum_balls6_omp(cmd, gd,
                                    p, nbody, histsincos, &hist1sincos);
    computeBodyProperties_sincos_kkk_sum_balls6_omp_N(cmd, gd,
                                    p, nbody, histsincos_N, &hist1sincos_N);

    histsincos->nbbcalcthread += histb->interact + histb->actlen - bptr;
    histsincos->nbccalcthread += cptr - histb->interact;

    ip = p - bodytable[gd->iCatalogs[0]] + 1;
    if (ip%cmd->stepState == 0) {
        verb_log_print(cmd->verbose_log, gd->outlog, " - Completed pivot: %ld\n", ip);
    }

    search_free_sincos_omp_kkk_N(cmd, gd, &hist1sincos_N);
    search_free_sincos_omp_kkk(cmd, gd, &hist1sincos);
    return SUCCESS;
}

local void sumnode_balls6_omp_N(struct cmdline_data* cmd,
                              struct  global_data* gd,
                              cellptr start, cellptr finish, bodyptr p0,
                              gdhistptr_sincos_omp_kkk hist,
                              gdhistptr_sincos_omp_kkk_N histN)
{
    cellptr p;
    compute_vector dr;
    real dr1;
    nodeptr pb;
    INTEGER ibodycount=0;
    int n;
    real xi;
    real xiN;
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
                    histN->histNNSubthread[n] = histN->histNNSubthread[n] + 1.;
#ifdef TWOPCF
                    balls4_accumulate_neighbor_2pcf(hist, pb, n);
#endif
                    xi = Kappa(pb);
                    xiN = 1.0;
                    //B Component of pb with respect to the axis of reference
                    if (!cballs_angular_phase(Pos(p0), dr, &cosphi, &sinphi))
                        continue;
                    CHEBYSHEVTUOMPKKK;
                    CHEBYSHEVTUOMPKKKN;
                } // ! 1 < n < HistN
            } // ! dr1 > rmin
        } // ! accept_body
    } // ! loop p
}

local void sumcell_balls6_omp_N(struct cmdline_data* cmd,
                                struct  global_data* gd,
                                cellptr start, cellptr finish, bodyptr p0,
                                gdhistptr_sincos_omp_kkk hist,
                                gdhistptr_sincos_omp_kkk_N histN)
{
    cellptr p;
    compute_vector dr;
    real dr1;
    nodeptr pb;
    INTEGER ibodycount=0;
    int n;
    real xi;
    real xiN;
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
                    histN->histNNSubthread[n] = histN->histNNSubthread[n] + 1.;
#ifdef TWOPCF
                    balls4_accumulate_neighbor_2pcf(hist, pb, n);
#endif
                    //B correction 2025-04-06
                    // needs to multiply xi by Weight(p) of the cell q
#ifdef NONORMHIST
                    if (cballs_opt_no_normalize_histzeta(cmd)) {
                        xi = Nb(pb)*Kappa(pb);
                        xiN = Nb(pb);
                    } else {
                        xi = Kappa(pb);
                        xiN = 1.0;
                    }
#else
                    xi = Kappa(pb);
                    xiN = 1.0;
#endif
//                    xi = Kappa(pb);
//                    xiN = 1.0;
                    //B Component of pb with respect to the axis of reference
                    if (!cballs_angular_phase(Pos(p0), dr, &cosphi, &sinphi))
                        continue;
                    CHEBYSHEVTUOMPKKK;
                    CHEBYSHEVTUOMPKKKN;
                } // ! 1 < n < HistN
            } // ! dr1 > rmin
        } // ! accept_body
    } // ! loop p
}

local void sumcellcell_balls6_omp_N(struct cmdline_data* cmd,
                                    struct  global_data* gd,
                                    cellptr start, cellptr finish, nodeptr p0,
                                    gdhistptr_sincos_omp_kkk hist,
                                    gdhistptr_sincos_omp_kkk_N histN)
{
    cellptr p;
    compute_vector dr;
    real dr1;
    nodeptr pb;
    INTEGER ibodycount=0;
    int n;
    real xi;
    real xiN;
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
                    histN->histNNSubthread[n] =
                                    histN->histNNSubthread[n] + Nb(pb);
#ifdef TWOPCF
                    balls4_accumulate_cell_pair_2pcf(hist, p0, pb, n);
#endif
                    //B correction 2025-04-06
                    // needs to multiply xi by Weight(p) of the cell q
#ifdef NONORMHIST
                    if (cballs_opt_no_normalize_histzeta(cmd)) {
                        xi = Nb(pb)*Kappa(pb);
                        xiN = Nb(pb);
                    } else {
                        xi = Kappa(pb);
                        xiN = 1.0;
                    }
#else
                    xi = Kappa(pb);
                    xiN = 1.0;
#endif
//                    xi = Kappa(pb);
//                    xiN = 1.0;
                    //B Component of pb with respect to the axis of reference
                    if (!cballs_angular_phase(Pos(p0), dr, &cosphi, &sinphi))
                        continue;
                    CHEBYSHEVTUOMPKKK;
                    CHEBYSHEVTUOMPKKKN;
                } // ! 1 < n < HistN
            } // ! dr1 > rmin
        } // ! accept_body
    } // ! loop p

    hist->ncccalcthread += 1;
}

local int computeBodyProperties_sincos_kkk_sum_balls6_omp_N(
                                struct  cmdline_data* cmd,
                                struct  global_data* gd,
                                bodyptr p, INTEGER nbody,
                                gdhistptr_sincos_omp_kkk_N histsincos,
                                gdhistptr_sincos_omp_kkk_N hist1sincos)
{
    int n;
    int m;
    real xi;

#ifdef NONORMHIST
    xi = 1.0;
#else
    xi = 1.0/nbody;
#endif

    for (m=1; m<=cmd->mChebyshev+1; m++) {
#ifndef NONORMHIST
        for (n=1; n<=cmd->sizeHistN; n++) {
            hist1sincos->histXithreadcos[m][n] /=
            MAX(hist1sincos->histNNSubthread[n],1.0);
            hist1sincos->histXithreadsin[m][n] /=
            MAX(hist1sincos->histNNSubthread[n],1.0);
        }
#endif
    }

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

//    for (n = 1; n <= cmd->sizeHistN; n++) {
//        histsincos->histNNSubthread[n] += hist1sincos->histNNSubthread[n];
//    }

    return SUCCESS;
}

#endif // ! NMultipoles


#ifdef TWOPCF
local int search_compute_Xi_kkk_balls4(struct cmdline_data* cmd,
                                      struct global_data* gd,
                                      INTEGER nbody,
                                      gdlptr_sincos_omp_kkk gdl)
{
    int k;
    int n;
    real volume = 1.0;
    real body_count = (real)nbody;

    DO_COORD(k)
        volume *= gd->Box[k];

    for (n = 1; n <= cmd->sizeHistN; n++) {
        real r0;
        real r1;

        if (gdl->histNN[n] == 0.0 || body_count == 0.0 || volume == 0.0) {
            gdl->histCF[n] = 0.0;
            continue;
        }

        if (cmd->useLogHist) {
            if (cmd->rminHist == 0.0) {
                r0 = rpow(10.0, ((real)(n-cmd->sizeHistN))
                    /cmd->logHistBinsPD + rlog10(cmd->rangeN));
                r1 = rpow(10.0, ((real)(n+1-cmd->sizeHistN))
                    /cmd->logHistBinsPD + rlog10(cmd->rangeN));
            } else {
                r0 = rpow(10.0, rlog10(cmd->rminHist) + (real)n*gd->deltaR);
                r1 = rpow(10.0, rlog10(cmd->rminHist)
                    + (real)(n+1)*gd->deltaR);
            }
        } else {
            r0 = cmd->rminHist + (real)(n-1)*gd->deltaR;
            r1 = cmd->rminHist + (real)n*gd->deltaR;
        }

#if NDIM == 3
        if (cballs_opt_cute_box(cmd)) {
            real shell_volume = 4.0*PI*(r1*r1*r1-r0*r0*r0)/3.0;
            real expected = body_count*body_count*shell_volume/volume;
            gdl->histCF[n] =
                cballs_normalize_or_zero(gdl->histNN[n], expected) - 1.0;
        } else {
            real norm = volume/(2.0*PI*rpow(gd->deltaR,3.0)
                                      *body_count*body_count);
            gdl->histCF[n] =
                gdl->histNN[n]*norm/rsqr((int)n-0.5) - 1.0;
        }
#else
        {
            real norm = volume/(PI*rpow(gd->deltaR,2.0)
                                    *body_count*body_count);
            gdl->histCF[n] = gdl->histNN[n]*norm/((real)n-0.5) - 1.0;
        }
#endif
    }

    return SUCCESS;
}

local int search_compute_HistN_kkk_balls4(struct cmdline_data* cmd,
                                         struct global_data* gd,
                                         INTEGER nbody,
                                         gdlptr_sincos_omp_kkk gdl)
{
    int n;

    if (!cballs_opt_asymmetric(cmd)) {
        for (n = 1; n <= cmd->sizeHistN; n++)
            gdl->histNN[n] *= 0.5;
    }

    if (cballs_opt_and_cf(cmd))
        return search_compute_Xi_kkk_balls4(cmd, gd, nbody, gdl);

    return SUCCESS;
}

local int PrintHistNN_kkk_balls4(struct cmdline_data* cmd,
                                struct global_data* gd,
                                gdlptr_sincos_omp_kkk gdl)
{
    string routineName = "PrintHistNN_kkk_balls4";
    real rBin, rbinlog;
    int n;
    stream outstr;

    OPEN_OUTPUT_OR_FAIL(outstr, gd->fpfnamehistNNFileName, "w!");
    for (n = 1; n <= cmd->sizeHistN; n++) {
        if (cmd->useLogHist) {
            if (cmd->rminHist == 0.0)
                rbinlog = ((real)(n-cmd->sizeHistN))/cmd->logHistBinsPD
                    + rlog10(cmd->rangeN);
            else
                rbinlog = rlog10(cmd->rminHist) + ((real)n-0.5)*gd->deltaR;
            rBin = rpow(10.0,rbinlog);
        } else {
            rBin = cmd->rminHist + ((real)n-0.5)*gd->deltaR;
        }
        WRITE_OUTPUT_OR_FAIL(outstr, gd->fpfnamehistNNFileName,
                             "%16.8e %16.8e\n", rBin, gdl->histNN[n]);
    }
    CLOSE_OUTPUT_OR_FAIL(outstr, gd->fpfnamehistNNFileName);

    if (cballs_opt_and_cf(cmd))
        return PrintHistCF_kkk_balls4(cmd, gd, gdl);
    return SUCCESS;
}

local int PrintHistCF_kkk_balls4(struct cmdline_data* cmd,
                                struct global_data* gd,
                                gdlptr_sincos_omp_kkk gdl)
{
    string routineName = "PrintHistCF_kkk_balls4";
    real rBin, rbinlog;
    int n;
    stream outstr;

    OPEN_OUTPUT_OR_FAIL(outstr, gd->fpfnamehistCFFileName, "w!");
    for (n = 1; n <= cmd->sizeHistN; n++) {
        if (cmd->useLogHist) {
            if (cmd->rminHist == 0.0)
                rbinlog = ((real)(n-cmd->sizeHistN))/cmd->logHistBinsPD
                    + rlog10(cmd->rangeN);
            else
                rbinlog = rlog10(cmd->rminHist) + ((real)n-0.5)*gd->deltaR;
            rBin = rpow(10.0,rbinlog);
        } else {
            rBin = cmd->rminHist + ((real)n-0.5)*gd->deltaR;
        }
        WRITE_OUTPUT_OR_FAIL(outstr, gd->fpfnamehistCFFileName,
                             "%16.8e %16.8e\n", rBin, gdl->histCF[n]);
    }
    CLOSE_OUTPUT_OR_FAIL(outstr, gd->fpfnamehistCFFileName);
    return SUCCESS;
}

local int PrintHistXi2pcf_kkk_balls4(struct cmdline_data* cmd,
                                    struct global_data* gd,
                                    gdlptr_sincos_omp_kkk gdl)
{
    string routineName = "PrintHistXi2pcf_kkk_balls4";
    real rBin, rbinlog;
    int n;
    stream outstr;
    char namebuf[256];

    if (format_checked(namebuf, sizeof(namebuf), "namebuf", "%s%s%s",
                       gd->fpfnamehistXi2pcfFileName,
                       cmd->suffixOutFiles, EXTFILES) != 0)
        return FAILURE;
    OPEN_OUTPUT_OR_FAIL(outstr, namebuf, "w!");
    for (n = 1; n <= cmd->sizeHistN; n++) {
        if (cmd->useLogHist) {
            if (cmd->rminHist == 0.0)
                rbinlog = ((real)(n-cmd->sizeHistN))/cmd->logHistBinsPD
                    + rlog10(cmd->rangeN);
            else
                rbinlog = rlog10(cmd->rminHist) + ((real)n-0.5)*gd->deltaR;
            rBin = rpow(10.0,rbinlog);
        } else {
            rBin = cmd->rminHist + ((real)n-0.5)*gd->deltaR;
        }
        if (cballs_opt_rbin_arcmin(cmd))
            rBin *= RADTOARCMIN;
        else if (cballs_opt_rbin_degree(cmd))
            rBin *= RADTOARCMIN/60.0;
        WRITE_OUTPUT_OR_FAIL(outstr, namebuf,
                             "%16.8e %16.8e\n", rBin, gdl->histXi2pcf[n]);
    }
    CLOSE_OUTPUT_OR_FAIL(outstr, namebuf);
    return SUCCESS;
}
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
    free(data);
#else
    free_dvector(data,1,NP);
#endif

    free_dmatrix3D(histZetaG_Im,1,NP,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix3D(histZetaG,1,NP,1,cmd->sizeHistN,1,cmd->sizeHistN);

    return SUCCESS;
}



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

/*
// Edge-correct each matrix ZetaM for each m multipole
local int HistZetaM_sincos_edge_effects(struct  cmdline_data* cmd,
                                        struct  global_data* gd,
                                        gdhistptr_sincos_omp_kkk hist,
                                        gdhistptr_sincos_omp_kkk_N histN)
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
                hist->histZetaMtmp[m][n1][n2] = mat5[m][n1][n2];
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
*/

// Saves matrix ZetaM for each m multipole


#endif // ! NONORMHIST

#endif // ! NMultipoles

#undef MHISTZETAHEADER
#undef MHISTZETA


//E Saving histograms section: case KKKCORRELATION:

