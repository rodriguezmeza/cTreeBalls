// Use:
//#include "datastruc_defs_balls_omp.h"

#ifndef _datastruc_defs_balls_omp_h
#define _datastruc_defs_balls_omp_h

//B CHEBYSHEVTUOMP definitions
// CLEAN THIS DEFINITIONS... ONLY ONE IS VALID, FIRST ONE, JUST BELOW
#ifdef TREENODEALLBODIES
#ifdef SINCOS
#define CHEBYSHEVTUOMP                                      \
{                                                           \
    histccsincos->ChebsT[1] = 1.0;                                        \
    xicosmphi = xi * histccsincos->ChebsT[1];                             \
    histccsincos->histXithreadcos[1][n] += xicosmphi;                        \
    histccsincos->ChebsT[2] = cosphi;                                     \
    xicosmphi = xi * histccsincos->ChebsT[2];                             \
    histccsincos->histXithreadcos[2][n] += xicosmphi;                        \
    histccsincos->ChebsT[3] = 2.0*(cosphi)*(cosphi) - (1.0);              \
    xicosmphi = xi * histccsincos->ChebsT[3];                             \
    histccsincos->histXithreadcos[3][n] += xicosmphi;                        \
    histccsincos->ChebsU[1] = 0.0;                                        \
    xisinmphi = xi * histccsincos->ChebsU[1] * sinphi;                    \
    histccsincos->histXithreadsin[1][n] += xisinmphi;                        \
    histccsincos->ChebsU[2] = 1.0;                                        \
    xisinmphi = xi * histccsincos->ChebsU[2] * sinphi;                    \
    histccsincos->histXithreadsin[2][n] += xisinmphi;                        \
    histccsincos->ChebsU[3] = 2.0*cosphi;                                 \
    xisinmphi = xi * histccsincos->ChebsU[3] * sinphi;                    \
    histccsincos->histXithreadsin[3][n] += xisinmphi;                        \
    for (m=4; m<=cmd->mChebyshev+1; m++){                    \
        histccsincos->ChebsT[m] = 2.0*(cosphi)*histccsincos->ChebsT[m-1] - histccsincos->ChebsT[m-2]; \
        xicosmphi = xi * histccsincos->ChebsT[m];                         \
        histccsincos->histXithreadcos[m][n] += xicosmphi;                    \
        histccsincos->ChebsU[m] = 2.0*(cosphi)*histccsincos->ChebsU[m-1] - histccsincos->ChebsU[m-2]; \
        xisinmphi = xi * histccsincos->ChebsU[m] * sinphi;                \
        histccsincos->histXithreadsin[m][n] += xisinmphi;                    \
    }}
#else // ! SINCOS
#define CHEBYSHEVTUOMP                                      \
{                                                           \
    hist->ChebsT[1] = 1.0;                                        \
    xicosmphi = xi * hist->ChebsT[1];                             \
    hist->histXithreadcos[1][n] += xicosmphi;                        \
    hist->ChebsT[2] = cosphi;                                     \
    xicosmphi = xi * hist->ChebsT[2];                             \
    hist->histXithreadcos[2][n] += xicosmphi;                        \
    hist->ChebsT[3] = 2.0*(cosphi)*(cosphi) - (1.0);              \
    xicosmphi = xi * hist->ChebsT[3];                             \
    hist->histXithreadcos[3][n] += xicosmphi;                        \
    hist->ChebsU[1] = 0.0;                                        \
    xisinmphi = xi * hist->ChebsU[1] * sinphi;                    \
    hist->histXithreadsin[1][n] += xisinmphi;                        \
    hist->ChebsU[2] = 1.0;                                        \
    xisinmphi = xi * hist->ChebsU[2] * sinphi;                    \
    hist->histXithreadsin[2][n] += xisinmphi;                        \
    hist->ChebsU[3] = 2.0*cosphi;                                 \
    xisinmphi = xi * hist->ChebsU[3] * sinphi;                    \
    hist->histXithreadsin[3][n] += xisinmphi;                        \
    for (m=4; m<=cmd->mChebyshev+1; m++){                    \
        hist->ChebsT[m] = 2.0*(cosphi)*hist->ChebsT[m-1] - hist->ChebsT[m-2]; \
        xicosmphi = xi * hist->ChebsT[m];                         \
        hist->histXithreadcos[m][n] += xicosmphi;                    \
        hist->ChebsU[m] = 2.0*(cosphi)*hist->ChebsU[m-1] - hist->ChebsU[m-2]; \
        xisinmphi = xi * hist->ChebsU[m] * sinphi;                \
        hist->histXithreadsin[m][n] += xisinmphi;                    \
    }}
#endif // ! SINCOS
#else // ! TREENODEALLBODIES

#ifdef TREENODEBALLS4

#ifdef SINCOS
#define CHEBYSHEVTUOMP                                      \
{                                                           \
    histccsincos->ChebsT[1] = 1.0;                                        \
    xicosmphi = xi * histccsincos->ChebsT[1];                             \
    histccsincos->histXithreadcos[1][n] += xicosmphi;                        \
    histccsincos->ChebsT[2] = cosphi;                                     \
    xicosmphi = xi * histccsincos->ChebsT[2];                             \
    histccsincos->histXithreadcos[2][n] += xicosmphi;                        \
    histccsincos->ChebsT[3] = 2.0*(cosphi)*(cosphi) - (1.0);              \
    xicosmphi = xi * histccsincos->ChebsT[3];                             \
    histccsincos->histXithreadcos[3][n] += xicosmphi;                        \
    histccsincos->ChebsU[1] = 0.0;                                        \
    xisinmphi = xi * histccsincos->ChebsU[1] * sinphi;                    \
    histccsincos->histXithreadsin[1][n] += xisinmphi;                        \
    histccsincos->ChebsU[2] = 1.0;                                        \
    xisinmphi = xi * histccsincos->ChebsU[2] * sinphi;                    \
    histccsincos->histXithreadsin[2][n] += xisinmphi;                        \
    histccsincos->ChebsU[3] = 2.0*cosphi;                                 \
    xisinmphi = xi * histccsincos->ChebsU[3] * sinphi;                    \
    histccsincos->histXithreadsin[3][n] += xisinmphi;                        \
    for (m=4; m<=cmd->mChebyshev+1; m++){                    \
        histccsincos->ChebsT[m] = 2.0*(cosphi)*histccsincos->ChebsT[m-1] - histccsincos->ChebsT[m-2]; \
        xicosmphi = xi * histccsincos->ChebsT[m];                         \
        histccsincos->histXithreadcos[m][n] += xicosmphi;                    \
        histccsincos->ChebsU[m] = 2.0*(cosphi)*histccsincos->ChebsU[m-1] - histccsincos->ChebsU[m-2]; \
        xisinmphi = xi * histccsincos->ChebsU[m] * sinphi;                \
        histccsincos->histXithreadsin[m][n] += xisinmphi;                    \
    }}
#else // ! SINCOS
#define CHEBYSHEVTUOMP                                      \
{                                                           \
    hist->ChebsT[1] = 1.0;                                        \
    xicosmphi = xi * hist->ChebsT[1];                             \
    hist->histXithreadcos[1][n] += xicosmphi;                        \
    hist->ChebsT[2] = cosphi;                                     \
    xicosmphi = xi * hist->ChebsT[2];                             \
    hist->histXithreadcos[2][n] += xicosmphi;                        \
    hist->ChebsT[3] = 2.0*(cosphi)*(cosphi) - (1.0);              \
    xicosmphi = xi * hist->ChebsT[3];                             \
    hist->histXithreadcos[3][n] += xicosmphi;                        \
    hist->ChebsU[1] = 0.0;                                        \
    xisinmphi = xi * hist->ChebsU[1] * sinphi;                    \
    hist->histXithreadsin[1][n] += xisinmphi;                        \
    hist->ChebsU[2] = 1.0;                                        \
    xisinmphi = xi * hist->ChebsU[2] * sinphi;                    \
    hist->histXithreadsin[2][n] += xisinmphi;                        \
    hist->ChebsU[3] = 2.0*cosphi;                                 \
    xisinmphi = xi * hist->ChebsU[3] * sinphi;                    \
    hist->histXithreadsin[3][n] += xisinmphi;                        \
    for (m=4; m<=cmd->mChebyshev+1; m++){                    \
        hist->ChebsT[m] = 2.0*(cosphi)*hist->ChebsT[m-1] - hist->ChebsT[m-2]; \
        xicosmphi = xi * hist->ChebsT[m];                         \
        hist->histXithreadcos[m][n] += xicosmphi;                    \
        hist->ChebsU[m] = 2.0*(cosphi)*hist->ChebsU[m-1] - hist->ChebsU[m-2]; \
        xisinmphi = xi * hist->ChebsU[m] * sinphi;                \
        hist->histXithreadsin[m][n] += xisinmphi;                    \
    }}
#endif // ! SINCOS

#else // ! TREENODEBALLS4

#define CHEBYSHEVTUOMP                                      \
{ real xicosmphi; int m;                                              \
    histccsincos->ChebsT[1] = 1.0;                                        \
    xicosmphi = xi * histccsincos->ChebsT[1];                             \
    histccsincos->histXithreadcos[1][n] += xj*xicosmphi;                        \
    histccsincos->ChebsT[2] = cosphi;                                     \
    xicosmphi = xi * histccsincos->ChebsT[2];                             \
    histccsincos->histXithreadcos[2][n] += xj*xicosmphi;                        \
    histccsincos->ChebsT[3] = 2.0*(cosphi)*(cosphi) - (1.0);              \
    xicosmphi = xi * histccsincos->ChebsT[3];                             \
    histccsincos->histXithreadcos[3][n] += xj*xicosmphi;                        \
    histccsincos->ChebsU[1] = 0.0;                                        \
    xisinmphi = xi * histccsincos->ChebsU[1] * sinphi;                    \
    histccsincos->histXithreadsin[1][n] += xj*xisinmphi;                        \
    histccsincos->ChebsU[2] = 1.0;                                        \
    xisinmphi = xi * histccsincos->ChebsU[2] * sinphi;                    \
    histccsincos->histXithreadsin[2][n] += xj*xisinmphi;                        \
    histccsincos->ChebsU[3] = 2.0*cosphi;                                 \
    xisinmphi = xi * histccsincos->ChebsU[3] * sinphi;                    \
    histccsincos->histXithreadsin[3][n] += xj*xisinmphi;                        \
    for (m=4; m<=cmd->mChebyshev+1; m++){                    \
        histccsincos->ChebsT[m] = 2.0*(cosphi)*histccsincos->ChebsT[m-1] - histccsincos->ChebsT[m-2]; \
        xicosmphi = xi * histccsincos->ChebsT[m];                         \
        histccsincos->histXithreadcos[m][n] += xj*xicosmphi;                    \
        histccsincos->ChebsU[m] = 2.0*(cosphi)*histccsincos->ChebsU[m-1] - histccsincos->ChebsU[m-2]; \
        xisinmphi = xi * histccsincos->ChebsU[m] * sinphi;                \
        histccsincos->histXithreadsin[m][n] += xj*xisinmphi;                    \
    }}
#endif // ! TREENODEBALLS4
#endif // ! TREENODEALLBODIES
//E CHEBYSHEVTUOMP definitions


//B Faster if commented declarations
// (comment and uncomment also lines in:
// sumnodes_bb_omp, sumnodes_bc_omp, sumnodes_cb_omp and sumnodes_cc_omp
#ifdef TREENODEALLBODIES
#define CHEBYSHEVOMPBALLSCC                                 \
{real xicosmphi; int m; \
      histcc->Chebs[1] = 1.0;                                    \
   xicosmphi = xi * histcc->Chebs[1];                         \
   histcc->histXithread[1][n] += xicosmphi;                \
   histcc->Chebs[2] = cosphi;                                 \
   xicosmphi = xi * histcc->Chebs[2];                         \
   histcc->histXithread[2][n] += xicosmphi;                \
   histcc->Chebs[3] = 2.0*(cosphi)*(cosphi) - (1.0);          \
   xicosmphi = xi * histcc->Chebs[3];                         \
   histcc->histXithread[3][n] += xicosmphi;                \
   for (m=4; m<=cmd->mChebyshev+1; m++){                     \
       histcc->Chebs[m] = 2.0*(cosphi)*histcc->Chebs[m-1] - histcc->Chebs[m-2];  \
       xicosmphi = xi * histcc->Chebs[m];                     \
       histcc->histXithread[m][n] += xicosmphi;            \
   }}
#else // ! TREENODEALLBODIES
#define CHEBYSHEVOMPBALLSCC                                   \
{real xicosmphi; int m;                                       \
      histcc->Chebs[1] = 1.0;                                 \
   xicosmphi = xi * histcc->Chebs[1];                         \
   histcc->histXithread[1][n] += xicosmphi;                \
   histcc->Chebs[2] = cosphi;                                 \
   xicosmphi = xi * histcc->Chebs[2];                         \
   histcc->histXithread[2][n] += xicosmphi;                \
   histcc->Chebs[3] = 2.0*(cosphi)*(cosphi) - (1.0);          \
   xicosmphi = xi * histcc->Chebs[3];                         \
   histcc->histXithread[3][n] += xicosmphi;                \
   for (m=4; m<=cmd->mChebyshev+1; m++){                     \
       histcc->Chebs[m] = 2.0*(cosphi)*histcc->Chebs[m-1] - histcc->Chebs[m-2];  \
       xicosmphi = xi * histcc->Chebs[m];                     \
       histcc->histXithread[m][n] += xicosmphi;            \
   }}
#endif // ! TREENODEALLBODIES


#endif	// ! _datastruc_defs_balls_omp_h
