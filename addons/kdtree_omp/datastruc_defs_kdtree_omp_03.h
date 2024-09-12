// Use:
//#include "datastruc_defs_kdtree_omp_03.h"

#ifndef _datastruc_defs_kdtree_omp_03_h
#define _datastruc_defs_kdtree_omp_03_h

#ifdef NMultipoles
//B kdtree-omp
#ifdef MANUALCHEBYSHEV
#define NCHEBYSHEVTUOMPSINCOS                                      \
{REAL xicosmphi,xisinmphi; int m;                                 \
    REAL cosphi2, cosphi3, cosphi4;                               \
    cosphi2 = cosphi*cosphi; cosphi3 = cosphi2*cosphi;            \
    cosphi4 = cosphi3*cosphi;                                     \
    histN->ChebsT[1] = 1.0;                                        \
    xicosmphi = xiN * histN->ChebsT[1];                             \
    histN->histXithreadcos[1][n] += xicosmphi;                     \
    histN->ChebsT[2] = cosphi;                                     \
    xicosmphi = xiN * histN->ChebsT[2];                             \
    histN->histXithreadcos[2][n] += xicosmphi;                     \
    histN->ChebsT[3] = 2.0*cosphi2 - (1.0);                        \
    xicosmphi = xiN * histN->ChebsT[3];                             \
    histN->histXithreadcos[3][n] += xicosmphi;                     \
    histN->ChebsT[4] = -3.0*cosphi + 4.0*cosphi3;                  \
    xicosmphi = xiN * histN->ChebsT[4];                             \
    histN->histXithreadcos[4][n] += xicosmphi;                     \
    histN->ChebsT[5] = 1.0 - 8.0*cosphi2 + 8.0*cosphi4;            \
    xicosmphi = xiN * histN->ChebsT[5];                             \
    histN->histXithreadcos[5][n] += xicosmphi;                     \
    histN->ChebsT[6] = 5.0*cosphi - 20.0*cosphi3 + 16.0*cosphi4*cosphi; \
    xicosmphi = xiN * histN->ChebsT[6];                             \
    histN->histXithreadcos[6][n] += xicosmphi;                     \
    histN->ChebsU[1] = 0.0;                                        \
    xisinmphi = xiN * histN->ChebsU[1] * sinphi;                    \
    histN->histXithreadsin[1][n] += xisinmphi;                     \
    histN->ChebsU[2] = 1.0;                                        \
    xisinmphi = xiN * histN->ChebsU[2] * sinphi;                    \
    histN->histXithreadsin[2][n] += xisinmphi;                     \
    histN->ChebsU[3] = 2.0*cosphi;                                 \
    xisinmphi = xiN * histN->ChebsU[3] * sinphi;                    \
    histN->histXithreadsin[3][n] += xisinmphi;                     \
    histN->ChebsU[4] = -1.0 + 4.0*cosphi2;                         \
    xisinmphi = xiN * histN->ChebsU[4] * sinphi;                    \
    histN->histXithreadsin[4][n] += xisinmphi;                     \
    histN->ChebsU[5] = -4.0*cosphi + 8.0*cosphi3;                  \
    xisinmphi = xiN * histN->ChebsU[5] * sinphi;                    \
    histN->histXithreadsin[5][n] += xisinmphi;                     \
    histN->ChebsU[6] = 1.0 -12.0*cosphi2 + 16.0*cosphi4;           \
    xisinmphi = xiN * histN->ChebsU[6] * sinphi;                    \
    histN->histXithreadsin[6][n] += xisinmphi;                     \
    for (m=7; m<=cmd->mChebyshev+1; m++){                          \
        histN->ChebsT[m] = 2.0*(cosphi)*histN->ChebsT[m-1] - histN->ChebsT[m-2]; \
        xicosmphi = xiN * histN->ChebsT[m];                         \
        histN->histXithreadcos[m][n] += xicosmphi;                 \
        histN->ChebsU[m] = 2.0*(cosphi)*histN->ChebsU[m-1] - histN->ChebsU[m-2]; \
        xisinmphi = xiN * histN->ChebsU[m] * sinphi;                \
        histN->histXithreadsin[m][n] += xisinmphi;                 \
    }}
#else // ! MANUALCHEBYSHEV
#define NCHEBYSHEVTUOMPSINCOS                                      \
{real xicosmphi,xisinmphi; int m;                                 \
    histN->ChebsT[1] = 1.0;                                        \
    xicosmphi = xiN * histN->ChebsT[1];                             \
    histN->histXithreadcos[1][n] += xicosmphi;                     \
    histN->ChebsT[2] = cosphi;                                     \
    xicosmphi = xiN * histN->ChebsT[2];                             \
    histN->histXithreadcos[2][n] += xicosmphi;                     \
    histN->ChebsT[3] = 2.0*(cosphi)*(cosphi) - (1.0);              \
    xicosmphi = xiN * histN->ChebsT[3];                             \
    histN->histXithreadcos[3][n] += xicosmphi;                     \
    histN->ChebsU[1] = 0.0;                                        \
    xisinmphi = xiN * histN->ChebsU[1] * sinphi;                    \
    histN->histXithreadsin[1][n] += xisinmphi;                     \
    histN->ChebsU[2] = 1.0;                                        \
    xisinmphi = xiN * histN->ChebsU[2] * sinphi;                    \
    histN->histXithreadsin[2][n] += xisinmphi;                     \
    histN->ChebsU[3] = 2.0*cosphi;                                 \
    xisinmphi = xiN * histN->ChebsU[3] * sinphi;                    \
    histN->histXithreadsin[3][n] += xisinmphi;                     \
    for (m=4; m<=cmd->mChebyshev+1; m++){                          \
        histN->ChebsT[m] = 2.0*(cosphi)*histN->ChebsT[m-1] - histN->ChebsT[m-2]; \
        xicosmphi = xiN * histN->ChebsT[m];                         \
        histN->histXithreadcos[m][n] += xicosmphi;                 \
        histN->ChebsU[m] = 2.0*(cosphi)*histN->ChebsU[m-1] - histN->ChebsU[m-2]; \
        xisinmphi = xiN * histN->ChebsU[m] * sinphi;                \
        histN->histXithreadsin[m][n] += xisinmphi;                 \
    }}
#endif // ! MANUALCHEBYSHEV
//E kdtree-omp
#endif // ! NMultipoles

#endif	// ! _datastruc_defs_kdtree_omp_03_h
