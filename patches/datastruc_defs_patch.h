// Use:
//NOLSST:
//#include "datastruc_defs_patch.h"

#ifndef _datastruc_defs_patch_h
#define _datastruc_defs_patch_h

#define NOCHEBYSHEVTU                                       \
{real ChebT, ChebU;for (m=1; m<=cmd.mchebyshev+1; m++){     \
    ChebT = rcos((real)(m-1) * racos(cosphi));              \
    ChebU = sinphi * rsin((real)(m-1) * racos(cosphi))/rsin(racos(cosphi)); \
    xicosmphi = xi * ChebT;                                 \
    xisinmphi = xi * ChebU;                                 \
    gd.histXicos[m][n] += xicosmphi;                        \
    gd.histXisin[m][n] += xisinmphi;                        \
    }}

#define CHEBYSHEVTU                                         \
{real ChebU;                                                \
    ChebsT[1] = 1.0;                                        \
    xicosmphi = xi * ChebsT[1];                             \
    gd.histXicos[1][n] += xicosmphi;                        \
    ChebsT[2] = cosphi;                                     \
    xicosmphi = xi * ChebsT[2];                             \
    gd.histXicos[2][n] += xicosmphi;                        \
    ChebsT[3] = 2.0*(cosphi)*(cosphi) - (1.0);              \
    xicosmphi = xi * ChebsT[3];                             \
    gd.histXicos[3][n] += xicosmphi;                        \
    ChebsU[1] = 0.0;                                        \
    xisinmphi = xi * ChebsU[1] * sinphi;                    \
    gd.histXisin[1][n] += xisinmphi;                        \
    ChebsU[2] = 1.0;                                        \
    xisinmphi = xi * ChebsU[2] * sinphi;                    \
    gd.histXisin[2][n] += xisinmphi;                        \
    ChebsU[3] = 2.0*cosphi;                                 \
    xisinmphi = xi * ChebsU[3] * sinphi;                    \
    gd.histXisin[3][n] += xisinmphi;                        \
    for (m=4; m<=cmd.mchebyshev+1; m++){                    \
        ChebsT[m] = 2.0*(cosphi)*ChebsT[m-1] - ChebsT[m-2]; \
        xicosmphi = xi * ChebsT[m];                         \
        gd.histXicos[m][n] += xicosmphi;                    \
        ChebsU[m] = 2.0*(cosphi)*ChebsU[m-1] - ChebsU[m-2]; \
        xisinmphi = xi * ChebsU[m] * sinphi;                \
        gd.histXisin[m][n] += xisinmphi;                    \
    }}


#define DIRECTSIMPLE            12
#define DIRECTSIMPLESINCOS      19

#endif	// ! _datastruc_defs_01_h
