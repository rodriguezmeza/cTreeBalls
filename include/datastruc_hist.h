/*==============================================================================
 HEADER: cmdline_data.h		[cTreeBalls]
 Written by: Mario A. Rodriguez-Meza
 Starting date: april 2023
 Purpose: Definitions of global variables and parameters
 Language: C
 Use: '#include "datastruc_hist.h"
 Major revisions:
 ==============================================================================*/
//        1          2          3          4        ^ 5          6          7

#ifndef _datastruc_hist_h
#define _datastruc_hist_h

//B Structure definitions for histograms
//

#ifdef USEGSL
typedef struct {
    int m;
    gsl_matrix_complex *histZetaM;
} mMatrix, *mMatrix_ptr;

global mMatrix_ptr histZetaMatrix;
#endif

typedef struct {
    real **xiOUTVPcos;
    real **xiOUTVPsin;
    real **xiOUTVPsincos;
    // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
    real **xiOUTVPcossin;
    real **histZetaMtmpcos;
    real **histZetaMtmpsin;
    real **histZetaMtmpsincos;
    // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
    real **histZetaMtmpcossin;
    real *ChebsT;
    real *ChebsU;
    real ***histZetaMthreadcos;
    real ***histZetaMthreadsin;
    real ***histZetaMthreadsincos;
    // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
    real ***histZetaMthreadcossin;
    realptr histNthread;
    realptr histNNSubthread;
// 2pcf
    realptr histNNSubXi2pcfthread;
//B kappa Avg Rmin
    realptr histNNSubXi2pcfthreadp;
    realptr histNNSubXi2pcfthreadtotal;
//E
//
    real **histXithreadcos;
    real **histXithreadsin;
    real *histXi2pcfthread;
    real *histXi2pcfthreadsub;

#ifdef SINGLEP
    float q0[NDIM];
    float drpq2, drpq;
    float dr0[NDIM];
#else
    vector q0;
    real drpq2, drpq;
    vector dr0;
#endif
    INTEGER ipcount;

    int actlen;
    int *activenb;
    int nblist;
} gdhist_sincos_omp, *gdhistptr_sincos_omp;

typedef struct {
    real **xiOUTVPcos;
    real **xiOUTVPsin;
    real **xiOUTVPsincos;
    real **histZetaMtmpcos;
    real **histZetaMtmpsin;
    real **histZetaMtmpsincos;
    real *ChebsT;
    real *ChebsU;

    real ***histZetaMthreadcos;
    real ***histZetaMthreadsin;
    real ***histZetaMthreadsincos;
    realptr histNthread;
    realptr histNNSubthread;
// 2pcf
    realptr histNNSubXi2pcfthread;
//B kappa Avg Rmin
    realptr histNNSubXi2pcfthreadp;
    realptr histNNSubXi2pcfthreadtotal;
//E
//

    real **histXithreadcos;
    real **histXithreadsin;
    real *histXi2pcfthread;
    real *histXi2pcfthreadsub;

    realptr histNNNthread;
    real ***histNNNSubthread;
    real ***histXi3pcfthread;

#ifdef SINGLEP
    float q0[NDIM];
    float drpq2, drpq;
    float dr0[NDIM];
#else
    vector q0;
    real drpq2, drpq;
    vector dr0;
#endif
    INTEGER ipcount;

    int actlen;
    int *activenb;
    int nblist;
} gdhist_sincos_omp_3pcfbf, *gdhistptr_sincos_omp_3pcfbf;

#ifdef ADDONS
#include "datastruc_hist_include.h"
#endif

//
//E Structure definitions for histograms

#endif // ! _datastruc_hist_h

