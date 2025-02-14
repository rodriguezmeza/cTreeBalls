// Use:
//#include "datastruc_hist_kdtree_omp.h"

// NMultipoles has been switched off for kdtree_omp
//  NMultipoles -> NMultipoles_kdtree

#ifndef _datastruc_hist_kdtree_omp_h
#define _datastruc_hist_kdtree_omp_h

/*
#ifdef NMultipoles_kdtree
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
} gdhistN_sincos_omp, *gdhistNptr_sincos_omp;
#endif // ! NMultipoles
*/

#endif	// ! _datastruc_hist_kdtree_omp_h
