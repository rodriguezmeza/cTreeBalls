/*==============================================================================
 HEADER: cmdline_data.h		[cTreeBalls]
 Written by: Mario A. Rodriguez-Meza
 Starting date: april 2023
 Purpose: Definitions of global variables and parameters
 Language: C
 Use: '#include "datastruc_hist.h"
 Major revisions:
 ==============================================================================*/
//        1          2          3          4          5          6          7

#ifndef _datastruc_hist_h
#define _datastruc_hist_h

//B Structure definitions for histograms
//

//#ifndef NOGSL
#ifdef USEGSL
typedef struct {
    int m;
    gsl_matrix_complex *histZetaM;
} mMatrix, *mMatrix_ptr;

global mMatrix_ptr histZetaMatrix;
#endif

typedef struct {
    real **xiOUTVP;
    real **histZetaMtmp;
    real *histXi2pcfsub;
    real *Chebs;

    vector q0;
    real drpq2, drpq;
    vector dr0;
    INTEGER ipcount;

} gdhist, *gdhistptr;

typedef struct {
// TPCF
    real **xiOUTVP;
    real **histZetaMtmp;
    real *histXi2pcfsub;
    real *Chebs;

    vector q0;
    real drpq2, drpq;
    vector dr0;
    INTEGER ipcount;

    int actlen;
    int *activenb;
    int nblist;
    
    nodeptr *active;
    cellptr interact;

} gdhist_balls, *gdhistptr_balls;

#ifdef MPICODE
typedef struct {
    real **xiOUTVPlocal;
    real **histZetaMtmplocal;
    real ***histZetaMlocal;
    real **histXilocal;
    real *histXi2pcflocalsub;
    realptr histNlocal;
    realptr histNSublocal;
    realptr histXi2pcflocal;
} gdhist_mpi, *gdhistptr_mpi;
#endif


typedef struct {
    real **xiOUTVP;
    real **histZetaMtmp;
    real *Chebs;
    real ***histZetaMthread;
    realptr histNthread;
    realptr histNSubthread;
// 2pcf
    realptr histNSubXi2pcfthread;
//
    real **histXithread;
    real *histXi2pcfthread;
    real *histXi2pcfthreadsub;

    vector q0;
    real drpq2, drpq;
    vector dr0;
    INTEGER ipcount;

    int actlen;
    int *activenb;
    int nblist;

//B BALLS4
    // PASAR ESTOS CAMBIOS A gdhist_omp_balls
    nodeptr *active;
    cellptr interact;
//E
} gdhist_omp, *gdhistptr_omp;


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
    realptr histNSubthread;
// 2pcf
    realptr histNSubXi2pcfthread;
//B kappa Avg Rmin
    realptr histNSubXi2pcfthreadp;
    realptr histNSubXi2pcfthreadtotal;
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
    real *Chebs;
    realptr histNthread;
    realptr histNSubthread;
    real *histXi2pcfthread;
    real *histXi2pcfthreadsub;
    
    realptr histNNNthread;
    real ***histNNNSubthread;
    real ***histXi3pcfthread;

    vector q0;
    real drpq2, drpq;
    vector dr0;
    INTEGER ipcount;

} gdhist_omp_3pcfbf, *gdhistptr_omp_3pcfbf;


typedef struct {
    real *Chebs;
    real *histXi2pcfsub;
    
    realptr histNNN;
    real ***histNNNSub;
    real ***histXi3pcf;

    vector q0;
    real drpq2, drpq;
    vector dr0;
    INTEGER ipcount;

} gdhist_3pcfbf, *gdhistptr_3pcfbf;
//
//E Structure definitions for histograms


#endif // ! _datastruc_hist_h

