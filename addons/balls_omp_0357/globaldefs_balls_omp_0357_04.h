// Use:
//#include "globaldefs_balls_omp_0357_04.h"
//
//
// included in file: addons/addons_include/include/datastruc_hist_include.h
//

#ifndef _globaldefs_balls_omp_0357_04_h
#define _globaldefs_balls_omp_0357_04_h

typedef struct {
    real **xiOUTVP;
    real **histZetaMtmp;
    real *Chebs;

    real ***histZetaMthread;
    realptr histNthread;
    realptr histNNSubthread;
// 2pcf
    realptr histNNSubXi2pcfthread;
    //B kappa Avg Rmin
    realptr histNNSubXi2pcfthreadp;
    realptr histNNSubXi2pcfthreadtotal;
    //E
//
    real **histXithread;
    real *histXi2pcfthread;
    real *histXi2pcfthreadsub;

    compute_vector q0;
    real drpq2, drpq;
    compute_vector dr0;
    INTEGER ipcount;

    int actlen;
    int *activenb;
    int nblist;

    nodeptr *active;
    cellptr interact;

} gdhist_omp_balls, *gdhistptr_omp_balls;

typedef struct {
    real **xiOUTVP;
    real **histZetaMtmp;
    real *Chebs;
    real ***histZetaMthread;
    realptr histNthread;
    realptr histNNSubthread;
// 2pcf
    realptr histNNSubXi2pcfthread;
//B kappa Avg Rmin
    realptr histNNSubXi2pcfthreadp;
    realptr histNNSubXi2pcfthreadtotal;
//E
//
    real **histXithread;
    real *histXi2pcfthread;
    real *histXi2pcfthreadsub;

    compute_vector q0;
    real drpq2, drpq;
    compute_vector dr0;
    INTEGER ipcount;

    int actlen;
    int *activenb;
    int nblist;
    
    nodeptr *active;
    cellptr interact;

} gdhist_omp_balls6, *gdhistptr_omp_balls6;

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

    compute_vector q0;
    real drpq2, drpq;
    compute_vector dr0;
    INTEGER ipcount;

    int actlen;
    int *activenb;
    int nblist;
    
    nodeptr *active;
    cellptr interact;

} gdhist_sincos_omp_balls6, *gdhistptr_sincos_omp_balls6;

#endif	// ! _globaldefs_balls_omp_0357_04_h
