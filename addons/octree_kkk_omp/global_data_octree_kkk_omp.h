// Use:
//#include "global_data_octree_kkk_omp.h"

#ifndef _global_data_octree_kkk_omp_h
#define _global_data_octree_kkk_omp_h

/*
// All of these are in addons_include/include/global_data_include_NMultipoles.h
//  CAN BE DELETED ALL THESE
#ifdef NMultipoles
// -----------------------------------
//B Histogram arrays
    realptr NhistNNSub;

    real ***NhistZetaMcos;
    real ***NhistZetaMsin;
    real ***NhistZetaMsincos;
    real ***NhistZetaMcossin;
    real ***NhistZetaM;
    real **NhistXi;
    real **NhistXicos;
    real **NhistXisin;
#ifdef USEGSL
    gsl_matrix_complex *NhistXi_gsl;
#endif

//B To save total 3pcf
    real **NhistZetaGcos;
    real **NhistZetaGsin;
    real ***NhistZetaGmRe;
    real ***NhistZetaGmIm;
//E
    
//E Histograms arrays
// -----------------------------------
#endif // ! NMultipoles
*/

// 2pcf
    realptr histNNSubN2pcf;
//B kappa Avg Rmin
    realptr histNNSubN2pcftotal;
//E
    real *histN2pcf;

    char fpfnamehistN2pcfFileName[MAXLENGTHOFFILES];


#endif	// ! _global_data_octree_kkk_omp_h
