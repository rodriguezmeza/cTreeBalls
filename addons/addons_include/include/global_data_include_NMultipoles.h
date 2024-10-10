// Use:
//#include "global_data_include_NMultipoles.h"

#ifndef _global_data_include_NMultipoles_h
#define _global_data_include_NMultipoles_h

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

//B Additional definitions to be use in kdtree_omp... delete ASAP
realptr NhistNN;
realptr NhistCF;

// 2pcf
    realptr NhistNNSubXi2pcf;
//B kappa Avg Rmin
    realptr NhistNNSubXi2pcftotal;
//E
    real *NhistXi2pcf;
//B TPCF

real *NhistNNN;
real ***NhistNNNSub;
real *NhistXi2pcf_omp;
real ***NhistXi3pcf;

//B To save total 3pcf shear
    real *NhistXitt;
    real *NhistXixx;
    real *NhistXitx;
//E

//E

#endif // ! NMultipoles

/*
 Add your addon item here
 */

#endif	// ! _global_data_include_NMultipoles_h
