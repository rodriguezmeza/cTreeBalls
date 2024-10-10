// Use:
//#include "global_data_kdtree_omp.h"

#ifndef _global_data_kdtree_omp_h
#define _global_data_kdtree_omp_h

#ifdef NMultipoles
// -----------------------------------
//B Histogram arrays
    realptr NhistNN;
    realptr NhistCF;
    realptr NhistNNSub;
// 2pcf
    realptr NhistNNSubXi2pcf;
//B kappa Avg Rmin
    realptr NhistNNSubXi2pcftotal;
//E
    real *NhistXi2pcf;
//B TPCF
    real ***NhistZetaMcos;
    real ***NhistZetaMsin;
    real ***NhistZetaMsincos;
// Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
    real ***NhistZetaMcossin;
//E
    real ***NhistZetaM;
//
    real *NhistNNN;
    real ***NhistNNNSub;
    real *NhistXi2pcf_omp;
    real ***NhistXi3pcf;
// TPCF
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
    
//B To save total 3pcf shear
    real *NhistXitt;
    real *NhistXixx;
    real *NhistXitx;
//E

//E Histograms arrays
// -----------------------------------
#endif // ! NMultipoles

#endif	// ! _global_data_kdtree_omp_h
