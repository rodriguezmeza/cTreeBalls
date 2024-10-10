// Use:
//#include "global_data_octree_kkk_omp.h"

#ifndef _global_data_octree_kkk_omp_h
#define _global_data_octree_kkk_omp_h

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

#endif	// ! _global_data_octree_kkk_omp_h
