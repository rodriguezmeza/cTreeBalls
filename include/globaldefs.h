/*==============================================================================
 HEADER: globaldefs.h		[cTreeBalls]
 Written by: Mario A. Rodriguez-Meza
 Starting date: april 2023
 Purpose: Definitions of global variables and parameters
 Language: C
 Use: '#include "globaldefs.h"
 Major revisions:
 ==============================================================================*/
//        1          2          3          4        ^ 5          6          7

#ifndef _globaldefs_h
#define _globaldefs_h

//B ===============================================
#include <string.h>                                 // Incluido para quitar
                                                    //  el warning:
                                                    // "Implicit declaration of
                                                    //  built-in function 
                                                    //  'strcpy' y 'strchr'"
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>                               // To get time of the day

#include "stdinc.h"
#include "vectdefs.h"
#include "vectmath.h"
#include "getparam.h"
#include "mathfns.h"
#include "mathutil.h"
#include "inout.h"
#include "constant.h"



//B GSL section
#ifdef USEGSL

#ifndef NOINTERNALGSL
#include <stdio.h>
#include <math.h>
#include "config.h"
#include "gsl_math.h"
#include "gsl_rng.h"
#include "gsl_chebyshev.h"
#include "gsl_types.h"
#include "gsl_blas_types.h"
#include "gsl_complex.h"
#include "gsl_blas.h"
#include "gsl_complex_math.h"
#include "gsl_matrix.h"
#include "gsl_eigen.h"
#include "gsl_matrix_complex_double.h"
#include "gsl_errno.h"
#include "gsl_fft_real.h"
#include "gsl_fft_halfcomplex.h"
#include "gsl_fft_complex.h"
#include "gsl_integration.h"
#else // ! NOINTERNALGSL
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_chebyshev.h>
#include <gsl/gsl_types.h>
#include <gsl/gsl_blas_types.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_matrix_complex_double.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_integration.h>
#endif // ! NOINTERNALGSL

#else // ! USEGSL
#include "numrec.h"
#endif // ! USEGSL
//E

#include "datastruc_defs.h"

#ifdef USEGSL
global gsl_rng * r_gsl;                             // seed for random generators
#else
global long idum;                                   // seed for random generators
#endif

//B OpenMP section
//  Timing variables
#ifdef OPENMPCODE
#include "omp.h"
global double relstart,relend,absstart,absend;
#else
#include <time.h>
global time_t relstart,relend,absstart,absend;
#endif
//E

#include "common_defs.h"

//E ===============================================

#include "cmdline_data.h"
#include "global_data.h"

global bodyptr bodytable[MAXITEMS];
global nodeptr *nodetablescanlev[MAXITEMS];
global nodeptr *nodetablescanlev_root[MAXITEMS];
global cellptr roottable[MAXITEMS];

//global bodyptr bodytab;
global bodyptr bodytabbf;
global bodyptr bodytabsm;
global bodyptr bodytabSel;

global nodeptr *nodetab;

// BALLS
global nodeptr *nodetabscanlev;
// Root nodes:
global nodeptr *nodetabscanlev_root;
//

global real *histXi2pcf_omp;                        // Auxiliary array.
                                                    //  Used in OMP segments

//B Tree
//global cellptr root;
// BALLS
global cellptr rootnode;                            // To make treenodes
global bodyptr nodetable;                           // To smooth minimum size 
                                                    //  cells
global bodyptr nodetable_root;
//E

#include "datastruc_hist.h"

// To use in inout and cballsio
global real *inout_xval;
global real *inout_yval;
global real *inout_zval;
global real *inout_wval;

#ifdef ADDONS
#include "globaldefs_include_04.h"
#endif


//B CLASSLIB section
// standard libraries from Julien Lesgourges CLASS
#ifdef CLASSLIB
#include "common.h"
global ErrorMsg errmsg;
#else // ! CLASSLIB

#if !defined(TRUE)
#define TRUE 1                                      // integer associated to true
                                                    //  statement
#endif
#if !defined(FALSE)
#define FALSE 0                                     // integer associated to false
                                                    //  statement
#endif
#define SUCCESS 0                                   // integer returned after successful
                                                    //  call of a function
#define FAILURE 1                                   // integer returned after failure in
                                                    //  a function
#define _ERRORMSGSIZE_ 2048                         // generic error messages are cut
                                                    //  beyond this number of characters
typedef char ErrorMsg[_ERRORMSGSIZE_];              // Generic error messages
#ifndef MIN
#define MIN(a,b) (((a)<(b)) ? (a) : (b) )           // the usual "min" function
#endif
#ifndef MAX
#define MAX(a,b) (((a)<(b)) ? (b) : (a) )           // the usual "max" function
#endif
#define SIGN(a) (((a)>0) ? 1. : -1. )
#define NRSIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
global ErrorMsg errmsg;
#endif // ! CLASSLIB
//E


#ifdef ADDONS
#include "globaldefs_include_05.h"
#endif

#include "protodefs.h"

#endif // ! _globaldefs_h

