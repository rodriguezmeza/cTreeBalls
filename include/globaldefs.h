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

//
// lines where there is a "//B socket:" string are places to include module files
//  that can be found in addons/addons_include folder
//

#ifndef _globaldefs_h
#define _globaldefs_h

//B ===============================================
//B general includes section
#include <string.h>                                 // Incluido para quitar
                                                    //  el warning:
                                                    // "Implicit declaration of
                                                    //  built-in function 
                                                    //  'strcpy' y 'strchr'"
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>                               // To get time of the day

//B
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
//E

#include "stdinc.h"
#include "vectdefs.h"
#include "vectmath.h"
#ifdef GETPARAM
#include "getparam.h"
#endif
#include "mathfns.h"
#include "mathutil.h"
#include "inout.h"
#include "constant.h"

#include "common_defs.h"


//B GSL section
//  so far only needed to compute edge-corrections and
//  testsmodels
//  however, generating test models can be done with NR randoms
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
#include "gsl_linalg.h"                             // Extract from sources...
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
#include <gsl/gsl_linalg.h>
#endif // ! NOINTERNALGSL

#else // ! USEGSL
#include "numrec.h"
#endif // ! USEGSL
//E GSL section

//B OpenMP section
#ifdef OPENMPCODE
#include "omp.h"
#endif
//E

#if !defined(global)
#  define global extern
#endif

#if !defined(global)
#  define global extern
#endif

#ifdef USEGSL
global gsl_rng * r_gsl;
#else
global long idum;
#endif

#include "datastruc_defs.h"                         // node, body, cell...
//E general includes section
//E ===============================================

#ifdef CLASSLIB
#include "common.h"
#endif
global ErrorMsg errmsg;

#define _ERRORMSGSIZE_ 2048
typedef char ErrorMsg[_ERRORMSGSIZE_];

#include "cmdline_data.h"
#include "options_cache.h"
#include "global_data.h"

global bodyptr bodytable[MAXITEMS];
global nodeptr *nodetablescanlev[MAXITEMS];
global nodeptr *nodetablescanlev_root[MAXITEMS];
global cellptr roottable[MAXITEMS];

#ifdef CBALLS_NEEDS_BALLS4_SCAN
global nodeptr *nodetablescanlevB4[MAXITEMS];
#endif

#ifndef MACONLY
//B celltable
global cellptr *celltable[MAXITEMS];
//E
#endif

// check this... it is repeated in global_data struct
global real *histXi2pcf_omp;                        // Auxiliary array.
                                                    //  Used in OMP segments
//B Tree
// BALLS
global cellptr rootnode;                            // To make treenodes
//E


#include "datastruc_hist.h"

// To use in inout and cballsio
global real *inout_xval;
global real *inout_yval;
global real *inout_zval;
global real *inout_uval;
global real *inout_vval;
global real *inout_wval;

//B socket:
#ifdef ADDONS
#include "globaldefs_include_04.h"
#endif
//E

//B socket:
#ifdef ADDONS
#include "globaldefs_include_05.h"
#endif
//E

/*
 * Mutable process globals are retained for the standalone executable and for
 * the existing OpenMP shared clauses.  The Python wrapper switches these
 * values through one runtime state per cballs object before entering C.
 */
#if ((defined(BALLS) && !defined(OCTREESMOOTHING)) || \
     defined(OCTREESMOOTHING) || defined(BALLS0357))
#define CBALLS_RUNTIME_BALLS_GLOBALS 1
#endif

global bool tree_is_threaded[MAXITEMS];

typedef struct cballs_runtime_state {
#ifdef USEGSL
    gsl_rng *r_gsl;
    mMatrix_ptr histZetaMatrix;
#else
    long idum;
#endif
    ErrorMsg errmsg;
    bodyptr bodytable[MAXITEMS];
    nodeptr *nodetablescanlev[MAXITEMS];
    nodeptr *nodetablescanlev_root[MAXITEMS];
    cellptr roottable[MAXITEMS];
#ifdef CBALLS_NEEDS_BALLS4_SCAN
    nodeptr *nodetablescanlevB4[MAXITEMS];
#endif
#ifndef MACONLY
    cellptr *celltable[MAXITEMS];
#endif
    bool tree_is_threaded[MAXITEMS];
    real *histXi2pcf_omp;
    cellptr rootnode;
    real *inout_xval;
    real *inout_yval;
    real *inout_zval;
    real *inout_uval;
    real *inout_vval;
    real *inout_wval;
#ifdef CBALLS_RUNTIME_BALLS_GLOBALS
    bodyptr bodytabbf;
    bodyptr bodytabsm;
    bodyptr bodytabSel;
    nodeptr *nodetab;
    nodeptr *nodetabscanlev;
    nodeptr *nodetabscanlev_root;
    bodyptr nodetable;
    bodyptr nodetable_root;
#endif
} cballs_runtime_state;

global cballs_runtime_state *cballs_runtime_create(void);
global int cballs_runtime_activate(cballs_runtime_state *state);
global const void *cballs_runtime_bodytable_at(
    const cballs_runtime_state *state, int ifile);
global void cballs_runtime_destroy(cballs_runtime_state *state);

#include "protodefs.h"
#ifdef CBALLS_MPI_ENABLED
#include "cballs_mpi_dispatch.h"
#endif

#endif // ! _globaldefs_h

