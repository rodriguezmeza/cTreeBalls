/*==============================================================================
 HEADER: globaldefs.h		[tpcf]
 Written by: Mario A. Rodriguez-Meza
 Starting date: april 2023
 Purpose: Definitions of global variables and parameters
 Language: C
 Use: '#include "global_defs.h"
 Major revisions:
 ==============================================================================*/

#ifndef _globaldefs_h
#define _globaldefs_h

//===============================================
#include <string.h>     // Incluido para quitar el warning:
                        // "Implicit declaration of built-in function 'strcpy' y 'strchr'"
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>       // To get time of the day

#include "stdinc.h"
#include "vectdefs.h"
#include "vectmath.h"
#include "getparam.h"
#include "mathfns.h"
#include "mathutil.h"
#include "inout.h"
#include "constant.h"

/*
//B CLASSLIB section
// Also in globaldefs_01.h, need to fix this part
// standard libraries from Julien Lesgourges CLASS
#ifdef CLASSLIB
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "string.h"
#include "float.h"
#ifdef _OPENMP
#include "omp.h"
#endif
#include "common.h"
#else
#define _SUCCESS_ 0 // integer returned after successful call of a function
#define _FAILURE_ 1 // integer returned after failure in a function
#define _ERRORMSGSIZE_ 2048 // generic error messages are cut beyond this number of characters
typedef char ErrorMsg[_ERRORMSGSIZE_]; // Generic error messages

#define MIN(a,b) (((a)<(b)) ? (a) : (b) ) // **< the usual "min" function
#define MAX(a,b) (((a)<(b)) ? (b) : (a) ) // **< the usual "max" function
#define SIGN(a) (((a)>0) ? 1. : -1. )
#define NRSIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define index_symmetric_matrix(i1,i2,N) (((i1)<=(i2)) ? ((i2)+N*(i1)-((i1)*((i1)+1))/2) : ((i1)+N*(i2)-((i2)*((i2)+1))/2)) //< assigns an index from 0 to [N(N+1)/2-1] to the coefficients M_{i1,i2} of an N*N symmetric matrix; useful for converting a symmetric matrix to a vector, without losing or double-counting any information

#endif
//E
*/

/*
//B OpenMP section
//  Timing variables
#ifdef OPENMPCODE
#include "omp.h"
//static double relstart,relend,absstart,absend;
global double relstart,relend,absstart,absend;
#else
#include <time.h>
//static time_t relstart,relend,absstart,absend;
global time_t relstart,relend,absstart,absend;
#endif
//E
*/

//B GSL section
#ifndef NOGSL

#include <stdio.h>
#include <math.h>
#ifdef GSLINTER
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
#include "gsl_integration.h"
#else
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
#include <gsl/gsl_integration.h>
#endif

#endif
//E

#include "datastruc_defs.h"

//B OpenMP section
//  Timing variables
#ifdef OPENMPCODE
#include "omp.h"
//static double relstart,relend,absstart,absend;
global double relstart,relend,absstart,absend;
#else
#include <time.h>
//static time_t relstart,relend,absstart,absend;
global time_t relstart,relend,absstart,absend;
#endif
//E

//B Usefule MACROS and other constants:
//#define CPUTIME         (cputime())
#define CPUTIME         (second()/60.0)         // Gives cputime in minutes
#define REALTIME        (gettimeofday(&gd.current_time, NULL)) // Gives time of the day in minutes
#define PRNUNITOFTIMEUSED   "min."

#define MAXITEMS    100

//#define MAXLENGTHOFFILES       200
#define MAXLENGTHOFFILES       300
#define MAXLENGTHOFSTRSCMD     200
#define EXTFILES            ".txt"
#define INMB                9.536743116E-7      // 1/(1024*1024)
// Already in Makefile_settings. Set it there
//#define NLOGBINPD   5       // bins per decade for log histograms
//E

//B
#ifdef MPICODE
#include <mpi.h>
// Get the rank of the process
//int world_rank;
global int ThisTask;
// Get the number of processes
//int world_size;
global int NTask;
global int PTask;

// Get the name of the processor
global char processor_name[MPI_MAX_PROCESSOR_NAME];
global int name_len;
#endif
//E

//===============================================

typedef struct {
// Every item in cmdline_defs.h must have an item here::

    string statefile;
    INTEGER stepState;
    string restorefile;
    string options;
    string version;
    short verbose;
    short verbose_log;
	string paramfile;

#ifndef GETPARAM
    char ParameterFile[MAXLENGTHOFFILES]; // May be you should incrase this number
#endif

    string searchMethod;

    real theta;

    int mchebyshev;
//    short dimension;

//B Parameters to set a region in the sky, for example for Takahasi data set.
    real thetaL;
    real thetaR;
    real phiL;
    real phiR;
//E

    int sizeHistN;
    real rangeN;
    real rminHist;
//    bool logHist;
    int sizeHistTheta;

    string infile;
    string infilefmt;
    string rootDir;
    string outfile;
    string outfilefmt;
    

    string histNFileName;
    string histXi2pcfFileName;
    string histZetaMFileName;
    string mhistZetaFileName;
    string suffixOutFiles;
    
    int numthreads;

//
//B NOLSST:
//
    int seed;
    string nsmooth;
#ifdef BALLS
    INTEGER ntosave;
    int scanLevel;
#endif
    INTEGER stepNodes;
    string ncritical;
//#endif
    string testmodel;
    INTEGER nbody;
    real lengthBox;
    int mToPlot;
//
//E
//

} cmdline_data, *cmdline_data_ptr;

typedef struct {
    real cpuinit;
    struct timeval current_time;
    long cpurealinit;       // get time of the day

	string headline0;
	string headline1;
	string headline2;
	string headline3;

    FILE *outlog;
//B To debug cells:
    FILE *outcells;
//E
	FILE *outstr_sols;

	char mode[2];

    int searchmethod_int;

//B Tree
    INTEGER ncell;
    int tdepth;
    INTEGER actmax;
// BALLS
    INTEGER ncccalc;
    FILE *outnodelev;
    FILE *outbodylev;
    INTEGER nsmoothcount;
//
    INTEGER nbccalc;
    INTEGER nbbcalc;
    real rSize;             // Maximum r of the box
    cellptr root;
//E

    real cputree;

// Tree:
     real Rcut;            // Cutoff radius
     real RcutSq;
    real cpusearch;
//
// Cell search
     vectorI cells;
    INTEGER *cellList;
//

    int infilefmt_int;
    
#ifndef NOGSL
    gsl_rng * r;
#endif

    vector Box;
    realptr histN;
    realptr histCF;
    realptr histNSub;
// 2pcf
    realptr histNSubXi2pcf;
//
    real *histNNN;
    real ***histNNNSub;
    real *histXi2pcf;
    real *histXi2pcf_omp;
    real ***histXi3pcf;
    real **histXi;
    real **histXicos;
    real **histXisin;
#ifndef NOGSL
    gsl_matrix_complex *histXi_gsl;
#endif
    real ***histZetaM;
    real ***histZetaMcos;
    real ***histZetaMsin;
    real ***histZetaMsincos;
    real rCutSq;
    int searchMethod_int;

    real deltaR;
#ifdef LOGHIST
    real *deltaRV;
#endif

    real rrRange;
    real deltaTheta;

// Activate later...
//    string histCFFileName;

    char logfilePath[MAXLENGTHOFFILES];
    char outputDir[MAXLENGTHOFFILES];
    char tmpDir[MAXLENGTHOFFILES];

    char fpfnameOutputFileName[MAXLENGTHOFFILES];
    char fpfnamehistNFileName[MAXLENGTHOFFILES];
    char fpfnamehistCFFileName[MAXLENGTHOFFILES];
    char fpfnamehistrBinsFileName[MAXLENGTHOFFILES];
    char fpfnamehistXi2pcfFileName[MAXLENGTHOFFILES];
    char fpfnamehistZetaMFileName[MAXLENGTHOFFILES];
    char fpfnamemhistZetaFileName[MAXLENGTHOFFILES];
    char fpfnameCPUFileName[MAXLENGTHOFFILES];

    string model_comment;
    
    int stopflag;
    INTEGER ip;
    real cputotalinout;
    real cputotal;

    INTEGER bytes_tot;
    INTEGER bytes_tot_cells;

//B To see the bodies belonging to a cell:
    INTEGER nbodySel;
//E
    INTEGER nbodysm;
    INTEGER nbodybf;

    bool bh86, sw94;

    real i_deltaR;


    char fnameData_kd[128];
    char fnameOut_kd[128];
    int input_format_kd;
    int use_tree_kd;
    int max_tree_order_kd;
    int max_tree_nparts_kd;
    int use_pm_kd;
    INTEGER n_objects_kd;
    float l_box_kd;
    float l_box_half_kd;

    int ninfiles;
    char *infilenames[MAXITEMS];
    char *infilefmtname[MAXITEMS];
    int nsmooth[MAXITEMS];
    //B To debug cells:
    char cellsfilePath[MAXLENGTHOFFILES];
    INTEGER nnode;
    INTEGER rnnode;
    //E
//#ifdef BALLS
    int ncritical[MAXITEMS];
    
    char nodesfilePath[MAXLENGTHOFFILES];
    int nnodescanlev;
    char bodiesfilePath[MAXLENGTHOFFILES];
//#endif
    
    bool flagSmooth;
    bool flagSetNbNoSel;

} global_data, *global_data_ptr;

global global_data gd;
global cmdline_data cmd;

global bodyptr bodytab;
global bodyptr bodytabbf;
global bodyptr bodytabsm;
global bodyptr bodytabSel;

//B To debug cells:
//global cellptr nodetab;
//global cellptr *nodetab;
global nodeptr *nodetab;
//local nodeptr *active;
//local cellptr interact;
//E

// BALLS
global nodeptr *nodetabscanlev;
//

global real *histXi2pcf_omp;        // Auxiliary array.
                                    // Used in OMP segments

//B Tree
global cellptr root;
// BALLS
global cellptr rootnode;
global bodyptr nodetable;
//E

#ifndef NOGSL
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

#ifdef OPENMPCODE

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
    
    nodeptr *active;
    cellptr interact;

} gdhist_omp_barnes, *gdhistptr_omp_barnes;

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
//
    real **histXithreadcos;
    real **histXithreadsin;
    real *histXi2pcfthread;
    real *histXi2pcfthreadsub;

    vector q0;
    real drpq2, drpq;
    vector dr0;
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

    nodeptr *active;
    cellptr interact;

} gdhist_omp_balls, *gdhistptr_omp_balls;

#endif

typedef struct {
    real *Chebs;
//    realptr histN;
//    realptr histNSub;
//    real *histXi2pcf;
    real *histXi2pcfsub;
    
    realptr histNNN;
    real ***histNNNSub;
    real ***histXi3pcf;

    vector q0;
    real drpq2, drpq;
    vector dr0;
    INTEGER ipcount;

} gdhist_3pcfbf, *gdhistptr_3pcfbf;

// From inout.h
global real *inout_xval;
global real *inout_yval;
global real *inout_zval;
global real *inout_wval;


//NOLSST:
// It is already in globaldef.h, need to fix this
//B CLASSLIB section
// standard libraries from Julien Lesgourges CLASS
#ifdef CLASSLIB
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "string.h"
#include "float.h"
#ifdef _OPENMP
#include "omp.h"
#endif
#include "common.h"
#else
#define _TRUE_ 1 //< integer associated to true statement
#define _FALSE_ 0 // **< integer associated to false statement
#define _SUCCESS_ 0 // integer returned after successful call of a function
#define _FAILURE_ 1 // integer returned after failure in a function
#define _ERRORMSGSIZE_ 2048 // generic error messages are cut beyond this number of characters
typedef char ErrorMsg[_ERRORMSGSIZE_]; // Generic error messages

#ifndef MIN
#define MIN(a,b) (((a)<(b)) ? (a) : (b) ) // **< the usual "min" function
#endif
#ifndef MAX
#define MAX(a,b) (((a)<(b)) ? (b) : (a) ) // **< the usual "max" function
#endif
#define SIGN(a) (((a)>0) ? 1. : -1. )
#define NRSIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define index_symmetric_matrix(i1,i2,N) (((i1)<=(i2)) ? ((i2)+N*(i1)-((i1)*((i1)+1))/2) : ((i1)+N*(i2)-((i2)*((i2)+1))/2)) // < assigns an index from 0 to [N(N+1)/2-1] to the coefficients M_{i1,i2} of an N*N symmetric matrix; useful for converting a symmetric matrix to a vector, without losing or double-counting any information
global ErrorMsg errmsg;
#endif
//E


#include "protodefs.h"

#endif // ! _globaldefs_h

