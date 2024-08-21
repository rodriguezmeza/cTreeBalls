# Set of declarations from C to python.

#B This is to use with USEGSL=1 and GSLINTERNAL=0
#cdef extern from "gsl/gsl_rng.h":
#    ctypedef struct gsl_rng_type
#    ctypedef struct gsl_rng

#    cdef gsl_rng_type *gsl_rng_default
#    unsigned long int gsl_rng_default_seed
#E

DEF MAXITEMS = 100
DEF MAXLENGTHOFFILES = 1024
DEF MAXLEVEL = 32
# Set it as is done in stdinc.h
DEF NDIM = 3

ctypedef double real
#ctypedef double* realptr
ctypedef long INTEGER
ctypedef char* string
# Correct to include FILE ctype
#ctypedef FILE *stream
ctypedef short int bool
ctypedef INTEGER vectorI[NDIM];
ctypedef double vector[NDIM];



cdef extern from "cballs.h":

    ctypedef char FileArg[40]

    ctypedef char* ErrorMsg


#B cmdline structure
    cdef struct cmdline_data:
# Every item in cmdline_defs.h must have an item here::

        INTEGER stepState

        string script
        string options
        string version
        short verbose
        short verbose_log
        string paramfile

#ifndef GETPARAM
        char ParameterFile[MAXLENGTHOFFILES]           # May be we should incrase
                                                    #  this number
#endif

        string searchMethod

        double theta

        int mchebyshev
#       short dimension;

#B Parameters to set a region in the sky, for example for Takahasi data set.
        real thetaL
        real thetaR
        real phiL
        real phiR
#E

        int sizeHistN
        real rangeN
        real rminHist
#       bool logHist

        string infile
        string infilefmt

        string iCatalogs

        string rootDir
        string outfile
        string outfilefmt


        string histNFileName
        string histXi2pcfFileName
        string histZetaMFileName
        string mhistZetaFileName
        string suffixOutFiles
    
#ifdef OPENMPCODE
        int numthreads
#endif

        int seed
        string nsmooth

        string rsmooth

#ifdef ADDONS
#include "globaldefs_include_01.h"
#endif

#    INTEGER stepNodes
#    string ncritical

        string testmodel
        INTEGER nbody
        real lengthBox

# Not used, delet it:
        int mToPlot


#E cmdline structure


#B globaldata structure
    cdef struct global_data:
        real cpuinit
# Correct to include this structure in cpython...
#        struct timeval current_time
        long cpurealinit                               # get time of the day

        string headline0
        string headline1
        string headline2
        string headline3

# Correct to include FILE ctype
#        FILE *outlog
# Is it used? (delete)
#        FILE *outstr_sols

        char mode[2]

        int searchmethod_int

#B Settings of the code
#       short dimension                                # Set dimension of the run
#       bool tpcfon                                    # Set 3pcf computation on
#E

#B Tree
        INTEGER ncell
        int tdepth
        INTEGER actmax
#B BALLS
        INTEGER ncccalc
# Correct to include FILE ctype
#        FILE *outnodelev
#        FILE *outbodylev
        INTEGER nsmoothcount
#E
        INTEGER nbccalc
        INTEGER nbbcalc
        real rSize                                     # Maximum r of the box

        real rSizeTable[MAXITEMS]                      # Maximum r of the box
        INTEGER ncellTable[MAXITEMS]
        INTEGER nbodyTable[MAXITEMS]
        int tdepthTable[MAXITEMS]
        INTEGER nnodescanlevTable[MAXITEMS]
        INTEGER nnodescanlev_rootTable[MAXITEMS]

# Set it using the struct node, body and cell definition
#        cellptr root
#E

        real cputree

# Tree:
        real Rcut                                     # Cutoff radius
        real RcutSq
        real cpusearch
#
# Cell search
        vectorI cells
        INTEGER *cellList
#

        int infilefmt_int

##ifndef NOGSL
##ifdef USEGSL
# Set it using cython version of GSL
#        gsl_rng * r
##endif

#ifdef SINGLEP
#        double Box[NDIM]
#else
        vector Box
#endif

#B Histogram arrays
        double* histN
        double* histCF
        double* histNSub
# 2pcf
        double* histNSubXi2pcf
#B kappa Avg Rmin
#    double* histNSubXi2pcfRmin
        double* histNSubXi2pcftotal
#E
        real *histXi2pcf
        real ***histZetaMcos
        real ***histZetaMsin
        real ***histZetaMsincos
#E Histogram arrays
        real ***histZetaM
#
        real *histNNN
        real ***histNNNSub
        real *histXi2pcf_omp
        real ***histXi3pcf
        real **histXi
        real **histXicos
        real **histXisin
##ifndef NOGSL
#ifdef USEGSL
# Set it using cython version of GSL
#        gsl_matrix_complex *histXi_gsl
#endif
        real rCutSq
        int searchMethod_int

        real deltaR
        real deltaRmin
        real deltaRmax
#ifdef LOGHIST
        real *deltaRV
        real *ddeltaRV
#endif

        real rrRange;
        real deltaTheta;

# Activate later...
#    string histCFFileName

        char logfilePath[MAXLENGTHOFFILES]
        char outputDir[MAXLENGTHOFFILES]
        char tmpDir[MAXLENGTHOFFILES]

        char fpfnameOutputFileName[MAXLENGTHOFFILES]
        char fpfnamehistNFileName[MAXLENGTHOFFILES]
        char fpfnamehistCFFileName[MAXLENGTHOFFILES]
        char fpfnamehistrBinsFileName[MAXLENGTHOFFILES]
        char fpfnamehistXi2pcfFileName[MAXLENGTHOFFILES]
        char fpfnamehistZetaMFileName[MAXLENGTHOFFILES]
        char fpfnamemhistZetaFileName[MAXLENGTHOFFILES]
        char fpfnameCPUFileName[MAXLENGTHOFFILES]

        string model_comment

        int stopflag
        INTEGER ip
        real cputotalinout
        real cputotal

        INTEGER bytes_tot
        INTEGER bytes_tot_cells

#B To see the bodies belonging to a cell:
        INTEGER nbodySel
#E
        INTEGER nbodysm
        INTEGER nbodybf

        bool bh86
        bool sw94

        real i_deltaR

#ifdef ADDONS
#include "globaldefs_include_02.h"
#endif

        char fnameData_kd[128]
        char fnameOut_kd[128]
        int input_format_kd
        int use_tree_kd
        int max_tree_order_kd
        int max_tree_nparts_kd
        int use_pm_kd
        INTEGER n_objects_kd
        float l_box_kd
        float l_box_half_kd

        int ninfiles
        char *infilenames[MAXITEMS]
        char *infilefmtname[MAXITEMS]
        int iCatalogs[MAXITEMS]

        int nsmooth[MAXITEMS]
        INTEGER nnode
        INTEGER rnnode
#E
##ifdef BALLS
#    real scanLevelMin[MAXITEMS]
        int scanLevelMin[MAXITEMS]
#    int ncritical[MAXITEMS]

        char nodesfilePath[MAXLENGTHOFFILES]
        int nnodescanlev
# Root nodes:
        int nnodescanlev_root
#define MAXLEVEL  32
        real Rcell[MAXLEVEL]
#undef MAXLEVEL
#
        char bodiesfilePath[MAXLENGTHOFFILES]
##endif

        bool flagSmoothCellMin
        bool flagSmooth
        bool flagSetNbNoSel
#B BUCKET
        real rminCell[2]
#E

        real rsmooth[MAXITEMS]
        bool rsmoothFlag

        int irsmooth

#ifdef ADDONS
#include "globaldefs_include_03.h"
#endif


#E globaldata structure



    cdef struct file_content:
        char * filename
        int size
        FileArg * name
        FileArg * value
        short * read

    cdef int FAILURE
    cdef int FALSE
    cdef int TRUE

    int input_read_from_file(void*, void*, char*)
    int StartRun_Common(void*, void*);
    int PrintParameterFile(void*, char *)

    int StartOutput(void*)
    int SetNumberThreads(void*)

    int MainLoop(void*, void*);
    int EndRun(void*, void*);

    int saveHistZetaM_sincos(void*, void*, int, int, double***);

