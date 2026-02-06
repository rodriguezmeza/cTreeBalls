# Set of declarations from C to python
# Only define quantities you will use,
# for input, output or intermediate manipulations
# by the python wrapper. The rest is internal in cBalls

cdef extern from "cballs.h":

    ctypedef char FileArg[40]
    ctypedef char* ErrorMsg

    cdef struct cmdline_data:
        char * searchMethod;
        int mChebyshev;
        int nsmooth;
        char * rsmooth;
        double theta;
        short computeTPCF;
        short usePeriodic;

        int sizeHistN
        int numberThreads
        int verbose
        int verbose_log
        double rangeN
        double rminHist
        char * preScript
        char * posScript
        char * rootDir
        char * iCatalogs
        char * columns
        char * options

    cdef struct global_data:
        short cmd_allocated
        short gd_allocated
        short gd_allocated_2
        short random_allocated
        short histograms_allocated
        short tree_allocated
        short bodytable_allocated
#B PXD functions
        double *vecPXD
        double **matPXD
        double *rBins
        double *histNN
        double *histNNPXD
        double *histZetaMFlatten
        double ***histZetaMcos
#E
        double Rcut
        double* histN
        double* histCF
        double *histXi2pcf
        double ***histZetaM
        double *histNNN
        double **histXi
        double **histXicos
        double **histXisin

# In common_defs.h #define MAXITEMS                100
        int columns[100]
        double rsmooth[100]
        short rsmoothFlag;

    cdef struct file_content:
        char * filename
        int size
        FileArg * name
        FileArg * value
        short * read

    cdef int FAILURE
    cdef int FALSE
    cdef int TRUE

    int input_read_from_file(void*, void*, void*, char*)
    int StartRun_Common(void*, void*);
    int PrintParameterFile(void*, void*, char*)

    int StartOutput(void*)
    int SetNumberThreads(void*)

    int MainLoop(void*, void*);
    int EndRun(void*, void*);
    int EndRun_FreeMemory(void*, void*);
    int freeTree(void*, void*);

#B parameters
    int get_cmdline(void*, void*)
    int get_nthreads(void*, int*)
    int get_nmonopoles(void*, int*)
    int get_theta(void*, double*)
    int get_rsmooth(void*, double*)
    int get_cputime(void*, double*)
    int getcputime(void*)
    int get_sizeHistN(void*, int*)
    int get_version(void*, char*)
#B parameters

#B free memory
    void EndRun_FreeMemory_cmd(void*, void*)
    void EndRun_FreeMemory_gd(void*, void*)
    void EndRun_FreeMemory_gd_2(void*, void*)
    void EndRun_FreeMemory_histograms(void*, void*)
    void EndRun_FreeMemory_tree(void*, void*)
    void EndRun_FreeMemory_bodytable(void*, void*)
#E

#B histograms
    int get_rBins(void*, void*)
    int get_HistNN(void*, void*)
    int get_HistCF(void*, void*)
    int get_HistXi2pcf(void*, void*)
    int get_HistZetaMsincos(void*, void*, int, int, char*)
#E histograms


