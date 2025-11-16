# Set of declarations from C to python
# Only define quantities you will use,
# for input, output or intermediate manipulations
# by the python wrapper. The rest is internal in cBalls

cdef extern from "cballs.h":

    ctypedef char FileArg[40]
    ctypedef char* ErrorMsg

    cdef struct cmdline_data:
        double theta
        int sizeHistN
        double rangeN
        double rminHist
        char * preScript
        char * posScript

    cdef struct global_data:
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
    int PrintParameterFile(void*, char*)

    int StartOutput(void*)
    int SetNumberThreads(void*)

    int MainLoop(void*, void*);
    int EndRun(void*, void*);

#B parameters
    int get_theta(void*, double*)
    int get_sizeHistN(void*, int*)
    int get_version(void*, char*)
#B parameters
#B histograms
    int get_rBins(void*, void*)
    int get_HistNN(void*, void*)
    int get_HistCF(void*, void*)
    int get_HistXi2pcf(void*, void*)
    int get_HistZetaMsincos(void*, void*, int, int, char*)
    int get_HistZetaM_sincos(void*, void*, int, int, char*)
#E histograms

