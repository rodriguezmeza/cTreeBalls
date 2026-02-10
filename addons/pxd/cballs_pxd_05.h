// Use:
//#include "cballs_pxd_05.h"

// included in addons/addons_include/source/cballs_include_05.h

#ifndef _cballs_pxd_05_h
#define _cballs_pxd_05_h

//B parameters section
/*
 paramfile=                      Parameter input file. Overwrite what follows
 searchMethod=balls-omp          Searching method to use [a: search]
 mChebyshev=7                    Number of Chebyshev polynomial to use (m+1) [a: mcheb]
 nsmooth=1                       Number of bodies to smooth out (or in a bucket) [a: nsm]
 rsmooth=                        Radius of the pivot smoothing neighbourhood. If empty a default is set [a: rsm]
 theta=1.0                       Control tree search parameter, can be used to increase speed
 computeTPCF=true                If true, compute 3pcf [a: tpcf]
 usePeriodic=false               If false, don't use periodic boundary condition [a: periodic]
 infile=                         File names with points to analyse [a: in]
 infileformat=columns-ascii      Data input files format (columns-ascii or binary) [a: infmt]
 iCatalogs=1                     index of point catalogs to analyse [a: icats]
 rootDir=Output                  Output dir, where output files will be written [a: root]
 outfile=                        Output file of the points analysed (default ext: .txt) [a: o]
 outfileformat=columns-ascii     Data output file format (columns-ascii or binary) [a: ofmt]
 thetaL=1.279928                 Angle theta left side of the region
 thetaR=1.861664                 Angle theta right side of the region
 phiL=1.280107                   Angle phi left side of the region
 phiR=1.861869                   Angle phi right side of the region
 useLogHist=true                 If true, use logaritmic scale for histograms [a: loghist]
 logHistBinsPD=5                 Bins per decades [a: binspd]
 sizeHistN=20                    array size for N histogram
 rangeN=0.0633205                Radial range for N histogram (217.68 arcmin)
 rminHist=0.00213811             Radial minimum value for histograms (7.35 arcmin) [a: rmin]
 sizeHistPhi=32                  array size for angular histogram
 histNNFileName=histNN           File name (without extension) to save histograms of NN [a: histNNfname]
 histXi2pcfFileName=histXi2pcf   File name (without extension) to save histograms of Xi2pcf [a: histXi2pcffname]
 histZetaFileName=histZeta       Prefix of file name to save histograms of matrix Zeta [a: histZfname]
 suffixOutFiles=                 Suffix to add to output filenames [a: suffix]
 seed=123                        Random number seed to test run or useful to change a random region in Takahasi simulations
 testmodel=unit-sphere-random    Test model name to analyse: simple-cubic-random, unit-sphere-random,... [a: tstmodel]
 nbody=262144                    Number of points to test (2^18)
 lengthBox=2                     Length of the box to test [a: lbox]
 preScript=                      Script in shell or python that can be run in pre-processing
 posScript=                      Script in shell or python that can be run in post-processing
 stepState=100000                number of steps to save a state-run info (pivot number completed in the log file)
 verbose=0                       Option to activate the amount of information sent to standard output [a: verb]
 verbose_log=0                   Option to activate the amount of information sent to log file [a: verblog]
 numberThreads=16                To set the number of threads to use (OpenMP) [a: nthreads]
 options=out-m-HistZeta          Various control options, i.e., no-one-ball (to use one-ball scheme),  compute-HistN, bh86, etc. [a: opt]
 scanLevel=6                     Scan level to start the search (look at tdepth value, will be the maximum for this parameter) [a: scl]
 scanLevelRoot=0                 Scan level of root cells to start the search (look at tdepth value, will be the maximum for this parameter) [a: sclroot]
 scanLevelMin=-0                 Scan level of size cells to stop the search. Integer negative values (look at tdepth value, will be tdepth-1+scanLevelMin+1) [a: sclmin]
// does not used... it was remove from all instances...
 ntosave=1000                    Number of found bodies to save; use in combination with 'bodyfound', balls4' method [a: ntsav]
//
 columns=1,2,3,4                 Columns to use as vector position and fields when using multi-columns-ascii or fits formats of point catalog [a: cols]
 Version=1.0.1                   Mario A. Rodríguez-Meza (2023-2025)

 */

//B flags
int get_tree_allocated(struct global_data* gd, short *value)
{
    *value = gd->tree_allocated;
    return SUCCESS;
}

int get_allocated_2(struct global_data* gd, short *value)
{
    *value = gd->gd_allocated_2;
    return SUCCESS;
}

int get_bodytable_allocated(struct global_data* gd, short *value)
{
    *value = gd->bodytable_allocated;
    return SUCCESS;
}

int get_histograms_allocated(struct global_data* gd, short *value)
{
    *value = gd->histograms_allocated;
    return SUCCESS;
}

int get_gd_allocated(struct global_data* gd, short *value)
{
    *value = gd->gd_allocated;
    return SUCCESS;
}

int get_cmd_allocated(struct global_data* gd, short *value)
{
    *value = gd->cmd_allocated;
    return SUCCESS;
}
//E

int get_nthreads(struct  cmdline_data* cmd, int *value)
{
    *value = cmd->numthreads;
    return SUCCESS;
}

//B version 1.0.1
int get_nmonopoles(struct  cmdline_data* cmd, int *value)
{
    *value = cmd->mChebyshev;
    return SUCCESS;
}
//E

int get_theta(struct  cmdline_data* cmd, real *theta)
{
    *theta = cmd->theta;
    return SUCCESS;
}

int get_rsmooth(struct  global_data* gd, real *value)
{
    *value = gd->rsmooth[0];
    return SUCCESS;
}

int get_cputime(struct  global_data* gd, real *cputime)
{
    *cputime = gd->cpusearch;
    return SUCCESS;
}

int get_sizeHistN(struct  cmdline_data* cmd, int *sizeHistN)
{
    *sizeHistN = cmd->sizeHistN;
    return SUCCESS;
}

// use same version as is in cmdline_defs.h and/or in addons/class_lib/common.h
//  see also setup.py
int get_version(struct  cmdline_data* cmd, char *param)
{
    sprintf(param,"%s","1.0.1");
    return SUCCESS;
}

int get_rootDir(struct  cmdline_data* cmd, char *value)
{
    sprintf(value,"%s",cmd->rootDir);
    return SUCCESS;
}
//E parameters section


//B histograms section
int get_rBins(struct  cmdline_data* cmd, struct  global_data* gd)
{
    real rBin, rbinlog;
    int n;

    for (n=1; n<=cmd->sizeHistN; n++) {
        if (cmd->useLogHist) {
            if (cmd->rminHist==0) {
                rbinlog = ((real)(n-cmd->sizeHistN))/cmd->logHistBinsPD + rlog10(cmd->rangeN);
            } else {
                rbinlog = rlog10(cmd->rminHist) + ((real)(n)-0.5)*gd->deltaR;
            }
            rBin=rpow(10.0,rbinlog);
        } else {
            rBin = cmd->rminHist + ((real)n-0.5)*gd->deltaR;
        }
        gd->rBins[n] = rBin;
    }

    return SUCCESS;
}

int get_HistNN(struct  cmdline_data* cmd, struct  global_data* gd)
{
    int n;

    for (n=1; n<=cmd->sizeHistN; n++) {
        gd->vecPXD[n] = gd->histNN[n];
    }

    return SUCCESS;
}

int get_HistCF(struct  cmdline_data* cmd, struct  global_data* gd)
{
    int n;

    for (n=1; n<=cmd->sizeHistN; n++) {
        gd->vecPXD[n] = gd->histCF[n];
    }

    return SUCCESS;
}

int get_HistXi2pcf(struct  cmdline_data* cmd, struct  global_data* gd)
{
    int n;

    for (n=1; n<=cmd->sizeHistN; n++) {
        gd->vecPXD[n] = gd->histXi2pcf[n];
    }

    return SUCCESS;
}


// get matrix ZetaM for each m multipole
#define _COS_         1
#define _SIN_         2
#define _SINCOS_      3
#define _COSSIN_      4

int get_HistZetaMsincos(struct  cmdline_data* cmd,
                     struct  global_data* gd,
                     int m, int type, ErrorMsg errmsg)
{
    int n1, n2;

    class_test((m <= 0 || m > cmd->mChebyshev + 1),
        errmsg,"\nget_HistZetaM_sincos: not allowed value of m = %d\n",m);
    for (n1=1; n1<=cmd->sizeHistN; n1++) {
        for (n2=1; n2<=cmd->sizeHistN; n2++) {
            gd->matPXD[n1][n2] = gd->histZetaMcos[m][n1][n2];
        }
    }

    switch(type) {
        case _COS_:
            for (n1=1; n1<=cmd->sizeHistN; n1++) {
                for (n2=1; n2<=cmd->sizeHistN; n2++) {
                    gd->matPXD[n1][n2] = gd->histZetaMcos[m][n1][n2];
                }
            }
            break;
        case _SIN_:
            for (n1=1; n1<=cmd->sizeHistN; n1++) {
                for (n2=1; n2<=cmd->sizeHistN; n2++) {
                    gd->matPXD[n1][n2] = gd->histZetaMsin[m][n1][n2];
                }
            }
            break;
        case _SINCOS_:
            for (n1=1; n1<=cmd->sizeHistN; n1++) {
                for (n2=1; n2<=cmd->sizeHistN; n2++) {
                    gd->matPXD[n1][n2] = gd->histZetaMsincos[m][n1][n2];
                }
            }
            break;
        case _COSSIN_:
            for (n1=1; n1<=cmd->sizeHistN; n1++) {
                for (n2=1; n2<=cmd->sizeHistN; n2++) {
                    gd->matPXD[n1][n2] = gd->histZetaMcossin[m][n1][n2];
                }
            }
            break;
        default:
            for (n1=1; n1<=cmd->sizeHistN; n1++) {
                for (n2=1; n2<=cmd->sizeHistN; n2++) {
                    gd->matPXD[n1][n2] = gd->histZetaMcos[m][n1][n2];
                }
            }
            break;
    }

    return SUCCESS;
}

#undef _COS_
#undef _SIN_
#undef _SINCOS_
#undef _COSSIN_
//E histograms section

//E PXD functions


#endif	// ! _cballs_pxd_05_h
