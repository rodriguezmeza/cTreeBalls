
#ifndef _cmdline_defs_02_test_unit_sphere_kappa_h
#define _cmdline_defs_02_test_unit_sphere_kappa_h

//
// set of default values of the parameters useful for testing
//  correlation of the convergence field on the unit sphere
//

string defv[] = {  ";"HEAD1": " HEAD2 "\n\t " HEAD3,
    "paramfile=",			            ";Parameter input file. Overwrite what follows",

    //B Parameters related to the searching method

//B socket:
#ifdef ADDONS
#include "cmdline_defs_include_00.h"
#else
    "searchMethod=tree-omp-sincos",     ";Searching method to use", ":search",
#endif
//E

    "mChebyshev=7",                     ";Number of Chebyshev polynomial to use (m+1)", ":mcheb",
    "nsmooth=1",                        ";Number of bodies to smooth out (or in a bucket)", ":nsm",
    "rsmooth=",                         ";Radius of the pivot smoothing neighbourhood. If empty a default is set", ":rsm",
    "theta=1.0",                        ";Control tree search parameter, can be used to increase speed",
    "computeTPCF=true",                 ";If true, compute 3pcf", ":tpcf",
    //B correction 2025-04-06
    // Move this to addon that computes shear correlations
//    "computeShearCF=true",              ";If true, compute shear cf", ":shearcf",
    //E
    "usePeriodic=false",                ";If false, don't use periodic boundary condition", ":periodic",
    //E

    //B Parameters to control the I/O file(s)
    // Input catalog parameters
    "infile=",                          ";File names with points to analyse", ":in",
    "infileformat=columns-ascii",       ";Data input files format (columns-ascii or binary)", ":infmt",
    "iCatalogs=1",                      ";index of point catalogs to analyse", ":icats",
    // Output parameters
    "rootDir=Output",                   ";Output dir, where output files will be written", ":root",
    "outfile=",                         ";Output file of the points analysed (default ext: .txt)", ":o",
    "outfileformat=columns-ascii",      ";Data output file format (columns-ascii or binary)", ":ofmt",
    // Parameters to set a region in the sky, for example for Takahasi data set
    "thetaL=1.279928",                  ";Angle theta left side of the region",
    "thetaR=1.861664",                  ";Angle theta right side of the region",
    "phiL=1.280107",                    ";Angle phi left side of the region",
    "phiR=1.861869",                    ";Angle phi right side of the region",
    //E

    //B Parameters to control histograms and their output files
    "useLogHist=true",                  ";If true, use logaritmic scale for histograms", ":loghist",
    "logHistBinsPD=5",                  ";Bins per decades", ":binspd",
    //
    "sizeHistN=20",                     ";array size for N histogram",
    // RADTOARCMIN   3437.74677
    "rangeN=0.0633205",                 ";Radial range for N histogram (217.68 arcmin)",
    "rminHist=0.00213811",              ";Radial minimum value for histograms (7.35 arcmin)", ":rmin",
    "sizeHistPhi=32",                   ";array size for angular histogram",
    //
    "histNNFileName=histNN",            ";File name (without extension) to save histograms of NN", ":histNNfname",
    "histXi2pcfFileName=histXi2pcf",    ";File name (without extension) to save histograms of Xi2pcf", ":histXi2pcffname",
    "histZetaFileName=histZeta",        ";Prefix of file name to save histograms of matrix Zeta", ":histZfname",
    "suffixOutFiles=",                  ";Suffix to add to output filenames", ":suffix",
    //E

    //B Set of parameters needed to construct a test model
    "seed=123",                         ";Random number seed to test run or useful to change a random region in Takahasi simulations",
    "testmodel=unit-sphere-random",     ";Test model name to analyse: simple-cubic-random, unit-sphere-random,...", ":tstmodel",
    "nbody=262144",                     ";Number of points to test (2^18)",
    "lengthBox=2",                      ";Length of the box to test", ":lbox",
    //E

    //B Miscellaneous parameters
//    "script=",                          ";Scripts in shell or python that can be run in pre-processing or post-processing",
    "preScript=",                          ";Script in shell or python that can be run in pre-processing",
    "posScript=",                          ";Script in shell or python that can be run in post-processing",
    "stepState=100000",                 ";number of steps to save a state-run info (pivot number completed in the log file)",
    "verbose=0",                        ";Option to activate the amount of information sent to standard output", ":verb",
    "verbose_log=0",                    ";Option to activate the amount of information sent to log file", ":verblog",
#ifdef OPENMPCODE
    "numberThreads=16",                 ";To set the number of threads to use (OpenMP)", ":nthreads",
#endif
    "options=out-m-HistZeta",           ";Various control options, i.e., no-one-ball (to use one-ball scheme),  compute-HistN, bh86, etc.", ":opt",
    //E

//B socket:
#ifdef ADDONS
#include "cmdline_defs_include.h"
#endif
//E

    "Version=1.0.0",			        ";Mario A. Rodríguez-Meza (2023-2025)",
    NULL,
};

#endif // ! _cmdline_defs_02_test_unit_sphere_kappa_h
