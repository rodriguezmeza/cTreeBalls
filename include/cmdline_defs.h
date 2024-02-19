/* ==============================================================================
!   HEADER: cmdline_defs.h		[cTreeBalls]                                    !
!   Written by: Mario A. Rodriguez-Meza                                         !
!   Starting date: april 2023                                                   !
!   Purpose: Definitions for importing arguments from the command line          !
!   Language: C                                                                 !
!   Use: '#include "...."                                                       !
!   Major revisions:                                                            !
 ==============================================================================*/
//        1          2          3          4          5          6          7

/*
If you need to add a parameter, the pattern is a line like the following:
 
 "parameter_name=default_value",   ";Parameter comment/or brief explanation", ":parameter_shortname",

Instert a line like above inside the defv array.
 
Note: it is not necessary to have the parameter_shortname item. Then in this case the line for a parameter is:
 
 "parameter_name=default_value",   ";Parameter comment/or brief explanation",

 */


#ifndef _cmdline_defs_h
#define _cmdline_defs_h

#define HEAD1	"LSST/CosmoININ"
#define HEAD2	"cballs Code for computing 3pcf, expected: O(N logN)"
#define HEAD3	"..."

string defv[] = {  ";"HEAD1": " HEAD2 "\n\t " HEAD3,
    "paramfile=",			            ";Parameter input file. Overwrite what follows",
//
    "searchMethod=tree-omp-sincos",     ";Searching method to use", ":search",
//
    "infile=",                          ";File names with points to analyse", ":in",
    "infileformat=columns-ascii",       ";Data input files format (columns-ascii or binary)", ":infmt",
    "iCatalogs=1",                      ";index of point catalogs to analyse", ":icats",
//E
//B Parameters to set a region in the sky, for example for Takahasi data set.
    "thetaL=1.279928",                  ";Angle theta left side of the region",
    "thetaR=1.861664",                  ";Angle theta right side of the region",
    "phiL=1.280107",                    ";Angle phi left side of the region",
    "phiR=1.861869",                    ";Angle phi right side of the region",
//E
    "mChebyshev=8",                     ";Number of Chebyshev polynomial to use (m+1)", ":mcheb",
//
    "rootDir=Output",                   ";Output dir, where output files will be written", ":root",
    "outfile=",                         ";Output file of the points analysed (default ext: .txt)", ":o",
    "outfileformat=columns-ascii",      ";Data output file format (columns-ascii or binary)", ":ofmt",
//B Parameters to control histograms
    "sizeHistN=40",                     ";array size for N histogram",
    "rangeN=100.0",                     ";Radial range for N histogram",
    "rminHist=1.0e-3",                  ";Radial minimum value for histograms", ":rmin",

    "histNFileName=histN",              ";File name (without extension) to save histograms of N", ":histNfname",
    "histXi2pcfFileName=histXi2pcf",    ";File name (without extension) to save histograms of Xi2pcf", ":histXi2pcffname",
    "histZetaMFileName=histZetaM",      ";File name (without extension) to save histograms of matrix ZetaM", ":histZMfname",
    "mhistZetaFileName=mhistZeta",      ";File name (without extension) to save histograms of matrix ZetaM in some directions", ":mhistZfname",
    "suffixOutFiles=",                  ";Suffix to add to output filenames", ":suffix",
//E
    "stepState=10000",                  ";number of steps to save a state-run info (pivot number completed in the log file)",

    "verbose=1",                        ";Option to activate the amount of information sent to standard output", ":verb",
    "verbose_log=1",                    ";Option to activate the amount of information sent to log file", ":verblog",

    "numberThreads=4",                  ";To set the number of threads to use (OpenMP)", ":nthreads",

    "seed=123",                         ";Random number seed to test run or useful to change a random region in Takahasi simulations",
    "nsmooth=2",                        ";Number of bodies to smooth out (or in a bucket)", ":nsm",
    "stepNodes=1",                      ";number of nodes to jump and make a list of nodes to search in balls4 method",
    "ncritical=100,10000",               ";Range of number of bodies select nodes to search in balls4 method (use also theta=0.5)", ":nc",

    "testmodel=simple-cubic-random",    ";Test model name to analyse", ":tstmodel",
    "nbody=16384",                      ";Number of points to test",
    "lengthBox=10000",                  ";Length of the box to test", ":lbox",

    "script=",                          ";Scripts in shell or python that can be run in pre-processing or post-processing.",

    "theta=1.0",                        ";Control tree search parameter, can be used to increase speed",

    "rsmooth=",                         ";Radius of the pivot smoothing neighbourhood", ":rsm",

#ifdef ADDONS
#include "cmdline_defs_include.h"
#endif

    "options=",                         ";Various control options, i.e., no-one-ball (to use one-ball scheme),  compute-HistN, bh86, etc.", ":opt",
    "Version=0.1",			            ";Mario A. Rodríguez-Meza (2023)",
    NULL,
};

#endif // ! _cmdline_defs_h
