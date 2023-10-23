/* ==============================================================================
!   HEADER: cmdline_defs.h		[tpcf]                                          !
!   Written by: Mario A. Rodriguez-Meza                                         !
!   Starting date: april 2023                                                   !
!   Purpose: Definitions for importing arguments from the command line          !
!   Language: C                                                                 !
!   Use: '#include "...."                                                       !
!   Major revisions:                                                            !
 ==============================================================================*/

#ifndef _cmdline_defs_h
#define _cmdline_defs_h

#define HEAD1	"LSST/CosmoININ"
#define HEAD2	"cballs Code for computing 3pcf, expected: O(N logN)"
#define HEAD3	"..."

string defv[] = {  ";"HEAD1": " HEAD2 "\n\t " HEAD3,
    "paramfile=",			            ";Parameter input file. Overwrite what follows",
//
#ifndef BALLSSINCOSEXEC
#ifdef BALLSEXEC
    "searchMethod=balls-omp",             ";Searching method to use", ":search",
#else
    "searchMethod=tree-omp",             ";Searching method to use", ":search",
#endif
#else
    "searchMethod=balls-omp",             ";Searching method to use", ":search",
#endif
#ifdef BALLS
    "theta=1.0",                        ";Control tree search parameter, can be used to increase speed (also use with options=BH86 or SW94 with search=balls4)",
    "scanLevel=6",                     ";Scan level to start the search (look at tdepth value, will be the maximum for this parameter)", ":scl",
#else
    "theta=1.0",                        ";Control tree search parameter, can be used to increase speed",
#endif
//
    "infile=",                          ";File name with points to analyse", ":in",
    "infileformat=columns-ascii",       ";Data input file format (columns-ascii or binary)", ":infmt",
//B Parameters to set a region in the sky, for example for Takahasi data set.
    "thetaL=1.279928",                      ";Angle theta left side of the region",
    "thetaR=1.861664",                      ";Angle theta right side of the region",
    "phiL=1.280107",                        ";Angle phi left side of the region",
    "phiR=1.861869",                        ";Angle phi right side of the region",
//E

//    "dimension=3",                      ";Number of dimensions of the run (3D or 2D systems); dummy by now, activate in Makefile_settings", ":dim",
    "mChebyshev=8",                     ";Number of Chebyshev polynomial to use (m+1)", ":mcheb",
//
    "rootDir=Output",                   ";Output dir, where output files will be written", ":root",
    "outfile=",                         ";Output file of the points analysed (default ext: .txt)", ":o",
    "outfileformat=columns-ascii",      ";Data output file format (columns-ascii or binary)", ":ofmt",
//B Parameters to control histograms
    "sizeHistN=40",                     ";array size for N histogram",
    "rangeN=100.0",                     ";Radial range for N histogram",
    "rminHist=1.0e-3",                     ";Radial minimum value for histograms", ":rmin",
//    "logHist=false",                    ";Log radial histograms; dummy by now, activate in Makefile_settings",
    "sizeHistTheta=40",                 ";array size for angular histogram (used when computing 3pcf brute force)",

    "histNFileName=histN",              ";File name (without extension) to save histograms of N", ":histNfname",
    "histXi2pcfFileName=histXi2pcf",    ";File name (without extension) to save histograms of Xi2pcf", ":histXi2pcffname",
    "histZetaMFileName=histZetaM",      ";File name (without extension) to save histograms of matrix ZetaM", ":histZMfname",
    "mhistZetaFileName=mhistZeta",      ";File name (without extension) to save histograms of matrix ZetaM in some directions", ":mhistZfname",
    "suffixOutFiles=",                  ";Suffix to add to output filenames", ":suffix",
//E
    "stepState=10000",                  ";number of steps to save a state-run info; or to save a state of the run which is not working yet",
//B Work in progress on this section
#ifdef SAVERESTORE
    "statefile=state",                  ";Write run state to a file; use in combination with stepState options above; not working yes", ":state",
    "restorefile=",                     ";Continue run from state file; not working yes", ":restore",
#endif
//E
    "verbose=1",                        ";Option to activate the amount of information sent to standard output", ":verb",
    "verbose_log=1",                    ";Option to activate the amount of information sent to log file", ":verblog",

    "numberThreads=8",                  ";To set the number of threads to use (OpenMP)", ":nthreads",

//B NOLSST:
    "seed=123",                         ";Random number seed to test run or useful to change a random region in Takahasi simulations",
    "nsmooth=2",                        ";Number of bodies to smooth out (or in a bucket)", ":nsm",
#ifdef BALLS
    "ntosave=1000",                     ";Number of found bodies to save; use in combination with 'bodyfound', balls4' method", ":ntsav",
#endif
    "stepNodes=1",                     ";number of nodes to jump and make a list of nodes to search in balls4 method",
    "ncritical=100,10000",               ";Range of number of bodies select nodes to search in balls4 method (use also theta=0.5)", ":nc",
//#endif
    "testmodel=simple-cubic-random",    ";Test model name to analyse", ":tstmodel",
    "nbody=16384",                      ";Number of points to test",
    "lengthBox=10000",                  ";Length of the box to test", ":lbox",
//    "mToPlot=1",                        ";m to plot histograms of histZetaM in some direcctions (only with search=direct-simple-exp)",
//E

    "options=",				            ";Various control options, i.e., no-one-ball (to use one-ball scheme),  compute-HistN, bh86, etc.", ":opt",
    "Version=0.1",			            ";Mario A. Rodríguez-Meza (2023)",
    NULL,
};

#endif // ! _cmdline_defs_h
