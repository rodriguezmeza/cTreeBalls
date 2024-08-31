/*==============================================================================
 HEADER: cmdline_data.h		[cTreeBalls]
 Written by: Mario A. Rodriguez-Meza
 Starting date: april 2023
 Purpose: Definitions of global variables and parameters
 Language: C
 Use: '#include "cmdline_data.h"
 Major revisions:
 ==============================================================================*/
//        1          2          3          4        ^ 5          6          7

#ifndef _cmdline_data_h
#define _cmdline_data_h

#include "common_defs.h"


struct cmdline_data{
// Every item in cmdline_defs.h must have an item here::

    string paramfile;

#ifndef GETPARAM
    char ParameterFile[MAXLENGTHOFFILES];           // May be we should incrase
                                                    //  this number
#endif

    //B Parameters related to the searching method
    string searchMethod;
    int mChebyshev;
    string nsmooth;
    string rsmooth;
    real theta;
    bool computeTPCF;
    bool computeShearCF;
    bool usePeriodic;
    //E

    //B Parameters about the I/O file(s)
    // Input catalog parameters
    string infile;
    string infilefmt;
    string iCatalogs;
    // Output parameters
    string rootDir;
    string outfile;
    string outfilefmt;
    // Parameters to set a region in the sky, for example for Takahasi data set
    real thetaL;
    real thetaR;
    real phiL;
    real phiR;
    //E

    //B Parameters to control histograms and their output files
    bool useLogHist;
    int logHistBinsPD;
    //
    int sizeHistN;
    real rangeN;
    real rminHist;
    int sizeHistTheta;
    //
    string histNNFileName;
    string histXi2pcfFileName;
    string histZetaFileName;
    string suffixOutFiles;
    //E

    //B Set of parameters needed to construct a test model
    int seed;
    string testmodel;
    INTEGER nbody;
    real lengthBox;
    //E

    //B Miscellaneous parameters
    string script;
    INTEGER stepState;
    short verbose;
    short verbose_log;
#ifdef OPENMPCODE
    int numthreads;
#endif
    string options;
    //E

    string version;

#ifdef ADDONS
#include "globaldefs_include_01.h"
#endif
};

#endif // ! _cmdline_data_h

