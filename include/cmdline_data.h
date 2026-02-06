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

//
// lines where there is a "//B socket:" string are places to include module files
//  that can be found in addons/addons_include folder
//

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
    int nsmooth;
    //E
    string rsmooth;                                 // use to smooth pivot
    real theta;
    bool usePeriodic;
    //E

    //B Parameters about the I/O file(s)
    // Input catalog parameters
    string infile;
    string infilefmt;
    //E
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
    int sizeHistPhi;
    //
    string histNNFileName;
    string histXi2pcfFileName;
    string histZetaFileName;
    string suffixOutFiles;                  // deprecate... in favor o no .txt
    //E

    //B Set of parameters needed to construct a test model
    int seed;
    string testmodel;
    INTEGER nbody;
    real lengthBox;
    //E

    //B Miscellaneous parameters
    string preScript;
    string posScript;
    INTEGER stepState;
    short verbose;
    short verbose_log;
#ifdef OPENMPCODE
    int numthreads;
#endif
    string options;
    //E

    string version;

//B socket:
#ifdef ADDONS
#include "globaldefs_include_01.h"
#endif
//E
};

#endif // ! _cmdline_data_h
