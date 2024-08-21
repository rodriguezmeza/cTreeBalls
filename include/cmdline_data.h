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

    INTEGER stepState;

    string script;
    string options;
    string version;
    short verbose;
    short verbose_log;
	string paramfile;

#ifndef GETPARAM
    char ParameterFile[MAXLENGTHOFFILES];           // May be we should incrase
                                                    //  this number
#endif

    string searchMethod;

    real theta;

//B TPCF
    int sizeHistTheta;
    int mChebyshev;
//E

//B Parameters to set a region in the sky, for example for Takahasi data set.
    real thetaL;
    real thetaR;
    real phiL;
    real phiR;
//E

    int sizeHistN;
    real rangeN;
    real rminHist;

    string infile;
    string infilefmt;

    string iCatalogs;

    string rootDir;
    string outfile;
    string outfilefmt;
    

    string histNNFileName;
    string histXi2pcfFileName;
    string histZetaFileName;
    string suffixOutFiles;
    
#ifdef OPENMPCODE
    int numthreads;
#endif

    int seed;
    string nsmooth;

    string rsmooth;

#ifdef ADDONS
#include "globaldefs_include_01.h"
#endif

    string testmodel;
    INTEGER nbody;
    real lengthBox;

    bool computeTPCF;
    bool useLogHist;
    int logHistBinsPD;
    bool usePeriodic;

};

#endif // ! _cmdline_data_h

