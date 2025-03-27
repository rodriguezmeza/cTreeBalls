// Use:
//NOLSST:
//#include "startrun_save_restore_00.h"

#ifndef _startrun_save_restore_00_h
#define _startrun_save_restore_00_h



#ifdef GETPARAM

local void startrun_ParamStat(struct  cmdline_data*,
                              struct  global_data*);

//B Section of parameter stat
local void startrun_ParamStat(struct  cmdline_data* cmd,
                              struct  global_data* gd)
{
// Every item in cmdline_defs.h must have an item here::

    //B Parameters related to the searching method
    if (GetParamStat("searchMethod") & ARGPARAM)
        cmd->searchMethod = GetParam("searchMethod");
    if (GetParamStat("mChebyshev") & ARGPARAM)
        cmd->mChebyshev = GetiParam("mChebyshev");
    if (GetParamStat("nsmooth") & ARGPARAM)
        cmd->nsmooth = GetParam("nsmooth");
    if (GetParamStat("rsmooth") & ARGPARAM)
        cmd->rsmooth = GetParam("rsmooth");
    if (GetParamStat("theta") & ARGPARAM)
        cmd->theta = GetdParam("theta");
    if (GetParamStat("computeTPCF") & ARGPARAM)
        cmd->computeTPCF = GetbParam("computeTPCF");
    if (GetParamStat("computeShearCF") & ARGPARAM)
        cmd->computeShearCF = GetbParam("computeShearCF");
    if (GetParamStat("usePeriodic") & ARGPARAM)
        cmd->usePeriodic = GetbParam("usePeriodic");
    //E

    //B Parameters to control the I/O file(s)
    // Input catalog parameters
    if (GetParamStat("infile") & ARGPARAM)
        cmd->infile = GetParam("infile");
    if (GetParamStat("infilefmt") & ARGPARAM)
        cmd->infilefmt = GetParam("infileformat");
    if (GetParamStat("iCatalogs") & ARGPARAM)
        cmd->iCatalogs = GetParam("iCatalogs");
    // Output parameters
    if (GetParamStat("rootDir") & ARGPARAM)
        cmd->rootDir = GetParam("rootDir");
    if (GetParamStat("outfile") & ARGPARAM)
        cmd->outfile = GetParam("outfile");
    if (GetParamStat("outfilefmt") & ARGPARAM)
        cmd->outfilefmt = GetParam("outfileformat");
    // Parameters to set a region in the sky, for example for Takahasi data set
    if (GetParamStat("thetaL") & ARGPARAM)
        cmd->thetaL = GetdParam("thetaL");
    if (GetParamStat("thetaR") & ARGPARAM)
        cmd->thetaR = GetdParam("thetaR");
    if (GetParamStat("phiL") & ARGPARAM)
        cmd->phiL = GetdParam("phiL");
    if (GetParamStat("thetaR") & ARGPARAM)
        cmd->phiR = GetdParam("phiR");
    //E

    //B Parameters to control histograms and their output files
    if (GetParamStat("useLogHist") & ARGPARAM)
        cmd->useLogHist = GetbParam("useLogHist");
    if (GetParamStat("logHistBinsPD") & ARGPARAM)
        cmd->logHistBinsPD = GetiParam("logHistBinsPD");
    //
    if (GetParamStat("sizeHistN") & ARGPARAM)
        cmd->sizeHistN = GetiParam("sizeHistN");
    if (GetParamStat("rangeN") & ARGPARAM)
        cmd->rangeN = GetdParam("rangeN");
    if (GetParamStat("rminHist") & ARGPARAM)
        cmd->rminHist = GetdParam("rminHist");
    if (GetParamStat("sizeHistPhi") & ARGPARAM)
        cmd->sizeHistPhi = GetiParam("sizeHistPhi");
    //
    if (GetParamStat("histNNFileName") & ARGPARAM)
        cmd->histNNFileName = GetParam("histNNFileName");
    if (GetParamStat("histXi2pcfFileName") & ARGPARAM)
        cmd->histXi2pcfFileName = GetParam("histXi2pcfFileName");
    if (GetParamStat("histZetaFileName") & ARGPARAM)
        cmd->histZetaFileName = GetParam("histZetaFileName");
    if (GetParamStat("suffixOutFiles") & ARGPARAM)
        cmd->suffixOutFiles = GetParam("suffixOutFiles");
    //E
    
    //B Set of parameters needed to construct a test model
    if (GetParamStat("seed") & ARGPARAM)
        cmd->seed = GetiParam("seed");
    if (GetParamStat("testmodel") & ARGPARAM)
        cmd->testmodel = GetParam("testmodel");
    if (GetParamStat("nbody") & ARGPARAM)
#ifdef LONGINT
        cmd->nbody = GetlParam("nbody");
#else
        cmd->nbody = GetiParam("nbody");
#endif
    if (GetParamStat("lengthBox") & ARGPARAM)
        cmd->lengthBox = GetdParam("lengthBox");
    //E

    //B Miscellaneous parameters
    if (GetParamStat("script") & ARGPARAM)
        cmd->script = GetParam("script");
    if (GetParamStat("verbose") & ARGPARAM)
        cmd->verbose = GetiParam("verbose");
    if (GetParamStat("verbose_log") & ARGPARAM)
        cmd->verbose_log = GetiParam("verbose_log");
#ifdef OPENMPCODE
    if (GetParamStat("numberThreads") & ARGPARAM)
        cmd->numthreads = GetiParam("numberThreads");
#endif
    if (GetParamStat("options") & ARGPARAM)
        cmd->options = GetParam("options");
    //E

#ifdef ADDONS
#include "startrun_include_06.h"
#endif
}
//E
#endif


#endif	// ! _startrun_save_restore_00_h
