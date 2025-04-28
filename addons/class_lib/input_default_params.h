
#ifndef _default_params_unit_sphere_h_
#define _default_params_unit_sphere_h_

// input_default_params_unit_sphere -> input_default_params
//  name does not change
int input_default_params(struct cmdline_data *cmd)
{
// Every item in cmdline_defs.h must have an item here::

    //B Parameters related to the searching method
#ifdef BALLS
    cmd->searchMethod = "balls-omp";
#else
    cmd->searchMethod = "tree-omp-sincos";
#endif
    cmd->mChebyshev = 7;
    cmd->nsmooth = "1";
    cmd->rsmooth = "";
    cmd->theta = 1.0;
    cmd->computeTPCF = 1;
    //B correction 2025-04-06
    // Move this to addon that computes shear correlations
//    cmd->computeShearCF = 0;
    //E
    cmd->usePeriodic = 0;
    //E

    //B Parameters about the I/O file(s)
    // Input catalog parameters
    cmd->infile = "";
    cmd->infilefmt = "columns-ascii";
    cmd->iCatalogs = "1";
    // Output parameters
    cmd->rootDir = "Output";
    cmd->outfile = "";
    cmd->outfilefmt = "columns-ascii";
    // Parameters to set a region in the sky, for example for Takahasi data set
    cmd->thetaL = 1.279928;
    cmd->thetaR = 1.861664;
    cmd->phiL = 1.280107;
    cmd->phiR = 1.861869;
    //E

    //B Parameters to control histograms and their output files
    cmd->useLogHist = 1;
    cmd->logHistBinsPD = 5;
    //
    cmd->sizeHistN = 20;
    cmd->rangeN = 0.0633205;                        // RADTOARCMIN   3437.74677
    cmd->rminHist = 0.00213811;
    cmd->sizeHistPhi = 32;
    //
    cmd->histNNFileName = "histNN";
    cmd->histXi2pcfFileName = "histXi2pcf";
    cmd->histZetaFileName = "histZeta";
    cmd->suffixOutFiles = "";
    //E

    //B Set of parameters needed to construct a test model
    cmd->seed=123;                                          // to always have
                                                            //  defaults Check in gsl
    cmd->testmodel = "unit-sphere-random";
    cmd->nbody = 262144;
    cmd->lengthBox = 2.0;
    //E

    //B Miscellaneous parameters
    cmd->script = "";
    cmd->stepState = 100000;
    cmd->verbose = 2;
    cmd->verbose_log = 2;
#ifdef OPENMPCODE
    cmd->numthreads = 16;
#endif
    cmd->options = "out-m-HistZeta";
    //E

//B socket:
#ifdef ADDONS
#include "class_lib_include_02.h"
#endif
//

  return SUCCESS;
}

#endif // ! _default_params_unit_sphere_h_
