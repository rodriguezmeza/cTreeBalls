
#ifndef _default_params_unit_sphere_h_
#define _default_params_unit_sphere_h_

// input_default_params_unit_sphere -> input_default_params
//  name does not change
int input_default_params(struct cmdline_data *cmd)
{
// Every item in cmdline_defs.h must have an item here::

    //B Parameters related to the searching method
    cmd->searchMethod = "tree-omp-sincos";
    cmd->mChebyshev = 7;
    cmd->nsmooth = 8;
    //E
    cmd->rsmooth = "\0";
    cmd->theta = 1.05;
    cmd->usePeriodic = 0;
    //E

    //B Parameters about the I/O file(s)
    // Input catalog parameters
    cmd->infile = "\0";
    //E
    cmd->infilefmt = "columns-ascii";
    cmd->iCatalogs = "1";
    // Output parameters
    cmd->rootDir = "Output";
    cmd->outfile = "\0";
    cmd->outfilefmt = "columns-ascii";
    // Parameters to set a region in the sky, for example for Takahashi data
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
    cmd->seed=123;                                  // to always have
                                                    //  defaults Check in gsl
    cmd->testmodel = "unit-sphere-random";
    cmd->nbody = 262144;
    cmd->lengthBox = 2.0;
    //E

    //B Miscellaneous parameters
    cmd->preScript = "";
    cmd->posScript = "";
    cmd->stepState = 100000;
    cmd->verbose = 0;
    cmd->verbose_log = 0;
#ifdef OPENMPCODE
    cmd->numthreads = 16;
#endif
    cmd->options = "\0";
    //E

//B socket:
#ifdef ADDONS
#include "class_lib_include_02.h"
#endif
//

  return SUCCESS;
}

#endif // ! _default_params_unit_sphere_h_
