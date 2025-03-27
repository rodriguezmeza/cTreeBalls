
#ifndef _cballsio_save_restore_11_h
#define _cballsio_save_restore_11_h

#define stopfilename    "stop"
#define stopsavestate    "stop-state"

/*
 Stop the run routine:
    It is activated with "echo > stop" in the terminal

 Global tructures used: gd, cmd
 Return:
    void
 */
global void checkstop(struct cmdline_data* cmd,
                      struct  global_data* gd)
{
    char   buf[200];
    FILE *fd;
    double cpudt;

    if ((fd=fopen(stopfilename,"r"))) {
        fclose(fd);
        gd->stopflag = 1;
        sprintf(buf,"rm -f %s", stopfilename);
        system(buf);
        printf("\nsaving a stop run state...\n\n");
        cpudt = CPUTIME-gd->cpuinit;
        gd->cputotal += cpudt;
        savestate(cmd, gd, stopsavestate);
        gd->cputotal -= cpudt;
    }
}

#undef stopfilename
#undef stopsavestate


#define savestatetmp    "savestate-tmp"

//B Save and restore state of the run

/*
 Save state of the run routine:

 Arguments:
    * `pattern`: name pattern of the file to save the state of the run
 Global tructures used: gd, cmd
 Counting encounters (in global gd): nbbcalc, nbccalc, ncccalc
 Return (the error status):
    int SUCCESS or FAILURE
 */
global void savestate(struct cmdline_data* cmd,
                      struct  global_data* gd,
                      string pattern)
{
// every item you add to cmdline_data, global_data structures
//  or to globaldefs.h file should you consider necessary to
//  recover the run, must be added to this routine.
// Also all items must be read by restorestate routine
//  in the same order as given here

    int m,n,k;
    double rval;
    char namebuf[256];
    stream str;
    int nchars;
    int ndim;
    int ifile;
    int i;

    double cpustart = CPUTIME;
    //B header: program name and version
    sprintf(namebuf, pattern);
    str = stropen(namebuf, "w!");
    nchars = strlen(getargv0()) + 1;
    safewrite(&nchars, sizeof(int), str);
    safewrite(getargv0(), nchars * sizeof(char), str);
    nchars = strlen(getversion()) + 1;
    safewrite(&nchars, sizeof(int), str);
    safewrite(getversion(), nchars * sizeof(char), str);
    //E

//=============================================
//B cmdline_data

    //B Parameters related to the searching method
    nchars = strlen(cmd->searchMethod) + 1;
    safewrite(&nchars, sizeof(int), str);
    safewrite(cmd->searchMethod, nchars * sizeof(char), str);

    safewrite(&cmd->mChebyshev, sizeof(int), str);

    nchars = strlen(cmd->nsmooth) + 1;
    safewrite(&nchars, sizeof(int), str);
    safewrite(cmd->nsmooth, nchars * sizeof(char), str);
    
    nchars = strlen(cmd->rsmooth) + 1;
    safewrite(&nchars, sizeof(int), str);
    safewrite(cmd->rsmooth, nchars * sizeof(char), str);

    safewrite(&cmd->theta, sizeof(real), str);
    safewrite(&cmd->computeTPCF, sizeof(bool), str);
    safewrite(&cmd->computeShearCF, sizeof(bool), str);
    safewrite(&cmd->usePeriodic, sizeof(bool), str);
    //E
    
    //B Parameters about the I/O file(s)
    // Input catalog parameters
    nchars = strlen(cmd->infile) + 1;
    safewrite(&nchars, sizeof(int), str);
    safewrite(cmd->infile, nchars * sizeof(char), str);

    nchars = strlen(cmd->infilefmt) + 1;
    safewrite(&nchars, sizeof(int), str);
    safewrite(cmd->infilefmt, nchars * sizeof(char), str);

    nchars = strlen(cmd->iCatalogs) + 1;
    safewrite(&nchars, sizeof(int), str);
    safewrite(cmd->iCatalogs, nchars * sizeof(char), str);

    // Output parameters
    nchars = strlen(cmd->rootDir) + 1;
    safewrite(&nchars, sizeof(int), str);
    safewrite(cmd->rootDir, nchars * sizeof(char), str);

    nchars = strlen(cmd->outfile) + 1;
    safewrite(&nchars, sizeof(int), str);
    safewrite(cmd->outfile, nchars * sizeof(char), str);

    nchars = strlen(cmd->outfilefmt) + 1;
    safewrite(&nchars, sizeof(int), str);
    safewrite(cmd->outfilefmt, nchars * sizeof(char), str);

    // Parameters to set a region in the sky, for example for Takahasi data set
    safewrite(&cmd->thetaL, sizeof(real), str);
    safewrite(&cmd->thetaR, sizeof(real), str);
    safewrite(&cmd->phiL, sizeof(real), str);
    safewrite(&cmd->phiR, sizeof(real), str);
    //E

    //B Parameters to control histograms and their output files
    safewrite(&cmd->useLogHist, sizeof(bool), str);
    safewrite(&cmd->logHistBinsPD, sizeof(int), str);
    //
    safewrite(&cmd->sizeHistN, sizeof(int), str);
    safewrite(&cmd->rangeN, sizeof(real), str);
    safewrite(&cmd->rminHist, sizeof(real), str);
    safewrite(&cmd->sizeHistPhi, sizeof(int), str);
    //
    nchars = strlen(cmd->histNNFileName) + 1;
    safewrite(&nchars, sizeof(int), str);
    safewrite(cmd->histNNFileName, nchars * sizeof(char), str);

    nchars = strlen(cmd->histXi2pcfFileName) + 1;
    safewrite(&nchars, sizeof(int), str);
    safewrite(cmd->histXi2pcfFileName, nchars * sizeof(char), str);

    nchars = strlen(cmd->histZetaFileName) + 1;
    safewrite(&nchars, sizeof(int), str);
    safewrite(cmd->histZetaFileName, nchars * sizeof(char), str);

    nchars = strlen(cmd->suffixOutFiles) + 1;
    safewrite(&nchars, sizeof(int), str);
    safewrite(cmd->suffixOutFiles, nchars * sizeof(char), str);
    //E

    //B Set of parameters needed to construct a test model
    safewrite(&cmd->seed, sizeof(int), str);

    nchars = strlen(cmd->testmodel) + 1;
    safewrite(&nchars, sizeof(int), str);
    safewrite(cmd->testmodel, nchars * sizeof(char), str);

#ifdef LONGINT
    safewrite(&cmd->nbody, sizeof(INTEGER), str);
#else
    safewrite(&cmd->nbody, sizeof(int), str);
#endif
    safewrite(&cmd->lengthBox, sizeof(real), str);
    //E

    //B Miscellaneous parameters
    nchars = strlen(cmd->script) + 1;
    safewrite(&nchars, sizeof(int), str);
    safewrite(cmd->script, nchars * sizeof(char), str);

#ifdef LONGINT
    safewrite(&cmd->stepState, sizeof(INTEGER), str);
#else
    safewrite(&cmd->stepState, sizeof(int), str);
#endif
    safewrite(&cmd->verbose, sizeof(short), str);
    safewrite(&cmd->verbose_log, sizeof(short), str);
#ifdef OPENMPCODE
    safewrite(&cmd->numthreads, sizeof(int), str);
#endif

    nchars = strlen(cmd->options) + 1;
    safewrite(&nchars, sizeof(int), str);
    safewrite(cmd->options, nchars * sizeof(char), str);
    //E

//B in cmdline_data.h there is an ADDONS ifdef
//  that include the file "globaldefs_include_01.h"
//  here it is the necessary lines...
#ifdef ADDONS
#ifdef BALLS
    safewrite(&cmd->ntosave, sizeof(INTEGER), str);
    safewrite(&cmd->scanLevel, sizeof(int), str);
    //B Root nodes:
    safewrite(&cmd->scanLevelRoot, sizeof(int), str);
    nchars = strlen(cmd->scanLevelMin) + 1;
    safewrite(&nchars, sizeof(int), str);
    safewrite(cmd->scanLevelMin, nchars * sizeof(char), str);
    //E
#endif
#ifdef IOLIB
    nchars = strlen(cmd->columns) + 1;
    safewrite(&nchars, sizeof(int), str);
    safewrite(cmd->columns, nchars * sizeof(char), str);
#endif
#endif // ! ADDONS

    nchars = strlen(cmd->statefile) + 1;
    safewrite(&nchars, sizeof(int), str);
    safewrite(cmd->statefile, nchars * sizeof(char), str);

//E

//E cmdline_data
//=============================================


//=============================================
//B global_data

    safewrite(&gd->cpuinit, sizeof(real), str);         // init at main
    safewrite(&gd->cpurealinit, sizeof(real), str);     // init at main

    //B this headlines are input at the begining of StartRun
    //  no needed to save...
    nchars = strlen(gd->headline0) + 1;
    safewrite(&nchars, sizeof(int), str);
    safewrite(gd->headline0, nchars * sizeof(char), str);

    nchars = strlen(gd->headline1) + 1;
    safewrite(&nchars, sizeof(int), str);
    safewrite(gd->headline1, nchars * sizeof(char), str);

    nchars = strlen(gd->headline2) + 1;
    safewrite(&nchars, sizeof(int), str);
    safewrite(gd->headline2, nchars * sizeof(char), str);

    nchars = strlen(gd->headline3) + 1;
    safewrite(&nchars, sizeof(int), str);
    safewrite(gd->headline3, nchars * sizeof(char), str);
    //E

    nchars = strlen(gd->mode) + 1;
    safewrite(&nchars, sizeof(int), str);
    safewrite(gd->mode, nchars * sizeof(char), str);

    //B Tree
    safewrite(&gd->ncell, sizeof(INTEGER), str);
    safewrite(&gd->tdepth, sizeof(int), str);
    safewrite(&gd->actmax, sizeof(INTEGER), str);
    safewrite(&gd->sameposcount, sizeof(INTEGER), str); // init at StartRun
    // BALLS
    safewrite(&gd->ncccalc, sizeof(INTEGER), str);
    safewrite(&gd->nsmoothcount, sizeof(INTEGER), str);
    //
    safewrite(&gd->nbccalc, sizeof(INTEGER), str);
    safewrite(&gd->nbbcalc, sizeof(INTEGER), str);
    safewrite(&gd->rSize, sizeof(real), str);

    safewrite(&gd->ninfiles, sizeof(int), str);
    for (ifile=0; ifile<gd->ninfiles; ifile++) {
        safewrite(&gd->rSizeTable[ifile], sizeof(real), str);
    }
    for (ifile=0; ifile<gd->ninfiles; ifile++) {
        safewrite(&gd->ncellTable[ifile], sizeof(INTEGER), str);
    }
    for (ifile=0; ifile<gd->ninfiles; ifile++) {
        safewrite(&gd->nbodyTable[ifile], sizeof(INTEGER), str);
    }
    for (ifile=0; ifile<gd->ninfiles; ifile++) {
        safewrite(&gd->tdepthTable[ifile], sizeof(int), str);
    }
    for (ifile=0; ifile<gd->ninfiles; ifile++) {
        safewrite(&gd->nnodescanlevTable[ifile], sizeof(INTEGER), str);
    }
    for (ifile=0; ifile<gd->ninfiles; ifile++) {
        safewrite(&gd->nnodescanlev_rootTable[ifile], sizeof(INTEGER), str);
    }

    safewrite(&gd->cputree, sizeof(real), str);
    //E

    safewrite(&gd->Rcut, sizeof(real), str);
    safewrite(&gd->RcutSq, sizeof(real), str);
    safewrite(&gd->cpusearch, sizeof(real), str);

#ifdef CELLMETHOD
// Cell search
//     vectorI cells;
//    INTEGER *cellList;
//
#endif

    safewrite(&gd->infilefmt_int, sizeof(int), str);


    //B save NDIM
    ndim = NDIM;
    safewrite(&ndim, sizeof(int), str);
    //E

//B here NDIM is used
#ifdef SINGLEP
    DO_COORD(k) {
        safewrite(&gd->Box[k], sizeof(double), str);
    }
#else
    DO_COORD(k) {
        safewrite(&gd->Box[k], sizeof(real), str);
    }
#endif
//E

    safewrite(&gd->searchMethod_int, sizeof(int), str);

    safewrite(&gd->deltaR, sizeof(real), str);
    safewrite(&gd->deltaRmin, sizeof(real), str);
    safewrite(&gd->deltaRmax, sizeof(real), str);
    if (cmd->useLogHist) {
        if (cmd->rminHist==0) {
            //B rminHist = 0 not allowed when useLogHist is true
            //  therefore this segment does not ocurr...
            for (n=1; n<=cmd->sizeHistN; n++) {
                safewrite(&gd->deltaRV[n], sizeof(real), str);
            }
            //E
        } else {
            for (n=1; n<=cmd->sizeHistN; n++) {
                safewrite(&gd->deltaRV[n], sizeof(real), str);
            }
            for (n=1; n<=cmd->sizeHistN-1; n++) {
                safewrite(&gd->ddeltaRV[n], sizeof(real), str);
            }
        }
    }

    safewrite(&gd->deltaPhi, sizeof(real), str);

    //B File pointers:
    nchars = strlen(gd->logfilePath) + 1;
    safewrite(&nchars, sizeof(int), str);
    safewrite(gd->logfilePath, nchars * sizeof(char), str);

    nchars = strlen(gd->outputDir) + 1;
    safewrite(&nchars, sizeof(int), str);
    safewrite(gd->outputDir, nchars * sizeof(char), str);

    nchars = strlen(gd->tmpDir) + 1;
    safewrite(&nchars, sizeof(int), str);
    safewrite(gd->tmpDir, nchars * sizeof(char), str);

    nchars = strlen(gd->fpfnameOutputFileName) + 1;
    safewrite(&nchars, sizeof(int), str);
    safewrite(gd->fpfnameOutputFileName, nchars * sizeof(char), str);

    nchars = strlen(gd->fpfnamehistNNFileName) + 1;
    safewrite(&nchars, sizeof(int), str);
    safewrite(gd->fpfnamehistNNFileName, nchars * sizeof(char), str);

    nchars = strlen(gd->fpfnamehistCFFileName) + 1;
    safewrite(&nchars, sizeof(int), str);
    safewrite(gd->fpfnamehistCFFileName, nchars * sizeof(char), str);

    nchars = strlen(gd->fpfnamehistrBinsFileName) + 1;
    safewrite(&nchars, sizeof(int), str);
    safewrite(gd->fpfnamehistrBinsFileName, nchars * sizeof(char), str);

    nchars = strlen(gd->fpfnamehistXi2pcfFileName) + 1;
    safewrite(&nchars, sizeof(int), str);
    safewrite(gd->fpfnamehistXi2pcfFileName, nchars * sizeof(char), str);

    nchars = strlen(gd->fpfnamehistZetaGFileName) + 1;
    safewrite(&nchars, sizeof(int), str);
    safewrite(gd->fpfnamehistZetaGFileName, nchars * sizeof(char), str);

    nchars = strlen(gd->fpfnamehistZetaGmFileName) + 1;
    safewrite(&nchars, sizeof(int), str);
    safewrite(gd->fpfnamehistZetaGmFileName, nchars * sizeof(char), str);

    nchars = strlen(gd->fpfnamehistZetaMFileName) + 1;
    safewrite(&nchars, sizeof(int), str);
    safewrite(gd->fpfnamehistZetaMFileName, nchars * sizeof(char), str);

    nchars = strlen(gd->fpfnamemhistZetaMFileName) + 1;
    safewrite(&nchars, sizeof(int), str);
    safewrite(gd->fpfnamemhistZetaMFileName, nchars * sizeof(char), str);

    nchars = strlen(gd->fpfnameCPUFileName) + 1;
    safewrite(&nchars, sizeof(int), str);
    safewrite(gd->fpfnameCPUFileName, nchars * sizeof(char), str);
    //E

    nchars = strlen(gd->model_comment) + 1;
    safewrite(&nchars, sizeof(int), str);
    safewrite(gd->model_comment, nchars * sizeof(char), str);

    nchars = strlen(gd->input_comment) + 1;
    safewrite(&nchars, sizeof(int), str);
    safewrite(gd->input_comment, nchars * sizeof(char), str);

    nchars = strlen(gd->output_comment) + 1;
    safewrite(&nchars, sizeof(int), str);
    safewrite(gd->output_comment, nchars * sizeof(char), str);

    //B save-restore
    safewrite(&gd->stopflag, sizeof(int), str);         // init at StartRun
    safewrite(&gd->ip, sizeof(INTEGER), str);
    //E

    safewrite(&gd->cputotalinout, sizeof(real), str);   // init at StartRun
    safewrite(&gd->cputotal, sizeof(real), str);        // init at StartRun

    safewrite(&gd->bytes_tot, sizeof(INTEGER), str);    // init at StartRun
    safewrite(&gd->bytes_tot_cells, sizeof(INTEGER), str);

    safewrite(&gd->nbodySel, sizeof(INTEGER), str);

    safewrite(&gd->nbodysm, sizeof(INTEGER), str);
    safewrite(&gd->nbodybf, sizeof(INTEGER), str);

    safewrite(&gd->bh86, sizeof(bool), str);
    safewrite(&gd->sw94, sizeof(bool), str);

    safewrite(&gd->i_deltaR, sizeof(real), str);

#ifdef ADDONS
//#include "global_data_include.h"
#ifdef NMultipoles
//B Already safewrite when KDTREEOMP is defined
//#ifndef KDTREEOMP             // NMultipoles Deactivated in KDTREEOMP
//#else // ! KDTREEOMP      // NMultipoles Deactivated in KDTREEOMP
//#endif // ! KDTREEOMP
#endif // ! NMultipoles

    nchars = strlen(gd->fpfnamehistN2pcfFileName) + 1;
    safewrite(&nchars, sizeof(int), str);
    safewrite(gd->fpfnamehistN2pcfFileName, nchars * sizeof(char), str);

#endif // ! ADDONS


    nchars = strlen(gd->fnameData_kd) + 1;
    safewrite(&nchars, sizeof(int), str);
    safewrite(gd->fnameData_kd, nchars * sizeof(char), str);

    nchars = strlen(gd->fnameOut_kd) + 1;
    safewrite(&nchars, sizeof(int), str);
    safewrite(gd->fnameOut_kd, nchars * sizeof(char), str);

    safewrite(&gd->input_format_kd, sizeof(int), str);
    safewrite(&gd->use_tree_kd, sizeof(int), str);
    safewrite(&gd->max_tree_order_kd, sizeof(int), str);
    safewrite(&gd->max_tree_nparts_kd, sizeof(int), str);
    safewrite(&gd->use_pm_kd, sizeof(int), str);
    safewrite(&gd->n_objects_kd, sizeof(INTEGER), str);
    safewrite(&gd->l_box_kd, sizeof(float), str);
    safewrite(&gd->l_box_half_kd, sizeof(float), str);

    for (ifile=0; ifile<gd->ninfiles; ifile++) {
        nchars = strlen(gd->infilenames[ifile]) + 1;
        safewrite(&nchars, sizeof(int), str);
        safewrite(gd->infilenames[ifile], nchars * sizeof(char), str);
    }
    for (ifile=0; ifile<gd->ninfiles; ifile++) {
        nchars = strlen(gd->infilefmtname[ifile]) + 1;
        safewrite(&nchars, sizeof(int), str);
        safewrite(gd->infilefmtname[ifile], nchars * sizeof(char), str);
    }
    for (ifile=0; ifile<gd->ninfiles; ifile++) {
        safewrite(&gd->iCatalogs[ifile], sizeof(int), str);
    }

    for (i=0; i<1; i++) {
        safewrite(&gd->nsmooth[i], sizeof(int), str);
    }
    safewrite(&gd->nnode, sizeof(INTEGER), str);
    safewrite(&gd->rnnode, sizeof(INTEGER), str);

    for (i=0; i<2; i++) {
        safewrite(&gd->scanLevelMin[i], sizeof(int), str);
    }

    nchars = strlen(gd->nodesfilePath) + 1;
    safewrite(&nchars, sizeof(int), str);
    safewrite(gd->nodesfilePath, nchars * sizeof(char), str);

    safewrite(&gd->nnodescanlev, sizeof(int), str);
    //B Root nodes:
    safewrite(&gd->nnodescanlev_root, sizeof(int), str);

//B must see MAXLEV definition in global_data.h
#define MAXLEVEL  32
    for (i=0; i<MAXLEVEL; i++) {
        safewrite(&gd->Rcell[i], sizeof(real), str);
    }
#undef MAXLEVEL
//E

    nchars = strlen(gd->bodiesfilePath) + 1;
    safewrite(&nchars, sizeof(int), str);
    safewrite(gd->bodiesfilePath, nchars * sizeof(char), str);

    safewrite(&gd->flagSmoothCellMin, sizeof(bool), str);
    safewrite(&gd->flagSmooth, sizeof(bool), str);
    safewrite(&gd->flagSetNbNoSel, sizeof(bool), str);

    //B BUCKET
    for (i=0; i<2; i++) {
        safewrite(&gd->rminCell[i], sizeof(real), str);
    }
    //E
    for (i=0; i<1; i++) {
        safewrite(&gd->rsmooth[i], sizeof(real), str);
    }

    safewrite(&gd->rsmoothFlag, sizeof(bool), str);

    safewrite(&gd->irsmooth, sizeof(int), str);

#ifdef ADDONS
//#include "globaldefs_include_03.h"
#ifdef BALLS
//#include "globaldefs_balls_omp_02.h"  // empty defs
#endif

#ifdef IOLIB
//#include "global_data_iolib.h"
    // pos, kappa, gamma1, gamma2, weight, seven places maximum
    // use in multi-columns-ascii and cfitsio
    // default is 3D: three pos and one convergence (kappa)
    for (i=0; i<6; i++) {
        safewrite(&gd->columns[i], sizeof(int), str);
    }
#endif
#endif // ! ADDONS

//E global_data
//=============================================


//=============================================
//B globaldefs.h

#ifdef USEGSL
//    gsl_rng * r;            // It is used r_gsl in globaldefs.h. Check!!!
#else
    safewrite(&idum, sizeof(long), str);
#endif

//B OpenMP section
//  Timing variables
#ifdef OPENMPCODE
    safewrite(&relstart, sizeof(double), str);
    safewrite(&relend, sizeof(double), str);
    safewrite(&absstart, sizeof(double), str);
    safewrite(&absend, sizeof(double), str);
#else
    safewrite(&relstart, sizeof(time_t), str);
    safewrite(&relend, sizeof(time_t), str);
    safewrite(&absstart, sizeof(time_t), str);
    safewrite(&absend, sizeof(time_t), str);
#endif
//E


//B saving catalogs
    for (ifile=0; ifile<gd->ninfiles; ifile++) {
        safewrite(bodytable[ifile], gd->nbodyTable[ifile] * sizeof(body), str);
    }
//E

// to avoid conflicts allocating memory...
//    tree must be reconstructed when restoring the run

#ifdef ADDONS
//#include "globaldefs_include_04.h"
//#include "cballsio_gadget_00.h"
//#include "globaldefs_balls_omp_04.h"
#endif

#ifdef ADDONS
//#include "globaldefs_include_05.h"
#endif

//E globaldefs.h
//=============================================



//=============================================
//B histograms

    for (k=1; k<=cmd->sizeHistN; k++){
        rval = gd->histNNSubXi2pcf[k];
        safewrite(&rval, sizeof(real), str);
    }

    for (k=1; k<=cmd->sizeHistN; k++){
        rval = gd->histNN[k];
        safewrite(&rval, sizeof(real), str);
    }

    for (k=1; k<=cmd->sizeHistN; k++){
        rval = gd->histXi2pcf[k];
        safewrite(&rval, sizeof(real), str);
    }

    for (m=1; m<=cmd->mChebyshev+1; m++)
        for (n=1; n<=cmd->sizeHistN; n++)
            for (k=1; k<=cmd->sizeHistN; k++) {
                rval = gd->histZetaMcos[m][n][k];
                safewrite(&rval, sizeof(real), str);
            }

    for (m=1; m<=cmd->mChebyshev+1; m++)
        for (n=1; n<=cmd->sizeHistN; n++)
            for (k=1; k<=cmd->sizeHistN; k++) {
                rval = gd->histZetaMsin[m][n][k];
                safewrite(&rval, sizeof(real), str);
            }

    for (m=1; m<=cmd->mChebyshev+1; m++)
        for (n=1; n<=cmd->sizeHistN; n++)
            for (k=1; k<=cmd->sizeHistN; k++) {
                rval = gd->histZetaMsincos[m][n][k];
                safewrite(&rval, sizeof(real), str);
            }

    for (m=1; m<=cmd->mChebyshev+1; m++)
        for (n=1; n<=cmd->sizeHistN; n++)
            for (k=1; k<=cmd->sizeHistN; k++) {
                rval = gd->histZetaMcossin[m][n][k];
                safewrite(&rval, sizeof(real), str);
            }

#ifdef ADDONS
//#include "global_data_include.h"
#ifdef NMultipoles
//B Already safewrite when KDTREEOMP is defined
//#ifndef KDTREEOMP             // NMultipoles Deactivated in KDTREEOMP
//#else // ! KDTREEOMP      // NMultipoles Deactivated in KDTREEOMP
//#endif // ! KDTREEOMP
#endif // ! NMultipoles
#endif // ! ADDONS

//E histograms
//=============================================

    fclose(str);
    gd->cputotalinout += CPUTIME - cpustart;

}

/*
 Restore state of the run routine:

 Arguments:
    * `pattern`: name pattern of the file to restore the state of the run
 Global tructures used: gd, cmd
 Counting encounters (in global gd): nbbcalc, nbccalc, ncccalc
 Return (the error status):
    int SUCCESS or FAILURE
 */
global int restorestate(struct cmdline_data* cmd,
                        struct  global_data* gd,
                        string file)
{
// every item you add to cmdline_data, global_data structures
//  or to globaldefs.h file should you consider necessary to
//  recover the run, must be added to this routine.
// Also all items must be read by restorestate routine
//  in the same order as was given in savestate routine

    int m,n,k;
    double rval;
    stream str;
    int nchars;
    int ndim;
    int ifile;
    int i;
    string program, version;

    //B must be after reading from restore file
    gd->model_comment = "start from a restore file";
    //E

    double cpustart = CPUTIME;
    //B header: program name and version
    str = stropen(file, "r");
    saferead(&nchars, sizeof(int), str);
    program = (string) allocate(nchars * sizeof(char));
    saferead(program, nchars * sizeof(char), str);
    saferead(&nchars, sizeof(int), str);
    version = (string) allocate(nchars * sizeof(char));
    saferead(version, nchars * sizeof(char), str);
    if (! streq(program, getargv0()) ||
        ! streq(version, getversion())) {
        printf("warning: state file may be outdated (%s %s)\n\n",
               program, version);
    }
    //E

//=============================================
//B cmdline_data

    //B Parameters related to the searching method
    saferead(&nchars, sizeof(int), str);
    cmd->searchMethod = (string) allocate(nchars * sizeof(char));
    saferead(cmd->searchMethod, nchars * sizeof(char), str);

    saferead(&cmd->mChebyshev, sizeof(int), str);

    saferead(&nchars, sizeof(int), str);
    cmd->nsmooth = (string) allocate(nchars * sizeof(char));
    saferead(cmd->nsmooth, nchars * sizeof(char), str);

    saferead(&nchars, sizeof(int), str);
    cmd->rsmooth = (string) allocate(nchars * sizeof(char));
    saferead(cmd->rsmooth, nchars * sizeof(char), str);

    saferead(&cmd->theta, sizeof(real), str);
    saferead(&cmd->computeTPCF, sizeof(bool), str);
    saferead(&cmd->computeShearCF, sizeof(bool), str);
    saferead(&cmd->usePeriodic, sizeof(bool), str);
    //E
    
    //B Parameters about the I/O file(s)
    // Input catalog parameters
    saferead(&nchars, sizeof(int), str);
    cmd->infile = (string) allocate(nchars * sizeof(char));
    saferead(cmd->infile, nchars * sizeof(char), str);

    saferead(&nchars, sizeof(int), str);
    cmd->infilefmt = (string) allocate(nchars * sizeof(char));
    saferead(cmd->infilefmt, nchars * sizeof(char), str);

    saferead(&nchars, sizeof(int), str);
    cmd->iCatalogs = (string) allocate(nchars * sizeof(char));
    saferead(cmd->iCatalogs, nchars * sizeof(char), str);

    // Output parameters
    saferead(&nchars, sizeof(int), str);
    cmd->rootDir = (string) allocate(nchars * sizeof(char));
    saferead(cmd->rootDir, nchars * sizeof(char), str);

    saferead(&nchars, sizeof(int), str);
    cmd->outfile = (string) allocate(nchars * sizeof(char));
    saferead(cmd->outfile, nchars * sizeof(char), str);

    saferead(&nchars, sizeof(int), str);
    cmd->outfilefmt = (string) allocate(nchars * sizeof(char));
    saferead(cmd->outfilefmt, nchars * sizeof(char), str);

    // Parameters to set a region in the sky, for example for Takahasi data set
    saferead(&cmd->thetaL, sizeof(real), str);
    saferead(&cmd->thetaR, sizeof(real), str);
    saferead(&cmd->phiL, sizeof(real), str);
    saferead(&cmd->phiR, sizeof(real), str);
    //E

    //B Parameters to control histograms and their output files
    saferead(&cmd->useLogHist, sizeof(bool), str);
    saferead(&cmd->logHistBinsPD, sizeof(int), str);
    //
    saferead(&cmd->sizeHistN, sizeof(int), str);
    saferead(&cmd->rangeN, sizeof(real), str);
    saferead(&cmd->rminHist, sizeof(real), str);
    saferead(&cmd->sizeHistPhi, sizeof(int), str);
    //
    saferead(&nchars, sizeof(int), str);
    cmd->histNNFileName = (string) allocate(nchars * sizeof(char));
    saferead(cmd->histNNFileName, nchars * sizeof(char), str);

    saferead(&nchars, sizeof(int), str);
    cmd->histXi2pcfFileName = (string) allocate(nchars * sizeof(char));
    saferead(cmd->histXi2pcfFileName, nchars * sizeof(char), str);

    saferead(&nchars, sizeof(int), str);
    cmd->histZetaFileName = (string) allocate(nchars * sizeof(char));
    saferead(cmd->histZetaFileName, nchars * sizeof(char), str);

    saferead(&nchars, sizeof(int), str);
    cmd->suffixOutFiles = (string) allocate(nchars * sizeof(char));
    saferead(cmd->suffixOutFiles, nchars * sizeof(char), str);
    //E

    //B Set of parameters needed to construct a test model
    saferead(&cmd->seed, sizeof(int), str);

    saferead(&nchars, sizeof(int), str);
    cmd->testmodel = (string) allocate(nchars * sizeof(char));
    saferead(cmd->testmodel, nchars * sizeof(char), str);

#ifdef LONGINT
    saferead(&cmd->nbody, sizeof(INTEGER), str);
#else
    saferead(&cmd->nbody, sizeof(int), str);
#endif
    saferead(&cmd->lengthBox, sizeof(real), str);
    //E

    //B Miscellaneous parameters
    saferead(&nchars, sizeof(int), str);
    cmd->script = (string) allocate(nchars * sizeof(char));
    saferead(cmd->script, nchars * sizeof(char), str);

#ifdef LONGINT
    saferead(&cmd->stepState, sizeof(INTEGER), str);
#else
    saferead(&cmd->stepState, sizeof(int), str);
#endif
    saferead(&cmd->verbose, sizeof(short), str);
    saferead(&cmd->verbose_log, sizeof(short), str);
#ifdef OPENMPCODE
    saferead(&cmd->numthreads, sizeof(int), str);
#endif

    saferead(&nchars, sizeof(int), str);
    cmd->options = (string) allocate(nchars * sizeof(char));
    saferead(cmd->options, nchars * sizeof(char), str);
    //E

//B in cmdline_data.h there is an ADDONS ifdef
//  that include the file "globaldefs_include_01.h"
//  here it is the necessary lines...
#ifdef ADDONS
#ifdef BALLS
    saferead(&cmd->ntosave, sizeof(INTEGER), str);
    saferead(&cmd->scanLevel, sizeof(int), str);
    //B Root nodes:
    saferead(&cmd->scanLevelRoot, sizeof(int), str);
    saferead(&nchars, sizeof(int), str);
    cmd->scanLevelMin = (string) allocate(nchars * sizeof(char));
    saferead(cmd->scanLevelMin, nchars * sizeof(char), str);
    //E
#endif
#ifdef IOLIB
    saferead(&nchars, sizeof(int), str);
    cmd->columns = (string) allocate(nchars * sizeof(char));
    saferead(cmd->columns, nchars * sizeof(char), str);
#endif
#endif // ! ADDONS

    saferead(&nchars, sizeof(int), str);
    cmd->statefile = (string) allocate(nchars * sizeof(char));
    saferead(cmd->statefile, nchars * sizeof(char), str);

//E

//E cmdline_data
//=============================================


//=============================================
//B global_data

    saferead(&gd->cpuinit, sizeof(real), str);
    saferead(&gd->cpurealinit, sizeof(real), str);

    saferead(&nchars, sizeof(int), str);
    gd->headline0 = (string) allocate(nchars * sizeof(char));
    saferead(gd->headline0, nchars * sizeof(char), str);

    saferead(&nchars, sizeof(int), str);
    gd->headline1 = (string) allocate(nchars * sizeof(char));
    saferead(gd->headline1, nchars * sizeof(char), str);

    saferead(&nchars, sizeof(int), str);
    gd->headline2 = (string) allocate(nchars * sizeof(char));
    saferead(gd->headline2, nchars * sizeof(char), str);

    saferead(&nchars, sizeof(int), str);
    gd->headline3 = (string) allocate(nchars * sizeof(char));
    saferead(gd->headline3, nchars * sizeof(char), str);

    saferead(&nchars, sizeof(int), str);
    saferead(gd->mode, nchars * sizeof(char), str);

    //B Tree
    saferead(&gd->ncell, sizeof(INTEGER), str);
    saferead(&gd->tdepth, sizeof(int), str);
    saferead(&gd->actmax, sizeof(INTEGER), str);
    saferead(&gd->sameposcount, sizeof(INTEGER), str);
    // BALLS
    saferead(&gd->ncccalc, sizeof(INTEGER), str);
    saferead(&gd->nsmoothcount, sizeof(INTEGER), str);
    //
    saferead(&gd->nbccalc, sizeof(INTEGER), str);
    saferead(&gd->nbbcalc, sizeof(INTEGER), str);
    saferead(&gd->rSize, sizeof(real), str);

    saferead(&gd->ninfiles, sizeof(int), str);
    for (ifile=0; ifile<gd->ninfiles; ifile++) {
        saferead(&gd->rSizeTable[ifile], sizeof(real), str);
    }
    for (ifile=0; ifile<gd->ninfiles; ifile++) {
        saferead(&gd->ncellTable[ifile], sizeof(INTEGER), str);
    }
    for (ifile=0; ifile<gd->ninfiles; ifile++) {
        saferead(&gd->nbodyTable[ifile], sizeof(INTEGER), str);
    }
    for (ifile=0; ifile<gd->ninfiles; ifile++) {
        saferead(&gd->tdepthTable[ifile], sizeof(int), str);
    }
    for (ifile=0; ifile<gd->ninfiles; ifile++) {
        saferead(&gd->nnodescanlevTable[ifile], sizeof(INTEGER), str);
    }
    for (ifile=0; ifile<gd->ninfiles; ifile++) {
        saferead(&gd->nnodescanlev_rootTable[ifile], sizeof(INTEGER), str);
    }

    saferead(&gd->cputree, sizeof(real), str);
    //E

    saferead(&gd->Rcut, sizeof(real), str);
    saferead(&gd->RcutSq, sizeof(real), str);
    saferead(&gd->cpusearch, sizeof(real), str);

#ifdef CELLMETHOD
// Cell search
//     vectorI cells;
//    INTEGER *cellList;
//
#endif

    saferead(&gd->infilefmt_int, sizeof(int), str);

    
    //B read NDIM
    saferead(&ndim, sizeof(int), str);
    if (ndim != NDIM)
        error(
              "\nrestorestate: DIM of restore file is different than NDIM\n\n"
              );
    //

//B here NDIM is used
#ifdef SINGLEP
    DO_COORD(k) {
        saferead(&gd->Box[k], sizeof(double), str);
    }
#else
    DO_COORD(k) {
        saferead(&gd->Box[k], sizeof(real), str);
    }
#endif
//E

    saferead(&gd->searchMethod_int, sizeof(int), str);

    saferead(&gd->deltaR, sizeof(real), str);
    saferead(&gd->deltaRmin, sizeof(real), str);
    saferead(&gd->deltaRmax, sizeof(real), str);
    if (cmd->useLogHist) {
        if (cmd->rminHist==0) {
            //B rminHist = 0 not allowed when useLogHist is true
            //  therefore this segment does not ocurr...
            gd->deltaRV = dvector(1,cmd->sizeHistN);
            for (n=1; n<=cmd->sizeHistN; n++) {
                saferead(&gd->deltaRV[n], sizeof(real), str);
            }
            //E
        } else {
            //B allocated after startrun_memoryAllocation
            //  deallocate before deallocate arrays
            //  in startrun_memoryAllocation
            gd->deltaRV = dvector(1,cmd->sizeHistN);
            gd->ddeltaRV = dvector(1,cmd->sizeHistN-1);
            //E
            for (n=1; n<=cmd->sizeHistN; n++) {
                saferead(&gd->deltaRV[n], sizeof(real), str);
            }
            for (n=1; n<=cmd->sizeHistN-1; n++) {
                saferead(&gd->ddeltaRV[n], sizeof(real), str);
            }
        }
    }

    saferead(&gd->deltaPhi, sizeof(real), str);

    //B File pointers:
    saferead(&nchars, sizeof(int), str);
    saferead(gd->logfilePath, nchars * sizeof(char), str);
    
    saferead(&nchars, sizeof(int), str);
    saferead(gd->outputDir, nchars * sizeof(char), str);

    saferead(&nchars, sizeof(int), str);
    saferead(gd->tmpDir, nchars * sizeof(char), str);

    saferead(&nchars, sizeof(int), str);
    saferead(gd->fpfnameOutputFileName, nchars * sizeof(char), str);

    saferead(&nchars, sizeof(int), str);
    saferead(gd->fpfnamehistNNFileName, nchars * sizeof(char), str);

    saferead(&nchars, sizeof(int), str);
    saferead(gd->fpfnamehistCFFileName, nchars * sizeof(char), str);

    saferead(&nchars, sizeof(int), str);
    saferead(gd->fpfnamehistrBinsFileName, nchars * sizeof(char), str);

    saferead(&nchars, sizeof(int), str);
    saferead(gd->fpfnamehistXi2pcfFileName, nchars * sizeof(char), str);

    saferead(&nchars, sizeof(int), str);
    saferead(gd->fpfnamehistZetaGFileName, nchars * sizeof(char), str);

    saferead(&nchars, sizeof(int), str);
    saferead(gd->fpfnamehistZetaGmFileName, nchars * sizeof(char), str);

    saferead(&nchars, sizeof(int), str);
    saferead(gd->fpfnamehistZetaMFileName, nchars * sizeof(char), str);

    saferead(&nchars, sizeof(int), str);
    saferead(gd->fpfnamemhistZetaMFileName, nchars * sizeof(char), str);

    saferead(&nchars, sizeof(int), str);
    saferead(gd->fpfnameCPUFileName, nchars * sizeof(char), str);
    //E

    saferead(&nchars, sizeof(int), str);
    gd->model_comment = (string) allocate(nchars * sizeof(char));
    saferead(gd->model_comment, nchars * sizeof(char), str);

    //B must be after reading from restore file
    verb_print(cmd->verbose,
            "\nrestorestate: model comment before reading restore file %s\n",
            gd->model_comment);
    gd->model_comment = "start from a restore file";
    verb_print(cmd->verbose,
            "restorestate: model comment after correcting %s\n",
            gd->model_comment);
    //E

    saferead(&nchars, sizeof(int), str);
    gd->input_comment = (string) allocate(nchars * sizeof(char));
    saferead(gd->input_comment, nchars * sizeof(char), str);

    saferead(&nchars, sizeof(int), str);
    gd->output_comment = (string) allocate(nchars * sizeof(char));
    saferead(gd->output_comment, nchars * sizeof(char), str);


    //B save-restore
    saferead(&gd->stopflag, sizeof(int), str);
    saferead(&gd->ip, sizeof(INTEGER), str);
    //E

    saferead(&gd->cputotalinout, sizeof(double), str);
    saferead(&gd->cputotal, sizeof(double), str);

    saferead(&gd->bytes_tot, sizeof(INTEGER), str);
    saferead(&gd->bytes_tot_cells, sizeof(INTEGER), str);

    saferead(&gd->nbodySel, sizeof(INTEGER), str);

    saferead(&gd->nbodysm, sizeof(INTEGER), str);
    saferead(&gd->nbodybf, sizeof(INTEGER), str);

    saferead(&gd->bh86, sizeof(bool), str);
    saferead(&gd->sw94, sizeof(bool), str);

    saferead(&gd->i_deltaR, sizeof(real), str);

#ifdef ADDONS
//#include "global_data_include.h"
#ifdef NMultipoles
//B Already safewrite when KDTREEOMP is defined
//#ifndef KDTREEOMP             // NMultipoles Deactivated in KDTREEOMP
//B histograms
//E histograms

#ifdef USEGSL
//    gsl_matrix_complex *NhistXi_gsl;
#endif

#endif // ! NMultipoles
    
    saferead(&nchars, sizeof(int), str);
    saferead(gd->fpfnamehistN2pcfFileName, nchars * sizeof(char), str);

#endif // ! ADDONS

    saferead(&nchars, sizeof(int), str);
    saferead(gd->fnameData_kd, nchars * sizeof(char), str);

    saferead(&nchars, sizeof(int), str);
    saferead(gd->fnameOut_kd, nchars * sizeof(char), str);

    saferead(&gd->input_format_kd, sizeof(int), str);
    saferead(&gd->use_tree_kd, sizeof(int), str);
    saferead(&gd->max_tree_order_kd, sizeof(int), str);
    saferead(&gd->max_tree_nparts_kd, sizeof(int), str);
    saferead(&gd->use_pm_kd, sizeof(int), str);
    saferead(&gd->n_objects_kd, sizeof(INTEGER), str);
    saferead(&gd->l_box_kd, sizeof(float), str);
    saferead(&gd->l_box_half_kd, sizeof(float), str);

    for (ifile=0; ifile<gd->ninfiles; ifile++) {
        saferead(&nchars, sizeof(int), str);
        gd->infilenames[ifile] = (string) allocate(nchars * sizeof(char));
        saferead(gd->infilenames[ifile], nchars * sizeof(char), str);
    }
    for (ifile=0; ifile<gd->ninfiles; ifile++) {
        saferead(&nchars, sizeof(int), str);
        gd->infilefmtname[ifile] = (string) allocate(nchars * sizeof(char));
        saferead(gd->infilefmtname[ifile], nchars * sizeof(char), str);
    }
    for (ifile=0; ifile<gd->ninfiles; ifile++) {
        saferead(&gd->iCatalogs[ifile], sizeof(int), str);
    }

    //B setting iCatalogs for ninfiles=1
    if (gd->ninfiles==1) {
        (gd->iCatalogs[0]) = 0;
        if (cmd->verbose>=2)
            verb_print(cmd->verbose,
                        "restorestate: iCatalogs final values: %d\n",
                        gd->iCatalogs[0]);
        (gd->iCatalogs[1]) = 0;
        if (cmd->verbose>=2)
            verb_print(cmd->verbose,
                        "restorestate: iCatalogs final values: %d\n",
                        gd->iCatalogs[1]);
    } else {
        for (i=0; i<gd->ninfiles; i++) {
            if (cmd->verbose>=2)
                verb_print(cmd->verbose,
                        "restorestate: iCatalogs final values: %d\n",
                        gd->iCatalogs[i]);
        }
    }
    //E

    for (i=0; i<1; i++) {
        saferead(&gd->nsmooth[i], sizeof(int), str);
    }
    saferead(&gd->nnode, sizeof(INTEGER), str);
    saferead(&gd->rnnode, sizeof(INTEGER), str);

    for (i=0; i<2; i++) {
        saferead(&gd->scanLevelMin[i], sizeof(int), str);
    }

    saferead(&nchars, sizeof(int), str);
    saferead(gd->nodesfilePath, nchars * sizeof(char), str);

    saferead(&gd->nnodescanlev, sizeof(int), str);
    //B Root nodes:
    saferead(&gd->nnodescanlev_root, sizeof(int), str);

//B must see MAXLEV definition in global_data.h
#define MAXLEVEL  32
    for (i=0; i<MAXLEVEL; i++) {
        saferead(&gd->Rcell[i], sizeof(real), str);
    }
#undef MAXLEVEL
//E

    saferead(&nchars, sizeof(int), str);
    saferead(gd->bodiesfilePath, nchars * sizeof(char), str);

    saferead(&gd->flagSmoothCellMin, sizeof(bool), str);
    saferead(&gd->flagSmooth, sizeof(bool), str);
    saferead(&gd->flagSetNbNoSel, sizeof(bool), str);

    //B BUCKET
    for (i=0; i<2; i++) {
        saferead(&gd->rminCell[i], sizeof(real), str);
    }
    //E
    for (i=0; i<1; i++) {
        saferead(&gd->rsmooth[i], sizeof(real), str);
    }

    saferead(&gd->rsmoothFlag, sizeof(bool), str);
    saferead(&gd->irsmooth, sizeof(int), str);

#ifdef ADDONS
//#include "globaldefs_include_03.h"
#ifdef BALLS
//#include "globaldefs_balls_omp_02.h"  // empty defs
#endif

#ifdef IOLIB
//#include "global_data_iolib.h"
    // pos, kappa, gamma1, gamma2, weight, seven places maximum
    // use in multi-columns-ascii and cfitsio
    // default is 3D: three pos and one convergence (kappa)
    for (i=0; i<6; i++) {
        saferead(&gd->columns[i], sizeof(int), str);
    }
#endif
#endif // ! ADDONS

//E global_data
//=============================================


//=============================================
//B globaldefs.h

#ifdef USEGSL
//    gsl_rng * r;            // It is used r_gsl in globaldefs.h. Check!!!
#else
    saferead(&idum, sizeof(long), str);
#endif

//B OpenMP section
//  Timing variables
#ifdef OPENMPCODE
    saferead(&relstart, sizeof(double), str);
    saferead(&relend, sizeof(double), str);
    saferead(&absstart, sizeof(double), str);
    saferead(&absend, sizeof(double), str);
#else
    saferead(&relstart, sizeof(time_t), str);
    saferead(&relend, sizeof(time_t), str);
    saferead(&absstart, sizeof(time_t), str);
    saferead(&absend, sizeof(time_t), str);
#endif
//E

    //B saving catalogs
    for (ifile=0; ifile<gd->ninfiles; ifile++) {
        bodytable[ifile] =
                    (bodyptr) allocate(gd->nbodyTable[ifile] * sizeof(body));
        saferead(bodytable[ifile], gd->nbodyTable[ifile] * sizeof(body), str);
    }
    //E
    
// to avoid conflicts allocating memory... tree must be reconstructed...
    
#ifdef ADDONS
//#include "globaldefs_include_04.h"
//#include "cballsio_gadget_00.h"
//#include "globaldefs_balls_omp_04.h"
#endif

#ifdef ADDONS
//#include "globaldefs_include_05.h"
#endif

//E globaldefs.h
//=============================================


//=============================================
//B histograms

// move this line to segment after restorestate has been
//  read cmd parameters needed to allocate histograms...
    class_call_cballs(startrun_memoryAllocation(cmd, gd), errmsg, errmsg);

    for (k=1; k<=cmd->sizeHistN; k++) {
        saferead(&rval, sizeof(double), str);
        gd->histNNSubXi2pcf[k] = rval;
    }
    for (k=1; k<=cmd->sizeHistN; k++) {
        saferead(&rval, sizeof(double), str);
        gd->histNN[k] = rval;
    }
    for (k=1; k<=cmd->sizeHistN; k++) {
        saferead(&rval, sizeof(double), str);
        gd->histXi2pcf[k] = rval;
    }

    for (m=1; m<=cmd->mChebyshev+1; m++)
        for (n=1; n<=cmd->sizeHistN; n++)
            for (k=1; k<=cmd->sizeHistN; k++) {
                saferead(&rval, sizeof(double), str);
                gd->histZetaMcos[m][n][k] = rval;
            }

    for (m=1; m<=cmd->mChebyshev+1; m++)
        for (n=1; n<=cmd->sizeHistN; n++)
            for (k=1; k<=cmd->sizeHistN; k++) {
                saferead(&rval, sizeof(double), str);
                gd->histZetaMsin[m][n][k] = rval;
            }

    for (m=1; m<=cmd->mChebyshev+1; m++)
        for (n=1; n<=cmd->sizeHistN; n++)
            for (k=1; k<=cmd->sizeHistN; k++) {
                saferead(&rval, sizeof(double), str);
                gd->histZetaMsincos[m][n][k] = rval;
            }

    for (m=1; m<=cmd->mChebyshev+1; m++)
        for (n=1; n<=cmd->sizeHistN; n++)
            for (k=1; k<=cmd->sizeHistN; k++) {
                saferead(&rval, sizeof(double), str);
                gd->histZetaMcossin[m][n][k] = rval;
            }


#ifdef ADDONS
//#include "global_data_include.h"
#ifdef NMultipoles
//B Already safewrite when KDTREEOMP is defined
//#ifndef KDTREEOMP             // NMultipoles Deactivated in KDTREEOMP

#ifdef USEGSL
//    gsl_matrix_complex *NhistXi_gsl;
#endif

//B To save total 3pcf... and savewrite if necessary
//B look for where these definitions are used...
//    real **NhistZetaGcos;
//    real **NhistZetaGsin;
//E
//E

//#else // ! KDTREEOMP      // NMultipoles Deactivated in KDTREEOMP
//#endif // ! KDTREEOMP
#endif // ! NMultipoles

#endif // ! ADDONS

//E histograms
//=============================================


    fclose(str);
    gd->cputotalinout += CPUTIME - cpustart;

    return SUCCESS;
}

#undef savestatetmp


#define FMTT    "%-35s = %s\n"
#define FMTTS    "%-35s = \"%s\"\n"
#define FMTI    "%-35s = %d\n"
#define FMTIL    "%-35s = %ld\n"
#define FMTR    "%-35s = %g\n"

/*
 Print parameters cmd, gd and global after
    restore the state of the run routine has been called:

 Arguments:
    * `fname`: name pattern of the file to printe the parameters
                of the state of the run
 Global tructures used: gd, cmd
 Counting encounters (in global gd): nbbcalc, nbccalc, ncccalc
 Return (the error status):
    int SUCCESS or FAILURE
 */
global int PrintState(struct  cmdline_data *cmd,
                      struct  global_data* gd,
                      char *fname)
{
// Every item you add to cmdline_data, global_data structures
//  or to globaldefs.h file should you consider necessary to
//  recover the run, must be added to this routine.
// Also all items should be printed
//  in the same order as was given in savestate routine

    FILE *fdout;
    char buf[200];
    int  errorFlag=0;
    int ifile;
    int i;
    int n;
    int k;

    //B Look for "/" if fname is composed: path and filename
    int ndefault = 0;
    int ipos;
    char *dp;
    for (int i; i< strlen(fname); i++) {
        if(fname[i] == '/') {
            ipos = i+1;
            ndefault++;
        }
    }

    if (ndefault == 0) {
        sprintf(buf,"%s/%s%s",cmd->rootDir,fname,"-restorevalues");
    } else {
        dp = (char*) malloc((strlen(fname)-ipos)*sizeof(char));
        strncpy(dp, fname + ipos, strlen(fname)-ipos);
        verb_print_q(3,cmd->verbose,
                    "PrintParameterFile: '/' counts %d pos %d and %s\n",
                    ndefault, ipos, dp);
        sprintf(buf,"%s/%s%s",cmd->rootDir,dp,"-restorevalues");
    }
    //E

    if(!(fdout=fopen(buf,"w"))) {
        fprintf(stdout,"error opening file '%s' \n",buf);
        errorFlag=1;
    } else {
        //=============================================
        //B cmdline_data

        fprintf(fdout,"\ncmdline data:\n\n");

        //B Parameters related to the searching method
        fprintf(fdout,FMTT,"searchMethod",cmd->searchMethod);
        fprintf(fdout,FMTI,"mChebyshev",cmd->mChebyshev);
        fprintf(fdout,FMTT,"nsmooth",cmd->nsmooth);
        fprintf(fdout,FMTT,"rsmooth",cmd->rsmooth);
        fprintf(fdout,FMTR,"theta",cmd->theta);
        fprintf(fdout,FMTT,"computeTPCF",cmd->computeTPCF ? "true" : "false");
        fprintf(fdout,FMTT,"computeShearCF",
                cmd->computeShearCF ? "true" : "false");
        fprintf(fdout,FMTT,"usePeriodic",cmd->usePeriodic ? "true" : "false");
        //E

        //B Parameters to control the I/O file(s)
        // Input catalog parameters
        fprintf(fdout,FMTT,"infile",cmd->infile);
        fprintf(fdout,FMTT,"infileformat",cmd->infilefmt);
        fprintf(fdout,FMTT,"iCatalogs",cmd->iCatalogs);
        // Output parameters
        fprintf(fdout,FMTT,"rootDir",cmd->rootDir);
        fprintf(fdout,FMTT,"outfile",cmd->outfile);
        fprintf(fdout,FMTT,"outfileformat",cmd->outfilefmt);
        // Parameters to set a region in the sky,
        //  for example for Takahasi data set
        fprintf(fdout,FMTR,"thetaL",cmd->thetaL);
        fprintf(fdout,FMTR,"thetaR",cmd->thetaR);
        fprintf(fdout,FMTR,"phiL",cmd->phiL);
        fprintf(fdout,FMTR,"phiR",cmd->phiR);
        //E

        //B Parameters to control histograms and their output files
        fprintf(fdout,FMTT,"useLogHist",cmd->useLogHist ? "true" : "false");
        fprintf(fdout,FMTI,"logHistBinsPD",cmd->logHistBinsPD);
        //
        fprintf(fdout,FMTI,"sizeHistN",cmd->sizeHistN);
        fprintf(fdout,FMTR,"rangeN",cmd->rangeN);
        fprintf(fdout,FMTR,"rminHist",cmd->rminHist);
        fprintf(fdout,FMTI,"sizeHistPhi",cmd->sizeHistPhi);
        //
        fprintf(fdout,FMTT,"histNNFileName",cmd->histNNFileName);
        fprintf(fdout,FMTT,"histXi2pcfFileName",cmd->histXi2pcfFileName);
        fprintf(fdout,FMTT,"histZetaFileName",cmd->histZetaFileName);
        fprintf(fdout,FMTT,"suffixOutFiles",cmd->suffixOutFiles);
        //E

        //B Set of parameters needed to construct a test model
        fprintf(fdout,FMTI,"seed",cmd->seed);
        fprintf(fdout,FMTT,"testmodel",cmd->testmodel);
#ifdef LONGINT
        fprintf(fdout,FMTIL,"nbody",cmd->nbody);
#else
        fprintf(fdout,FMTI,"nbody",cmd->nbody);
#endif
        fprintf(fdout,FMTR,"lengthBox",cmd->lengthBox);
        //E

        //B Miscellaneous parameters
        fprintf(fdout,FMTTS,"script",cmd->script);
#ifdef LONGINT
        fprintf(fdout,FMTIL,"stepState",cmd->stepState);
#else
        fprintf(fdout,FMTI,"stepState",cmd->stepState);
#endif
        fprintf(fdout,FMTI,"verbose",cmd->verbose);
        fprintf(fdout,FMTI,"verbose_log",cmd->verbose_log);
#ifdef OPENMPCODE
        fprintf(fdout,FMTI,"numberThreads",cmd->numthreads);
#endif
        if (cmd->verbose>=VERBOSEDEBUGINFO)
            verb_print(cmd->verbose,
                       "\nPrintParamterFile: script: %s\n", cmd->script);

        fprintf(fdout,FMTT,"options",cmd->options);
        //E

//B this lines as in startrun, do the job...
#ifdef ADDONS
#include "startrun_include_08.h"
#endif
//E

        //E cmdline_data
        //=============================================

        //=============================================
        //B global_data

        fprintf(fdout,"\n\nglobal data:\n\n");

        //B global_data
        fprintf(fdout,FMTR,"cpuinit",gd->cpuinit);
        fprintf(fdout,FMTR,"cpurealinit",gd->cpurealinit);

        fprintf(fdout,FMTT,"headline0",gd->headline0);
        fprintf(fdout,FMTT,"headline1",gd->headline1);
        fprintf(fdout,FMTT,"headline2",gd->headline2);
        fprintf(fdout,FMTT,"headline3",gd->headline3);
        fprintf(fdout,FMTT,"mode",gd->mode);

        //B Tree
        fprintf(fdout,"\n\n\ttree data (begin):\n\n");

        fprintf(fdout,FMTIL,"ncell",gd->ncell);
        fprintf(fdout,FMTI,"tdepth",gd->tdepth);
        fprintf(fdout,FMTIL,"actmax",gd->actmax);
        fprintf(fdout,FMTIL,"sameposcount",gd->sameposcount);
        // BALLS
        fprintf(fdout,FMTIL,"ncccalc",gd->ncccalc);
        fprintf(fdout,FMTIL,"nsmoothcount",gd->nsmoothcount);
        //
        fprintf(fdout,FMTIL,"nbccalc",gd->nbccalc);
        fprintf(fdout,FMTIL,"nbbcalc",gd->nbbcalc);
        fprintf(fdout,FMTR,"rSize",gd->rSize);

        fprintf(fdout,FMTI,"ninfiles",gd->ninfiles);
        for (ifile=0; ifile<gd->ninfiles; ifile++) {
            fprintf(fdout,FMTR,"rSizeTable",gd->rSizeTable[ifile]);
        }
        for (ifile=0; ifile<gd->ninfiles; ifile++) {
            fprintf(fdout,FMTIL,"ncellTable",gd->ncellTable[ifile]);
        }
        for (ifile=0; ifile<gd->ninfiles; ifile++) {
            fprintf(fdout,FMTIL,"nbodyTable",gd->nbodyTable[ifile]);
        }
        for (ifile=0; ifile<gd->ninfiles; ifile++) {
            fprintf(fdout,FMTI,"tdepthTable",gd->tdepthTable[ifile]);
        }
        for (ifile=0; ifile<gd->ninfiles; ifile++) {
            fprintf(fdout,FMTIL,"nnodescanlevTable",
                    gd->nnodescanlevTable[ifile]);
        }
        for (ifile=0; ifile<gd->ninfiles; ifile++) {
            fprintf(fdout,FMTIL,"nnodescanlev_rootTable",
                    gd->nnodescanlev_rootTable[ifile]);
        }

        fprintf(fdout,FMTR,"cputree",gd->cputree);

        fprintf(fdout,"\n\n\ttree data (end):\n\n");
        //E

        fprintf(fdout,FMTR,"Rcut",gd->Rcut);
        fprintf(fdout,FMTR,"RcutSq",gd->RcutSq);
        fprintf(fdout,FMTR,"cpusearch",gd->cpusearch);

#ifdef CELLMETHOD
// Cell search
//     vectorI cells;
//    INTEGER *cellList;
//
#endif

        fprintf(fdout,FMTI,"infilefmt_int",gd->infilefmt_int);

        int ndim = NDIM;
        fprintf(fdout,FMTI,"NDIM",ndim);

//B here NDIM is used
#ifdef SINGLEP
        DO_COORD(k) {
            fprintf(fdout,FMTR,"Box",gd->Box[k]);
        }
#else
        DO_COORD(k) {
            fprintf(fdout,FMTR,"Box",gd->Box[k]);
        }
#endif
//E

        fprintf(fdout,FMTI,"searchMethod_int",gd->searchMethod_int);

        fprintf(fdout,FMTR,"deltaR",gd->deltaR);
        fprintf(fdout,FMTR,"deltaRmin",gd->deltaRmin);
        fprintf(fdout,FMTR,"deltaRmax",gd->deltaRmax);
        if (cmd->useLogHist) {
            for (n=1; n<=cmd->sizeHistN; n++) {
                fprintf(fdout,FMTR,"deltaRV",gd->deltaRV[n]);
            }
            for (n=1; n<=cmd->sizeHistN-1; n++) {
                fprintf(fdout,FMTR,"ddeltaRV",gd->ddeltaRV[n]);
            }
        }

        fprintf(fdout,FMTR,"deltaPhi",gd->deltaPhi);

        //B File pointers:
        fprintf(fdout,FMTT,"logfilePath",gd->logfilePath);
        fprintf(fdout,FMTT,"outputDir",gd->outputDir);
        fprintf(fdout,FMTT,"tmpDir",gd->tmpDir);
        fprintf(fdout,FMTT,"fpfnameOutputFileName",gd->fpfnameOutputFileName);
        fprintf(fdout,FMTT,"fpfnamehistNNFileName",gd->fpfnamehistNNFileName);
        fprintf(fdout,FMTT,"fpfnamehistCFFileName",gd->fpfnamehistCFFileName);
        fprintf(fdout,FMTT,"fpfnamehistrBinsFileName",
                gd->fpfnamehistrBinsFileName);
        fprintf(fdout,FMTT,"fpfnamehistXi2pcfFileName",
                gd->fpfnamehistXi2pcfFileName);
        fprintf(fdout,FMTT,"fpfnamehistZetaGFileName",
                gd->fpfnamehistZetaGFileName);
        fprintf(fdout,FMTT,"fpfnamehistZetaGmFileName",
                gd->fpfnamehistZetaGmFileName);
        fprintf(fdout,FMTT,"fpfnamehistZetaMFileName",
                gd->fpfnamehistZetaMFileName);
        fprintf(fdout,FMTT,"fpfnamemhistZetaMFileName",
                gd->fpfnamemhistZetaMFileName);
        fprintf(fdout,FMTT,"fpfnameCPUFileName",gd->fpfnameCPUFileName);
        //E

        fprintf(fdout,FMTT,"model_comment",gd->model_comment);
        fprintf(fdout,FMTT,"input_comment",gd->input_comment);
        fprintf(fdout,FMTT,"output_comment",gd->output_comment);

        //B save-restore
        fprintf(fdout,FMTI,"stopflag",gd->stopflag);
        fprintf(fdout,FMTIL,"ip",gd->ip);
        //E

        fprintf(fdout,FMTR,"cputotalinout",gd->cputotalinout);
        fprintf(fdout,FMTR,"cputotal",gd->cputotal);

        fprintf(fdout,FMTIL,"bytes_tot",gd->bytes_tot);

        fprintf(fdout,"\n\n\ttree data (begin):\n\n");

        fprintf(fdout,FMTIL,"bytes_tot_cells",gd->bytes_tot_cells);
        fprintf(fdout,FMTIL,"nbodySel",gd->nbodySel);
        fprintf(fdout,FMTIL,"nbodysm",gd->nbodysm);
        fprintf(fdout,FMTIL,"nbodybf",gd->nbodybf);

        fprintf(fdout,"\n\n\ttree data (end):\n\n");

        fprintf(fdout,FMTT,"bh86",gd->bh86 ? "true" : "false");
        fprintf(fdout,FMTT,"sw94",gd->sw94 ? "true" : "false");

        fprintf(fdout,FMTR,"i_deltaR",gd->i_deltaR);

#ifdef ADDONS
//#include "global_data_include.h"
#ifdef NMultipoles
//B Already safewrite when KDTREEOMP is defined
//#ifndef KDTREEOMP             // NMultipoles Deactivated in KDTREEOMP
//B histograms
//E histograms
    
    
#ifdef USEGSL
//    gsl_matrix_complex *NhistXi_gsl;
#endif

#endif // ! NMultipoles

        fprintf(fdout,FMTT,"fpfnamehistN2pcfFileName",
                gd->fpfnamehistN2pcfFileName);

#endif // ! ADDONS

        fprintf(fdout,"\n\n\tkd data (begin):\n\n");

        fprintf(fdout,FMTT,"fnameData_kd", gd->fnameData_kd);
        fprintf(fdout,FMTT,"fnameOut_kd", gd->fnameOut_kd);
        fprintf(fdout,FMTI,"input_format_kd",gd->input_format_kd);
        fprintf(fdout,FMTI,"use_tree_kd",gd->use_tree_kd);
        fprintf(fdout,FMTI,"max_tree_order_kd",gd->max_tree_order_kd);
        fprintf(fdout,FMTI,"max_tree_nparts_kd",gd->max_tree_nparts_kd);
        fprintf(fdout,FMTI,"use_pm_kd",gd->use_pm_kd);
        fprintf(fdout,FMTIL,"n_objects_kd",gd->n_objects_kd);
        fprintf(fdout,FMTR,"l_box_kd",gd->l_box_kd);
        fprintf(fdout,FMTR,"l_box_half_kd",gd->l_box_half_kd);

        fprintf(fdout,"\n\n\tkd data (end):\n\n");


        for (ifile=0; ifile<gd->ninfiles; ifile++) {
            fprintf(fdout,FMTT,"infilenames",gd->infilenames[ifile]);
        }
        for (ifile=0; ifile<gd->ninfiles; ifile++) {
            fprintf(fdout,FMTT,"infilefmtname",gd->infilefmtname[ifile]);
        }
        for (ifile=0; ifile<gd->ninfiles; ifile++) {
            fprintf(fdout,FMTI,"iCatalogs",gd->iCatalogs[ifile]);
        }

        for (i=0; i<1; i++) {
            fprintf(fdout,FMTI,"nsmooth",gd->nsmooth[i]);
        }


        fprintf(fdout,"\n\n\ttree data (begin):\n\n");

        fprintf(fdout,FMTIL,"nnode",gd->nnode);
        fprintf(fdout,FMTIL,"rnnode",gd->rnnode);
        for (i=0; i<2; i++) {
            fprintf(fdout,FMTI,"scanLevelMin",gd->scanLevelMin[i]);
        }
        fprintf(fdout,FMTT,"nodesfilePath",gd->nodesfilePath);
        fprintf(fdout,FMTI,"nnodescanlev",gd->nnodescanlev);
        //B Root nodes:
        fprintf(fdout,FMTI,"nnodescanlev_root",gd->nnodescanlev_root);

//B must see MAXLEV definition in global_data.h
#define MAXLEVEL  32
        for (i=0; i<MAXLEVEL; i++) {
            fprintf(fdout,FMTR,"Rcell",gd->Rcell[i]);
        }
#undef MAXLEVEL
//E

        fprintf(fdout,FMTT,"bodiesfilePath",gd->bodiesfilePath);

        fprintf(fdout,FMTT,"flagSmoothCellMin",
                gd->flagSmoothCellMin ? "true" : "false");
        fprintf(fdout,FMTT,"flagSmooth",gd->flagSmooth ? "true" : "false");
        fprintf(fdout,FMTT,"flagSetNbNoSel",
                gd->flagSetNbNoSel ? "true" : "false");

        //B BUCKET
        for (i=0; i<2; i++) {
            fprintf(fdout,FMTR,"rminCell",gd->rminCell[i]);
        }
        //E

        fprintf(fdout,"\n\n\ttree data (end):\n\n");



        for (i=0; i<1; i++) {
            fprintf(fdout,FMTR,"rsmooth",gd->rsmooth[i]);
        }

        fprintf(fdout,FMTT,"rsmoothFlag",gd->rsmoothFlag ? "true" : "false");
        fprintf(fdout,FMTI,"irsmooth",gd->irsmooth);

#ifdef ADDONS
//#include "globaldefs_include_03.h"
#ifdef BALLS
//#include "globaldefs_balls_omp_02.h"  // empty defs
#endif

#ifdef IOLIB
//#include "global_data_iolib.h"
    // pos, kappa, gamma1, gamma2, weight, seven places maximum
    // use in multi-columns-ascii and cfitsio
    // default is 3D: three pos and one convergence (kappa)
    for (i=0; i<6; i++) {
        fprintf(fdout,FMTI,"columns",gd->columns[i]);
    }
#endif
#endif // ! ADDONS

        //E global_data
        //=============================================

        //B histograms

        //E histograms


        //=============================================
        //B globaldefs.h

#ifdef USEGSL
//    gsl_rng * r;            // It is used r_gsl in globaldefs.h. Check!!!
#else
        fprintf(fdout,FMTIL,"idum",idum);
#endif

//B OpenMP section
//  Timing variables
#ifdef OPENMPCODE
        fprintf(fdout,FMTR,"relstart",relstart);
        fprintf(fdout,FMTR,"relend",relend);
        fprintf(fdout,FMTR,"absstart",absstart);
        fprintf(fdout,FMTR,"absend",absend);
#else
        fprintf(fdout,FMTR,"relstart",relstart);
        fprintf(fdout,FMTR,"relend",relend);
        fprintf(fdout,FMTR,"absstart",absstart);
        fprintf(fdout,FMTR,"absend",absend);
#endif
//E

        bodyptr p;
        for (ifile=0; ifile<gd->ninfiles; ifile++) {
            fprintf(fdout,FMTI,"bodytable ifile:",ifile);
            for (p = bodytable[ifile];
                 p < bodytable[ifile] + gd->nbodyTable[ifile]; p++) {
                // DUMMY so far...
            }
        }

#ifdef ADDONS
//#include "globaldefs_include_04.h"
//#include "cballsio_gadget_00.h"
//#include "globaldefs_balls_omp_04.h"
#endif

#ifdef ADDONS
//#include "globaldefs_include_05.h"
#endif

//E globaldefs.h
//=============================================

        fprintf(fdout,"\n\n");

    }
    fclose(fdout);

    if(errorFlag) {
        exit(0);
    }

    return SUCCESS;

}

#undef FMTT
#undef FMTTS
#undef FMTI
#undef FMTIL
#undef FMTR

//E
//#endif // SAVERESTORE

#endif	// ! _cballsio_save_restore_11_h
