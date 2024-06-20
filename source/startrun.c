/*==============================================================================
 MODULE: startrun.c				[cTreeBalls]
 Written by: Mario A. Rodriguez-Meza
 Starting date: april 2023
 Purpose: routines to initialize the main code
 Language: C
 Use: 'StartRun();'
 
 Mayor revisions:
 ==============================================================================*/
//        1          2          3          4          5          6          7

//
// We must check the order of memory allocation and dealocation!!!
// Here and in EndRun in tpfc_io.c
//

#include "globaldefs.h"

local void ReadParameterFile(struct  cmdline_data* cmd, struct  global_data* gd, char *);

local int startrun_parameterfile(struct  cmdline_data* cmd, struct  global_data* gd);
local int startrun_cmdline(struct  cmdline_data* cmd, struct  global_data* gd);
local void ReadParametersCmdline(struct  cmdline_data* cmd, struct  global_data* gd);
local void ReadParametersCmdline_short(struct  cmdline_data* cmd, struct  global_data* gd);
local void startrun_ParamStat(struct  cmdline_data* cmd, struct  global_data* gd);
local int CheckParameters(struct  cmdline_data* cmd, struct  global_data* gd);
local int random_init(struct  cmdline_data* cmd, struct  global_data* gd, int seed);
local void search_method_string_to_int(string method_str,int *method_int);
local int expandbox(struct  cmdline_data* cmd, struct  global_data* gd, bodyptr btab, int nbody, int);

local int infilefmt_string_to_int(string infmt_str,int *infmt_int);

local int startrun_getParamsSpecial(struct  cmdline_data* cmd, struct  global_data* gd);
local int scaniOption(struct  cmdline_data* cmd, struct  global_data* gd,
                      string, int *, int *, int, int, string);
local int scanrOption(struct  cmdline_data* cmd, struct  global_data* gd,
                      string, double *, int *, int, int, string);

#ifndef USEGSL
local long saveidum;
#endif

/*
 StartRun routine:

 To be called in main:
 StartRun(argv[0], HEAD1, HEAD2, HEAD3);
 
 This routine is in charge of setting all global structures in order
    the comutation process run smoothly with all parameters given
    by the user set and checked.

 Arguments:
    * `head0`: Input: string
    * `head1`: Input: string
    * `head2`: Input: string
    * `head3`: Input: string
 Return:
    The error status: int
 */
#ifndef CLASSLIB
int StartRun(struct  cmdline_data* cmd, struct  global_data* gd, string head0, string head1, string head2, string head3)
{
    double cpustart = CPUTIME;

    gd->headline0 = head0; gd->headline1 = head1;
    gd->headline2 = head2; gd->headline3 = head3;
#ifdef MPICODE
    if (ThisTask==0) {
#endif
        printf("\n%s\n%s: %s\n\t %s\n",
                   gd->headline0, gd->headline1, gd->headline2, gd->headline3);
        printf("Version = %s\n", getversion());
#ifdef MPICODE
    }
#endif

    gd->stopflag = 0;
    gd->cputotalinout = 0.;
    gd->cputotal = 0.;
    gd->bytes_tot = 0;

#ifdef GETPARAM
    cmd->paramfile = GetParam("paramfile");
    if (*(cmd->paramfile)=='-')
        error("bad parameter %s\n", cmd->paramfile);
    if (!strnull(cmd->paramfile))
		startrun_parameterfile(cmd, gd);
    else
		startrun_cmdline(cmd, gd);
#else
    startrun_parameterfile(cmd, gd);
#endif

    class_call_cballs(StartOutput(cmd), errmsg, errmsg);

    gd->bytes_tot += sizeof(struct  global_data);
    gd->bytes_tot += sizeof(struct cmdline_data);
    verb_print(cmd->verbose,
               "\nStartRun: Total allocated %g MByte storage so far.",
               gd->bytes_tot*INMB);

#ifdef OPENMPCODE
    class_call_cballs(set_number_threads(cmd), errmsg, errmsg);
#endif

    gd->cputotalinout += CPUTIME - cpustart;
    verb_print(cmd->verbose, "\nStartRun: elapsed time: %g\n\n",
               CPUTIME - cpustart);

    return _SUCCESS_;
}

#else // ! CLASSLIB

#include "input.h"

int StartRun(struct  cmdline_data* cmd, struct  global_data* gd,
             string head0, string head1, string head2, string head3)
{
    struct file_content fc;

    double cpustart = CPUTIME;
    
    gd->headline0 = head0; gd->headline1 = head1;
    gd->headline2 = head2; gd->headline3 = head3;
#ifdef MPICODE
    if (ThisTask==0) {
#endif
        printf("\n%s\n%s: %s\n\t %s\n",
               gd->headline0, gd->headline1, gd->headline2, gd->headline3);
        printf("Version = %s\n", getversion());
#ifdef MPICODE
    }
#endif
    
    gd->stopflag = 0;
    gd->cputotalinout = 0.;
    gd->cputotal = 0.;
    gd->bytes_tot = 0;
    
#ifdef GETPARAM
    cmd->paramfile = GetParam("paramfile");
    if (*(cmd->paramfile)=='-')
        error("bad parameter %s\n", cmd->paramfile);
    if (!strnull(cmd->paramfile)) {
        class_call_cballs(input_find_file(cmd->paramfile, &fc, errmsg), errmsg, errmsg);
        class_call_cballs(input_read_from_file(cmd, &fc, errmsg), errmsg, errmsg);
        class_call_cballs(parser_free(&fc), errmsg, errmsg);
    } else
        startrun_cmdline(cmd, gd);
#else
    class_call_cballs(input_find_file(cmd->ParameterFile, &fc, errmsg), errmsg, errmsg);
    class_call_cballs(input_read_from_file(cmd, &fc, errmsg), errmsg, errmsg);
    class_call_cballs(parser_free(&fc), errmsg, errmsg);
#endif

    class_call_cballs(startrun_Common(cmd, gd), errmsg, errmsg);

#ifdef GETPARAM
    if (!strnull(cmd->paramfile))
        PrintParameterFile(cmd, cmd->paramfile);
#else
    PrintParameterFile(cmd, cmd->ParameterFile);
#endif

    class_call_cballs(StartOutput(cmd), errmsg, errmsg);

    gd->bytes_tot += sizeof(struct  global_data);
    gd->bytes_tot += sizeof(struct cmdline_data);
    verb_print(cmd->verbose,
               "\nStartRun: Total allocated %g MByte storage so far.",
               gd->bytes_tot*INMB);

#ifdef OPENMPCODE
    class_call_cballs(set_number_threads(cmd), errmsg, errmsg);
#endif

    gd->cputotalinout += CPUTIME - cpustart;
    verb_print(cmd->verbose, "\nStartRun: elapsed time: %g\n\n",
               CPUTIME - cpustart);

    return _SUCCESS_;
}
#endif // ! CLASSLIB


local int startrun_parameterfile(struct  cmdline_data* cmd, struct  global_data* gd)
{
#ifdef GETPARAM
	ReadParameterFile(cmd, gd, cmd->paramfile);
    ReadParametersCmdline_short(cmd, gd);

#ifdef ADDONS
#include "startrun_include_01.h"
#endif

#else // ! GETPARAM
    ReadParameterFile(cmd, gd, cmd->ParameterFile);
#endif  // ! GETPARAM

	startrun_Common(cmd, gd);
#ifdef GETPARAM
	PrintParameterFile(cmd, cmd->paramfile);
#else
    PrintParameterFile(cmd, cmd->ParameterFile);
#endif

    return _SUCCESS_;
}


#ifdef GETPARAM
#define parameter_null	"parameters_null-cballs"

//B Section for reading parameters from the command line

local int startrun_cmdline(struct  cmdline_data* cmd, struct  global_data* gd)
{
	ReadParametersCmdline(cmd, gd);
	startrun_Common(cmd, gd);
	PrintParameterFile(cmd, parameter_null);

    return _SUCCESS_;
}

local void ReadParametersCmdline(struct  cmdline_data* cmd, struct  global_data* gd)
{
// Every item in cmdline_defs.h must have an item here::

#ifdef MPICODE
    if (ThisTask==0) {                              // Input all parameters on
                                                    //  proccess 0
#endif
    cmd->searchMethod = GetParam("searchMethod");
    cmd->theta = GetdParam("theta");
    cmd->mchebyshev = GetiParam("mChebyshev");
//B Parameters to set a region in the sky, for example for Takahasi data set.
    cmd->thetaL = GetdParam("thetaL");
    cmd->thetaR = GetdParam("thetaR");
    cmd->phiL = GetdParam("phiL");
    cmd->phiR = GetdParam("phiR");
//E
    cmd->sizeHistN = GetiParam("sizeHistN");
    cmd->rangeN = GetdParam("rangeN");
    cmd->rminHist = GetdParam("rminHist");
    cmd->infile = GetParam("infile");
    cmd->infilefmt = GetParam("infileformat");
    cmd->iCatalogs = GetParam("iCatalogs");
    cmd->rootDir = GetParam("rootDir");
    cmd->outfile = GetParam("outfile");
    cmd->outfilefmt = GetParam("outfileformat");

    cmd->histNFileName = GetParam("histNFileName");
    cmd->histXi2pcfFileName = GetParam("histXi2pcfFileName");
    cmd->histZetaMFileName = GetParam("histZetaMFileName");
    cmd->mhistZetaFileName = GetParam("mhistZetaFileName");
    cmd->suffixOutFiles = GetParam("suffixOutFiles");

    cmd->stepState = GetlParam("stepState");

    cmd->verbose = GetiParam("verbose");
    cmd->verbose_log = GetiParam("verbose_log");

#ifdef OPENMPCODE
    cmd->numthreads = GetiParam("numberThreads");
#endif

    cmd->script = GetParam("script");
    cmd->options = GetParam("options");

//
//B NOLSST:
//
    cmd->seed=GetiParam("seed");    // to always have defaults // Check in gsl
    cmd->nsmooth = GetParam("nsmooth");
    cmd->testmodel = GetParam("testmodel");
    cmd->nbody = GetiParam("nbody");
    cmd->lengthBox = GetdParam("lengthBox");
//
//E
//

    cmd->rsmooth = GetParam("rsmooth");

        cmd->computeTPCF = GetbParam("computeTPCF");
        cmd->useLogHist = GetbParam("useLogHist");
        cmd->logHistBinsPD = GetiParam("logHistBinsPD");
        cmd->usePeriodic = GetbParam("usePeriodic");

#ifdef ADDONS
#include "startrun_include_02.h"
#endif

#ifdef MPICODE
    }
    MPI_Bcast(&cmd, sizeof(cmdline_data), MPI_BYTE, 0, MPI_COMM_WORLD);
#endif


}

local void ReadParametersCmdline_short(struct  cmdline_data* cmd, struct  global_data* gd)
{
#ifdef MPICODE
    if (ThisTask==0) {                              // Input all parameters on
                                                    //  proccess 0
#endif
//B Here add parameters needed to be read after reading parameter file
//    cmd->script = GetParam("script");
//E
#ifdef MPICODE
    }
    MPI_Bcast(cmd, sizeof(cmdline_data), MPI_BYTE, 0, MPI_COMM_WORLD);
#endif
}

//E

#undef parameter_null
#endif // end of GETPARAM

//local int startrun_Common(void)
 int startrun_Common(struct  cmdline_data* cmd, struct  global_data* gd)
{
    int ifile;

#ifdef ADDONS
#include "startrun_include_03.h"
#endif

     setFilesDirs(cmd, gd);
    setFilesDirs_log(cmd, gd);
    strcpy(gd->mode,"w");
#ifdef MPICODE
        if(ThisTask==0) {
#endif
            if(!(gd->outlog=fopen(gd->logfilePath, gd->mode)))
                error("\nstart_Common: error opening file '%s' \n",
                      gd->logfilePath);
#ifdef MPICODE
        }
#endif

     class_call_cballs(startrun_getParamsSpecial(cmd, gd), errmsg, errmsg);
     class_call_cballs(random_init(cmd, gd, cmd->seed), errmsg, errmsg);
     class_call_cballs(CheckParameters(cmd, gd), errmsg, errmsg);
     class_call_cballs(startrun_memoryAllocation(cmd, gd), errmsg, errmsg);

//B Pre-processing necessary for reading data files:
        double cpustart;
        char buf[200];
        if (scanopt(cmd->options, "pre-processing")) {
            cpustart = CPUTIME;
            sprintf(buf,"%s",cmd->script);
            verb_print(cmd->verbose,
                       "\npre-processing: executing %s...",cmd->script);
            system(buf);
            verb_print(cmd->verbose, " done.\n");
            gd->cputotalinout += CPUTIME - cpustart;
            verb_print(cmd->verbose, "cpu time expended in this script %g\n\n",
                       CPUTIME - cpustart);
            if (scanopt(cmd->options, "stop")) {
                verb_print(cmd->verbose, "\n\tMainLoop: stopping...\n\n");
                exit(1);
            }
        }
//E

//B In this section update computation of rSize and center-of-mass if necessary
//      so we have a common root size and c-of-m
        for (ifile=0; ifile<gd->ninfiles; ifile++) {
            if (!strnull(cmd->infile)) {
                class_call_cballs(infilefmt_string_to_int(gd->infilefmtname[ifile],
                                &gd->infilefmt_int), errmsg, errmsg);
                class_call_cballs(inputdata(cmd, gd, gd->infilenames[ifile],ifile),
                                  errmsg, errmsg);
            } else {
                testdata(cmd, gd);
            }
        }
//E

#ifdef MPICODE
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Bcast(&gd, sizeof(struct  global_data), MPI_BYTE, 0, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);

        if (ThisTask!=0)
            bodytab = (bodyptr) allocate(cmd->nbody * sizeof(body));

        MPI_Bcast(bodytab, cmd->nbody*sizeof(body), MPI_BYTE, 0, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
#endif

        for (ifile=0; ifile<gd->ninfiles; ifile++) {
            gd->bytes_tot += gd->nbodyTable[ifile]*sizeof(body);
            verb_print(cmd->verbose,
                "\n\nAllocated %g MByte for particle storage (file %d).\n",
                       gd->nbodyTable[ifile]*sizeof(body)*INMB, ifile);
        }

        search_method_string_to_int(cmd->searchMethod, &gd->searchMethod_int);
        for (ifile=0; ifile<gd->ninfiles; ifile++) {
            gd->rSizeTable[ifile] = 1.0;
            expandbox(cmd, gd, bodytable[ifile], gd->nbodyTable[ifile], ifile);
            if (cmd->rangeN > gd->rSizeTable[ifile])
                verb_print(cmd->verbose,
        "\nstartrun_Common: warning! rangeN (%g) is greather than rSize (%g) of the system...\n",
                cmd->rangeN, gd->rSizeTable[ifile]);
        }

//B Tree search:
        gd->Rcut = cmd->rangeN;                       // Maximum search radius
        gd->RcutSq = gd->Rcut*gd->Rcut;
//E
     if (cmd->useLogHist) {
         real rBin, rbinlog;
         if (cmd->rminHist==0) {
             gd->deltaRV = dvector(1,cmd->sizeHistN);
             verb_log_print(cmd->verbose_log, gd->outlog, "\ndeltaRV:\n");
             int n;
             for (n=1; n<=cmd->sizeHistN; n++) {
                 gd->deltaRV[n] =
                 cmd->rangeN*rpow( 10.0, ( (real)(n-cmd->sizeHistN) )/cmd->logHistBinsPD );
                 verb_log_print(cmd->verbose_log, gd->outlog,
                                " %d %lg\n",n,gd->deltaRV[n]);
             }
         } else {
             gd->deltaR = rlog10(cmd->rangeN/cmd->rminHist)/cmd->sizeHistN;
             gd->deltaRV = dvector(1,cmd->sizeHistN);
             gd->ddeltaRV = dvector(1,cmd->sizeHistN-1);
             verb_log_print(cmd->verbose_log, gd->outlog,
                            "deltaRV (deltaR=%lf logscale):\n", gd->deltaR);
             int n;
             rbinlog = rlog10(cmd->rminHist) + ((real)(1))*gd->deltaR;
             rBin=rpow(10.0,rbinlog);
             gd->deltaRV[1] = rBin;
             verb_log_print(cmd->verbose_log, gd->outlog,
                            " %d %lg\n",n,gd->deltaRV[1]);
             gd->deltaRmin=rBin;
             gd->deltaRmax=rBin;
             for (n=2; n<=cmd->sizeHistN; n++) {
                 rbinlog = rlog10(cmd->rminHist) + ((real)(n))*gd->deltaR;
                 rBin=rpow(10.0,rbinlog);
                 gd->deltaRV[n] = rBin;
                 verb_log_print(cmd->verbose_log, gd->outlog,
                                " %d %lg %lg\n",
                                n,gd->deltaRV[n],gd->deltaRV[n]-gd->deltaRV[n-1]);
                 //B Not working MIN and MAX macros... Test these macros
                 gd->ddeltaRV[n-1] = gd->deltaRV[n]-gd->deltaRV[n-1];
                 if (gd->deltaRmax < gd->deltaRV[n]-gd->deltaRV[n-1])
                     gd->deltaRmax = gd->deltaRV[n]-gd->deltaRV[n-1];
                 if (gd->deltaRmin > gd->deltaRV[n]-gd->deltaRV[n-1])
                     gd->deltaRmin = gd->deltaRV[n]-gd->deltaRV[n-1];
                 //E
             }
             verb_log_print(cmd->verbose_log, gd->outlog,
                            "deltaRV min and max: %lg %lg\n",
                            gd->deltaRmin, gd->deltaRmax);
         }
     } else { // ! useLogHist
         gd->deltaR = (cmd->rangeN-cmd->rminHist)/cmd->sizeHistN;
         verb_log_print(cmd->verbose_log, gd->outlog,
                        "deltaR=%lf normal scale):\n",gd->deltaR);
     } // ! useLogHist

        MULVS(gd->cells, gd->Box, 1.0/gd->Rcut);       // Only needed for cellmethod
        AllocMem(gd->cellList, VProd (gd->cells)      // Only needed for cellmethod
                + cmd->nbody, INTEGER);
        gd->bytes_tot += (VProd(gd->cells)+cmd->nbody)*sizeof(INTEGER);
        verb_print(cmd->verbose,
                   "\n\nAllocated %g MByte for cells storage...\n",
                   (VProd(gd->cells)+cmd->nbody)*sizeof(INTEGER)*INMB);
//E
        real Vol = 1.0;
        int k;
        DO_COORD(k)
            Vol = Vol*gd->Box[k];
        
        gd->i_deltaR = 1.0/gd->deltaR;                // This is gd->i_r_max 
                                                    //  change...
        verb_log_print(cmd->verbose_log, gd->outlog,
                       "\nRcut, deltaR: %g %g\n",gd->Rcut,gd->deltaR);
#if NDIM == 3
// CHECK!!! gd->Box[1]
        verb_log_print(cmd->verbose_log, gd->outlog,
                       "lbox: %g %g %g\n",gd->Box[0],gd->Box[1],gd->Box[2]);
#else
        verb_log_print(cmd->verbose_log, gd->outlog,
                       "lbox: %g %g\n",gd->Box[0],gd->Box[1]);
#endif
        verb_log_print(cmd->verbose_log, gd->outlog, "Box volume = %e\n",Vol);
    verb_log_print(cmd->verbose_log, gd->outlog,
                   "(V/N)^(1/3): %g\n\n",rpow(Vol/cmd->nbody,1.0/3.0));
    verb_log_print(cmd->verbose_log, gd->outlog,
                   "Unit sphere (Takahasi): (S/N)^(1/2): %g\n\n",
                   rpow(2.0*TWOPI/gd->nbodyTable[gd->iCatalogs[0]],1.0/2.0));
    real rSizeTmp;
    int i, idepth=64;
    rSizeTmp = gd->rSizeTable[gd->iCatalogs[0]];
    for (i = 1; i <= idepth; i++) {
        rSizeTmp = rSizeTmp/2.0;
        verb_log_print(cmd->verbose_log, gd->outlog, "Cell size = %e\n",rSizeTmp);
        if (rSizeTmp < rpow(2.0*TWOPI/gd->nbodyTable[gd->iCatalogs[0]],1.0/2.0)) {
            verb_log_print(cmd->verbose_log, gd->outlog,
                           "Cell size threshold = %d\n",i);
            gd->irsmooth = i;
            verb_print(cmd->verbose,
            "startrun_Common: i threshold, cell size and rsmooth = %d %e %e\n",
                       gd->irsmooth, rSizeTmp, gd->rsmooth[0]);
            break;
        }
    }


#ifdef ADDONS
#include "startrun_include_04.h"
#endif

    return _SUCCESS_;
}

local int startrun_getParamsSpecial(struct  cmdline_data* cmd, struct  global_data* gd)
{
    char *pch;
    int nitems, ndummy=1;
    char inputnametmp[MAXLENGTHOFSTRSCMD];
    int i;

    if (strnull(cmd->infile)) {
        if (cmd->verbose_log>=3)
       verb_log_print(cmd->verbose_log, gd->outlog,
                       "\nstartrun_getParamsSpecial: no inputfile was given, making data ...\n");
        gd->ninfiles=1;                              // To test data...
    } else {
        strcpy(inputnametmp,cmd->infile);
        if (cmd->verbose_log>=3)
        verb_log_print(cmd->verbose_log, gd->outlog,
                       "\nSplitting string \"%s\" in tokens:\n",inputnametmp);
        gd->ninfiles=0;
        pch = strtok(inputnametmp," ,");
        while (pch != NULL) {
            gd->infilenames[gd->ninfiles] = (string) malloc(MAXLENGTHOFFILES);
            strcpy(gd->infilenames[gd->ninfiles],pch);
            ++gd->ninfiles;
            if (cmd->verbose_log>=3)
            verb_log_print(cmd->verbose_log, gd->outlog, "%s\n",gd->infilenames[gd->ninfiles-1]);
            pch = strtok (NULL, " ,");
        }
        if (cmd->verbose_log>=3)
        verb_log_print(cmd->verbose_log, gd->outlog,
                       "num. of files in infile %s =%d\n",cmd->infile,gd->ninfiles);
    }

    if (!strnull(cmd->infilefmt)) {
        strcpy(inputnametmp,cmd->infilefmt);
        if (cmd->verbose_log>=3)
        verb_log_print(cmd->verbose_log, gd->outlog,
                       "\nSplitting string \"%s\" in tokens:\n",inputnametmp);
        nitems=0;
        pch = strtok(inputnametmp," ,");
        while (pch != NULL) {
            gd->infilefmtname[nitems] = (string) malloc(30);
            strcpy(gd->infilefmtname[nitems],pch);
            ++nitems;
            if (cmd->verbose_log>=3)
            verb_log_print(cmd->verbose_log, gd->outlog, "%s\n",gd->infilefmtname[nitems-1]);
            pch = strtok (NULL, " ,");
        }
        if (cmd->verbose_log>=3)
        verb_log_print(cmd->verbose_log, gd->outlog,
                       "num. of items in infilefmt %s =%d\n",cmd->infilefmt,nitems);
        if (nitems != gd->ninfiles)
            error("\nstartrun_Common: nitems must be equal to number of files\n\n");
    }

    scaniOption(cmd, gd, cmd->nsmooth, gd->nsmooth, &nitems, 1, 1, "nsmooth");

#ifdef BALLS
    verb_log_print(cmd->verbose_log, gd->outlog, "\nSplitting string \"%s\" in tokens:\n",cmd->scanLevelMin);
    scaniOption(cmd, gd, cmd->scanLevelMin, gd->scanLevelMin, &nitems, ndummy,
                2, "scanLevelMin");

    if (nitems==1)
        gd->scanLevelMin[1]=gd->scanLevelMin[0]-1;
    if (strnull(cmd->scanLevelMin)) {
        gd->scanLevelMin[0]=0;
        if (cmd->verbose_log>=3)
            verb_log_print(cmd->verbose_log, gd->outlog, "option: %d\n",
                           gd->scanLevelMin[0]);
        gd->scanLevelMin[1]=gd->scanLevelMin[0]-1;
        if (cmd->verbose_log>=3)
            verb_log_print(cmd->verbose_log, gd->outlog, "option: %d\n",
                           gd->scanLevelMin[1]);
    }
#endif

    scanrOption(cmd, gd, cmd->rsmooth, gd->rsmooth, &nitems, ndummy,
                2, "rsmooth");
    gd->rsmoothFlag = TRUE;
    if (nitems!=1 && !strnull(cmd->rsmooth)) {
        if (cmd->verbose_log>=3)
            verb_log_print(cmd->verbose_log, gd->outlog,
                           "option: rsmooth=%s is not valid... going out\n",
                           cmd->rsmooth);
        gd->rsmoothFlag = FALSE;
    }
//#endif

    scaniOption(cmd, gd, cmd->iCatalogs, gd->iCatalogs,
                &nitems, gd->ninfiles, 0, "iCatalogs");
    if (gd->ninfiles==1) {
        (gd->iCatalogs[0]) = 0;
        if (cmd->verbose_log>=2)
            verb_log_print(cmd->verbose_log, gd->outlog,
                        "option: iCatalogs final values: %d\n",
                        gd->iCatalogs[0]);
        (gd->iCatalogs[1]) = 0;
        if (cmd->verbose_log>=2)
            verb_log_print(cmd->verbose_log, gd->outlog,
                        "option: iCatalogs final values: %d\n",
                        gd->iCatalogs[1]);
    } else {
        for (i=0; i<gd->ninfiles; i++) {
            (gd->iCatalogs[i])--;
            if (cmd->verbose_log>=2)
                verb_log_print(cmd->verbose_log, gd->outlog,
                        "option: iCatalogs final values: %d\n",
                        gd->iCatalogs[i]);
        }
    }


    return _SUCCESS_;
}

local int scaniOption(struct  cmdline_data* cmd, struct  global_data* gd,
                      string optionstr, int *option, int *noption,
    int nfiles, int flag, string message)
{
    char *pch;
    char *poptionstr[30],  optiontmp[100];
    int i;
//
// Decide what is better: DEBUG or if (cmd->verbose_log>=3)
//
#ifdef DEBUG
    if (cmd->verbose_log>=3)
    verb_log_print(cmd->verbose_log, gd->outlog, "\nProcessing '%s' option:\n", message);
#endif

    verb_log_print(cmd->verbose_log, gd->outlog, "\nSplitting string \"%s\" in tokens:\n",message);

    if (!strnull(optionstr)) {
        strcpy(optiontmp,optionstr);
        if (cmd->verbose_log>=3)
        verb_log_print(cmd->verbose_log, gd->outlog, "\nSplitting string \"%s\" in tokens:\n",optiontmp);
        *noption=0;
        pch = strtok(optiontmp," ,");
        while (pch != NULL) {
            poptionstr[*noption] = (string) malloc(10);
            strcpy(poptionstr[*noption],pch);
            ++(*noption);
            if (cmd->verbose_log>=3)
            verb_log_print(cmd->verbose_log, gd->outlog, "%s\n",poptionstr[*noption-1]);
            pch = strtok (NULL, " ,");
        }

        if (cmd->verbose_log>=3)
        verb_log_print(cmd->verbose_log, gd->outlog, "num. of tokens in option %s =%d\n", optionstr,*noption);

        if (flag == 0)
            if (*noption != nfiles)
                error("\nscanOption: noption = %d must be equal to number of infiles\n\n",*noption);
        if (*noption > MAXITEMS)
            error("\nscanOption: noption = %d must be less than the maximum num. of lines\n\n",*noption);

        for (i=0; i<*noption; i++) {
            option[i]=atoi(poptionstr[i]);
            if (cmd->verbose_log>=3)
            verb_log_print(cmd->verbose_log, gd->outlog, "option: %d\n",option[i]);
        }
        if (cmd->verbose_log>=3)
        verb_log_print(cmd->verbose_log, gd->outlog, "\nnoptions, nfiles: %d %d\n",*noption,nfiles);
        if (flag == 1) {
            if (*noption > nfiles)
                error("\nscanOption: noption = %d must be less or equal to number of files\n\n",*noption);
            else {
                for (i=*noption; i<nfiles; i++) {
                    option[i]=option[i-1]+1;
                    if (cmd->verbose_log>=3)
                    verb_log_print(cmd->verbose_log, gd->outlog, "option: %d\n",option[i]);
                }
            }
        }
    } else {
        for (i=0; i<nfiles; i++) {
            option[i]=1;
            if (cmd->verbose_log>=3)
            verb_log_print(cmd->verbose_log, gd->outlog, "option: %d\n",option[i]);
        }
    }

    return _SUCCESS_;
}

local int scanrOption(struct  cmdline_data* cmd, struct  global_data* gd,
                      string optionstr, double *option, int *noption,
    int nfiles, int flag, string message)
{
    char *pch;
    char *poptionstr[30],  optiontmp[100];
    int i;
//
// Decide what is better: DEBUG or if (cmd->verbose_log>=3)
//
#ifdef DEBUG
    if (cmd->verbose_log>=3)
    verb_log_print(cmd->verbose_log, gd->outlog, "\nProcessing '%s' option:\n", message);
#endif

    verb_log_print(cmd->verbose_log, gd->outlog, "\nSplitting string \"%s\" in tokens:\n",message);

    if (!strnull(optionstr)) {
        strcpy(optiontmp,optionstr);
        if (cmd->verbose_log>=3)
        verb_log_print(cmd->verbose_log, gd->outlog, "\nSplitting string \"%s\" in tokens:\n",optiontmp);
        *noption=0;
        pch = strtok(optiontmp," ,");
        while (pch != NULL) {
            poptionstr[*noption] = (string) malloc(10);
            strcpy(poptionstr[*noption],pch);
            ++(*noption);
            if (cmd->verbose_log>=3)
            verb_log_print(cmd->verbose_log, gd->outlog, "%s\n",poptionstr[*noption-1]);
            pch = strtok (NULL, " ,");
        }
        if (cmd->verbose_log>=3)
        verb_log_print(cmd->verbose_log, gd->outlog, "num. of tokens in option %s =%d\n", optionstr,*noption);

        if (flag == 0)
            if (*noption != nfiles)
                error("\nscanOption: noption = %d must be equal to number of files\n\n",*noption);
        if (*noption > MAXITEMS)
            error("\nscanOption: noption = %d must be less than the maximum num. of lines\n\n",*noption);

        for (i=0; i<*noption; i++) {
            option[i]=atof(poptionstr[i]);
            if (cmd->verbose_log>=3)
            verb_log_print(cmd->verbose_log, gd->outlog, "option: %g\n",option[i]);
        }

        if (cmd->verbose_log>=3)
        verb_log_print(cmd->verbose_log, gd->outlog, "\nnoptions, nfiles: %d %d\n",*noption,nfiles);
        if (flag == 1) {
            if (*noption > nfiles)
                error("\nscanOption: noption = %d must be less or equal to number of files\n\n",*noption);
            else {
                for (i=*noption; i<nfiles; i++) {
                    option[i]=option[i-1]+1;
                    if (cmd->verbose_log>=3)
                    verb_log_print(cmd->verbose_log, gd->outlog, "option: %g\n",option[i]);
                }
            }
        }
    } else {
        for (i=0; i<nfiles; i++) {
//B Be aware of this values
            option[i]=0;
//E
            if (cmd->verbose_log>=3)
            verb_log_print(cmd->verbose_log, gd->outlog, "option: %d\n",option[i]);
        }
    }

    return _SUCCESS_;
}


#ifdef GETPARAM
//B Section of parameter stat
local void startrun_ParamStat(struct  cmdline_data* cmd, struct  global_data* gd)
{
// Every item in cmdline_defs.h must have an item here::

#ifdef MPICODE
    if (ThisTask==0) {                                  // Input all parameters
                                                        //  on proccess 0
#endif
    if (GetParamStat("searchMethod") & ARGPARAM)
        cmd->searchMethod = GetParam("searchMethod");

    if (GetParamStat("theta") & ARGPARAM)
        cmd->theta = GetdParam("theta");
    if (GetParamStat("infile") & ARGPARAM)
        cmd->infile = GetParam("infile");
    if (GetParamStat("infilefmt") & ARGPARAM)
        cmd->infilefmt = GetParam("infileformat");
    if (GetParamStat("iCatalogs") & ARGPARAM)
        cmd->iCatalogs = GetParam("iCatalogs");
    if (GetParamStat("rootDir") & ARGPARAM)
        cmd->rootDir = GetParam("rootDir");
    if (GetParamStat("outfile") & ARGPARAM)
        cmd->outfile = GetParam("outfile");
    if (GetParamStat("outfilefmt") & ARGPARAM)
        cmd->outfilefmt = GetParam("outfileformat");

    if (GetParamStat("mChebyshev") & ARGPARAM)
        cmd->mchebyshev = GetiParam("mChebyshev");
//B Parameters to set a region in the sky, for example for Takahasi data set.
    if (GetParamStat("thetaL") & ARGPARAM)
        cmd->thetaL = GetdParam("thetaL");
    if (GetParamStat("thetaR") & ARGPARAM)
        cmd->thetaR = GetdParam("thetaR");
    if (GetParamStat("phiL") & ARGPARAM)
        cmd->phiL = GetdParam("phiL");
    if (GetParamStat("thetaR") & ARGPARAM)
        cmd->phiR = GetdParam("phiR");
//E
    if (GetParamStat("sizeHistN") & ARGPARAM)
        cmd->sizeHistN = GetiParam("sizeHistN");
    if (GetParamStat("rangeN") & ARGPARAM)
        cmd->rangeN = GetdParam("rangeN");
    if (GetParamStat("rminHist") & ARGPARAM)
        cmd->rminHist = GetdParam("rminHist");

    if (GetParamStat("histNFileName") & ARGPARAM)
        cmd->histNFileName = GetParam("histNFileName");
    if (GetParamStat("histXi2pcfFileName") & ARGPARAM)
        cmd->histXi2pcfFileName = GetParam("histXi2pcfFileName");
    if (GetParamStat("histZetaMFileName") & ARGPARAM)
        cmd->histZetaMFileName = GetParam("histZetaMFileName");
    if (GetParamStat("mhistZetaFileName") & ARGPARAM)
        cmd->mhistZetaFileName = GetParam("mhistZetaFileName");
    if (GetParamStat("suffixOutFiles") & ARGPARAM)
        cmd->suffixOutFiles = GetParam("suffixOutFiles");

    if (GetParamStat("verbose") & ARGPARAM)
        cmd->verbose = GetiParam("verbose");
    if (GetParamStat("verbose_log") & ARGPARAM)
        cmd->verbose_log = GetiParam("verbose_log");

#ifdef OPENMPCODE
    if (GetParamStat("numberThreads") & ARGPARAM)
        cmd->numthreads = GetiParam("numberThreads");
#endif

    if (GetParamStat("script") & ARGPARAM)
        cmd->script = GetParam("script");
    if (GetParamStat("options") & ARGPARAM)
        cmd->options = GetParam("options");

//
    if (GetParamStat("seed") & ARGPARAM)
        cmd->seed = GetiParam("seed");
    if (GetParamStat("nsmooth") & ARGPARAM)
        cmd->nsmooth = GetParam("nsmooth");
    if (GetParamStat("testmodel") & ARGPARAM)
        cmd->testmodel = GetParam("testmodel");
    if (GetParamStat("nbody") & ARGPARAM)
        cmd->nbody = GetiParam("nbody");
    if (GetParamStat("lengthBox") & ARGPARAM)
        cmd->lengthBox = GetdParam("lengthBox");
//

    if (GetParamStat("rsmooth") & ARGPARAM)
        cmd->rsmooth = GetParam("rsmooth");

        if (GetParamStat("computeTPCF") & ARGPARAM)
            cmd->computeTPCF = GetbParam("computeTPCF");
        if (GetParamStat("useLogHist") & ARGPARAM)
            cmd->useLogHist = GetbParam("useLogHist");
        if (GetParamStat("logHistBinsPD") & ARGPARAM)
            cmd->logHistBinsPD = GetiParam("logHistBinsPD");
        if (GetParamStat("usePeriodic") & ARGPARAM)
            cmd->usePeriodic = GetbParam("usePeriodic");

#ifdef ADDONS
#include "startrun_include_05.h"
#endif

#ifdef MPICODE
    }
    MPI_Bcast(&cmd, sizeof(cmdline_data), MPI_BYTE, 0, MPI_COMM_WORLD);
#endif

}
//E
#endif


//B Section of parameter check
local int CheckParameters(struct  cmdline_data* cmd, struct  global_data* gd)
{
// If it is necessary: an item in cmdline_defs.h must have an item here::


//    if (!strnull(cmd->restorefile) && !strnull(cmd->infile))
//        fprintf(stdout,"\nCheckParameters: Warning! : %s\n\n",
//            "You are using options restorefile and infile at the same time");

    if (cmd->theta < 0)
        error("CheckParameters: absurd value for theta\n");
    if (cmd->mchebyshev < 1)
        error("CheckParameters: absurd value for mchebyshev\n");
//B Parameters to set a region in the sky, for example for Takahasi data set.
    if (cmd->thetaL < 0 || cmd->thetaL > PI)
        error("CheckParameters: absurd value for thetaL (must be in the range 0--pi)\n");
    if (cmd->thetaR < 0 || cmd->thetaR > PI)
        error("CheckParameters: absurd value for thetaR (must be in the range 0--pi)\n");
    if (cmd->phiL < 0 || cmd->phiL > PI)
        error("CheckParameters: absurd value for phiL (must be in the range 0--2pi)\n");
    if (cmd->phiR < 0 || cmd->phiR > PI)
        error("CheckParameters: absurd value for phiR (must be in the range 0--2pi)\n");
    if (cmd->thetaL > cmd->thetaR)
        error("CheckParameters: absurd value for thetaL (must be greater thatn thetaR)\n");
    if (cmd->phiL > cmd->phiR)
        error("CheckParameters: absurd value for phiL (must be greater thatn phiR)\n");
//E
    if (cmd->sizeHistN < 2)
        error("CheckParameters: absurd value for sizeHistN\n");
    if (cmd->rangeN < 0)
        error("CheckParameters: absurd value for rangeN\n");
    if (cmd->rminHist < 0 || cmd->rminHist > cmd->rangeN)
        error("CheckParameters: absurd value for rminHist\n");

#ifdef OPENMPCODE
    if (cmd->numthreads <= 0)
        error("CheckParameters: absurd value for numberThreads must be an integer >= 0\n");
#endif

//
    if (gd->nsmooth[0] < 1)
        error("CheckParameters: absurd value for nsmooth\n");
    if (cmd->nbody < 3)
        error("CheckParameters: absurd value for nbody\n");
    if (cmd->lengthBox <= 0)
        error("CheckParameters: absurd value for lengthBox\n");

    gd->bh86 = scanopt(cmd->options, "bh86");             // Barnes, J. & Hut, P.
                                                        //  1986. Nature 324,
                                                        //  446.
    gd->sw94 = scanopt(cmd->options, "sw94");             // Salmon, J.K. &
                                                        //  Warren, M.S. 1994.
                                                        //  J. Comp. Phys. 111,
                                                        //  136
    if (gd->bh86 && gd->sw94)
        error("CheckParameters: incompatible options bh86 and sw94\n");
//

    if (gd->rsmooth[0] < 0 || gd->rsmoothFlag==FALSE)
        error("CheckParameters: absurd value for rsmooth (%s)\n",cmd->rsmooth);

    if (cmd->useLogHist==FALSE && (strcmp(cmd->searchMethod,"balls-omp-sincos") == 0))
        error("CheckParameters: can´t have loghist and balls-omp-sincos (%d %s)\n",
              cmd->useLogHist, cmd->searchMethod);

#ifdef ADDONS
#include "startrun_include_06.h"
#endif

    return _SUCCESS_;
}
//E

//B Section of parameter reading from a file
local void ReadParameterFile(struct  cmdline_data* cmd, struct  global_data* gd, char *fname)
{
// Every item in cmdline_defs.h must have an item here::
#define DOUBLE 1
#define STRING 2
#define INT 3
#define LONG 6
#define BOOLEAN 4
#define MAXTAGS 300
#define MAXCHARBUF 300

    FILE *fd;

  char buf[MAXCHARBUF],buf1[MAXCHARBUF],buf2[MAXCHARBUF],buf3[MAXCHARBUF];
  int  i,j,nt;
  int  id[MAXTAGS];
  void *addr[MAXTAGS];
  char tag[MAXTAGS][50];
  int  errorFlag=0;

#ifdef MPICODE
    if (ThisTask==0) {                                  // Input all parameters
                                                        //  on proccess 0
#endif
  nt=0;

    SPName(cmd->searchMethod,"searchMethod",MAXLENGTHOFSTRSCMD);
    RPName(cmd->theta,"theta");
    SPName(cmd->infile,"infile",MAXLENGTHOFSTRSCMD);
    SPName(cmd->infilefmt,"infileformat",MAXLENGTHOFSTRSCMD);
    SPName(cmd->iCatalogs,"iCatalogs",MAXLENGTHOFSTRSCMD);
    SPName(cmd->rootDir,"rootDir",MAXLENGTHOFSTRSCMD);
    SPName(cmd->outfile,"outfile",MAXLENGTHOFSTRSCMD);
    SPName(cmd->outfilefmt,"outfileformat",MAXLENGTHOFSTRSCMD);
    IPName(cmd->mchebyshev,"mChebyshev");
//B Parameters to set a region in the sky, for example for Takahasi data set.
    RPName(cmd->thetaL,"thetaL");
    RPName(cmd->thetaR,"thetaR");
    RPName(cmd->phiL,"phiL");
    RPName(cmd->phiR,"phiR");
//E
    IPName(cmd->sizeHistN,"sizeHistN");
    RPName(cmd->rangeN,"rangeN");
    RPName(cmd->rminHist,"rminHist");
    SPName(cmd->histNFileName,"histNFileName",MAXLENGTHOFSTRSCMD);
    SPName(cmd->histXi2pcfFileName,"histXi2pcfFileName",MAXLENGTHOFSTRSCMD);
    SPName(cmd->histZetaMFileName,"histZetaMFileName",MAXLENGTHOFSTRSCMD);
    SPName(cmd->mhistZetaFileName,"mhistZetaFileName",MAXLENGTHOFSTRSCMD);
    SPName(cmd->suffixOutFiles,"suffixOutFiles",MAXLENGTHOFSTRSCMD);
    LPName(cmd->stepState,"stepState");

    IPName(cmd->verbose,"verbose");
    IPName(cmd->verbose_log,"verbose_log");
#ifdef OPENMPCODE
    IPName(cmd->numthreads,"numberThreads");
#endif
    SPName(cmd->script,"script",MAXLENGTHOFSTRSCMD);
    SPName(cmd->options,"options",MAXLENGTHOFSTRSCMD);

//
    IPName(cmd->seed,"seed");                            // to always have
                                                        //  defaults
    SPName(cmd->nsmooth,"nsmooth",MAXLENGTHOFSTRSCMD);
    SPName(cmd->testmodel,"testmodel",MAXLENGTHOFSTRSCMD);
    IPName(cmd->nbody,"nbody");
    RPName(cmd->lengthBox,"lengthBox");
//

    SPName(cmd->rsmooth,"rsmooth",MAXLENGTHOFSTRSCMD);

    BPName(cmd->computeTPCF,"computeTPCF");
    BPName(cmd->useLogHist,"useLogHist");
    IPName(cmd->logHistBinsPD,"logHistBinsPD");
    BPName(cmd->usePeriodic,"usePeriodic");

#ifdef ADDONS
#include "startrun_include_07.h"
#endif

    size_t slen;
    char *script1;
    char *script2;
    char *script3;
    char *script4;
    char buf4[MAXCHARBUF];
    char buf5[MAXCHARBUF];

	if((fd=fopen(fname,"r"))) {
		while(!feof(fd)) {
			fgets(buf,MAXCHARBUF,fd);
            if(sscanf(buf,"%s%s%s%s",buf1,buf2,buf3,buf4)<1)
                continue;
            if(sscanf(buf,"%s%s%s%s",buf1,buf2,buf3,buf4)<3)
                *buf3='\0';
            if(buf1[0]=='%' || buf1[0]=='#')
                continue;
            if(buf2[0]!='=')
                continue;
            for(i=0,j=-1;i<nt;i++)
                if(strcmp(buf1,tag[i])==0) {
                    j=i;
					tag[i][0]=0;
					break;
				}
			if(j>=0) {
                switch(id[j]) {
                    case DOUBLE:
                        *((double*)addr[j])=atof(buf3);
                        break;
                    case STRING:
                        if (strcmp(buf1,"script") == 0){ // To remove both '"'
                            printf("ReadParamterFile: script: %s %s\n", buf3, buf4);
                            slen = strlen(buf3);
                            cmd->script = (char*) malloc((slen-2)*sizeof(char));
                            script1 = (char*) malloc(slen*sizeof(char));
                            memcpy(script1,buf3,slen);
                            script2 = strchr(script1, '"');
                            slen = strlen(script2);
                            script3 = (char*) malloc(slen*sizeof(char));
                            memcpy(script3,script2+1,slen);
                            printf("ReadParamterFile: script3: %s\n", script3);
                            slen = strlen(buf4);
                            script4 = (char*) malloc(slen*sizeof(char));
                            memcpy(script4,buf4,slen-1);
                            printf("ReadParamterFile: script4: %s\n", script4);
                            sprintf(buf5,"%s %s",script3,script4);
                            printf("ReadParamterFile: buf5: %s\n", buf5);
                            strcpy(cmd->script,buf5);
                            printf("ReadParamterFile: addr: %s\n", cmd->script);
                        } else {
                            strcpy(addr[j],buf3);
                        }
						break;
                    case INT:
                        *((int*)addr[j])=atoi(buf3);
                        break;
                    case LONG:
                        *((long*)addr[j])=atol(buf3);
                        break;
					case BOOLEAN:
						if (strchr("tTyY1", *buf3) != NULL) {
							*((bool*)addr[j])=TRUE;
                        } else 
                            if (strchr("fFnN0", *buf3) != NULL)  {
                                *((bool*)addr[j])=FALSE;
                            } else {
                                error("getbparam: %s=%s not bool\n",buf1,buf3);
                            }
						break;
                }
            } else {
                fprintf(stdout, "Error in file %s: Tag '%s' %s.\n",
					fname, buf1, "not allowed or multiple defined");
                errorFlag=1;
//                exit(0);
            }
        }
        fclose(fd);
    } else {
        fprintf(stdout,"Parameter file %s not found.\n", fname);
#ifdef MPICODE
        errorFlag=2;
        MPI_Bcast(&errorFlag, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Finalize();
        exit(0);
#else
        errorFlag=2;
        exit(0);
#endif
    }
  
    for(i=0;i<nt;i++) {
        if(*tag[i]) {
//            if (cmd->verbose>=2)
                if (cmd->verbose>2)
                fprintf(stdout,
                "Warning! I miss a value for tag '%s' in parameter file '%s'.\n",
                tag[i],fname);
            switch(id[i]) {
                case DOUBLE:
                    *((double*)addr[i])=GetdParam(tag[i]);
                    break;
                case STRING:
                    strcpy(addr[i],GetParam(tag[i]));
                    break;
                case INT:
                    *((int*)addr[i])=GetiParam(tag[i]);
                    break;
                case LONG:
                    *((long*)addr[i])=GetlParam(tag[i]);
                    break;
                case BOOLEAN:
                    *((bool*)addr[i])=GetbParam(tag[i]);
                    /*
                    if (strchr("tTyY1", *buf3) != NULL) {
                        *((bool*)addr[i])=TRUE;
                    } else
                        if (strchr("fFnN0", *buf3) != NULL)  {
                            *((bool*)addr[i])=FALSE;
                        } else {
                            error("getbparam: %s=%s not bool\n",buf1,buf3);
                        }
                    */
                    break;
            }
            errorFlag=1;
        }
    }

//        printf("ReadParamterFile: cmd->script: %s\n", cmd->script);
#ifdef MPICODE
//        MPI_Bcast(&errorFlag, 1, MPI_INT, 0, MPI_COMM_WORLD);
//        MPI_Finalize();
//        exit(0);
#else
//        if (errorFlag==3)
//            error("getbparam: %s=%s not bool\n",);
#endif

#ifdef MPICODE
    }
    MPI_Bcast(&errorFlag, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&cmd, sizeof(cmdline_data), MPI_BYTE, 0, MPI_COMM_WORLD);
#endif


#undef DOUBLE 
#undef STRING 
#undef INT 
#undef BOOLEAN
#undef MAXTAGS
#undef MAXCHARBUF
}
//E

#define FMTT    "%-35s = %s\n"
#define FMTTS    "%-35s = \"%s\"\n"
#define FMTI    "%-35s = %d\n"
#define FMTIL    "%-35s = %ld\n"
#define FMTR	"%-35s = %g\n"

//B Section of parameter writing to a file
int PrintParameterFile(struct  cmdline_data *cmd, char *fname)
{
// Every item in cmdline_defs.h must have an item here::
    FILE *fdout;
    char buf[200];
    int  errorFlag=0;

#ifdef MPICODE
    if (ThisTask==0) {                                  // Output only on 
                                                        //  proccess 0
#endif
    sprintf(buf,"%s/%s%s",cmd->rootDir,fname,"-usedvalues");
    if(!(fdout=fopen(buf,"w"))) {
        fprintf(stdout,"error opening file '%s' \n",buf);
        errorFlag=1;
    } else {
        fprintf(fdout,FMTT,"searchMethod",cmd->searchMethod);
        fprintf(fdout,FMTR,"theta",cmd->theta);
        fprintf(fdout,FMTT,"infile",cmd->infile);
        fprintf(fdout,FMTT,"infileformat",cmd->infilefmt);
        fprintf(fdout,FMTT,"iCatalogs",cmd->iCatalogs);
        fprintf(fdout,FMTT,"rootDir",cmd->rootDir);
        fprintf(fdout,FMTT,"outfile",cmd->outfile);
        fprintf(fdout,FMTT,"outfileformat",cmd->outfilefmt);
        fprintf(fdout,FMTI,"mChebyshev",cmd->mchebyshev);
//B Parameters to set a region in the sky, for example for Takahasi data set.
        fprintf(fdout,FMTR,"thetaL",cmd->thetaL);
        fprintf(fdout,FMTR,"thetaR",cmd->thetaR);
        fprintf(fdout,FMTR,"phiL",cmd->phiL);
        fprintf(fdout,FMTR,"phiR",cmd->phiR);
//E
        fprintf(fdout,FMTI,"sizeHistN",cmd->sizeHistN);
        fprintf(fdout,FMTR,"rangeN",cmd->rangeN);
        fprintf(fdout,FMTR,"rminHist",cmd->rminHist);
        fprintf(fdout,FMTT,"histNFileName",cmd->histNFileName);
        fprintf(fdout,FMTT,"histXi2pcfFileName",cmd->histXi2pcfFileName);
        fprintf(fdout,FMTT,"histZetaMFileName",cmd->histZetaMFileName);
        fprintf(fdout,FMTT,"mhistZetaFileName",cmd->mhistZetaFileName);
        fprintf(fdout,FMTT,"suffixOutFiles",cmd->suffixOutFiles);
        fprintf(fdout,FMTIL,"stepState",cmd->stepState);

        fprintf(fdout,FMTI,"verbose",cmd->verbose);
        fprintf(fdout,FMTI,"verbose_log",cmd->verbose_log);
#ifdef OPENMPCODE
        fprintf(fdout,FMTI,"numberThreads",cmd->numthreads);
#endif
        printf("PrintParamterFile: script: %s\n", cmd->script);
        fprintf(fdout,FMTTS,"script",cmd->script);
        fprintf(fdout,FMTT,"options",cmd->options);
//
        fprintf(fdout,FMTI,"seed",cmd->seed);
        fprintf(fdout,FMTT,"nsmooth",cmd->nsmooth);
        fprintf(fdout,FMTT,"testmodel",cmd->testmodel);
        fprintf(fdout,FMTIL,"nbody",cmd->nbody);
        fprintf(fdout,FMTR,"lengthBox",cmd->lengthBox);
//

        fprintf(fdout,FMTT,"rsmooth",cmd->rsmooth);

        fprintf(fdout,FMTT,"computeTPCF",cmd->computeTPCF ? "true" : "false");
        fprintf(fdout,FMTT,"useLogHist",cmd->useLogHist ? "true" : "false");
        fprintf(fdout,FMTI,"logHistBinsPD",cmd->logHistBinsPD);
        fprintf(fdout,FMTT,"usePeriodic",cmd->usePeriodic ? "true" : "false");

#ifdef ADDONS
#include "startrun_include_08.h"
#endif

        fprintf(fdout,"\n\n");
    }
    fclose(fdout);

#ifdef MPICODE
    }
    if(errorFlag) {
        MPI_Finalize();
        exit(0);
    }
#else
    if(errorFlag) {
        exit(0);
    }
#endif

    return _SUCCESS_;
}
//E

#undef FMTT
#undef FMTTS
#undef FMTI
#undef FMTR


local int infilefmt_string_to_int(string infmt_str,int *infmt_int)
{
    *infmt_int=-1;
    if (strcmp(infmt_str,"columns-ascii") == 0)         *infmt_int = INCOLUMNS;
    if (strcmp(infmt_str,"columns-ascii-pos") == 0)     *infmt_int = INCOLUMNSPOS;
    if (strnull(infmt_str))                             *infmt_int = INNULL;
    if (strcmp(infmt_str,"binary") == 0)                *infmt_int = INCOLUMNSBIN;
    if (strcmp(infmt_str,"takahasi") == 0)              *infmt_int = INTAKAHASI;
    if (strcmp(infmt_str,"columns-ascii-2d-to-3d") == 0)
                                                    *infmt_int = INCOLUMNS2DTO3D;

#ifdef ADDONS
#include "startrun_include_09.h"
#endif

    return _SUCCESS_;
}

local int random_init(struct  cmdline_data* cmd, struct  global_data* gd, int seed)
{
#ifdef USEGSL
    const gsl_rng_type * T;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    gd->r = gsl_rng_alloc (T);
    gsl_rng_set(gd->r, seed);
    r_gsl = gd->r;
    gd->bytes_tot += (1)*sizeof(gsl_rng_type);
#else
    saveidum=idum;
    xsrandom(idum);
#endif

    return _SUCCESS_;
}

global int startrun_memoryAllocation(struct  cmdline_data *cmd, struct  global_data* gd)
{
    // Free allocated memory in reverse order as were allocated
    INTEGER bytes_tot_local=0;
    
    gd->histN = dvector(1,cmd->sizeHistN);
    gd->histCF = dvector(1,cmd->sizeHistN);
    bytes_tot_local += cmd->sizeHistN*sizeof(real);
    gd->histNSub = dvector(1,cmd->sizeHistN);
    // 2pcf
    gd->histNSubXi2pcf = dvector(1,cmd->sizeHistN);
    //B kappa Avg Rmin
    gd->histNSubXi2pcftotal = dvector(1,cmd->sizeHistN);
    //E
    //
    bytes_tot_local += cmd->sizeHistN*sizeof(real);
    gd->histNNN = dvector(1,cmd->sizeHistN);
    bytes_tot_local += cmd->sizeHistN*sizeof(real);
#ifdef ADDONS
#include "startrun_include_10.h"
#endif
    gd->histXi2pcf = dvector(1,cmd->sizeHistN);
    bytes_tot_local += cmd->sizeHistN*sizeof(real);
    if (cmd->computeTPCF) {
    gd->histXi = dmatrix(1,cmd->mchebyshev+1,1,cmd->sizeHistN);
    bytes_tot_local += (cmd->mchebyshev+1)*cmd->sizeHistN*sizeof(real);
    gd->histXicos = dmatrix(1,cmd->mchebyshev+1,1,cmd->sizeHistN);
    bytes_tot_local += (cmd->mchebyshev+1)*cmd->sizeHistN*sizeof(real);
    gd->histXisin = dmatrix(1,cmd->mchebyshev+1,1,cmd->sizeHistN);
    bytes_tot_local += (cmd->mchebyshev+1)*cmd->sizeHistN*sizeof(real);
    gd->histZetaM = dmatrix3D(1,cmd->mchebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    bytes_tot_local += (cmd->mchebyshev+1)*cmd->sizeHistN*cmd->sizeHistN*sizeof(real);
    gd->histZetaMcos =
    dmatrix3D(1,cmd->mchebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    bytes_tot_local +=
    (cmd->mchebyshev+1)*cmd->sizeHistN*cmd->sizeHistN*sizeof(real);
    gd->histZetaMsin =
    dmatrix3D(1,cmd->mchebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    bytes_tot_local += (cmd->mchebyshev+1)*cmd->sizeHistN*cmd->sizeHistN*sizeof(real);
    gd->histZetaMsincos =
    dmatrix3D(1,cmd->mchebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    bytes_tot_local +=
    (cmd->mchebyshev+1)*cmd->sizeHistN*cmd->sizeHistN*sizeof(real);
    
#ifdef USEGSL
    gd->histXi_gsl = gsl_matrix_complex_calloc(cmd->mchebyshev+1,cmd->sizeHistN);
    bytes_tot_local += (cmd->mchebyshev+1)*cmd->sizeHistN*2*sizeof(real);
    
    histZetaMatrix = (mMatrix_ptr) allocate((cmd->mchebyshev+1) * sizeof(mMatrix));
    int m;
    for (m=0; m<=cmd->mchebyshev; m++){
        histZetaMatrix[m].histZetaM =
        gsl_matrix_complex_calloc(cmd->sizeHistN,cmd->sizeHistN);
    }
    bytes_tot_local += (cmd->mchebyshev+1)*cmd->sizeHistN*cmd->sizeHistN*2*sizeof(real);
#endif
    }

    gd->bytes_tot += bytes_tot_local;
    verb_print(cmd->verbose,
    "\n\nstartrun_memoryAllocation: Allocated %g MByte for histograms storage.\n",
    bytes_tot_local*INMB);

    return _SUCCESS_;
}

local void search_method_string_to_int(string method_str,int *method_int)
{
// Every search method must have an item here::
    *method_int=-1;
    if (strnull(method_str))
                *method_int = SEARCHNULL;
    if (strcmp(method_str,"tree-omp-sincos") == 0)
                *method_int = TREEOMPMETHODSINCOS;

#ifdef ADDONS
#include "startrun_include_11.h"
#endif
}

local int expandbox(struct  cmdline_data* cmd, struct  global_data* gd, bodyptr btab, int nbody, int ifile)
{
    real dmax, d;
    bodyptr p;
    int k;
    bodyptr root;
    
    root = (bodyptr) allocate(1 * sizeof(body));
    gd->bytes_tot += (1)*sizeof(body);

    CLRV(Pos(root));        // Assuming that (0,0,...) is the center of the Box

    dmax = 0.0;
    DO_BODY(p, btab, btab+nbody)
        DO_COORD(k) {
            d = rabs(Pos(p)[k] - Pos(root)[k]);
            if (d > dmax)
                dmax = d;
        }
    while (gd->rSizeTable[ifile] < 2 * dmax)
      gd->rSizeTable[ifile] = 2 * gd->rSizeTable[ifile];

#ifdef DEBUG
    if (cmd->verbose_log>=3)
    verb_log_print(cmd->verbose_log, gd->outlog,"\nexpandbox: rSize = %g\n", gd->rSizeTable[ifile]);
#endif

    return _SUCCESS_;
}

#ifdef OPENMPCODE
int set_number_threads(struct  cmdline_data *cmd)
{
    omp_set_num_threads(cmd->numthreads);

    return _SUCCESS_;
}
#endif
