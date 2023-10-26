/*==============================================================================
 MODULE: startrun.c				[cTreeBalls]
 Written by: Mario A. Rodriguez-Meza
 Starting date: april 2023
 Purpose: routines to initialize the main code
 Language: C
 Use: 'StartRun();'
 
 Mayor revisions:
 ==============================================================================*/

//
// We must check the order of memory allocation and dealocation!!!
// Here and in EndRun in tpfc_io.c
//

#include "globaldefs.h"

local void ReadParameterFile(char *);
local void PrintParameterFile(char *);

local int startrun_parameterfile(void);
local int startrun_cmdline(void);
local void ReadParametersCmdline(void);
local int startrun_Common(void);
local void startrun_ParamStat(void);
local void CheckParameters(void);
local int random_init(int seed);
local void search_method_string_to_int(string method_str,int *method_int);
local int expandbox(bodyptr btab, int nbody);

local int infilefmt_string_to_int(string infmt_str,int *infmt_int);

local int startrun_getParamsSpecial(void);
local int scaniOption(string, int *, int *, int, int, string);

#ifdef OPENMPCODE
int set_number_threads(void);
#endif

// To debug including MPICODE
//verb_print_debug(1, "\nAqui voy (0)\n");
//

int StartRun(string head0, string head1, string head2, string head3)
{
    double cpustart = CPUTIME;

    gd.headline0 = head0; gd.headline1 = head1;
    gd.headline2 = head2; gd.headline3 = head3;
#ifdef MPICODE
    if (ThisTask==0) {
#endif
        printf("\n%s\n%s: %s\n\t %s\n",
                   gd.headline0, gd.headline1, gd.headline2, gd.headline3);
        printf("Version = %s\n", getversion());
#ifdef MPICODE
    }
#endif

    gd.stopflag = 0;
    gd.cputotalinout = 0.;
    gd.cputotal = 0.;
    gd.bytes_tot = 0;

#ifdef GETPARAM
    cmd.paramfile = GetParam("paramfile");
    if (!strnull(cmd.paramfile))
		startrun_parameterfile();
    else
		startrun_cmdline();
#else
    startrun_parameterfile();
#endif

	StartOutput();

    gd.bytes_tot += sizeof(global_data);
    gd.bytes_tot += sizeof(cmdline_data);
    verb_print(cmd.verbose, "\nStartRun: Total allocated %g MByte storage so far.",
               gd.bytes_tot*INMB);

#ifdef OPENMPCODE
    set_number_threads();
#endif

    gd.cputotalinout += CPUTIME - cpustart;
//    verb_print(cmd.verbose, "\nStartRun: elapsed time: %g sec.\n\n",CPUTIME - cpustart);
    verb_print(cmd.verbose, "\nStartRun: elapsed time: %g\n\n",CPUTIME - cpustart);

    return _SUCCESS_;
}

local int startrun_parameterfile(void)
{
#ifdef GETPARAM
	ReadParameterFile(cmd.paramfile);

#ifdef SAVERESTORE
    if (strnull(cmd.restorefile))
        startrun_ParamStat();
#endif
#else
    ReadParameterFile(cmd.ParameterFile);
#endif

	startrun_Common();
#ifdef GETPARAM
	PrintParameterFile(cmd.paramfile);
#else
    PrintParameterFile(cmd.ParameterFile);
#endif

    return _SUCCESS_;
}


#ifdef GETPARAM
#define parameter_null	"parameters_null-cballs"

//B Section of read parameters from the command line

local int startrun_cmdline(void)
{
	ReadParametersCmdline();
	startrun_Common();
	PrintParameterFile(parameter_null);

    return _SUCCESS_;
}

local void ReadParametersCmdline(void)
{
// Every item in cmdline_defs.h must have an item here::

#ifdef MPICODE
    if (ThisTask==0) {                              // Input all parameters on proccess 0
#endif
    cmd.searchMethod = GetParam("searchMethod");
    cmd.theta = GetdParam("theta");
    cmd.mchebyshev = GetiParam("mChebyshev");
//B Parameters to set a region in the sky, for example for Takahasi data set.
        cmd.thetaL = GetdParam("thetaL");
        cmd.thetaR = GetdParam("thetaR");
        cmd.phiL = GetdParam("phiL");
        cmd.phiR = GetdParam("phiR");
//E
//    cmd.dimension = GetiParam("dimension");
    cmd.sizeHistN = GetiParam("sizeHistN");
    cmd.rangeN = GetdParam("rangeN");
    cmd.rminHist = GetdParam("rminHist");
//    cmd.logHist = GetbParam("logHist");
    cmd.sizeHistTheta = GetiParam("sizeHistTheta");
    cmd.infile = GetParam("infile");
    cmd.infilefmt = GetParam("infileformat");
    cmd.rootDir = GetParam("rootDir");
    cmd.outfile = GetParam("outfile");
    cmd.outfilefmt = GetParam("outfileformat");

    cmd.histNFileName = GetParam("histNFileName");
    cmd.histXi2pcfFileName = GetParam("histXi2pcfFileName");
    cmd.histZetaMFileName = GetParam("histZetaMFileName");
    cmd.mhistZetaFileName = GetParam("mhistZetaFileName");
    cmd.suffixOutFiles = GetParam("suffixOutFiles");

    cmd.stepState = GetiParam("stepState");
#ifdef SAVERESTORE
    cmd.statefile = GetParam("statefile");
    cmd.restorefile = GetParam("restorefile");
#endif
    cmd.verbose = GetiParam("verbose");
    cmd.verbose_log = GetiParam("verbose_log");

    cmd.numthreads = GetiParam("numberThreads");

    cmd.options = GetParam("options");

//
//B NOLSST:
//
    cmd.seed=GetiParam("seed");    // to always have defaults // Check in gsl
    cmd.nsmooth = GetParam("nsmooth");
#ifdef BALLS
    cmd.ntosave = GetiParam("ntosave");
        cmd.scanLevel = GetiParam("scanLevel");
#endif
    cmd.stepNodes = GetiParam("stepNodes");
    cmd.ncritical = GetParam("ncritical");
//#endif
    cmd.testmodel = GetParam("testmodel");
    cmd.nbody = GetiParam("nbody");
    cmd.lengthBox = GetdParam("lengthBox");
//    cmd.mToPlot = GetiParam("mToPlot");
//
//E
//

#ifdef MPICODE
    }
    MPI_Bcast(&cmd, sizeof(cmdline_data), MPI_BYTE, 0, MPI_COMM_WORLD);
#endif


}

//E

#undef parameter_null
#endif // end of GETPARAM

local int startrun_Common(void)
{
#ifdef SAVERESTORE
    if (strnull(cmd.restorefile)) {
#endif

        setFilesDirs();
            setFilesDirs_log();
        strcpy(gd.mode,"w");
#ifdef MPICODE
        if(ThisTask==0) {
#endif
            if(!(gd.outlog=fopen(gd.logfilePath, gd.mode)))
                error("\nstart_Common: error opening file '%s' \n",gd.logfilePath);
#ifdef MPICODE
        }
#endif


        startrun_getParamsSpecial();
        

#ifdef CLASSLIB
        class_call(random_init(cmd.seed), errmsg, errmsg);
#else
        random_init(cmd.seed);
#endif

        CheckParameters();
#ifdef CLASSLIB
        class_call(startrun_memoryAllocation(), errmsg, errmsg);
#else
        startrun_memoryAllocation();
#endif

        if (! strnull(cmd.infile)) {
#ifdef CLASSLIB
            class_call(infilefmt_string_to_int(cmd.infilefmt, &gd.infilefmt_int),
                   errmsg, errmsg);
            class_call(inputdata(), errmsg, errmsg);
#else
            infilefmt_string_to_int(cmd.infilefmt, &gd.infilefmt_int);
            inputdata();
#endif
        } else {
            testdata();
        }

#ifdef MPICODE
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Bcast(&gd, sizeof(global_data), MPI_BYTE, 0, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);

        if (ThisTask!=0)
            bodytab = (bodyptr) allocate(cmd.nbody * sizeof(body));

        MPI_Bcast(bodytab, cmd.nbody*sizeof(body), MPI_BYTE, 0, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
#endif

        gd.bytes_tot += cmd.nbody*sizeof(body);
        verb_print(cmd.verbose, "\n\nAllocated %g MByte for particle storage.\n",
                   cmd.nbody*sizeof(body)*INMB);

        search_method_string_to_int(cmd.searchMethod, &gd.searchMethod_int);
        gd.rSize = 1.0;
        expandbox(bodytab, cmd.nbody);
        if (cmd.rangeN > gd.rSize)
            verb_print(cmd.verbose,
                "\nstartrun_Common: warning! rangeN (%g) is greather than rSize (%g) of the system...\n",
                cmd.rangeN, gd.rSize);

//        setFilesDirs();

//B Tree search:
        gd.Rcut = cmd.rangeN;                       // Maximum search radius
        gd.RcutSq = gd.Rcut*gd.Rcut;
//B This is how to use logscale bins
//        rlog = rlog10(cmd.rmin) + gd.deltaR*((real)(i));
//        rv = rpow(10.0,rlog);
//E
#ifdef LOGHIST      // This segment is ignored:
//        if (cmd.logHist) {
/*
        gd.deltaRV = dvector(1,cmd.sizeHistN);
        verb_log_print(cmd.verbose_log, gd.outlog, "deltaRV:\n");
        int n;
        for (n=1; n<=cmd.sizeHistN; n++) {
//            gd.deltaRV[n] = cmd.rangeN*rpow(10.0, ((real)(n-cmd.sizeHistN))/NLOGBINPD)
//                            *(rpow(10.0,1.0/NLOGBINPD) - 1.0);
            gd.deltaRV[n] = cmd.rangeN*rpow( 10.0, ( (real)(n-cmd.sizeHistN) )/NLOGBINPD );
            verb_log_print(cmd.verbose_log, gd.outlog, " %d %lg\n",n,gd.deltaRV[n]);
        }
*/
        real rBin, rbinlog;
        if (cmd.rminHist==0) {
//                gd.deltaR = rlog10(cmd.rangeN)/cmd.sizeHistN;
            gd.deltaRV = dvector(1,cmd.sizeHistN);
            verb_log_print(cmd.verbose_log, gd.outlog, "\ndeltaRV:\n");
            int n;
            for (n=1; n<=cmd.sizeHistN; n++) {
                gd.deltaRV[n] = cmd.rangeN*rpow( 10.0, ( (real)(n-cmd.sizeHistN) )/NLOGBINPD );
                verb_log_print(cmd.verbose_log, gd.outlog, " %d %lg\n",n,gd.deltaRV[n]);
            }
        } else {
            gd.deltaR = rlog10(cmd.rangeN/cmd.rminHist)/cmd.sizeHistN;
            gd.deltaRV = dvector(1,cmd.sizeHistN);
            verb_log_print(cmd.verbose_log, gd.outlog, "deltaRV (deltaR=%lf logscale):\n", gd.deltaR);
            int n;
            for (n=1; n<=cmd.sizeHistN; n++) {
                rbinlog = rlog10(cmd.rminHist) + ((real)(n))*gd.deltaR;
                rBin=rpow(10.0,rbinlog);
                gd.deltaRV[n] = rBin;
                verb_log_print(cmd.verbose_log, gd.outlog, " %d %lg\n",n,gd.deltaRV[n]);
            }
        }
//        } else
#else // ! LOGHIST
            gd.deltaR = (cmd.rangeN-cmd.rminHist)/cmd.sizeHistN;
        verb_log_print(cmd.verbose_log, gd.outlog, "deltaR=%lf normal scale):\n",gd.deltaR);
#endif

//
//B For 3pcf brute force:
//        gd.deltaTheta = TWOPI/cmd.sizeHistTheta;
        gd.deltaTheta = PI/cmd.sizeHistTheta;
//E
        MULVS(gd.cells, gd.Box, 1.0/gd.Rcut);       // Only needed for cellmethod
        AllocMem(gd.cellList, VProd (gd.cells)      // Only needed for cellmethod
                + cmd.nbody, INTEGER);
        gd.bytes_tot += (VProd(gd.cells)+cmd.nbody)*sizeof(INTEGER);
//E
        real Vol = 1.0;
        int k;
        DO_COORD(k)
            Vol = Vol*gd.Box[k];
        
        gd.i_deltaR = 1.0/gd.deltaR;     // This is gd.i_r_max change...

        verb_log_print(cmd.verbose_log, gd.outlog,
                       "\nRcut, deltaR: %g %g\n",gd.Rcut,gd.deltaR);
#if NDIM == 3
        verb_log_print(cmd.verbose_log, gd.outlog,
                       "lbox: %g %g %g\n",gd.Box[0],gd.Box[2],gd.Box[2]);
#else
        verb_log_print(cmd.verbose_log, gd.outlog,
                       "lbox: %g %g\n",gd.Box[0],gd.Box[1]);
#endif
        verb_log_print(cmd.verbose_log, gd.outlog, "Box volume = %e\n",Vol);
        verb_log_print(cmd.verbose_log, gd.outlog,
                       "(V/N)^(1/3): %g\n\n",rpow(Vol/cmd.nbody,1.0/3.0));



#ifdef SAVERESTORE
    }    else {                              // if !strnull(cmd.restorefile)

// We must check the order of memory allocation and dealocation
//NOLSST
        restorestate(cmd.restorefile);

#ifdef CLASSLIB
        class_call(random_init(cmd.seed), errmsg, errmsg);
#else
        random_init(cmd.seed);
#endif

        setFilesDirs_log();
        strcpy(gd.mode,"a");
#ifdef MPICODE
        if(ThisTask==0) {
#endif
        if(!(gd.outlog=fopen(gd.logfilePath, gd.mode)))
            error("\nstart_Common: error opening file '%s' \n",gd.logfilePath);
#ifdef MPICODE
        }
#endif
        verb_log_print(cmd.verbose_log, gd.outlog, "\n\nAdded after restart from restart file\n");
        verb_log_print(cmd.verbose_log, gd.outlog, "\nnbody=%d\n",cmd.nbody);

        setFilesDirs();

// We must check the order of memory allocation and dealocation
        AllocMem(gd.cellList, VProd (gd.cells)      // Only needed for cellmethod
                + cmd.nbody, INTEGER);
        gd.bytes_tot += (VProd(gd.cells)+cmd.nbody)*sizeof(INTEGER);

    }
#endif // ! SAVERESTORE


    return _SUCCESS_;
}

local int startrun_getParamsSpecial(void)
{
    char *pch;
    int nitems, ndummy=1;
    char inputnametmp[MAXLENGTHOFSTRSCMD];

    if (strnull(cmd.infile)) {
        if (cmd.verbose_log>=3)
       verb_log_print(cmd.verbose_log, gd.outlog,
                       "\nstartrun_getParamsSpecial: no inputfile was given making data ...\n");
        gd.ninfiles=1;                            // To test data...
    } else {
        strcpy(inputnametmp,cmd.infile);
        if (cmd.verbose_log>=3)
        verb_log_print(cmd.verbose_log, gd.outlog,
                       "\nSplitting string \"%s\" in tokens:\n",inputnametmp);
        gd.ninfiles=0;
        pch = strtok(inputnametmp," ,");
        while (pch != NULL) {
            gd.infilenames[gd.ninfiles] = (string) malloc(MAXLENGTHOFFILES);
            strcpy(gd.infilenames[gd.ninfiles],pch);
            ++gd.ninfiles;
            if (cmd.verbose_log>=3)
            verb_log_print(cmd.verbose_log, gd.outlog, "%s\n",gd.infilenames[gd.ninfiles-1]);
            pch = strtok (NULL, " ,");
        }
        if (cmd.verbose_log>=3)
        verb_log_print(cmd.verbose_log, gd.outlog,
                       "num. of files in infile %s =%d\n",cmd.infile,gd.ninfiles);
    }

    if (!strnull(cmd.infilefmt)) {
        strcpy(inputnametmp,cmd.infilefmt);
        if (cmd.verbose_log>=3)
        verb_log_print(cmd.verbose_log, gd.outlog,
                       "\nSplitting string \"%s\" in tokens:\n",inputnametmp);
        nitems=0;
        pch = strtok(inputnametmp," ,");
        while (pch != NULL) {
            gd.infilefmtname[nitems] = (string) malloc(30);
            strcpy(gd.infilefmtname[nitems],pch);
            ++nitems;
            if (cmd.verbose_log>=3)
            verb_log_print(cmd.verbose_log, gd.outlog, "%s\n",gd.infilefmtname[nitems-1]);
            pch = strtok (NULL, " ,");
        }
        if (cmd.verbose_log>=3)
        verb_log_print(cmd.verbose_log, gd.outlog,
                       "num. of items in infilefmt %s =%d\n",cmd.infilefmt,nitems);
        if (nitems != gd.ninfiles)
            error("\nstartrun_Common: nitems must be equal to number of files\n\n");
    }

//    scaniOption(cmd.nsmooth, gd.nsmooth, &nitems, ndummy, 2, "nsmooth");
    scaniOption(cmd.nsmooth, gd.nsmooth, &nitems, 1, 1, "nsmooth");
//#ifdef BALLS
    scaniOption(cmd.ncritical, gd.ncritical, &nitems, ndummy, 2, "ncritical");
//#endif

    return _SUCCESS_;
}

local int scaniOption(string optionstr, int *option, int *noption,
    int nfiles, int flag, string message)
{
    char *pch;
    char *poptionstr[30],  optiontmp[100];
    int i;
//
// Decide what is better: DEBUG or if (cmd.verbose_log>=3)
//
#ifdef DEBUG
    if (cmd.verbose_log>=3)
    verb_log_print(cmd.verbose_log, gd.outlog, "\nProcessing '%s' option:\n", message);
#endif

    if (!strnull(optionstr)) {
        strcpy(optiontmp,optionstr);
        if (cmd.verbose_log>=3)
        verb_log_print(cmd.verbose_log, gd.outlog, "\nSplitting string \"%s\" in tokens:\n",optiontmp);
        *noption=0;
        pch = strtok(optiontmp," ,");
        while (pch != NULL) {
            poptionstr[*noption] = (string) malloc(10);
            strcpy(poptionstr[*noption],pch);
            ++(*noption);
            if (cmd.verbose_log>=3)
            verb_log_print(cmd.verbose_log, gd.outlog, "%s\n",poptionstr[*noption-1]);
            pch = strtok (NULL, " ,");
        }
        if (cmd.verbose_log>=3)
        verb_log_print(cmd.verbose_log, gd.outlog, "num. of tokens in option %s =%d\n", optionstr,*noption);

        if (flag == 0)
            if (*noption != nfiles)
                error("\nscanOption: noption = %d must be equal to number of files\n\n",*noption);
        if (*noption > MAXITEMS)
            error("\nscanOption: noption = %d must be less than the maximum num. of lines\n\n",*noption);

        for (i=0; i<*noption; i++) {
            option[i]=atoi(poptionstr[i]);
            if (cmd.verbose_log>=3)
            verb_log_print(cmd.verbose_log, gd.outlog, "option: %d\n",option[i]);
        }
        if (cmd.verbose_log>=3)
        verb_log_print(cmd.verbose_log, gd.outlog, "\nnoptions, nfiles: %d %d\n",*noption,nfiles);
        if (flag == 1) {
            if (*noption > nfiles)
                error("\nscanOption: noption = %d must be less or equal to number of files\n\n",*noption);
            else {
                for (i=*noption; i<nfiles; i++) {
                    option[i]=option[i-1]+1;
                    if (cmd.verbose_log>=3)
                    verb_log_print(cmd.verbose_log, gd.outlog, "option: %d\n",option[i]);
                }
            }
        }
    } else {
        for (i=0; i<nfiles; i++) {
            option[i]=1;
            if (cmd.verbose_log>=3)
            verb_log_print(cmd.verbose_log, gd.outlog, "option: %d\n",option[i]);
        }
    }

    return _SUCCESS_;
}

#ifdef GETPARAM
//B Section of parameter stat
local void startrun_ParamStat(void)
{
// Every item in cmdline_defs.h must have an item here::

#ifdef MPICODE
    if (ThisTask==0) {                              // Input all parameters on proccess 0
#endif
    if (GetParamStat("searchMethod") & ARGPARAM)
        cmd.searchMethod = GetParam("searchMethod");

    if (GetParamStat("theta") & ARGPARAM)
        cmd.theta = GetdParam("theta");
    if (GetParamStat("infile") & ARGPARAM)
        cmd.infile = GetParam("infile");
    if (GetParamStat("infilefmt") & ARGPARAM)
        cmd.infilefmt = GetParam("infileformat");
    if (GetParamStat("rootDir") & ARGPARAM)
        cmd.rootDir = GetParam("rootDir");
    if (GetParamStat("outfile") & ARGPARAM)
        cmd.outfile = GetParam("outfile");
    if (GetParamStat("outfilefmt") & ARGPARAM)
        cmd.outfilefmt = GetParam("outfileformat");

    if (GetParamStat("mChebyshev") & ARGPARAM)
        cmd.mchebyshev = GetiParam("mChebyshev");
//B Parameters to set a region in the sky, for example for Takahasi data set.
        if (GetParamStat("thetaL") & ARGPARAM)
            cmd.thetaL = GetdParam("thetaL");
        if (GetParamStat("thetaR") & ARGPARAM)
            cmd.thetaR = GetdParam("thetaR");
        if (GetParamStat("phiL") & ARGPARAM)
            cmd.phiL = GetdParam("phiL");
        if (GetParamStat("thetaR") & ARGPARAM)
            cmd.phiR = GetdParam("phiR");
//E
//    if (GetParamStat("dimension") & ARGPARAM)
//        cmd.dimension = GetiParam("dimension");

    if (GetParamStat("sizeHistN") & ARGPARAM)
        cmd.sizeHistN = GetiParam("sizeHistN");
    if (GetParamStat("rangeN") & ARGPARAM)
        cmd.rangeN = GetdParam("rangeN");
    if (GetParamStat("rminHist") & ARGPARAM)
        cmd.rminHist = GetdParam("rminHist");
 //   if (GetParamStat("logHist") & ARGPARAM)
 //       cmd.logHist = GetdParam("logHist");
    if (GetParamStat("sizeHistTheta") & ARGPARAM)
        cmd.sizeHistTheta = GetiParam("sizeHistTheta");

    if (GetParamStat("histNFileName") & ARGPARAM)
        cmd.histNFileName = GetParam("histNFileName");
    if (GetParamStat("histXi2pcfFileName") & ARGPARAM)
        cmd.histXi2pcfFileName = GetParam("histXi2pcfFileName");
    if (GetParamStat("histZetaMFileName") & ARGPARAM)
        cmd.histZetaMFileName = GetParam("histZetaMFileName");
    if (GetParamStat("mhistZetaFileName") & ARGPARAM)
        cmd.mhistZetaFileName = GetParam("mhistZetaFileName");
    if (GetParamStat("suffixOutFiles") & ARGPARAM)
        cmd.suffixOutFiles = GetParam("suffixOutFiles");

    if (GetParamStat("verbose") & ARGPARAM)
        cmd.verbose = GetiParam("verbose");
    if (GetParamStat("verbose_log") & ARGPARAM)
        cmd.verbose_log = GetiParam("verbose_log");

    if (GetParamStat("numberThreads") & ARGPARAM)
        cmd.numthreads = GetiParam("numberThreads");

    if (GetParamStat("options") & ARGPARAM)
        cmd.options = GetParam("options");

//
//B NOLSST:
//
    if (GetParamStat("seed") & ARGPARAM)
        cmd.seed = GetiParam("seed");
    if (GetParamStat("nsmooth") & ARGPARAM)
        cmd.nsmooth = GetParam("nsmooth");
//        cmd.nsmooth = GetiParam("nsmooth");
#ifdef BALLS
    if (GetParamStat("ntosave") & ARGPARAM)
        cmd.ntosave = GetiParam("ntosave");
    if (GetParamStat("scanLevel") & ARGPARAM)
            cmd.scanLevel = GetiParam("scanLevel");
#endif
    if (GetParamStat("stepNodes") & ARGPARAM)
        cmd.stepNodes = GetiParam("stepNodes");
    if (GetParamStat("ncritical") & ARGPARAM)
        cmd.ncritical = GetParam("ncritical");
//#endif
    if (GetParamStat("testmodel") & ARGPARAM)
        cmd.testmodel = GetParam("testmodel");
    if (GetParamStat("nbody") & ARGPARAM)
        cmd.nbody = GetiParam("nbody");
    if (GetParamStat("lengthBox") & ARGPARAM)
        cmd.lengthBox = GetdParam("lengthBox");
//    if (GetParamStat("mToPlot") & ARGPARAM)
//        cmd.mToPlot = GetiParam("mToPlot");
//
//E
//

#ifdef MPICODE
    }
    MPI_Bcast(&cmd, sizeof(cmdline_data), MPI_BYTE, 0, MPI_COMM_WORLD);
#endif

}
//E
#endif


//B Section of parameter check
local void CheckParameters(void)
{
// If it is necessary: an item in cmdline_defs.h must have an item here::


//    if (!strnull(cmd.restorefile) && !strnull(cmd.infile))
//        fprintf(stdout,"\nCheckParameters: Warning! : %s\n\n",
//            "You are using options restorefile and infile at the same time");

    if (cmd.theta < 0)
        error("CheckParameters: absurd value for theta\n");
    if (cmd.mchebyshev < 1)
        error("CheckParameters: absurd value for mchebyshev\n");
//B Parameters to set a region in the sky, for example for Takahasi data set.
    if (cmd.thetaL < 0 || cmd.thetaL > PI)
        error("CheckParameters: absurd value for thetaL (must be in the range 0--pi)\n");
    if (cmd.thetaR < 0 || cmd.thetaR > PI)
        error("CheckParameters: absurd value for thetaR (must be in the range 0--pi)\n");
    if (cmd.phiL < 0 || cmd.phiL > PI)
        error("CheckParameters: absurd value for phiL (must be in the range 0--2pi)\n");
    if (cmd.phiR < 0 || cmd.phiR > PI)
        error("CheckParameters: absurd value for phiR (must be in the range 0--2pi)\n");
    if (cmd.thetaL > cmd.thetaR)
        error("CheckParameters: absurd value for thetaL (must be greater thatn thetaR)\n");
    if (cmd.phiL > cmd.phiR)
        error("CheckParameters: absurd value for phiL (must be greater thatn phiR)\n");
//E
//    if (cmd.dimension != 2 && cmd.dimension != 3)
//        error("CheckParameters: absurd value for dimension\n");
    if (cmd.sizeHistN < 2)
        error("CheckParameters: absurd value for sizeHistN\n");
    if (cmd.rangeN < 0)
        error("CheckParameters: absurd value for rangeN\n");
    if (cmd.rminHist < 0 || cmd.rminHist > cmd.rangeN)
        error("CheckParameters: absurd value for rminHist\n");
    if (cmd.sizeHistTheta < 2)
        error("CheckParameters: absurd value for sizeHistTheta\n");
    if (cmd.numthreads <= 0)
        error("CheckParameters: absurd value for numberThreads must be an integer >= 0\n");

//
//B NOLSST:
//
    if (gd.nsmooth[0] < 1)
        error("CheckParameters: absurd value for nsmooth\n");
#ifdef BALLS
    if (scanopt(cmd.options, "bodyfound"))
        if (cmd.ntosave < 1 || cmd.ntosave > cmd.nbody)
            error("CheckParameters: absurd value for ntosave\n");
    if (cmd.scanLevel < 0)
        error("CheckParameters: absurd value for scanLevel (%d)\n",cmd.scanLevel);
//    if (cmd.scanLevel > 12)
//        error("CheckParameters: too big value for scanLevel (%d)\n",cmd.scanLevel);
#endif
    if (cmd.stepNodes < 1)
        error("CheckParameters: absurd value for stepNodes\n");
    if (gd.ncritical[0] < 2)
        error("CheckParameters: absurd value for ncritical one\n");
    if (gd.ncritical[1] < 2)
        error("CheckParameters: absurd value for ncritical two\n");
//#endif
    if (cmd.nbody < 3)
        error("CheckParameters: absurd value for nbody\n");
    if (cmd.lengthBox <= 0)
        error("CheckParameters: absurd value for lengthBox\n");
//    if (cmd.mToPlot <= 0)
//        error("CheckParameters: absurd value for mToPlot must be an integer >= 0\n");

    gd.bh86 = scanopt(cmd.options, "bh86"); // Barnes, J. & Hut, P. 1986. Nature 324, 446.
    gd.sw94 = scanopt(cmd.options, "sw94"); // Salmon, J.K. & Warren, M.S. 1994. J. Comp. Phys. 111, 136
    if (gd.bh86 && gd.sw94)
        error("CheckParameters: incompatible options bh86 and sw94\n");
//
//E
//

}
//E

//B Section of parameter reading from a file
local void ReadParameterFile(char *fname)
{
// Every item in cmdline_defs.h must have an item here::
#define DOUBLE 1
#define STRING 2
#define INT 3
#define BOOLEAN 4
#define MAXTAGS 300

    FILE *fd;

  char buf[200],buf1[200],buf2[200],buf3[200];
  int  i,j,nt;
  int  id[MAXTAGS];
  void *addr[MAXTAGS];
  char tag[MAXTAGS][50];
  int  errorFlag=0;

#ifdef MPICODE
    if (ThisTask==0) {                              // Input all parameters on proccess 0
#endif
  nt=0;

    SPName(cmd.searchMethod,"searchMethod",MAXLENGTHOFSTRSCMD);
    RPName(cmd.theta,"theta");
    SPName(cmd.infile,"infile",MAXLENGTHOFSTRSCMD);
    SPName(cmd.infilefmt,"infileformat",MAXLENGTHOFSTRSCMD);
    SPName(cmd.rootDir,"rootDir",MAXLENGTHOFSTRSCMD);
    SPName(cmd.outfile,"outfile",MAXLENGTHOFSTRSCMD);
    SPName(cmd.outfilefmt,"outfileformat",MAXLENGTHOFSTRSCMD);
    IPName(cmd.mchebyshev,"mChebyshev");
//B Parameters to set a region in the sky, for example for Takahasi data set.
        RPName(cmd.thetaL,"thetaL");
        RPName(cmd.thetaR,"thetaR");
        RPName(cmd.phiL,"phiL");
        RPName(cmd.phiR,"phiR");
//E
//    IPName(cmd.dimension,"dimension");
    IPName(cmd.sizeHistN,"sizeHistN");
    RPName(cmd.rangeN,"rangeN");
    RPName(cmd.rminHist,"rminHist");
//    BPName(cmd.logHist,"logHist");
    IPName(cmd.sizeHistTheta,"sizeHistTheta");
    SPName(cmd.histNFileName,"histNFileName",MAXLENGTHOFSTRSCMD);
    SPName(cmd.histXi2pcfFileName,"histXi2pcfFileName",MAXLENGTHOFSTRSCMD);
    SPName(cmd.histZetaMFileName,"histZetaMFileName",MAXLENGTHOFSTRSCMD);
    SPName(cmd.mhistZetaFileName,"mhistZetaFileName",MAXLENGTHOFSTRSCMD);
    SPName(cmd.suffixOutFiles,"suffixOutFiles",MAXLENGTHOFSTRSCMD);
    IPName(cmd.stepState,"stepState");
#ifdef SAVERESTORE
    SPName(cmd.statefile,"statefile",MAXLENGTHOFSTRSCMD);
    SPName(cmd.restorefile,"restorefile",MAXLENGTHOFSTRSCMD);
#endif
    IPName(cmd.verbose,"verbose");
    IPName(cmd.verbose_log,"verbose_log");
    IPName(cmd.numthreads,"numberThreads");
	SPName(cmd.options,"options",MAXLENGTHOFSTRSCMD);

//
//B NOLSST:
//
    IPName(cmd.seed,"seed");                        // to always have defaults
    SPName(cmd.nsmooth,"nsmooth",MAXLENGTHOFSTRSCMD);
#ifdef BALLS
    IPName(cmd.ntosave,"ntosave");
    IPName(cmd.scanLevel,"scanLevel");
#endif
    IPName(cmd.stepNodes,"stepNodes");
    SPName(cmd.ncritical,"ncritical",MAXLENGTHOFSTRSCMD);
//#endif
    SPName(cmd.testmodel,"testmodel",MAXLENGTHOFSTRSCMD);
    IPName(cmd.nbody,"nbody");
    RPName(cmd.lengthBox,"lengthBox");
//    IPName(cmd.mToPlot,"mToPlot");
//
//E
//


	if((fd=fopen(fname,"r"))) {
		while(!feof(fd)) {
			fgets(buf,200,fd);
            if(sscanf(buf,"%s%s%s",buf1,buf2,buf3)<1)
                continue;
            if(sscanf(buf,"%s%s%s",buf1,buf2,buf3)<2)
                *buf2='\0';
            if(buf1[0]=='%')
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
						*((double*)addr[j])=atof(buf2); 
						break;
					case STRING:
						strcpy(addr[j],buf2);
						break;
					case INT:
						*((int*)addr[j])=atoi(buf2);
						break;
					case BOOLEAN:
						if (strchr("tTyY1", *buf2) != NULL) {          
							*((bool*)addr[j])=TRUE;
                        } else 
                            if (strchr("fFnN0", *buf2) != NULL)  {
                                *((bool*)addr[j])=FALSE;
                            } else {
                                error("getbparam: %s=%s not bool\n",buf1,buf2);
                            }
						break;
                }
            } else {
                fprintf(stdout, "Error in file %s: Tag '%s' %s.\n",
					fname, buf1, "not allowed or multiple defined");
                errorFlag=1;
            }
        }
        fclose(fd);
    } else {
        fprintf(stdout,"Parameter file %s not found.\n", fname);
        errorFlag=1;
//        exit(1);
    }
  
    for(i=0;i<nt;i++) {
        if(*tag[i]) {
            fprintf(stdout,
                "Error. I miss a value for tag '%s' in parameter file '%s'.\n",
                tag[i],fname);
            errorFlag=1;
//            exit(0);
        }
    }

#ifdef MPICODE
    }
    MPI_Bcast(&errorFlag, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if(errorFlag) {
        MPI_Finalize();
        exit(0);
    }
    MPI_Bcast(&cmd, sizeof(cmdline_data), MPI_BYTE, 0, MPI_COMM_WORLD);
#else
    if(errorFlag) {
        exit(0);
    }
#endif


#undef DOUBLE 
#undef STRING 
#undef INT 
#undef BOOLEAN
#undef MAXTAGS
}
//E

#define FMTT	"%-35s%s\n"
#define FMTI    "%-35s%d\n"
#define FMTIL    "%-35s%ld\n"
#define FMTR	"%-35s%g\n"

//B Section of parameter writing to a file
local void PrintParameterFile(char *fname)
{
// Every item in cmdline_defs.h must have an item here::
    FILE *fdout;
    char buf[200];
    int  errorFlag=0;

#ifdef MPICODE
    if (ThisTask==0) {                              // Output only on proccess 0
#endif
    sprintf(buf,"%s/%s%s",cmd.rootDir,fname,"-usedvalues");
    if(!(fdout=fopen(buf,"w"))) {
        fprintf(stdout,"error opening file '%s' \n",buf);
        errorFlag=1;
    } else {
        fprintf(fdout,FMTT,"searchMethod",cmd.searchMethod);
        fprintf(fdout,FMTR,"theta",cmd.theta);
        fprintf(fdout,FMTT,"infile",cmd.infile);
        fprintf(fdout,FMTT,"infileformat",cmd.infilefmt);
        fprintf(fdout,FMTT,"rootDir",cmd.rootDir);
        fprintf(fdout,FMTT,"outfile",cmd.outfile);
        fprintf(fdout,FMTT,"outfileformat",cmd.outfilefmt);
        fprintf(fdout,FMTI,"mChebyshev",cmd.mchebyshev);
//B Parameters to set a region in the sky, for example for Takahasi data set.
        fprintf(fdout,FMTR,"thetaL",cmd.thetaL);
        fprintf(fdout,FMTR,"thetaR",cmd.thetaR);
        fprintf(fdout,FMTR,"phiL",cmd.phiL);
        fprintf(fdout,FMTR,"phiR",cmd.phiR);
//E
//        fprintf(fdout,FMTI,"dimension",cmd.dimension);
        fprintf(fdout,FMTI,"sizeHistN",cmd.sizeHistN);
        fprintf(fdout,FMTI,"sizeHistTheta",cmd.sizeHistTheta);
        fprintf(fdout,FMTR,"rangeN",cmd.rangeN);
        fprintf(fdout,FMTR,"rminHist",cmd.rminHist);
//        fprintf(fdout,FMTT,"logHist",cmd.logHist ? "true" : "false");
        fprintf(fdout,FMTT,"histNFileName",cmd.histNFileName);
        fprintf(fdout,FMTT,"histXi2pcfFileName",cmd.histXi2pcfFileName);
        fprintf(fdout,FMTT,"histZetaMFileName",cmd.histZetaMFileName);
        fprintf(fdout,FMTT,"mhistZetaFileName",cmd.mhistZetaFileName);
        fprintf(fdout,FMTT,"suffixOutFiles",cmd.suffixOutFiles);
        fprintf(fdout,FMTIL,"stepState",cmd.stepState);
#ifdef SAVERESTORE
        fprintf(fdout,FMTT,"statefile",cmd.statefile);
        fprintf(fdout,FMTT,"restorefile",cmd.restorefile);
#endif
        fprintf(fdout,FMTI,"verbose",cmd.verbose);
        fprintf(fdout,FMTI,"verbose_log",cmd.verbose_log);
        fprintf(fdout,FMTI,"numberThreads",cmd.numthreads);
        fprintf(fdout,FMTT,"options",cmd.options);
//
//B NOLSST:
//
        fprintf(fdout,FMTI,"seed",cmd.seed);
        fprintf(fdout,FMTT,"nsmooth",cmd.nsmooth);
#ifdef BALLS
        fprintf(fdout,FMTIL,"ntosave",cmd.ntosave);
        fprintf(fdout,FMTI,"scanLevel",cmd.scanLevel);
#endif
        fprintf(fdout,FMTIL,"stepNodes",cmd.stepNodes);
        fprintf(fdout,FMTT,"ncritical",cmd.ncritical);
//#endif
        fprintf(fdout,FMTT,"testmodel",cmd.testmodel);
        fprintf(fdout,FMTIL,"nbody",cmd.nbody);
        fprintf(fdout,FMTR,"lengthBox",cmd.lengthBox);
//        fprintf(fdout,FMTI,"mToPlot",cmd.mToPlot);
//
//E
//
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

}
//E

#undef FMTT
#undef FMTI
#undef FMTR


local int infilefmt_string_to_int(string infmt_str,int *infmt_int)
{
    *infmt_int=-1;
    if (strcmp(infmt_str,"columns-ascii") == 0)             *infmt_int = INCOLUMNS;
    if (strnull(infmt_str))                                 *infmt_int = INNULL;
    if (strcmp(infmt_str,"binary") == 0)                    *infmt_int = INCOLUMNSBIN;
    if (strcmp(infmt_str,"takahasi") == 0)                  *infmt_int = INTAKAHASI;
    if (strcmp(infmt_str,"columns-ascii-2d-to-3d") == 0)    *infmt_int = INCOLUMNS2DTO3D;


    return _SUCCESS_;
}

local int random_init(int seed)
{
#ifndef NOGSL
    const gsl_rng_type * T;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    gd.r = gsl_rng_alloc (T);
    gsl_rng_set(gd.r, seed);
    gd.bytes_tot += (1)*sizeof(gsl_rng_type);
#endif

    return _SUCCESS_;
}

global int startrun_memoryAllocation(void)
{
// Free allocated memory in reverse order as were allocated
    INTEGER bytes_tot_local=0;

    gd.histN = dvector(1,cmd.sizeHistN);
    gd.histCF = dvector(1,cmd.sizeHistN);
    bytes_tot_local += cmd.sizeHistN*sizeof(real);
    gd.histNSub = dvector(1,cmd.sizeHistN);
// 2pcf
    gd.histNSubXi2pcf = dvector(1,cmd.sizeHistN);
//
    bytes_tot_local += cmd.sizeHistN*sizeof(real);
    gd.histNNN = dvector(1,cmd.sizeHistN);
    bytes_tot_local += cmd.sizeHistN*sizeof(real);
    gd.histNNNSub = dmatrix3D(1,cmd.sizeHistN,1,cmd.sizeHistN,1,cmd.sizeHistTheta);
    bytes_tot_local += (cmd.sizeHistN*cmd.sizeHistN*cmd.sizeHistTheta)*sizeof(real);
    gd.histXi2pcf = dvector(1,cmd.sizeHistN);
    bytes_tot_local += cmd.sizeHistN*sizeof(real);
    gd.histXi3pcf = dmatrix3D(1,cmd.sizeHistN,1,cmd.sizeHistN,1,cmd.sizeHistTheta);
    bytes_tot_local += (cmd.sizeHistN*cmd.sizeHistN*cmd.sizeHistTheta)*sizeof(real);
#ifdef TPCF
    gd.histXi = dmatrix(1,cmd.mchebyshev+1,1,cmd.sizeHistN);
    bytes_tot_local += (cmd.mchebyshev+1)*cmd.sizeHistN*sizeof(real);
    gd.histXicos = dmatrix(1,cmd.mchebyshev+1,1,cmd.sizeHistN);
    bytes_tot_local += (cmd.mchebyshev+1)*cmd.sizeHistN*sizeof(real);
    gd.histXisin = dmatrix(1,cmd.mchebyshev+1,1,cmd.sizeHistN);
    bytes_tot_local += (cmd.mchebyshev+1)*cmd.sizeHistN*sizeof(real);
    gd.histZetaM = dmatrix3D(1,cmd.mchebyshev+1,1,cmd.sizeHistN,1,cmd.sizeHistN);
    bytes_tot_local += (cmd.mchebyshev+1)*cmd.sizeHistN*cmd.sizeHistN*sizeof(real);
    gd.histZetaMcos = dmatrix3D(1,cmd.mchebyshev+1,1,cmd.sizeHistN,1,cmd.sizeHistN);
    bytes_tot_local += (cmd.mchebyshev+1)*cmd.sizeHistN*cmd.sizeHistN*sizeof(real);
    gd.histZetaMsin = dmatrix3D(1,cmd.mchebyshev+1,1,cmd.sizeHistN,1,cmd.sizeHistN);
    bytes_tot_local += (cmd.mchebyshev+1)*cmd.sizeHistN*cmd.sizeHistN*sizeof(real);
    gd.histZetaMsincos = dmatrix3D(1,cmd.mchebyshev+1,1,cmd.sizeHistN,1,cmd.sizeHistN);
    bytes_tot_local += (cmd.mchebyshev+1)*cmd.sizeHistN*cmd.sizeHistN*sizeof(real);

#ifndef NOGSL
    gd.histXi_gsl = gsl_matrix_complex_calloc(cmd.mchebyshev+1,cmd.sizeHistN);
    bytes_tot_local += (cmd.mchebyshev+1)*cmd.sizeHistN*2*sizeof(real);

    histZetaMatrix = (mMatrix_ptr) allocate((cmd.mchebyshev+1) * sizeof(mMatrix));
    int m;
    for (m=0; m<=cmd.mchebyshev; m++){
        histZetaMatrix[m].histZetaM = gsl_matrix_complex_calloc(cmd.sizeHistN,cmd.sizeHistN);
    }
    bytes_tot_local += (cmd.mchebyshev+1)*cmd.sizeHistN*cmd.sizeHistN*2*sizeof(real);
#endif

#endif // ! TPCF

    gd.bytes_tot += bytes_tot_local;
    verb_print(cmd.verbose,
               "\n\nstartrun_memoryAllocation: Allocated %g MByte for histograms storage.\n",
               bytes_tot_local*INMB);

    return _SUCCESS_;
}

local void search_method_string_to_int(string method_str,int *method_int)
{
// Every search method must have an item here::
    *method_int=-1;
    if (strnull(method_str))                                *method_int = SEARCHNULL;
    if (strcmp(method_str,"tree-omp") == 0)                 *method_int = TREEOMPMETHOD;
    if (strcmp(method_str,"tree-3pcf-direct-omp") == 0)     *method_int = TREE3PCFBFOMPMETHOD;
    if (strcmp(method_str,"tree-omp-sincos") == 0)          *method_int = TREEOMPMETHODSINCOS;
    if (strcmp(method_str,"balls-omp") == 0)               *method_int = BALLSOMPMETHOD;



}

local int expandbox(bodyptr btab, int nbody)
{
    real dmax, d;
    bodyptr p;
    int k;
    bodyptr root;
    
    root = (bodyptr) allocate(1 * sizeof(body));
    gd.bytes_tot += (1)*sizeof(body);

    CLRV(Pos(root));        // Assuming that (0,0,...) is the center of the Box

    dmax = 0.0;
    DO_BODY(p, btab, btab+nbody)
        DO_COORD(k) {
            d = rabs(Pos(p)[k] - Pos(root)[k]);
            if (d > dmax)
                dmax = d;
        }
    while (gd.rSize < 2 * dmax)
      gd.rSize = 2 * gd.rSize;

#ifdef DEBUG
    if (cmd.verbose_log>=3)
    verb_log_print(cmd.verbose_log, gd.outlog,"\nexpandbox: rSize = %g\n", gd.rSize);
#endif

    return _SUCCESS_;
}

#ifdef OPENMPCODE
int set_number_threads(void)
{
    omp_set_num_threads(cmd.numthreads);

    return _SUCCESS_;
}
#endif
