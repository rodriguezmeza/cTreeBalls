/*==============================================================================
 MODULE: startrun.c				[cTreeBalls]
 Written by: Mario A. Rodriguez-Meza
 Starting date: april 2023
 Purpose: routines to initialize the main code
 Language: C
 Use: 'StartRun();'
 
 Mayor revisions:
 ==============================================================================*/
//        1          2          3          4        ^ 5          6          7

//
// We must check the order of memory allocation and deallocation!!!
// Here and in EndRun in cballsio.c
//

#include "globaldefs.h"

local void ReadParameterFile(struct  cmdline_data*, struct  global_data*, char *);
local int startrun_parameterfile(struct  cmdline_data*, struct  global_data*);
local int startrun_cmdline(struct  cmdline_data*, struct  global_data*);
local void ReadParametersCmdline(struct  cmdline_data*, struct  global_data*);
local void ReadParametersCmdline_short(struct  cmdline_data*, 
                                       struct  global_data*);
local int CheckParameters(struct  cmdline_data*, struct  global_data*);
local int random_init(struct  cmdline_data*, struct  global_data*, int);
local void search_method_string_to_int(string method_str,int *method_int);

//B Special routines to scan command line
local int startrun_getParamsSpecial(struct  cmdline_data*, struct  global_data*);
local int scaniOption(struct  cmdline_data*, struct  global_data*,
                      string, int *, int *, int, int, string);
local int scanrOption(struct  cmdline_data*, struct  global_data*,
                      string, double *, int *, int, int, string);
//E

#ifndef USEGSL
local long saveidum;
#endif

#ifdef ADDONS
#include "startrun_include_00.h"
#endif


/*
 StartRun routine:

 To be called in main:
 StartRun(&cmd, &gd, argv[0], HEAD1, HEAD2, HEAD3);
 
 This routine is in charge of setting all global structures in order to
    the comutation process run smoothly with all parameters given
    by the user set and checked.

 Arguments:
    * `cmd`: Input: structure cmdline_data pointer
    * `gd`: Input: structure global_data pointer
    * `head0`: Input: string
    * `head1`: Input: string
    * `head2`: Input: string
    * `head3`: Input: string
 Return (the error status):
    int SUCCESS or FAILURE
 */
#ifndef CLASSLIB
int StartRun(struct  cmdline_data* cmd, struct  global_data* gd, 
             string head0, string head1, string head2, string head3)
{
    double cpustart = CPUTIME;

    gd->headline0 = head0; gd->headline1 = head1;
    gd->headline2 = head2; gd->headline3 = head3;

    printf("\n%s\n%s: %s\n\t %s\n",
           gd->headline0, gd->headline1, gd->headline2, gd->headline3);
    printf("Version = %s\n", getversion());

    gd->stopflag = 0;
    gd->cputotalinout = 0.;
    gd->cputotal = 0.;
    gd->bytes_tot = 0;
    gd->sameposcount = 0;

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

    gd->bytes_tot += sizeof(struct  global_data);
    gd->bytes_tot += sizeof(struct cmdline_data);
    verb_print(cmd->verbose,
               "\nStartRun: Total allocated %g MByte storage so far.\n",
               gd->bytes_tot*INMB);

//B If uncommented there will be a warning in the setup.py process
//#ifdef OPENMPCODE
    class_call_cballs(SetNumberThreads(cmd), errmsg, errmsg);
//#endif
//E
    gd->cputotalinout += CPUTIME - cpustart;
    verb_print(cmd->verbose, "\nStartRun CPU time: %g %s\n",
               CPUTIME - cpustart, PRNUNITOFTIMEUSED);

    return SUCCESS;
}

#else // ! CLASSLIB

#ifdef USEGSL
#error `USEGSL` and `CLASSLIB` can not used at the same time. Switched off one of them
#endif

#include "input.h"

int StartRun(struct  cmdline_data* cmd, struct  global_data* gd,
             string head0, string head1, string head2, string head3)
{
    struct file_content fc;

    double cpustart = CPUTIME;
    
    gd->headline0 = head0; gd->headline1 = head1;
    gd->headline2 = head2; gd->headline3 = head3;
    printf("\n%s\n%s: %s\n\t %s\n",
           gd->headline0, gd->headline1, gd->headline2, gd->headline3);
    printf("Version = %s\n", getversion());
    
    gd->stopflag = 0;
    gd->cputotalinout = 0.;
    gd->cputotal = 0.;
    gd->bytes_tot = 0;
    
#ifdef GETPARAM
    cmd->paramfile = GetParam("paramfile");
    if (*(cmd->paramfile)=='-')
        error("bad parameter %s\n", cmd->paramfile);
    if (!strnull(cmd->paramfile)) {
        class_call_cballs(input_find_file(cmd->paramfile, &fc, errmsg), 
                          errmsg, errmsg);
        class_call_cballs(input_read_from_file(cmd, &fc, errmsg), errmsg, errmsg);
        class_call_cballs(parser_free(&fc), errmsg, errmsg);
    } else {
        startrun_cmdline(cmd, gd);
    }

#else
    class_call_cballs(input_find_file(cmd->ParameterFile, &fc, errmsg), 
                      errmsg, errmsg);
    class_call_cballs(input_read_from_file(cmd, &fc, errmsg), errmsg, errmsg);
    class_call_cballs(parser_free(&fc), errmsg, errmsg);
#endif

    if (!strnull(cmd->paramfile))
        class_call_cballs(StartRun_Common(cmd, gd), errmsg, errmsg);

#ifdef GETPARAM
    if (!strnull(cmd->paramfile))
        PrintParameterFile(cmd, cmd->paramfile);
#else
    PrintParameterFile(cmd, cmd->ParameterFile);
#endif

    gd->bytes_tot += sizeof(struct  global_data);
    gd->bytes_tot += sizeof(struct cmdline_data);
    verb_print(cmd->verbose,
               "\nStartRun: Total allocated %g MByte storage so far.\n",
               gd->bytes_tot*INMB);

#ifdef OPENMPCODE
    class_call_cballs(SetNumberThreads(cmd), errmsg, errmsg);
#endif

    gd->cputotalinout += CPUTIME - cpustart;
    verb_print(cmd->verbose, "\nStartRun CPU time: %g %s\n",
               CPUTIME - cpustart, PRNUNITOFTIMEUSED);

    return SUCCESS;
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
#endif // ! GETPARAM

	StartRun_Common(cmd, gd);
#ifdef GETPARAM
	PrintParameterFile(cmd, cmd->paramfile);
#else
    PrintParameterFile(cmd, cmd->ParameterFile);
#endif

    return SUCCESS;
}


#ifdef GETPARAM
#define parameter_null	"parameters_null-cballs"

//B Section for reading parameters from the command line

local int startrun_cmdline(struct  cmdline_data* cmd, struct  global_data* gd)
{
	ReadParametersCmdline(cmd, gd);
	StartRun_Common(cmd, gd);
	PrintParameterFile(cmd, parameter_null);

    return SUCCESS;
}

local void ReadParametersCmdline(struct  cmdline_data* cmd, 
                                 struct  global_data* gd)
{
// Every item in cmdline_defs.h must have an item here::

    //B Parameters related to the searching method
    cmd->searchMethod = GetParam("searchMethod");
    cmd->mChebyshev = GetiParam("mChebyshev");
    cmd->nsmooth = GetParam("nsmooth");
    cmd->rsmooth = GetParam("rsmooth");
    cmd->theta = GetdParam("theta");
    cmd->computeTPCF = GetbParam("computeTPCF");
    cmd->computeShearCF = GetbParam("computeShearCF");
    cmd->usePeriodic = GetbParam("usePeriodic");
    //E

    //B Parameters about the I/O file(s)
    // Input catalog parameters
    cmd->infile = GetParam("infile");
    cmd->infilefmt = GetParam("infileformat");
    cmd->iCatalogs = GetParam("iCatalogs");
    // Output parameters
    cmd->rootDir = GetParam("rootDir");
    cmd->outfile = GetParam("outfile");
    cmd->outfilefmt = GetParam("outfileformat");
    // Parameters to set a region in the sky, for example for Takahasi data set
    cmd->thetaL = GetdParam("thetaL");
    cmd->thetaR = GetdParam("thetaR");
    cmd->phiL = GetdParam("phiL");
    cmd->phiR = GetdParam("phiR");
    //E

    //B Parameters to control histograms and their output files
    cmd->useLogHist = GetbParam("useLogHist");
    cmd->logHistBinsPD = GetiParam("logHistBinsPD");
    //
    cmd->sizeHistN = GetiParam("sizeHistN");
    cmd->rangeN = GetdParam("rangeN");
    cmd->rminHist = GetdParam("rminHist");
    cmd->sizeHistPhi = GetiParam("sizeHistPhi");
    //
    cmd->histNNFileName = GetParam("histNNFileName");
    cmd->histXi2pcfFileName = GetParam("histXi2pcfFileName");
    cmd->histZetaFileName = GetParam("histZetaFileName");
    cmd->suffixOutFiles = GetParam("suffixOutFiles");
    //E

    //B Set of parameters needed to construct a test model
    cmd->seed=GetiParam("seed");                    // to always have defaults
                                                    //  Check in gsl
    cmd->testmodel = GetParam("testmodel");
#ifdef LONGINT
    cmd->nbody = GetlParam("nbody");
#else
    cmd->nbody = GetiParam("nbody");
#endif
    cmd->lengthBox = GetdParam("lengthBox");
    //E

    //B Miscellaneous parameters
    cmd->script = GetParam("script");
#ifdef LONGINT
        cmd->stepState = GetlParam("stepState");
#else
        cmd->stepState = GetiParam("stepState");
#endif
    cmd->verbose = GetiParam("verbose");
    cmd->verbose_log = GetiParam("verbose_log");
#ifdef OPENMPCODE
    cmd->numthreads = GetiParam("numberThreads");
#endif
    cmd->options = GetParam("options");
    //E

#ifdef ADDONS
#include "startrun_include_02.h"
#endif
}

local void ReadParametersCmdline_short(struct  cmdline_data* cmd, struct  global_data* gd)
{
//B Here add parameters needed to be read after reading parameter file
//    cmd->script = GetParam("script");
//E
}

//E

#undef parameter_null
#endif // end of GETPARAM

//B Section of parameter reading from a file
local void ReadParameterFile(struct  cmdline_data* cmd, 
                             struct  global_data* gd, char *fname)
{
// Every item in cmdline_defs.h must have an item here::
#define DOUBLE 1
#define STRING 2
#define INT 3
#define LONG 6
#define BOOLEAN 4
#define MAXTAGS 300
#define MAXCHARBUF 1024

    FILE *fd;

  int  i,j,nt;
  int  id[MAXTAGS];
  void *addr[MAXTAGS];
  char tag[MAXTAGS][50];
  int  errorFlag=0;

  nt=0;

    //B Parameters related to the searching method
    SPName(cmd->searchMethod,"searchMethod",MAXLENGTHOFSTRSCMD);
    IPName(cmd->mChebyshev,"mChebyshev");
    SPName(cmd->nsmooth,"nsmooth",MAXLENGTHOFSTRSCMD);
    SPName(cmd->rsmooth,"rsmooth",MAXLENGTHOFSTRSCMD);
    RPName(cmd->theta,"theta");
    BPName(cmd->computeTPCF,"computeTPCF");
    BPName(cmd->computeShearCF,"computeShearCF");
    BPName(cmd->usePeriodic,"usePeriodic");
    //E

    //B Parameters to control the I/O file(s)
    // Input catalog parameters
    SPName(cmd->infile,"infile",MAXLENGTHOFSTRSCMD);
    SPName(cmd->infilefmt,"infileformat",MAXLENGTHOFSTRSCMD);
    SPName(cmd->iCatalogs,"iCatalogs",MAXLENGTHOFSTRSCMD);
    // Output parameters
    SPName(cmd->rootDir,"rootDir",MAXLENGTHOFSTRSCMD);
    SPName(cmd->outfile,"outfile",MAXLENGTHOFSTRSCMD);
    SPName(cmd->outfilefmt,"outfileformat",MAXLENGTHOFSTRSCMD);
    //B Parameters to set a region in the sky, for example for Takahasi data set.
    RPName(cmd->thetaL,"thetaL");
    RPName(cmd->thetaR,"thetaR");
    RPName(cmd->phiL,"phiL");
    RPName(cmd->phiR,"phiR");
    //E

    //B Parameters to control histograms and their output files
    BPName(cmd->useLogHist,"useLogHist");
    IPName(cmd->logHistBinsPD,"logHistBinsPD");
    //
    IPName(cmd->sizeHistN,"sizeHistN");
    RPName(cmd->rangeN,"rangeN");
    RPName(cmd->rminHist,"rminHist");
    IPName(cmd->sizeHistPhi,"sizeHistPhi");
    //
    SPName(cmd->histNNFileName,"histNNFileName",MAXLENGTHOFSTRSCMD);
    SPName(cmd->histXi2pcfFileName,"histXi2pcfFileName",MAXLENGTHOFSTRSCMD);
    SPName(cmd->histZetaFileName,"histZetaFileName",MAXLENGTHOFSTRSCMD);
    SPName(cmd->suffixOutFiles,"suffixOutFiles",MAXLENGTHOFSTRSCMD);
    //E

    //B Set of parameters needed to construct a test model
    IPName(cmd->seed,"seed");                       // to always have
                                                    //  defaults
    SPName(cmd->testmodel,"testmodel",MAXLENGTHOFSTRSCMD);
#ifdef LONGINT
    LPName(cmd->nbody,"nbody");
#else
    IPName(cmd->nbody,"nbody");
#endif
    RPName(cmd->lengthBox,"lengthBox");
    //E

    //B Miscellaneous parameters
    SPName(cmd->script,"script",MAXLENGTHOFSTRSCMD);
#ifdef LONGINT
    LPName(cmd->stepState,"stepState");
#else
    IPName(cmd->stepState,"stepState");
#endif
    IPName(cmd->verbose,"verbose");
    IPName(cmd->verbose_log,"verbose_log");
#ifdef OPENMPCODE
    IPName(cmd->numthreads,"numberThreads");
#endif
    SPName(cmd->options,"options",MAXLENGTHOFSTRSCMD);
    //E
    
#ifdef ADDONS
#include "startrun_include_03.h"
#endif

    size_t slen;
    char *script1;
    char *script2;
    char *script3;
    char *script4;
    char buf4[MAXCHARBUF];
    char buf5[MAXCHARBUF];

//B
#ifndef _LINE_LENGTH_MAX_
#define _LINE_LENGTH_MAX_ 1024
#endif
#define _ARGUMENT_LENGTH_MAX_ 1024
        char line[_LINE_LENGTH_MAX_];
        char name[_ARGUMENT_LENGTH_MAX_];
        char value[_ARGUMENT_LENGTH_MAX_];
        char * phash;
        char * pequal;
        char * left;
        char * right;
//E

    if((fd=fopen(fname,"r"))) {
        while(!feof(fd)) {
//B
            fgets(line,MAXCHARBUF,fd);

            pequal=strchr(line,'=');
            if (pequal == NULL)
                continue;
            phash=strchr(line,'#');
            if ((phash != NULL) && (phash-pequal<2))
                continue;

            left=line;
            while (left[0]==' ') {
              left++;
            }
            if(left[0]=='\'' || left[0]=='\"'){
              left++;
            }
            right=pequal-1;
            while (right[0]==' ') {
              right--;
            }
            if(right[0]=='\'' || right[0]=='\"'){
              right--;
            }

            if (right-left < 0) {
                fprintf(stdout,
        "Error in file %s: there is no variable name before '=' in line: '%s'\n",
                    fname, line);
                errorFlag=1;
                continue;
            }

            strncpy(name,left,right-left+1);
            name[right-left+1]='\0';

            left = pequal+1;
            while (left[0]==' ') {
              left++;
            }

            if (phash == NULL)
              right = line+strlen(line)-1;
            else
              right = phash-1;

            while (right[0]<=' ') {
              right--;
            }

            if (right-left < 0)
                continue;

            strncpy(value,left,right-left+1);
            value[right-left+1]='\0';
//E

            for(i=0,j=-1;i<nt;i++)
                if(strcmp(name,tag[i])==0) {
                    j=i;
                    tag[i][0]=0;
                    break;
                }
            if(j>=0) {
                switch(id[j]) {
                    case DOUBLE:
                        *((double*)addr[j])=atof(value);
                        break;
                    case STRING:
                        if (strcmp(name,"script") == 0){ // To remove both '"'
                            int index;
                            size_t slen;
                                  slen = strlen(value);
                                  cmd->script = (char*) malloc((slen-2)*sizeof(char));
                                  script1 = (char*) malloc(slen*sizeof(char));
                                  memcpy(script1,value,slen);
                                  script2 = strchr(script1, '"');
                                  memcpy(cmd->script,script2+1,slen-2);
                        } else {
                            strcpy(addr[j],value);
                        }
                        break;
                    case INT:
                        *((int*)addr[j])=atoi(value);
                        break;
                    case LONG:
                        *((long*)addr[j])=atol(value);
                        break;
                    case BOOLEAN:
                        if (strchr("tTyY1", *value) != NULL) {
                            *((bool*)addr[j])=TRUE;
                        } else
                            if (strchr("fFnN0", *value) != NULL)  {
                                *((bool*)addr[j])=FALSE;
                            } else {
                                error("getbparam: %s=%s not bool\n",name,value);
                            }
                        break;
                }
            } else {
                fprintf(stdout, "Error in file %s: Tag '%s' %s %s.\n",
                        fname, name,
                        "not allowed or multiple defined...\n",
                        "look at saved parameter file which value was used");
                errorFlag=1;
            }
        } // ! while loop
        fclose(fd);
    } else {
        fprintf(stdout,"Parameter file %s not found.\n", fname);
        errorFlag=2;
        exit(0);
    }

    for(i=0;i<nt;i++) {
        if(*tag[i]) {
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
                    break;
            }
            errorFlag=1;
        }
    }

#undef DOUBLE
#undef STRING
#undef INT
#undef BOOLEAN
#undef MAXTAGS
#undef MAXCHARBUF
}
//E


int StartRun_Common(struct  cmdline_data* cmd, struct  global_data* gd)
{
    int ifile;

#ifdef ADDONS
#include "startrun_include_04.h"
#endif

    class_call_cballs(StartOutput(cmd, gd), errmsg, errmsg);

    setFilesDirs(cmd, gd);
    setFilesDirs_log(cmd, gd);
    strcpy(gd->mode,"w");
            if(!(gd->outlog=fopen(gd->logfilePath, gd->mode)))
                error("\nstart_Common: error opening file '%s' \n",
                      gd->logfilePath);

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
    if (scanopt(cmd->options, "statistics-histograms")) {
        statHistogram(cmd, gd);
        verb_print(cmd->verbose, "\n\tpre-processing: stopping...\n\n");
        exit(1);
    }
//E

//B In this section update computation of rSize and center-of-mass if necessary
//      so we have a common root size and c-of-m
    for (ifile=0; ifile<gd->ninfiles; ifile++) {
        if (!strnull(cmd->infile)) {
            class_call_cballs(infilefmt_string_to_int(gd->infilefmtname[ifile],
                            &gd->infilefmt_int), errmsg, errmsg);
            class_call_cballs(InputData(cmd, gd, gd->infilenames[ifile],ifile),
                            errmsg, errmsg);
            gd->model_comment = "input data file";
        } else {
            verb_print(cmd->verbose,"\nNo data catalog was given...");
            verb_print(cmd->verbose,"creating a test model...\n");
            TestData(cmd, gd);
            gd->input_comment = "no data file given";
        }
    }
//E
    for (ifile=0; ifile<gd->ninfiles; ifile++) {
        gd->bytes_tot += gd->nbodyTable[ifile]*sizeof(body);
        verb_print(cmd->verbose,
                   "\nAllocated %g MByte for particle storage (file %d).\n\n",
                   gd->nbodyTable[ifile]*sizeof(body)*INMB, ifile);
    }

    search_method_string_to_int(cmd->searchMethod, &gd->searchMethod_int);
    for (ifile=0; ifile<gd->ninfiles; ifile++) {
//B
        cellptr root;                           // Set it up a temporal root
        root = (cellptr) allocate(1 * sizeof(body));
        findRootCenter(cmd, gd,
                       bodytable[ifile], gd->nbodyTable[ifile], ifile, root);
        centerBodies(bodytable[ifile], gd->nbodyTable[ifile], ifile, root);
        findRootCenter(cmd, gd,
                           bodytable[ifile], gd->nbodyTable[ifile], ifile, root);
//E
        gd->rSizeTable[ifile] = 1.0;
        expandbox(cmd, gd, bodytable[ifile], gd->nbodyTable[ifile], ifile, root);
        free(root);
        if (cmd->rangeN > gd->rSizeTable[ifile])
            verb_print(cmd->verbose,
                       "\nstartrun_Common: warning! rangeN (%g) is greather than rSize (%g) of the system...\n",
                        cmd->rangeN, gd->rSizeTable[ifile]);
    }

    //B set output comment
    if (! strnull(cmd->outfile)) {
        gd->output_comment = "output data file was given";
    } else {
        gd->output_comment = "no output data file was given";
    }
    //E

//B Tree search:
    gd->Rcut = cmd->rangeN;                     // Maximum search radius
    gd->RcutSq = gd->Rcut*gd->Rcut;
//E
    if (cmd->useLogHist) {
        real rBin, rbinlog;
        if (cmd->rminHist==0) {
            //B rminHist = 0 not allowed when useLogHist is true
            //  therefore this segment does not ocurr...
            gd->deltaRV = dvector(1,cmd->sizeHistN);
            verb_log_print(cmd->verbose_log, gd->outlog, "\ndeltaRV:\n");
            int n;
            for (n=1; n<=cmd->sizeHistN; n++) {
                gd->deltaRV[n] =
                cmd->rangeN*rpow( 10.0,
                                ( (real)(n-cmd->sizeHistN) )/cmd->logHistBinsPD );
                verb_log_print(cmd->verbose_log, gd->outlog,
                            " %d %lg\n",n,gd->deltaRV[n]);
            }
            //E
        } else {
            gd->deltaR = rlog10(cmd->rangeN/cmd->rminHist)/cmd->sizeHistN;
            //B allocated after startrun_memoryAllocation
            //  deallocate before deallocate arrays in startrun_memoryAllocation
            gd->deltaRV = dvector(1,cmd->sizeHistN);
            gd->ddeltaRV = dvector(1,cmd->sizeHistN-1);
            //E
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

#ifdef CELLMETHOD
        MULVS(gd->cells, gd->Box, 1.0/gd->Rcut);    // Only needed for cellmethod
        AllocMem(gd->cellList, VProd (gd->cells)    // Only needed for cellmethod
                + cmd->nbody, INTEGER);
        gd->bytes_tot += (VProd(gd->cells)+cmd->nbody)*sizeof(INTEGER);
        verb_print(cmd->verbose,
                   "\n\nAllocated %g MByte for cells storage...\n",
                   (VProd(gd->cells)+cmd->nbody)*sizeof(INTEGER)*INMB);
#endif
        real Vol = 1.0;
        int k;
        DO_COORD(k)
            Vol = Vol*gd->Box[k];
        
        gd->i_deltaR = 1.0/gd->deltaR;              // This is gd->i_r_max
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
    real avgDistance = rpow(2.0*TWOPI/gd->nbodyTable[gd->iCatalogs[0]],1.0/2.0);
    verb_log_print(cmd->verbose_log, gd->outlog,
                   "Unit sphere (Takahasi): (S/N)^(1/2): %g\n",
                   avgDistance);
    verb_log_print(cmd->verbose_log, gd->outlog,
                   "and Nsmooth (Takahasi): (N*rs^2)/4: %g\n\n",
                   gd->nbodyTable[gd->iCatalogs[0]]*rpow(avgDistance,2.0)*0.25);


#define RADTOARCMIN   3437.74677
    //B cell size threshold computed for gd->iCatalogs[0] only...
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
            "\nstartrun_Common: i threshold, cell size and rsmooth = %d %e %e\n",
                       gd->irsmooth, rSizeTmp, gd->rsmooth[0]);
            verb_print(cmd->verbose,
                       "\t\t same in arcmin (useful for unit sphere)= %d %e %e\n",
                       gd->irsmooth, rSizeTmp*RADTOARCMIN,
                       gd->rsmooth[0]*RADTOARCMIN);
            break;
        }
    }
    //E
#undef RADTOARCMIN

    gd->deltaPhi = TWOPI/cmd->sizeHistPhi;
    // Nyquist frecuency = 1/(2 deltaPhi)
    verb_log_print(cmd->verbose_log, gd->outlog,
                   "\nNyquist frequency in phi bins = %g\n",0.5/gd->deltaPhi);

#ifdef ADDONS
#include "startrun_include_05.h"
#endif

    return SUCCESS;
}




//B Section of parameter check
local int CheckParameters(struct  cmdline_data* cmd, struct  global_data* gd)
{
// If it is necessary: an item in cmdline_defs.h must have an item here::

    //B Parameters related to the searching method
    if (cmd->useLogHist==FALSE && (strcmp(cmd->searchMethod,"balls-omp") == 0))
        error("CheckParameters: can´t have loghist false and balls-omp (%d %s)\n",
              cmd->useLogHist, cmd->searchMethod);
    if (cmd->computeTPCF) {
#ifndef USEGSL
        if (cmd->mChebyshev + 1 < 2 || (cmd->mChebyshev + 1)&(cmd->mChebyshev))
            error("CheckParameters: absurd value for mChebyshev + 1 (=%d)\n\t\t\tmust be positive and a power of 2\n", cmd->mChebyshev+1);
#else
        if (cmd->mChebyshev + 1 < 2)
            error("CheckParameters: absurd value for mChebyshev + 1 (=%d)\n\t\t\tmust be positive\n", cmd->mChebyshev+1);
#endif
    }
    if (cmd->mChebyshev + 1 > cmd->sizeHistPhi)
        error("CheckParameters: mChebyshev + 1 must be less than sizeHistPhi\n");
    if (gd->nsmooth[0] < 1)
        error("CheckParameters: absurd value for nsmooth\n");
    if (gd->rsmooth[0] < 0 || gd->rsmoothFlag==FALSE)
        error("CheckParameters: absurd value for rsmooth (%s)\n",cmd->rsmooth);
    if (cmd->theta < 0)
        error("CheckParameters: absurd value for theta\n");
    // We can get ride off one parameter if do sizeHistPhi = 2(mChebyshev + 1)
    //E

    //B Parameters to control the I/O file(s)
    // Input catalog parameters
    // Output parameters
    // Parameters to set a region in the sky, for example for Takahasi data set
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

    //B Parameters to control histograms and their output files
    if (cmd->useLogHist==TRUE && cmd->rminHist == 0)
        error("CheckParameters: can´t have useLogHist=true and rminHist=0 (%d %g)\n",
              cmd->useLogHist, cmd->rminHist);
    //
    if (cmd->sizeHistN < 2)
        error("CheckParameters: absurd value for sizeHistN\n");
    if (cmd->rangeN < 0)
        error("CheckParameters: absurd value for rangeN\n");
    if (cmd->rminHist < 0 || cmd->rminHist > cmd->rangeN)
        error("CheckParameters: absurd value for rminHist\n");
    if (cmd->sizeHistPhi < 2 || cmd->sizeHistPhi&(cmd->sizeHistPhi-1))
    error("CheckParameters: absurd value for sizeHistPhi\n\tmust be a power of 2\n");
    //
    //E
    
    //B Set of parameters needed to construct a test model
    if (cmd->nbody < 3) {
        error("CheckParameters: absurd value for nbody: %d\n", cmd->nbody);
    }
    if (cmd->lengthBox <= 0)
        error("CheckParameters: absurd value for lengthBox\n");

    //E

    //B Miscellaneous parameters
    if (cmd->stepState <= 0)
        error("CheckParameters: absurd value for stepState must be an integer > 0\n");
    if (cmd->verbose <= 0)
        error("CheckParameters: absurd value for stepState must be an integer > 0\n");
    if (cmd->verbose_log <= 0)
        error("CheckParameters: absurd value for stepState must be an integer > 0\n");
#ifdef OPENMPCODE
    if (cmd->numthreads <= 0)
        error("CheckParameters: absurd value for numberThreads must be an integer >= 0\n");
#endif
    //E

    //
    gd->bh86 = scanopt(cmd->options, "bh86");       // Barnes, J. & Hut, P.
                                                    //  1986. Nature 324,
                                                    //  446.
    gd->sw94 = scanopt(cmd->options, "sw94");       // Salmon, J.K. &
                                                    //  Warren, M.S. 1994.
                                                    //  J. Comp. Phys. 111,
                                                    //  136
    if (gd->bh86 && gd->sw94)
        error("CheckParameters: incompatible options bh86 and sw94\n");
    //

#ifdef ADDONS
#include "startrun_include_07.h"
#endif

    return SUCCESS;
}
//E


#define FMTT    "%-35s = %s\n"
#define FMTTS    "%-35s = \"%s\"\n"
#define FMTI    "%-35s = %d\n"
#define FMTIL    "%-35s = %ld\n"
#define FMTR	"%-35s = %g\n"

//B Section of parameter writing to a file
//local int PrintParameterFile(struct  cmdline_data *cmd, char *fname)
int PrintParameterFile(struct  cmdline_data *cmd, char *fname)
{
// Every item in cmdline_defs.h must have an item here::

    FILE *fdout;
    char buf[200];
    int  errorFlag=0;

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
            sprintf(buf,"%s/%s%s",cmd->rootDir,fname,"-usedvalues");
        } else {
            dp = (char*) malloc((strlen(fname)-ipos)*sizeof(char));
            strncpy(dp, fname + ipos, strlen(fname)-ipos);
            verb_print_q(3,cmd->verbose,
                         "PrintParameterFile: '/' counts %d pos %d and %s\n",
                         ndefault, ipos, dp);
            sprintf(buf,"%s/%s%s",cmd->rootDir,dp,"-usedvalues");
        }
//E

    if(!(fdout=fopen(buf,"w"))) {
        fprintf(stdout,"error opening file '%s' \n",buf);
        errorFlag=1;
    } else {
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
        // Parameters to set a region in the sky, for example for Takahasi data set
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
            verb_print(cmd->verbose, "\nPrintParamterFile: script: %s\n", cmd->script);
        fprintf(fdout,FMTT,"options",cmd->options);
        //E

#ifdef ADDONS
#include "startrun_include_08.h"
#endif

        fprintf(fdout,"\n\n");
    }
    fclose(fdout);

    if(errorFlag) {
        exit(0);
    }

    return SUCCESS;
}
//E

#undef FMTT
#undef FMTTS
#undef FMTI
#undef FMTR


local int random_init(struct  cmdline_data* cmd, struct  global_data* gd, int seed)
{
#ifdef USEGSL
    const gsl_rng_type * T;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    gd->r = gsl_rng_alloc (T);                      // Deallocate properly at EndRun...
    if (cmd->verbose>=3)
        verb_print(cmd->verbose, "\nrandom_init: gd->r and seed = %d %d\n",
                   *(gd->r), seed);
    gsl_rng_set(gd->r, seed);
    if (cmd->verbose>=3)
        verb_print(cmd->verbose, "\nrandom_init: gd->r and seed = %d %d\n",
                   *(gd->r), seed);
    r_gsl = gd->r;
    gd->bytes_tot += (1)*sizeof(gsl_rng_type);
#else
    idum = (long)seed;
    saveidum=idum;
    if (cmd->verbose>=3)
        verb_print(cmd->verbose, "\nrandom_init: idum and seed = %d %d\n",
                   idum, seed);
    xsrandom(idum);
    if (cmd->verbose>=3)
        verb_print(cmd->verbose, "\nrandom_init: idum = %d\n",idum);
#endif

    return SUCCESS;
}

global int startrun_memoryAllocation(struct  cmdline_data *cmd, 
                                     struct  global_data* gd)
{
    // Free allocated memory in reverse order as were allocated
    //  First is allocated above gsl structure gd->r

    INTEGER bytes_tot_local=0;
    //B PXD functions
    gd->rBins = dvector(1,cmd->sizeHistN);
    gd->histZetaMFlatten = dvector(1,cmd->sizeHistN*cmd->sizeHistN);
    //E PXD functions
    gd->histNN = dvector(1,cmd->sizeHistN);
    gd->histCF = dvector(1,cmd->sizeHistN);
    gd->histNNSub = dvector(1,cmd->sizeHistN);
    // 2pcf
    gd->histNNSubXi2pcf = dvector(1,cmd->sizeHistN);
    //B kappa Avg Rmin
    gd->histNNSubXi2pcftotal = dvector(1,cmd->sizeHistN);
    //E
    //
    gd->histNNN = dvector(1,cmd->sizeHistN);
    gd->histXi2pcf = dvector(1,cmd->sizeHistN);

    bytes_tot_local += 7*cmd->sizeHistN*sizeof(real);

    if (cmd->computeTPCF) {
        gd->histXi = dmatrix(1,cmd->mChebyshev+1,1,cmd->sizeHistN);
        gd->histXicos = dmatrix(1,cmd->mChebyshev+1,1,cmd->sizeHistN);
        gd->histXisin = dmatrix(1,cmd->mChebyshev+1,1,cmd->sizeHistN);
        bytes_tot_local += 3*(cmd->mChebyshev+1)*cmd->sizeHistN*sizeof(real);
        gd->histZetaM = dmatrix3D(1,cmd->mChebyshev+1,1,cmd->sizeHistN,
                                  1,cmd->sizeHistN);
        bytes_tot_local += 
                (cmd->mChebyshev+1)*cmd->sizeHistN*cmd->sizeHistN*sizeof(real);

        gd->histZetaMcos =
                dmatrix3D(1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
        gd->histZetaMsin =
                dmatrix3D(1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
        gd->histZetaMsincos =
                dmatrix3D(1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
        // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
        gd->histZetaMcossin =
                dmatrix3D(1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
        bytes_tot_local +=
                4*(cmd->mChebyshev+1)*cmd->sizeHistN*cmd->sizeHistN*sizeof(real);
#ifdef USEGSL
        //B Do not use anymore...
        /*
        gd->histXi_gsl = gsl_matrix_complex_calloc(cmd->mChebyshev+1,cmd->sizeHistN);
        bytes_tot_local += (cmd->mChebyshev+1)*cmd->sizeHistN*2*sizeof(real);

        histZetaMatrix = (mMatrix_ptr) allocate((cmd->mChebyshev+1) * sizeof(mMatrix));
        int m;
        for (m=0; m<=cmd->mChebyshev; m++){
            histZetaMatrix[m].histZetaM =
            gsl_matrix_complex_calloc(cmd->sizeHistN,cmd->sizeHistN);
        }
        bytes_tot_local += 
                (cmd->mChebyshev+1)*cmd->sizeHistN*cmd->sizeHistN*2*sizeof(real);
        */
        //E
#endif
        gd->histZetaGmRe =
                    dmatrix3D(1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
        gd->histZetaGmIm =
                    dmatrix3D(1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
        bytes_tot_local +=
                2*(cmd->mChebyshev+1)*cmd->sizeHistN*cmd->sizeHistN*sizeof(real);
    } // ! computeTPCF

    if (cmd->computeShearCF) {
        gd->histXitt = dvector(1,cmd->sizeHistN);
        gd->histXixx = dvector(1,cmd->sizeHistN);
        gd->histXitx = dvector(1,cmd->sizeHistN);
    }

#ifdef ADDONS
#include "startrun_include_10.h"                    // should be sync with
                                                    //  "cballsio_include_10.h"
#endif

    gd->bytes_tot += bytes_tot_local;
    verb_print(cmd->verbose,
    "\nstartrun_memoryAllocation: Allocated %g MByte for histograms storage.\n",
    bytes_tot_local*INMB);

    return SUCCESS;
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
#include "startrun_include_11.h"                    // See in this file
                                                    //  the last tag number used
#endif
}


#ifdef OPENMPCODE
int SetNumberThreads(struct  cmdline_data *cmd)
{
    omp_set_num_threads(cmd->numthreads);

    return SUCCESS;
}
#endif


//B Special routines to scan command line

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
            verb_log_print(cmd->verbose_log, gd->outlog,
                           "%s\n",gd->infilefmtname[nitems-1]);
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
    verb_log_print(cmd->verbose_log, gd->outlog, 
                   "\nSplitting string \"%s\" in tokens:\n",cmd->scanLevelMin);
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

#ifdef ADDONS
#include "startrun_include_12.h"
#endif

    return SUCCESS;
}

local int scaniOption(struct  cmdline_data* cmd, struct  global_data* gd,
                      string optionstr, int *option, int *noption,
                      int nfiles, int flag, string message)
{
    char *pch;
    char *poptionstr[30],  optiontmp[100];
    int i;

    if (cmd->verbose_log>=3)
    verb_log_print(cmd->verbose_log, gd->outlog, 
                   "\nProcessing '%s' option:\n", message);

    verb_log_print(cmd->verbose_log, gd->outlog,
                   "\nSplitting string \"%s=%s\" in tokens:\n",
                   message, optionstr);

    if (!strnull(optionstr)) {
        strcpy(optiontmp,optionstr);
        if (cmd->verbose_log>=3)
        verb_log_print(cmd->verbose_log, gd->outlog, 
                       "\nSplitting string \"%s\" in tokens:\n",
                       optiontmp);
        *noption=0;
        pch = strtok(optiontmp," ,");
        while (pch != NULL) {
            poptionstr[*noption] = (string) malloc(10);
            strcpy(poptionstr[*noption],pch);
            ++(*noption);
            if (cmd->verbose_log>=3)
            verb_log_print(cmd->verbose_log, gd->outlog,
                           "%s\n",poptionstr[*noption-1]);
            pch = strtok (NULL, " ,");
        }

        if (cmd->verbose_log>=3)
        verb_log_print(cmd->verbose_log, gd->outlog, 
                       "num. of tokens in option %s =%d\n", 
                       optionstr,*noption);

        if (flag == 0)
            if (*noption != nfiles)
                error("\nscanOption: noption = %d %s",
                      *noption,
                      "must be equal to number of infiles\n\n");
        if (*noption > MAXITEMS)
            error("\nscaniOption: noption = %d %s",
                  *noption,
                  "must be less than the maximum num. of lines\n\n");

        for (i=0; i<*noption; i++) {
            option[i]=atoi(poptionstr[i]);
            if (cmd->verbose_log>=3)
            verb_log_print(cmd->verbose_log, gd->outlog, 
                           "option: %d\n",option[i]);
        }
        if (cmd->verbose_log>=3)
        verb_log_print(cmd->verbose_log, gd->outlog, 
                       "\nnoptions, nfiles: %d %d\n",*noption,nfiles);
        if (flag == 1) {
            if (*noption > nfiles)
                error("\nscanOption: noption = %d %s",
                      *noption,
                      "must be less or equal to number of files\n\n");
            else {
                for (i=*noption; i<nfiles; i++) {
                    option[i]=option[i-1]+1;
                    if (cmd->verbose_log>=3)
                    verb_log_print(cmd->verbose_log, gd->outlog, 
                                   "option: %d\n",option[i]);
                }
            }
        }
    } else {
        for (i=0; i<nfiles; i++) {
            option[i]=1;
            if (cmd->verbose_log>=3)
            verb_log_print(cmd->verbose_log, gd->outlog, 
                           "option: %d\n",option[i]);
        }
    }

    return SUCCESS;
}

local int scanrOption(struct  cmdline_data* cmd, struct  global_data* gd,
                      string optionstr, double *option, int *noption,
    int nfiles, int flag, string message)
{
    char *pch;
    char *poptionstr[30],  optiontmp[100];
    int i;

    if (cmd->verbose_log>=3)
    verb_log_print(cmd->verbose_log, gd->outlog, 
                   "\nProcessing '%s' option:\n", message);

    verb_log_print(cmd->verbose_log, gd->outlog, 
                   "\nSplitting string \"%s\" in tokens:\n",message);

    if (!strnull(optionstr)) {
        strcpy(optiontmp,optionstr);
        if (cmd->verbose_log>=3)
        verb_log_print(cmd->verbose_log, gd->outlog, 
                       "\nSplitting string \"%s\" in tokens:\n",optiontmp);
        *noption=0;
        pch = strtok(optiontmp," ,");
        while (pch != NULL) {
            poptionstr[*noption] = (string) malloc(10);
            strcpy(poptionstr[*noption],pch);
            ++(*noption);
            if (cmd->verbose_log>=3)
            verb_log_print(cmd->verbose_log, gd->outlog,
                           "%s\n",poptionstr[*noption-1]);
            pch = strtok (NULL, " ,");
        }
        if (cmd->verbose_log>=3)
        verb_log_print(cmd->verbose_log, gd->outlog, 
                       "num. of tokens in option %s =%d\n", optionstr,*noption);

        if (flag == 0)
            if (*noption != nfiles)
            error("\nscanOption: noption = %d must be equal to number of files\n\n",
                      *noption);
        if (*noption > MAXITEMS)
    error("\nscanOption: noption = %d must be less than the maximum num. of lines\n\n",
                  *noption);

        for (i=0; i<*noption; i++) {
            option[i]=atof(poptionstr[i]);
            if (cmd->verbose_log>=3)
            verb_log_print(cmd->verbose_log, gd->outlog, "option: %g\n",option[i]);
        }

        if (cmd->verbose_log>=3)
        verb_log_print(cmd->verbose_log, gd->outlog, 
                       "\nnoptions, nfiles: %d %d\n",*noption,nfiles);
        if (flag == 1) {
            if (*noption > nfiles)
    error("\nscanOption: noption = %d must be less or equal to number of files\n\n",
                  *noption);
            else {
                for (i=*noption; i<nfiles; i++) {
                    option[i]=option[i-1]+1;
                    if (cmd->verbose_log>=3)
                    verb_log_print(cmd->verbose_log, gd->outlog, 
                                   "option: %g\n",option[i]);
                }
            }
        }
    } else {
        for (i=0; i<nfiles; i++) {
            option[i]=0;                            // Be aware of this values
            if (cmd->verbose_log>=3)
            verb_log_print(cmd->verbose_log, gd->outlog, "option: %d\n",option[i]);
        }
    }

    return SUCCESS;
}

//E Special routines to scan command line
