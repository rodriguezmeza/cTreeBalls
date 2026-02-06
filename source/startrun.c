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

//
// lines where there is a "//B socket:" string are places to include module files
//  that can be found in addons/addons_include folder
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

local int print_make_info(struct cmdline_data* cmd,
                     struct  global_data* gd);

#ifndef USEGSL
local long saveidum;
#endif

//B socket:
#ifdef ADDONS
#include "startrun_include_00.h"
#endif
//E

/*
 StartRun routine:

 To be called in main:
 StartRun(&cmd, &gd, argv[0], HEAD1, HEAD2, HEAD3);
 
 This routine is in charge of setting all global structures in order to
    the comutation process run smoothly with all parameters given
    by the user, set and checked.

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
    string routineName = "StartRun";
    double cpustart = CPUTIME;

    gd->headline0 = head0; gd->headline1 = head1;
    gd->headline2 = head2; gd->headline3 = head3;

    printf("\n%s\n%s: %s\n\t %s\n",
           gd->headline0, gd->headline1, gd->headline2, gd->headline3);
    printf("Version = %s\n", getversion());

    //B move all these to Startrun_Common... or make an appropriate change
    gd->cmd_allocated = FALSE;
    gd->stopflag = 0;
    gd->cputotalinout = 0.;
    gd->cputotal = 0.;
    gd->bytes_tot = 0;
    gd->sameposcount = 0;
    //E

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
    verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                "\n%s: Total allocated %g MByte storage so far.\n",
                        routineName, gd->bytes_tot*INMB);

//B If uncommented there will be a warning in the setup.py process
//#ifdef OPENMPCODE
    class_call_cballs(SetNumberThreads(cmd), errmsg, errmsg);
//#endif
//E
    gd->cputotalinout += CPUTIME - cpustart;
    verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                        "\n%s CPU time: %g %s\n",
                        routineName, CPUTIME - cpustart, PRNUNITOFTIMEUSED);

    return SUCCESS;
}

#else // ! CLASSLIB

#include "input.h"

int StartRun(struct  cmdline_data* cmd, struct  global_data* gd,
             string head0, string head1, string head2, string head3)
{
    string routineName = "StartRun";
    struct file_content fc;

    double cpustart = CPUTIME;
    
    gd->headline0 = head0; gd->headline1 = head1;
    gd->headline2 = head2; gd->headline3 = head3;
    printf("\n%s\n%s: %s\n\t %s\n",
           gd->headline0, gd->headline1, gd->headline2, gd->headline3);
    printf("Version = %s\n", getversion());

    //B move all these to Startrun_Common... or make an appropriate change
    gd->cmd_allocated = FALSE;
    gd->stopflag = 0;
    gd->cputotalinout = 0.;
    gd->cputotal = 0.;
    gd->bytes_tot = 0;
    gd->sameposcount = 0;
    //E

#ifdef GETPARAM
    cmd->paramfile = GetParam("paramfile");
    if (*(cmd->paramfile)=='-')
        error("bad parameter %s\n", cmd->paramfile);
    if (!strnull(cmd->paramfile)) {
        class_call_cballs(input_find_file(cmd, gd, cmd->paramfile, &fc, errmsg),
                          errmsg, errmsg);
        class_call_cballs(input_read_from_file(cmd, gd, &fc, errmsg),
                          errmsg, errmsg);
        class_call_cballs(parser_free(&fc), errmsg, errmsg);
    } else {
        startrun_cmdline(cmd, gd);
    }

#else
    // this segment must be checked!!! (used when GETPARMON = 0)
    class_call_cballs(input_find_file(cmd->ParameterFile, &fc, errmsg),
                      errmsg, errmsg);
    class_call_cballs(input_read_from_file(cmd, gd, &fc, errmsg), errmsg, errmsg);
    class_call_cballs(parser_free(&fc), errmsg, errmsg);
    //E
#endif

    if (!strnull(cmd->paramfile))
        class_call_cballs(StartRun_Common(cmd, gd), errmsg, errmsg);

//    if (gd->flagPrint==TRUE && gd->rootDirFlag==TRUE) {
#ifdef GETPARAM
        if (!strnull(cmd->paramfile))
            PrintParameterFile(cmd, gd, cmd->paramfile);
#else
        PrintParameterFile(cmd, gd, cmd->ParameterFile);
#endif
//    }

    gd->bytes_tot += sizeof(struct  global_data);
    gd->bytes_tot += sizeof(struct cmdline_data);
    verb_print(cmd->verbose,
               "\n%s: Total allocated %g MByte storage so far.\n",
               routineName, gd->bytes_tot*INMB);

#ifdef OPENMPCODE
    class_call_cballs(SetNumberThreads(cmd), errmsg, errmsg);
#endif

    gd->cputotalinout += CPUTIME - cpustart;
    verb_print(cmd->verbose, "\n%s CPU time: %g %s\n",
               routineName, CPUTIME - cpustart, PRNUNITOFTIMEUSED);

    return SUCCESS;
}
#endif // ! CLASSLIB


local int startrun_parameterfile(struct  cmdline_data* cmd,
                                 struct  global_data* gd)
{
#ifdef GETPARAM
	ReadParameterFile(cmd, gd, cmd->paramfile);
    ReadParametersCmdline_short(cmd, gd);

//B socket:
#ifdef ADDONS
#include "startrun_include_01.h"
#endif
//E

#else // ! GETPARAM
    ReadParameterFile(cmd, gd, cmd->ParameterFile);
#endif // ! GETPARAM

	StartRun_Common(cmd, gd);
//    if (gd->flagPrint==TRUE && gd->rootDirFlag==TRUE) {
#ifdef GETPARAM
        PrintParameterFile(cmd, gd, cmd->paramfile);
#else
        PrintParameterFile(cmd, gd, cmd->ParameterFile);
#endif
//    }

    return SUCCESS;
}


#ifdef GETPARAM
#define parameter_null	"parameters_null-cballs"

//B Section for reading parameters from the command line

local int startrun_cmdline(struct  cmdline_data* cmd, struct  global_data* gd)
{
	ReadParametersCmdline(cmd, gd);
	StartRun_Common(cmd, gd);
//    if (gd->flagPrint==TRUE && gd->rootDirFlag==TRUE) {
        PrintParameterFile(cmd, gd, parameter_null);
//    }

    return SUCCESS;
}

local void ReadParametersCmdline(struct  cmdline_data* cmd, 
                                 struct  global_data* gd)
{
// Every item in cmdline_defs.h must have an item here::

    //B Parameters related to the searching method
    cmd->searchMethod = GetParam("searchMethod");
    cmd->mChebyshev = GetiParam("mChebyshev");
    cmd->nsmooth = GetiParam("nsmooth");
    //E
    cmd->rsmooth = GetParam("rsmooth");
    cmd->theta = GetdParam("theta");
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
    cmd->preScript = GetParam("preScript");
    cmd->posScript = GetParam("posScript");
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

//B socket:
#ifdef ADDONS
#include "startrun_include_02.h"
#endif
//E
}

local void ReadParametersCmdline_short(struct  cmdline_data* cmd, struct  global_data* gd)
{
//B Here add parameters needed to be read after reading parameter file
    //B Miscellaneous parameters
//    cmd->preScript = GetParam("preScript");
//    cmd->posScript = GetParam("posScript");
    //E
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

    string routineName = "ReadParameterFile";
    FILE *fd;

  int  i,j,nt;
  int  id[MAXTAGS];
  void *addr[MAXTAGS];
  char tag[MAXTAGS][50];
  int  errorFlag=0;

    int input_verbose = 2;
    verb_print(input_verbose, "\nparsing input parameters...\n");

  nt=0;

    //B Parameters related to the searching method
    SPName(cmd->searchMethod,"searchMethod",MAXLENGTHOFSTRSCMD);
    IPName(cmd->mChebyshev,"mChebyshev");
    IPName(cmd->nsmooth,"nsmooth");
    //E
    SPName(cmd->rsmooth,"rsmooth",MAXLENGTHOFSTRSCMD);
    RPName(cmd->theta,"theta");
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
    SPName(cmd->preScript,"preScript",MAXLENGTHOFSTRSCMD);
    SPName(cmd->posScript,"posScript",MAXLENGTHOFSTRSCMD);
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

//B socket:
#ifdef ADDONS
#include "startrun_include_03.h"
#endif
//E

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
            phash=strchr(line,'%');
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
                        if (strcmp(name,"preScript") == 0){ // To remove both '"'
                            int index;
                            size_t slen;
                            slen = strlen(value);
                            //B
                            script1 = (char*) malloc(1*sizeof(char));
                            memcpy(script1,value,1);
                            script2 = strchr(script1, '"');
                            if (script2 == NULL)
                                error("preScript parameter needs enclosing script with \"\"!! (%s)\n\n",
                                      value);
                            free(script1);
                            //E
                            //B
                            script1 = (char*) malloc((slen-1)*sizeof(char));
                            memcpy(script1,value+1,slen-1);
                            script2 = strchr(script1, '"');
                            if (script2 == NULL)
                                error("preScript parameter needs enclosing script with \"\"!! (%s)\n\n",
                                      value);
                            free(script1);
                            //E
                            cmd->preScript = (char*) malloc((slen-2)*sizeof(char));
                            script1 = (char*) malloc(slen*sizeof(char));
                            memcpy(script1,value,slen);
                            script2 = strchr(script1, '"');
                            if (script2 == NULL)
                                error("preScript parameter needs enclosing script with \"\"!! (%s)\n\n",
                                      script1);
                            memcpy(cmd->preScript,script2+1,slen-2);
                            free(script1);
                        } else {
                            if (strcmp(name,"posScript")==0){// To remove both '"'
                                int index;
                                size_t slen;
                                slen = strlen(value);
                                //B
                                script1 = (char*) malloc(1*sizeof(char));
                                memcpy(script1,value,1);
                                script2 = strchr(script1, '"');
                                if (script2 == NULL)
                                    error("posScript parameter needs enclosing script with \"\"!! (%s)\n\n",
                                          value);
                                free(script1);
                                //E
                                //B
                                script1 = (char*) malloc((slen-1)*sizeof(char));
                                memcpy(script1,value+1,slen-1);
                                script2 = strchr(script1, '"');
                                if (script2 == NULL)
                                    error("posScript parameter needs enclosing script with \"\"!! (%s)\n\n",
                                          value);
                                free(script1);
                                //E
                                cmd->posScript=(char*) malloc((slen-2)*sizeof(char));
                                script1 = (char*) malloc(slen*sizeof(char));
                                memcpy(script1,value,slen);
                                script2 = strchr(script1, '"');
                                if (script2 == NULL)
                                    error("posScript parameter needs enclosing script with \"\"!! (%s)\n\n",
                                          script1);
                                memcpy(cmd->posScript,script2+1,slen-2);
                                free(script1);
                            } else
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
                fprintf(stdout, "\n%s: Error in file %s: Tag '%s' %s\n",
                        routineName, fname, name,
                        "not allowed or multiple defined...");
//                        "look at saved parameter file which value was used");
                errorFlag=1;
            }
        } // ! while loop
        fclose(fd);
    } else {
        fprintf(stdout,"Parameter file %s not found.\n", fname);
        errorFlag=2;
        exit(0);
    }

    if (errorFlag==1)
        error("%s: going out\n", routineName);

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
    string routineName = "StartRun_Common";
    int ifile;
    double cpustart;
    double cpustartMiddle;

    gd->cmd_allocated = TRUE;
    gd->histograms_allocated = FALSE;
    gd->random_allocated = FALSE;
    gd->gd_allocated = FALSE;
    gd->gd_allocated_2 = FALSE;
    gd->tree_allocated = FALSE;
    gd->bodytable_allocated = FALSE;

    if (strlen(cmd->rootDir)==0 || strnull(cmd->rootDir))
        gd->rootDirFlag = FALSE;
    else
        gd->rootDirFlag = TRUE;

    gd->flagPrint = TRUE;

    if (scanopt(cmd->options, "make-info"))
        print_make_info(cmd, gd);

//B socket:
#ifdef ADDONS
#include "startrun_include_04.h"
#endif
//E

//B correction 2025-05-03 :: look for edge-effects
#if defined(NMultipoles) && defined(NONORMHIST)
    if (scanopt(cmd->options, "patch-with-all")) {
        gd->pivotCount = 0;
    }
#endif
//E
    gd->pivotNumber = cmd->nbody;

    class_call_cballs(StartOutput(cmd, gd), errmsg, errmsg);

    debug_tracking_s("001", routineName);

    setFilesDirs(cmd, gd);
    debug_tracking("002");

    setFilesDirs_log(cmd, gd);
    debug_tracking("003");

    strcpy(gd->mode,"w");
    if (cmd->verbose_log>0) {               // gd->outlog is defined
        if(!(gd->outlog=fopen(gd->logfilePath, gd->mode)))
            error("\n%s: error opening file '%s' \n",
                  routineName, gd->logfilePath);
    }

     class_call_cballs(startrun_getParamsSpecial(cmd, gd), errmsg, errmsg);
    debug_tracking("004");

     class_call_cballs(random_init(cmd, gd, cmd->seed), errmsg, errmsg);

     class_call_cballs(CheckParameters(cmd, gd), errmsg, errmsg);

     class_call_cballs(startrun_memoryAllocation(cmd, gd), errmsg, errmsg);

    coordinate_string_to_int(cmd, gd);              // set coordTag
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
            "\n%s: coordTag: %d\n", routineName, gd->coordTag);

    debug_tracking("005");

//B Pre-processing necessary for reading data files:
    char buf[BUFFERSIZE];
    if (scanopt(cmd->options, "pre-processing")) {
        cpustartMiddle = CPUTIME;
        sprintf(buf,"%s",cmd->preScript);
        verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                               "\n%s: pre-processing: executing %s...\n",
                               routineName, cmd->preScript);
        system(buf);
        verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                               " done.\n");
        gd->cputotalinout += CPUTIME - cpustart;
        verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                               "%s: cpu time expended in this script %g\n\n",
                               routineName, CPUTIME - cpustartMiddle,
                               PRNUNITOFTIMEUSED);
        if (scanopt(cmd->options, "stop")) {
            verb_print_normal_info(cmd->verbose,
                                   cmd->verbose_log, gd->outlog,
                                   "\n\tMainLoop: stopping...\n\n");
            exit(1);
        }
    }
    if (scanopt(cmd->options, "statistics-histograms")) {
        statHistogram(cmd, gd);
        verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                               "\n\tpre-processing: stopping...\n\n");
        exit(1);
    }

    if (scanopt(cmd->options, "edge-corrections-from-files")) {
        computeEdgeCorrections(cmd, gd);
        verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                               "\n\tpre-processing: stopping...\n\n");
        exit(1);
    }
//E

    debug_tracking("006");

//B In this section update computation of rSize
//      and center-of-mass if necessary
//      so we have a common root size and c-of-m
    if (scanopt(cmd->options, "read-mask")) {
        if (gd->ninfiles < 2)
            error("\tevalHist:: read-mask ninfiles = %d is absurd\n", gd->ninfiles);
        ifile=0;
        class_call_cballs(infilefmt_string_to_int(gd->infilefmtname[ifile],
                    &gd->infilefmt_int), errmsg, errmsg);
        class_call_cballs(InputData(cmd, gd, gd->infilenames[ifile],ifile),
                    errmsg, errmsg);
        gd->model_comment = "input data file with mask";
        ifile=1;
        class_call_cballs(infilefmt_string_to_int(gd->infilefmtname[ifile],
                    &gd->infilefmt_int), errmsg, errmsg);
        class_call_cballs(InputData(cmd, gd, gd->infilenames[ifile],ifile),
                    errmsg, errmsg);
    } else {
        for (ifile=0; ifile<gd->ninfiles; ifile++) {
            debug_tracking_s("007",gd->infilenames[ifile]);

            if (!strnull(cmd->infile)) {
                debug_tracking_i("008",ifile);
                debug_tracking_s("009", gd->infilefmtname[ifile]);
                debug_tracking("010");
                class_call_cballs(
                            infilefmt_string_to_int(gd->infilefmtname[ifile],
                            &gd->infilefmt_int), errmsg, errmsg);
                debug_tracking_i("011",gd->infilefmt_int);
                class_call_cballs(InputData(cmd, gd,
                                            gd->infilenames[ifile],ifile),
                                            errmsg, errmsg);
                gd->model_comment = "input data file";
                debug_tracking("012");
            } else {
                verb_print_normal_info(cmd->verbose,
                                       cmd->verbose_log, gd->outlog,
                                       "\nNo data catalog was given...");
                verb_print_normal_info(cmd->verbose,
                                       cmd->verbose_log, gd->outlog,
                                       "creating a test model...\n");
                TestData(cmd, gd);
                gd->input_comment = "no data file given";
            }
        }
    }
    debug_tracking("013");

    //B consider moving below after computing rsize
    if (scanopt(cmd->options, "all-in-one")) {
        class_call_cballs(InputData_all_in_one(cmd, gd),
                          errmsg, errmsg);
    }
    //E
//E

    if (scanopt(cmd->options, "all-in-one")) {
        ifile = 0;
        gd->ninfiles = 1;
        verb_print_warning(cmd->verbose,
                           "%s: Warning! ninfile has been set to 1.\n\n",
                           routineName);
    }

    search_method_string_to_int(cmd->searchMethod, &gd->searchMethod_int);

    if (scanopt(cmd->options, "read-mask")) {
        ifile=0;
        cellptr root;                            // Set it up a temporal root
        root = (cellptr) allocate(1 * sizeof(body));
        FindRootCenter(cmd, gd,
                       bodytable[ifile], gd->nbodyTable[ifile], ifile, root);
        centerBodies(bodytable[ifile], gd->nbodyTable[ifile], ifile, root);
        FindRootCenter(cmd, gd,
                       bodytable[ifile], gd->nbodyTable[ifile], ifile, root);
        gd->rSizeTable[ifile] = 1.0;
        expandbox(cmd, gd, bodytable[ifile],
                  gd->nbodyTable[ifile], ifile, root);
        free(root);
        if (cmd->rangeN > gd->rSizeTable[ifile])
        verb_print_warning(cmd->verbose,
        "\n%s: warning! rangeN (%g) is greather than rSize (%g) of the system...\n",
                           routineName, cmd->rangeN, gd->rSizeTable[ifile]);
    } else {
        for (ifile=0; ifile<gd->ninfiles; ifile++) {
            debug_tracking("014");
            //B
            cellptr root;                        // Set it up a temporal root
            root = (cellptr) allocate(1 * sizeof(body));
            FindRootCenter(cmd, gd,
                           bodytable[ifile],
                           gd->nbodyTable[ifile], ifile, root);
            centerBodies(bodytable[ifile],
                         gd->nbodyTable[ifile], ifile, root);
            FindRootCenter(cmd, gd,
                           bodytable[ifile],
                           gd->nbodyTable[ifile], ifile, root);
            //E
            gd->rSizeTable[ifile] = 1.0;
            expandbox(cmd, gd, bodytable[ifile], gd->nbodyTable[ifile],
                      ifile, root);
            free(root);
            if (cmd->rangeN > gd->rSizeTable[ifile])
                verb_print_warning(cmd->verbose,
    "\n%s: warning! rangeN (%g) is greather than rSize (%g) of the system...\n",
                           routineName, cmd->rangeN, gd->rSizeTable[ifile]);
            debug_tracking("015");
        }
    }
    debug_tracking("016");

    //B set output comment
    if (! strnull(cmd->outfile)) {
        gd->output_comment = "output data file was given";
    } else {
        gd->output_comment = "no output data file was given";
    }
    //E
    debug_tracking("017");

//B Tree search:
    gd->Rcut = cmd->rangeN;                         // Maximum search radius
    //B correction 2025-05-03 :: look for edge-effects
    if (scanopt(cmd->options, "Rcut/theta")) {
        if (cmd->theta>0)
            gd->Rcut /= cmd->theta;                 // Maximum search radius
    }
    //E
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
            gd->bytes_tot += (cmd->sizeHistN)*sizeof(real);
            gd->bytes_tot += (cmd->sizeHistN-1)*sizeof(real);
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

    real Vol = 1.0;
    int k;
    DO_COORD(k)
        Vol = Vol*gd->Box[k];

    gd->i_deltaR = 1.0/gd->deltaR;                  // This is gd->i_r_max
                                                    //  change...
    verb_log_print(cmd->verbose_log, gd->outlog,
                       "\nRcut, deltaR: %g %g\n",gd->Rcut,gd->deltaR);
#if NDIM == 3
// CHECK!!! gd->Box[1]
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                           "%s: lbox: %g %g %g\n",
                            routineName, gd->Box[0],gd->Box[1],gd->Box[2]);
    DO_COORD(k)
        gd->Box[k] = cmd->lengthBox;
    verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                        "%s: Warning!! lbox size was changed!! %s\n",
                        routineName, "\n... make sure lengthBox was given right!!");
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                           "%s: lbox: %g %g %g\n",
                           routineName, gd->Box[0],gd->Box[1],gd->Box[2]);
#else
    verb_log_print(cmd->verbose_log, gd->outlog,
                   "lbox: %g %g\n",gd->Box[0],gd->Box[1]);
#endif
    verb_log_print(cmd->verbose_log, gd->outlog, "Box volume = %e\n",Vol);
    verb_log_print(cmd->verbose_log, gd->outlog,
                   "(V/N)^(1/3): %g\n\n",rpow(Vol/cmd->nbody,1.0/3.0));
    real avgDistance =
                rpow(2.0*TWOPI/gd->nbodyTable[gd->iCatalogs[0]],1.0/2.0);
    verb_log_print(cmd->verbose_log, gd->outlog,
                   "Unit sphere (Takahasi): (S/N)^(1/2): %g\n",
                   avgDistance);
    verb_log_print(cmd->verbose_log, gd->outlog,
                   "and Nsmooth (Takahasi): (N*rs^2)/4: %g\n\n",
                gd->nbodyTable[gd->iCatalogs[0]]*rpow(avgDistance,2.0)*0.25);

//B kappa Avg Rmin
#ifdef SMOOTHPIVOT
// 180*60/Pi
#define RADTOARCMIN   3437.74677
    //B cell size threshold computed for gd->iCatalogs[0] only...
    real rSizeTmp;
    int i, idepth=64;
    rSizeTmp = gd->rSizeTable[gd->iCatalogs[0]];
    for (i = 1; i <= idepth; i++) {
        rSizeTmp = rSizeTmp/2.0;
        verb_log_print(cmd->verbose_log, gd->outlog, "Cell size = %e\n",
                       rSizeTmp);
        if (rSizeTmp < rpow(2.0*TWOPI/gd->nbodyTable[gd->iCatalogs[0]],1.0/2.0)) {
            verb_log_print(cmd->verbose_log, gd->outlog,
                           "Cell size threshold = %d\n",i);
            gd->irsmooth = i;
            verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
            "\nstartrun_Common: i threshold, cell size and rsmooth = %d %e %e\n",
                        gd->irsmooth, rSizeTmp, gd->rsmooth[0]);
            verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                        "\t\t same in arcmin (useful for unit sphere)= %d %e %e\n",
                        gd->irsmooth, rSizeTmp*RADTOARCMIN,
                        gd->rsmooth[0]*RADTOARCMIN);
            break;
        }
    }
    //E
#undef RADTOARCMIN
#endif
//E

    gd->stepState = (INTEGER)(((real)cmd->stepState)/10.0);

    gd->deltaPhi = TWOPI/cmd->sizeHistPhi;
    // Nyquist frecuency = 1/(2 deltaPhi)
    verb_log_print(cmd->verbose_log, gd->outlog,
                   "\nNyquist frequency in phi bins = %g\n",
                   0.5/gd->deltaPhi);

//B socket:
#ifdef ADDONS
#include "startrun_include_05.h"
#endif
//E

    gd->gd_allocated = TRUE;
    gd->gd_allocated_2 = TRUE;
    gd->bodytable_allocated = TRUE;

    debug_tracking("018... final");

    return SUCCESS;
}


//B Section of parameter check
local int CheckParameters(struct  cmdline_data* cmd, struct  global_data* gd)
{
// If it is necessary: an item in cmdline_defs.h must have an item here::
    string routineName = "CheckParameters";

    debug_tracking_s("001", routineName);

    //B Parameters related to the searching method
    if (cmd->useLogHist==FALSE &&
        (strcmp(cmd->searchMethod,"balls-omp") == 0))
        error("%s: can´t have loghist false and balls-omp (%d %s)\n",
              routineName, cmd->useLogHist, cmd->searchMethod);
#ifdef TPCF
#ifndef USEGSL
        //  for recursivity needs that at least 3 multipoles be evaluated
        if (scanopt(cmd->options, "out-HistZetaG")) {
            if (cmd->mChebyshev + 1 < 3
                || (cmd->mChebyshev + 1)&(cmd->mChebyshev)) {
                verb_print(cmd->verbose,
                    "\n%s: using option out-HistZetaG...\n", routineName);
                error("%s: %s (=%d)\n\t\t\t%s\n",
                      routineName, "absurd value for mChebyshev + 1",
                      cmd->mChebyshev+1, "must be positive and a power of 2");
            }
        } else {
            if (cmd->mChebyshev + 1 < 3)
            error("%s: %s (=%d)\n\t\t\tmust be positive\n",
                  "absurd value for mChebyshev + 1",
                  routineName, cmd->mChebyshev+1);
        }
#else
        //  for recursivity needs that at least 3 multipoles be evaluated
        if (cmd->mChebyshev + 1 < 3)
            error("CheckParameters: absurd value for mChebyshev + 1 (=%d)\n\t\t\tmust be positive\n", cmd->mChebyshev+1);
#endif
#endif
    if (cmd->mChebyshev + 1 > cmd->sizeHistPhi)
        error("CheckParameters: mChebyshev + 1 must be less than sizeHistPhi\n");
    if (cmd->nsmooth < 1)
        error("%s: absurd value for nsmooth: %d\n",
              routineName, cmd->nsmooth);
    debug_tracking("002");

    if (!strnull(cmd->rsmooth)) {
        if (gd->rsmooth[0] < 0.0 || gd->rsmoothFlag==FALSE)
            error("CheckParameters: absurd value for rsmooth (%s, %g, %d)\n",
                  cmd->rsmooth, gd->rsmooth[0], gd->rsmoothFlag);
    } else {
        gd->rsmooth[0] = 0.0;
    }
    if (cmd->theta < 0)
        error("CheckParameters: absurd value for theta\n");
    // We can get ride off one parameter if do sizeHistPhi = 2(mChebyshev+1)
    //E

    debug_tracking("003");

    //B Parameters to control the I/O file(s)
    // Input catalog parameters
    // Output parameters
    // Parameters to set a region in the sky, for example for Takahashi data
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

    debug_tracking("004");

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
        error("%s: absurd value for nbody: %d\n", routineName, cmd->nbody);
    }
    if (cmd->lengthBox <= 0)
        error("%s: absurd value for lengthBox\n", routineName);

    //E

    debug_tracking("005");

    //B Miscellaneous parameters
    if (cmd->stepState <= 0)
        error("%s: absurd value for stepState must be an integer > 0\n",
              routineName);
    if (cmd->verbose < 0)
        error("%s: absurd value for verbose must be an integer >= 0\n",
              routineName);
    if (cmd->verbose_log < 0)
        error("%s: absurd value for verbose_log must be an integer >= 0\n",
              routineName);
#ifdef OPENMPCODE
    if (cmd->numthreads <= 0)
        error("%s: absurd value for numberThreads must be an integer >= 0\n",
              routineName);
#endif
    //E
    debug_tracking("006");

    //
    gd->bh86 = scanopt(cmd->options, "bh86");       // Barnes, J. & Hut, P.
                                                    //  1986. Nature 324,
                                                    //  446.
    gd->sw94 = scanopt(cmd->options, "sw94");       // Salmon, J.K. &
                                                    //  Warren, M.S. 1994.
                                                    //  J. Comp. Phys. 111,
                                                    //  136
    debug_tracking("007");

    if (gd->bh86 && gd->sw94)
        error("%s: incompatible options bh86 and sw94\n", routineName);
    //

    debug_tracking("008... follows OCTREESMOOTHING and IOLIB...");

//B socket:
#ifdef ADDONS
#include "startrun_include_07.h"
#endif
//E
    debug_tracking("013... final");

    return SUCCESS;
}
//E


#define FMTT    "%-35s = %s\n"
#define FMTTS    "%-35s = \"%s\"\n"
#define FMTI    "%-35s = %d\n"
#define FMTIL    "%-35s = %ld\n"
#define FMTR	"%-35s = %g\n"

//B Section of parameter writing to a file
int PrintParameterFile(struct  cmdline_data *cmd,
                       struct  global_data* gd, char *fname)
{
// Every item in cmdline_defs.h must have an item here::
    string routineName = "PrintParameterFile";

    FILE *fdout;
    char buf[200];
    int  errorFlag=0;

    debug_tracking_s("001", routineName);

    if (gd->flagPrint==TRUE && gd->rootDirFlag==TRUE) {
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
            fprintf(fdout,FMTI,"nsmooth",cmd->nsmooth);
            fprintf(fdout,FMTT,"rsmooth",cmd->rsmooth);
            fprintf(fdout,FMTR,"theta",cmd->theta);
            fprintf(fdout,FMTT,"usePeriodic",
                    cmd->usePeriodic ? "true" : "false");
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
            //  for example for Takahashi data set
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
            fprintf(fdout,FMTTS,"preScript",cmd->preScript);
            fprintf(fdout,FMTTS,"posScript",cmd->posScript);
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
            if (cmd->verbose>=VERBOSEDEBUGINFO) {
                verb_print(cmd->verbose, "\n%s: PrintParamterFile: preScript: %s\n",
                           routineName, cmd->preScript);
                verb_print(cmd->verbose, "\n%s: PrintParamterFile: posScript: %s\n",
                           routineName, cmd->posScript);
            }
            fprintf(fdout,FMTT,"options",cmd->options);
            //E
            
            //B socket:
#ifdef ADDONS
#include "startrun_include_08.h"
#endif
            //E
            
            fprintf(fdout,"\n\n");
        }
        fclose(fdout);
        
        if(errorFlag) {
            exit(0);
        }
        
        if (ndefault != 0)
            free(dp);
        
    } // ! gd->flagPrint==TRUE && gd->rootDirFlag==TRUE

    debug_tracking("002... final");


    return SUCCESS;
}
//E

#undef FMTT
#undef FMTTS
#undef FMTI
#undef FMTR


local int random_init(struct  cmdline_data* cmd, struct  global_data* gd, int seed)
{
    string routineName = "random_init";

    gd->random_allocated = FALSE;

#ifdef USEGSL
    const gsl_rng_type * T;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    gd->r = gsl_rng_alloc (T);                      // Deallocate at EndRun...
    verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                          "\n%s: gd->r and seed = %d %d\n",
                          routineName, *(gd->r), seed);
    gsl_rng_set(gd->r, seed);
    verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                          "\n%s: gd->r and seed = %d %d\n",
                          routineName, *(gd->r), seed);
    r_gsl = gd->r;
    gd->bytes_tot += (1)*sizeof(gsl_rng_type);

    gd->random_allocated = TRUE;

#else
    idum = (long)seed;
    saveidum=idum;
    verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                          "\n%s: idum and seed = %d %d\n",
                          routineName, idum, seed);
    xsrandom(idum);
    verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                          "\n%s: idum = %d\n",
                          routineName, idum);
#endif

    return SUCCESS;
}

global int startrun_memoryAllocation(struct  cmdline_data *cmd, 
                                     struct  global_data* gd)
{
    string routineName = "startrun_memoryAllocation";
    // Free allocated memory in reverse order as were allocated
    //  First is allocated above gsl structure gd->r

    INTEGER bytes_tot_local=0;
    //B PXD functions
#ifdef PXD
    gd->vecPXD = dvector(1,cmd->sizeHistN);
    gd->matPXD = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
    bytes_tot_local += cmd->sizeHistN*sizeof(real);
    bytes_tot_local += cmd->sizeHistN*cmd->sizeHistN*sizeof(real);
    gd->rBins = dvector(1,cmd->sizeHistN);
    gd->histZetaMFlatten = dvector(1,cmd->sizeHistN*cmd->sizeHistN);
#endif
    //E PXD functions
    gd->histNN = dvector(1,cmd->sizeHistN);
    gd->histCF = dvector(1,cmd->sizeHistN);
    gd->histNNSub = dvector(1,cmd->sizeHistN);
    // 2pcf
    gd->histNNSubXi2pcf = dvector(1,cmd->sizeHistN);
#ifdef SMOOTHPIVOT
    gd->histNNSubXi2pcftotal = dvector(1,cmd->sizeHistN);
#endif
    //
    gd->histNNN = dvector(1,cmd->sizeHistN);
    gd->histXi2pcf = dvector(1,cmd->sizeHistN);

    bytes_tot_local += 8*cmd->sizeHistN*sizeof(real);
    bytes_tot_local += cmd->sizeHistN*cmd->sizeHistN*sizeof(real);

#ifdef TPCF
        gd->histXicos = dmatrix(1,cmd->mChebyshev+1,1,cmd->sizeHistN);
        gd->histXisin = dmatrix(1,cmd->mChebyshev+1,1,cmd->sizeHistN);
        bytes_tot_local += 2*(cmd->mChebyshev+1)*cmd->sizeHistN*sizeof(real);
        gd->histZetaM = dmatrix3D(1,cmd->mChebyshev+1,1,cmd->sizeHistN,
                                  1,cmd->sizeHistN);
        bytes_tot_local += 
                (cmd->mChebyshev+1)*cmd->sizeHistN*cmd->sizeHistN*sizeof(real);

        gd->histZetaMcos =
                dmatrix3D(1,cmd->mChebyshev+1,1,
                          cmd->sizeHistN,1,cmd->sizeHistN);
        gd->histZetaMsin =
                dmatrix3D(1,cmd->mChebyshev+1,1,
                          cmd->sizeHistN,1,cmd->sizeHistN);
        gd->histZetaMsincos =
                dmatrix3D(1,cmd->mChebyshev+1,1,
                          cmd->sizeHistN,1,cmd->sizeHistN);
        // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
        gd->histZetaMcossin =
                dmatrix3D(1,cmd->mChebyshev+1,1,
                          cmd->sizeHistN,1,cmd->sizeHistN);
        bytes_tot_local +=
            4*(cmd->mChebyshev+1)*cmd->sizeHistN*cmd->sizeHistN*sizeof(real);
        gd->histZetaGmRe =
                    dmatrix3D(1,cmd->mChebyshev+1,1,
                              cmd->sizeHistN,1,cmd->sizeHistN);
        gd->histZetaGmIm =
                    dmatrix3D(1,cmd->mChebyshev+1,1,
                              cmd->sizeHistN,1,cmd->sizeHistN);
        bytes_tot_local +=
            2*(cmd->mChebyshev+1)*cmd->sizeHistN*cmd->sizeHistN*sizeof(real);
#endif

//B socket:
#ifdef ADDONS
    // this is empty and can be remove these 3 lines
#include "startrun_include_10.h"                    // should be sync with
                                                //  "cballsio_include_10.h"
#endif
//E

    //B this were in startrun_include_10.h -> startrun_octree_kkk_omp_10.h above...
    // problems with OCTREEKKKOMPON = 0
    // 2pcf
    gd->histNNSubN2pcf = dvector(1,cmd->sizeHistN);
#ifdef SMOOTHPIVOT
    gd->histNNSubN2pcftotal = dvector(1,cmd->sizeHistN);
#endif
    gd->histN2pcf = dvector(1,cmd->sizeHistN);
    bytes_tot_local += 3*cmd->sizeHistN*sizeof(real);
    //E

    gd->bytes_tot += bytes_tot_local;
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                           "\n%s: Allocated %g MByte for histograms storage.\n",
                           routineName, bytes_tot_local*INMB);

    gd->histograms_allocated = TRUE;

    return SUCCESS;
}

local void search_method_string_to_int(string method_str,int *method_int)
{
// Every search method must have an item here::
    *method_int=-1;
    if (strnull(method_str))
                *method_int = SEARCHNULL;
    if (strcmp(method_str,"octree-sincos-omp") == 0)
                *method_int = OCTREESINCOSOMPMETHOD;

//B socket:
#ifdef ADDONS
#include "startrun_include_11.h"                    // See in this file
                                                    //  the last tag number used
#endif
//E
}


#ifdef OPENMPCODE
int SetNumberThreads(struct  cmdline_data *cmd)
{
    omp_set_num_threads(cmd->numthreads);

    return SUCCESS;
}
#endif


//B Special routines to scan command line

local int startrun_getParamsSpecial(struct  cmdline_data* cmd,
                                    struct  global_data* gd)
{
    string routineName = "startrun_getParamsSpecial";
    char *pch;
    int nitems, ndummy=1;
    char inputnametmp[MAXLENGTHOFSTRSCMD];
    int i;
    size_t slen;

    debug_tracking_s("001", routineName);

    if (strnull(cmd->infile)) {
        verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                        "\n%s: no inputfile was given, making data ...\n",
                        routineName);
        gd->ninfiles=1;                              // To test data...
    } else {
        strcpy(inputnametmp,cmd->infile);
        verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                              "\n%s: Splitting string \"%s\" in tokens:\n",
                              routineName, inputnametmp);
        gd->ninfiles=0;
        pch = strtok(inputnametmp," ,");
        while (pch != NULL) {
            slen = strlen(pch);
            snprintf(gd->infilenames[gd->ninfiles],
                     slen+1, "%s", pch);
            ++gd->ninfiles;
            verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                                  "%s: %s\n",
                            routineName, gd->infilenames[gd->ninfiles-1]);
            pch = strtok (NULL, " ,");
        }
        verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                              "%s: num. of files in infile %s =%d\n",
                            routineName, cmd->infile, gd->ninfiles);
    }
    debug_tracking("002");

    if (!strnull(cmd->infilefmt)) {
        strcpy(inputnametmp,cmd->infilefmt);
        verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                              "\n%s: Splitting string \"%s\" in tokens:\n",
                            routineName, inputnametmp);
        nitems=0;
        pch = strtok(inputnametmp," ,");
        while (pch != NULL) {
            slen = strlen(pch);
            debug_tracking_s("003: pch",pch);
            snprintf(gd->infilefmtname[nitems],
                     slen+1, "%s", pch);
            ++nitems;
            verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                                  "%s: %s\n",
                                routineName, gd->infilefmtname[nitems-1]);
            pch = strtok (NULL, " ,");
            debug_tracking_i("004: not null infilefmt", nitems);
            debug_tracking_s("...", gd->infilefmtname[nitems-1]);
        }
        verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                              "%s: num. of items in infilefmt %s =%d\n",
                            routineName, cmd->infilefmt, nitems);
        if (nitems != gd->ninfiles)
    error("\nstartrun_Common: nitems must be equal to number of files\n\n");
        debug_tracking("005: not null infilefmt");

    } else {
        nitems=1;
        pch = "columns-ascii";
        strcpy(gd->infilefmtname[nitems],pch);
        verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                              "%s: warning!! found null infilefmt:\n",
                            routineName);
        verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                              "%s: %s\n",
                            routineName, gd->infilefmtname[nitems-1]);
        verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                              "%s: num. of items in infilefmt %s =%d\n",
                            routineName, cmd->infilefmt, nitems);
        if (nitems != gd->ninfiles)
            error("\nstartrun_Common: nitems must be equal to number of files\n\n");
        debug_tracking("006: null infilefmt");
    }

    debug_tracking("007");

    scanrOption(cmd, gd, cmd->rsmooth, gd->rsmooth, &nitems, ndummy,
                2, "rsmooth");
    gd->rsmoothFlag = TRUE;
    if (nitems!=1 && !strnull(cmd->rsmooth)) {
        verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                              "%s: option: rsmooth=%s is not valid... going out\n",
                              routineName, cmd->rsmooth);
        gd->rsmoothFlag = FALSE;
    }

    debug_tracking("008");

    scaniOption(cmd, gd, cmd->iCatalogs, gd->iCatalogs,
                &nitems, gd->ninfiles, 0, "iCatalogs");
    if (gd->ninfiles==1) {
        (gd->iCatalogs[0]) = 0;
        if (cmd->verbose_log>=2)
        verb_log_print(cmd->verbose_log, gd->outlog,
                       "%s: option: iCatalogs final values: %d\n",
                       routineName, gd->iCatalogs[0]);
        (gd->iCatalogs[1]) = 0;
        if (cmd->verbose_log>=2)
        verb_log_print(cmd->verbose_log, gd->outlog,
                       "%s: option: iCatalogs final values: %d\n",
                       routineName, gd->iCatalogs[1]);
    } else {
        for (i=0; i<gd->ninfiles; i++) {
            (gd->iCatalogs[i])--;
            if (cmd->verbose_log>=2)
                verb_log_print(cmd->verbose_log, gd->outlog,
                               "%s:option: iCatalogs final values: %d\n",
                               routineName, gd->iCatalogs[i]);
        }
    }

    debug_tracking("009");

//B socket:
#ifdef ADDONS
#include "startrun_include_12.h"
#endif
//E
    debug_tracking("010... final");

    return SUCCESS;
}

local int scaniOption(struct  cmdline_data* cmd, struct  global_data* gd,
                      string optionstr, int *option, int *noption,
                      int nfiles, int flag, string message)
{
    string routineName = "scaniOption";
    char *pch;
    char *poptionstr[30],  optiontmp[100];
    int i;

    verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                           "\n%s: Processing '%s' option:\n",
                        routineName, message);

    verb_log_print(cmd->verbose_log, gd->outlog,
                           "%s: Splitting string \"%s=%s\" in tokens:\n",
                        routineName, message, optionstr);

    if (!strnull(optionstr)) {
        strcpy(optiontmp,optionstr);
        verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                               "%s: Splitting string \"%s\" in tokens:\n",
                            routineName, optiontmp);
        *noption=0;
        pch = strtok(optiontmp," ,");
        while (pch != NULL) {
            poptionstr[*noption] = (string) malloc(10);
            strcpy(poptionstr[*noption],pch);
            ++(*noption);
            verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                                  "%s: %s\n",
                                routineName, poptionstr[*noption-1]);
            pch = strtok (NULL, " ,");
        }

        verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                               "%s: num. of tokens in option %s =%d\n",
                            routineName, optionstr, *noption);

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
            verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                                   "%s: option: %d\n",
                                routineName, option[i]);
        }
        verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                               "%s: noptions, nfiles: %d %d\n\n",
                            routineName, *noption,nfiles);
        if (flag == 1) {
            if (*noption > nfiles)
                error("\nscanOption: noption = %d %s",
                      *noption,
                      "must be less or equal to number of files\n\n");
            else {
                for (i=*noption; i<nfiles; i++) {
                    option[i]=option[i-1]+1;
                    verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                                           "%s: option: %d\n",
                                        routineName, option[i]);
                }
            }
        }
    } else {
        //B ask a question if optionstr is columns then fix columns to default values..
        if (scanopt(message, "columns")) {
            gd->columns[0] = 1;
            gd->columns[1] = 2;
            gd->columns[2] = 3;
            gd->columns[3] = 4;
            gd->columns[4] = 5;
            gd->columns[5] = 6;
            gd->columns[6] = 7;
        } else {
            for (i=0; i<nfiles; i++) {
                option[i]=1;
                verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                                       "%s: option: %d\n",
                                       routineName, option[i]);
            }
        }
        //E
    }

    return SUCCESS;
}

local int scanrOption(struct  cmdline_data* cmd, struct  global_data* gd,
                      string optionstr, double *option, int *noption,
    int nfiles, int flag, string message)
{
    string routineName = "scanrOption";

    char *pch;
    char *poptionstr[30],  optiontmp[100];
    int i;

    verb_log_print(cmd->verbose_log, gd->outlog,
                           "\n%s: Processing '%s' option:\n",
                        routineName, message);

    verb_log_print(cmd->verbose_log, gd->outlog,
                           "%s: Splitting string \"%s\" in tokens:\n",
                        routineName, message);

    if (!strnull(optionstr)) {
        strcpy(optiontmp,optionstr);
        verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                               "\n%s: Splitting string \"%s\" in tokens:\n",
                            routineName, optiontmp);
        *noption=0;
        pch = strtok(optiontmp," ,");
        while (pch != NULL) {
            poptionstr[*noption] = (string) malloc(MAXLENGTHOFREAL);
            strcpy(poptionstr[*noption],pch);
            ++(*noption);
            verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                                  "%s: %s\n",
                                routineName, poptionstr[*noption-1]);
            pch = strtok (NULL, " ,");
        }
        verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                              "%s: num. of tokens in option %s =%d\n",
                            routineName, optionstr, *noption);

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

        verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                              "\n%s: noptions, nfiles: %d %d\n",
                            routineName, *noption, nfiles);
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

        for (i=0; i<*noption; i++) {
            free(poptionstr[*noption]);
        }

    } else {
        for (i=0; i<nfiles; i++) {
            option[i]=0.0;                          // Be aware of this values
            verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                                   "%s: option: %d\n",
                                routineName, option[i]);
        }
    }

    return SUCCESS;
}

//E Special routines to scan command line

local int print_make_info(struct cmdline_data* cmd,
                                  struct  global_data* gd)
{
    verb_print(cmd->verbose,
               "\nprint_make_info:\n");

#ifdef TWODIMCODE
    verb_print(cmd->verbose, "TWODIMCODE\n");
#endif
#ifdef THREEDIMCODE
    verb_print(cmd->verbose, "THREEDIMCODE\n");
#endif

#ifdef OPENMPCODE
    verb_print(cmd->verbose, "using OpenMP\n");
#endif

#ifdef SINGLEP
    verb_print(cmd->verbose, "SINGLEP\n");
#endif

#ifdef SMOOTHPIVOT
        verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                            "with option smooth-pivot...\n");
#endif

#ifdef BALLS4SCANLEV
    verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                        "with BALLS4SCANLEV... \n");
#endif

#ifdef TWOPCF
    verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                        "with 2pcf convergence computation... \n");
#endif

#ifdef THREEPCFCONVERGENCE
    verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                        "with 3pcf convergence computation... \n");
#endif

#ifdef THREEPCFSHEAR
    verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                        "with 3pcf shear computation... \n");
#endif

#ifdef KappaAvgON
    verb_print(cmd->verbose, "KappaAvgON\n");
#endif

#ifdef LONGINT
    verb_print(cmd->verbose, "LONGINT\n");
#endif

#ifdef PTOPIVOTROTATION
    verb_print(cmd->verbose, "PTOPIVOTROTATION\n");
#endif

#ifdef DEBUG
    verb_print(cmd->verbose, "DEBUG\n");
#endif

#ifdef GETPARAM
    verb_print(cmd->verbose, "GETPARAM\n");
#endif

#ifdef USEGSL
#ifdef GSLINTERNAL
    verb_print(cmd->verbose, "using internal GSL\n");
#else
    verb_print(cmd->verbose, "using GSL\n");
#endif
#endif

#ifdef ADDONS
    verb_print(cmd->verbose, "with ADDONS\n");
#endif

#ifdef BALLS
    verb_print(cmd->verbose, "with BALLS\n");
#endif

#ifdef OCTREESMOOTHING
    verb_print(cmd->verbose, "with OCTREESMOOTHING\n");
#endif

#ifdef SINCOS
    verb_print(cmd->verbose, "with SINCOS\n");
#endif

#ifdef TREENODEALLBODIES
    verb_print(cmd->verbose, "with TREENODEALLBODIES\n");
#endif

#ifdef TREENODEBALLS4
    verb_print(cmd->verbose, "with TREENODEBALLS4\n");
#endif

#ifdef PIVOTEXTERNAL
    verb_print(cmd->verbose, "with PIVOTEXTERNAL\n");
#endif

#ifdef GADGETIO
    verb_print(cmd->verbose, "with GADGETIO\n");
#endif

#ifdef CLASSLIB
    verb_print(cmd->verbose, "with CLASSLIB\n");
#endif

#ifdef PXD
    verb_print(cmd->verbose, "with PXD\n");
#endif

#ifdef IOLIB
    verb_print(cmd->verbose, "with IOLIB\n");
#endif

#ifdef CFITSIO
    verb_print(cmd->verbose, "with CFITSIO\n");
#endif

#ifdef KDTREEOMP
    verb_print(cmd->verbose, "with KDTREEOMP\n");
#endif

#ifdef OCTREEKKKOMP
    verb_print(cmd->verbose, "with OCTREEKKKOMP\n");
#endif

#ifdef OCTREEGGGOMP
    verb_print(cmd->verbose, "with OCTREEGGGOMP\n");
#endif

#ifdef NMultipoles
    verb_print(cmd->verbose, "with NMultipoles\n");
#else
    verb_print(cmd->verbose, "without NMultipoles\n");
#endif
#ifdef NONORMHIST
    verb_print(cmd->verbose, "with NONORMHIST\n");
#else
    verb_print(cmd->verbose, "without NONORMHIST\n");
#endif

#ifdef POLARAXIS
    verb_print(cmd->verbose, "with POLARAXIS\n");
#endif

#ifdef NOLIMBER
    verb_print(cmd->verbose, "with NOLIMBER\n");
#endif

#ifdef ADDPIVOTNEIGHBOURS
        verb_print(cmd->verbose, "with ADDPIVOTNEIGHBOURS\n");
#endif

#ifdef OVERCOUNTING
        verb_print(cmd->verbose, "with OVERCOUNTING\n");
#endif

#ifdef SAVERESTORE
        verb_print(cmd->verbose, "with SAVERESTORE\n");
#endif

#ifdef DIRECTMETHOD
        verb_print(cmd->verbose, "with DIRECTMETHOD\n");
#endif

#ifdef DIRECTMETHODSIMPLE
        verb_print(cmd->verbose, "with DIRECTMETHODSIMPLE\n");
#endif

    return SUCCESS;
}

