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
#include "cballs_make_info.h"
//B
#include <limits.h>
#include <errno.h>

#ifdef CLASSLIB
#define cBALLS_FAIL(cmd, ...)                                      \
    do {                                                          \
        snprintf((cmd)->error_message, _ERRORMSGSIZE_, __VA_ARGS__); \
        return FAILURE;                                           \
    } while (0)
#else
#define cBALLS_FAIL(cmd, ...) error(__VA_ARGS__)
#endif
//E

//B
local int ReadParameterFile(struct  cmdline_data*,
                             struct  global_data*, char *);
//E
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
local int print_options(struct cmdline_data* cmd,
                        struct  global_data* gd);
local int print_search_methods(struct cmdline_data* cmd,
                               struct global_data* gd);

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
    string routineName = "StartRun (no CLASSLIB)";
    double cpustart = CPUTIME;

    gd->headline0 = head0; gd->headline1 = head1;
    gd->headline2 = head2; gd->headline3 = head3;

    printf("\n%s\n%s: %s\n\t %s\n",
           gd->headline0, gd->headline1, gd->headline2, gd->headline3);
#ifdef GETPARAM
    printf("Version = %s\n", getversion());
#endif

    //B move all these to Startrun_Common... or make an appropriate change
    gd->startrun_cputime = TRUE;
    gd->cmd_allocated = FALSE;
    gd->stopflag = 0;
    gd->cputotalinout = 0.;
    gd->cputotal = 0.;
    gd->bytes_tot = 0;
    gd->sameposcount = 0;
    //E

#ifdef GETPARAM
    cmd->paramfile = GetParam("paramfile");
    if (*(cmd->paramfile) == '-') {
        cBALLS_FAIL(cmd, "bad parameter %s\n", cmd->paramfile);
    }

    if (!strnull(cmd->paramfile)) {
        if (startrun_parameterfile(cmd, gd) == FAILURE)
            return FAILURE;
    } else {
        if (startrun_cmdline(cmd, gd) == FAILURE)
            return FAILURE;
    }
#else
    if (startrun_parameterfile(cmd, gd) == FAILURE)
        return FAILURE;
#endif

    gd->bytes_tot += sizeof(struct  global_data);
    gd->bytes_tot += sizeof(struct cmdline_data);
    verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                "\n%s: Total allocated %g MByte storage so far.\n",
                        routineName, gd->bytes_tot*INMB);

    class_call_cballs(SetNumberThreads(cmd), errmsg, errmsg);
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
    string routineName = "StartRun (CLASSLIB)";
    struct file_content fc;
    double cpustart = CPUTIME;

    gd->headline0 = head0; gd->headline1 = head1;
    gd->headline2 = head2; gd->headline3 = head3;
    printf("\n%s\n%s: %s\n\t %s\n",
           gd->headline0, gd->headline1, gd->headline2, gd->headline3);
#ifdef GETPARAM
    printf("Version = %s\n", getversion());
#endif

    //B move all these to Startrun_Common... or make an appropriate change
    gd->startrun_cputime = TRUE;
    gd->cmd_allocated = FALSE;
    gd->stopflag = 0;
    gd->cputotalinout = 0.;
    gd->cputotal = 0.;
    gd->bytes_tot = 0;
    gd->sameposcount = 0;
    //E

#ifdef GETPARAM
    cmd->paramfile = GetParam("paramfile");
    if (*(cmd->paramfile) == '-') {
        cBALLS_FAIL(cmd, "bad parameter %s\n", cmd->paramfile);
    }

    if (!strnull(cmd->paramfile)) {
        //B
        class_call_cballs(input_find_file(cmd, gd, cmd->paramfile, &fc, errmsg),
                          errmsg, cmd->error_message);

        if (input_read_from_file(cmd, gd, &fc, errmsg) == FAILURE) {
            parser_free(&fc);
            class_call_cballs(FAILURE, errmsg, cmd->error_message);
        }

        class_call_cballs(parser_free(&fc), errmsg, errmsg);
        
        if (StartRun_Common(cmd, gd) == FAILURE)
            return FAILURE;

        if (PrintParameterFile(cmd, gd, cmd->paramfile) == FAILURE)
            return FAILURE;
        //E
        
    } else {
        if (startrun_cmdline(cmd, gd) == FAILURE)
            return FAILURE;
    }

#else
    //B
    // this segment must be checked!!! (used when GETPARMON = 0)
    class_call_cballs(input_find_file(cmd, gd, cmd->paramfile, &fc, errmsg),
                      errmsg, cmd->error_message);

    if (input_read_from_file(cmd, gd, &fc, errmsg) == FAILURE) {
        parser_free(&fc);
        class_call_cballs(FAILURE, errmsg, cmd->error_message);
    }

    class_call_cballs(parser_free(&fc), errmsg, errmsg);
    
    if (StartRun_Common(cmd, gd) == FAILURE)
        return FAILURE;

    if (PrintParameterFile(cmd, gd, cmd->ParameterFile) == FAILURE)
        return FAILURE;
    //E
#endif


    gd->bytes_tot += sizeof(struct  global_data);
    gd->bytes_tot += sizeof(struct cmdline_data);
    verb_print(cmd->verbose,
               "\n%s: Total allocated %g MByte storage so far.\n",
               routineName, gd->bytes_tot*INMB);

    class_call_cballs(SetNumberThreads(cmd), errmsg, errmsg);

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
    //B
    if (ReadParameterFile(cmd, gd, cmd->paramfile) == FAILURE) {
        if (cmd->error_message[0] == '\0') {
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     "ReadParameterFile: failed parsing parameter file '%s'\n",
                     cmd->paramfile);
        }
        return FAILURE;
    }

    ReadParametersCmdline_short(cmd, gd);
    //E

//B socket:
#ifdef ADDONS
#include "startrun_include_01.h"
#endif
//E

#else // ! GETPARAM
    if (ReadParameterFile(cmd, gd, cmd->ParameterFile) == FAILURE)
        return FAILURE;
#endif // ! GETPARAM

    if (StartRun_Common(cmd, gd) == FAILURE)
        return FAILURE;
    
#ifdef GETPARAM
    if (PrintParameterFile(cmd, gd, cmd->paramfile) == FAILURE)
        return FAILURE;
#else
    if (PrintParameterFile(cmd, gd, cmd->ParameterFile) == FAILURE)
        return FAILURE;
#endif

    return SUCCESS;
}


#ifdef GETPARAM
#define parameter_null	"parameters_null-cballs"

//B Section for reading parameters from the command line

local int startrun_cmdline(struct  cmdline_data* cmd, struct  global_data* gd)
{
    ReadParametersCmdline(cmd, gd);

    if (StartRun_Common(cmd, gd) == FAILURE)
        return FAILURE;

    if (PrintParameterFile(cmd, gd, parameter_null) == FAILURE)
        return FAILURE;

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
    cmd->numthreads = GetiParam("numberThreads");
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
local int ReadParameterFile(struct  cmdline_data* cmd,
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
    FILE *fd = NULL;

  int  i,j,nt;
  int  id[MAXTAGS];
  void *addr[MAXTAGS];
  string *string_slot[MAXTAGS];
  char tag[MAXTAGS][50];

    size_t str_size[MAXTAGS];

#undef SPName
#define SPName(param,paramtext,n)                                       \
  do {                                                                  \
      SET_TAG_NAME(paramtext);                                          \
      string_slot[nt] = &(param);                                       \
      (param) = (string)malloc(n);                                      \
      addr[nt] = (param);                                               \
      str_size[nt] = (n);                                               \
      id[nt] = STRING;                                                   \
      nt++;                                                              \
      if ((param) == NULL) {                                            \
          snprintf(cmd->error_message, _ERRORMSGSIZE_,                  \
                   "%s: not enough memory for parameter '%s'",         \
                   routineName, paramtext);                             \
          goto fail;                                                     \
      }                                                                 \
  } while (0)

    int input_verbose = 2;
    verb_print(input_verbose, "\nparsing input parameters...\n");

    cmd->error_message[0] = '\0';
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
    IPName(cmd->numthreads,"numberThreads");
    SPName(cmd->options,"options",MAXLENGTHOFSTRSCMD);
    //E

//B socket:
#ifdef ADDONS
#include "startrun_include_03.h"
#endif
//E

//B
#ifndef _LINE_LENGTH_MAX_
#define _LINE_LENGTH_MAX_ 1024
#endif
#define _ARGUMENT_LENGTH_MAX_ 1024
        char line[_LINE_LENGTH_MAX_];
        char name[_ARGUMENT_LENGTH_MAX_];
        char value[_ARGUMENT_LENGTH_MAX_];
        int has_line;
        int is_data;
        unsigned long line_number = 0;
//E

    fd = fopen(fname, "r");
    if (fd == NULL) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: parameter file '%s' not found", routineName, fname);
        goto fail;
    }

    while (TRUE) {
        if (read_parameter_line_checked(fd, line, sizeof(line),
                                        &has_line, &line_number, fname,
                                        cmd->error_message,
                                        _ERRORMSGSIZE_) == FAILURE)
            goto fail;
        if (has_line == FALSE)
            break;

        if (parse_parameter_line_checked(line,
                                         name, sizeof(name),
                                         value, sizeof(value),
                                         &is_data, fname, line_number,
                                         cmd->error_message,
                                         _ERRORMSGSIZE_) == FAILURE)
            goto fail;
        if (is_data == FALSE)
            continue;

        for (i = 0, j = -1; i < nt; i++) {
            if (strcmp(name, tag[i]) == 0) {
                j = i;
                tag[i][0] = '\0';
                break;
            }
        }

        if (j < 0) {
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     "%s:%lu: parameter '%s' is unknown or duplicated",
                     fname, line_number, name);
            goto fail;
        }

        switch (id[j]) {
            case DOUBLE:
                if (parse_double_checked(value, (double *)addr[j],
                                         cmd->error_message,
                                         _ERRORMSGSIZE_, name) == FAILURE)
                    goto fail;
                break;
            case STRING:
                if (strcmp(name, "preScript") == 0
                    || strcmp(name, "posScript") == 0) {
                    size_t slen = strlen(value);

                    if (slen < 2 || value[0] != '"'
                        || value[slen - 1] != '"') {
                        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                                 "%s:%lu: %s needs enclosing double quotes",
                                 fname, line_number, name);
                        goto fail;
                    }
                    if (slen - 2 >= str_size[j]) {
                        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                                 "%s:%lu: %s is too long",
                                 fname, line_number, name);
                        goto fail;
                    }

                    memcpy((char *)addr[j], value + 1, slen - 2);
                    ((char *)addr[j])[slen - 2] = '\0';
                } else if (copy_checked((char *)addr[j], str_size[j],
                                        value, name) != 0) {
                    snprintf(cmd->error_message, _ERRORMSGSIZE_,
                             "%s:%lu: string parameter '%s' is too long",
                             fname, line_number, name);
                    goto fail;
                }
                break;
            case INT:
                if (parse_int_checked(value, (int *)addr[j],
                                      cmd->error_message,
                                      _ERRORMSGSIZE_, name) == FAILURE)
                    goto fail;
                break;
            case LONG:
                if (parse_long_checked(value, (long *)addr[j],
                                       cmd->error_message,
                                       _ERRORMSGSIZE_, name) == FAILURE)
                    goto fail;
                break;
            case BOOLEAN:
                if (parse_bool_checked(value, (bool *)addr[j],
                                       cmd->error_message,
                                       _ERRORMSGSIZE_, name) == FAILURE)
                    goto fail;
                break;
        }
    }

    if (fclose(fd) != 0) {
        fd = NULL;
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: could not close parameter file '%s'",
                 routineName, fname);
        goto fail;
    }
    fd = NULL;
    
#ifdef GETPARAM
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
                    if (copy_checked((char *)addr[i], str_size[i],
                                     GetParam(tag[i]), tag[i]) != 0) {
                        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                                 "%s: default string parameter '%s' is too long",
                                 routineName, tag[i]);
                        goto fail;
                    }
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
        }
    }
#endif

    return SUCCESS;

fail:
    if (fd != NULL)
        fclose(fd);
    for (i = 0; i < nt; i++) {
        if (id[i] == STRING && string_slot[i] != NULL
            && *string_slot[i] != NULL) {
            free(*string_slot[i]);
            *string_slot[i] = NULL;
        }
    }
    return FAILURE;

#undef DOUBLE
#undef STRING
#undef INT
#undef LONG
#undef BOOLEAN
#undef MAXTAGS
#undef MAXCHARBUF
#undef SPName
}
//E

int StartRun_Common(struct  cmdline_data* cmd, struct  global_data* gd)
{
    string routineName = "StartRun_Common";
    int ifile;
    double cpustartMiddle;
    const int memory_ninfiles = gd->bodytable_allocated ? gd->ninfiles : 0;
    const bool input_catalogs_in_memory = memory_ninfiles > 0;

    cballs_refresh_option_cache(cmd);

#ifdef THREEPCFCONVERGENCE
    gd->computeTPCF = TRUE;
#else
    gd->computeTPCF = FALSE;
#endif
    gd->inputHeaderFlag = FALSE;
    gd->cmd_allocated = TRUE;
    gd->histograms_allocated = FALSE;
    gd->random_allocated = FALSE;
    gd->gd_allocated = FALSE;
    gd->gd_allocated_2 = FALSE;
    gd->tree_allocated = FALSE;
    gd->bodytable_allocated = input_catalogs_in_memory;

    //B
#if defined(DEBUG) && defined(BODYTABBF_ON)
    bodytabbf = NULL;
#endif
    //E

    gd->outlog = NULL;
    gd->outlogFlagFree = FALSE;
    
    if (strlen(cmd->rootDir)==0 || strnull(cmd->rootDir))
        gd->rootDirFlag = FALSE;
    else
        gd->rootDirFlag = TRUE;

    gd->flagPrint = TRUE;

#ifdef CBALLS_MPI_ENABLED
    if (cballs_mpi_prepare(cmd, gd) == FAILURE)
        return FAILURE;
#endif

    if (scanopt(cmd->options, "make-info")) {
        print_make_info(cmd, gd);
        if (!scanopt(cmd->options, "no-stop")) {        // make a better way out
            gd->stopflag = TRUE;
            return FAILURE;
        }
    }

    if (scanopt(cmd->options, "print-options")) {
        print_options(cmd, gd);
        if (!scanopt(cmd->options, "no-stop")) {        // make a better way out
            gd->stopflag = TRUE;
            return FAILURE;
        }
    }

    if (scanopt(cmd->options, "print-search-methods")) {
        print_search_methods(cmd, gd);
        if (!scanopt(cmd->options, "no-stop")) {
            gd->stopflag = TRUE;
            return FAILURE;
        }
    }

//B socket:
#ifdef ADDONS
#include "startrun_include_04.h"
#endif
//E

#if defined(NMultipoles) && defined(NONORMHIST)
    if (scanopt(cmd->options, "patch-with-all")) {
        gd->pivotCount = 0;
    }
#endif
    gd->pivotNumber = cmd->nbody;

    int output_setup_status = StartOutput(cmd, gd);
#ifdef OCTREE3PCF3DMPI
    output_setup_status = cb3d_mpi_consensus(
        cmd, output_setup_status, "3D output initialization");
#endif
#ifdef LYAFORESTMPI
    output_setup_status = lya_forest_mpi_consensus(
        cmd, output_setup_status, "Ly-alpha output initialization");
#endif
    class_call_cballs(output_setup_status, errmsg, errmsg);
    output_setup_status = SUCCESS;
#ifdef CBALLS_MPI_ENABLED
    if (cballs_mpi_output_enabled(cmd)) {
#endif
        output_setup_status = setFilesDirs(cmd, gd);
        if (output_setup_status == SUCCESS)
            output_setup_status = setFilesDirs_log(cmd, gd);

        if (output_setup_status == SUCCESS) {
            strcpy(gd->mode,"w");
            if (cmd->verbose_log > 0) {
                gd->outlog = fopen(gd->logfilePath, gd->mode);
                if (gd->outlog == NULL) {
                    snprintf(cmd->error_message, _ERRORMSGSIZE_,
                             "\n%s: error opening file '%s' \n",
                             routineName, gd->logfilePath);
                    output_setup_status = FAILURE;
                } else {
                    gd->outlogFlagFree = TRUE;
                }
            }
        }
#ifdef CBALLS_MPI_ENABLED
    }
    output_setup_status = cballs_mpi_consensus(
        cmd, output_setup_status, "MPI output setup");
#endif
    if (output_setup_status == FAILURE) return FAILURE;

     int catalog_setup_status = startrun_getParamsSpecial(cmd, gd);
#ifdef OCTREE3PCF3DMPI
     catalog_setup_status = cb3d_mpi_consensus(
         cmd, catalog_setup_status, "3D startup options");
#endif
#ifdef LYAFORESTMPI
     catalog_setup_status = lya_forest_mpi_consensus(
         cmd, catalog_setup_status, "Ly-alpha startup options");
#endif
     class_call_cballs(catalog_setup_status, errmsg, errmsg);

     search_method_string_to_int(cmd->searchMethod, &gd->searchMethod_int);
     catalog_setup_status = CheckParameters(cmd, gd);
#ifdef OCTREE3PCF3DMPI
     catalog_setup_status = cb3d_mpi_consensus(
         cmd, catalog_setup_status, "3D parameter validation");
#endif
#ifdef LYAFORESTMPI
     catalog_setup_status = lya_forest_mpi_consensus(
         cmd, catalog_setup_status, "Ly-alpha parameter validation");
#endif
     class_call_cballs(catalog_setup_status, errmsg, errmsg);

     catalog_setup_status = random_init(cmd, gd, cmd->seed);
#ifdef OCTREE3PCF3DMPI
     catalog_setup_status = cb3d_mpi_consensus(
         cmd, catalog_setup_status, "3D random initialization");
#endif
#ifdef LYAFORESTMPI
     catalog_setup_status = lya_forest_mpi_consensus(
         cmd, catalog_setup_status, "Ly-alpha random initialization");
#endif
     class_call_cballs(catalog_setup_status, errmsg, errmsg);

     catalog_setup_status = startrun_memoryAllocation(cmd, gd);
#ifdef OCTREE3PCF3DMPI
     catalog_setup_status = cb3d_mpi_consensus(
         cmd, catalog_setup_status, "3D startup allocation");
#endif
#ifdef LYAFORESTMPI
     catalog_setup_status = lya_forest_mpi_consensus(
         cmd, catalog_setup_status, "Ly-alpha startup allocation");
#endif
     class_call_cballs(catalog_setup_status, errmsg, errmsg);

    coordinate_string_to_int(cmd, gd);              // set coordTag
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
            "\n%s: coordTag: %d\n", routineName, gd->coordTag);

//B Pre-processing necessary for reading data files:
    char buf[BUFFERSIZE];
    int preprocessing_owner = TRUE;
    int preprocessing_status = SUCCESS;
#ifdef CBALLS_MPI_ENABLED
    preprocessing_owner = cballs_mpi_output_enabled(cmd);
#endif
    if (scanopt(cmd->options, "pre-processing")) {
        if (preprocessing_owner) {
            cpustartMiddle = CPUTIME;
            if (copy_checked(buf, sizeof(buf), cmd->preScript,
                             "preScript") != 0) {
                snprintf(cmd->error_message, _ERRORMSGSIZE_,
                         "%s: preScript command too long", routineName);
                preprocessing_status = FAILURE;
            } else {
                verb_print_normal_info(cmd->verbose, cmd->verbose_log,
                                       gd->outlog,
                                       "\n%s: pre-processing: executing %s...\n",
                                       routineName, cmd->preScript);
                preprocessing_status =
                    cballs_system_checked(cmd, routineName, buf);
                if (preprocessing_status == SUCCESS) {
                    verb_print_normal_info(cmd->verbose, cmd->verbose_log,
                                           gd->outlog, " done.\n");
                    gd->cputotalinout += CPUTIME - cpustartMiddle;
                    verb_print_normal_info(cmd->verbose, cmd->verbose_log,
                                           gd->outlog,
                                           "%s: cpu time expended in this script %g\n\n",
                                           routineName,
                                           CPUTIME - cpustartMiddle,
                                           PRNUNITOFTIMEUSED);
                }
            }
        }
#ifdef CBALLS_MPI_ENABLED
        preprocessing_status = cballs_mpi_consensus(
            cmd, preprocessing_status, "MPI pre-processing");
#endif
        if (preprocessing_status == FAILURE) return FAILURE;
        if (scanopt(cmd->options, "stop")) {
            verb_print_normal_info(cmd->verbose,
                                   cmd->verbose_log, gd->outlog,
                                   "\n\tMainLoop: stopping...\n\n");
            gd->stopflag = TRUE;
            return SUCCESS;
        }
    }
    if (scanopt(cmd->options, "statistics-histograms")) {
        preprocessing_status = preprocessing_owner
            ? statHistogram(cmd, gd) : SUCCESS;
#ifdef CBALLS_MPI_ENABLED
        preprocessing_status = cballs_mpi_consensus(
            cmd, preprocessing_status, "MPI histogram statistics");
#endif
        if (preprocessing_status == FAILURE) return FAILURE;

        verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                               "\n\tpre-processing: stopping...\n\n");

//        gd->stopflag = TRUE;
//        return SUCCESS;
        if (scanopt(cmd->options, "stop")) {        // make a better way out
            gd->stopflag = TRUE;
            return FAILURE;
        } else {
            gd->stopflag = TRUE;
            return FAILURE;
        }

    }

    if (scanopt(cmd->options, "edge-corrections-from-files")) {
        preprocessing_status = preprocessing_owner
            ? computeEdgeCorrections(cmd, gd) : SUCCESS;
#ifdef CBALLS_MPI_ENABLED
        preprocessing_status = cballs_mpi_consensus(
            cmd, preprocessing_status, "MPI edge-correction preprocessing");
#endif
        if (preprocessing_status == FAILURE) return FAILURE;
        verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                               "\n\tpre-processing: stopping...\n\n");
        gd->stopflag = TRUE;
        return SUCCESS;
    }
//E

//B In this section update computation of rSize
//      and center-of-mass if necessary
//      so we have a common root size and c-of-m
    if (input_catalogs_in_memory) {
        if (gd->ninfiles != memory_ninfiles)
            cBALLS_FAIL(cmd,
                        "%s: parsed %d in-memory descriptors for %d loaded catalogs\n",
                        routineName, gd->ninfiles, memory_ninfiles);
        for (ifile = 0; ifile < memory_ninfiles; ifile++) {
            if (bodytable[ifile] == NULL || gd->nbodyTable[ifile] < 3)
                cBALLS_FAIL(cmd,
                            "%s: in-memory catalog %d is missing or invalid\n",
                            routineName, ifile);
        }
        gd->model_comment = "Python in-memory catalog";
        gd->input_comment = "catalog copied from Python memory";
        verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                               "\n%s: using %d in-memory catalog(s).\n",
                               routineName, memory_ninfiles);
    } else if (scanopt(cmd->options, "read-mask")) {
        bool embedded_octree_mask = FALSE;
#if defined(OCTREE2BALLSOMP) || defined(OCTREE2BALLSMPI) || defined(OCTREEBALLS4OMP) || defined(OCTREEBALLS4MPI)
        bool octree_two_balls = FALSE;
#ifdef OCTREE2BALLSOMP
        octree_two_balls |= gd->searchMethod_int == OCTREE2BALLSMETHOD;
#endif
#ifdef OCTREE2BALLSMPI
        octree_two_balls |= gd->searchMethod_int == OCTREE2BALLSMPIMETHOD;
#endif
#ifdef OCTREEBALLS4OMP
        octree_two_balls |= gd->searchMethod_int == OCTREEBALLS4OMPMETHOD;
#endif
#ifdef OCTREEBALLS4MPI
        octree_two_balls |= gd->searchMethod_int == OCTREEBALLS4MPIMETHOD;
#endif
        embedded_octree_mask = octree_two_balls && gd->ninfiles == 1
            && (strcmp(gd->infilefmtname[0], "columns-ascii-all") == 0
                || strcmp(gd->infilefmtname[0], "binary-all") == 0);
#endif
        if (gd->ninfiles < 2 && !embedded_octree_mask)
            cBALLS_FAIL(cmd, "\t%s:: read-mask ninfiles = %d is absurd\n", routineName, gd->ninfiles);
        ifile=0;
        class_call_cballs(infilefmt_string_to_int(gd->infilefmtname[ifile],
                    &gd->infilefmt_int), errmsg, errmsg);
        class_call_cballs(InputData(cmd, gd, gd->infilenames[ifile],ifile),
                    errmsg, errmsg);
        if (gd->stopflag)
            return SUCCESS;
        
        gd->model_comment = "input data file with mask";
        if (!embedded_octree_mask) {
            ifile=1;
            class_call_cballs(infilefmt_string_to_int(gd->infilefmtname[ifile],
                    &gd->infilefmt_int), errmsg, errmsg);
            class_call_cballs(InputData(cmd, gd, gd->infilenames[ifile],ifile),
                    errmsg, errmsg);
            if (gd->stopflag)
                return SUCCESS;
        }
        
    } else {
        for (ifile=0; ifile<gd->ninfiles; ifile++) {
            if (!strnull(cmd->infile)) {
                class_call_cballs(
                            infilefmt_string_to_int(gd->infilefmtname[ifile],
                            &gd->infilefmt_int), errmsg, errmsg);
                class_call_cballs(InputData(cmd, gd,
                                            gd->infilenames[ifile],ifile),
                                            errmsg, errmsg);
                if (gd->stopflag)
                    return SUCCESS;
                
                gd->model_comment = "input data file";
            } else {
                verb_print_normal_info(cmd->verbose,
                                       cmd->verbose_log, gd->outlog,
                                       "\nNo data catalog was given...");
                verb_print_normal_info(cmd->verbose,
                                       cmd->verbose_log, gd->outlog,
                                       "creating a test model...\n");
                class_call_cballs(TestData(cmd, gd), errmsg, errmsg);
                gd->input_comment = "no data file given";
            }
        }
    }

#if defined(OCTREE3PCF3DOMP) || defined(OCTREE3PCF3DMPI)
    /* Readers with a real LOS column publish it through the flag below.  All
     * other inputs receive unique IDs, so LOS exclusion is deterministic and
     * cannot inspect an uninitialized addon field. */
    for (ifile = 0; ifile < gd->ninfiles; ifile++) {
        if (!gd->octree3pcf3d_los_ids[ifile]) {
            bodyptr p_los;
            DO_BODY(p_los, bodytable[ifile],
                    bodytable[ifile] + gd->nbodyTable[ifile])
                Octree3pcf3dLosId(p_los) = Id(p_los);
            gd->octree3pcf3d_los_ids[ifile] = TRUE;
        }
    }
#endif

    if (gd->inputHeaderFlag==TRUE) return SUCCESS;

    //B consider moving below after computing rsize
    if (scanopt(cmd->options, "all-in-one")) {
        class_call_cballs(InputData_all_in_one(cmd, gd),
                          errmsg, errmsg);
        if (gd->stopflag)
            return SUCCESS;
        
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

    bool preserve_common_catalog_frame = cballs_observer_frame(cmd);
#ifdef OCTREESHEAROMP
    preserve_common_catalog_frame |=
        gd->searchMethod_int == OCTREESHEARMETHOD;
    if (gd->searchMethod_int == OCTREESHEARMETHOD
        && prepare_octree_shear_catalogs(cmd, gd, bodytable,
                                         gd->nbodyTable) == FAILURE)
        return FAILURE;
#endif
#if defined(OCTREE3PCF3DOMP) || defined(OCTREE3PCF3DMPI)
    if (cb3d_is_method(cmd->searchMethod)
        && (scanopt(cmd->options, "survey-estimator-3d")
            || scanopt(cmd->options, "encore-survey-estimator")
            || scanopt(cmd->options, "survey-edge-correction"))) {
        preserve_common_catalog_frame = TRUE;
        if (cb3d_parallel_consensus(cmd,
                cb3d_prepare_common_frame(cmd, gd, bodytable, gd->nbodyTable),
                "3D common catalog frame") == FAILURE)
            return FAILURE;
    }
#endif
#ifdef BALLTREE2BALLSOMP
    if (gd->searchMethod_int == BALLTREE2BALLSMETHOD
        && gd->iCatalogs[0] != gd->iCatalogs[1])
        preserve_common_catalog_frame = TRUE;
#endif
#ifdef BALLTREE2BALLSMPI
    if (gd->searchMethod_int == BALLTREE2BALLSMPIMETHOD
        && gd->iCatalogs[0] != gd->iCatalogs[1])
        preserve_common_catalog_frame = TRUE;
#endif
#ifdef OCTREE2BALLSOMP
    if (gd->searchMethod_int == OCTREE2BALLSMETHOD
        && gd->iCatalogs[0] != gd->iCatalogs[1])
        preserve_common_catalog_frame = TRUE;
#endif
#ifdef BALLTREE2BALLSOMP3PCF
    if (gd->searchMethod_int == BALLTREE2BALLSOMP3PCFMETHOD
        && gd->iCatalogs[0] != gd->iCatalogs[1])
        preserve_common_catalog_frame = TRUE;
#endif
#ifdef BALLTREE2BALLSMPI3PCF
    if (gd->searchMethod_int == BALLTREE2BALLSMPI3PCFMETHOD
        && gd->iCatalogs[0] != gd->iCatalogs[1])
        preserve_common_catalog_frame = TRUE;
#endif
    if (scanopt(cmd->options, "read-mask")) {
        ifile=0;
        //B
        cellptr root = NULL;                // Set it up a temporal root
        int rc = FAILURE;

        if (cballs_calloc_checked((void **)&root, 1, sizeof(body),
                                  "startup root", cmd->error_message,
                                  _ERRORMSGSIZE_) == FAILURE)
            goto cleanup_root_read_mask;
        if (FindRootCenter(cmd, gd, bodytable[ifile], gd->nbodyTable[ifile], ifile, root) == FAILURE)
            goto cleanup_root_read_mask;
        if (!preserve_common_catalog_frame
            && centerBodies(bodytable[ifile], gd->nbodyTable[ifile],
                            ifile, root) == FAILURE)
            goto cleanup_root_read_mask;
        if (FindRootCenter(cmd, gd, bodytable[ifile], gd->nbodyTable[ifile], ifile, root) == FAILURE)
            goto cleanup_root_read_mask;
        gd->rSizeTable[ifile] = 1.0;
        if (expandbox(cmd, gd, bodytable[ifile], gd->nbodyTable[ifile], ifile, root) == FAILURE)
            goto cleanup_root_read_mask;

        rc = SUCCESS;

        cleanup_root_read_mask:
        free(root);
#ifdef OCTREE3PCF3DMPI
        rc = cb3d_mpi_consensus(cmd, rc, "3D masked root sizing");
#endif
        if (rc == FAILURE)
            return FAILURE;
        //E

        if (cmd->rangeN > gd->rSizeTable[ifile])
        verb_print_warning(cmd->verbose,
        "\n%s: warning! rangeN (%g) is greather than rSize (%g) of the system...\n",
                           routineName, cmd->rangeN, gd->rSizeTable[ifile]);
    } else {
        for (ifile=0; ifile<gd->ninfiles; ifile++) {
            //B
            cellptr root = NULL;                // Set it up a temporal root
            int rc = FAILURE;

            if (cballs_calloc_checked((void **)&root, 1, sizeof(body),
                                      "startup root", cmd->error_message,
                                      _ERRORMSGSIZE_) == FAILURE)
                goto cleanup_root;
            if (FindRootCenter(cmd, gd, bodytable[ifile], gd->nbodyTable[ifile], ifile, root) == FAILURE)
                goto cleanup_root;
            if (!preserve_common_catalog_frame
                && centerBodies(bodytable[ifile], gd->nbodyTable[ifile],
                                ifile, root) == FAILURE)
                goto cleanup_root;
            if (FindRootCenter(cmd, gd, bodytable[ifile], gd->nbodyTable[ifile], ifile, root) == FAILURE)
                goto cleanup_root;
            gd->rSizeTable[ifile] = 1.0;
            if (expandbox(cmd, gd, bodytable[ifile], gd->nbodyTable[ifile], ifile, root) == FAILURE)
                goto cleanup_root;

            rc = SUCCESS;

            cleanup_root:
            free(root);
#ifdef OCTREE3PCF3DMPI
            rc = cb3d_mpi_consensus(cmd, rc, "3D root sizing");
#endif
            if (rc == FAILURE)
                return FAILURE;
            //E
            
            if (cmd->rangeN > gd->rSizeTable[ifile])
                verb_print_warning(cmd->verbose,
    "\n%s: warning! rangeN (%g) is greather than rSize (%g) of the system...\n",
                           routineName, cmd->rangeN, gd->rSizeTable[ifile]);
        }
    }

    //B set output comment
    if (! strnull(cmd->outfile)) {
        gd->output_comment = "output data file was given";
    } else {
        gd->output_comment = "no output data file was given";
    }
    //E

//B Tree search:
    gd->Rcut = cmd->rangeN;                         // Maximum search radius
    if (scanopt(cmd->options, "Rcut/theta")) {
        if (cmd->theta>0)
            gd->Rcut /= cmd->theta;                 // Maximum search radius
    }
    gd->RcutSq = gd->Rcut*gd->Rcut;
//E
    if (cmd->useLogHist) {
        real rBin, rbinlog;
        if (cmd->rminHist==0) {
            //B rminHist = 0 not allowed when useLogHist is true
            //  therefore this segment does not ocurr...
            gd->deltaRV = dvector(1,cmd->sizeHistN);
            gd->gd_allocated_2 = TRUE;
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
            gd->deltaR = (rlog10(cmd->rangeN) - rlog10(cmd->rminHist))
                       / cmd->sizeHistN;
            if (!isfinite(gd->deltaR) || gd->deltaR <= 0.0)
                cBALLS_FAIL(cmd,
                            "%s: logarithmic histogram bin width must be finite and positive (deltaR=%g)\n",
                            routineName, gd->deltaR);
            //B allocated after startrun_memoryAllocation
            //  deallocate before deallocate arrays in startrun_memoryAllocation
            gd->deltaRV = dvector(1,cmd->sizeHistN);
            gd->gd_allocated_2 = TRUE;
            gd->ddeltaRV = dvector(1,cmd->sizeHistN-1);
            gd->bytes_tot += (cmd->sizeHistN)*sizeof(real);
            gd->bytes_tot += (cmd->sizeHistN-1)*sizeof(real);
            //E
            verb_log_print(cmd->verbose_log, gd->outlog,
                           "deltaRV (deltaR=%lf logscale):\n", gd->deltaR);
            int n=1;                                // check this value of 1
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

    if (!isfinite(gd->deltaR) || gd->deltaR <= 0.0)
        cBALLS_FAIL(cmd,
                    "%s: histogram bin width must be finite and positive (deltaR=%g)\n",
                    routineName, gd->deltaR);

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
// in common_defs.h
// 180*60/Pi
//#define RADTOARCMIN   3437.74677
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
//#undef RADTOARCMIN
#endif
//E

    gd->stepState = cmd->stepState/10;
    if (gd->stepState < 1)
        gd->stepState = 1;

    gd->deltaPhi = TWOPI/cmd->sizeHistPhi;
    if (!isfinite(gd->deltaPhi) || gd->deltaPhi <= 0.0)
        cBALLS_FAIL(cmd,
                    "%s: angular bin width must be finite and positive (deltaPhi=%g)\n",
                    routineName, gd->deltaPhi);
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

    return SUCCESS;
}


//B Section of parameter check
local int CheckParameters(struct  cmdline_data* cmd, struct  global_data* gd)
{
// If it is necessary: an item in cmdline_defs.h must have an item here::
    string routineName = "CheckParameters";

    cballs_refresh_option_cache(cmd);

#if defined(OCTREEBALLS4OMP) || defined(OCTREEBALLS4MPI)
    if (cmd->searchMethod && (strcmp(cmd->searchMethod, "octree-balls4-omp") == 0
        || strcmp(cmd->searchMethod, "octree-balls4-mpi") == 0)
        && cballs_opt_smooth_pivot(cmd))
        cBALLS_FAIL(cmd, "%s does not support options=smooth-pivot\n", cmd->searchMethod);
#endif

    //B Parameters related to the searching method
    if (cmd->useLogHist==FALSE &&
        (strcmp(cmd->searchMethod,"balls-omp") == 0))
        cBALLS_FAIL(cmd, "%s: can´t have loghist false and balls-omp (%d %s)\n",
                    routineName, cmd->useLogHist, cmd->searchMethod);
#ifdef TPCF
#ifndef USEGSL
        //  for recursivity needs that at least 3 multipoles be evaluated
        if (scanopt(cmd->options, "out-HistZetaG")) {
            if (cmd->mChebyshev + 1 < 3
                || (cmd->mChebyshev + 1)&(cmd->mChebyshev)) {
                verb_print(cmd->verbose,
                    "\n%s: using option out-HistZetaG...\n", routineName);
                cBALLS_FAIL(cmd, "%s: %s (=%d)\n\t\t\t%s\n",
                            routineName, "absurd value for mChebyshev + 1",
                            cmd->mChebyshev+1, "must be positive and a power of 2");
            }
        } else {
            if (cmd->mChebyshev + 1 < 3)
                cBALLS_FAIL(cmd, "%s: %s (=%d)\n\t\t\tmust be positive\n",
                            "absurd value for mChebyshev + 1",
                            routineName, cmd->mChebyshev+1);
        }
#else
        //  for recursivity needs that at least 3 multipoles be evaluated
        if (cmd->mChebyshev + 1 < 3)
            cBALLS_FAIL(cmd, "CheckParameters: absurd value for mChebyshev + 1 (=%d)\n\t\t\tmust be positive\n", cmd->mChebyshev+1);
#endif
#endif
    if (cmd->mChebyshev + 1 > cmd->sizeHistPhi)
        cBALLS_FAIL(cmd, "CheckParameters: mChebyshev + 1 must be less than sizeHistPhi\n");
    if (cmd->nsmooth < 1)
        cBALLS_FAIL(cmd, "%s: absurd value for nsmooth: %d\n",
                    routineName, cmd->nsmooth);

    if (!strnull(cmd->rsmooth)) {
        if (gd->rsmooth[0] < 0.0 || gd->rsmoothFlag==FALSE)
            cBALLS_FAIL(cmd, "CheckParameters: absurd value for rsmooth (%s, %g, %d)\n",
                        cmd->rsmooth, gd->rsmooth[0], gd->rsmoothFlag);
    } else {
        gd->rsmooth[0] = 0.0;
    }
    if (!isfinite(cmd->theta) || cmd->theta < 0.0)
        cBALLS_FAIL(cmd,
                    "%s: theta must be finite and non-negative (got %g)\n",
                    routineName, cmd->theta);
    //E

    //B Parameters to control the I/O file(s)
    // Input catalog parameters
    // Output parameters
    // Parameters to set a region in the sky, for example for Takahashi data
    if (cmd->thetaL < 0 || cmd->thetaL > PI)
        cBALLS_FAIL(cmd, "CheckParameters: absurd value for thetaL (must be in the range 0--pi)\n");
    if (cmd->thetaR < 0 || cmd->thetaR > PI)
        cBALLS_FAIL(cmd, "CheckParameters: absurd value for thetaR (must be in the range 0--pi)\n");
    if (cmd->phiL < 0 || cmd->phiL > TWOPI)
        cBALLS_FAIL(cmd, "CheckParameters: absurd value for phiL (must be in the range 0--2pi)\n");
    if (cmd->phiR < 0 || cmd->phiR > TWOPI)
        cBALLS_FAIL(cmd, "CheckParameters: absurd value for phiR (must be in the range 0--2pi)\n");
    if (cmd->thetaL > cmd->thetaR)
        cBALLS_FAIL(cmd, "CheckParameters: absurd value for thetaL (must not be greater than thetaR)\n");
    if (cmd->phiL > cmd->phiR)
        cBALLS_FAIL(cmd, "CheckParameters: absurd value for phiL (must not be greater than phiR)\n");
    //E

    //B Parameters to control histograms and their output files
    if (cmd->sizeHistN < 2)
        cBALLS_FAIL(cmd, "%s: sizeHistN must be at least 2 (got %d)\n",
                    routineName, cmd->sizeHistN);
    if (scanopt(cmd->options, "out-m-HistZeta") && cmd->sizeHistN < 4)
        cBALLS_FAIL(cmd,
                    "%s: out-m-HistZeta requires sizeHistN >= 4 (got %d)\n",
                    routineName, cmd->sizeHistN);

    if (!isfinite(cmd->rangeN) || cmd->rangeN <= 0.0)
        cBALLS_FAIL(cmd, "%s: rangeN must be finite and positive (got %g)\n",
                    routineName, cmd->rangeN);
    if (!isfinite(cmd->rminHist) || cmd->rminHist < 0.0
        || cmd->rminHist >= cmd->rangeN)
        cBALLS_FAIL(cmd,
                    "%s: rminHist must be finite and satisfy 0 <= rminHist < rangeN (got %g, rangeN=%g)\n",
                    routineName, cmd->rminHist, cmd->rangeN);
    if (cmd->useLogHist == TRUE && cmd->rminHist <= 0.0)
        cBALLS_FAIL(cmd,
                    "%s: logarithmic histograms require rminHist > 0 (got %g)\n",
                    routineName, cmd->rminHist);
    if (cmd->useLogHist == TRUE && cmd->logHistBinsPD <= 0)
        cBALLS_FAIL(cmd,
                    "%s: logarithmic histograms require logHistBinsPD > 0 (got %d)\n",
                    routineName, cmd->logHistBinsPD);

    if (cmd->sizeHistPhi < 2 || cmd->sizeHistPhi&(cmd->sizeHistPhi-1))
        cBALLS_FAIL(cmd, "CheckParameters: absurd value for sizeHistPhi\n\tmust be a power of 2\n");
    //
    //E
    
    //B Set of parameters needed to construct a test model
    if (cmd->nbody < 3) {
        cBALLS_FAIL(cmd, "%s: absurd value for nbody: %" INTEGER_FMT "\n",
                    routineName, cmd->nbody);
    }
    if (!isfinite(cmd->lengthBox) || cmd->lengthBox <= 0.0)
        cBALLS_FAIL(cmd,
                    "%s: lengthBox must be finite and positive (got %g)\n",
                    routineName, cmd->lengthBox);
    //E

    //B Miscellaneous parameters
    if (cmd->stepState <= 0)
        cBALLS_FAIL(cmd, "%s: absurd value for stepState must be an integer > 0\n",
                    routineName);
    if (cmd->verbose < 0)
        cBALLS_FAIL(cmd, "%s: absurd value for verbose must be an integer >= 0\n",
                    routineName);
    if (cmd->verbose_log < 0)
        cBALLS_FAIL(cmd, "%s: absurd value for verbose_log must be an integer >= 0\n",
                    routineName);
    if (cmd->numthreads <= 0)
        cBALLS_FAIL(cmd, "%s: absurd value for numberThreads must be an integer >= 0\n",
                    routineName);
    //
    gd->bh86 = scanopt(cmd->options, "bh86");       // Barnes, J. & Hut, P.
                                                    //  1986. Nature 324,
                                                    //  446.
    gd->sw94 = scanopt(cmd->options, "sw94");       // Salmon, J.K. &
                                                    //  Warren, M.S. 1994.
                                                    //  J. Comp. Phys. 111,
                                                    //  136

    if (gd->bh86 && gd->sw94)
        cBALLS_FAIL(cmd, "%s: incompatible options bh86 and sw94\n", routineName);
    //

    //B options section
    if (scanopt(cmd->options, "and-CF")) {
        if (!scanopt(cmd->options, "compute-HistN")) {
            cBALLS_FAIL(cmd, "\n%s: 'and-CF' in options must be %s\n %s\n",
            routineName, "accompanied by 'compute-HistN'.",
                        "Useful only for periodic box catalogs.");
        }
    }

    if (scanopt(cmd->options, "edge-corrections")) {
        if (!gd->computeTPCF)
            cBALLS_FAIL(cmd, "%s: edge-corrections requires TPCFON=1\n", routineName);
        if (scanopt(cmd->options, "only-2pcf"))
            cBALLS_FAIL(cmd, "%s: edge-corrections requires 3PCF; remove only-2pcf\n",
                        routineName);
        if (cmd->mChebyshev < 0 || cmd->mChebyshev > (INT_MAX - 1) / 2)
            cBALLS_FAIL(cmd, "%s: edge-correction window order is out of range\n",
                        routineName);
        if (!scanopt(cmd->options, "no-normalize-HistZeta")) {
            cBALLS_FAIL(cmd, "\n%s: 'edge-corrections' in options must be %s\n",
                        routineName, "accompanied by 'no-normalize-HistZeta'.\n");
        }
    }

    if (scanopt(cmd->options, "read-mask")) {
        verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
            "\n\t%s: warning! 'read-mask' in options must be %s\n\t\t %s\n\n",
            routineName,
            "accompanied by properly setting infile, infileformat and iCatalogs parameters'.",
                            "For a FITS companion mask, use the same Nside as the data catalog.");
    }

    if (scanopt(cmd->options, "no-normalize-HistZeta")) {
        if (cballs_raw_legacy_multipoles(cmd) && cballs_opt_smooth_pivot(cmd))
            cBALLS_FAIL(cmd, "%s: raw distinct legacy multipoles require smooth-pivot disabled\n",
                        routineName);
        if (!scanopt(cmd->options, "edge-corrections")) {
            verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                    "\n\t%s: warning! 'no-normalize-HistZeta' in options %s\n", routineName,
                                "will give no-normalized 3pcf.");

        }
    }

    if (scanopt(cmd->options, "patch")) {
        verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                "\n\t%s: warning! 'patch' in options %s\n\t\t %s\n\n", routineName,
                               "is active only for takahashi and fits-healpix full-sky catalogs.",
                               "Use in combination with parameters thetaL, thetaR, phiL and phiR");
    }
    //E

//B socket:
#ifdef ADDONS
#include "startrun_include_07.h"
#endif
//E

    if (gd->searchMethod_int < 0) {
        cBALLS_FAIL(cmd,
            "%s: searchMethod '%s' is not available in this executable; "
            "enable the corresponding addon in addons/Makefile_addons_settings "
            "or choose a compiled search method\n",
            routineName, cmd->searchMethod);
    }

#if defined(BALLS) && defined(OCTREESMOOTHING) && defined(TREENODEALLBODIES) \
    && !defined(ALLOW_BALLS_WITH_SMOOTHING)
    if (strcmp(cmd->searchMethod, "balls-omp") == 0) {
        cBALLS_FAIL(cmd,
            "%s: balls-omp is not validated with OCTREESMOOTHING and "
            "TREENODEALLBODIES enabled; rebuild without smoothing or set "
            "ALLOW_BALLS_WITH_SMOOTHING=1 for experimental use\n",
            routineName);
    }
#endif

    return SUCCESS;
}
//E


//B Section of parameter writing to a file
#define FMTT    "%-35s = %s\n"
#define FMTTS    "%-35s = \"%s\"\n"
#define FMTI    "%-35s = %d\n"
#define FMTIL    "%-35s = %ld\n"
#define FMTR	"%-35s = %g\n"

int PrintParameterFile(struct  cmdline_data *cmd,
                       struct  global_data* gd, char *fname)
{
// Every item in cmdline_defs.h must have an item here::
    string routineName = "PrintParameterFile";

    FILE *fdout = NULL;
    char buf[MAXLENGTHOFSTRSCMD + MAXLENGTHOFFILES + 32];

    if (gd->flagPrint==TRUE && gd->rootDirFlag==TRUE) {
        //B Look for "/" if fname is composed: path and filename
        int ndefault = 0;
        int ipos = 0;
        const char *dp = fname;
        int nwritten;
        for (size_t i = 0; fname[i] != '\0'; i++) {
            if(fname[i] == '/') {
                ipos = (int)i + 1;
                ndefault++;
            }
        }

        if (ndefault != 0) {
            dp = fname + ipos;
            verb_print_q(3,cmd->verbose,
                         "PrintParameterFile: '/' counts %d pos %d and %s\n",
                         ndefault, ipos, dp);
        }

        nwritten = snprintf(buf, sizeof(buf), "%s/%s%s",
                            cmd->rootDir, dp, "-usedvalues");
        if (nwritten < 0 || (size_t)nwritten >= sizeof(buf))
            cBALLS_FAIL(cmd, "%s: used-values path too long\n", routineName);
        //E
        
        OPEN_OUTPUT_OR_FAIL(fdout, buf, "w!");
        
            //B Parameters related to the searching method
        WRITE_OUTPUT_OR_FAIL(fdout, buf, FMTT, "searchMethod",
                             cmd->searchMethod);
        WRITE_OUTPUT_OR_FAIL(fdout, buf, FMTI,"mChebyshev",cmd->mChebyshev);
        WRITE_OUTPUT_OR_FAIL(fdout, buf, FMTI,"nsmooth",cmd->nsmooth);
        WRITE_OUTPUT_OR_FAIL(fdout, buf, FMTT,"rsmooth",cmd->rsmooth);
        WRITE_OUTPUT_OR_FAIL(fdout, buf, FMTR,"theta",cmd->theta);
        WRITE_OUTPUT_OR_FAIL(fdout, buf, FMTT,"usePeriodic",
                             cmd->usePeriodic ? "true" : "false");
            //E
            
            //B Parameters to control the I/O file(s)
            // Input catalog parameters
        WRITE_OUTPUT_OR_FAIL(fdout, buf, FMTT,"infile",cmd->infile);
        WRITE_OUTPUT_OR_FAIL(fdout, buf, FMTT,"infileformat",cmd->infilefmt);
        WRITE_OUTPUT_OR_FAIL(fdout, buf, FMTT,"iCatalogs",cmd->iCatalogs);
            // Output parameters
        WRITE_OUTPUT_OR_FAIL(fdout, buf, FMTT,"rootDir",cmd->rootDir);
        WRITE_OUTPUT_OR_FAIL(fdout, buf, FMTT,"outfile",cmd->outfile);
        WRITE_OUTPUT_OR_FAIL(fdout, buf, FMTT,"outfileformat",cmd->outfilefmt);
            // Parameters to set a region in the sky,
            //  for example for Takahashi data set
        WRITE_OUTPUT_OR_FAIL(fdout, buf, FMTR,"thetaL",cmd->thetaL);
        WRITE_OUTPUT_OR_FAIL(fdout, buf, FMTR,"thetaR",cmd->thetaR);
        WRITE_OUTPUT_OR_FAIL(fdout, buf, FMTR,"phiL",cmd->phiL);
        WRITE_OUTPUT_OR_FAIL(fdout, buf, FMTR,"phiR",cmd->phiR);
            //E
            
            //B Parameters to control histograms and their output files
        WRITE_OUTPUT_OR_FAIL(fdout, buf, FMTT,"useLogHist",
                             cmd->useLogHist ? "true" : "false");
        WRITE_OUTPUT_OR_FAIL(fdout, buf, FMTI,
                             "logHistBinsPD",cmd->logHistBinsPD);
            //
        WRITE_OUTPUT_OR_FAIL(fdout, buf, FMTI,"sizeHistN",cmd->sizeHistN);
        WRITE_OUTPUT_OR_FAIL(fdout, buf, FMTR,"rangeN",cmd->rangeN);
        WRITE_OUTPUT_OR_FAIL(fdout, buf, FMTR,"rminHist",cmd->rminHist);
        WRITE_OUTPUT_OR_FAIL(fdout, buf, FMTI,"sizeHistPhi",cmd->sizeHistPhi);
            //
        WRITE_OUTPUT_OR_FAIL(fdout, buf, FMTT,
                             "histNNFileName",cmd->histNNFileName);
        WRITE_OUTPUT_OR_FAIL(fdout, buf, FMTT,
                             "histXi2pcfFileName",cmd->histXi2pcfFileName);
        WRITE_OUTPUT_OR_FAIL(fdout, buf, FMTT,
                             "histZetaFileName",cmd->histZetaFileName);
        WRITE_OUTPUT_OR_FAIL(fdout, buf, FMTT,
                             "suffixOutFiles",cmd->suffixOutFiles);
            //E
            
            //B Set of parameters needed to construct a test model
        WRITE_OUTPUT_OR_FAIL(fdout, buf, FMTI,"seed",cmd->seed);
        WRITE_OUTPUT_OR_FAIL(fdout, buf, FMTT,"testmodel",cmd->testmodel);
#ifdef LONGINT
        WRITE_OUTPUT_OR_FAIL(fdout, buf, FMTIL,"nbody",cmd->nbody);
#else
        WRITE_OUTPUT_OR_FAIL(fdout, buf, FMTI,"nbody",cmd->nbody);
#endif
        WRITE_OUTPUT_OR_FAIL(fdout, buf, FMTR,"lengthBox",cmd->lengthBox);
            //E

            //B Miscellaneous parameters
        WRITE_OUTPUT_OR_FAIL(fdout, buf, FMTTS,"preScript",cmd->preScript);
        WRITE_OUTPUT_OR_FAIL(fdout, buf, FMTTS,"posScript",cmd->posScript);
#ifdef LONGINT
        WRITE_OUTPUT_OR_FAIL(fdout, buf, FMTIL,"stepState",cmd->stepState);
#else
        WRITE_OUTPUT_OR_FAIL(fdout, buf, FMTI,"stepState",cmd->stepState);
#endif
        WRITE_OUTPUT_OR_FAIL(fdout, buf, FMTI,"verbose",cmd->verbose);
        WRITE_OUTPUT_OR_FAIL(fdout, buf, FMTI,"verbose_log",cmd->verbose_log);
        WRITE_OUTPUT_OR_FAIL(fdout, buf, FMTI,"numberThreads",cmd->numthreads);
        if (cmd->verbose>=VERBOSEDEBUGINFO) {
            verb_print(cmd->verbose, "\n%s: PrintParamterFile: preScript: %s\n",
                       routineName, cmd->preScript);
            verb_print(cmd->verbose, "\n%s: PrintParamterFile: posScript: %s\n",
                           routineName, cmd->posScript);
        }
        WRITE_OUTPUT_OR_FAIL(fdout, buf, FMTT,"options",cmd->options);
            //E
            
//B socket:
#ifdef ADDONS
#include "startrun_include_08.h"
#endif
//E
            
        WRITE_OUTPUT_OR_FAIL(fdout, buf, "\n\n");
        CLOSE_OUTPUT_OR_FAIL(fdout, buf);
    } // ! gd->flagPrint==TRUE && gd->rootDirFlag==TRUE

    return SUCCESS;
}

#undef FMTT
#undef FMTTS
#undef FMTI
#undef FMTR
//E


local int random_init(struct  cmdline_data* cmd,
                      struct  global_data* gd, int seed)
{
    string routineName = "random_init";

    gd->random_allocated = FALSE;

#ifdef USEGSL
    const gsl_rng_type * T;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    gd->r = gsl_rng_alloc (T);                      // Deallocate at EndRun...
    if (gd->r == NULL) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: memory allocation failed for GSL random generator",
                 routineName);
        return FAILURE;
    }
    gd->random_allocated = TRUE;
    gd->gd_allocated = TRUE;
    verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                          "\n%s: gd->r and seed = %d %d\n",
                          routineName, *(gd->r), seed);
    gsl_rng_set(gd->r, seed);
    verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                          "\n%s: gd->r and seed = %d %d\n",
                          routineName, *(gd->r), seed);
    r_gsl = gd->r;
    gd->bytes_tot += (1)*sizeof(gsl_rng_type);

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
    gd->histograms_allocated = TRUE;
    //B PXD functions
#ifdef PXD
    gd->vecPXD = dvector(1,cmd->sizeHistN);
    //B offset at 0 in order to work with Cython
    gd->matPXD = dmatrix(0,cmd->sizeHistN-1,0,cmd->sizeHistN-1);
    //E
    bytes_tot_local += cmd->sizeHistN*sizeof(real);
    bytes_tot_local += cmd->sizeHistN*cmd->sizeHistN*sizeof(real);
    gd->rBins = dvector(1,cmd->sizeHistN);
    gd->histZetaMFlatten = dvector(1,cmd->sizeHistN*cmd->sizeHistN);
#endif
    //E PXD functions
    gd->histNN = dvector(1,cmd->sizeHistN);
    gd->histCF = dvector(1,cmd->sizeHistN);
    gd->histNNSub = dvector(1,cmd->sizeHistN);
    //B only in search_direct_method_simple
//    gd->histW = dvector(1,cmd->sizeHistN);
//    gd->histWW = dvector(1,cmd->sizeHistN);
    //E
    // 2pcf
    gd->histNNSubXi2pcf = dvector(1,cmd->sizeHistN);
#ifdef SMOOTHPIVOT
    gd->histNNSubXi2pcftotal = dvector(1,cmd->sizeHistN);
#endif
    //
    gd->histNNN = dvector(1,cmd->sizeHistN);
    gd->histXi2pcf = dvector(1,cmd->sizeHistN);
    //B cross
    gd->histXi2pcf12 = dvector(1,cmd->sizeHistN);
    gd->histXi2pcf13 = dvector(1,cmd->sizeHistN);
    //E

    bytes_tot_local += 8*cmd->sizeHistN*sizeof(real);
    bytes_tot_local += cmd->sizeHistN*cmd->sizeHistN*sizeof(real);

#ifdef TPCF
        gd->histXicos = dmatrix(1,cmd->mChebyshev+1,1,cmd->sizeHistN);
        gd->histXisin = dmatrix(1,cmd->mChebyshev+1,1,cmd->sizeHistN);
        bytes_tot_local += 2*(cmd->mChebyshev+1)*cmd->sizeHistN*sizeof(real);
        gd->histZetaM = dmatrix3D(1,cmd->mChebyshev+1,1,cmd->sizeHistN,
                                  1,cmd->sizeHistN);
        gd->histZetaM_EE = dmatrix3D(1,cmd->mChebyshev+1,1,cmd->sizeHistN,
                                  1,cmd->sizeHistN);
        gd->histZetaM_EE_Im = dmatrix3D(1,cmd->mChebyshev+1,
                                       1,cmd->sizeHistN,1,cmd->sizeHistN);
        bytes_tot_local +=
            3.0*(cmd->mChebyshev+1)*cmd->sizeHistN*cmd->sizeHistN*sizeof(real);

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

    //B this was in startrun_include_10.h -> startrun_octree_kkk_omp_10.h above...
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

global int cballs_search_method_id(const char *method)
{
    int method_id = -1;

    if (method == NULL) return -1;
    search_method_string_to_int((string)method, &method_id);
    return method_id;
}


//B
#ifdef OPENMPCODE
int SetNumberThreads(struct  cmdline_data *cmd)
{
    if (cmd == NULL || cmd->numthreads < 1) return FAILURE;
    
    omp_set_num_threads(cmd->numthreads);

    return SUCCESS;
}
#else // dummy for no OPENMCODE
int SetNumberThreads(struct  cmdline_data *cmd)
{
    return SUCCESS;
}
#endif
//E

//B Special routines to scan command line

local int startrun_getParamsSpecial(struct  cmdline_data* cmd,
                                    struct  global_data* gd)
{
    string routineName = "startrun_getParamsSpecial";
    char *pch;
    int nitems, ndummy=1;
    char inputnametmp[MAXLENGTHOFSTRSCMD];
    int i;

    if (strnull(cmd->infile)) {
        verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                        "\n%s: no inputfile was given, making data ...\n",
                        routineName);
        gd->ninfiles=1;                              // To test data...
    } else {
        if (copy_checked(inputnametmp, sizeof(inputnametmp),
                         cmd->infile, "infile") != 0)
            cBALLS_FAIL(cmd, "%s: infile option too long\n", routineName);
        verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                              "\n%s: Splitting string \"%s\" in tokens:\n",
                              routineName, inputnametmp);
        gd->ninfiles=0;
        pch = strtok(inputnametmp," ,");
        while (pch != NULL) {
            if (gd->ninfiles >= MAXITEMS)
                cBALLS_FAIL(cmd, "%s: too many input files\n", routineName);
            if (copy_checked(gd->infilenames[gd->ninfiles],
                             sizeof(gd->infilenames[gd->ninfiles]),
                             pch, "infile token") != 0)
                cBALLS_FAIL(cmd, "%s: infile token too long\n", routineName);
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

    if (!strnull(cmd->infilefmt)) {
        if (copy_checked(inputnametmp, sizeof(inputnametmp),
                         cmd->infilefmt, "infileformat") != 0)
            cBALLS_FAIL(cmd, "%s: infileformat option too long\n", routineName);
        verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                              "\n%s: Splitting string \"%s\" in tokens:\n",
                            routineName, inputnametmp);
        nitems=0;
        pch = strtok(inputnametmp," ,");
        while (pch != NULL) {
            if (nitems >= MAXITEMS)
                cBALLS_FAIL(cmd, "%s: too many infileformat items\n", routineName);
            if (copy_checked(gd->infilefmtname[nitems],
                             sizeof(gd->infilefmtname[nitems]),
                             pch, "infileformat token") != 0)
                cBALLS_FAIL(cmd, "%s: infileformat token too long\n",
                            routineName);
            ++nitems;
            verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                                  "%s: %s\n",
                                routineName, gd->infilefmtname[nitems-1]);
            pch = strtok (NULL, " ,");
        }
        verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                              "%s: num. of items in infilefmt %s =%d\n",
                            routineName, cmd->infilefmt, nitems);
        if (nitems != gd->ninfiles)
            cBALLS_FAIL(cmd,
            "\nstartrun_Common: nitems must be equal to number of files\n\n");

    } else {
        nitems=0;
        pch = "columns-ascii";
        if (copy_checked(gd->infilefmtname[nitems],
                         sizeof(gd->infilefmtname[nitems]),
                         pch, "default infileformat") != 0)
            cBALLS_FAIL(cmd, "%s: default infileformat too long\n",
                        routineName);
        ++nitems;
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
            cBALLS_FAIL(cmd,
            "\nstartrun_Common: nitems must be equal to number of files\n\n");
    }

    class_call_cballs(scanrOption(cmd, gd, cmd->rsmooth, gd->rsmooth,
                                  &nitems, ndummy, 2, "rsmooth"),
                      errmsg, errmsg);

    gd->rsmoothFlag = TRUE;
    if (nitems!=1 && !strnull(cmd->rsmooth)) {
        verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                            "%s: option: rsmooth=%s is not valid... going out\n",
                            routineName, cmd->rsmooth);
        gd->rsmoothFlag = FALSE;
    }

    class_call_cballs(scaniOption(cmd, gd, cmd->iCatalogs, gd->iCatalogs,
                                  &nitems, gd->ninfiles, 0, "iCatalogs"),
                      errmsg, errmsg);

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

//B socket:
#ifdef ADDONS
#include "startrun_include_12.h"
#endif
//E

    return SUCCESS;
}

local int scaniOption(struct  cmdline_data* cmd, struct  global_data* gd,
                      string optionstr, int *option, int *noption,
                      int nfiles, int flag, string message)
{
    string routineName = "scaniOption";
    char *pch;
    char *poptionstr[MAXITEMS], optiontmp[MAXLENGTHOFSTRSCMD];
    int i;

    verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                           "\n%s: Processing '%s' option:\n",
                        routineName, message);

    verb_log_print(cmd->verbose_log, gd->outlog,
                           "%s: Splitting string \"%s=%s\" in tokens:\n",
                        routineName, message, optionstr);

    if (!strnull(optionstr)) {
        if (copy_checked(optiontmp, sizeof(optiontmp),
                         optionstr, message) != 0)
            cBALLS_FAIL(cmd, "%s: option '%s' too long\n",
                        routineName, message);
        verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                               "%s: Splitting string \"%s\" in tokens:\n",
                            routineName, optiontmp);
        *noption=0;
        pch = strtok(optiontmp," ,");
        while (pch != NULL) {
            if (*noption >= MAXITEMS)
                cBALLS_FAIL(cmd, "\n%s: noption = %d must be less than "
                            "the maximum num. of lines\n\n",
                            routineName, *noption + 1);
            poptionstr[*noption] = pch;
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
                cBALLS_FAIL(cmd, "\nscanOption: noption = %d %s",
                            *noption,
                            "must be equal to number of infiles\n\n");
        if (*noption > MAXITEMS)
            cBALLS_FAIL(cmd,
                        "\nscaniOption: noption = %d %s",
                              *noption,
                              "must be less than the maximum num. of lines\n\n");

        for (i=0; i<*noption; i++) {
            if (parse_int_checked(poptionstr[i], &option[i],
                                  cmd->error_message, _ERRORMSGSIZE_,
                                  message) == FAILURE)
                return FAILURE;
            verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                                   "%s: option: %d\n",
                                routineName, option[i]);
        }
        verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                               "%s: noptions, nfiles: %d %d\n\n",
                            routineName, *noption,nfiles);
        if (flag == 1) {
            if (*noption > nfiles)
                cBALLS_FAIL(cmd,
                            "\nscanOption: noption = %d %s",
                            *noption,
                            "must be less or equal to number of files\n\n");
            else {
                for (i=*noption; i<nfiles; i++) {
                    option[i]=option[i-1]+1;
                    verb_print_debug_info(cmd->verbose, cmd->verbose_log,
                                          gd->outlog, "%s: option: %d\n",
                                          routineName, option[i]);
                }
            }
        }
    } else {
        //B ask a question if optionstr is columns then fix columns to default values..
        if (scanopt(message, "columns")) {
#if defined(ADDONS) && defined(IOLIB)
            gd->columns[0] = 1;
            gd->columns[1] = 2;
            gd->columns[2] = 3;
            gd->columns[3] = 4;
            gd->columns[4] = 5;
            gd->columns[5] = 6;
            gd->columns[6] = 7;
#endif
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
    char *poptionstr[MAXITEMS], optiontmp[MAXLENGTHOFSTRSCMD];
    int i;

    verb_log_print(cmd->verbose_log, gd->outlog,
                           "\n%s: Processing '%s' option:\n",
                        routineName, message);

    verb_log_print(cmd->verbose_log, gd->outlog,
                           "%s: Splitting string \"%s\" in tokens:\n",
                        routineName, message);

    if (!strnull(optionstr)) {
        if (copy_checked(optiontmp, sizeof(optiontmp),
                         optionstr, message) != 0)
            cBALLS_FAIL(cmd, "%s: option '%s' too long\n",
                        routineName, message);
        verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                               "\n%s: Splitting string \"%s\" in tokens:\n",
                            routineName, optiontmp);
        *noption=0;
        pch = strtok(optiontmp," ,");
        while (pch != NULL) {
            if (*noption >= MAXITEMS)
                cBALLS_FAIL(cmd, "\n%s: noption = %d must be less than "
                            "the maximum num. of lines\n\n",
                            routineName, *noption + 1);
            poptionstr[*noption] = pch;
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
                cBALLS_FAIL(cmd,
                            "\nscanOption: noption = %d must be equal to number of files\n\n",
                                      *noption);
        if (*noption > MAXITEMS)
            cBALLS_FAIL(cmd, "\nscanOption: noption = %d must be less than the maximum num. of lines\n\n",
                        *noption);

        for (i=0; i<*noption; i++) {
            if (parse_double_checked(poptionstr[i], &option[i],
                                     cmd->error_message, _ERRORMSGSIZE_,
                                     message) == FAILURE)
                return FAILURE;
            if (cmd->verbose_log>=3)
            verb_log_print(cmd->verbose_log, gd->outlog, "option: %g\n",option[i]);
        }

        verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                              "\n%s: noptions, nfiles: %d %d\n",
                            routineName, *noption, nfiles);
        if (flag == 1) {
            if (*noption > nfiles)
                cBALLS_FAIL(cmd,
                            "\nscanOption: noption = %d must be less or equal to number of files\n\n",
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
    typedef struct {
        const char *name;
        const char *value;
    } make_setting;

#define MAKE_SETTING(name) {#name, CBALLS_MAKE_##name}
    static const make_setting code_settings[] = {
        MAKE_SETTING(DEFDIMENSION),
        MAKE_SETTING(USEGSL),
        MAKE_SETTING(GSLINTERNAL),
        MAKE_SETTING(OPENMPMACHINE),
        MAKE_SETTING(ADDONSON),
    };
    static const make_setting machine_settings[] = {
        MAKE_SETTING(SINGLEPON),
        MAKE_SETTING(LONGINTON),
        MAKE_SETTING(MACONLYON),
        MAKE_SETTING(DEBUGON),
        MAKE_SETTING(DEBUGCOMPILINGON),
        MAKE_SETTING(DEBUGTREEON),
        MAKE_SETTING(DEBUGTRACKINGON),
        MAKE_SETTING(GENERALLIBS_ON),
        MAKE_SETTING(GETPARAM_ON),
    };
    static const make_setting addon_settings[] = {
        MAKE_SETTING(OCTREEGGGOMPON),
        MAKE_SETTING(KDTREEOMPON),
        MAKE_SETTING(BALLTREEOMPON),
        MAKE_SETTING(BALLTREE2BALLSOMPON),
        MAKE_SETTING(BALLTREE2BALLSMPION),
        MAKE_SETTING(OCTREE2BALLSOMPON),
        MAKE_SETTING(OCTREE2BALLSMPION),
        MAKE_SETTING(BALLTREE2BALLSOMP3PCFON),
        MAKE_SETTING(BALLTREE2BALLSMPI3PCFON),
        MAKE_SETTING(BALLTREEMPION),
        MAKE_SETTING(OCTREEGGGMPION),
        MAKE_SETTING(GGG_OMP_PIVOT_CHUNK_SIZE),
        MAKE_SETTING(GGG_MPI_PIVOT_CLAIM_SIZE),
        MAKE_SETTING(LYAFORESTOMPON),
        MAKE_SETTING(LYAFORESTMPION),
        MAKE_SETTING(LYA1D_OMP_PIVOT_BLOCK_SIZE),
        MAKE_SETTING(OCTREESHEAROMPON),
        MAKE_SETTING(KDTREEBOXOMPON),
        MAKE_SETTING(NEIGHBORBOXESOMPON),
        MAKE_SETTING(GADGETIOON),
        MAKE_SETTING(CLASSLIBON),
        MAKE_SETTING(PXDON),
        MAKE_SETTING(IOLIBON),
        MAKE_SETTING(CFITSIOLIBON),
        MAKE_SETTING(CFITSIOON),
        MAKE_SETTING(SMOOTHPIVOTON),
        MAKE_SETTING(NOWKAvgON),
        MAKE_SETTING(NMultipolesON),
        MAKE_SETTING(NONORMHISTON),
        MAKE_SETTING(POLARAXISON),
        MAKE_SETTING(NOLIMBERON),
        MAKE_SETTING(OVERCOUNTINGON),
        MAKE_SETTING(NOSTANDARNORMHISTON),
        MAKE_SETTING(CMDLINE_DEFS_UNITSPHERE_ON),
        MAKE_SETTING(TWOPCFON),
        MAKE_SETTING(TPCFON),
        MAKE_SETTING(TPCFSHEARON),
        MAKE_SETTING(PRUNEON),
        MAKE_SETTING(BALLS4SCANLEVON),
        MAKE_SETTING(THETA),
        MAKE_SETTING(NORMALHISTSCALEON),
        MAKE_SETTING(_DEBUGON_),
        MAKE_SETTING(LOGBINON),
        MAKE_SETTING(ORIGINALCBON),
        MAKE_SETTING(ADDPIVOTNEIGHBOURSON),
        MAKE_SETTING(OCTREE3PCF3DOMPON),
        MAKE_SETTING(OCTREE3PCF3DMPION),
        MAKE_SETTING(CB3D_OMP_PIVOT_BLOCK_SIZE),
        MAKE_SETTING(OCTREESMOOTHINGON),
        MAKE_SETTING(OCTREEKKKOMPON),
        MAKE_SETTING(BALLSON),
        MAKE_SETTING(BALLS0357ON),
        MAKE_SETTING(ALLOW_BALLS_WITH_SMOOTHING),
        MAKE_SETTING(OCTREEBALLS4OMPON),
        MAKE_SETTING(OCTREEBALLS4MPION),
        MAKE_SETTING(OCTREEGGGCROSSOMPON),
        MAKE_SETTING(DIRECTMETHODON),
        MAKE_SETTING(OCTREEGGGON),
        MAKE_SETTING(DIRECTMETHODSIMPLEON),
        MAKE_SETTING(DIRECTMETHODSIMPLELOOPIDON),
        MAKE_SETTING(COMPLEXLIBON),
        MAKE_SETTING(OCTREEGGGOMPTRIANGLESON),
        MAKE_SETTING(OCTREEBOXOMPON),
        MAKE_SETTING(KDTREECUTEBOXON),
        MAKE_SETTING(OCTREEKKKBALLS4OMPTRIANGLESON),
        MAKE_SETTING(COSMOLIBON),
        MAKE_SETTING(SAVERESTOREON),
        MAKE_SETTING(TREEOMPSINCOSON),
        MAKE_SETTING(OCTREESINCOSOMPON),
    };
    static const make_setting toolchain_settings[] = {
        MAKE_SETTING(MDIR),
        MAKE_SETTING(WRKDIR),
        MAKE_SETTING(MACHINES_DIR),
        MAKE_SETTING(CC),
        MAKE_SETTING(PYTHON),
        MAKE_SETTING(OPTFLAG),
        MAKE_SETTING(OMPFLAG),
        MAKE_SETTING(OPT1),
        MAKE_SETTING(OPT2),
        MAKE_SETTING(OPT2LIB),
        MAKE_SETTING(CCFLAG),
        MAKE_SETTING(LDFLAG),
        MAKE_SETTING(INCLUDES),
        MAKE_SETTING(MLIBS),
        MAKE_SETTING(FITSIOLIBS),
        MAKE_SETTING(EXTERNAL),
        MAKE_SETTING(AR),
        MAKE_SETTING(SANITIZER),
        MAKE_SETTING(MACOSX_DEPLOYMENT_TARGET),
        MAKE_SETTING(PKG_CONFIG),
        MAKE_SETTING(CFITSIO_PKG),
        MAKE_SETTING(CFITSIO_INCLUDE),
        MAKE_SETTING(CFITSIO_LIB),
        MAKE_SETTING(CFITSIO_INCL),
        MAKE_SETTING(CFITSIO_LDFLAGS),
        MAKE_SETTING(CFITSIO_LIBS),
        MAKE_SETTING(CFITSIO_RPATH),
        MAKE_SETTING(GSL_CONFIG),
        MAKE_SETTING(GSL_INCLUDE),
        MAKE_SETTING(GSL_LIB),
        MAKE_SETTING(GSL_CFLAGS),
        MAKE_SETTING(GSL_LDFLAGS),
        MAKE_SETTING(GSL_LIBS),
        MAKE_SETTING(BUILD_ADDONS),
        MAKE_SETTING(BUILD_CLASSLIB),
        MAKE_SETTING(BUILD_OPENMP),
        MAKE_SETTING(BUILD_MPI),
        MAKE_SETTING(BUILD_LYA),
        MAKE_SETTING(BUILD_SHEAR),
        MAKE_SETTING(BUILD_PRECISION),
        MAKE_SETTING(DEFINEFLAGS),
        MAKE_SETTING(SANITIZER_FLAGS),
        MAKE_SETTING(PROJECT_WARNING_FLAGS),
        MAKE_SETTING(DIMCODE),
        MAKE_SETTING(OMPCODE),
        MAKE_SETTING(SINGLEP),
        MAKE_SETTING(LONG),
        MAKE_SETTING(PYTHON_MACOSX_DEPLOYMENT_TARGET),
        MAKE_SETTING(OPENMP_MACOSX_DEPLOYMENT_TARGET),
        MAKE_SETTING(MACOS_DEPLOYMENT_FLAG),
        MAKE_SETTING(GENERALLIBS),
        MAKE_SETTING(GETPARAM),
        MAKE_SETTING(MPICC),
        MAKE_SETTING(UNAME_S),
    };
    const make_setting *groups[] = {
        code_settings,
        machine_settings,
        addon_settings,
        toolchain_settings,
    };
    const size_t group_sizes[] = {
        sizeof(code_settings) / sizeof(code_settings[0]),
        sizeof(machine_settings) / sizeof(machine_settings[0]),
        sizeof(addon_settings) / sizeof(addon_settings[0]),
        sizeof(toolchain_settings) / sizeof(toolchain_settings[0]),
    };
    const char *group_names[] = {
        "Makefile_settings",
        "Makefile_machine switches",
        "addons/Makefile_addons_settings",
        "Makefile_machine toolchain and resolved paths",
    };
    size_t group;
    size_t i;

#undef MAKE_SETTING

    (void) gd;
    verb_print_zero(cmd->verbose,
               "\nprint_make_info:\n\n");

    for (group = 0; group < sizeof(groups) / sizeof(groups[0]); group++) {
        verb_print_zero(cmd->verbose, "[%s]\n", group_names[group]);
        for (i = 0; i < group_sizes[group]; i++) {
            verb_print_zero(cmd->verbose, "  %-38s = %s\n",
                            groups[group][i].name,
                            groups[group][i].value);
        }
        verb_print_zero(cmd->verbose, "\n");
    }

    verb_print_zero(cmd->verbose,
                    "[effective compile-time features (enabled)]\n");

#ifdef ADDONS
    verb_print_zero(cmd->verbose, "with ADDONS\n");
#endif

#ifdef TWODIMCODE
    verb_print_zero(cmd->verbose, "TWODIMCODE\n");
#endif
#ifdef THREEDIMCODE
    verb_print_zero(cmd->verbose, "THREEDIMCODE\n");
#endif

#ifdef OPENMPCODE
    verb_print_zero(cmd->verbose, "using OpenMP\n");
#endif

    verb_print_zero(cmd->verbose, "precision profile: %s\n",
                    CballsStoragePrecision);

#ifdef LONGINT
    verb_print_zero(cmd->verbose, "LONGINT\n");
#endif

//B engines
#ifdef BALLS
    verb_print_zero(cmd->verbose, "with BALLS engine\n");
#endif

#ifdef OCTREESMOOTHING
    verb_print_zero(cmd->verbose, "with OCTREESMOOTHING\n");
#endif

#ifdef SINCOS
    verb_print_zero(cmd->verbose, "with SINCOS\n");
#endif

#ifdef TREENODEALLBODIES
    verb_print_zero(cmd->verbose, "with TREENODEALLBODIES\n");
#endif

#ifdef TREENODEBALLS4
    verb_print_zero(cmd->verbose, "with TREENODEBALLS4 engine\n");
#endif

#ifdef OCTREEGGGOMP
    verb_print_zero(cmd->verbose, "with OCTREEGGGOMP engine\n");
#endif

#ifdef NMultipoles
    verb_print_zero(cmd->verbose, "with NMultipoles\n");
#else
    verb_print_zero(cmd->verbose, "without NMultipoles\n");
#endif
#ifdef NONORMHIST
    verb_print_zero(cmd->verbose, "with NONORMHIST\n");
#else
    verb_print_zero(cmd->verbose, "without NONORMHIST\n");
#endif

#ifdef KDTREEOMP
    verb_print_zero(cmd->verbose, "with KDTREEOMP engine\n");
#endif

#ifdef BALLTREEOMP
    verb_print_zero(cmd->verbose, "with BALLTREEOMP (FCFC-style) engine\n");
#endif

#ifdef BALLTREE2BALLSOMP
    verb_print_zero(cmd->verbose,
                    "with BALLTREE2BALLSOMP (TreeCorr-style dual/triple-node) engine\n");
#endif

#ifdef BALLTREE2BALLSMPI
    verb_print_zero(cmd->verbose,
                    "with BALLTREE2BALLSMPI (distributed dual/triple-node) engine\n");
#endif

#ifdef OCTREE2BALLSOMP
    verb_print_zero(cmd->verbose,
                    "with OCTREE2BALLSOMP (TreeCorr scan over native octree) engine\n");
#endif

#ifdef OCTREE2BALLSMPI
    verb_print_zero(cmd->verbose,
                    "with OCTREE2BALLSMPI (distributed LogMultipole) engine\n");
#endif
#ifdef OCTREEBALLS4OMP
    verb_print_zero(cmd->verbose, "with OCTREEBALLS4OMP (masked B4/complex edge) engine\n");
#endif
#ifdef OCTREEBALLS4MPI
    verb_print_zero(cmd->verbose, "with OCTREEBALLS4MPI (distributed B4/complex edge) engine\n");
#endif

#ifdef BALLTREE2BALLSOMP3PCF
    verb_print_zero(cmd->verbose,
                    "with BALLTREE2BALLSOMP3PCF (LogMultipole frontier, 3PCF only) engine\n");
#endif

#ifdef BALLTREE2BALLSMPI3PCF
    verb_print_zero(cmd->verbose,
                    "with BALLTREE2BALLSMPI3PCF (distributed FCFC ball-tree LogMultipole) engine\n");
#endif

#ifdef BALLTREEMPI
    verb_print_zero(cmd->verbose,
                    "with BALLTREEMPI (FCFC-style MPI) engine\n");
#endif
#ifdef OCTREEGGGMPI
    verb_print_zero(cmd->verbose,
                    "with OCTREEGGGMPI (FCFC-style MPI) engine\n");
#endif

#ifdef OCTREEKKKOMP
    verb_print_zero(cmd->verbose, "with OCTREEKKKOMP engine\n");
#endif

#ifdef LYAFORESTOMP
    verb_print_zero(cmd->verbose, "with LYAFORESTOMP (seven radial/3D forest engines)\n");
#endif
#ifdef LYAFORESTMPI
    verb_print_zero(cmd->verbose, "with LYAFORESTMPI (seven MPI+OpenMP forest engines)\n");
#endif
#ifdef OCTREE3PCF3DOMP
    verb_print_zero(cmd->verbose, "with OCTREE3PCF3DOMP (scalar/survey Legendre estimator)\n");
#endif
#ifdef OCTREE3PCF3DMPI
    verb_print_zero(cmd->verbose, "with OCTREE3PCF3DMPI (distributed scalar/survey estimator)\n");
#endif
#ifdef OCTREESHEAROMP
    verb_print_zero(cmd->verbose, "with OCTREESHEAROMP (flat-sky spin-2 estimator)\n");
#endif

#ifdef DIRECTMETHOD
        verb_print_zero(cmd->verbose, "with DIRECTMETHOD engine\n");
#endif

#ifdef DIRECTMETHODSIMPLE
        verb_print_zero(cmd->verbose, "with DIRECTMETHODSIMPLE engine\n");
#endif
//E

//B correlations
#ifdef TWOPCF
    verb_print_zero(cmd->verbose,
                    "with 2pcf convergence computation\n");
#endif

#ifdef THREEPCFCONVERGENCE
    verb_print_zero(cmd->verbose,
                    "with 3pcf convergence computation\n");
#endif

#ifdef THREEPCFSHEAR
    verb_print_zero(cmd->verbose,
                    "with 3pcf shear computation \n");
#endif
//E

//B speed up
#ifdef SMOOTHPIVOT
        verb_print_zero(cmd->verbose,
                            "with SMOOTHPIVOT compiled; enable options=smooth-pivot only on supported engines\n");
#endif

#ifdef BALLS4SCANLEV
    verb_print_zero(cmd->verbose,
                        "with BALLS4SCANLEV\n");
#endif
//E

//B Cython interface
#ifdef CLASSLIB
    verb_print_zero(cmd->verbose, "with CLASSLIB for cyballs module\n");
#endif

#ifdef PXD
    verb_print_zero(cmd->verbose, "with PXD Cython getters\n");
#endif
//E

//B I/O modules
#ifdef IOLIB
    verb_print_zero(cmd->verbose, "with IOLIB\n");
#endif

#ifdef CFITSIO
    verb_print_zero(cmd->verbose, "with CFITSIO\n");
#endif

#ifdef GADGETIO
    verb_print_zero(cmd->verbose, "with GADGETIO\n");
#endif
//E

#ifdef POLARAXIS
    verb_print_zero(cmd->verbose, "with POLARAXIS\n");
#endif

#ifdef NOLIMBER
    verb_print_zero(cmd->verbose, "with NOLIMBER\n");
#endif

#ifdef NOWKAvg
    verb_print_zero(cmd->verbose, "NOWKAvg\n");
#endif

#ifdef PTOPIVOTROTATION
    verb_print_zero(cmd->verbose, "PTOPIVOTROTATION\n");
#endif

#ifdef PIVOTEXTERNAL
    verb_print_zero(cmd->verbose, "with PIVOTEXTERNAL\n");
#endif

#ifdef ADDPIVOTNEIGHBOURS
        verb_print_zero(cmd->verbose, "with ADDPIVOTNEIGHBOURS\n");
#endif

#ifdef OVERCOUNTING
        verb_print_zero(cmd->verbose, "with OVERCOUNTING\n");
#endif

#ifdef SAVERESTORE
        verb_print_zero(cmd->verbose, "with SAVERESTORE\n");
#endif

//B libraries
#ifdef GETPARAM
    verb_print_zero(cmd->verbose, "using GETPARAM\n");
#endif

#ifdef USEGSL
#ifdef GSLINTERNAL
    verb_print_zero(cmd->verbose, "using internal GSL\n");
#else
    verb_print_zero(cmd->verbose, "using GSL\n");
#endif
#endif
//E

#ifdef DEBUG
    verb_print_zero(cmd->verbose, "DEBUG\n");
#endif

    verb_print_zero(cmd->verbose,
               "\n");


    return SUCCESS;
}

local int print_search_methods(struct cmdline_data* cmd,
                               struct global_data* gd)
{
    typedef struct {
        const char *name;
        const char *geometry;
        const char *correlations;
        const char *usage;
    } search_method_help;

    static const search_method_help known_methods[] = {
        {"octree-sincos-omp", "2D/3D scalar; octree; OpenMP",
         "standard 2PCF and sine/cosine 3PCF multipoles",
         "Core method. Select KKKCorrelation for convergence 3PCF output."},
        {"balls-omp", "scalar; octree ball aggregation; OpenMP",
         "2PCF/3PCF multipoles",
         "Use no-one-ball for exact body traversal; theta controls cell acceptance."},
        {"kdtree-omp", "scalar; balanced k-d tree; OpenMP",
         "2PCF/3PCF multipoles",
         "Exact traversal is the default. Add behavior-ball for aggregate-node acceptance; no-one-ball disables it; aggregation requires useLogHist=true. no-normalize-HistZeta selects weighted raw distinct-triplet multipoles; smooth-pivot must be disabled for this mode."},
        {"balltree-omp", "scalar; FCFC-style PCA ball tree; OpenMP",
         "2PCF/3PCF multipoles",
         "nsmooth sets leaf capacity. Exact traversal is the default; behavior-ball enables aggregate-node acceptance, no-one-ball disables it, and aggregation requires useLogHist=true. no-normalize-HistZeta selects weighted raw distinct-triplet multipoles; smooth-pivot must be disabled for this mode."},
        {"balltree-2balls-omp", "scalar; FCFC PCA ball tree; TreeCorr dual/triple-node OpenMP scan",
         "2PCF and angular-multipole 3PCF with auto- and cross-catalog support",
         "TWOPCFON and TPCFON switch 2PCF and 3PCF independently at build time; only-2pcf and only-3pcf select one at runtime when both are built. The direct distinct-triplet 3PCF has cubic worst-case work and is not intended for full-sky NSIDE=1024 catalogs. TreeCorr-style cell aggregation is the default; nsmooth sets leaf capacity, theta scales radial and angular tolerances, weights-norm uses the distinct-triplet weight denominator, no-normalize-HistZeta returns raw sums, and no-two-balls forces exact body pairs and triplets."},
        {"balltree-2balls-mpi", "scalar; FCFC PCA ball tree; deterministic MPI+OpenMP dual/triple-node scan",
         "distributed 2PCF and angular-multipole 3PCF with auto- and cross-catalog support",
         "Build with BALLTREE2BALLSMPION=1 and run with mpiexec. TWOPCFON and TPCFON select the compiled orders; only-2pcf and only-3pcf select one at runtime. numberThreads sets OpenMP threads per rank. Runtime acceptance and normalization options match balltree-2balls-omp. The explicit triple-node 3PCF remains cubic even when distributed."},
        {"octree-2balls-omp", "scalar; native octree; dual-node 2PCF and LogMultipole 3PCF",
         "2PCF and angular-multipole 3PCF with auto- and cross-catalog support",
         "Build with OCTREE2BALLSOMPON=1. TWOPCFON and TPCFON switch correlation orders; only-2pcf and only-3pcf select one at runtime. The 2PCF uses TreeCorr-style two-ball acceptance. The 3PCF scans pivot-neighbor pairs, accumulates radial angular moments, removes q==r with second moments, and forms all radial-bin pairs without enumerating triples. theta controls radial and angular slop, no-two-balls selects exact body moments, and treecorr-direct-triples selects the slower cubic validation traversal. read-mask selects masked autocorrelation: use a FITS companion mask, one columns-ascii-all/binary-all catalog, or set_catalog(mask=...) in cyballs. Only valid bodies enter pivots, neighbors, and normalization; an empty mask is an error."},
        {"octree-2balls-mpi", "scalar; native octree; deterministic MPI+OpenMP frontier",
         "distributed 2PCF and LogMultipole angular-multipole 3PCF",
         "Build with OCTREE2BALLSMPION=1 and run with mpiexec. numberThreads sets OpenMP threads per rank. It accepts the octree-2balls-omp runtime options, including read-mask with valid-body-only normalization, except treecorr-direct-triples, which remains a serial validation oracle. All ranks build the same masked tree and reduce its frontier tasks deterministically."},
        {"balltree-2balls-omp_3pcf", "scalar; FCFC PCA ball tree; TreeCorr LogMultipole pair factorization",
         "angular-multipole 3PCF with auto- and cross-catalog support",
         "TPCFON=1 is required. Pivot-neighbor tree scans accumulate radial field multipoles and form 3PCF bins from pair products, excluding q==r with node second moments. TreeCorr-style radial-bin slop and angular phase-error criteria control aggregation; completed large-radius moments are inherited while unresolved bins split the pivot. nsmooth sets the default leaf capacity; treecorr-singleton-leaves forces one body per leaf. theta=1 is the production setting, theta approaching zero converges to exact. Multipoles are normalized by distinct-triplet count by default; weights-norm uses the distinct-triplet weight denominator and no-normalize-HistZeta returns raw sums. no-two-balls selects exact body moments, and treecorr-direct-triples selects the slower validation traversal."},
        {"balltree-2balls-mpi_3pcf", "scalar; FCFC PCA ball tree; deterministic MPI+OpenMP LogMultipole frontier",
         "distributed angular-multipole 3PCF with auto- and cross-catalog support",
         "Build with BALLTREE2BALLSMPI3PCFON=1 and TPCFON=1, then run with mpiexec. numberThreads sets OpenMP threads per rank. It accepts the production balltree-2balls-omp_3pcf options except treecorr-direct-triples, which remains an OpenMP validation oracle."},
        {"balltree-mpi", "scalar; FCFC-style ball tree; MPI+OpenMP",
         "distributed 2PCF/3PCF multipoles",
         "Run with mpiexec; numberThreads sets OpenMP threads per MPI rank."},
        {"octree-kkk-omp", "scalar convergence; octree; OpenMP",
         "KK, KKK, NN, and configured count estimators",
         "Select the correlation with KKCorrelation, KKKCorrelation, NNCorrelation, or NNEstimator."},
        {"octree-ggg-omp", "projected scalar; octree; OpenMP",
         "2PCF and angular-multipole 3PCF",
         "Common options include KKKCorrelation, compute-HistN, read-mask, and edge-corrections. Add only-2pcf to skip all angular multipoles and 3PCF storage. Window modes extend through mChebyshev normally, or 2*mChebyshev with edge-corrections or ggg-full-window. ggg-profile reports work and ordered-reduction timings."},
        {"octree-ggg-mpi", "projected scalar; octree; MPI+OpenMP",
         "distributed octree-ggg-omp estimator",
         "Run with mpiexec; numberThreads sets OpenMP threads per MPI rank."},
        {"octree-shear-omp", "flat-sky spin-2 shear; exact octree leaves; OpenMP",
         "xi+/xi- and four Gamma-x shear 3PCF components",
         "Use Cartesian tangent-plane catalogs with pos-and-shear; periodic geometry is not supported."},
        {"octree-sincos-omp-addons", "scalar; addon octree; OpenMP",
         "sine/cosine 3PCF multipoles",
         "Addon variant of the core octree-sincos method."},
        {"tree-omp-sincos", "scalar; legacy tree; OpenMP",
         "sine/cosine 3PCF multipoles",
         "Legacy comparison method using the common histogram parameters."},
        {"octree-ggg-cross-omp", "projected scalar cross-catalog; octree; OpenMP",
         "cross-catalog angular-multipole 3PCF",
         "Provide the required distinct catalogs through infile and iCatalogs."},
        {"direct-sincos", "scalar; direct enumeration",
         "reference sine/cosine correlations",
         "Brute-force validation method; intended for small catalogs."},
        {"kdtree-cute-box", "periodic Cartesian box; k-d tree",
         "CUTE-compatible box correlations",
         "Use usePeriodic=true and the cute-box/cute-box-fmt options as needed."},
        {"direct-simple-sincos", "scalar; direct simple loops",
         "reference sine/cosine correlations",
         "Development/reference method for small catalogs."},
        {"octree-ggg-omp-triangles", "projected scalar; octree; OpenMP",
         "explicit triangle-oriented 3PCF output",
         "Triangle-output variant of octree-ggg-omp."},
        {"octree-balls4-omp", "scalar convergence; normal octree; OpenMP",
         "weighted convergence 2PCF and balls4 KKK 3PCF multipoles",
         "Requires 3D, usePeriodic=false, and useLogHist=true. Enable OCTREEBALLS4OMPON=1. Build with TWOPCFON=1 for 2PCF and TPCFON=1 for 3PCF; compute-HistN,and-CF writes pair-count/CF files. Works with SMOOTHPIVOTON=0 or 1; options=smooth-pivot is not supported. read-mask accepts an embedded columns-ascii-all/binary-all mask or an in-memory mask. edge-corrections,no-normalize-HistZeta uses distinct-neighbor body-pivot multipoles, window modes through 2*mChebyshev, and a complex solve; weights-norm enables weighted windows. Either no-one-ball or no-two-balls makes edge moments exact."},
        {"octree-balls4-mpi", "scalar convergence; normal octree; MPI+OpenMP",
         "distributed BALLS4 2PCF and 3PCF with masks and complex edge corrections",
         "Enable OCTREEBALLS4MPION=1 and run with mpiexec; numberThreads is threads per rank. The runtime options match octree-balls4-omp, including read-mask and edge-corrections. Ranks own disjoint B4 work; raw sums are reduced before normalization and only rank 0 writes output. Edge pivot blocks are reduced in a fixed order. SMOOTHPIVOT is unsupported."},
        {"octree-kkk-balls4-omp-triangles", "scalar convergence; octree; OpenMP",
         "balls4 triangle-oriented KKK 3PCF",
         "Triangle-output variant of octree-balls4-omp."},
        {"kdtree-box-omp", "periodic Cartesian box; k-d tree; OpenMP",
         "box 2PCF/3PCF",
         "Use usePeriodic=true; cute-box-fmt writes CUTE-compatible output."},
        {"octree-box-omp", "periodic Cartesian box; octree; OpenMP",
         "box 2PCF/3PCF",
         "Use usePeriodic=true and configure rangeN/rminHist/sizeHistN."},
        {"neighbor-boxes-omp", "periodic Cartesian box; linked boxes; OpenMP",
         "local-neighbor box 2PCF/3PCF",
         "Use usePeriodic=true; add only-pos for catalogs without a convergence field."},
        {"octree-ggg", "projected scalar; serial octree",
         "legacy angular-multipole 3PCF",
         "Serial comparison variant of octree-ggg-omp."},
        {"direct-simple-sincos-loopId", "scalar; direct loops with IDs",
         "reference sine/cosine correlations",
         "Development/reference method that retains loop identifiers."},
        {"balls-omp-0357", "scalar; legacy balls profile; OpenMP",
         "legacy 2PCF/3PCF multipoles",
         "Optional compatibility profile; prefer balls-omp for new runs."},
        {"octree-3pcf-3d-omp", "3D scalar; exact octree leaves; OpenMP",
         "spherical-harmonic 2PCF/3PCF multipoles",
         "Use x,y,z,delta,weight input. compute-2pcf-3d and compute-3pcf-3d select outputs; octree-ggg-3d-omp is an alias."},
        {"octree-3pcf-3d-mpi", "3D scalar; exact octree leaves; MPI+OpenMP pivot blocks",
         "spherical-harmonic 2PCF/3PCF; data/random survey estimator and edge correction",
         "Enable OCTREE3PCF3DMPION=1 and run with mpiexec. compute-2pcf-3d and compute-3pcf-3d select outputs. survey-estimator-3d uses data,random catalogs with iCatalogs=1,2. Catalogs are replicated; numberThreads is per rank. octree-ggg-3d-mpi is an alias."},
        {"lya-2pcf-omp", "3D Lyman-alpha forest pixels; exact octree leaves; OpenMP",
         "weighted anisotropic 2PCF",
         "Requires one forest catalog, non-periodic geometry, and forest IDs; see input below."},
        {"lya-3pcf-omp", "3D Lyman-alpha forest pixels; exact octree leaves; OpenMP",
         "weighted five-dimensional 3PCF",
         "Requires one forest catalog, non-periodic geometry, and forest IDs; see input below."},
        {"lya-2pcf-3pcf-omp", "3D Lyman-alpha forest pixels; exact octree leaves; OpenMP",
         "weighted 2PCF and 3PCF in one traversal",
         "Requires one forest catalog, non-periodic geometry, and forest IDs; see input below."},
        {"lya-1d-2pcf-omp", "radial Lyman-alpha pixels; sorted 1D range scan; OpenMP",
         "weighted radial-only 2PCF",
         "Ignores transverse separation; lya2RpMax and lya2RpBins define |delta radial distance| bins."},
        {"lya-1d-3pcf-omp", "radial Lyman-alpha pixels; sorted 1D range scan; OpenMP",
         "weighted radial-only 3PCF",
         "Ignores transverse separation; lya3RMax and lya3RBins define signed pivot-lag bins."},
        {"lya-1d-2pcf-3pcf-omp", "radial Lyman-alpha pixels; sorted 1D range scan; OpenMP",
         "weighted radial-only 2PCF and 3PCF",
         "Requires forest IDs; all correlated pixels must belong to distinct quasars."},
        {"lya-1d-tree-2pcf-omp", "radial Lyman-alpha pixels; exact 1D interval tree; OpenMP",
         "weighted radial-only 2PCF",
         "Uses exact same-bin node aggregation and subtracts within-quasar pairs; lya2RpMax and lya2RpBins set the domain."},
        {"lya-2pcf-mpi", "Lyman-alpha pixels; 3D octree; MPI+OpenMP",
         "weighted anisotropic 2PCF",
         "Enable LYAFORESTMPION=1; run with mpiexec and numberThreads per rank. Requires one forest catalog (file or in-memory, see below), DEFDIMENSION=3, and usePeriodic=false. Options and output columns match lya-2pcf-omp. Raw sums are reduced before normalization; only rank 0 writes."},
        {"lya-3pcf-mpi", "Lyman-alpha pixels; 3D octree; MPI+OpenMP",
         "weighted five-dimensional 3PCF",
         "Enable LYAFORESTMPION=1; run with mpiexec and numberThreads per rank. Requires one forest catalog (file or in-memory, see below), DEFDIMENSION=3, and usePeriodic=false. Options and output columns match lya-3pcf-omp. Raw sums are reduced before normalization; only rank 0 writes."},
        {"lya-2pcf-3pcf-mpi", "Lyman-alpha pixels; 3D octree; MPI+OpenMP",
         "weighted 2PCF and 3PCF",
         "Enable LYAFORESTMPION=1; run with mpiexec and numberThreads per rank. Requires one forest catalog (file or in-memory, see below), DEFDIMENSION=3, and usePeriodic=false. Options and output columns match lya-2pcf-3pcf-omp. Raw sums are reduced before normalization; only rank 0 writes."},
        {"lya-1d-2pcf-mpi", "Lyman-alpha pixels; radial range scan; MPI+OpenMP",
         "weighted radial-only 2PCF",
         "Enable LYAFORESTMPION=1; run with mpiexec and numberThreads per rank. Requires one forest catalog (file or in-memory, see below), DEFDIMENSION=3, and usePeriodic=false. Options and output columns match lya-1d-2pcf-omp. Raw sums are reduced before normalization; only rank 0 writes."},
        {"lya-1d-3pcf-mpi", "Lyman-alpha pixels; radial range scan; MPI+OpenMP",
         "weighted radial-only 3PCF",
         "Enable LYAFORESTMPION=1; run with mpiexec and numberThreads per rank. Requires one forest catalog (file or in-memory, see below), DEFDIMENSION=3, and usePeriodic=false. Options and output columns match lya-1d-3pcf-omp. Raw sums are reduced before normalization; only rank 0 writes."},
        {"lya-1d-2pcf-3pcf-mpi", "Lyman-alpha pixels; radial range scan; MPI+OpenMP",
         "weighted radial-only 2PCF and 3PCF",
         "Enable LYAFORESTMPION=1; run with mpiexec and numberThreads per rank. Requires one forest catalog (file or in-memory, see below), DEFDIMENSION=3, and usePeriodic=false. Options and output columns match lya-1d-2pcf-3pcf-omp. Raw sums are reduced before normalization; only rank 0 writes."},
        {"lya-1d-tree-2pcf-mpi", "Lyman-alpha pixels; radial interval tree; MPI+OpenMP",
         "weighted radial-only 2PCF with same-forest subtraction",
         "Enable LYAFORESTMPION=1; run with mpiexec and numberThreads per rank. Requires one forest catalog (file or in-memory, see below), DEFDIMENSION=3, and usePeriodic=false. Options and output columns match lya-1d-tree-2pcf-omp. Raw sums are reduced before normalization; only rank 0 writes."}
    };
    const size_t method_count =
        sizeof(known_methods) / sizeof(known_methods[0]);
    size_t available_count = 0;
    size_t i;
    int method_id;

    (void)gd;

    for (i = 0; i < method_count; ++i) {
        search_method_string_to_int((string)known_methods[i].name, &method_id);
        if (method_id >= 0) ++available_count;
    }

    verb_print_zero(cmd->verbose,
                    "\nSearching methods registered in this executable (%zu):\n",
                    available_count);
    verb_print_zero(cmd->verbose,
                    "Select one with search=<name> or searchMethod=<name>.\n");
    verb_print_zero(cmd->verbose,
                    "Common controls: rangeN, rminHist, sizeHistN, useLogHist, "
                    "theta, numberThreads, and options.\n");
    verb_print_zero(cmd->verbose,
                    "Only methods enabled by this executable's build profile are shown.\n");

    for (i = 0; i < method_count; ++i) {
        search_method_string_to_int((string)known_methods[i].name, &method_id);
        if (method_id < 0) continue;
        verb_print_zero(cmd->verbose, "\n- %s (id=%d)\n",
                        known_methods[i].name, method_id);
        verb_print_zero(cmd->verbose, "  geometry: %s\n",
                        known_methods[i].geometry);
        verb_print_zero(cmd->verbose, "  computes: %s\n",
                        known_methods[i].correlations);
        verb_print_zero(cmd->verbose, "  use: %s\n",
                        known_methods[i].usage);
        if (strncmp(known_methods[i].name, "lya-", 4) == 0)
            verb_print_zero(cmd->verbose,
                "  input: one x y z delta weight forest_id catalog (lya-ascii), or "
                "cyballs.set_forest_catalog(positions, delta, weights, forest_ids). "
                "DEFDIMENSION=3, observer-centered comoving coordinates and "
                "usePeriodic=false are required even for radial searches. "
                "Pairs exclude the same quasar; triplets require three distinct "
                "quasars. Histograms are weight-normalized; empty bins are zero. "
                "Use python/lya_corr_all_engines.py for DESI FITS/NPZ/ASCII, "
                "one-time loading and MPI broadcasting. No smooth-pivot.\n");
        if (strncmp(known_methods[i].name, "octree-3pcf-3d-", 15) == 0)
            verb_print_zero(cmd->verbose,
                "  modes: only-2pcf-3d, only-3pcf-3d, or "
                "compute-2pcf-3d,compute-3pcf-3d. survey-estimator-3d "
                "uses data/random catalogs and window correction. "
                "exclude-same-los excludes pivot LOS matches, not neighbor-neighbor "
                "LOS matches: this is not the three-distinct-forest estimator.\n");
        if (strstr(known_methods[i].name, "2balls") != NULL)
            verb_print_zero(cmd->verbose,
                "  edge correction: add edge-corrections,no-normalize-HistZeta "
                "for complex scalar 3PCF window deconvolution. Window modes "
                "extend through 2*mChebyshev. weights-norm weights both signal "
                "and window; empty/singular bins are zero. 2PCF is unchanged.\n");
    }

    verb_print_zero(cmd->verbose,
                    "\nScalar angular engines preserve the observer frame in 3D and use "
                    "tangent-plane angles with Euclidean chord bins. Coincident/radial/antipodal "
                    "legs are excluded from angular multipoles, not from ordinary pair counts.\n"
                    "Use options=make-info to inspect the build profile and "
                    "options=print-options for the full option list.\n\n");
    return SUCCESS;
}

local int print_options(struct cmdline_data* cmd,
                        struct global_data* gd)
{
    typedef struct {
        const char *name;
        const char *scope;
        const char *description;
    } option_help;

    static const option_help registered_options[] = {
        {"all", "test models",
         "show selected and requested body counts for random unit-sphere data"},
        {"all-in-one", "input",
         "merge all input catalogs into one catalog before the search"},
        {"and-CF", "histograms",
         "also compute correlation functions; requires compute-HistN"},
        {"arfken", "coordinates",
         "use the Arfken spherical-coordinate convention (the default)"},
        {"asymmetric", "search",
         "use asymmetric pivot/neighbor traversal for cross-catalog searches"},
        {"behavior-ball", "tree search",
         "select the ball-oriented tree acceptance behavior"},
        {"behavior-tree-omp", "tree search addon",
         "select the OpenMP tree behavior in KKK/balls searches"},
        {"bh86", "tree search",
         "use the Barnes-Hut 1986 cell acceptance criterion"},
        {"celestial", "coordinates",
         "interpret angular input as celestial/equatorial coordinates"},
        {"center-of-mass", "tree search",
         "place aggregate cells at their center of mass"},
        {"check-eq-pos", "input",
         "reject input catalogs containing bodies at identical positions"},
        {"compute-HistN", "histograms",
         "compute pair-count histograms used for normalization and 2PCF output"},
        {"compute-2pcf-3d", "3D scalar addon",
         "enable the scalar three-dimensional 2PCF"},
        {"compute-3pcf-3d", "3D scalar addon",
         "enable the spherical-harmonic three-dimensional 3PCF"},
        {"compute-j-no-eq-i", "legacy/development",
         "reserved legacy option; no effect in the production search kernels"},
        {"pivot-loop", "legacy/development",
         "reserved legacy option; its old KKK branch is disabled"},
        {"cute-box", "box output",
         "apply the CUTE-compatible periodic-box coordinate convention"},
        {"cute-box-fmt", "box output addon",
         "write correlation output in CUTE box format"},
        {"cute-box-rmin", "box output",
         "use the CUTE-compatible radial-bin coordinate starting at rminHist"},
        {"default-rsmooth", "smoothing addon",
         "derive smoothing radii from the normal default rule"},
        {"ecliptic", "coordinates",
         "interpret angular input as ecliptic coordinates"},
        {"edge-corrections", "multipoles",
         "apply edge corrections to unnormalized multipole histograms"},
        {"edge-corrections-from-files", "multipoles",
         "compute edge corrections from previously generated histogram files"},
        {"edge-effects", "octree KKK addon",
         "enable edge-effect count normalization in the KKK search"},
        {"ggg-full-window", "octree GGG addon",
         "retain window modes through 2*mChebyshev for diagnostics or later edge correction"},
        {"ggg-profile", "octree GGG addon",
         "report per-rank work, ordered-wait, merge, and wall timings"},
        {"exclude-los", "3D scalar addon",
         "alias for exclude-same-los"},
        {"exclude-pivot-los", "3D scalar addon",
         "alias for exclude-same-los"},
        {"exclude-same-los", "3D scalar addon",
         "exclude pivot-neighbor pairs sharing a line-of-sight identifier"},
        {"fits-type-file", "CFITSIO",
         "open FITS input as a generic FITS file instead of a data unit"},
        {"fix-center", "sky selection",
         "select a fixed angular center using lengthBox"},
        {"fix-rsmooth", "smoothing addon",
         "keep the configured smoothing radius instead of deriving it"},
        {"full-sky", "sky output",
         "use full-sky angular normalization when writing correlations"},
        {"galactic", "coordinates",
         "interpret angular input as galactic coordinates"},
        {"GGGCorrelation", "correlation",
         "select the shear three-point correlation output"},
        {"header-info", "input",
         "print input-file header and column information"},
        {"in-degrees", "I/O addon",
         "convert angular input columns from degrees to radians"},
        {"kappa", "CFITSIO output",
         "write the convergence field in FITS output"},
        {"kappa-constant", "input/test models",
         "replace input convergence values with a constant"},
        {"kappa-constant-one", "input/test models",
         "use one, rather than two, for kappa-constant"},
        {"lya-output-empty-bins", "Lyman-alpha addon",
         "write zero-denominator bins in radial and five-dimensional 3PCF output; supported by OpenMP and MPI"},
        {"compute-2pcf-3d", "3D scalar OpenMP/MPI",
         "enable the scalar 2PCF; combine with compute-3pcf-3d for both orders"},
        {"compute-3pcf-3d", "3D scalar OpenMP/MPI",
         "enable spherical-harmonic 3PCF multipoles"},
        {"survey-estimator-3d", "3D scalar OpenMP/MPI",
         "use data/random catalogs (iCatalogs=1,2), form D-alpha*R and random multipoles, then solve the window-coupling system"},
        {"encore-survey-estimator", "3D scalar OpenMP/MPI",
         "alias of survey-estimator-3d"},
        {"survey-edge-correction", "3D scalar OpenMP/MPI",
         "alias of survey-estimator-3d, not angular edge-corrections"},
        {"survey-keep-top-multipole", "3D scalar OpenMP/MPI",
         "also report the highest measured survey multipole instead of reserving it for window correction"},
        {"exclude-same-los", "3D scalar OpenMP/MPI",
         "exclude pivot/neighbor LOS_ID matches; does not exclude neighbor/neighbor LOS_ID matches"},
        {"exclude-los", "3D scalar OpenMP/MPI",
         "alias of exclude-same-los"},
        {"exclude-pivot-los", "3D scalar OpenMP/MPI",
         "alias of exclude-same-los"},
        {"KKCorrelation", "correlation",
         "select the convergence two-point correlation output"},
        {"KKKCorrelation", "correlation",
         "select the convergence three-point correlation output"},
        {"make-info", "startup",
         "print the compile-time Makefile configuration"},
        {"make-tree", "tree search",
         "build the tree and return without evaluating histograms"},
        {"mask-inside", "CFITSIO/NumPy input",
         "retain HEALPix pixels inside, rather than outside, the mask"},
        {"measure-cputime", "I/O",
         "measure and report CPU time spent reading data"},
        {"NNCorrelation", "correlation",
         "select number-count two-point correlation output"},
        {"NNEstimator", "correlation",
         "select the configured number-count correlation estimator"},
        {"NNLandySzalay1", "octree NN addon",
         "use the first Landy-Szalay number-count estimator"},
        {"NNLandySzalay2", "octree NN addon",
         "use the second Landy-Szalay number-count estimator"},
        {"NNStandard", "octree NN addon",
         "use the standard number-count estimator"},
        {"no-arfken", "coordinates",
         "use longitude/latitude ordering instead of the Arfken convention"},
        {"no-check-two-bodies-eq-pos", "tree search",
         "skip the tree-build check for bodies at identical positions"},
        {"no-normalize-HistZeta", "multipoles",
         "return raw three-point multipoles; KD-tree, balltree and BALLS4 exclude repeated neighbors in this mode"},
        {"no-one-ball", "tree search",
         "disable the one-ball aggregation path"},
        {"no-out-Hist", "output",
         "suppress histogram files"},
        {"no-stop", "startup",
         "continue after make-info, print-options, or print-search-methods"},
        {"no-two-ball", "balls addon",
         "disable the singular two-ball traversal path"},
        {"no-two-balls", "balls addon",
         "disable two-ball cell aggregation"},
        {"only-2pcf", "octree-ggg-omp and two-ball methods",
         "run only the 2PCF when both correlation orders are built"},
        {"only-2pcf-3d", "3D scalar addon",
         "compute only the scalar three-dimensional 2PCF"},
        {"only-3pcf", "TreeCorr-style two-ball methods",
         "run only the triple-node 3PCF when both correlation orders are built"},
        {"only-3pcf-3d", "3D scalar addon",
         "compute only the scalar three-dimensional 3PCF"},
        {"only-pos", "I/O addon",
         "read or write positions without convergence/shear fields"},
        {"out-HistZetaG", "output",
         "write the full three-point histogram grid"},
        {"out-m-HistZeta", "output",
         "write three-point histograms separated by multipole order"},
        {"patch", "sky input/test models",
         "restrict data to the configured angular patch"},
        {"patch-with-all", "multipoles",
         "combine patch selection with all-pivot edge-correction counting"},
        {"pivot-number", "GGG addon",
         "limit the search to the configured number of pivots"},
        {"plot-map-gif", "CFITSIO",
         "write the temporary HEALPix map used by the map-plot workflow"},
        {"pos-and-convergence", "I/O addon",
         "read position columns followed by convergence"},
        {"pos-and-convergence-weight", "I/O addon",
         "read position columns followed by convergence and weight"},
        {"pos-and-shear", "I/O addon",
         "read position columns followed by the two shear components"},
        {"post-processing", "workflow",
         "run posScript after the main computation"},
        {"pre-processing", "workflow",
         "run preScript before reading data"},
        {"print-options", "startup",
         "print this list of recognized options"},
        {"print-search-methods", "startup",
         "print searching methods registered in this executable and brief usage"},
        {"ra-dec-only", "CFITSIO",
         "ignore radial distance and place RA/DEC input on the unit sphere"},
        {"ra-reversed", "coordinates",
         "reverse right ascension as 2*pi minus RA"},
        {"random-point", "sky selection",
         "choose a random center for the selected sky region"},
        {"rbin-arcmin", "GGG output",
         "write radial-bin coordinates in arcminutes"},
        {"rbin-degree", "GGG output",
         "write radial-bin coordinates in degrees"},
        {"Rcut/theta", "tree search",
         "set the tree cutoff radius to rangeN/theta"},
        {"read-mask", "input/search",
         "use a supported companion or embedded mask, including set_catalog(mask=...); engine-specific support is listed by print-search-methods"},
        {"remove-mean", "tree search",
         "subtract the catalog mean convergence before tree construction"},
        {"rotation", "sky selection",
         "rotate a selected sky region to the reference center"},
        {"same-infiles", "input",
         "reuse the same input filename for multiple catalog entries"},
        {"save-ra-dec", "output",
         "preserve angular RA/DEC coordinates when saving converted data"},
        {"set-default-param", "smoothing addon",
         "replace smoothing parameters with addon defaults"},
        {"set-Nb-noSel", "smoothing addon",
         "set aggregate body counts without the selection pass"},
        {"smooth", "smoothing addon",
         "enable tree smoothing"},
        {"smooth-min-cell", "smoothing addon",
         "smooth using the minimum accepted cell scale"},
        {"smooth-pivot", "OpenMP tree addons",
         "opt into supported pivot smoothing; not supported by forest, exact scalar-3D, or octree-balls4 engines; compiling SMOOTHPIVOT does not activate it"},
        {"statistics-histograms", "workflow",
         "compute statistics from a set of histogram realizations"},
        {"stop", "workflow",
         "stop after the requested input, conversion, or preprocessing phase"},
        {"stop-fits", "CFITSIO",
         "stop after reading and reporting FITS metadata"},
        {"stop-numpy", "CFITSIO/NumPy input",
         "stop after reading the NumPy/HEALPix input"},
        {"sw94", "tree search",
         "use the Salmon-Warren 1994 cell acceptance criterion"},
        {"treecorr-direct-triples", "balltree-2balls-omp_3pcf",
         "use the slower genuine triple-node validation traversal; also supported by octree-2balls-omp, not its MPI variants"},
        {"treecorr-bucket-leaves", "balltree-2balls-omp",
         "use nsmooth bodies per leaf in the explicit triple-node algorithm"},
        {"treecorr-singleton-leaves", "balltree-2balls-omp_3pcf",
         "force one body per leaf instead of the nsmooth leaf capacity"},
        {"weights-norm", "correlation addons",
         "normalize correlations by pair weights instead of pair counts"},
        {"with-weight", "CFITSIO",
         "read and use the configured weight column"},
        {"xr1r2", "correlation output",
         "write the alternate x-r1-r2 three-point representation"}
    };
    size_t i;
    const size_t option_count =
        sizeof(registered_options) / sizeof(registered_options[0]);

    (void)gd;

    verb_print_zero(cmd->verbose,
                    "\nRegistered options (%zu):\n", option_count);
    verb_print_zero(cmd->verbose,
                    "Options are comma-separated and case-sensitive.\n");

    for (i = 0; i < option_count; ++i) {
        verb_print_zero(cmd->verbose, "\n- %s [%s]: %s.\n",
                        registered_options[i].name,
                        registered_options[i].scope,
                        registered_options[i].description);
    }

    verb_print_zero(cmd->verbose, "\n");

    return SUCCESS;
}
