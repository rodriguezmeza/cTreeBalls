/*==============================================================================
 MODULE: input.c                [cTreeBalls]
 Written by: Mario A. Rodriguez-Meza
 based in: input.c module by Julien Lesgourgues
 Starting date:    april 2023
 Purpose:
 Language: C
 Use:
 Major revisions:
 ==============================================================================*/
//        1          2          3          4        ^ 5          6          7

//
// lines where there is a "//B socket:" string are places to include module files
//  that can be found in addons/addons_include folder
//

#include "globaldefs.h"
#include "input.h"

typedef struct {
    struct cmdline_data *cmd;
    struct global_data *gd;
    struct file_content *pfc;
    char *errmsg;
} input_guard_context;

static int input_read_from_file_callback(void *argument)
{
    input_guard_context *context = argument;
    return input_read_from_file(context->cmd, context->gd,
                                context->pfc, context->errmsg);
}

int input_read_from_file_guarded(struct cmdline_data *cmd,
                                 struct global_data *gd,
                                 struct file_content *pfc,
                                 ErrorMsg errmsg)
{
    input_guard_context context;

    context.cmd = cmd;
    context.gd = gd;
    context.pfc = pfc;
    context.errmsg = errmsg;
    return cballs_allocation_guard(input_read_from_file_callback,
                                   &context, errmsg, _ERRORMSGSIZE_);
}

//B helper routine
static char *copy_param_string(const char *src)
{
    size_t len = strlen(src);
    char *dst = (char *) malloc((len + 1) * sizeof(char));
    if (dst != NULL)
        memcpy(dst, src, len + 1);
    return dst;
}

static char *copy_script_param_string(const char *src)
{
    size_t len;
    char *dst;

    len = strlen(src);

    if (len < 2 || src[0] != '"' || src[len - 1] != '"')
        return NULL;

    dst = (char *) malloc((len - 1) * sizeof(char));
    if (dst == NULL)
        return NULL;

    memcpy(dst, src + 1, len - 2);
    dst[len - 2] = '\0';

    return dst;
}
//E

local int testParameterFile(struct cmdline_data*,
                            struct global_data*,
                            char *,
                            ErrorMsg);

int input_find_file(struct  cmdline_data* cmd, struct  global_data* gd,
                    char *fname,
                    struct file_content * fc,
                    ErrorMsg errmsg)
{
    string routineName = "input_find_file";
    stream outstr = NULL;
  struct file_content fc_input;
  struct file_content fc_precision;
  struct file_content * pfc_input;
  struct file_content fc_setroot;

  int i;
  char extension[5];
  char input_file[_ARGUMENT_LENGTH_MAX_];
  char precision_file[_ARGUMENT_LENGTH_MAX_];

    //B test if cmd->paramfile exist...
    if (!strnull(fname)) {
        class_call(testParameterFile(cmd, gd, cmd->paramfile, errmsg),
                   errmsg,
                   errmsg);
    }
    //E

  pfc_input = &fc_input;

    fc->size = 0;
    fc->filename = NULL;
    fc->name = NULL;
    fc->value = NULL;
    fc->read = NULL;

    fc_input.size = 0;
    fc_input.filename = NULL;
    fc_input.name = NULL;
    fc_input.value = NULL;
    fc_input.read = NULL;

    fc_precision.size = 0;
    fc_precision.filename = NULL;
    fc_precision.name = NULL;
    fc_precision.value = NULL;
    fc_precision.read = NULL;

    fc_setroot.size = 0;
    fc_setroot.filename = NULL;
    fc_setroot.name = NULL;
    fc_setroot.value = NULL;
    fc_setroot.read = NULL;
    
    input_file[0]='\0';
    precision_file[0]='\0';

    if (strnull(fname)) {
        verb_print(1,
        "If you intend to use a parameter file call with <ParameterFileName>\n");
        return FAILURE;
    } else {
        //B
        if (stropen_checked(fname, "r", &outstr, errmsg, _ERRORMSGSIZE_) == FAILURE)
            return FAILURE;

        if (cballs_stream_close_checked(cmd, routineName, &outstr, fname) == FAILURE) {
            snprintf(errmsg, _ERRORMSGSIZE_, "%s", cmd->error_message);
            return FAILURE;
        }
        //E
    }

    if (copy_checked(input_file, sizeof(input_file), fname, "input file") != 0)
        return FAILURE;
    
    if (!strnull(input_file)) {
      class_call_except(parser_read_file(input_file, &fc_input, errmsg),
                        errmsg,
                        errmsg,
                        parser_free(&fc_input););

      class_call_except(input_set_root(input_file, &pfc_input, &fc_setroot, errmsg),
                        errmsg,
                        errmsg,
                        parser_free(pfc_input);
                        if (pfc_input != &fc_setroot)
                          parser_free(&fc_setroot);
                        parser_free(fc););
    }

    if ((input_file[0] != '\0') || (precision_file[0] != '\0')) {
      class_call_except(parser_cat(pfc_input, &fc_precision, fc, errmsg),
                        errmsg,
                        errmsg,
                        parser_free(pfc_input);
                        if (pfc_input != &fc_setroot)
                          parser_free(&fc_setroot);
                        parser_free(&fc_precision);
                        parser_free(fc););
    }

    class_call(parser_free(pfc_input), errmsg, errmsg);

    if (pfc_input != &fc_setroot)
      class_call(parser_free(&fc_setroot), errmsg, errmsg);

    class_call(parser_free(&fc_precision), errmsg, errmsg);
    
  return SUCCESS;
}

int input_set_root(char* input_file,
                   struct file_content** ppfc_input,
                   struct file_content * pfc_setroot,
                   ErrorMsg errmsg) {

    int flag1;
  int index_root_in_fc_input = -1;
  int overwrite_root;

  FileArg outfname;

  struct file_content fc_root;                      // Temporary structure with
                                                    //  only the root name
  FileArg string1;                                  // Is ignored
  struct file_content * pfc = *ppfc_input;

    //B
    fc_root.size = 0;
    fc_root.filename = NULL;
    fc_root.name = NULL;
    fc_root.value = NULL;
    fc_root.read = NULL;

    pfc_setroot->size = 0;
    pfc_setroot->filename = NULL;
    pfc_setroot->name = NULL;
    pfc_setroot->value = NULL;
    pfc_setroot->read = NULL;
    //E
    
    class_call_except(parser_read_string(pfc, "rootDir", &string1, &flag1, errmsg),
                      errmsg,
                      errmsg,
                      parser_free(&fc_root);
                      parser_free(pfc_setroot););

//B To behave as not class_lib parameter file
    overwrite_root = TRUE;
    class_read_flag("overwrite_root",overwrite_root);
    overwrite_root = TRUE;
//E

    if (flag1 == FALSE) {
        if (copy_checked(outfname, sizeof(outfname), "Output", "rootDir") != 0)
            return FAILURE;
    } else {
        for (index_root_in_fc_input=0; index_root_in_fc_input<pfc->size;
             ++index_root_in_fc_input) {
            if (strcmp(pfc->name[index_root_in_fc_input], "rootDir") == 0) {
                if (copy_checked(outfname, sizeof(outfname),
                                 pfc->value[index_root_in_fc_input],
                                 "rootDir") != 0)
                    return FAILURE;
                break;
            }
        }
    }

    if(flag1 == FALSE) {
        class_call(parser_init(&fc_root, 1, pfc->filename, errmsg),
                   errmsg,errmsg);
        if (copy_checked(fc_root.name[0], sizeof(FileArg), "rootDir", "rootDir name") != 0) {
            parser_free(&fc_root);
            return FAILURE;
        }

        if (copy_checked(fc_root.value[0], sizeof(FileArg),
                         outfname, "rootDir value") != 0) {
            parser_free(&fc_root);
            return FAILURE;
        }

        fc_root.read[0] = FALSE;

        //B
        class_call_except(parser_cat(pfc, &fc_root, pfc_setroot, errmsg),
                          errmsg,
                          errmsg,
                          parser_free(&fc_root);
                          parser_free(pfc_setroot););

        parser_free(pfc);
        parser_free(&fc_root);

        (*ppfc_input) = pfc_setroot;
        //E

    } else {
        if (copy_checked(pfc->value[index_root_in_fc_input], sizeof(FileArg),
                         outfname, "rootDir value") != 0)
            return FAILURE;

        (*ppfc_input) = pfc;
    }

  return SUCCESS;
}


int input_read_from_file(struct cmdline_data *cmd, struct  global_data* gd,
                         struct file_content * pfc,
                         ErrorMsg errmsg)
{
    int input_verbose = 0;

    if (gd->startrun_cputime==FALSE) {
        gd->cpuinit = CPUTIME;                       // init of cpu time
        gd->cpurealinit = rcpu_time();               // init of real time
    }

    gd->cmd_allocated = TRUE;

    class_read_int("verbose",input_verbose);
    verb_print(input_verbose, "\nReading input parameters...\n");

    class_call(input_read_parameters(cmd, gd, pfc, errmsg),errmsg,errmsg);

    return SUCCESS;
}

int input_read_parameters(struct cmdline_data *cmd,
                          struct  global_data* gd,
                          struct file_content * pfc,
                          ErrorMsg errmsg)
{
    int input_verbose=0;

    class_call(input_default_params(cmd),errmsg,errmsg);
    class_read_int("input_verbose",input_verbose);
    class_call(input_read_parameters_general(cmd, gd,  pfc,errmsg),errmsg,errmsg);

    return SUCCESS;
}

int input_read_parameters_general(struct cmdline_data *cmd,
                                  struct  global_data* gd,
                                  struct file_content * pfc, ErrorMsg errmsg)
{
#define FREE_CMD_STRING(flag, ptr)        \
    do {                                  \
        if ((flag) == TRUE && (ptr) != NULL) { \
            free(ptr);                    \
            (ptr) = NULL;                 \
            (flag) = FALSE;               \
        }                                 \
    } while (0)

#ifdef SAVERESTORE
#define SAVERESTORE_FREE_STRINGS_ON_FAILURE()              \
    do {                                                   \
        FREE_CMD_STRING(gd->statefileFlag, cmd->statefile); \
        FREE_CMD_STRING(gd->restorefileFlag, cmd->restorefile); \
    } while (0)
#else
#define SAVERESTORE_FREE_STRINGS_ON_FAILURE() do { } while (0)
#endif

#ifdef IOLIB
#define IOLIB_FREE_STRINGS_ON_FAILURE()                 \
    do {                                                \
        FREE_CMD_STRING(gd->columnsFlag, cmd->columns); \
    } while (0)
#else
#define IOLIB_FREE_STRINGS_ON_FAILURE() do { } while (0)
#endif
    
#define BASE_FREE_STRINGS_ON_FAILURE()                 \
    do {                                               \
        FREE_CMD_STRING(gd->searchMethodFlag, cmd->searchMethod); \
        FREE_CMD_STRING(gd->rsmoothFlagFree, cmd->rsmooth);       \
        FREE_CMD_STRING(gd->infileFlag, cmd->infile);             \
        FREE_CMD_STRING(gd->infilefmtFlag, cmd->infilefmt);       \
        FREE_CMD_STRING(gd->iCatalogsFlag, cmd->iCatalogs);       \
        FREE_CMD_STRING(gd->rootDirFlagFree, cmd->rootDir);       \
        FREE_CMD_STRING(gd->outfileFlag, cmd->outfile);           \
        FREE_CMD_STRING(gd->outfilefmtFlag, cmd->outfilefmt);     \
        FREE_CMD_STRING(gd->histNNFileNameFlag, cmd->histNNFileName); \
        FREE_CMD_STRING(gd->histXi2pcfFileNameFlag, cmd->histXi2pcfFileName); \
        FREE_CMD_STRING(gd->histZetaFileNameFlag, cmd->histZetaFileName); \
        FREE_CMD_STRING(gd->suffixOutFilesFlag, cmd->suffixOutFiles); \
        FREE_CMD_STRING(gd->testmodelFlag, cmd->testmodel);       \
        FREE_CMD_STRING(gd->preScriptFlag, cmd->preScript);       \
        FREE_CMD_STRING(gd->posScriptFlag, cmd->posScript);       \
        FREE_CMD_STRING(gd->optionsFlag, cmd->options);           \
        SAVERESTORE_FREE_STRINGS_ON_FAILURE();                    \
        IOLIB_FREE_STRINGS_ON_FAILURE();                \
    } while (0)

#define PARSER_READ(call) \
    class_call_except(call, errmsg, errmsg, BASE_FREE_STRINGS_ON_FAILURE();)

    string routineName = "input_read_parameters_general";
    int flag;
    int flag1,flag2;
    int param;
    int index;
    size_t slen;
    double param1,param2;
    char string1[_ARGUMENT_LENGTH_MAX_];

    // All malloc have to be freed at the end of the run (EndRun)

    //B Parameters related to the searching method
    class_call_except(parser_read_string(pfc, "searchMethod", &string1, &flag1, errmsg),
                      errmsg,
                      errmsg,
                      BASE_FREE_STRINGS_ON_FAILURE(););
    
    gd->searchMethodFlag=FALSE;
    if (flag1 == TRUE) {
        for (index=0;index<pfc->size;++index){
          if (strcmp(pfc->name[index],"searchMethod") == 0){
              cmd->searchMethod = copy_param_string(pfc->value[index]);
              if (cmd->searchMethod == NULL) {
                  BASE_FREE_STRINGS_ON_FAILURE();
                  return FAILURE;
              }
              gd->searchMethodFlag=TRUE;
            break;
          }
        }
    }
    
    PARSER_READ(parser_read_int(pfc, "mChebyshev", &param, &flag, errmsg));
    
    if (flag == TRUE)
      cmd->mChebyshev = param;

    PARSER_READ(parser_read_int(pfc, "nsmooth", &param, &flag, errmsg));
    
    if (flag == TRUE)
      cmd->nsmooth = param;

    PARSER_READ(parser_read_string(pfc, "rsmooth", &string1, &flag1, errmsg));

    gd->rsmoothFlagFree=FALSE;
    if (flag1 == TRUE) {
        for (index=0;index<pfc->size;++index){
          if (strcmp(pfc->name[index],"rsmooth") == 0){
              cmd->rsmooth = copy_param_string(pfc->value[index]);
              if (cmd->rsmooth == NULL) {
                  BASE_FREE_STRINGS_ON_FAILURE();
                  return FAILURE;
              }
              gd->rsmoothFlagFree = TRUE;
            break;
          }
        }
    }

    PARSER_READ(parser_read_double(pfc, "theta", &param1, &flag1, errmsg));

    if (flag1 == TRUE){
      cmd->theta = param1;
    }

    PARSER_READ(parser_read_string(pfc, "usePeriodic", &string1, &flag1, errmsg));

    if (flag1 == TRUE) {
        if (strchr("tTyY1", *string1) != NULL)
            cmd->usePeriodic=1;
        if (strchr("fFnN0", *string1) != NULL)
            cmd->usePeriodic=0;
    }
    //E

    //B Parameters about the I/O file(s)
    // Input catalog parameters
    class_call_except(parser_read_string(pfc, "infile",
                                         &string1, &flag1, errmsg),
                      errmsg, errmsg,
                      BASE_FREE_STRINGS_ON_FAILURE(););
    
    gd->infileFlag=FALSE;
    if (flag1 == TRUE) {
        for (index=0;index<pfc->size;++index){
          if (strcmp(pfc->name[index],"infile") == 0){
              cmd->infile = copy_param_string(pfc->value[index]);
              if (cmd->infile == NULL) {
                  BASE_FREE_STRINGS_ON_FAILURE();
                  return FAILURE;
              }
              gd->infileFlag = TRUE;
            break;
          }
        }
    }

    class_call_except(parser_read_string(pfc, "infileformat",
                                         &string1, &flag1, errmsg),
                      errmsg,
                      errmsg,
                      BASE_FREE_STRINGS_ON_FAILURE(););

    gd->infilefmtFlag=FALSE;
    if (flag1 == TRUE) {
        for (index=0;index<pfc->size;++index){
          if (strcmp(pfc->name[index],"infileformat") == 0){
              cmd->infilefmt = copy_param_string(pfc->value[index]);
              if (cmd->infilefmt == NULL) {
                  BASE_FREE_STRINGS_ON_FAILURE();
                  return FAILURE;
              }
              gd->infilefmtFlag = TRUE;
            break;
          }
        }
    }

    class_call_except(parser_read_string(pfc, "iCatalogs",
                                         &string1, &flag1, errmsg),
                      errmsg,
                      errmsg,
                      BASE_FREE_STRINGS_ON_FAILURE(););

    gd->iCatalogsFlag=FALSE;
    if (flag1 == TRUE) {
        for (index=0;index<pfc->size;++index){
          if (strcmp(pfc->name[index],"iCatalogs") == 0){
              cmd->iCatalogs = copy_param_string(pfc->value[index]);
              if (cmd->iCatalogs == NULL) {
                  BASE_FREE_STRINGS_ON_FAILURE();
                  return FAILURE;
              }
              gd->iCatalogsFlag = TRUE;
            break;
          }
        }
    }

    // Output parameters
    class_call_except(parser_read_string(pfc, "rootDir",
                                         &string1, &flag1, errmsg),
                      errmsg,
                      errmsg,
                      BASE_FREE_STRINGS_ON_FAILURE(););

    gd->rootDirFlagFree=FALSE;
    if (flag1 == TRUE) {
        for (index=0;index<pfc->size;++index){
          if (strcmp(pfc->name[index],"rootDir") == 0){
              cmd->rootDir = copy_param_string(pfc->value[index]);
              if (cmd->rootDir == NULL) {
                  BASE_FREE_STRINGS_ON_FAILURE();
                  return FAILURE;
              }
              gd->rootDirFlagFree = TRUE;
            break;
          }
        }
    }

    class_call_except(parser_read_string(pfc, "outfile",
                                         &string1, &flag1, errmsg),
                      errmsg,
                      errmsg,
                      BASE_FREE_STRINGS_ON_FAILURE(););

    gd->outfileFlag=FALSE;
    if (flag1 == TRUE) {
        for (index=0;index<pfc->size;++index){
          if (strcmp(pfc->name[index],"outfile") == 0){
              cmd->outfile = copy_param_string(pfc->value[index]);
              if (cmd->outfile == NULL) {
                  BASE_FREE_STRINGS_ON_FAILURE();
                  return FAILURE;
              }
              gd->outfileFlag = TRUE;
            break;
          }
        }
    }

    class_call_except(parser_read_string(pfc, "outfileformat",
                                         &string1, &flag1, errmsg),
                      errmsg,
                      errmsg,
                      BASE_FREE_STRINGS_ON_FAILURE(););

    gd->outfilefmtFlag=FALSE;
    if (flag1 == TRUE) {
        for (index=0;index<pfc->size;++index){
          if (strcmp(pfc->name[index],"outfileformat") == 0){
              cmd->outfilefmt = copy_param_string(pfc->value[index]);
              if (cmd->outfilefmt == NULL) {
                  BASE_FREE_STRINGS_ON_FAILURE();
                  return FAILURE;
              }
              gd->outfilefmtFlag=TRUE;
            break;
          }
        }
    }

    //B Parameters to set a region in the sky, for example for Takahashi data
    PARSER_READ(parser_read_double(pfc, "thetaL", &param1, &flag1, errmsg));
    if (flag1 == TRUE){
      cmd->thetaL = param1;
    }
    PARSER_READ(parser_read_double(pfc, "thetaR", &param1, &flag1, errmsg));
    if (flag1 == TRUE){
      cmd->thetaR = param1;
    }
    PARSER_READ(parser_read_double(pfc, "phiL", &param1, &flag1, errmsg));
    if (flag1 == TRUE){
      cmd->phiL = param1;
    }
    PARSER_READ(parser_read_double(pfc, "phiR", &param1, &flag1, errmsg));
    if (flag1 == TRUE){
      cmd->phiR = param1;
    }
    //E

    //B Parameters to control histograms and their output files
    class_call_except(parser_read_string(pfc, "useLogHist",
                                         &string1, &flag1, errmsg),
                      errmsg,
                      errmsg,
                      BASE_FREE_STRINGS_ON_FAILURE(););
    if (flag1 == TRUE) {
        if (strchr("tTyY1", *string1) != NULL)
            cmd->useLogHist=1;
        if (strchr("fFnN0", *string1) != NULL)
            cmd->useLogHist=0;
    }
    PARSER_READ(parser_read_int(pfc, "logHistBinsPD", &param, &flag, errmsg));
    if (flag == TRUE)
      cmd->logHistBinsPD = param;
    //
    PARSER_READ(parser_read_int(pfc, "sizeHistN", &param, &flag, errmsg));
    if (flag == TRUE)
      cmd->sizeHistN = param;
    PARSER_READ(parser_read_double(pfc, "rangeN", &param1, &flag1, errmsg));
    if (flag1 == TRUE){
      cmd->rangeN = param1;
    }
    PARSER_READ(parser_read_double(pfc, "rminHist", &param1, &flag1, errmsg));
    if (flag1 == TRUE){
      cmd->rminHist = param1;
    }

    PARSER_READ(parser_read_int(pfc, "sizeHistPhi", &param, &flag, errmsg));
    if (flag == TRUE)
      cmd->sizeHistPhi = param;
    //

    class_call_except(parser_read_string(pfc, "histNNFileName",
                                         &string1, &flag1, errmsg),
                      errmsg,
                      errmsg,
                      BASE_FREE_STRINGS_ON_FAILURE(););
    gd->histNNFileNameFlag=FALSE;
    if (flag1 == TRUE) {
        for (index=0;index<pfc->size;++index){
          if (strcmp(pfc->name[index],"histNNFileName") == 0){
              cmd->histNNFileName = copy_param_string(pfc->value[index]);
              if (cmd->histNNFileName == NULL) {
                  BASE_FREE_STRINGS_ON_FAILURE();
                  return FAILURE;
              }
              gd->histNNFileNameFlag=TRUE;
            break;
          }
        }
    }

    class_call_except(parser_read_string(pfc, "histXi2pcfFileName",
                                         &string1, &flag1, errmsg),
                      errmsg,
                      errmsg,
                      BASE_FREE_STRINGS_ON_FAILURE(););
    gd->histXi2pcfFileNameFlag=FALSE;
    if (flag1 == TRUE) {
        for (index=0;index<pfc->size;++index){
          if (strcmp(pfc->name[index],"histXi2pcfFileName") == 0){
              cmd->histXi2pcfFileName = copy_param_string(pfc->value[index]);
              if (cmd->histXi2pcfFileName == NULL) {
                  BASE_FREE_STRINGS_ON_FAILURE();
                  return FAILURE;
              }
              gd->histXi2pcfFileNameFlag=TRUE;
            break;
          }
        }
    }
    class_call_except(parser_read_string(pfc, "histZetaFileName",
                                         &string1, &flag1, errmsg),
                      errmsg,
                      errmsg,
                      BASE_FREE_STRINGS_ON_FAILURE(););
    gd->histZetaFileNameFlag=FALSE;
    if (flag1 == TRUE) {
        for (index=0;index<pfc->size;++index){
          if (strcmp(pfc->name[index],"histZetaFileName") == 0){
              cmd->histZetaFileName = copy_param_string(pfc->value[index]);
              if (cmd->histZetaFileName == NULL) {
                  BASE_FREE_STRINGS_ON_FAILURE();
                  return FAILURE;
              }
              gd->histZetaFileNameFlag=TRUE;
            break;
          }
        }
    }

    class_call_except(parser_read_string(pfc, "suffixOutFiles",
                                         &string1, &flag1, errmsg),
                      errmsg,
                      errmsg,
                      BASE_FREE_STRINGS_ON_FAILURE(););
    gd->suffixOutFilesFlag=FALSE;
    if (flag1 == TRUE) {
        for (index=0;index<pfc->size;++index){
          if (strcmp(pfc->name[index],"suffixOutFiles") == 0){
              cmd->suffixOutFiles = copy_param_string(pfc->value[index]);
              if (cmd->suffixOutFiles == NULL) {
                  BASE_FREE_STRINGS_ON_FAILURE();
                  return FAILURE;
              }
              gd->suffixOutFilesFlag=TRUE;
            break;
          }
        }
    }
    //E

    //B Set of parameters needed to construct a test model
    PARSER_READ(parser_read_int(pfc, "seed", &param, &flag, errmsg));
    if (flag == TRUE)
      cmd->seed = param;

    PARSER_READ(parser_read_string(pfc, "testmodel", &string1, &flag1, errmsg));
    gd->testmodelFlag=FALSE;
    if (flag1 == TRUE) {
        for (index=0;index<pfc->size;++index){
          if (strcmp(pfc->name[index],"testmodel") == 0){
              cmd->testmodel = copy_param_string(pfc->value[index]);
              if (cmd->testmodel == NULL) {
                  BASE_FREE_STRINGS_ON_FAILURE();
                  return FAILURE;
              }
              gd->testmodelFlag=TRUE;
            break;
          }
        }
    }

    PARSER_READ(parser_read_int(pfc, "nbody", &param, &flag, errmsg));
    if (flag == TRUE)
      cmd->nbody = param;
    PARSER_READ(parser_read_double(pfc, "lengthBox", &param1, &flag1, errmsg));
    if (flag1 == TRUE){
      cmd->lengthBox = param1;
    }
    //E

    //B Miscellaneous parameters
    class_call_except(parser_read_string(pfc, "preScript",
                                         &string1, &flag1, errmsg),
                      errmsg,
                      errmsg,
                      BASE_FREE_STRINGS_ON_FAILURE(););
    char *script1;
    char *script2;
    gd->preScriptFlag=FALSE;
    if (flag1 == TRUE) {
        for (index=0;index<pfc->size;++index){
          if (strcmp(pfc->name[index],"preScript") == 0){
              slen = strlen(pfc->value[index]);
              cmd->preScript = copy_script_param_string(pfc->value[index]);
              if (cmd->preScript == NULL) {
                  BASE_FREE_STRINGS_ON_FAILURE();
                  return FAILURE;
              }
              gd->preScriptFlag = TRUE;
              
            break;
          }
        }
    }

    class_call_except(parser_read_string(pfc, "posScript",
                                         &string1, &flag1, errmsg),
                      errmsg,
                      errmsg,
                      BASE_FREE_STRINGS_ON_FAILURE(););
    gd->posScriptFlag=FALSE;
    if (flag1 == TRUE) {
        for (index=0;index<pfc->size;++index){
          if (strcmp(pfc->name[index],"posScript") == 0){
              slen = strlen(pfc->value[index]);
              cmd->posScript = copy_script_param_string(pfc->value[index]);
              if (cmd->posScript == NULL) {
                  BASE_FREE_STRINGS_ON_FAILURE();
                  return FAILURE;
              }
              gd->posScriptFlag = TRUE;
            break;
          }
        }
    }

    PARSER_READ(parser_read_int(pfc, "stepState", &param, &flag, errmsg));
    if (flag == TRUE)
      cmd->stepState = param;

    PARSER_READ(parser_read_int(pfc, "verbose", &param, &flag, errmsg));
    if (flag == TRUE)
      cmd->verbose = param;
    PARSER_READ(parser_read_int(pfc, "verbose_log", &param, &flag, errmsg));
    if (flag == TRUE)
      cmd->verbose_log = param;

#ifdef OPENMPCODE
    PARSER_READ(parser_read_int(pfc, "numberThreads", &param, &flag, errmsg));
    if (flag == TRUE)
      cmd->numthreads = param;
#endif

    class_call_except(parser_read_string(pfc, "options",
                                         &string1, &flag1, errmsg),
                      errmsg,
                      errmsg,
                      BASE_FREE_STRINGS_ON_FAILURE(););
    gd->optionsFlag=FALSE;
    if (flag1 == TRUE) {
        for (index=0;index<pfc->size;++index){
          if (strcmp(pfc->name[index],"options") == 0){
              cmd->options = copy_param_string(pfc->value[index]);
              if (cmd->options == NULL) {
                  BASE_FREE_STRINGS_ON_FAILURE();
                  return FAILURE;
              }
              gd->optionsFlag=TRUE;
            break;
          }
        }
    }
    //E

    //B base_path is needed in cyballs.px... but not needed in C cballs code
    FileArg base_path_tmp;
    class_call_except(parser_read_string(pfc,"base_path",&base_path_tmp,&flag1,errmsg),
                      errmsg,
                      errmsg,
                      BASE_FREE_STRINGS_ON_FAILURE(););
    if (flag1 == TRUE) {
        
        if (copy_checked(cmd->base_path, sizeof(cmd->base_path),
                         base_path_tmp, "base_path") != 0) {
            BASE_FREE_STRINGS_ON_FAILURE();
            return FAILURE;
        }
        
    }
    //E

//B socket:
#ifdef ADDONS
#include "class_lib_include_01.h"
#endif
//

    PARSER_READ(parser_read_int(pfc, "sizeHistPhi", &param, &flag, errmsg));
    if (flag == TRUE)
      cmd->sizeHistPhi = param;

#undef PARSER_READ
#undef IOLIB_FREE_STRINGS_ON_FAILURE
#undef BASE_FREE_STRINGS_ON_FAILURE
#undef FREE_CMD_STRING
#undef SAVERESTORE_FREE_STRINGS_ON_FAILURE

  return SUCCESS;
}

//B cTreeBalls default values
#ifndef CMDLINE_DEFS_UNITSPHERE

int input_default_params(struct cmdline_data *cmd)
{
// Every item in cmdline_defs.h must have an item here::

    //B Parameters related to the searching method
    cmd->searchMethod = "octree-ggg-omp";
    cmd->mChebyshev = 7;
    cmd->nsmooth = 8;
    cmd->rsmooth = "\0";
    cmd->theta = 1.05;
    cmd->usePeriodic = 0;
    //E

    //B Parameters about the I/O file(s)
    // Input catalog parameters
    cmd->infile = "\0";
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
    cmd->sizeHistN = 40;
    cmd->rangeN = 100.0;
    cmd->rminHist = 1.0e-3;
    cmd->sizeHistPhi = 32;
    //
    cmd->histNNFileName = "histNN";
    cmd->histXi2pcfFileName = "histXi2pcf";
    cmd->histZetaFileName = "histZeta";
    cmd->suffixOutFiles = "";
    //E

    //B Set of parameters needed to construct a test model
    cmd->seed=123;                                          // to always have
                                                        // defaults Check in gsl
    cmd->testmodel = "simple-cubic-random";
    cmd->nbody = 16384;
    cmd->lengthBox = 10000.0;
    //E

    //B Miscellaneous parameters
    cmd->preScript = "";
    cmd->posScript = "";
    cmd->stepState = 10000;
    cmd->verbose = 0;
    cmd->verbose_log = 0;
#ifdef OPENMPCODE
    cmd->numthreads = 4;
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

#else
#include "input_default_params.h"
#endif

//E


//B parameter reading/testing from a file
local int testParameterFile(struct cmdline_data* cmd,
                            struct global_data* gd,
                            char *fname,
                            ErrorMsg errmsg)
{
// Every item in cmdline_defs.h must have an item here::
#define DOUBLE 1
#define STRING 2
#define INT 3
#define LONG 6
#define BOOLEAN 4
#define MAXTAGS 300
#define MAXCHARBUF 1024

    string routineName = "testParameterFile";
    FILE *fd;

  int  i,j,nt;
  int  id[MAXTAGS];
  void *addr[MAXTAGS];
  char tag[MAXTAGS][50];

    size_t str_size[MAXTAGS];

#undef SPName
#define SPName(param,paramtext,n)                                       \
  {SET_TAG_NAME(paramtext);                                             \
  addr[nt]=malloc(n);                                                   \
  if (addr[nt] == NULL) {                                               \
      snprintf(errmsg, _ERRORMSGSIZE_,                                  \
               "%s: not enough memory for parameter '%s'\n",            \
               routineName, paramtext);                                 \
      FREE_TEST_STRINGS();                                              \
      return FAILURE;                                                   \
  }                                                                     \
  str_size[nt]=(n);                                                     \
  id[nt++]=STRING;}

#define FREE_TEST_STRINGS()                                             \
  do {                                                                  \
      int free_i;                                                       \
      for (free_i = 0; free_i < nt; free_i++) {                         \
          if (id[free_i] == STRING && addr[free_i] != NULL) {           \
              free(addr[free_i]);                                       \
              addr[free_i] = NULL;                                      \
          }                                                             \
      }                                                                 \
  } while (0)

    int input_verbose = 2;
    verb_print(input_verbose, "\nparsing input parameters...\n");

    nt=0;

    //B Parameters related to the searching method
    SPName(cmd->searchMethod,"searchMethod",MAXLENGTHOFSTRSCMD);
    IPName(cmd->mChebyshev,"mChebyshev");
    IPName(cmd->nsmooth,"nsmooth");
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
    //B Parameters to set a region in the sky, for example for Takahashi data
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
        int has_line;
        int is_data;
        unsigned long line_number = 0;
//E

    if((fd=fopen(fname,"r"))) {
        while (TRUE) {
            if (read_parameter_line_checked(fd, line, sizeof(line),
                                            &has_line, &line_number, fname,
                                            errmsg,
                                            _ERRORMSGSIZE_) == FAILURE) {
                fclose(fd);
                FREE_TEST_STRINGS();
                return FAILURE;
            }
            if (has_line == FALSE)
                break;

            if (parse_parameter_line_checked(line,
                                             name, sizeof(name),
                                             value, sizeof(value),
                                             &is_data, fname, line_number,
                                             errmsg,
                                             _ERRORMSGSIZE_) == FAILURE) {
                fclose(fd);
                FREE_TEST_STRINGS();
                return FAILURE;
            }
            if (is_data == FALSE)
                continue;

            for(i=0,j=-1;i<nt;i++)
                if(strcmp(name,tag[i])==0) {
                    j=i;
                    tag[i][0]=0;
                    break;
                }
            if(j>=0) {
                switch(id[j]) {
                    case DOUBLE:
                        if (parse_double_checked(value, (double *)addr[j],
                                                 errmsg, _ERRORMSGSIZE_,
                                                 name) == FAILURE) {
                            fclose(fd);
                            FREE_TEST_STRINGS();
                            return FAILURE;
                        }
                        break;
                    case STRING:
                        if (strcmp(name,"preScript") == 0){ // To remove both '"'
                                size_t slen = strlen(value);
                            if (slen < 2 || value[0] != '"' || value[slen - 1] != '"') {
                                snprintf(errmsg, _ERRORMSGSIZE_,
                "preScript parameter needs enclosing script with \"\"!! (%s)\n\n",
                                value);
                                fclose(fd);
                                FREE_TEST_STRINGS();
                                return FAILURE;
                            }

                            if (slen - 2 >= str_size[j]) {
                                snprintf(errmsg, _ERRORMSGSIZE_,
                                         "preScript parameter too long!! (%s)\n\n", value);
                                        fclose(fd);
                                        FREE_TEST_STRINGS();
                                        return FAILURE;
                            }

                                memcpy(addr[j], value + 1, slen - 2);
                                ((char *)addr[j])[slen - 2] = '\0';
                        } else {
                            if (strcmp(name,"posScript") == 0){ // To remove both '"'
                                    size_t slen = strlen(value);

                                if (slen < 2 || value[0] != '"' || value[slen - 1] != '"') {
                                    snprintf(errmsg, _ERRORMSGSIZE_,
                "posScript parameter needs enclosing script with \"\"!! (%s)\n\n",
                                             value);
                                            fclose(fd);
                                            FREE_TEST_STRINGS();
                                            return FAILURE;
                                }

                                if (slen - 2 >= str_size[j]) {
                                    snprintf(errmsg, _ERRORMSGSIZE_,
                                             "posScript parameter too long!! (%s)\n\n",
                                             value);
                                    fclose(fd);
                                    FREE_TEST_STRINGS();
                                    return FAILURE;
                                }

                                    memcpy(addr[j], value + 1, slen - 2);
                                    ((char *)addr[j])[slen - 2] = '\0';
                            } else {
                                if (copy_checked((char *)addr[j], str_size[j], value, name) != 0) {
                                    snprintf(errmsg, _ERRORMSGSIZE_,
                                    "%s: string parameter '%s' too long\n",
                                             routineName, name);
                                    fclose(fd);
                                    FREE_TEST_STRINGS();
                                    return FAILURE;
                                }
                            }
                        }
                        break;
                    case INT:
                        if (parse_int_checked(value, (int *)addr[j],
                                              errmsg, _ERRORMSGSIZE_,
                                              name) == FAILURE) {
                            fclose(fd);
                            FREE_TEST_STRINGS();
                            return FAILURE;
                        }
                        break;
                    case LONG:
                        if (parse_long_checked(value, (long *)addr[j],
                                               errmsg, _ERRORMSGSIZE_,
                                               name) == FAILURE) {
                            fclose(fd);
                            FREE_TEST_STRINGS();
                            return FAILURE;
                        }
                        break;
                    case BOOLEAN:
                        if (parse_bool_checked(value, (bool *)addr[j],
                                               errmsg, _ERRORMSGSIZE_,
                                               name) == FAILURE) {
                            fclose(fd);
                            FREE_TEST_STRINGS();
                            return FAILURE;
                        }
                        break;
                }
            } else {
                snprintf(errmsg, _ERRORMSGSIZE_,
                         "%s:%lu: parameter '%s' is unknown or duplicated",
                         fname, line_number, name);
                fclose(fd);
                FREE_TEST_STRINGS();
                return FAILURE;
            }
        } // ! while loop
        if (fclose(fd) != 0) {
            snprintf(errmsg, _ERRORMSGSIZE_,
                     "%s: could not close parameter file '%s'",
                     routineName, fname);
            FREE_TEST_STRINGS();
            return FAILURE;
        }
    } else {
        snprintf(errmsg, _ERRORMSGSIZE_, "Cannot open parameter file %s\n", fname);
        FREE_TEST_STRINGS();
        return FAILURE;
    }

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
                        snprintf(errmsg, _ERRORMSGSIZE_,
                                 "%s: default string parameter '%s' too long\n",
                                 routineName, tag[i]);
                        FREE_TEST_STRINGS();
                        return FAILURE;
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

    FREE_TEST_STRINGS();

    return SUCCESS;

#undef DOUBLE
#undef STRING
#undef INT
#undef BOOLEAN
#undef MAXTAGS
#undef MAXCHARBUF
#undef FREE_TEST_STRINGS
#undef SPName
}
//E
